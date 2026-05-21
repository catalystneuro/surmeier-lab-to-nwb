"""Segmentation extractor for Bruker Prairie View Brightness Over Time (BOT) recordings.

Bruker BOT acquisitions produce two parallel data streams: full-frame OME-TIFFs (one per
channel per frame) and a per-region brightness CSV (``<dataset>_Cycle{N}-botData.csv``)
computed live by Prairie View. The regions are user-defined rectangles in the imaging
field, declared in the master ``<dataset>.xml`` under a ``<PVBOTs>`` element.

This extractor reads the CSV and XML and exposes the regions as ROIs with rectangular
image masks. The OME-TIFF imaging frames are not loaded; use ``BrukerTiffImagingExtractor``
on the same folder if you need them.

See https://github.com/catalystneuro/surmeier-lab-to-nwb obsidian docs
``source_formats/bruker/bot_format.md`` for the format reference.
"""

from datetime import datetime
from pathlib import Path
from zoneinfo import ZoneInfo

import numpy as np
import pandas as pd
import xmltodict

# Imports use the absolute `roiextractors` path while this module lives in
# surmeier-lab-to-nwb. When this is lifted into
# `roiextractors/extractors/tiffimagingextractors/`, switch to:
#   from ...extraction_tools import PathType
#   from ...segmentationextractor import SegmentationExtractor, _ROIMasks, _RoiResponse
from roiextractors.extraction_tools import PathType
from roiextractors.segmentationextractor import (
    SegmentationExtractor,
    _ROIMasks,
    _RoiResponse,
)


class BrukerBOTSegmentationExtractor(SegmentationExtractor):
    """A segmentation extractor for Bruker Prairie View BOT (Brightness Over Time) data.

    Each ``<Region>`` declared in the master XML becomes one ROI. ROI traces come from the
    matching ``Region N`` column of the BOT CSV. Spatial masks are rectangular image masks
    derived from the region geometry (x, y, width, height).

    Parameters
    ----------
    folder_path : PathType
        Path to a folder containing a Prairie View BOT recording. The folder must contain
        a ``<folder_name>.xml`` master file and at least one
        ``<folder_name>_Cycle{NNNNN}-botData.csv`` companion CSV.
    cycle : int, default: 1
        Cycle index to read. Most BOT recordings have a single cycle.
    frame_shape : tuple of int, optional
        ``(height, width)`` of the imaging frame. If not provided, inferred from the
        ``linesPerFrame``/``pixelsPerLine`` ``PVStateValue`` entries, falling back to
        ``(512, 512)`` if those are not present.
    channel_name : str, optional
        Restrict the segmentation to regions with this ``channelName`` (e.g. ``"Ch2"``).
        When ``None``, all regions across all channels are exposed as ROIs.

    Notes
    -----
    Bruker BOT regions can be replicated across channels (the same spatial rectangle is
    monitored independently per channel). When ``channel_name`` is ``None``, each replica
    becomes its own ROI. Use ``channel_name`` to filter to one channel's regions.
    """

    extractor_name = "BrukerBOTSegmentation"
    mode = "folder"

    def __init__(
        self,
        folder_path: PathType,
        cycle: int = 1,
        frame_shape: tuple[int, int] | None = None,
        channel_name: str | None = None,
    ):
        SegmentationExtractor.__init__(self)
        self.folder_path = Path(folder_path)
        self.cycle = cycle
        self.channel_name = channel_name

        xml_path = self.folder_path / f"{self.folder_path.name}.xml"
        if not xml_path.is_file():
            raise FileNotFoundError(f"Bruker XML configuration file not found at '{xml_path}'.")
        with open(xml_path, "r", encoding="utf-8") as f:
            self._xml_metadata = xmltodict.parse(f.read())

        sequence_metadata = self._xml_metadata["PVScan"]["Sequence"]
        if isinstance(sequence_metadata, list):
            try:
                sequence_metadata = next(s for s in sequence_metadata if int(s.get("@cycle", 0)) == cycle)
            except StopIteration:
                raise ValueError(f"Cycle {cycle} not found in {xml_path}.") from None

        if sequence_metadata.get("@type", "").lower() not in {
            "brightnessovertime",
            "tseries brightness over time element",
        }:
            raise ValueError(
                f"Sequence in {xml_path} has type {sequence_metadata.get('@type')!r}; "
                "expected 'BrightnessOverTime' or 'TSeries Brightness Over Time Element'."
            )

        # Map channel index ("2", "3", ...) → channel name ("Ch2", "Dodt", ...)
        first_frame = sequence_metadata["Frame"]
        if isinstance(first_frame, list):
            first_frame = first_frame[0]
        file_metadata = first_frame["File"]
        if not isinstance(file_metadata, list):
            file_metadata = [file_metadata]
        channel_index_to_name = {entry["@channel"]: entry["@channelName"] for entry in file_metadata}
        self._available_channels = sorted(set(channel_index_to_name.values()))
        self._channel_names = list(self._available_channels)

        if channel_name is not None and channel_name not in self._available_channels:
            raise ValueError(f"channel_name={channel_name!r} not in available channels {self._available_channels}.")

        # Parse <PVBOTs>/<Region>
        bot_metadata = sequence_metadata["PVBOTs"]
        bot_csv_filename = bot_metadata["@botData"]
        regions = bot_metadata["Region"]
        if not isinstance(regions, list):
            regions = [regions]

        # Each region's CSV column is "Region N" using 1-based XML order.
        all_regions: list[dict] = []
        for index, region in enumerate(regions, start=1):
            region_channel = channel_index_to_name[region["@channel"]]
            all_regions.append(
                {
                    "csv_column": f"Region {index}",
                    "channel_name": region_channel,
                    "x": int(region["@x"]),
                    "y": int(region["@y"]),
                    "width": int(region["@width"]),
                    "height": int(region["@height"]),
                }
            )

        if channel_name is not None:
            self._regions = [r for r in all_regions if r["channel_name"] == channel_name]
        else:
            self._regions = all_regions

        if not self._regions:
            raise ValueError(
                f"No regions found in BOT XML for channel_name={channel_name!r}. "
                f"Available regions: {[(r['channel_name'], r['csv_column']) for r in all_regions]}"
            )

        # Frame shape (height, width)
        if frame_shape is not None:
            self._frame_shape = tuple(frame_shape)
        else:
            self._frame_shape = self._infer_frame_shape()

        # Load BOT CSV traces
        bot_csv_path = self.folder_path / bot_csv_filename
        if not bot_csv_path.is_file():
            raise FileNotFoundError(f"BOT CSV not found at '{bot_csv_path}'.")
        self._bot_csv_path = bot_csv_path
        bot_df = pd.read_csv(bot_csv_path)
        if "Timestamp" not in bot_df.columns:
            raise ValueError(f"'Timestamp' column not found in {bot_csv_path}.")
        self._timestamps = bot_df["Timestamp"].to_numpy(dtype=float)

        missing_cols = [r["csv_column"] for r in self._regions if r["csv_column"] not in bot_df.columns]
        if missing_cols:
            raise ValueError(
                f"BOT CSV {bot_csv_path} is missing expected region columns: {missing_cols}. "
                f"Available columns: {list(bot_df.columns)}"
            )

        traces = np.stack(
            [bot_df[r["csv_column"]].to_numpy(dtype=float) for r in self._regions],
            axis=1,
        )  # shape: (num_samples, num_rois)

        # ROI IDs: use the channel + region label for uniqueness
        roi_ids = [f"{r['channel_name']}_{r['csv_column'].replace(' ', '')}" for r in self._regions]
        self._roi_ids = roi_ids
        self._roi_responses.append(_RoiResponse("raw", traces, roi_ids))

        # Image masks: dense bool array (height, width, num_rois)
        height, width = self._frame_shape
        image_masks = np.zeros((height, width, len(self._regions)), dtype=bool)
        for index, region in enumerate(self._regions):
            x, y, w, h = region["x"], region["y"], region["width"], region["height"]
            image_masks[y : y + h, x : x + w, index] = True
        self._roi_masks = _ROIMasks(
            data=image_masks,
            mask_tpe="nwb-image_mask",
            field_of_view_shape=(height, width),
            roi_id_map={roi_id: index for index, roi_id in enumerate(roi_ids)},
        )

        # Sampling frequency from CSV timestamps (CSV is the canonical timing source)
        if len(self._timestamps) > 1:
            intervals = np.diff(self._timestamps)
            mean_interval = float(np.mean(intervals))
            if mean_interval > 0 and np.allclose(intervals, mean_interval, rtol=1e-3):
                self._sampling_frequency = 1.0 / mean_interval
            else:
                self._times = self._timestamps

    def _infer_frame_shape(self) -> tuple[int, int]:
        """Read ``linesPerFrame`` / ``pixelsPerLine`` from PVStateShard, falling back to (512, 512)."""
        try:
            state_values = self._xml_metadata["PVScan"]["PVStateShard"]["PVStateValue"]
        except KeyError:
            return (512, 512)
        lines = next((s for s in state_values if s.get("@key") == "linesPerFrame"), None)
        pixels = next((s for s in state_values if s.get("@key") == "pixelsPerLine"), None)
        if lines is None or pixels is None:
            return (512, 512)
        return (int(lines["@value"]), int(pixels["@value"]))

    def get_session_start_time(self, tz: ZoneInfo | None = None) -> datetime:
        """Return ``<PVScan @date>`` parsed to a tz-aware datetime.

        The Bruker XML records local time at the rig with no timezone marker. Pass ``tz``
        to attach a timezone; defaults to a naive datetime when omitted.
        """
        date_str = self._xml_metadata["PVScan"]["@date"]
        dt = datetime.strptime(date_str, "%m/%d/%Y %I:%M:%S %p")
        if tz is not None:
            dt = dt.replace(tzinfo=tz)
        return dt

    def get_available_channels(self) -> list[str]:
        """Return the channel names declared in this BOT recording's XML."""
        return list(self._available_channels)

    def get_native_timestamps(self, start_sample: int | None = None, end_sample: int | None = None) -> np.ndarray:
        """Return per-sample timestamps from the BOT CSV ``Timestamp`` column."""
        start = start_sample or 0
        end = end_sample if end_sample is not None else len(self._timestamps)
        return self._timestamps[start:end]

    def get_frame_shape(self) -> tuple[int, int]:
        """Return ``(height, width)`` of the imaging frame the regions are defined in."""
        return self._frame_shape
