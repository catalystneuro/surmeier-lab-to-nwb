"""Segmentation extractor for Bruker Prairie View line scan recordings.

Bruker line-scan (``<Sequence type="Linescan">``) acquisitions repeatedly scan a single
user-drawn line through the imaging field, producing a 2D kymograph per channel
(``<folder>_Cycle{N}_Ch{C}_000001.ome.tif``: rows = line repetitions, columns = pixels
along the line) plus a processed CSV ``<folder>_Cycle{N}_LineProfileData.csv`` with the
per-channel 1D fluorescence profile averaged along the line versus time.

This extractor reads the LineProfileData CSV and the XML ``<PVLinescanDefinition>``
element to expose the line as a single pixel-mask ROI with its time trace. The raw
kymograph TIFF is not loaded here; it is best handled by a separate imaging extractor.

See https://github.com/catalystneuro/surmeier-lab-to-nwb obsidian docs
``source_formats/bruker/line_scan_format.md`` for the format reference.
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

# Channel index in the XML (1-based) maps to "Prof N" column in the CSV.
_PROFILE_COLUMN_FOR_CHANNEL = {"Ch1": (" Prof 1", "Prof 1 Time(ms)"), "Ch2": (" Prof 2", " Prof 2 Time(ms)")}


def _find_csv_column(available_columns: list[str], template: str) -> str:
    """Match a column name tolerating leading/trailing whitespace differences.

    Prairie View's LineProfileData CSV exports have inconsistent header whitespace
    across versions (e.g. ``" Prof 2"`` vs ``"Prof 2"``).
    """
    template_stripped = template.strip()
    for col in available_columns:
        if col == template or col.strip() == template_stripped:
            return col
    raise ValueError(
        f"Could not find column matching {template!r} (stripped: {template_stripped!r}). "
        f"Available: {available_columns}"
    )


class BrukerLineScanSegmentationExtractor(SegmentationExtractor):
    """A segmentation extractor for one channel of a Bruker line-scan recording.

    The line itself is the single ROI. Its trace is the averaged fluorescence profile
    along the line at each time point, read from ``LineProfileData.csv``.

    Parameters
    ----------
    folder_path : PathType
        Folder containing a Prairie View line-scan recording with a master
        ``<folder>.xml`` and a ``<folder>_Cycle{N}_LineProfileData.csv``.
    cycle : int, default: 1
        Cycle index to read.
    channel_name : str, default: "Ch1"
        Prairie View channel name to extract. Must be one of ``"Ch1"`` or ``"Ch2"``
        (the channels Prairie View emits per-channel ``Prof N`` columns for).
    frame_shape : tuple of int, optional
        ``(height, width)`` of the source frame the scan line is drawn in. Defaults
        to the ``SourcePixelsPerLine``/``SourceLinesPerFrame`` attributes on
        ``<PVLinescanDefinition>``.
    """

    extractor_name = "BrukerLineScanSegmentation"
    mode = "folder"

    def __init__(
        self,
        folder_path: PathType,
        cycle: int = 1,
        channel_name: str = "Ch1",
        frame_shape: tuple[int, int] | None = None,
    ):
        SegmentationExtractor.__init__(self)
        self.folder_path = Path(folder_path)
        self.cycle = cycle
        self.channel_name = channel_name

        if channel_name not in _PROFILE_COLUMN_FOR_CHANNEL:
            raise ValueError(
                f"channel_name={channel_name!r} not supported. " f"Available: {sorted(_PROFILE_COLUMN_FOR_CHANNEL)}"
            )

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

        sequence_type = sequence_metadata.get("@type", "")
        if sequence_type.lower() != "linescan":
            raise ValueError(f"Sequence in {xml_path} has type {sequence_type!r}; expected 'Linescan'.")

        line_definition = sequence_metadata["PVLinescanDefinition"]
        line_metadata = line_definition["Line"]
        self._line_start_x = int(line_metadata["@startPixelX"])
        self._line_start_y = int(line_metadata["@startPixelY"])
        self._line_stop_x = int(line_metadata.get("@stopPixelX", line_metadata["@startPixelX"]))
        self._line_stop_y = int(line_metadata.get("@stopPixelY", line_metadata["@startPixelY"]))
        self._line_length_microns = float(line_metadata.get("@lineLength", 0.0))

        if frame_shape is None:
            source_h = int(line_definition.get("@SourceLinesPerFrame", 512))
            source_w = int(line_definition.get("@SourcePixelsPerLine", 512))
            self._frame_shape = (source_h, source_w)
        else:
            self._frame_shape = tuple(frame_shape)

        # Load profile CSV
        profile_csv_filename = line_definition["LineScanProfiles"]["@DataFile"]
        profile_csv_path = self.folder_path / profile_csv_filename
        if not profile_csv_path.is_file():
            raise FileNotFoundError(f"Line scan profile CSV not found at '{profile_csv_path}'.")
        self._profile_csv_path = profile_csv_path
        profile_df = pd.read_csv(profile_csv_path)
        available_cols = list(profile_df.columns)
        profile_template, time_template = _PROFILE_COLUMN_FOR_CHANNEL[channel_name]
        profile_col = _find_csv_column(available_cols, profile_template)
        time_col = _find_csv_column(available_cols, time_template)

        # Drop NaN rows that show up at the tail of some CSV exports
        valid_mask = profile_df[profile_col].notna() & profile_df[time_col].notna()
        profile_data = profile_df.loc[valid_mask, profile_col].to_numpy(dtype=float)
        timestamps_ms = profile_df.loc[valid_mask, time_col].to_numpy(dtype=float)
        self._timestamps = timestamps_ms / 1000.0  # seconds

        # Single ROI: the line, identified by channel
        roi_id = f"{channel_name}_Line"
        self._roi_ids = [roi_id]
        # Shape (num_samples, num_rois=1)
        traces = profile_data.reshape(-1, 1)
        self._roi_responses.append(_RoiResponse("raw", traces, self._roi_ids))

        # Pixel mask: line endpoints' integer pixels along the line, weight=1.0 each
        height, width = self._frame_shape
        if self._line_start_y == self._line_stop_y:
            xs = np.arange(min(self._line_start_x, self._line_stop_x), max(self._line_start_x, self._line_stop_x) + 1)
            ys = np.full_like(xs, self._line_start_y)
        elif self._line_start_x == self._line_stop_x:
            ys = np.arange(min(self._line_start_y, self._line_stop_y), max(self._line_start_y, self._line_stop_y) + 1)
            xs = np.full_like(ys, self._line_start_x)
        else:
            # Diagonal — Bresenham via linspace + rounding
            num_steps = (
                max(abs(self._line_stop_x - self._line_start_x), abs(self._line_stop_y - self._line_start_y)) + 1
            )
            xs = np.round(np.linspace(self._line_start_x, self._line_stop_x, num_steps)).astype(int)
            ys = np.round(np.linspace(self._line_start_y, self._line_stop_y, num_steps)).astype(int)

        image_mask = np.zeros((height, width, 1), dtype=bool)
        ys_clipped = np.clip(ys, 0, height - 1)
        xs_clipped = np.clip(xs, 0, width - 1)
        image_mask[ys_clipped, xs_clipped, 0] = True
        self._roi_masks = _ROIMasks(
            data=image_mask,
            mask_tpe="nwb-image_mask",
            field_of_view_shape=(height, width),
            roi_id_map={roi_id: 0},
        )

        # Sampling frequency: derive from line scan period if regular
        if len(self._timestamps) > 1:
            intervals = np.diff(self._timestamps)
            mean_interval = float(np.mean(intervals))
            if mean_interval > 0 and np.allclose(intervals, mean_interval, rtol=1e-3):
                self._sampling_frequency = 1.0 / mean_interval
            else:
                self._times = self._timestamps

    def get_session_start_time(self, tz: ZoneInfo | None = None) -> datetime:
        """Return ``<PVScan @date>`` parsed to a datetime, optionally tz-aware."""
        date_str = self._xml_metadata["PVScan"]["@date"]
        dt = datetime.strptime(date_str, "%m/%d/%Y %I:%M:%S %p")
        if tz is not None:
            dt = dt.replace(tzinfo=tz)
        return dt

    def get_line_length_microns(self) -> float:
        """Return the physical line length in microns from ``<Line lineLength=...>``."""
        return self._line_length_microns

    def get_native_timestamps(self, start_sample: int | None = None, end_sample: int | None = None) -> np.ndarray:
        """Return per-sample timestamps in seconds from ``Prof N Time(ms)`` column."""
        start = start_sample or 0
        end = end_sample if end_sample is not None else len(self._timestamps)
        return self._timestamps[start:end]

    def get_frame_shape(self) -> tuple[int, int]:
        """Return ``(height, width)`` of the source frame the scan line was drawn in."""
        return self._frame_shape
