"""Segmentation extractor for Bruker PrairieView line scans.

A PrairieView line scan (``<Sequence type="Linescan">``) scans one line, drawn on a snapshot of the field of view,
over and over. Each trial is its own folder, holding the master ``<folder>.xml``, one kymograph TIFF per channel
(one row per scan, one column per pixel along the line), a snapshot of the field of view per channel
(``...Source.tif``) and a ``LineProfileData`` CSV. For each channel the CSV holds the time of every scan and the
mean of the kymograph over pixels ``x0`` to ``x1`` of the line, the window of the channel's ``<Profile>``.

This extractor exposes those pixels as one ROI and the CSV as its trace, for a list of trials of the same line
placed on one timeline by the start of each trial. The kymographs and snapshots are read by the interface.

See the vault note ``source_formats/bruker/line_scan_format.md`` for the format reference.
"""

from dataclasses import dataclass
from datetime import datetime, time
from pathlib import Path

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


@dataclass(frozen=True)
class LineScanTrial:
    """One trial of a line scan, read from the master XML of its folder."""

    name: str  # the folder name, which PrairieView also puts at the start of every file in it
    start_time: datetime  # the date of `PVScan` with the time of day of its `Sequence`, which has the microseconds
    line: tuple[int, int, int, int]  # start x, start y, stop x and stop y, in pixels of the snapshot
    profile_pixels: tuple[int, int]  # x0 and x1, the first and last pixel along the line that the profile averages
    frame_shape: tuple[int, int]  # height and width of the snapshot
    scan_line_period: float  # seconds
    num_scans: int
    profile_file_path: Path
    profile_column: str  # "Prof N", the CSV column of the channel's profile
    kymograph_file_path: Path
    source_image_file_path: Path

    @property
    def num_pixels(self) -> int:
        """The number of pixels along the line, which is the width of the kymograph."""
        start_x, start_y, stop_x, stop_y = self.line
        return max(abs(stop_x - start_x), abs(stop_y - start_y)) + 1


def _as_list(value) -> list:
    return value if isinstance(value, list) else [value]


def _read_master_xml(folder_path: Path) -> dict:
    """Return the ``PVScan`` element of a trial folder's master XML, which must hold one line scan."""
    xml_path = folder_path / f"{folder_path.name}.xml"
    pv_scan = xmltodict.parse(xml_path.read_text(encoding="utf-8-sig"))["PVScan"]
    sequence = pv_scan["Sequence"]
    if isinstance(sequence, list) or sequence["@type"] != "Linescan":
        raise ValueError(f"'{xml_path}' does not hold one line scan.")
    return pv_scan


def _read_scan_line(pv_scan: dict) -> tuple[int, int, int, int]:
    definition = pv_scan["Sequence"]["PVLinescanDefinition"]
    if definition["@mode"] != "straightLine":
        raise ValueError(f"Only straight lines are read, and this line scan is a {definition['@mode']!r}.")
    line = definition["Line"]
    start_x, start_y = int(line["@startPixelX"]), int(line["@startPixelY"])
    # PrairieView leaves out a stop coordinate the line shares with its start, as stopPixelY on a horizontal line
    return start_x, start_y, int(line.get("@stopPixelX", start_x)), int(line.get("@stopPixelY", start_y))


def _read_line_scan_trial(folder_path: Path, channel_name: str) -> LineScanTrial:
    pv_scan = _read_master_xml(folder_path)
    sequence = pv_scan["Sequence"]
    definition = sequence["PVLinescanDefinition"]
    frames = _as_list(sequence["Frame"])
    if len(frames) != 1:
        raise ValueError(f"The line scan in '{folder_path}' has {len(frames)} frames, where a trial has one.")
    files = _as_list(frames[0]["File"])
    channel_file = next((entry for entry in files if entry["@channelName"] == channel_name), None)
    if channel_file is None:
        channel_names = [entry["@channelName"] for entry in files]
        raise ValueError(f"The line scan in '{folder_path}' has no channel {channel_name!r}, only {channel_names}.")

    profiles = [
        profile
        for profile in _as_list(definition["LineScanProfiles"]["Profile"])
        if profile["@channel"] == channel_file["@channel"]
    ]
    windows = {(int(profile["@x0"]), int(profile["@x1"])) for profile in profiles}
    if len(windows) != 1:
        raise ValueError(f"The line scan in '{folder_path}' has {len(windows)} profile windows for {channel_name}.")

    pv_scan_date = datetime.strptime(pv_scan["@date"], "%m/%d/%Y %I:%M:%S %p")
    # `time.fromisoformat` takes up to six decimals and PrairieView writes seven
    start_time = datetime.combine(pv_scan_date.date(), time.fromisoformat(sequence["@time"][:15]))
    if not 0 <= (start_time - pv_scan_date).total_seconds() < 2:
        raise ValueError(f"The line scan in '{folder_path}' starts at {start_time}, but its date is {pv_scan_date}.")

    # Settings with one value per laser or detector hold them in child elements instead of a value attribute
    frame_state = {
        value["@key"]: value["@value"]
        for value in _as_list(frames[0]["PVStateShard"]["PVStateValue"])
        if "@value" in value
    }
    return LineScanTrial(
        name=folder_path.name,
        start_time=start_time,
        line=_read_scan_line(pv_scan),
        profile_pixels=windows.pop(),
        frame_shape=(int(definition["@SourceLinesPerFrame"]), int(definition["@SourcePixelsPerLine"])),
        scan_line_period=float(frame_state["scanLinePeriod"]),
        num_scans=int(frame_state["linesPerFrame"]),
        profile_file_path=folder_path / definition["LineScanProfiles"]["@DataFile"],
        profile_column=f"Prof {profiles[0]['@number']}",
        kymograph_file_path=folder_path / channel_file["@filename"],
        source_image_file_path=folder_path / channel_file["@source"],
    )


class BrukerLineScanSegmentationExtractor(SegmentationExtractor):
    """A segmentation extractor for one channel of a Bruker line scan, over the trials of one scan line.

    The ROI is the pixels of the line that PrairieView averages into its line profile, and its trace is that
    profile, read from the ``LineProfileData`` CSV of each trial. The trials are put in the order they were
    recorded, and the time of each scan is the start of its trial, relative to the first trial, plus its time
    in the CSV.

    Parameters
    ----------
    folder_paths : list of PathType
        The trial folders, each with its master XML and ``LineProfileData`` CSV. Their scan line, profile window
        and snapshot size must be the same.
    channel_name : str, default: "Ch1"
        PrairieView channel, as the XML names it.
    """

    extractor_name = "BrukerLineScanSegmentation"

    def __init__(self, folder_paths: list[PathType], channel_name: str = "Ch1"):
        SegmentationExtractor.__init__(self)
        if not folder_paths:
            raise ValueError("folder_paths is empty; give at least one trial folder.")
        trials = sorted(
            (_read_line_scan_trial(Path(folder_path), channel_name) for folder_path in folder_paths),
            key=lambda trial: trial.start_time,
        )
        first = trials[0]
        for trial in trials[1:]:
            for field in ("line", "profile_pixels", "frame_shape"):
                if getattr(trial, field) != getattr(first, field):
                    raise ValueError(
                        f"Trial '{trial.name}' has {field} {getattr(trial, field)} and trial '{first.name}' has "
                        f"{getattr(first, field)}. The trials of one line scan share their line, profile window and "
                        "snapshot; group the folders by scan line."
                    )
        self._trials = trials

        profiles, timestamps = [], []
        for trial in trials:
            profile = pd.read_csv(trial.profile_file_path)
            profile.columns = profile.columns.str.strip()
            values = profile[trial.profile_column].to_numpy(dtype=float)
            if len(values) != trial.num_scans or np.isnan(values).any():
                raise ValueError(
                    f"'{trial.profile_file_path}' has {len(values)} rows for {trial.num_scans} scans, or empty values."
                )
            profiles.append(values)
            start = (trial.start_time - first.start_time).total_seconds()
            timestamps.append(start + profile[f"{trial.profile_column} Time(ms)"].to_numpy(dtype=float) / 1000)
        self._timestamps = np.concatenate(timestamps)
        self._sampling_frequency = 1 / first.scan_line_period

        roi_id = "LineProfile"
        self._roi_ids = [roi_id]
        self._roi_responses.append(_RoiResponse("raw", np.concatenate(profiles)[:, np.newaxis], self._roi_ids))

        start_x, start_y, stop_x, stop_y = first.line
        x0, x1 = first.profile_pixels
        if x1 >= first.num_pixels:
            raise ValueError(f"Trial '{first.name}' profiles pixels {x0} to {x1} of a line of {first.num_pixels}.")
        line_xs = np.round(np.linspace(start_x, stop_x, first.num_pixels)).astype(int)
        line_ys = np.round(np.linspace(start_y, stop_y, first.num_pixels)).astype(int)
        image_mask = np.zeros((*first.frame_shape, 1), dtype=bool)
        image_mask[line_ys[x0 : x1 + 1], line_xs[x0 : x1 + 1], 0] = True
        self._roi_masks = _ROIMasks(
            data=image_mask,
            mask_tpe="nwb-image_mask",
            field_of_view_shape=first.frame_shape,
            roi_id_map={roi_id: 0},
        )

    def get_trials(self) -> list[LineScanTrial]:
        """Return the trials in the order they were recorded, which is the order of their samples."""
        return list(self._trials)

    def get_native_timestamps(self, start_sample: int | None = None, end_sample: int | None = None) -> np.ndarray:
        """Return the time of each scan, in seconds from the start of the first trial."""
        return self._timestamps[start_sample:end_sample]

    def _has_native_timestamps(self) -> bool:
        return True

    def get_frame_shape(self) -> tuple[int, int]:
        """Return ``(height, width)`` of the snapshot the scan line is drawn on."""
        return self._trials[0].frame_shape
