"""Interface for the raw kymograph TIFF of a Bruker line-scan recording.

A line-scan acquisition writes a 2D OME-TIFF per channel
(``<folder>_Cycle{N}_Ch{C}_000001.ome.tif``): rows are line repetitions over time,
columns are pixels along the scan line. This interface reads that TIFF and writes
it as a ``TimeSeries`` in ``nwbfile.acquisition``. The processed 1D fluorescence
profile (averaged along the line per time point) is handled separately by
:class:`BrukerLineScanInterface`.

Kept local to this repo so it can be lifted upstream (likely as a Bruker-line-scan
imaging extractor backed by `MultiTIFFMultiPageExtractor`) when neuroconv defines
that pattern.
"""

from pathlib import Path
from typing import Optional

import numpy as np
import tifffile
from neuroconv.basedatainterface import BaseDataInterface
from neuroconv.utils import DeepDict
from pydantic import DirectoryPath
from pynwb.base import TimeSeries

from .prairie_view_utils import (
    get_pv_scalar,
    get_session_start_time,
    parse_prairie_view_xml,
)


class BrukerLineScanKymographInterface(BaseDataInterface):
    """Write the raw 2D kymograph of one channel of a Bruker line-scan recording.

    Parameters
    ----------
    folder_path : DirectoryPath
        Folder containing the Prairie View line-scan recording.
    cycle : int, default: 1
    channel_name : str, default: "Ch1"
        Prairie View channel name (e.g. ``"Ch1"`` or ``"Ch2"``).
    metadata_key : str, optional
        Defaults to ``"bruker_line_scan_kymograph_{channel_name}"``. Used as the key
        for the interface's TimeSeries metadata entry.
    """

    display_name = "Bruker Line Scan Raw Kymograph"
    associated_suffixes = (".ome.tif", ".xml")
    info = "Raw 2D kymograph TIFF from a Bruker line-scan recording."

    def __init__(
        self,
        folder_path: DirectoryPath,
        *,
        cycle: int = 1,
        channel_name: str = "Ch1",
        metadata_key: Optional[str] = None,
    ):
        super().__init__()
        self.folder_path = Path(folder_path)
        self.cycle = cycle
        self.channel_name = channel_name
        self.metadata_key = metadata_key or f"bruker_line_scan_kymograph_{channel_name}"
        self._t_start: float = 0.0

        xml_path = self.folder_path / f"{self.folder_path.name}.xml"
        if not xml_path.is_file():
            raise FileNotFoundError(f"Bruker XML configuration file not found at '{xml_path}'.")
        self._xml_metadata = parse_prairie_view_xml(xml_path)

        sequence_metadata = self._xml_metadata["PVScan"]["Sequence"]
        if isinstance(sequence_metadata, list):
            try:
                sequence_metadata = next(s for s in sequence_metadata if int(s.get("@cycle", 0)) == cycle)
            except StopIteration:
                raise ValueError(f"Cycle {cycle} not found in {xml_path}.") from None
        self._sequence_metadata = sequence_metadata

        frame_metadata = sequence_metadata["Frame"]
        if isinstance(frame_metadata, list):
            frame_metadata = frame_metadata[0]
        file_entries = frame_metadata["File"]
        if not isinstance(file_entries, list):
            file_entries = [file_entries]
        filenames = {entry["@channelName"]: entry["@filename"] for entry in file_entries}
        if channel_name not in filenames:
            raise ValueError(f"channel_name={channel_name!r} not in available channels {sorted(filenames)}.")
        self._raw_tiff_path = self.folder_path / filenames[channel_name]
        if not self._raw_tiff_path.is_file():
            raise FileNotFoundError(f"Raw line-scan TIFF not found at '{self._raw_tiff_path}'.")

        scan_line_period = get_pv_scalar(self._xml_metadata, "scanLinePeriod")
        self._rate = 1.0 / float(scan_line_period) if scan_line_period else None

    def get_metadata(self) -> DeepDict:
        metadata = super().get_metadata()
        metadata["NWBFile"]["session_start_time"] = get_session_start_time(self._xml_metadata)
        metadata["TimeSeries"][self.metadata_key] = {
            "name": f"TimeSeriesLineScanRawData{self.channel_name}",
            "description": (
                f"Raw line-scan kymograph for {self.channel_name}: rows are line "
                "repetitions over time, columns are pixels along the scan line."
            ),
            "unit": "a.u.",
        }
        return metadata

    def set_aligned_starting_time(self, aligned_starting_time: float) -> None:
        self._t_start = aligned_starting_time

    def add_to_nwbfile(self, nwbfile, metadata: Optional[dict] = None):
        if metadata is None:
            metadata = self.get_metadata()
        ts_metadata = metadata["TimeSeries"][self.metadata_key]

        raw = tifffile.imread(self._raw_tiff_path)
        # Some recordings package multiple channels in one 3D TIFF: (channels, lines, pixels)
        if raw.ndim == 3:
            channel_index_map = {"Ch1": 0, "Ch2": 1}
            if self.channel_name not in channel_index_map:
                raise ValueError(
                    f"Unsupported channel name {self.channel_name!r} for 3D kymograph; expected Ch1 or Ch2."
                )
            raw = raw[channel_index_map[self.channel_name]]
        elif raw.ndim != 2:
            raise ValueError(f"Unexpected kymograph TIFF shape {raw.shape}; expected 2D or 3D.")

        kwargs = {
            "name": ts_metadata["name"],
            "description": ts_metadata["description"],
            "data": raw,
            "unit": ts_metadata["unit"],
        }
        if self._rate is not None:
            kwargs["rate"] = float(self._rate)
            kwargs["starting_time"] = float(self._t_start)
        else:
            # Fall back to uniform spacing implied by line count if scanLinePeriod is missing.
            kwargs["timestamps"] = np.arange(raw.shape[0], dtype=float) + self._t_start

        nwbfile.add_acquisition(TimeSeries(**kwargs))
