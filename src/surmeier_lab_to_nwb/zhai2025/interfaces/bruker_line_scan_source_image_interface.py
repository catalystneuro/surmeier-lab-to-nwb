"""Interface for the field-of-view source image of a Bruker line-scan recording.

Each line-scan acquisition writes a per-channel snapshot of the imaging field with
the scan-line drawn on top (``<folder>-Cycle{N}_Ch{C}Source.tif``). This interface
adds those snapshots to a shared ``Images`` container under ``nwbfile.acquisition``
so the recorded region can be visually identified later.
"""

from pathlib import Path
from typing import Optional

import numpy as np
import tifffile
from neuroconv.basedatainterface import BaseDataInterface
from neuroconv.utils import DeepDict
from pydantic import DirectoryPath
from pynwb.image import Image, Images

from .prairie_view_utils import parse_prairie_view_xml

_DEFAULT_CONTAINER_NAME = "ImageLineScanSource"


class BrukerLineScanSourceImageInterface(BaseDataInterface):
    """Write the line-scan field-of-view snapshot of one channel into NWB.

    Parameters
    ----------
    folder_path : DirectoryPath
        Folder containing the Prairie View line-scan recording.
    cycle : int, default: 1
    channel_name : str, default: "Ch1"
    metadata_key : str, optional
        Defaults to ``"bruker_line_scan_source_image_{channel_name}"``. Used as the
        key for the interface's image metadata entry.
    container_name : str, default: ``"ImageLineScanSource"``
        Name of the shared ``Images`` container the snapshot is appended to.
    """

    display_name = "Bruker Line Scan Source Image"
    associated_suffixes = (".tif", ".xml")
    info = "Field-of-view snapshot of a Bruker line-scan recording, with scan-line overlay."

    def __init__(
        self,
        folder_path: DirectoryPath,
        *,
        cycle: int = 1,
        channel_name: str = "Ch1",
        metadata_key: Optional[str] = None,
        container_name: str = _DEFAULT_CONTAINER_NAME,
    ):
        super().__init__()
        self.folder_path = Path(folder_path)
        self.cycle = cycle
        self.channel_name = channel_name
        self.metadata_key = metadata_key or f"bruker_line_scan_source_image_{channel_name}"
        self.container_name = container_name

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

        frame_metadata = sequence_metadata["Frame"]
        if isinstance(frame_metadata, list):
            frame_metadata = frame_metadata[0]
        file_entries = frame_metadata["File"]
        if not isinstance(file_entries, list):
            file_entries = [file_entries]
        source_filenames = {entry["@channelName"]: entry["@source"] for entry in file_entries if "@source" in entry}
        if channel_name not in source_filenames:
            raise ValueError(
                f"channel_name={channel_name!r} has no @source TIFF declared in XML. "
                f"Available with source: {sorted(source_filenames)}"
            )
        self._source_image_path = self.folder_path / source_filenames[channel_name]
        if not self._source_image_path.is_file():
            raise FileNotFoundError(f"Source image TIFF not found at '{self._source_image_path}'.")

    def get_metadata(self) -> DeepDict:
        metadata = super().get_metadata()
        metadata["Acquisition"]["SourceImages"][self.metadata_key] = {
            "name": f"ImageSourceChannel{self.channel_name}",
            "description": f"Source image for {self.channel_name}: field of view with scan-line overlay.",
        }
        return metadata

    def add_to_nwbfile(self, nwbfile, metadata: Optional[dict] = None):
        if metadata is None:
            metadata = self.get_metadata()
        image_metadata = metadata["Acquisition"]["SourceImages"][self.metadata_key]

        if self.container_name in nwbfile.acquisition:
            container = nwbfile.acquisition[self.container_name]
        else:
            container = Images(
                name=self.container_name,
                description=(
                    "Source images for line scan experiments: field of view snapshots with "
                    "scan-line overlay, used to identify the recorded region."
                ),
            )
            nwbfile.add_acquisition(container)

        data = tifffile.imread(self._source_image_path)
        # Convert RGB(A) snapshots to grayscale to match the bespoke implementation.
        if data.ndim == 3:
            if data.shape[2] == 3:
                data = np.mean(data, axis=2).astype(np.float32)
            elif data.shape[2] == 4:
                data = np.mean(data[:, :, :3], axis=2).astype(np.float32)

        image = Image(
            name=image_metadata["name"],
            data=data,
            description=image_metadata["description"],
        )
        container.add_image(image)
