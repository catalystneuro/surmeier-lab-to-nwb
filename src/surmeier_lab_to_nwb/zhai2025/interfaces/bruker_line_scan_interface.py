"""NeuroConv-portable interface for Bruker Prairie View line scan recordings.

Wraps :class:`roiextractors.BrukerLineScanSegmentationExtractor` and emits metadata
in the new dict-based ophys shape. Lives in this repo for now; structured so it can
be lifted into ``neuroconv.datainterfaces.ophys.brukertiff`` with minimal changes.

The raw kymograph TIFF and the line-overlay source image are out of scope for this
interface (use a separate imaging extractor for the kymograph).
"""

from pathlib import Path

from neuroconv.datainterfaces.ophys.basesegmentationextractorinterface import (
    BaseSegmentationExtractorInterface,
)
from neuroconv.utils import DeepDict
from pydantic import DirectoryPath


class BrukerLineScanInterface(BaseSegmentationExtractorInterface):
    """Conversion interface for one channel of a Bruker line-scan recording.

    Reads the LineProfileData CSV via the matching roiextractors extractor and
    exposes the line as a single ROI with its 1D fluorescence trace.

    Parameters
    ----------
    folder_path : DirectoryPath
        Folder containing the Prairie View line-scan recording.
    cycle : int, default: 1
    channel_name : str, default: "Ch1"
        Prairie View channel. Must be ``"Ch1"`` or ``"Ch2"``.
    metadata_key : str, optional
        Defaults to ``"bruker_line_scan_{channel_name}"``.
    verbose : bool, default: False
    """

    display_name = "Bruker Line Scan"
    associated_suffixes = (".xml", ".csv")
    info = "Interface for Bruker Prairie View Linescan recordings."

    @classmethod
    def get_extractor_class(cls):
        from ..extractors.bruker_line_scan_segmentation_extractor import (
            BrukerLineScanSegmentationExtractor,
        )

        return BrukerLineScanSegmentationExtractor

    def __init__(
        self,
        folder_path: DirectoryPath,
        *,
        cycle: int = 1,
        channel_name: str = "Ch1",
        metadata_key: str | None = None,
        verbose: bool = False,
    ):
        if metadata_key is None:
            metadata_key = f"bruker_line_scan_{channel_name}"

        self.folder_path = Path(folder_path)
        self.channel_name = channel_name

        super().__init__(
            folder_path=folder_path,
            cycle=cycle,
            channel_name=channel_name,
            verbose=verbose,
            metadata_key=metadata_key,
        )

    def add_to_nwbfile(self, nwbfile, metadata=None, **kwargs):
        """Add segmentation to NWB, defaulting to new-shape metadata when not provided.

        Local override so callers in this repo don't need to remember the
        ``use_new_metadata_format=True`` opt-in. Also injects the shared Bruker
        DeviceModel into the Device metadata so the resulting NWB device links to
        a DeviceModel rather than using the deprecated ``manufacturer`` field.
        Safe to delete when lifted upstream and the upstream default is decided.
        """
        from surmeier_lab_to_nwb.zhai2025.devices import (
            get_or_create_bruker_ultima_model,
        )

        if metadata is None:
            metadata = self.get_metadata(use_new_metadata_format=True)
        metadata["Devices"][self.metadata_key]["model"] = get_or_create_bruker_ultima_model(nwbfile)
        return super().add_to_nwbfile(nwbfile=nwbfile, metadata=metadata, **kwargs)

    def get_metadata(self, *, use_new_metadata_format: bool = False) -> DeepDict:
        """Return metadata.

        New-format keys populated: ``Devices[metadata_key]``,
        ``Ophys.ImagingPlanes[metadata_key]``, ``Ophys.PlaneSegmentations[metadata_key]``,
        ``Ophys.RoiResponses[metadata_key]``. Excitation/indicator/location are left as
        placeholders; set them at the call site (typically per-experiment).
        """
        metadata = (
            super().get_metadata()
            if not use_new_metadata_format
            else super().get_metadata(use_new_metadata_format=True)
        )

        extractor = self.segmentation_extractor
        try:
            session_start_time = extractor.get_session_start_time()
        except Exception:
            session_start_time = None
        if session_start_time is not None:
            metadata["NWBFile"]["session_start_time"] = session_start_time

        if not use_new_metadata_format:
            return metadata

        line_length = float(extractor.get_line_length_microns())

        metadata["Devices"] = {
            self.metadata_key: {
                "name": "BrukerFluorescenceMicroscope",
                "description": "Bruker Prairie View two-photon microscope (line-scan acquisition).",
            }
        }
        # Placeholders for NWB-required imaging plane fields not in the Bruker XML;
        # callers should override per experiment (e.g. Fluo-4 dendritic excitability:
        # excitation_lambda=810 nm, indicator="Fluo-4").
        import numpy as np

        metadata["Ophys"] = {
            "ImagingPlanes": {
                self.metadata_key: {
                    "name": f"ImagingPlaneLineScan{self.channel_name}",
                    "description": (
                        f"Line-scan imaging plane for {self.channel_name}; line length {line_length:.3f} um."
                    ),
                    "device_metadata_key": self.metadata_key,
                    "excitation_lambda": np.nan,
                    "indicator": "unknown",
                    "location": "Caudoputamen",
                    "optical_channel": [
                        {
                            "name": f"OpticalChannel{self.channel_name}",
                            "description": "An optical channel of the microscope.",
                            "emission_lambda": np.nan,
                        }
                    ],
                }
            },
            "PlaneSegmentations": {
                self.metadata_key: {
                    "name": f"PlaneSegmentationLineScan{self.channel_name}",
                    "description": f"Line-scan ROI segmentation ({self.channel_name}).",
                    "imaging_plane_metadata_key": self.metadata_key,
                }
            },
            "RoiResponses": {
                self.metadata_key: {
                    "raw": {
                        "name": f"RoiResponseSeriesLineScan{self.channel_name}",
                        "description": (
                            f"Line-scan fluorescence profile averaged along the scan line " f"({self.channel_name})."
                        ),
                        "unit": "a.u.",
                    },
                }
            },
        }

        return metadata
