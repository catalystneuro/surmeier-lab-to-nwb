"""NeuroConv-portable interface for Bruker Prairie View BOT segmentation data.

Wraps :class:`roiextractors.BrukerBOTSegmentationExtractor` and emits metadata in the
new dict-based ophys shape (PR catalystneuro/neuroconv#1653). Lives in this repo for
now; structured so it can be lifted into ``neuroconv.datainterfaces.ophys.brukertiff``
with minimal changes.
"""

from pathlib import Path

from neuroconv.datainterfaces.ophys.basesegmentationextractorinterface import (
    BaseSegmentationExtractorInterface,
)
from neuroconv.utils import DeepDict
from pydantic import DirectoryPath


class BrukerBOTSegmentationInterface(BaseSegmentationExtractorInterface):
    """Conversion interface for Bruker Prairie View Brightness Over Time (BOT) data.

    BOT acquisitions monitor fluorescence in user-defined regions during continuous
    imaging. Prairie View writes a per-region brightness CSV computed live (one column
    per ``<Region>`` in the master XML). This interface wraps the corresponding
    roiextractors extractor and exposes the regions as ROIs in NWB.

    Parameters
    ----------
    folder_path : DirectoryPath
        Folder containing the Prairie View BOT recording: ``<folder>.xml`` master file
        and ``<folder>_Cycle{NNNNN}-botData.csv``.
    cycle : int, default: 1
        Cycle index to read.
    channel_name : str, optional
        Restrict the segmentation to one Prairie View channel (e.g. ``"Ch2"``). When
        ``None``, all regions across all channels are exposed.
    metadata_key : str, optional
        Metadata key. Defaults to ``"bruker_bot_segmentation"`` (with a
        ``_{channel_name}`` suffix when ``channel_name`` is provided).
    verbose : bool, default: False
    """

    display_name = "Bruker BOT Segmentation"
    associated_suffixes = (".xml", ".csv")
    info = "Interface for Bruker Prairie View Brightness Over Time (BOT) recordings."

    @classmethod
    def get_extractor_class(cls):
        from ..extractors.bruker_bot_segmentation_extractor import (
            BrukerBOTSegmentationExtractor,
        )

        return BrukerBOTSegmentationExtractor

    def __init__(
        self,
        folder_path: DirectoryPath,
        *,
        cycle: int = 1,
        channel_name: str | None = None,
        metadata_key: str | None = None,
        verbose: bool = False,
    ):
        if metadata_key is None:
            metadata_key = "bruker_bot_segmentation"
            if channel_name is not None:
                metadata_key = f"{metadata_key}_{channel_name}"

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
        ``use_new_metadata_format=True`` opt-in. Safe to delete when this is lifted
        upstream and the upstream interface's default is decided.
        """
        if metadata is None:
            metadata = self.get_metadata(use_new_metadata_format=True)
        return super().add_to_nwbfile(nwbfile=nwbfile, metadata=metadata, **kwargs)

    def get_metadata(self, *, use_new_metadata_format: bool = False) -> DeepDict:
        """Return metadata.

        When ``use_new_metadata_format`` is True (or when called via the default path
        below) populates ``Devices[metadata_key]``, ``Ophys.ImagingPlanes[metadata_key]``,
        ``Ophys.PlaneSegmentations[metadata_key]``, and ``Ophys.RoiResponses[metadata_key]``
        keyed by ``self.metadata_key``. Device defaults match Bruker Prairie View; imaging
        plane defaults are conservative placeholders since the BOT XML does not record
        excitation/emission wavelengths or indicator identity (set those at the call site).
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

        # New dict-based ophys shape.
        device_description = "Bruker Prairie View two-photon microscope (BOT acquisition)."
        metadata["Devices"] = {
            self.metadata_key: {
                "name": "BrukerFluorescenceMicroscope",
                "description": device_description,
            }
        }

        channel_label = self.channel_name or "AllChannels"
        imaging_plane_description = f"Imaging plane for Bruker BOT recording ({channel_label})."
        roi_responses_description = "BOT region brightness traces from Prairie View live computation."

        # Placeholders for the NWB-required imaging-plane fields the Bruker XML does not
        # declare (excitation/emission/indicator/location). Callers should override these
        # for their specific protocol (e.g. GRABACh3.0 = 920 nm excitation, 520 nm emission).
        import numpy as np

        metadata["Ophys"] = {
            "ImagingPlanes": {
                self.metadata_key: {
                    "name": f"ImagingPlane{channel_label}",
                    "description": imaging_plane_description,
                    "device_metadata_key": self.metadata_key,
                    "excitation_lambda": np.nan,
                    "indicator": "unknown",
                    "location": "unknown",
                    "optical_channel": [
                        {
                            "name": f"OpticalChannel{channel_label}",
                            "description": "An optical channel of the microscope.",
                            "emission_lambda": np.nan,
                        }
                    ],
                }
            },
            "PlaneSegmentations": {
                self.metadata_key: {
                    "name": f"PlaneSegmentation{channel_label}",
                    "description": f"BOT region segmentation ({channel_label}).",
                    "imaging_plane_metadata_key": self.metadata_key,
                }
            },
            "RoiResponses": {
                self.metadata_key: {
                    "raw": {
                        "name": f"RoiResponseSeries{channel_label}",
                        "description": roi_responses_description,
                        "unit": "a.u.",
                    },
                }
            },
        }

        return metadata
