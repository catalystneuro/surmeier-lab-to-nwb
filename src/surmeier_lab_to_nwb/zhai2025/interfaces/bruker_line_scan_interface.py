"""NeuroConv-portable interface for Bruker PrairieView line scans.

Wraps :class:`BrukerLineScanSegmentationExtractor` for the ROI and its trace, and writes the kymographs and the
snapshots of the field of view itself. Lives in this repo for now; structured so it can be lifted into
``neuroconv.datainterfaces.ophys.brukertiff`` with minimal changes.
"""

from functools import partial
from pathlib import Path

import numpy as np
import tifffile
from neuroconv._temporal_alignment import _TemporalAlignment
from neuroconv.datainterfaces.ophys._metadata_schema import _keyed_registry
from neuroconv.datainterfaces.ophys.basesegmentationextractorinterface import (
    BaseSegmentationExtractorInterface,
)
from neuroconv.tools.roiextractors.roiextractors import _add_segmentation_to_nwbfile
from neuroconv.utils import DeepDict
from pydantic import DirectoryPath, validate_call
from pynwb import NWBFile, TimeSeries
from pynwb.image import GrayscaleImage, Images

from ..extractors.bruker_line_scan_segmentation_extractor import (
    BrukerLineScanSegmentationExtractor,
    _read_master_xml,
    _read_scan_line,
)

_SOURCE_IMAGES_NAME = "ImageLineScanSource"


class BrukerLineScanInterface(BaseSegmentationExtractorInterface):
    """
    Interface for one channel of a Bruker PrairieView line scan, over the trials of one scan line.

    PrairieView saves each trial of a line scan in its own folder. The interface writes the trials it is given in
    the order they were recorded, each at its own start: the kymographs (one row per scan, one column per pixel
    along the line) as one ``TimeSeries`` in acquisition, the pixels PrairieView averages into its line profile as
    the one row of a ``PlaneSegmentation``, with the profile as its ``RoiResponseSeries``, and the snapshot of the
    field of view of each trial as an image in the ``ImageLineScanSource`` container in acquisition. The two
    series share one timestamps array, so the gaps between trials are kept.

    Each trial is its own key in ``temporal_alignment``, its folder name, so that it can be placed on its own.
    The trials must share their scan line; ``get_scan_line`` reads it so that folders can be grouped first.
    """

    display_name = "Bruker Line Scan"
    associated_suffixes = (".xml", ".csv", ".tif")
    info = "Interface for Bruker PrairieView line scans."

    @classmethod
    def get_extractor_class(cls):
        return BrukerLineScanSegmentationExtractor

    @classmethod
    def get_scan_line(cls, folder_path: DirectoryPath) -> tuple[int, int, int, int]:
        """Return the line a trial scanned: start x, start y, stop x and stop y in pixels of its snapshot."""
        return _read_scan_line(_read_master_xml(Path(folder_path)))

    @validate_call
    def __init__(
        self,
        folder_paths: list[DirectoryPath],
        *,
        channel_name: str = "Ch1",
        roi_id: str = "LineProfile",
        metadata_key: str | None = None,
        verbose: bool = False,
    ):
        """
        Parameters
        ----------
        folder_paths : list of DirectoryPath
            The trial folders, each holding a PrairieView line scan of the same line.
        channel_name : str, default: "Ch1"
            PrairieView channel, as the XML names it (``"Ch1"``, ``"Ch2"``).
        roi_id : str, default: "LineProfile"
            ID of the ROI, written as its ``roi_name``. Line scans whose ROIs share a ``PlaneSegmentation`` need
            different IDs.
        metadata_key : str, optional
            Key of this line scan's entries in the metadata. Defaults to ``"bruker_line_scan_{channel_name}"``.
        verbose : bool, default: False
        """
        super().__init__(
            folder_paths=folder_paths,
            channel_name=channel_name,
            roi_id=roi_id,
            metadata_key=metadata_key or f"bruker_line_scan_{channel_name}",
            verbose=verbose,
        )
        self.channel_name = channel_name

        # One time-bearing object per trial, instead of the base's one for the whole series
        self._alignment = _TemporalAlignment()
        start_sample = 0
        for trial in self.segmentation_extractor.get_trials():
            end_sample = start_sample + trial.num_scans
            self._alignment._register(
                key=trial.name,
                get_native_times=partial(
                    self.segmentation_extractor.get_native_timestamps, start_sample=start_sample, end_sample=end_sample
                ),
                get_number_of_samples=partial(int, trial.num_scans),
            )
            start_sample = end_sample

    def get_trial_sample_ranges(self) -> dict[str, tuple[int, int]]:
        """Return the first sample and the number of samples of each trial in the written series, by trial name."""
        ranges = {}
        start_sample = 0
        for trial in self.segmentation_extractor.get_trials():
            ranges[trial.name] = (start_sample, trial.num_scans)
            start_sample += trial.num_scans
        return ranges

    def get_metadata_schema(self) -> dict:
        metadata_schema = super().get_metadata_schema()
        ophys_schema = metadata_schema["properties"]["Ophys"]
        ophys_schema["definitions"]["LineScanKymographEntry"] = dict(
            type="object",
            additionalProperties=True,
            properties=dict(name=dict(type="string"), description=dict(type="string"), unit=dict(type="string")),
        )
        ophys_schema["definitions"]["LineScanSourceImageEntry"] = dict(
            type="object",
            additionalProperties=True,
            properties=dict(name=dict(type="string"), description=dict(type="string")),
        )
        ophys_schema["properties"]["LineScanKymographs"] = _keyed_registry(
            "#/properties/Ophys/definitions/LineScanKymographEntry"
        )
        # Keyed by metadata_key and then by trial name
        ophys_schema["properties"]["LineScanSourceImages"] = dict(
            type="object",
            additionalProperties=_keyed_registry("#/properties/Ophys/definitions/LineScanSourceImageEntry"),
        )
        return metadata_schema

    def get_metadata(self) -> DeepDict:
        metadata = super().get_metadata()
        trials = self.segmentation_extractor.get_trials()
        x0, x1 = trials[0].profile_pixels
        channel_name = self.channel_name
        metadata["NWBFile"]["session_start_time"] = trials[0].start_time
        metadata["Devices"] = {"bruker_device": {"name": "BrukerFluorescenceMicroscope"}}
        metadata["Ophys"] = {
            "ImagingPlanes": {
                self.metadata_key: {
                    "name": f"ImagingPlane{channel_name}",
                    "description": f"Line-scan imaging plane for {channel_name}.",
                    "device_metadata_key": "bruker_device",
                }
            },
            "PlaneSegmentations": {
                self.metadata_key: {
                    "name": f"PlaneSegmentationLineScan{channel_name}",
                    "description": (
                        f"Pixels {x0} to {x1} along the scan line, counted from its start at 0, which PrairieView "
                        "averages into its line profile."
                    ),
                    "imaging_plane_metadata_key": self.metadata_key,
                }
            },
            "RoiResponses": {
                self.metadata_key: {
                    "raw": {
                        "name": f"RoiResponseSeriesLineScan{channel_name}",
                        "description": (
                            f"Mean over pixels {x0} to {x1} along the scan line at each scan, from PrairieView's "
                            "LineProfileData file."
                        ),
                        "unit": "a.u.",
                    }
                }
            },
            "LineScanKymographs": {
                self.metadata_key: {
                    "name": f"LineScan{channel_name}",
                    "description": "Line scan: one row per scan of the line, one column per pixel along it.",
                    "unit": "a.u.",
                }
            },
            "LineScanSourceImages": {
                self.metadata_key: {
                    trial.name: {
                        "name": f"Image{channel_name}Trial{trial_number}",
                        "description": f"Snapshot of the field of view PrairieView saved with trial {trial.name}.",
                    }
                    for trial_number, trial in enumerate(trials, start=1)
                }
            },
        }
        return metadata

    def add_to_nwbfile(self, nwbfile: NWBFile, metadata: dict | None = None):
        metadata = metadata or self.get_metadata()
        trials = self.segmentation_extractor.get_trials()

        timestamps = np.concatenate([self.temporal_alignment[trial.name].get_times() for trial in trials])
        if np.any(np.diff(timestamps) <= 0):
            raise ValueError("Once aligned, the trials overlap or are out of the order they were recorded in.")
        _add_segmentation_to_nwbfile(
            segmentation_extractor=self.segmentation_extractor,
            nwbfile=nwbfile,
            metadata=metadata,
            metadata_key=self.metadata_key,
            plane_segmentation_name=None,
            timestamps=timestamps,
            starting_time=None,
        )

        # The kymograph shares the times of the line profile, as both have one sample per scan
        roi_response_name = metadata["Ophys"]["RoiResponses"][self.metadata_key]["raw"]["name"]
        roi_response_series = nwbfile.processing["ophys"]["Fluorescence"][roi_response_name]
        if roi_response_series.timestamps is None:
            times = dict(starting_time=roi_response_series.starting_time, rate=roi_response_series.rate)
        else:
            times = dict(timestamps=roi_response_series)
        kymographs = []
        for trial in trials:
            # The kymograph of the first channel is an OME-TIFF that can list the files of the other channels too
            kymograph = tifffile.imread(trial.kymograph_file_path, is_ome=False)
            if kymograph.shape != (trial.num_scans, trial.num_pixels):
                raise ValueError(
                    f"'{trial.kymograph_file_path}' has shape {kymograph.shape}, where the XML gives "
                    f"{trial.num_scans} scans of {trial.num_pixels} pixels."
                )
            kymographs.append(kymograph)
        kymograph_metadata = metadata["Ophys"]["LineScanKymographs"][self.metadata_key]
        nwbfile.add_acquisition(TimeSeries(**kymograph_metadata, data=np.concatenate(kymographs), **times))

        if _SOURCE_IMAGES_NAME not in nwbfile.acquisition:
            nwbfile.add_acquisition(
                Images(
                    name=_SOURCE_IMAGES_NAME,
                    description=(
                        "Snapshots of the field of view that PrairieView saves with each line scan, one per trial "
                        "and channel. The scan lines and their ROIs are given in pixels of these images."
                    ),
                )
            )
        source_images = nwbfile.acquisition[_SOURCE_IMAGES_NAME]
        images_metadata = metadata["Ophys"]["LineScanSourceImages"][self.metadata_key]
        for trial in trials:
            image = tifffile.imread(trial.source_image_file_path)
            image_metadata = dict(images_metadata[trial.name])
            if image.shape != trial.frame_shape:
                height, width = trial.frame_shape
                image_metadata["description"] += (
                    f" This snapshot is {image.shape[1]} by {image.shape[0]} pixels, while the scan line and its ROI "
                    f"are given on a grid of {width} by {height}."
                )
            # NWB images are indexed (x, y), as the ROI's image mask is, so the (row, column) snapshot is transposed
            source_images.add_image(GrayscaleImage(**image_metadata, data=image.T))
