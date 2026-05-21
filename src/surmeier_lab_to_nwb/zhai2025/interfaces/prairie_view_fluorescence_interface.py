from datetime import datetime
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
from neuroconv.basedatainterface import BaseDataInterface
from neuroconv.utils import DeepDict, calculate_regular_series_rate
from pynwb.device import Device
from pynwb.ophys import (
    Fluorescence,
    ImageSegmentation,
    OpticalChannel,
    RoiResponseSeries,
)

from .prairie_view_utils import get_session_start_time, parse_prairie_view_xml

# Channel-name to user-facing display mapping used in NWB object names.
# Surmeier-specific: Ch2 carries the GRABACh3.0 sensor, Dodt is the DIC reference channel.
_CHANNEL_DISPLAY_MAP = {"Ch2": "GRABCh", "Dodt": "DoDT"}


def _channel_display(channel_name: str) -> str:
    return _CHANNEL_DISPLAY_MAP.get(channel_name, channel_name)


class PrairieViewFluorescenceInterface(BaseDataInterface):
    """Interface for Prairie View Brightness Over Time (BOT) acetylcholine fluorescence data.

    Reads the per-region BOT CSV produced by Prairie View together with the master
    `.xml` describing region geometry and channel assignment, and writes
    `RoiResponseSeries` traces under one `Fluorescence` container in the NWB file.

    Emits metadata under the dict-based ophys shape (PR catalystneuro/neuroconv#1653):
    `Devices`, `Ophys.ImagingPlanes`, `Ophys.PlaneSegmentations`, `Ophys.RoiResponses`,
    cross-referenced by `device_metadata_key` / `imaging_plane_metadata_key`. Per-channel
    keys are formed as `f"{metadata_key}_{channel_name}"`.

    See `obsidian_docs/source_formats/bruker/bot_format.md` for the XML/CSV structure.
    """

    def __init__(
        self,
        bot_csv_data_file_path,
        xml_metadata_file_path,
        image_dimensions=None,
        *,
        metadata_key: str = "bot_acetylcholine",
    ):
        """Initialize the interface.

        Parameters
        ----------
        bot_csv_data_file_path : str or Path
            Path to the BOT CSV file with brightness measurements.
        xml_metadata_file_path : str or Path
            Path to the Prairie View master XML with `<PVBOTs>` region metadata.
        image_dimensions : tuple, optional
            Image dimensions (height, width) for creating ROI masks. Defaults to (512, 512).
        metadata_key : str, default: "bot_acetylcholine"
            Top-level metadata key shared across all channels of this acquisition. Per-channel
            keys are derived as `f"{metadata_key}_{channel_name}"`.
        """
        super().__init__()

        self.bot_csv_data_file_path = Path(bot_csv_data_file_path)
        self.xml_metadata_file_path = Path(xml_metadata_file_path)
        self.image_dimensions = image_dimensions or (512, 512)
        self.metadata_key = metadata_key
        self._t_start: float = 0.0

        self._load_roi_region_data()

    def get_session_start_time(self) -> datetime:
        """Return the session start time (Central Time) from the parsed XML."""
        return get_session_start_time(self.xml_metadata)

    def set_aligned_starting_time(self, aligned_starting_time: float) -> None:
        """Set the offset (seconds, relative to session start) to apply when adding traces."""
        self._t_start = aligned_starting_time

    def _channel_metadata_key(self, channel_name: str) -> str:
        return f"{self.metadata_key}_{channel_name}"

    def _load_roi_region_data(self):
        """Parse `<PVBOTs>` region geometry and the channel index → name mapping from XML."""
        self.xml_metadata = parse_prairie_view_xml(self.xml_metadata_file_path)

        sequence_metadata = self.xml_metadata["PVScan"]["Sequence"]
        first_frame = sequence_metadata["Frame"][0]
        file_metadata = first_frame["File"]
        channel_index_to_channel_name = {m["@channel"]: m["@channelName"] for m in file_metadata}
        region_metadata_list = sequence_metadata["PVBOTs"]["Region"]

        self.regions_info = {}
        for index, region in enumerate(region_metadata_list):
            channel_name = channel_index_to_channel_name[region["@channel"]]
            self.regions_info[channel_name] = {
                "x": int(region["@x"]),
                "y": int(region["@y"]),
                "width": int(region["@width"]),
                "height": int(region["@height"]),
                "region": f"Region {index + 1}",
            }

    def get_metadata(self) -> DeepDict:
        """Return metadata in the dict-based ophys format.

        Populates `NWBFile.session_start_time`, `Devices[metadata_key]`, and per-channel
        entries under `Ophys.ImagingPlanes`, `Ophys.PlaneSegmentations`, and
        `Ophys.RoiResponses`. Default values reflect the Surmeier GRABACh3.0 protocol
        (920 nm excitation, 520 nm emission, 21.26 Hz, "striatum"); override by mutating
        the returned dict before calling :meth:`add_to_nwbfile`.
        """
        metadata = super().get_metadata()
        metadata["NWBFile"]["session_start_time"] = self.get_session_start_time()

        metadata["Devices"][self.metadata_key] = {
            "name": "BrukerFluorescenceMicroscope",
            "description": "Bruker two-photon microscope for acetylcholine GRAB biosensor imaging",
        }

        imaging_planes = metadata["Ophys"]["ImagingPlanes"]
        plane_segmentations = metadata["Ophys"]["PlaneSegmentations"]
        roi_responses = metadata["Ophys"]["RoiResponses"]

        for channel_name in self.regions_info:
            channel_key = self._channel_metadata_key(channel_name)
            channel_display = _channel_display(channel_name)

            imaging_planes[channel_key] = {
                "name": f"ImagingPlane{channel_name}",
                "description": f"Imaging plane for acetylcholine fluorescence measurements in {channel_name}",
                "device_metadata_key": self.metadata_key,
                "excitation_lambda": 920.0,
                "imaging_rate": 21.26,
                "indicator": "GRABACh3.0",
                "location": "striatum",
                "optical_channel": [
                    {
                        "name": f"OpticalChannel{channel_name}",
                        "description": f"Optical channel for {channel_name} fluorescence detection",
                        "emission_lambda": 520.0,
                    }
                ],
            }

            plane_segmentations[channel_key] = {
                "name": f"PlaneSegmentation{channel_display}",
                "description": (f"Segmentation of regions for brightness over time analysis in {channel_display}."),
                "imaging_plane_metadata_key": channel_key,
            }

            roi_responses[channel_key] = {
                "raw": {
                    "name": f"RoiResponseSeries{channel_display}",
                    "description": (
                        f"Acetylcholine fluorescence brightness over time measurements "
                        f"for {self.regions_info[channel_name]['region']}"
                    ),
                    "unit": "a.u.",
                },
            }

        return metadata

    def add_to_nwbfile(
        self,
        nwbfile,
        metadata: Optional[dict] = None,
        fluorescence_series_name: Optional[str] = None,
        repetition_number: Optional[int] = None,
    ):
        """Add brightness over time data to the NWB file.

        Reads from `metadata` (defaulting to :meth:`get_metadata`) to determine device,
        imaging plane, segmentation, and series names. Devices and imaging planes are
        shared across calls (created on first use, reused thereafter). Plane segmentations
        are created once per channel; ROI response series are created per trial.

        Parameters
        ----------
        nwbfile : NWBFile
        metadata : dict, optional
            Dict-based metadata; uses :meth:`get_metadata` when None.
        fluorescence_series_name : str, optional
            Trial-specific suffix appended to the RoiResponseSeries name from metadata.
        repetition_number : int, optional
            Repetition number appended to the series name as `RepetitionNNN`.
        """
        if metadata is None:
            metadata = self.get_metadata()

        bot_df = pd.read_csv(self.bot_csv_data_file_path)
        if "Timestamp" not in bot_df.columns:
            raise ValueError(
                f"'Timestamp' column not found in {self.bot_csv_data_file_path}. "
                "Ensure the CSV file contains a 'Timestamp' column."
            )
        timestamps = bot_df["Timestamp"].to_numpy()

        if "ophys" not in nwbfile.processing:
            ophys_module = nwbfile.create_processing_module(
                name="ophys", description="optical physiology processed data"
            )
        else:
            ophys_module = nwbfile.processing["ophys"]

        if "ImageSegmentation" in ophys_module.data_interfaces:
            image_segmentation = ophys_module["ImageSegmentation"]
        else:
            image_segmentation = ImageSegmentation()
            ophys_module.add(image_segmentation)

        if "Fluorescence" in ophys_module.data_interfaces:
            fluorescence = ophys_module["Fluorescence"]
        else:
            fluorescence = Fluorescence(name="Fluorescence")
            ophys_module.add(fluorescence)

        device = self._get_or_create_device(nwbfile, metadata)

        for channel_name, region_metadata in self.regions_info.items():
            channel_key = self._channel_metadata_key(channel_name)
            imaging_plane = self._get_or_create_imaging_plane(nwbfile, metadata, channel_key, device)
            plane_segmentation = self._get_or_create_plane_segmentation(
                metadata, channel_key, image_segmentation, imaging_plane, region_metadata
            )
            self._add_roi_response_series(
                metadata,
                channel_key,
                plane_segmentation,
                fluorescence,
                bot_df,
                timestamps,
                region_metadata,
                fluorescence_series_name,
                repetition_number,
            )

    def _get_or_create_device(self, nwbfile, metadata):
        device_entry = metadata["Devices"][self.metadata_key]
        device_name = device_entry["name"]
        if device_name in nwbfile.devices:
            return nwbfile.devices[device_name]
        device = Device(name=device_name, description=device_entry.get("description"))
        nwbfile.add_device(device)
        return device

    def _get_or_create_imaging_plane(self, nwbfile, metadata, channel_key, device):
        plane_entry = metadata["Ophys"]["ImagingPlanes"][channel_key]
        plane_name = plane_entry["name"]
        if plane_name in nwbfile.imaging_planes:
            return nwbfile.imaging_planes[plane_name]

        optical_channels = [OpticalChannel(**oc) for oc in plane_entry["optical_channel"]]
        return nwbfile.create_imaging_plane(
            name=plane_name,
            description=plane_entry["description"],
            optical_channel=optical_channels,
            device=device,
            excitation_lambda=plane_entry["excitation_lambda"],
            imaging_rate=plane_entry["imaging_rate"],
            indicator=plane_entry["indicator"],
            location=plane_entry["location"],
        )

    def _get_or_create_plane_segmentation(
        self, metadata, channel_key, image_segmentation, imaging_plane, region_metadata
    ):
        seg_entry = metadata["Ophys"]["PlaneSegmentations"][channel_key]
        seg_name = seg_entry["name"]
        if seg_name in image_segmentation.plane_segmentations:
            return image_segmentation.plane_segmentations[seg_name]

        plane_segmentation = image_segmentation.create_plane_segmentation(
            name=seg_name,
            description=seg_entry["description"],
            imaging_plane=imaging_plane,
        )
        image_mask = np.zeros(shape=self.image_dimensions, dtype=bool)
        x, y = region_metadata["x"], region_metadata["y"]
        w, h = region_metadata["width"], region_metadata["height"]
        image_mask[y : y + h, x : x + w] = True
        plane_segmentation.add_roi(image_mask=image_mask)
        return plane_segmentation

    def _add_roi_response_series(
        self,
        metadata,
        channel_key,
        plane_segmentation,
        fluorescence,
        bot_df,
        timestamps,
        region_metadata,
        fluorescence_series_name,
        repetition_number,
    ):
        roi_entry = metadata["Ophys"]["RoiResponses"][channel_key]["raw"]
        base_name = roi_entry["name"]
        if fluorescence_series_name:
            series_name = f"{base_name}{fluorescence_series_name}"
        else:
            series_name = base_name
        if repetition_number is not None:
            series_name += f"Repetition{repetition_number:03d}"

        roi_table_region = plane_segmentation.create_roi_table_region(
            region=[0],
            description=roi_entry.get("description", "ROI"),
        )

        region_data = bot_df[region_metadata["region"]].to_numpy()
        rate = calculate_regular_series_rate(series=timestamps)
        common_kwargs = dict(
            name=series_name,
            description=roi_entry["description"],
            data=region_data,
            rois=roi_table_region,
            unit=roi_entry["unit"],
        )
        if rate is not None:
            roi_response_series = RoiResponseSeries(
                **common_kwargs,
                starting_time=float(self._t_start + timestamps[0]),
                rate=rate,
            )
        else:
            roi_response_series = RoiResponseSeries(
                **common_kwargs,
                timestamps=timestamps + self._t_start,
            )

        fluorescence.add_roi_response_series(roi_response_series)
