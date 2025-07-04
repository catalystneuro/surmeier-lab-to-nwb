"""Optical physiology interfaces for converting Prairie View data."""

from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
import tifffile
import xmltodict
from neuroconv.basedatainterface import BaseDataInterface
from neuroconv.utils import DeepDict
from pynwb.base import TimeSeries
from pynwb.file import NWBFile
from pynwb.image import Image, Images
from pynwb.ophys import (
    Fluorescence,
    ImageSegmentation,
    OpticalChannel,
    RoiResponseSeries,
)


class PrairieViewLineScanInterface(BaseDataInterface):
    """Interface for Prairie View line scan data."""

    keywords = ("ophys", "imaging", "line scan", "prairie view")

    def __init__(
        self,
        xml_metadata_file_path: str | Path,
        channel_name: str,
        ophys_metadata_key: str = "PrairieView",
    ):
        """
        Initialize the interface with path to the XML metadata file.

        Parameters
        ----------
        xml_metadata_file_path : str | Path
            Path to the XML metadata file containing line scan parameters
        channel_name : str
            Channel name to process ("Ch1" for structural/Alexa568, "Ch2" for calcium/Fluo4)
        ophys_metadata_key : str, default="PrairieView"
            Key to use for organizing metadata in the NWB file
        """
        self.xml_metadata_file_path = Path(xml_metadata_file_path)
        self.channel_name = channel_name
        self.ophys_metadata_key = ophys_metadata_key

        parent_folder = self.xml_metadata_file_path.parent
        # Load XML metadata
        with open(self.xml_metadata_file_path, "r") as xml_file:
            xml_general_metadata_file = xml_file.read()
            self.xml_general_metadata_dict = xmltodict.parse(xml_general_metadata_file)

        # Extract line scan metadata
        prairie_view_metadata = self.xml_general_metadata_dict["PVScan"]
        sequence_metadata = prairie_view_metadata["Sequence"]
        line_scan_definition_metadata = sequence_metadata["PVLinescanDefinition"]
        profiles_metadata = line_scan_definition_metadata["LineScanProfiles"]
        profile_csv_data_file_name = profiles_metadata["@DataFile"]

        self.profile_csv_data_file_path = parent_folder / f"{profile_csv_data_file_name}"
        if not self.profile_csv_data_file_path.is_file():
            raise FileNotFoundError(f"Line scan profile data file not found: {self.profile_csv_data_file_path}")

        # Fetch source and raw line scan data file paths for the specific channel
        frame_metadata = sequence_metadata["Frame"]
        file_metadata = frame_metadata["File"]
        source_data_file_names_dict = {file["@channelName"]: file["@source"] for file in file_metadata}
        line_scan_data_raw_data_file_names_dict = {file["@channelName"]: file["@filename"] for file in file_metadata}

        if self.channel_name not in source_data_file_names_dict:
            raise ValueError(
                f"Channel {self.channel_name} not found in metadata. Available channels: {list(source_data_file_names_dict.keys())}"
            )

        # Get file paths for this specific channel only
        source_data_file_name = source_data_file_names_dict[self.channel_name]
        self.source_data_file_path = parent_folder / source_data_file_name
        if not self.source_data_file_path.is_file():
            raise FileNotFoundError(f"Source image file not found: {self.source_data_file_path}")

        raw_line_scan_data_file_name = line_scan_data_raw_data_file_names_dict[self.channel_name]
        self.line_scan_raw_data_file_path = parent_folder / raw_line_scan_data_file_name
        if not self.line_scan_raw_data_file_path.is_file():
            raise FileNotFoundError(f"Raw line scan data file not found: {self.line_scan_raw_data_file_path}")

    def get_metadata(self) -> DeepDict:
        """
        Extract metadata from the XML files.

        Returns
        -------
        DeepDict
            Dictionary containing metadata
        """
        metadata = super().get_metadata()

        # Add session start time if available
        session_start_time = self.get_session_start_time()
        if session_start_time is not None:
            metadata["NWBFile"]["session_start_time"] = session_start_time

        # Extract microscope metadata from XML
        prairie_view_metadata = self.xml_general_metadata_dict["PVScan"]
        prairie_view_state_metadata_list = prairie_view_metadata["PVStateShard"]["PVStateValue"]

        # Get objective lens info
        objective_lens = next(
            (item["@value"] for item in prairie_view_state_metadata_list if item["@key"] == "objectiveLens"),
            "Olympus 60X",
        )

        # Device metadata
        self.device_name = "BrukerUltima"
        metadata["Devices"][self.ophys_metadata_key] = {
            "name": self.device_name,
            "description": "Bruker two-photon microscope",
        }

        # Extract spatial resolution
        microns_per_pixel = next(
            (item for item in prairie_view_state_metadata_list if item["@key"] == "micronsPerPixel"), None
        )
        microns_per_pixel_dict = {value["@index"]: value["@value"] for value in microns_per_pixel["IndexedValue"]}
        grid_spacing = [float(microns_per_pixel_dict["XAxis"]), float(microns_per_pixel_dict["YAxis"])]

        # Ophys-specific metadata
        ophys_metadata = metadata["Ophys"]

        self.imaging_plane_name = "ImagingPlaneLineScan"
        ophys_metadata["ImagingPlanes"][self.ophys_metadata_key] = {
            "name": self.imaging_plane_name,
            "description": f"Line scan imaging using {objective_lens} objective",
            "device": self.ophys_metadata_key,
            "excitation_lambda": 810.0,
            "indicator": "Fluo-4",
            "location": "dorsolateral striatum",
            "grid_spacing": grid_spacing,
            "grid_spacing_unit": "microns",
        }

        self.plane_segmentation_name = "PlaneSegmentationLineScan"
        ophys_metadata["PlaneSegmentation"][self.ophys_metadata_key] = {
            "name": self.plane_segmentation_name,
            "description": "Line scan ROI segmentation for dendritic calcium imaging",
            "imaging_plane": self.ophys_metadata_key,
        }

        # Use generic channel-based naming - fluorophore info comes from metadata updates
        indicator_name = f"Channel{self.channel_name}"
        channel_type = "imaging"
        indicator_description = f"Fluorescence data from {self.channel_name}"
        timeseries_description = f"Line scan raw data for {self.channel_name}"
        image_description = f"Source image for {self.channel_name}"

        # Single ROI response series for this channel
        ophys_metadata["RoiResponseSeries"][self.ophys_metadata_key] = {
            "name": f"RoiResponseSeriesLineScan{indicator_name}",
            "description": f"{channel_type.capitalize()} fluorescence from {indicator_description}",
            "unit": "a.u.",
        }

        # Single TimeSeries for this channel (top level)
        metadata["TimeSeries"][self.ophys_metadata_key] = {
            "name": f"TimeSeriesLineScanRawData{indicator_name}",
            "description": timeseries_description,
            "unit": "a.u.",
        }

        # Single source image in acquisition metadata
        acquisition_metadata = metadata["Acquisition"]
        self.source_images_name = "ImagesLineScan"
        acquisition_metadata["SourceImages"][self.ophys_metadata_key] = {
            "name": f"ImageSource{indicator_name}",
            "description": image_description,
        }

        return metadata

    @staticmethod
    def get_session_start_time_from_file(xml_metadata_file_path: str | Path) -> datetime:
        """
        Extract session start time from XML metadata file without creating interface instance.

        Parameters
        ----------
        xml_metadata_file_path : str | Path
            Path to the XML metadata file

        Returns
        -------
        datetime
            Session start time from XML metadata
        """
        from datetime import datetime
        from pathlib import Path
        from zoneinfo import ZoneInfo

        import xmltodict

        xml_metadata_file_path = Path(xml_metadata_file_path)

        # Load XML metadata
        with open(xml_metadata_file_path, "r") as xml_file:
            xml_general_metadata_file = xml_file.read()
            xml_general_metadata_dict = xmltodict.parse(xml_general_metadata_file)

        # Extract session start time from PVScan element
        prairie_view_metadata = xml_general_metadata_dict["PVScan"]

        # Get the date attribute from PVScan element
        # Format is like "5/23/2016 1:58:21 PM"
        if "@date" in prairie_view_metadata:
            date_str = prairie_view_metadata["@date"]
        else:
            raise ValueError(f"Could not find date attribute in PVScan element")

        # Parse the date string in format "M/D/YYYY H:MM:SS AM/PM"
        # Illinois is in Central Time Zone
        central_tz = ZoneInfo("America/Chicago")

        # Parse using strptime - this handles various month/day formats automatically
        session_start_time = datetime.strptime(date_str, "%m/%d/%Y %I:%M:%S %p")

        # Add timezone info
        session_start_time = session_start_time.replace(tzinfo=central_tz)

        return session_start_time

    def get_session_start_time(self) -> datetime:
        """
        Extract session start time from the XML metadata.

        Returns
        -------
        datetime
            Session start time from XML metadata
        """
        from datetime import datetime
        from zoneinfo import ZoneInfo

        # Extract session start time from PVScan element
        prairie_view_metadata = self.xml_general_metadata_dict["PVScan"]

        # Get the date attribute from PVScan element
        # Format is like "5/23/2016 1:58:21 PM"
        if "@date" in prairie_view_metadata:
            date_str = prairie_view_metadata["@date"]
        else:
            raise ValueError(f"Could not find date attribute in PVScan element")

        # Parse the date string in format "M/D/YYYY H:MM:SS AM/PM"
        # Illinois is in Central Time Zone
        central_tz = ZoneInfo("America/Chicago")

        # Parse using strptime - this handles various month/day formats automatically
        session_start_time = datetime.strptime(date_str, "%m/%d/%Y %I:%M:%S %p")

        # Add timezone info
        session_start_time = session_start_time.replace(tzinfo=central_tz)

        return session_start_time

    @staticmethod
    def get_available_channels(xml_metadata_file_path: str | Path) -> list[str]:
        """
        Get available channels from XML metadata file without creating interface instance.

        Parameters
        ----------
        xml_metadata_file_path : str | Path
            Path to the XML metadata file

        Returns
        -------
        list[str]
            List of available channel names from the XML metadata
        """
        from pathlib import Path

        import xmltodict

        xml_metadata_file_path = Path(xml_metadata_file_path)

        # Load XML metadata
        with open(xml_metadata_file_path, "r") as xml_file:
            xml_general_metadata_file = xml_file.read()
            xml_general_metadata_dict = xmltodict.parse(xml_general_metadata_file)

        # Extract channel information
        prairie_view_metadata = xml_general_metadata_dict["PVScan"]
        sequence_metadata = prairie_view_metadata["Sequence"]
        frame_metadata = sequence_metadata["Frame"]
        file_metadata = frame_metadata["File"]

        # Get all channel names
        available_channels = [file["@channelName"] for file in file_metadata]

        return available_channels

    def add_to_nwbfile(self, nwbfile: NWBFile, metadata: dict | None = None):
        """
        Add the line scan data to the NWB file.

        Parameters
        ----------
        nwbfile : NWBFile
            NWB file to add the data to
        metadata : dict | None
            Metadata dictionary containing naming parameters like trial_id, position, cell_info
        """
        # Get organized metadata
        metadata = metadata or self.get_metadata()

        # Create or get device
        device_metadata = metadata["Devices"][self.ophys_metadata_key]
        device_name = device_metadata["name"]
        if device_name in nwbfile.devices:
            device = nwbfile.devices[device_name]
        else:
            device_description = device_metadata["description"]
            device = nwbfile.create_device(name=device_name, description=device_description)

        # Create or get imaging plane
        imaging_plane_metadata = metadata["Ophys"]["ImagingPlanes"][self.ophys_metadata_key]
        imaging_plane_name = imaging_plane_metadata["name"]
        available_imaging_planes = nwbfile.imaging_planes
        if imaging_plane_name in available_imaging_planes:
            imaging_plane = available_imaging_planes[imaging_plane_name]
        else:
            # Create optical channels
            optical_channel_fluo4 = OpticalChannel(
                name="GaAsPFluo4",
                description="GaAsP PMT for Fluo-4 calcium indicator",
                emission_lambda=525.0,  # Center of 490-560 nm range
            )

            optical_channel_alexa568 = OpticalChannel(
                name="PMTAlexa568",
                description="Side-on PMT for Alexa Fluor 568 structural dye",
                emission_lambda=600.0,  # Center of 580-620 nm range
            )

            optical_channels = [optical_channel_fluo4, optical_channel_alexa568]

            imaging_plane = nwbfile.create_imaging_plane(
                name=imaging_plane_name,
                description=imaging_plane_metadata["description"],
                optical_channel=optical_channels,
                excitation_lambda=imaging_plane_metadata["excitation_lambda"],
                indicator=imaging_plane_metadata["indicator"],
                location=imaging_plane_metadata["location"],
                grid_spacing=imaging_plane_metadata["grid_spacing"],
                grid_spacing_unit=imaging_plane_metadata["grid_spacing_unit"],
                device=device,
            )

        # Create ophys processing module if it doesn't exist
        if "ophys" not in nwbfile.processing:
            ophys_module = nwbfile.create_processing_module(
                name="ophys",
                description="optical physiology processed data",
            )
        else:
            ophys_module = nwbfile.processing["ophys"]

        # Create image segmentation (fixed name)
        image_segmentation_name = "ImageSegmentation"
        if image_segmentation_name in ophys_module.data_interfaces:
            image_segmentation = ophys_module.data_interfaces[image_segmentation_name]
        else:
            image_segmentation = ImageSegmentation(name=image_segmentation_name)
            ophys_module.add(image_segmentation)

        # Create plane segmentation
        plane_segmentation_metadata = metadata["Ophys"]["PlaneSegmentation"][self.ophys_metadata_key]
        plane_segmentation_name = plane_segmentation_metadata["name"]
        available_plane_segmentations = image_segmentation.plane_segmentations
        if plane_segmentation_name in available_plane_segmentations:
            plane_segmentation = available_plane_segmentations[plane_segmentation_name]
        else:
            plane_segmentation = image_segmentation.create_plane_segmentation(
                name=plane_segmentation_name,
                description=plane_segmentation_metadata["description"],
                imaging_plane=imaging_plane,
            )

        # Load line scan profile data for this specific channel
        line_scan_data_df = pd.read_csv(self.profile_csv_data_file_path)

        # Extract profile data for the specific channel
        # Map channel names to profile columns
        profile_mapping = {"Ch1": (" Prof 1", "Prof 1 Time(ms)"), "Ch2": (" Prof 2", " Prof 2 Time(ms)")}

        if self.channel_name not in profile_mapping:
            raise ValueError(
                f"Unsupported channel: {self.channel_name}. Supported channels: {list(profile_mapping.keys())}"
            )

        profile_col, time_col = profile_mapping[self.channel_name]
        profile_data = line_scan_data_df[profile_col].to_numpy()
        timestamps = line_scan_data_df[time_col].to_numpy(dtype=float) / 1000.0  # Convert to seconds

        # Extract line scan geometry from XML for ROI creation
        prairie_view_metadata = self.xml_general_metadata_dict["PVScan"]
        sequence_metadata = prairie_view_metadata["Sequence"]
        line_scan_definition_metadata = sequence_metadata["PVLinescanDefinition"]
        line_metadata = line_scan_definition_metadata["Line"]
        start_pixel_y = int(line_metadata["@startPixelY"])
        start_pixel_x = int(line_metadata["@startPixelX"])
        stop_pixel_x = int(line_metadata["@stopPixelX"])
        line_length = float(line_metadata["@lineLength"])

        # Create pixel mask (list of [x, y, weight] for each pixel in the ROI)
        pixel_mask = []
        for x in range(start_pixel_x, stop_pixel_x + 1):
            pixel_mask.append([x, start_pixel_y, 1.0])
        pixel_mask = np.array(pixel_mask)

        # Add ROI to plane segmentation
        plane_segmentation.add_roi(
            id=0,
            pixel_mask=pixel_mask,
        )

        # Create ROI table region with better description
        roi_table_region = plane_segmentation.create_roi_table_region(
            region=[0], description="Dendritic segment for calcium imaging"
        )

        # Create ROI response series using metadata structure
        roi_metadata = metadata["Ophys"]["RoiResponseSeries"][self.ophys_metadata_key]

        roi_response_series = RoiResponseSeries(
            name=roi_metadata["name"],
            description=roi_metadata["description"],
            data=profile_data,
            rois=roi_table_region,
            unit=roi_metadata["unit"],
            timestamps=timestamps,
        )

        # Get ophys module (created by helper method)
        ophys_module = nwbfile.processing["ophys"]

        # Create fluorescence container (fixed name)
        name = "Fluorescence"
        if name in ophys_module.data_interfaces:
            fluorescence = ophys_module.data_interfaces[name]
        else:
            fluorescence = Fluorescence(name=name)
            ophys_module.add(fluorescence)

        fluorescence.add_roi_response_series(roi_response_series)

        # Add fluorescence to ophys module

        # Add source image using metadata structure
        source_image_metadata = metadata["Acquisition"]["SourceImages"][self.ophys_metadata_key]

        # Skip Dodt channel if it's the current channel
        if self.channel_name == "Dodt":
            print(f"Skipping unsupported channel: {self.channel_name}")
            return

        # Create or get SourceImagesLineScan Images container
        source_images_container_name = "SourceImagesLineScan"
        if source_images_container_name in nwbfile.acquisition:
            source_images_container = nwbfile.acquisition[source_images_container_name]
        else:
            source_images_container = Images(
                name=source_images_container_name,
                description="Source images for line scan experiments, can be used to identify the recorded region",
            )
            nwbfile.add_acquisition(source_images_container)

        # Create and add source image to the container
        source_data = tifffile.imread(self.source_data_file_path)
        source_image = Image(
            name=source_image_metadata["name"], data=source_data, description=source_image_metadata["description"]
        )
        source_images_container.add_image(source_image)

        # Add raw line scan data
        frame_metadata = sequence_metadata["Frame"]
        configuration_metadata = frame_metadata["PVStateShard"]["PVStateValue"]
        scan_line_period = next((item for item in configuration_metadata if item["@key"] == "scanLinePeriod"), None)
        scan_line_period = float(scan_line_period["@value"]) if scan_line_period else None
        lines_per_frame = next((item for item in configuration_metadata if item["@key"] == "linesPerFrame"), None)
        lines_per_frame = int(lines_per_frame["@value"]) if lines_per_frame else None
        pixels_per_line = next((item for item in configuration_metadata if item["@key"] == "pixelsPerLine"), None)
        pixels_per_line = int(pixels_per_line["@value"]) if pixels_per_line else None

        # Calculate scanning parameters
        rate = 1.0 / scan_line_period if scan_line_period else None

        # Add raw line scan data using metadata structure
        timeseries_metadata = metadata["TimeSeries"][self.ophys_metadata_key]

        raw_line_scan_data = tifffile.imread(self.line_scan_raw_data_file_path)
        raw_line_scan_data_series = TimeSeries(
            name=timeseries_metadata["name"],
            description=timeseries_metadata["description"],
            data=raw_line_scan_data,
            unit=timeseries_metadata["unit"],
            rate=rate,
        )
        nwbfile.add_acquisition(raw_line_scan_data_series)
