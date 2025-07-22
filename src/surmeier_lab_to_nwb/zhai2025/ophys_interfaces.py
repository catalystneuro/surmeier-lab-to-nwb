"""Optical physiology interfaces for converting Prairie View data."""

from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
import tifffile
import xmltodict
from neuroconv.basedatainterface import BaseDataInterface
from neuroconv.utils import DeepDict, calculate_regular_series_rate
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
    """
    Interface for Prairie View line scan data from combined patch clamp and two-photon imaging experiments.

    This interface handles line scan recordings acquired during dendritic excitability experiments where
    brief current steps are delivered to neurons while simultaneously recording calcium transients using
    two-photon laser scanning microscopy.

    Session Structure
    -----------------
    Each recording session corresponds to a single experimental trial and contains:

    ```
    session_folder/
    ├── session_name.xml                                    # Master experiment file (this interface's input)
    ├── session_name_Cycle00001_Ch1_000001.ome.tif         # Fluo-4 calcium channel raw kymograph data
    ├── session_name-Cycle00001_Ch1Source.tif               # Field of view with scan line overlay (Ch1)
    ├── session_name_Cycle00001_Ch2_000001.ome.tif         # Alexa Fluor 568 structural channel raw data
    ├── session_name-Cycle00001_Ch2Source.tif               # Field of view with scan line overlay (Ch2)
    ├── session_name_Cycle00001_LineProfileData.csv         # Processed fluorescence profiles over time
    ├── session_name_Cycle00001_VoltageOutput_001.xml       # Electrical stimulus protocol definition
    ├── session_name_Cycle00001_VoltageRecording_001.csv    # Intracellular electrophysiology data
    ├── session_name_Cycle00001_VoltageRecording_001.xml    # Electrophysiology metadata
    ├── session_name.env                                    # Environment configuration
    └── References/                                         # Calibration and reference files
    ```

    XML Metadata File Contents
    --------------------------
    The master XML file (session_name.xml) contains comprehensive experimental metadata including:

    - **Timing Information**: Session start time with precise timestamps for synchronization
    - **Optical Configuration**: Two-photon laser parameters (810 nm excitation), objective lens details
    - **Channel Definitions**:
        - Ch1: Fluo-4 calcium indicator (GaAsP PMT, 490-560 nm detection)
        - Ch2: Alexa Fluor 568 structural dye (side-on PMT, 580-620 nm detection)
    - **Line Scan Geometry**: Start/stop coordinates, line length, pixel dimensions
    - **Acquisition Parameters**: Pixel dwell time (10 μs), lines per frame, scanning speed
    - **Spatial Calibration**: Microns per pixel conversion factors for accurate measurements
    - **File References**: Paths to raw TIFF data, processed CSV profiles, and source images

    The interface processes one channel at a time to enable modular handling of structural vs.
    calcium imaging data, allowing separate metadata organization and temporal alignment.

    Experimental Context
    -------------------
    These recordings are part of dendritic excitability experiments where:
    1. Neurons are patch clamped and filled with calcium-sensitive (Fluo-4) and structural (Alexa 568) dyes
    2. Brief current steps (typically 3x 2nA injections, 2ms each, at 50Hz) trigger action potentials
    3. Line scans capture calcium transients in dendrites as action potentials back-propagate
    4. Simultaneous structural imaging provides anatomical reference for ROI placement
    """

    keywords = ("ophys", "imaging", "line scan", "prairie view")

    def __init__(
        self,
        file_path: str | Path,
        channel_name: str,
        ophys_metadata_key: str = "PrairieView",
    ):
        """
        Initialize the interface with path to the master XML metadata file.

        Parameters
        ----------
        file_path : str | Path
            Path to the master XML metadata file (session_name.xml) containing complete line scan
            experimental metadata including timing, optical configuration, channel definitions,
            scan geometry, and file references. This is the main experiment file that coordinates
            all data streams for the recording session.
        channel_name : str
            Channel name to process. Must be one of:
            - "Ch1": Fluo-4 calcium indicator channel (GaAsP PMT, 490-560 nm)
            - "Ch2": Alexa Fluor 568 structural reference channel (PMT, 580-620 nm)
        ophys_metadata_key : str, default="PrairieView"
            Unique key to use for organizing this channel's metadata in the NWB file.
            Should be descriptive and include recording identifier for multi-channel experiments.
        """
        self.xml_metadata_file_path = Path(file_path)
        self.channel_name = channel_name
        self.ophys_metadata_key = ophys_metadata_key
        self._t_start = 0.0  # Default starting time

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
    def get_session_start_time_from_file(file_path: str | Path) -> datetime:
        """
        Extract session start time from master XML metadata file without creating interface instance.

        Parameters
        ----------
        file_path : str | Path
            Path to the master XML metadata file (session_name.xml)

        Returns
        -------
        datetime
            Session start time from XML metadata
        """
        from datetime import datetime
        from pathlib import Path
        from zoneinfo import ZoneInfo

        import xmltodict

        xml_metadata_file_path = Path(file_path)

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
    def get_available_channels(file_path: str | Path) -> list[str]:
        """
        Get available channels from master XML metadata file without creating interface instance.

        Parameters
        ----------
        file_path : str | Path
            Path to the master XML metadata file (session_name.xml)

        Returns
        -------
        list[str]
            List of available channel names from the XML metadata
        """
        from pathlib import Path

        import xmltodict

        xml_metadata_file_path = Path(file_path)

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

    def set_aligned_starting_time(self, aligned_starting_time: float) -> None:
        """
        Set the aligned starting time for all time series data.

        Parameters
        ----------
        aligned_starting_time : float
            Starting time in seconds relative to session start for temporal alignment
        """
        self._t_start = aligned_starting_time

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
        # Map channel names to profile columns with flexible matching
        profile_mapping = {"Ch1": (" Prof 1", "Prof 1 Time(ms)"), "Ch2": (" Prof 2", " Prof 2 Time(ms)")}

        if self.channel_name not in profile_mapping:
            raise ValueError(
                f"Unsupported channel: {self.channel_name}. Supported channels: {list(profile_mapping.keys())}"
            )

        profile_col_template, time_col_template = profile_mapping[self.channel_name]

        # Find actual column names in the dataframe with robust matching
        # This handles inconsistent whitespace in CSV column headers that can occur
        # across different experimental sessions or Prairie View software versions
        available_cols = line_scan_data_df.columns.tolist()

        # Look for profile column with flexible matching strategy:
        # 1. Try exact match first (fastest)
        # 2. Try stripped match (handles leading/trailing whitespace differences)
        # This is necessary because Prairie View CSV exports can have inconsistent
        # whitespace in column headers (e.g., " Prof 2" vs "Prof 2")
        profile_col = None
        for col in available_cols:
            if col == profile_col_template or col.strip() == profile_col_template.strip():
                profile_col = col
                break

        if profile_col is None:
            raise ValueError(
                f"Could not find profile column for {self.channel_name}. "
                f"Looking for '{profile_col_template}' (or stripped variant) "
                f"but available columns are: {available_cols}"
            )

        # Look for time column with the same flexible matching strategy
        # Time columns can also have inconsistent whitespace in headers
        time_col = None
        for col in available_cols:
            if col == time_col_template or col.strip() == time_col_template.strip():
                time_col = col
                break

        if time_col is None:
            raise ValueError(
                f"Could not find time column for {self.channel_name}. "
                f"Looking for '{time_col_template}' (or stripped variant) "
                f"but available columns are: {available_cols}"
            )
        profile_data = line_scan_data_df[profile_col].to_numpy()
        timestamps = line_scan_data_df[time_col].to_numpy(dtype=float) / 1000.0  # Convert to seconds

        # Apply aligned starting time offset
        timestamps = timestamps + self._t_start

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
            pixel_mask=pixel_mask,
        )

        # Create ROI table region with better description
        roi_table_region = plane_segmentation.create_roi_table_region(
            region=[0], description="Dendritic segment for calcium imaging"
        )

        # Create ROI response series using metadata structure
        roi_metadata = metadata["Ophys"]["RoiResponseSeries"][self.ophys_metadata_key]

        # Use neuroconv pattern: check if timestamps are regular, use rate+starting_time or timestamps
        rate = calculate_regular_series_rate(series=timestamps)
        recording_t_start = timestamps[0]

        if rate is not None:
            # Regular timestamps - use starting_time + rate
            roi_response_series = RoiResponseSeries(
                name=roi_metadata["name"],
                description=roi_metadata["description"],
                data=profile_data,
                rois=roi_table_region,
                unit=roi_metadata["unit"],
                starting_time=float(recording_t_start),
                rate=rate,
            )
        else:
            # Irregular timestamps - use explicit timestamps
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

        # Handle both 2D (timestamps, pixels) and 3D (channels, timestamps, pixels) data formats
        if len(raw_line_scan_data.shape) == 3:
            # 3D data: (channels, timestamps, pixels) - extract the appropriate channel
            num_channels = raw_line_scan_data.shape[0]
            if self.channel_name == "Ch1":
                raw_line_scan_data = raw_line_scan_data[0]  # First channel
            elif self.channel_name == "Ch2":
                if num_channels < 2:
                    raise ValueError(
                        f"Requested Ch2 but only {num_channels} channels available in data with shape {raw_line_scan_data.shape}"
                    )
                raw_line_scan_data = raw_line_scan_data[1]  # Second channel
            else:
                raise ValueError(
                    f"Unknown channel name '{self.channel_name}' for 3D data with shape {raw_line_scan_data.shape}"
                )
        elif len(raw_line_scan_data.shape) == 2:
            # 2D data: (timestamps, pixels) - use as is for single-channel files
            pass
        else:
            raise ValueError(f"Unexpected raw line scan data shape: {raw_line_scan_data.shape}. Expected 2D or 3D.")

        # Use aligned starting time if set, otherwise use rate-based timing
        if self._t_start != 0.0:
            raw_line_scan_data_series = TimeSeries(
                name=timeseries_metadata["name"],
                description=timeseries_metadata["description"],
                data=raw_line_scan_data,
                unit=timeseries_metadata["unit"],
                starting_time=self._t_start,
                rate=rate,
            )
        else:
            raw_line_scan_data_series = TimeSeries(
                name=timeseries_metadata["name"],
                description=timeseries_metadata["description"],
                data=raw_line_scan_data,
                unit=timeseries_metadata["unit"],
                rate=rate,
            )
        nwbfile.add_acquisition(raw_line_scan_data_series)
