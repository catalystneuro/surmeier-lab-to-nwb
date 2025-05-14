"""Primary classes for converting Prairie View data."""

from datetime import datetime
from pathlib import Path
from typing import Optional, Union

import numpy as np
import pandas as pd
import tifffile
import xmltodict
from neuroconv.basedatainterface import BaseDataInterface
from neuroconv.utils import DeepDict
from pynwb.file import NWBFile
from pynwb.icephys import CurrentClampSeries
from pynwb.image import Image, Images
from pynwb.ophys import (
    Fluorescence,
    ImageSegmentation,
    OpticalChannel,
    RoiResponseSeries,
    TwoPhotonSeries,
)


class PrairieViewLineScanInterface(BaseDataInterface):
    """Interface for Prairie View line scan data."""

    keywords = ("ophys", "imaging", "line scan", "prairie view")

    def __init__(self, xml_metadata_file_path: Union[str, Path]):
        """
        Initialize the interface with path to the XML metadata file.

        Parameters
        ----------
        xml_metadata_file_path : Union[str, Path]
            Path to the XML metadata file containing line scan parameters
        """
        self.xml_metadata_file_path = Path(xml_metadata_file_path)
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

        # Fetch source and raw line scan data file paths
        file_metadata = frame_metadata["File"]
        source_data_file_names_dict = {file["@channelName"]: file["@source"] for file in file_metadata}
        line_scan_data_raw_data_file_names_dict = {file["@channelName"]: file["@source"] for file in file_metadata}
        self.source_data_file_paths_dict = {}
        self.line_scan_raw_data_file_paths_dict = {}
        for channel_name in source_data_file_names_dict:
            source_data_file_name = source_data_file_names_dict[channel_name]
            source_data_file_path = parent_folder / source_data_file_name
            if not source_data_file_path.is_file():
                raise FileNotFoundError(f"Source image file not found: {source_data_file_path}")
            self.source_data_file_paths_dict[channel_name] = source_data_file_path

            raw_line_scan_data_file_name = line_scan_data_raw_data_file_names_dict[channel_name]
            raw_line_scan_data_file_path = parent_folder / raw_line_scan_data_file_name
            if not raw_line_scan_data_file_path.is_file():
                raise FileNotFoundError(f"Raw line scan data file not found: {raw_line_scan_data_file_path}")
            self.line_scan_raw_data_file_paths_dict

    def get_metadata(self) -> DeepDict:
        """
        Extract metadata from the XML files.

        Returns
        -------
        DeepDict
            Dictionary containing metadata
        """
        metadata = super().get_metadata()

        return metadata

    def add_to_nwbfile(self, nwbfile: NWBFile, metadata: dict | None = None):
        """
        Add the line scan data to the NWB file.

        Parameters
        ----------
        nwbfile : NWBFile
            NWB file to add the data to
        metadata : dict
            Metadata dictionary
        """

        # Create two-photon imaging device
        if "Bruker Ultima" in nwbfile.devices:
            two_photon_device = nwbfile.devices["Bruker Ultima"]
        else:
            two_photon_device = nwbfile.create_device(
                name="Bruker Ultima", description="Bruker two-photon microscope", manufacturer="Bruker"
            )

        # Create optical channel
        optical_channel = OpticalChannel(name="GaAsP", description="GaAsP PMT", emission_lambda=525.0)

        # Create imaging plane
        imaging_plane = nwbfile.create_imaging_plane(
            name="dSPN_dendrite_plane",
            description="",
            optical_channel=optical_channel,
            excitation_lambda=920.0,
            indicator="Fluo-4",
            location="dorsolateral striatum",
            device=two_photon_device,
        )

        # Load line scan profile data
        line_scan_data_df = pd.read_csv(self.profile_csv_data_file_path)

        # Extract line scan profiles
        line_scan_profile_channel1 = line_scan_data_df[" Prof 1"].to_numpy()
        line_scan_profile_channel2 = line_scan_data_df[" Prof 2"].to_numpy()
        timestamps = line_scan_data_df["Prof 1 Time(ms)"].to_numpy(dtype=float) / 1000.0  # Convert to seconds

        # Create ophys processing module if it doesn't exist
        if "ophys" not in nwbfile.processing:
            ophys_module = nwbfile.create_processing_module(
                name="ophys", description="optical physiology processed data"
            )
        else:
            ophys_module = nwbfile.processing["ophys"]

        # Create image segmentation
        image_segmentation = ImageSegmentation()
        ophys_module.add(image_segmentation)

        # Create plane segmentation
        plane_segmentation = image_segmentation.create_plane_segmentation(
            name="PlaneSegmentation",
            description="output from segmenting my favorite imaging plane",
            imaging_plane=imaging_plane,
        )

        # Create pixel mask (list of [x, y, weight] for each pixel in the ROI)
        # Extract line coordinates
        # Extract line scan metadata
        prairie_view_metadata = self.xml_general_metadata_dict["PVScan"]
        sequence_metadata = prairie_view_metadata["Sequence"]
        line_scan_definition_metadata = sequence_metadata["PVLinescanDefinition"]

        # Extract source image dimensions
        source_pixels_per_line = int(line_scan_definition_metadata["@SourcePixelsPerLine"])
        source_lines_per_frame = int(line_scan_definition_metadata["@SourceLinesPerFrame"])
        line_metadata = line_scan_definition_metadata["Line"]
        start_pixel_y = int(line_metadata["@startPixelY"])
        start_pixel_x = int(line_metadata["@startPixelX"])
        stop_pixel_x = int(line_metadata["@stopPixelX"])
        line_length = float(line_metadata["@lineLength"])

        pixel_mask = []
        for x in range(start_pixel_x, stop_pixel_x + 1):
            pixel_mask.append([x, start_pixel_y, 1.0])
        pixel_mask = np.array(pixel_mask)

        # Add ROI to plane segmentation
        plane_segmentation.add_roi(
            id=0,
            pixel_mask=pixel_mask,
        )

        # Create ROI table region
        roi_table_region = plane_segmentation.create_roi_table_region(region=[0], description="the first of two ROIs")

        # Create ROI response series
        roi_resp_series1 = RoiResponseSeries(
            name="LineScanChannel1RoiResponseSeries",
            description="Fluorescence responses for ROI",
            data=line_scan_profile_channel1,
            rois=roi_table_region,
            unit="n.a",
            timestamps=timestamps,
        )

        roi_resp_series2 = RoiResponseSeries(
            name="LineScanChannel2RoiResponseSeries",
            description="Fluorescence responses for ROI",
            data=line_scan_profile_channel2,
            rois=roi_table_region,
            unit="n.a",
            timestamps=timestamps,
        )

        # Create fluorescence container
        fluorescence = Fluorescence()
        fluorescence.add_roi_response_series(roi_resp_series1)
        fluorescence.add_roi_response_series(roi_resp_series2)

        # Add fluorescence to ophys module
        ophys_module.add(fluorescence)

        # Load source image files
        frame_metadata = sequence_metadata["Frame"]
        file_metadata = frame_metadata["File"]

        # Load source images
        source_images = []
        for channel_name, source_data_file_path in self.source_data_file_paths_dict.items():
            source_data = tifffile.imread(source_data_file_path)
            source_image = Image(name=f"LineScanChannel{channel_name}SourceImage", data=source_data)
            source_images.append(source_image)

        # Create image container and add to ophys module
        image_container = Images(name="SourceImages", images=source_images)
        ophys_module.add(image_container)

        frame_metadata = sequence_metadata["Frame"]
        configuration_metadata = frame_metadata["PVStateShard"]["PVStateValue"]
        configuration_metadata
        scan_line_period = next((item for item in configuration_metadata if item["@key"] == "scanLinePeriod"), None)
        scan_line_period = float(scan_line_period["@value"]) if scan_line_period else None

        rate = 1.0 / scan_line_period

        # Add raw line scan data
        for channel_name, raw_line_scan_data_file_path in self.line_scan_raw_data_file_paths_dict.items():
            raw_line_scan_data = tifffile.imread(raw_line_scan_data_file_path)
            raw_line_scan_data_series = TwoPhotonSeries(
                name=f"LineScanRawDataChannel{channel_name}TimeSeries",
                data=raw_line_scan_data,
                imaging_plane=imaging_plane,
                unit="a.u.",
                rate=rate,
                starting_frame=0,
                format="raw",
                dimension=[source_pixels_per_line, source_lines_per_frame],
                field_of_view=[line_length, source_lines_per_frame],
                timestamps=timestamps,
            )
            nwbfile.add_acquisition(raw_line_scan_data_series)


class PraireViewIntracellularRecordingInterface(BaseDataInterface):
    """Interface for Prairie View intracellular recording data."""

    keywords = ["intracellular", "recording", "prairie view"]

    def __init__(self, xml_file_path: Union[str, Path], recording_csv_file_path: Optional[Union[str, Path]] = None):
        """
        Initialize the interface with paths to the XML and CSV files.

        Parameters
        ----------
        xml_file_path : Union[str, Path]
            Path to the XML file containing recording metadata
        recording_csv_file_path : Optional[Union[str, Path]], default=None
            Path to the CSV file containing recording data. If None, it will be inferred from the XML file.
        """
        self.xml_file_path = Path(xml_file_path)

        # Load XML data
        with open(self.xml_file_path, "r") as xml_file:
            xml_content = xml_file.read()
            self.xml_recording_dict = xmltodict.parse(xml_content)

        # Extract basic metadata to determine CSV file path
        vrec_session_entry = self.xml_recording_dict["VRecSessionEntry"]

        # Determine recording CSV file path if not provided
        if recording_csv_file_path is None:
            recording_data_file_path_stem = vrec_session_entry["DataFile"]
            self.recording_csv_file_path = self.xml_file_path.parent / f"{recording_data_file_path_stem}.csv"
        else:
            self.recording_csv_file_path = Path(recording_csv_file_path)

        if not self.recording_csv_file_path.is_file():
            raise FileNotFoundError(f"Recording data file not found: {self.recording_csv_file_path}")

    def get_metadata(self) -> DeepDict:
        """
        Extract metadata from the XML file.

        Returns
        -------
        DeepDict
            Dictionary containing metadata
        """
        metadata = super().get_metadata()

        # Extract basic metadata
        vrec_session_entry = self.xml_recording_dict["VRecSessionEntry"]
        experiment_metadata = vrec_session_entry["Experiment"]
        acquisition_time = experiment_metadata["AcquisitionTime"]
        sampling_rate = float(experiment_metadata["Rate"])

        # Extract signal information
        signal_list = experiment_metadata["SignalList"]
        assert len(signal_list) == 1, "Only one signal is expected in the XML file"
        signal_list = signal_list["VRecSignal"]
        enabled_signals = [signal for signal in signal_list if signal["Enabled"] == "true"]
        expected_signal = enabled_signals[0]

        # Extract unit metadata
        unit_metadata = expected_signal["Unit"]
        unit = unit_metadata["UnitName"]
        device = unit_metadata["PatchclampDevice"]

        # Extract session start time
        session_start_time = None
        if "DateTime" in vrec_session_entry:
            datetime_str = vrec_session_entry["DateTime"]
            try:
                # Parse ISO format datetime string
                session_start_time = datetime.fromisoformat(datetime_str.replace("Z", "+00:00"))
            except (ValueError, AttributeError):
                # If parsing fails, leave as None
                pass

        # Add recording-specific metadata
        metadata["Ecephys"] = {
            "Device": [{"name": device, "description": f"Prairie View intracellular recording device: {device}"}],
            "IntracelllularSeries": [
                {
                    "name": "PrairieViewIntracellularRecording",
                    "description": f"Intracellular recording from {device}",
                    "unit": unit,
                    "rate": sampling_rate,
                }
            ],
        }

        # Add session start time if available
        if session_start_time is not None:
            metadata["NWBFile"] = {"session_start_time": session_start_time}

        return metadata

    def add_to_nwbfile(self, nwbfile: NWBFile, metadata: dict, electrode=None):
        """
        Add the recording data to the NWB file.

        Parameters
        ----------
        nwbfile : NWBFile
            NWB file to add the data to
        metadata : dict
            Metadata dictionary
        electrode : optional
            Electrode to use for the recording. If None, a new electrode will be created.
        """
        # Extract metadata from XML
        vrec_session_entry = self.xml_recording_dict["VRecSessionEntry"]
        experiment_metadata = vrec_session_entry["Experiment"]
        signal_list = experiment_metadata["SignalList"]["VRecSignal"]
        enabled_signals = [signal for signal in signal_list if signal["Enabled"] == "true"]
        expected_signal = enabled_signals[0]
        gain = float(expected_signal["Gain"])
        unit_metadata = expected_signal["Unit"]
        multiplier = float(unit_metadata["Multiplier"])
        divisor = float(unit_metadata["Divisor"])
        unit = unit_metadata["UnitName"]
        device_name = unit_metadata["PatchclampDevice"]

        # Load recording data
        recording_data_df = pd.read_csv(
            self.recording_csv_file_path,
            header=0,
        )
        timestamps = recording_data_df["Time(ms)"] / 1_000.0  # Convert to seconds
        voltage_mV = recording_data_df[" Primary"] * multiplier / divisor

        # Create or get device
        if device_name in nwbfile.devices:
            device = nwbfile.devices[device_name]
        else:
            device = nwbfile.create_device(
                name=device_name, description=f"Prairie View intracellular recording device: {device_name}"
            )

        # Create or use provided electrode
        if electrode is None:
            # Create electrode
            electrode = nwbfile.create_icephys_electrode(
                name="PrairieViewElectrode",
                description=f"Prairie View electrode for {device_name}",
                device=device,
            )

        # Create current clamp series
        current_clamp_series = CurrentClampSeries(
            name="PrairieViewIntracellularRecording",
            data=voltage_mV.values,
            timestamps=timestamps.values,
            electrode=electrode,
            unit=unit,
            description=f"Intracellular recording from {device_name}",
        )

        # Add current clamp series to acquisition
        nwbfile.add_acquisition(current_clamp_series)
