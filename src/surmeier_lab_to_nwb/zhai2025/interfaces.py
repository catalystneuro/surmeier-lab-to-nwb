"""Primary classes for converting Prairie View data."""

from datetime import datetime
from pathlib import Path
from typing import Literal, Optional

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

    def __init__(
        self,
        xml_metadata_file_path: str | Path,
        cell_id: str = "CellID",
        trial_id: str = "",
        position: Literal["", "Proximal", "Distal"] = "",
        condition: Literal["", "LIDOffState", "LIDOnState", "LIDOnStateWithSCH"] = "",
    ):
        """
        Initialize the interface with path to the XML metadata file.

        Parameters
        ----------
        xml_metadata_file_path : str | Path
            Path to the XML metadata file containing line scan parameters
        """
        self.xml_metadata_file_path = Path(xml_metadata_file_path)
        self.cell_id = cell_id
        self.trial_id = trial_id
        self.position = position
        self.condition = condition

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
        frame_metadata = sequence_metadata["Frame"]
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
        microscope_name = "BrukerUltima"
        if microscope_name in nwbfile.devices:
            microscope_device = nwbfile.devices[microscope_name]
        else:
            microscope_device = nwbfile.create_device(
                name=microscope_name,
                description="Bruker two-photon microscope",
            )

        # Create optical channels - separate channels for Fluo-4 and Alexa 568
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

        # Extract laser wavelength from metadata if available
        prairie_view_metadata = self.xml_general_metadata_dict["PVScan"]
        prairie_view_state_metadata_list = prairie_view_metadata["PVStateShard"]["PVStateValue"]

        # Create imaging plane
        microns_per_pixel = next(
            (item for item in prairie_view_state_metadata_list if item["@key"] == "micronsPerPixel"), None
        )
        microns_per_pixel_dict = {value["@index"]: value["@value"] for value in microns_per_pixel["IndexedValue"]}
        grid_spacing = [float(microns_per_pixel_dict["XAxis"]), float(microns_per_pixel_dict["YAxis"])]
        grid_spacing_unit = "microns"

        # Get objective lens info
        objective_lens = next(
            (item["@value"] for item in prairie_view_state_metadata_list if item["@key"] == "objectiveLens"),
            "Olympus 60X",
        )
        excitation_lambda = 810.0  # Chameleon Ultra II
        imaging_plane_name = f"dSPNImagingPlaneCondition{self.condition}Cell{self.cell_id}{self.position}"
        available_imaging_planes = nwbfile.imaging_planes
        if imaging_plane_name in available_imaging_planes:
            imaging_plane = available_imaging_planes[imaging_plane_name]
        else:
            imaging_plane = nwbfile.create_imaging_plane(
                name=imaging_plane_name,
                description=f"Line scan imaging of dendrite using {objective_lens} objective",
                optical_channel=optical_channels,
                excitation_lambda=excitation_lambda,
                indicator="Fluo-4",
                location="dorsolateral striatum",
                grid_spacing=grid_spacing,
                grid_spacing_unit=grid_spacing_unit,
                device=microscope_device,
            )

        # Load line scan profile data
        line_scan_data_df = pd.read_csv(self.profile_csv_data_file_path)

        # Extract and clearly label line scan profiles
        # Based on our analysis: Channel 1 is Alexa 568, Channel 2 is Fluo-4
        alexa568_profile = line_scan_data_df[" Prof 1"].to_numpy()  # Structural reference
        fluo4_profile = line_scan_data_df[" Prof 2"].to_numpy()  # Calcium indicator
        timestamps = line_scan_data_df["Prof 1 Time(ms)"].to_numpy(dtype=float) / 1000.0  # Convert to seconds

        # Create ophys processing module if it doesn't exist
        if "ophys" not in nwbfile.processing:
            ophys_module = nwbfile.create_processing_module(
                name="ophys",
                description="optical physiology processed data",
            )
        else:
            ophys_module = nwbfile.processing["ophys"]

        # Create image segmentation
        image_segmentation_name = f"ImageSegmentationCondition{self.condition}Cell{self.cell_id}{self.position}"
        if image_segmentation_name in ophys_module.data_interfaces:
            image_segmentation = ophys_module.data_interfaces[image_segmentation_name]
        else:
            image_segmentation = ImageSegmentation(name=image_segmentation_name)
            ophys_module.add(image_segmentation)

        # Add the line scan data as metadata
        sequence_metadata = prairie_view_metadata["Sequence"]
        line_scan_definition_metadata = sequence_metadata["PVLinescanDefinition"]
        line_metadata = line_scan_definition_metadata["Line"]
        start_pixel_y = int(line_metadata["@startPixelY"])
        start_pixel_x = int(line_metadata["@startPixelX"])
        stop_pixel_x = int(line_metadata["@stopPixelX"])
        line_length = float(line_metadata["@lineLength"])

        # Add additional metadata to the plane segmentation description
        segmentation_description = (
            f"Line scan ROI from pixel ({start_pixel_x},{start_pixel_y}) to ({stop_pixel_x},{start_pixel_y}), "
            f"length: {line_length} microns. Line scan used to measure calcium dynamics during "
            f"back-propagating action potentials in striatal spiny projection neurons."
        )

        # Create plane segmentation with better description
        plane_segmentation_name = (
            f"PlaneSegmentationCondition{self.condition}Cell{self.cell_id}{self.position}Trial{self.trial_id}"
        )
        available_plane_segmentations = image_segmentation.plane_segmentations
        if plane_segmentation_name in available_plane_segmentations:
            plane_segmentation = available_plane_segmentations[plane_segmentation_name]
        else:
            plane_segmentation = image_segmentation.create_plane_segmentation(
                name=plane_segmentation_name,
                description=segmentation_description,
                imaging_plane=imaging_plane,
            )

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

        # Create ROI response series with better names and descriptions
        alexa568_series = RoiResponseSeries(
            name=f"Alexa568StructuralDyeResponseSeriesCondition{self.condition}Cell{self.cell_id}{self.position}Trial{self.trial_id}",
            description="Structural reference fluorescence from Alexa Fluor 568 hydrazide (50 μM)",
            data=alexa568_profile,
            rois=roi_table_region,
            unit="a.u.",  # Arbitrary units for fluorescence
            timestamps=timestamps,
        )

        fluo4_series = RoiResponseSeries(
            name=f"Fluo4CalciumIndicatorResponseSeriesCondition{self.condition}Cell{self.cell_id}{self.position}Trial{self.trial_id}",
            description="Calcium indicator fluorescence from Fluo-4 (100 μM)",
            data=fluo4_profile,
            rois=roi_table_region,
            unit="a.u.",  # Arbitrary units for fluorescence
            timestamps=timestamps,
        )

        # Create fluorescence container
        name = f"FluorescenceContainerCondition{self.condition}Cell{self.cell_id}{self.position}"
        if name in ophys_module.data_interfaces:
            fluorescence = ophys_module.data_interfaces[name]
        else:
            fluorescence = Fluorescence(name=name)
            ophys_module.add(fluorescence)

        fluorescence.add_roi_response_series(alexa568_series)
        fluorescence.add_roi_response_series(fluo4_series)

        # Add fluorescence to ophys module

        # Add source images with better labels
        # Define source image metadata dictionary
        source_image_metadata = {
            "Ch1": {
                "name": f"Alexa568StructuralDyeSourceImageCondition{self.condition}Cell{self.cell_id}{self.position}Trial{self.trial_id}",
                "description": "Source image for Alexa Fluor 568 structural reference",
            },
            "Ch2": {
                "name": f"Fluo4CalciumIndicatorSourceImageCondition{self.condition}Cell{self.cell_id}{self.position}Trial{self.trial_id}",
                "description": "Source image for Fluo-4 calcium indicator",
            },
        }

        name = f"SourceImagesContainerCondition{self.condition}Cell{self.cell_id}"
        if name in nwbfile.acquisition:
            image_container = nwbfile.acquisition[name]
        else:
            image_container = Images(name=name)
            nwbfile.add_acquisition(image_container)

        invalid_channels = ["Dodt"]  # TODO: figure out what this channel is
        for channel_name, source_data_file_path in self.source_data_file_paths_dict.items():
            if channel_name in invalid_channels:
                continue
            source_data = tifffile.imread(source_data_file_path)

            name = source_image_metadata[channel_name]["name"]
            description = source_image_metadata[channel_name]["description"]
            source_image = Image(name=name, data=source_data, description=description)
            image_container.add_image(source_image)

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

        # Define channel metadata dictionary
        channel_metadata = {
            "Ch1": {
                "name": f"Alexa568StructuralDyeLineScanCondition{self.condition}Cell{self.cell_id}{self.position}Trial{self.trial_id}",
                "description": "Line scan raw data for Alexa Fluor 568 structural reference",
            },
            "Ch2": {
                "name": f"Fluo4CalciumIndicatorLineScanCondition{self.condition}Cell{self.cell_id}{self.position}Trial{self.trial_id}",
                "description": "Line scan raw data for Fluo-4 calcium indicator",
            },
        }

        # Add raw line scan data with better labels
        for channel_name, raw_line_scan_data_file_path in self.line_scan_raw_data_file_paths_dict.items():
            raw_line_scan_data = tifffile.imread(raw_line_scan_data_file_path)

            name = channel_metadata[channel_name]["name"]
            description = channel_metadata[channel_name]["description"]
            raw_line_scan_data_series = TwoPhotonSeries(
                name=metadata["name"],
                description=metadata["description"],
                data=raw_line_scan_data,
                imaging_plane=imaging_plane,
                unit="a.u.",
                rate=rate,
                dimension=[pixels_per_line, lines_per_frame],
            )
            nwbfile.add_acquisition(raw_line_scan_data_series)


class PrairieViewIntracellularRecordingInterface(BaseDataInterface):
    """Interface for Prairie View intracellular recording data."""

    keywords = ("intracellular", "recording", "prairie view")

    def __init__(
        self,
        xml_file_path: str | Path,
        recording_csv_file_path: Optional[str | Path] = None,
        cell_id: str = "CellID",
        trial_id: str = "",
        position: Literal["", "Proximal", "Distal"] = "",
        condition: Literal["", "LIDOffState", "LIDOnState", "LIDOnStateWithSCH"] = "",
    ):
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
        self.cell_id = cell_id
        self.trial_id = trial_id
        self.position = position
        self.condition = condition

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

    def add_to_nwbfile(self, nwbfile: NWBFile, metadata: Optional[dict] = None, electrode=None):
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
        # device_name = unit_metadata["PatchclampDevice"]

        # Load recording data
        recording_data_df = pd.read_csv(
            self.recording_csv_file_path,
            header=0,
        )
        timestamps = recording_data_df["Time(ms)"] / 1_000.0  # Convert to seconds
        voltage_mV = recording_data_df[" Primary"] * multiplier / divisor
        voltage_volts = voltage_mV.values / 1_000.0  # Convert mV to volts

        # Create or get device
        device_name = "MultiClamp 700B"
        device_description = f"{device_name} amplifier (Axon Instruments/Molecular Devices) with signals filtered at 2 kHz and digitized at 10 kHz"
        if device_name in nwbfile.devices:
            device = nwbfile.devices[device_name]
        else:
            device = nwbfile.create_device(name=device_name, description=device_description)

        # Create or use provided electrode
        electrode_name = f"IntracellularElectrodeCondition{self.condition}Cell{self.cell_id}{self.position}"
        available_electrodes = nwbfile.icephys_electrodes
        if electrode_name in available_electrodes:
            electrode = available_electrodes[electrode_name]
        else:
            electrode = nwbfile.create_icephys_electrode(
                name=electrode_name,
                description=f"Prairie View electrode for {device_name}",
                device=device,
            )

        # Create current clamp series
        name = f"CurrentClampSeriesCondition{self.condition}Cell{self.cell_id}{self.position}Trial{self.trial_id}"
        current_clamp_series = CurrentClampSeries(
            name=name,
            description=f"Intracellular recording from {device_name}",
            data=voltage_volts,
            timestamps=timestamps.values,
            electrode=electrode,
        )

        # Add current clamp series to acquisition
        nwbfile.add_acquisition(current_clamp_series)
