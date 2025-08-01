import re
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
import xmltodict
from neuroconv.basedatainterface import BaseDataInterface
from neuroconv.utils import calculate_regular_series_rate
from pynwb.ophys import (
    Fluorescence,
    ImageSegmentation,
    RoiResponseSeries,
)


class PrairieViewFluorescenceInterface(BaseDataInterface):
    """Interface for Prairie View Brightness Over Time (BOT) acetylcholine fluorescence data.

    This interface handles brightness over time measurements from Prairie View for
    acetylcholine GRAB sensor experiments. It processes:

    1. BOT CSV file with brightness measurements
    2. XML metadata with region information
    3. Creates ROI segmentation and fluorescence time series for each region

    The BOT data contains timestamps and brightness values for multiple
    regions of interest, which are added as ROI response series to the NWB file.
    """

    def __init__(self, bot_csv_data_file_path, xml_metadata_file_path, image_dimensions=None):
        """Initialize the interface.

        Parameters
        ----------
        bot_csv_data_file_path : str or Path
            Path to the BOT CSV file with brightness measurements
        xml_metadata_file_path : str or Path
            Path to the XML metadata file with region information
        image_dimensions : tuple, optional
            Image dimensions (height, width) for creating ROI masks.
            If not provided, will use (512, 512) as default.
        """
        super().__init__()

        self.bot_csv_data_file_path = Path(bot_csv_data_file_path)
        self.xml_metadata_file_path = Path(xml_metadata_file_path)
        self.image_dimensions = image_dimensions or (512, 512)

        # Load XML metadata to get region information
        self._load_roi_region_data()

        # Find background reference files
        self._find_reference_files()

    def _find_reference_files(self):
        """Find background reference TIFF files in the References directory."""
        self.reference_files = {}

        # Look for References directory in the same folder as the BOT CSV file
        references_dir = self.bot_csv_data_file_path.parent / "References"

        if references_dir.exists():
            # Find all reference TIFF files, prioritizing uncompressed 16-bit files
            reference_tiffs = list(references_dir.glob("*Reference.tif"))

            for ref_file in reference_tiffs:
                # Parse filename to extract channel/window information
                filename = ref_file.stem

                # Include all reference files, handling both compressed and uncompressed
                if "Ch2" in filename:
                    if "16bit" in filename:
                        self.reference_files["Ch2_16bit"] = ref_file
                    elif "8bit" in filename and "Window1" in filename:
                        self.reference_files["Ch2_8bit_Window1"] = ref_file
                elif "Ch3" in filename and "16bit" in filename:
                    self.reference_files["Ch3_16bit"] = ref_file
                elif "Dodt" in filename and "Window2" in filename:
                    self.reference_files["Dodt_Window2"] = ref_file
        else:
            # No references directory found
            self.reference_files = {}

    def get_session_start_time(self):
        """Get session start time from XML metadata.

        Returns
        -------
        datetime
            Session start time extracted from XML metadata
        """
        with open(self.xml_metadata_file_path, "r", encoding="utf-8") as f:
            content = f.read()
            match = re.search(r'<PVScan.*?date="(.*?)"', content)
            if match:
                date_str = match.group(1)
                try:
                    dt_object = datetime.strptime(date_str, "%m/%d/%Y %I:%M:%S %p")
                    return dt_object
                except ValueError:
                    raise ValueError(f"Could not parse date string: {date_str} in file: {self.xml_metadata_file_path}")
            else:
                raise ValueError(f"No date found in XML metadata file: {self.xml_metadata_file_path}")

    def _load_roi_region_data(self):
        """Load region information from XML metadata."""
        with open(self.xml_metadata_file_path, "r") as f:
            xml_content = f.read()
            self.xml_metadata = xmltodict.parse(xml_content)

        sequence_metadata = self.xml_metadata["PVScan"]["Sequence"]
        first_frame = sequence_metadata["Frame"][0]
        file_metadata = first_frame["File"]
        channel_index_to_channel_name = {m["@channel"]: m["@channelName"] for m in file_metadata}
        brightness_over_time_metadata = sequence_metadata["PVBOTs"]
        region_metadata_list = brightness_over_time_metadata["Region"]

        # Extract region information from XML
        self.regions_info = {}
        for index, region in enumerate(region_metadata_list):
            channel_index = region["@channel"]
            channel_name = channel_index_to_channel_name[channel_index]
            self.regions_info[channel_name] = {
                "x": int(region["@x"]),
                "y": int(region["@y"]),
                "width": int(region["@width"]),
                "height": int(region["@height"]),
                "region": f"Region {index + 1}",
            }

    def add_to_nwbfile(self, nwbfile, metadata=None, trial_id=None):
        """Add brightness over time data to the NWB file.

        Parameters
        ----------
        nwbfile : NWBFile
            The NWB file to add data to
        metadata : dict, optional
            Additional metadata (not used by this interface)
        trial_id : str, optional
            Unique trial identifier for multi-trial sessions
        """
        # Load brightness over time data
        bot_df = pd.read_csv(self.bot_csv_data_file_path)

        # Validate required columns
        if "Timestamp" not in bot_df.columns:
            raise ValueError(
                f"'Timestamp' column not found in {self.bot_csv_data_file_path}. "
                "Ensure the CSV file contains a 'Timestamp' column."
            )

        # Extract timestamps and region data
        timestamps = bot_df["Timestamp"].to_numpy()

        # Get available region columns
        region_columns = [col for col in bot_df.columns if col.startswith("Region")]

        if not region_columns:
            raise ValueError(
                f"No region columns found in {self.bot_csv_data_file_path}. "
                "Ensure the CSV file contains columns starting with 'Region'."
            )

        # Create or get ophys processing module
        if "ophys" not in nwbfile.processing:
            ophys_module = nwbfile.create_processing_module(
                name="ophys",
                description="optical physiology processed data",
            )
        else:
            ophys_module = nwbfile.processing["ophys"]

        # Create or get image segmentation
        if "ImageSegmentation" in ophys_module.data_interfaces:
            image_segmentation = ophys_module["ImageSegmentation"]
        else:
            image_segmentation = ImageSegmentation()
            ophys_module.add(image_segmentation)

        # Create unique naming suffix if trial_id provided
        name_suffix = trial_id if trial_id else ""

        # Process each channel/region
        for channel_name, region_metadata in self.regions_info.items():
            # For fluorescence data, we need an imaging plane reference
            # Since we don't have raw imaging data in this interface,
            # we'll create a placeholder imaging plane
            if f"ImagingPlane{channel_name}" not in nwbfile.imaging_planes:
                # Get device from NWB file
                if "default" in nwbfile.devices:
                    device = nwbfile.devices["default"]
                else:
                    # Create default device if not exists
                    from pynwb.device import Device

                    device = Device(
                        name="default",
                        description="Bruker two-photon microscope for acetylcholine GRAB biosensor imaging",
                    )
                    nwbfile.add_device(device)

                # Create optical channel
                from pynwb.ophys import OpticalChannel

                optical_channel = OpticalChannel(
                    name=f"OpticalChannel{channel_name}",
                    description=f"Optical channel for {channel_name} fluorescence detection",
                    emission_lambda=520.0,  # GRABACh3.0 emission wavelength
                )

                # Create imaging plane for fluorescence data
                imaging_plane = nwbfile.create_imaging_plane(
                    name=f"ImagingPlane{channel_name}",
                    description=f"Imaging plane for acetylcholine fluorescence measurements in {channel_name}",
                    optical_channel=[optical_channel],
                    device=device,
                    excitation_lambda=920.0,  # Two-photon excitation wavelength
                    imaging_rate=21.26,  # From protocol specifications
                    indicator="GRABACh3.0",
                    location="striatum",
                )
            else:
                imaging_plane = nwbfile.imaging_planes[f"ImagingPlane{channel_name}"]

            # Better channel names for consistency
            if channel_name == "Ch2":
                channel_display = "GRABCh"
            elif channel_name == "Dodt":
                channel_display = "DoDT"
            else:
                channel_display = channel_name

            # Create plane segmentation with unique naming
            # Follow same pattern as TwoPhotonSeries: PlaneSegmentation + name_suffix + channel + Trial
            if "Trial" in name_suffix:
                base_part, trial_part = name_suffix.rsplit("Trial", 1)
                plane_seg_name = f"PlaneSegmentation{base_part}{channel_display}Trial{trial_part}"
            else:
                plane_seg_name = f"PlaneSegmentation{name_suffix}{channel_display}"
            plane_segmentation = image_segmentation.create_plane_segmentation(
                name=plane_seg_name,
                description="Segmentation of regions for acetylcholine brightness over time analysis.",
                imaging_plane=imaging_plane,
            )

            # Create ROI mask for the region
            region_metadata = self.regions_info[channel_name]
            x = int(region_metadata["x"])
            y = int(region_metadata["y"])
            width = int(region_metadata["width"])
            height = int(region_metadata["height"])

            # Create image mask for the specific region
            image_mask = np.zeros(shape=self.image_dimensions, dtype=bool)
            image_mask[y : y + height, x : x + width] = True

            # Add ROI to plane segmentation
            plane_segmentation.add_roi(image_mask=image_mask)

            # Create ROI table region for this ROI
            roi_table_region = plane_segmentation.create_roi_table_region(
                region=[0], description=f"{region_metadata['region']} for acetylcholine brightness over time analysis"
            )

            # Create fluorescence container
            fluorescence_name = f"Fluorescence"
            if fluorescence_name in ophys_module.data_interfaces:
                fluorescence = ophys_module[fluorescence_name]
            else:
                fluorescence = Fluorescence(name=fluorescence_name)
                ophys_module.add(fluorescence)

            region_data = bot_df[region_metadata["region"]].to_numpy()

            # Use neuroconv pattern: check if timestamps are regular, use rate+starting_time or timestamps
            rate = calculate_regular_series_rate(series=timestamps)
            recording_t_start = timestamps[0]

            # Create ROI response series with unique naming
            # Follow same pattern as TwoPhotonSeries: Fluorescence + name_suffix + channel + Trial
            # name_suffix already contains the base pattern like "ULControlControlSinglePulseTrial004"
            # We need to insert the channel before "Trial"
            if "Trial" in name_suffix:
                base_part, trial_part = name_suffix.rsplit("Trial", 1)
                series_name = f"AcetylcholineFluorescence{base_part}{channel_display}Trial{trial_part}"
            else:
                series_name = f"AcetylcholineFluorescence{name_suffix}{channel_display}"

            if rate is not None:
                # Regular timestamps - use starting_time + rate
                roi_response_series = RoiResponseSeries(
                    name=series_name,
                    description=f"Acetylcholine fluorescence brightness over time measurements for {region_metadata['region']}",
                    data=region_data,
                    rois=roi_table_region,
                    unit="a.u.",  # Arbitrary units for fluorescence
                    starting_time=float(recording_t_start),
                    rate=rate,
                )
            else:
                # Irregular timestamps - use explicit timestamps
                roi_response_series = RoiResponseSeries(
                    name=series_name,
                    description=f"Acetylcholine fluorescence brightness over time measurements for {region_metadata['region']}",
                    data=region_data,
                    rois=roi_table_region,
                    unit="a.u.",  # Arbitrary units for fluorescence
                    timestamps=timestamps,
                )

            fluorescence.add_roi_response_series(roi_response_series)
