"""Complete implementation for PrairieViewBrightnessOverTimeInterface."""

from pathlib import Path

import numpy as np
import pandas as pd
import xmltodict
from neuroconv.converters import BrukerTiffSinglePlaneConverter
from pynwb.ophys import (
    Fluorescence,
    ImageSegmentation,
    RoiResponseSeries,
)


class PrairieViewBrightnessOverTimeInterface(BrukerTiffSinglePlaneConverter):
    """Interface for Prairie View Brightness Over Time (BOT) data.

    This interface extends the BrukerTiffSinglePlaneConverter to handle
    brightness over time measurements from Prairie View. It automatically:

    1. Loads the imaging data (TIFF files) using the parent class
    2. Finds and loads the corresponding BOT CSV file with brightness measurements
    3. Extracts region information from XML metadata
    4. Creates ROI segmentation and fluorescence time series for each region

    The BOT data typically contains timestamps and brightness values for multiple
    regions of interest, which are added as ROI response series to the NWB file.
    """

    ExtractorName = "BrukerTiffSinglePlaneImagingExtractor"

    def __init__(self, folder_path, bot_csv_data_file_path, xml_metadata_file_path):
        """Initialize the interface."""
        super().__init__(folder_path=folder_path)

        self.folder_path = Path(folder_path)

        self.bot_csv_data_file_path = bot_csv_data_file_path
        self.xml_metadata_file_path = xml_metadata_file_path

        # Load XML metadata to get region information
        self._load_roi_region_data()

    def _load_roi_region_data(self):

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

    def add_to_nwbfile(self, nwbfile, metadata=None, conversion_options=None):
        """Add brightness over time data to the NWB file."""
        # First add the imaging data using parent class

        metadata = metadata or self.get_metadata()
        conversion_options = conversion_options or {}
        super().add_to_nwbfile(nwbfile, metadata, **conversion_options)

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

        # Create image segmentation
        image_segmentation = ImageSegmentation()
        ophys_module.add(image_segmentation)

        for interface_name, interface in self.data_interface_objects.items():
            channel_name = interface._stream_name

            # Get the imaging plane created by parent class
            imaging_plane = nwbfile.imaging_planes[f"ImagingPlane{channel_name}"]

            # Create plane segmentation
            plane_segmentation = image_segmentation.create_plane_segmentation(
                name=f"PlaneSegmentation{channel_name}",
                description="Segmentation of regions for brightness over time analysis.",
                imaging_plane=imaging_plane,
            )

            # Get image dimensions from the TwoPhotonSeries
            two_photon_series = nwbfile.acquisition[f"TwoPhotonSeries{channel_name}"]
            dimension = two_photon_series.dimension

            # Create Image Mask for the region

            region_metadata = self.regions_info[channel_name]
            x = int(region_metadata["x"])
            y = int(region_metadata["y"])
            width = int(region_metadata["width"])
            height = int(region_metadata["height"])

            # Create image mask for the specific region
            image_mask = np.zeros(shape=dimension, dtype=bool)
            image_mask[y : y + height, x : x + width] = True

            # image_mask = np.ones(shape=dimension, dtype=bool)

            # Add ROI to plane segmentation
            plane_segmentation.add_roi(image_mask=image_mask)

            # Create ROI table region for this ROI
            roi_table_region = plane_segmentation.create_roi_table_region(
                region=[0], description=f"{region_metadata['region']} for brightness over time analysis"
            )

            # Create fluorescence container
            fluoresence_name = f"Fluorescence"
            if fluoresence_name in ophys_module.data_interfaces:
                fluorescence = ophys_module[fluoresence_name]
            else:
                fluorescence = Fluorescence(name=fluoresence_name)
                ophys_module.add(fluorescence)

            region_data = bot_df[region_metadata["region"]].to_numpy()

            roi_response_series = RoiResponseSeries(
                name=f"BrightnessOverTime{channel_name}",
                description=f"Brightness over time measurements for {region_metadata['region']}",
                data=region_data,
                rois=roi_table_region,
                unit="a.u.",  # Arbitrary units for fluorescence
                timestamps=timestamps,
            )

            fluorescence.add_roi_response_series(roi_response_series)
