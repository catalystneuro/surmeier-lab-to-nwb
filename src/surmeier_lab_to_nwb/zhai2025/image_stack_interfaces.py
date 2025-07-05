from datetime import datetime
from pathlib import Path
from typing import Optional
from zoneinfo import ZoneInfo

from astropy.time import Time
from neuroconv.basedatainterface import BaseDataInterface
from neuroconv.utils import DeepDict
from pydantic import FilePath, validate_call
from pynwb import NWBFile
from pynwb.base import Images
from pynwb.image import GrayscaleImage


class NikonImageStackInterface(BaseDataInterface):
    """Data interface for Nikon image stacks."""

    @validate_call
    def __init__(self, file_path: FilePath):
        """
        Initialize the NikonImageStackInterface.

        Parameters
        ----------
        file_path : FilePath
            Path to the image stack file.
        """
        super().__init__()
        self.file_path = Path(file_path)

        import nd2

        self.file_handle = nd2.ND2File(file_path)

        number_of_experiments = len(self.file_handle.experiment)
        assert number_of_experiments == 1, "Only single-experiment ND2 files are supported."

        self.experiment_info = self.file_handle.experiment[0]
        assert (
            self.experiment_info.type == "ZStackLoop"
        ), "This interface is meant for ZStackLoop image stack acquisition."

        number_of_channels = self.file_handle.attributes.channelCount
        assert number_of_channels == 1, "Only single-channel image stacks are supported."
        self.channel_index = 0

        self.width_pixels = self.file_handle.attributes.widthPx
        self.height_pixels = self.file_handle.attributes.heightPx
        self.number_of_z_planes = self.file_handle.attributes.sequenceCount

        self.dimension_size = {dimension: axis_length for dimension, axis_length in self.file_handle.sizes.items()}
        self.dimension_to_array_axis = {dimension: axis for axis, dimension in enumerate(self.dimension_size)}

    def get_metadata(self) -> DeepDict:

        metadata = super().get_metadata()

        session_start_time = self._get_session_start_time()
        metadata["NWBFile"]["session_start_time"] = session_start_time

        return metadata

    def add_to_nwbfile(self, nwbfile: NWBFile, metadata: Optional[dict] = None):

        image_array = self.file_handle.asarray()

        depth_axis = self.dimension_to_array_axis["Z"]
        assert depth_axis == 0
        number_of_images = image_array.shape[0]

        # Create Images container
        container_name = "ZStackImagesContainer"
        description = (
            "High-resolution confocal images are captured from fixed brain sections "
            "using a 60× oil immersion objective (NA = 1.49), with z-steps of 0.125 μm "
            "and a pixel size of 0.09 μm."
        )
        images_container = Images(
            name=container_name,
            description=description,
        )
        for image_frame_index in range(number_of_images):
            image_data = image_array[image_frame_index, ...]

            nwb_image = GrayscaleImage(
                name=f"ZStackImage{image_frame_index}",
                data=image_data,
            )
            images_container.add_image(nwb_image)

        nwbfile.add_acquisition(images_container)

    def get_session_start_time(self) -> datetime:

        frame_index = 0
        first_frame_metadata = self.file_handle.frame_metadata(frame_index)
        time_info_first_frame = first_frame_metadata.channels[self.channel_index].time

        julian_day_number = time_info_first_frame.absoluteJulianDayNumber
        relative_time_milliseconds = time_info_first_frame.relativeTimeMs

        seconds_in_a_day = 86400.0
        milliseconds_in_a_day = seconds_in_a_day * 1000.0
        first_frame_offset = relative_time_milliseconds / milliseconds_in_a_day  # Convert milliseconds to days

        dt = Time(julian_day_number + first_frame_offset, format="jd").to_datetime(ZoneInfo("UTC"))
        return dt

    def __del__(self):
        """Clean up resources."""
        if hasattr(self, "file_handle"):
            self.file_handle.close()
