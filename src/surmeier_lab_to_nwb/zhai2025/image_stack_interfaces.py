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
        """
        Extract the session start time from the first frame of the ND2 image stack.

        This method retrieves timing metadata from the first frame of the Nikon ND2 image file
        and converts it to a standard datetime object. The timing information is extracted from
        the frame metadata which contains both an absolute Julian day number and a relative
        time offset in milliseconds within that day.

        The conversion process involves:
        1. Retrieving metadata from the first frame (index 0)
        2. Extracting timing information from the specified channel
        3. Getting the Julian day number (absolute reference point)
        4. Getting the relative time offset in milliseconds within that day
        5. Converting milliseconds to fractional days
        6. Combining Julian day + fractional day offset
        7. Converting from Julian date format to UTC datetime

        Returns
        -------
        datetime
            The session start time as a timezone-aware datetime object in UTC.
            This represents the exact timestamp when the first frame of the
            image stack was acquired.

        Notes
        -----
        The method assumes that:
        - The ND2 file contains at least one frame
        - The specified channel index exists in the frame metadata
        - The timing metadata is available and valid

        The Julian date system is used as an intermediate format because it provides
        a continuous count of days since a fixed epoch, making it suitable for
        high-precision astronomical and scientific timing calculations.

        This implementation uses the astropy.time.Time module to handle the Julian
        date to datetime conversion, avoiding potential errors in manual date/time

        Examples
        --------
        >>> interface = NikonImageStackInterface("path/to/image.nd2")
        >>> start_time = interface.get_session_start_time()
        >>> print(start_time)
        2023-07-15 14:30:25.123000+00:00
        """
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
