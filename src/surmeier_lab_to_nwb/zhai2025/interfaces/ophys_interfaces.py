"""Optical physiology interfaces for converting Prairie View data.

Only ``BrukerReferenceImagesInterface`` remains here. The Bruker BOT and line-scan
flows have moved to dedicated portable interface modules in this package
(:mod:`bruker_bot_segmentation_interface`, :mod:`bruker_line_scan_interface`,
:mod:`bruker_line_scan_kymograph_interface`, :mod:`bruker_line_scan_source_image_interface`).
"""

from pathlib import Path
from typing import Any, Optional

import numpy as np
from neuroconv.datainterfaces import ImageInterface
from neuroconv.utils import DeepDict
from pynwb.file import NWBFile


class BrukerReferenceImagesInterface(ImageInterface):
    """
    Interface for handling Bruker reference images used for background subtraction.

    This interface extends ImageInterface to handle reference images found in
    "References/" folders within experimental data directories. These images are
    PMT background measurements taken with zero laser power for background subtraction.

    The interface automatically handles:
    - Multi-channel image conversion (RGBA/RGB to grayscale)
    - Metadata generation with proper descriptions
    - Consistent naming and organization in NWB files

    Used in both:
    - Figure 5 acetylcholine biosensor experiments
    - Dendritic excitability experiments

    Parameters
    ----------
    references_folder_path : Path or str
        Path to the "References/" folder containing reference image files
    container_name : str, optional
        Name for the Images container in the NWB file.
        Default: "BackgroundReferences"
    images_location : str, optional
        Where to place images in NWB file ('acquisition' or 'stimulus').
        Default: 'acquisition'
    verbose : bool, optional
        Whether to print detailed information during processing.
        Default: False
    """

    def __init__(
        self,
        references_folder_path: str | Path,
        container_name: str = "BackgroundReferences",
        images_location: str = "acquisition",
        verbose: bool = False,
    ):
        # Convert to Path object
        self.references_folder = Path(references_folder_path)

        if not self.references_folder.exists():
            raise FileNotFoundError(f"References folder does not exist: {self.references_folder}")

        if not self.references_folder.is_dir():
            raise ValueError(f"References path is not a directory: {self.references_folder}")

        # Find all TIFF files in the references folder
        tiff_extensions = [".tif", ".tiff"]
        reference_files = []

        for ext in tiff_extensions:
            reference_files.extend(self.references_folder.glob(f"*{ext}"))
            reference_files.extend(self.references_folder.glob(f"*{ext.upper()}"))

        if not reference_files:
            raise FileNotFoundError(f"No TIFF files found in references folder: {self.references_folder}")

        # Sort files for consistent ordering
        reference_files.sort()

        # Store reference file information for metadata generation
        self.reference_file_info = self._parse_reference_files(reference_files)

        # Initialize parent ImageInterface with the reference files
        super().__init__(
            file_paths=reference_files,
            metadata_key=container_name,
            images_location=images_location,
            verbose=verbose,
        )

        self.container_name = container_name

    def _parse_reference_files(self, reference_files: list[Path]) -> dict[str, dict[str, Any]]:
        """
        Parse reference file names to extract metadata information.

        Expected formats:
        - Dendritic excitability: *-Ch1-16bit-Reference.tif, *-Ch2-16bit-Reference.tif, etc.
        - Figure 5 (if applicable): Similar channel-based naming

        Parameters
        ----------
        reference_files : list[Path]
            List of reference file paths

        Returns
        -------
        dict[str, dict[str, Any]]
            Dictionary mapping file paths to metadata information
        """
        file_info = {}

        for file_path in reference_files:
            filename = file_path.name

            # Parse filename to extract channel, bit depth, and type information
            info = {
                "original_filename": filename,
                "channel": "unknown",
                "bit_depth": "unknown",
                "image_type": "reference",
                "description": f"Background reference image from {filename}. PMT background measured with zero laser power for background subtraction.",
            }

            # Extract channel information
            if "-Ch1-" in filename:
                info["channel"] = "Ch1"
            elif "-Ch2-" in filename:
                info["channel"] = "Ch2"
            elif "-Ch3-" in filename:
                info["channel"] = "Ch3"

            # Extract bit depth information
            if "16bit" in filename:
                info["bit_depth"] = "16bit"
            elif "8bit" in filename:
                info["bit_depth"] = "8bit"

            # Create descriptive name following the pattern: Image{Ch1|Ch2}{16|8}BitReference{Window1|Window2}
            if info["channel"] != "unknown" and info["bit_depth"] != "unknown":
                # Add window info if present
                window_suffix = ""
                if "Window" in filename:
                    import re

                    window_match = re.search(r"Window\d+", filename)
                    if window_match:
                        window_suffix = window_match.group()

                # Create name in format: Image{Ch1|Ch2}{16|8}BitReference{Window1|Window2}
                bit_depth_clean = info["bit_depth"].replace("bit", "")  # Remove 'bit' suffix
                info["display_name"] = f"Image{info['channel']}{bit_depth_clean}BitReference{window_suffix}"
                info["description"] = (
                    f"Background reference image for {info['channel']} ({info['bit_depth']}). PMT background measured with zero laser power for background subtraction."
                )
            else:
                # Fallback for non-standard naming - use original filename to ensure uniqueness
                base_name = filename.replace(".tif", "").replace(".tiff", "")
                info["display_name"] = f"Image{base_name}Reference"

            # Handle special cases like window references
            if "Window" in filename:
                info["description"] += " (Window-specific reference image)"
                info["image_type"] = "window_reference"

            file_info[str(file_path)] = info

        return file_info

    def get_metadata(self) -> DeepDict:
        """
        Get metadata for the reference images with custom image information.

        Returns
        -------
        DeepDict
            Metadata dictionary with reference image information
        """
        # Get base metadata from parent ImageInterface
        metadata = super().get_metadata()

        # Create custom metadata for individual reference images
        images_metadata = {}

        for file_path_str, file_info in self.reference_file_info.items():
            images_metadata[file_path_str] = {
                "name": file_info["display_name"],
                "description": file_info["description"],
            }

        # Update the Images container metadata with custom image information
        if "Images" in metadata and self.container_name in metadata["Images"]:
            metadata["Images"][self.container_name].update(
                {
                    "name": self.container_name,
                    "description": (
                        "Background reference images used for PMT background subtraction. "
                        "These images were acquired with zero laser power to measure "
                        "baseline PMT response for background correction as described in methods."
                    ),
                    "images": images_metadata,
                }
            )

        return metadata

    def add_to_nwbfile(
        self,
        nwbfile: NWBFile,
        metadata: Optional[dict[str, Any]] = None,
        qualifier: str = "",
        repetition_number: Optional[int] = None,
    ) -> None:
        """
        Add reference images to the NWB file in a shared container.

        Parameters
        ----------
        nwbfile : NWBFile
            The NWB file to add the reference images to
        metadata : dict[str, Any], optional
            Metadata for the images. If None, will use get_metadata()
        qualifier : str, optional
            Qualifier to prepend to image names (e.g., condition+treatment+stimulation)
        repetition_number : int, optional
            Repetition number to append to image names in 3-digit format
        """
        import tifffile
        from pynwb.image import Image, Images

        # Use provided metadata or get default
        if metadata is None:
            metadata = self.get_metadata()

        # Check if shared container already exists
        container_name = "ImagesBackgroundReferences"
        if container_name in nwbfile.acquisition:
            images_container = nwbfile.acquisition[container_name]
        else:
            # Create shared container
            images_container = Images(
                name=container_name,
                description=(
                    "Background reference images used for PMT background subtraction across all trials. "
                    "These images were acquired with zero laser power to measure "
                    "baseline PMT response for background correction as described in methods."
                ),
            )
            nwbfile.add_acquisition(images_container)

        # Process and add each reference image with new naming convention
        for file_path_str, file_info in self.reference_file_info.items():
            file_path = Path(file_path_str)

            # Load image data
            image_data = tifffile.imread(file_path)

            # Handle multi-channel images by converting to grayscale
            if len(image_data.shape) > 2:
                # Convert to grayscale if multi-channel (RGB or RGBA)
                if image_data.shape[2] == 3:  # RGB
                    image_data = np.mean(image_data, axis=2)
                elif image_data.shape[2] == 4:  # RGBA
                    image_data = np.mean(image_data[:, :, :3], axis=2)  # Ignore alpha channel
                image_data = image_data.astype(np.float32)

            # Extract channel information for display name
            channel_info = file_info.get("channel", "unknown")
            channel_display = channel_info  # Ch1, Ch2, etc.

            # Determine image type
            bit_depth = file_info.get("bit_depth", "unknown")
            if bit_depth == "16bit":
                type_suffix = "16BitReference"
            elif bit_depth == "8bit":
                type_suffix = "8BitReference"
            else:
                type_suffix = "Reference"

            # Build new image name following pattern: Image{qualifier}{channel_to_display}{type_of_image}{repetition_if_exists}
            image_name = f"Image{qualifier}{channel_display}{type_suffix}"
            if repetition_number is not None:
                image_name += f"Repetition{repetition_number:03d}"

            # Create image with new name
            image = Image(
                name=image_name,
                data=image_data,
                description=file_info["description"],
            )

            # Add to shared container
            images_container.add_image(image)

    @classmethod
    def from_trial_folder(
        cls, trial_folder: str | Path, container_name: Optional[str] = None, **kwargs
    ) -> "BrukerReferenceImagesInterface":
        """
        Create BrukerReferenceImagesInterface from a trial folder containing a References subfolder.

        Parameters
        ----------
        trial_folder : Path or str
            Path to the trial folder containing a "References/" subfolder
        container_name : str, optional
            Name for the Images container. If None, will generate from trial folder name
        **kwargs
            Additional arguments passed to the constructor

        Returns
        -------
        BrukerReferenceImagesInterface
            Interface instance for the reference images
        """
        trial_path = Path(trial_folder)
        references_path = trial_path / "References"

        if not references_path.exists():
            raise FileNotFoundError(f"No References folder found in trial folder: {trial_path}")

        # Generate container name if not provided
        if container_name is None:
            trial_name = trial_path.name
            container_name = f"BackgroundReferences_{trial_name}"

        return cls(references_folder_path=references_path, container_name=container_name, **kwargs)
