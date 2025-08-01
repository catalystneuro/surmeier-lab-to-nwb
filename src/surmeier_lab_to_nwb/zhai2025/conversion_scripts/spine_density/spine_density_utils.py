"""
Utility functions for spine density conversion scripts - Zhai et al. 2025
========================================================================

This module contains shared utility functions used across all spine density
conversion scripts for parsing XML metadata files and other common operations.
"""

import re
import uuid
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional
from zoneinfo import ZoneInfo

import xmltodict
from neuroconv.datainterfaces import ImageInterface
from neuroconv.tools.nwb_helpers import make_nwbfile_from_metadata
from neuroconv.utils import dict_deep_update, load_dict_from_file
from pynwb import NWBFile
from pynwb.device import Device

from surmeier_lab_to_nwb.zhai2025.conversion_scripts.conversion_utils import (
    PHARMACOLOGY_ADDITIONS,
    format_condition,
    generate_canonical_session_id,
)


def parse_xml_metadata(xml_file: Path, verbose: bool = False) -> Dict[str, Any]:
    """
    Parse XML metadata file to extract acquisition parameters.

    Parameters
    ----------
    xml_file : Path
        Path to XML metadata file
    verbose : bool, default=False
        Enable verbose output

    Returns
    -------
    Dict[str, Any]
        Dictionary containing acquisition parameters
    """
    if not xml_file.exists():
        raise FileNotFoundError(f"XML metadata file does not exist: {xml_file}")

    # Try xmltodict first, fall back to pure XML parsing
    with open(xml_file, "r", encoding="utf-8", errors="ignore") as f:
        xml_content = f.read()
        # Remove BOM if present
        if xml_content.startswith("\ufeff"):
            xml_content = xml_content[1:]
        # Clean invalid XML characters and HTML entities that can cause parsing errors
        xml_content = re.sub(r"[\x00-\x08\x0B\x0C\x0E-\x1F\x7F]", "", xml_content)
        # Also clean up HTML encoded null characters
        xml_content = re.sub(r"&#x0;", "", xml_content)
        xml_dict = xmltodict.parse(xml_content)

    metadata = {}

    # Check if this is PVScan or AutoQuant format
    if "PVScan" in xml_dict:
        # PVScan format
        pv_scan = xml_dict.get("PVScan", {})
        version = pv_scan.get("@version") or pv_scan.get("version")
        if version:
            metadata["system_version"] = version

        # Navigate to PVStateValue entries
        pv_state_shard = xml_dict.get("PVScan", {}).get("PVStateShard", {})
        pv_state_values = pv_state_shard.get("PVStateValue", [])

        # Handle case where PVStateValue is a single dict or list of dicts
        if isinstance(pv_state_values, dict):
            pv_state_values = [pv_state_values]

        # Extract key acquisition parameters
        for state_value in pv_state_values:
            # Handle both @key and key formats
            key = state_value.get("@key") or state_value.get("key")

            if key == "micronsPerPixel":
                # Extract X and Y pixel sizes from IndexedValue entries
                indexed_values = state_value.get("IndexedValue", [])
                if isinstance(indexed_values, dict):
                    indexed_values = [indexed_values]

                for indexed_value in indexed_values:
                    index = indexed_value.get("@index") or indexed_value.get("index")
                    value = indexed_value.get("@value") or indexed_value.get("value")
                    if index == "XAxis":
                        metadata["pixel_size_x"] = float(value)
                    elif index == "YAxis":
                        metadata["pixel_size_y"] = float(value)

            elif key == "dwellTime":
                value = state_value.get("@value") or state_value.get("value")
                metadata["dwell_time"] = float(value)

            elif key == "objectiveLens":
                value = state_value.get("@value") or state_value.get("value")
                metadata["objective_lens"] = value

            elif key == "objectiveLensNA":
                value = state_value.get("@value") or state_value.get("value")
                metadata["numerical_aperture"] = float(value)

            elif key == "objectiveLensMag":
                value = state_value.get("@value") or state_value.get("value")
                metadata["magnification"] = float(value)

            elif key == "opticalZoom":
                value = state_value.get("@value") or state_value.get("value")
                metadata["optical_zoom"] = float(value)

            elif key == "pixelsPerLine":
                value = state_value.get("@value") or state_value.get("value")
                metadata["pixels_per_line"] = int(value)

            elif key == "linesPerFrame":
                value = state_value.get("@value") or state_value.get("value")
                metadata["lines_per_frame"] = int(value)

            elif key == "bitDepth":
                value = state_value.get("@value") or state_value.get("value")
                metadata["bit_depth"] = int(value)

            elif key == "activeMode":
                value = state_value.get("@value") or state_value.get("value")
                metadata["scan_mode"] = value

    elif "cDatasetXSD" in xml_dict:
        # AutoQuant format
        metadata["system_version"] = "AutoQuant X3.0.4"

        # Extract parameters from ParameterTable entries
        parameter_tables = xml_dict.get("cDatasetXSD", {}).get("ParameterTable", [])
        if isinstance(parameter_tables, dict):
            parameter_tables = [parameter_tables]

        for param_table in parameter_tables:
            name = param_table.get("NameString", "")
            value = param_table.get("ValueString", "")

            if name == "XSpacing" and value:
                try:
                    metadata["pixel_size_x"] = float(value)
                except (ValueError, TypeError):
                    pass
            elif name == "YSpacing" and value:
                try:
                    metadata["pixel_size_y"] = float(value)
                except (ValueError, TypeError):
                    pass
            elif name == "BitDepth" and value:
                try:
                    metadata["bit_depth"] = int(value)
                except (ValueError, TypeError):
                    pass
            elif name == "DWidth" and value:
                try:
                    metadata["pixels_per_line"] = int(value)
                except (ValueError, TypeError):
                    pass
            elif name == "DHeight" and value:
                try:
                    metadata["lines_per_frame"] = int(value)
                except (ValueError, TypeError):
                    pass
            elif name == "DDepth" and value:
                try:
                    metadata["z_slices"] = int(value)
                except (ValueError, TypeError):
                    pass

    return metadata


def find_xml_metadata_file(subfolder: Path) -> Optional[Path]:
    """
    Find XML metadata file in subfolder.

    Parameters
    ----------
    subfolder : Path
        Path to image stack subfolder

    Returns
    -------
    Optional[Path]
        Path to XML metadata file or None if not found
    """
    # Each folder contains exactly one XML file - find it by suffix
    xml_files = [f for f in subfolder.iterdir() if f.suffix == ".xml"]

    if xml_files:
        return xml_files[0]  # Return the single XML file

    return None


def extract_z_slice_info(tiff_file: Path) -> Dict[str, Any]:
    """
    Extract Z-slice information from TIFF filename.

    Example filenames:
    - "20_ZSeries-20170707_Cell2_prox12-001_Cycle00001_Ch1_#.ome_Z026.tif"
    - "20_ZSeries-20190411_Cell1_dist1-001_Cycle00001_Ch1_#.ome_Z01.tif"
    - "ZSeries-20160812_Cell3_dist12-001_Cycle00001_Ch1_000001.ome.tif"

    Parameters
    ----------
    tiff_file : Path
        Path to TIFF file

    Returns
    -------
    Dict[str, Any]
        Dictionary containing Z-slice number, name, and description
    """
    filename = tiff_file.name

    # Extract Z-slice number from filename - handle multiple patterns
    z_slice_num = None

    # Pattern 1: _Z## format (e.g., _Z026.tif)
    z_match = re.search(r"_Z(\d+)\.tif", filename)
    if z_match:
        z_slice_num = int(z_match.group(1))
    else:
        # Pattern 2: ######.ome.tif format (e.g., 000001.ome.tif)
        z_match = re.search(r"_(\d{6})\.ome\.tif", filename)
        if z_match:
            z_slice_num = int(z_match.group(1))
        else:
            raise ValueError(f"Could not extract Z-slice number from filename: {filename}")

    # Extract cell and location info for descriptive name
    cell_match = re.search(r"Cell(\d+)", filename)
    cell_number = cell_match.group(1) if cell_match else "unknown"

    if "dist" in filename:
        location_match = re.search(r"dist([a-z0-9]+)", filename, re.IGNORECASE)
        location = "distal"
        dendrite_num = location_match.group(1) if location_match else "1"
    elif "prox" in filename:
        location_match = re.search(r"prox([a-z0-9]+)", filename, re.IGNORECASE)
        location = "proximal"
        dendrite_num = location_match.group(1) if location_match else "1"
    else:
        location = "unknown"
        dendrite_num = "1"

    # Create descriptive name and description
    z_depth_um = z_slice_num * 0.3  # 0.3 μm z-steps
    image_name = f"ImageZ{z_slice_num:03d}"

    description = (
        f"Z-slice {z_slice_num} at depth {z_depth_um:.1f} μm from cell {cell_number} "
        f"{location} dendrite {dendrite_num}. Acquired with two-photon microscopy "
        f"using 0.15 μm pixels and 0.3 μm z-steps for spine density analysis."
    )

    return {
        "z_slice_number": z_slice_num,
        "z_depth_um": z_depth_um,
        "image_name": image_name,
        "description": description,
        "cell_number": cell_number,
        "location": location,
        "dendrite_number": dendrite_num,
    }


def create_image_metadata(
    tiff_files: List[Path], container_info: Dict[str, str], xml_metadata: Dict[str, Any], verbose: bool = False
) -> Dict[str, Any]:
    """
    Create metadata for ImageInterface with individual image names and descriptions.

    Parameters
    ----------
    tiff_files : List[Path]
        List of TIFF file paths
    container_info : Dict[str, str]
        Container information with name and description
    xml_metadata : Dict[str, Any]
        Metadata extracted from XML file
    verbose : bool, default=False
        Enable verbose output

    Returns
    -------
    Dict[str, Any]
        Metadata dictionary for ImageInterface
    """
    # Sort files by Z-slice number
    tiff_files_with_z = []
    for tiff_file in tiff_files:
        z_info = extract_z_slice_info(tiff_file)
        tiff_files_with_z.append((tiff_file, z_info))

    # Sort by Z-slice number
    tiff_files_with_z.sort(key=lambda x: x[1]["z_slice_number"])

    # Calculate resolution in pixels per cm (converted from μm per pixel)
    pixel_size_x_um = xml_metadata.get("pixel_size_x", 0.15)
    pixel_size_y_um = xml_metadata.get("pixel_size_y", 0.15)
    resolution_x = 10000.0 / pixel_size_x_um  # Convert μm/pixel to pixels/cm
    resolution_y = 10000.0 / pixel_size_y_um
    resolution = (resolution_x + resolution_y) / 2  # Average resolution

    # Create comprehensive description with acquisition parameters
    detailed_description = (
        f"{container_info['description']} "
        f"This Z-stack contains {len(tiff_files_with_z)} individual images acquired at 0.3 μm intervals. "
        f"Acquisition parameters: pixel_size_x={xml_metadata.get('pixel_size_x', 0.15):.6f} μm, "
        f"pixel_size_y={xml_metadata.get('pixel_size_y', 0.15):.6f} μm, "
        f"pixel_size_z=0.3 μm, "
        f"dwell_time={xml_metadata.get('dwell_time', 10.0)} μs, "
        f"optical_zoom={xml_metadata.get('optical_zoom', 5.2)}, "
        f"scan_mode={xml_metadata.get('scan_mode', 'Galvo')}, "
        f"bit_depth={xml_metadata.get('bit_depth', 12)}, "
        f"objective_lens='{xml_metadata.get('objective_lens', 'Olympus 60X')}', "
        f"numerical_aperture={xml_metadata.get('numerical_aperture', 1.0)}, "
        f"magnification={xml_metadata.get('magnification', 60)}x, "
        f"system_version={xml_metadata.get('system_version', 'PVScan')}."
    )

    # Create metadata structure for ImageInterface
    images_metadata = {"description": detailed_description, "images": {}}

    # Add metadata for each individual image
    for tiff_file, z_info in tiff_files_with_z:
        file_path_str = str(tiff_file)
        images_metadata["images"][file_path_str] = {
            "name": z_info["image_name"],
            "description": z_info["description"],
            "resolution": resolution,
        }

    return images_metadata


def create_microscope_device(nwbfile: NWBFile, xml_metadata: Dict[str, Any]) -> Device:
    """
    Create microscope device with specifications from XML metadata.

    Parameters
    ----------
    nwbfile : NWBFile
        NWB file to add device to
    xml_metadata : Dict[str, Any]
        Metadata extracted from XML file

    Returns
    -------
    Device
        Device object for microscope
    """
    # Use metadata from XML or defaults from paper
    objective_lens = xml_metadata.get("objective_lens", "Olympus 60X")
    numerical_aperture = xml_metadata.get("numerical_aperture", 1.0)
    magnification = xml_metadata.get("magnification", 60)
    system_version = xml_metadata.get("system_version", "PVScan")

    device_description = (
        f"Two-photon laser scanning microscope with {objective_lens} objective "
        f"(NA={numerical_aperture}, {magnification}x magnification). "
        f"System: {system_version} (Bruker/Prairie View). "
        f"Used for spine density imaging with 0.15 μm pixels and 0.3 μm z-steps."
    )

    microscope_device = Device(name="TwoPhotonMicroscope", description=device_description, manufacturer="Bruker")

    nwbfile.add_device(microscope_device)
    return microscope_device


class TiffImageStackInterface(ImageInterface):
    """
    Custom ImageInterface for spine density TIFF image stacks.

    This interface extends ImageInterface to:
    1. Automatically filter out projection files
    2. Create microscope device from XML metadata
    3. Generate custom image metadata with Z-slice information
    """

    def __init__(self, subfolder: Path, container_info: Dict[str, str], verbose: bool = False):
        """
        Initialize TiffImageStackInterface.

        Parameters
        ----------
        subfolder : Path
            Path to subfolder containing TIFF files and XML metadata
        container_info : Dict[str, str]
            Container information with name and description
        verbose : bool, default=False
            Enable verbose output
        """
        self.subfolder = subfolder
        self.container_info = container_info
        self.verbose = verbose

        # Get all TIFF files in the folder, excluding projection images
        all_tiff_files = [f for f in subfolder.iterdir() if f.suffix.lower() == ".tif"]

        # Filter out projection images and keep only Z-stack images
        tiff_files = sorted(
            [
                f
                for f in all_tiff_files
                if not (
                    f.name.startswith("_")
                    or "project" in f.name.lower()
                    or f.name.endswith("_xy.tif")
                    or f.name.endswith("_xz.tif")
                    or f.name.endswith("_zy.tif")
                )
            ]
        )

        if len(tiff_files) == 0:
            raise ValueError(f"No Z-stack TIFF files found in {subfolder.name}")

        # Find and parse XML metadata file
        xml_file = find_xml_metadata_file(subfolder)
        if not xml_file:
            raise FileNotFoundError(f"No XML metadata file found in {subfolder.name}")

        self.xml_metadata = parse_xml_metadata(xml_file, verbose=verbose)
        if not self.xml_metadata:
            raise ValueError(f"Failed to parse XML metadata from {xml_file.name}")

        # Initialize parent ImageInterface
        super().__init__(file_paths=tiff_files, images_container_metadata_key=container_info["container_name"])

    def get_metadata(self) -> Dict[str, Any]:
        """
        Get metadata with custom image information.

        Returns
        -------
        Dict[str, Any]
            Metadata dictionary with custom image information
        """
        # Get base metadata from parent
        metadata = super().get_metadata()

        # Create custom metadata for individual images
        images_metadata = create_image_metadata(
            tiff_files=self.file_paths,
            container_info=self.container_info,
            xml_metadata=self.xml_metadata,
            verbose=self.verbose,
        )

        # Update the metadata with custom image information
        metadata["Images"][self.container_info["container_name"]].update(images_metadata)

        return metadata

    def add_to_nwbfile(self, nwbfile: NWBFile, metadata: Optional[Dict[str, Any]] = None) -> None:
        """
        Add image stack to NWB file and create microscope device.

        Parameters
        ----------
        nwbfile : NWBFile
            NWB file to add data to
        metadata : Optional[Dict[str, Any]], default=None
            Metadata dictionary
        """
        # Create microscope device if it doesn't exist
        if "TwoPhotonMicroscope" not in nwbfile.devices:
            create_microscope_device(nwbfile, self.xml_metadata)

        # Use provided metadata or get default
        if metadata is None:
            metadata = self.get_metadata()

        # Add to NWB file using parent method
        super().add_to_nwbfile(nwbfile=nwbfile, metadata=metadata)


def extract_date_from_tiff_filename(tiff_file: Path) -> datetime:
    """
    Extract date from TIFF filename using standard patterns.

    Example filenames:
    - "20_ZSeries-20170707_Cell2_prox12-001_Cycle00001_Ch1_#.ome_Z026.tif"
    - "20_ZSeries-20190411_Cell1_dist1-001_Cycle00001_Ch1_#.ome_Z01.tif"

    Parameters
    ----------
    tiff_file : Path
        Path to TIFF file

    Returns
    -------
    datetime
        datetime object with extracted date
    """
    filename = tiff_file.name

    # Look for date pattern YYYYMMDD or MMDDYYYY in filename
    date_match = re.search(r"(\d{8})", filename)
    if date_match:
        date_str = date_match.group(1)

        # Try YYYYMMDD format first
        year = int(date_str[:4])
        month = int(date_str[4:6])
        day = int(date_str[6:8])

        # If month is invalid, try MMDDYYYY format
        if month < 1 or month > 12:
            month = int(date_str[:2])
            day = int(date_str[2:4])
            year = int(date_str[4:8])

            # Validate the corrected date
            if month < 1 or month > 12 or day < 1 or day > 31:
                raise ValueError(f"Could not parse valid date from: {date_str} in filename: {filename}")
    else:
        raise ValueError(f"Could not extract date from filename: {filename}")

    # Illinois is in Central Time Zone
    central_tz = ZoneInfo("America/Chicago")

    return datetime(year, month, day, 0, 0, 0, tzinfo=central_tz)


def parse_session_info(session_folder: Path) -> Dict[str, Any]:
    """
    Parse session information, preferring TIFF filename dates over folder names.

    Parameters
    ----------
    session_folder : Path
        Path to session folder containing subfolders with TIFF files

    Returns
    -------
    Dict[str, Any]
        Dictionary containing session information including start time, animal ID, and session ID
    """
    # Look for TIFF files in subfolders
    session_start_time = None
    for subfolder in session_folder.iterdir():
        if subfolder.is_dir():
            tiff_files = [f for f in subfolder.iterdir() if f.suffix.lower() == ".tif"]
            if tiff_files:
                session_start_time = extract_date_from_tiff_filename(tiff_files[0])
                break

    if not session_start_time:
        raise ValueError(f"Could not extract date from TIFF files in folder: {session_folder.name}")

    # Extract animal ID from folder name
    folder_name = session_folder.name
    if "2019" in folder_name:
        animal_id = folder_name[8:] if len(folder_name) > 8 else "unknown"
    else:
        animal_id = folder_name[4:] if len(folder_name) > 4 else "unknown"

    return {
        "session_start_time": session_start_time,
        "animal_id": animal_id,
        "date_str": f"{session_start_time.year}-{session_start_time.month:02d}-{session_start_time.day:02d}",
    }


def parse_container_info(subfolder_name: str) -> Dict[str, str]:
    """
    Parse container information from subfolder names.

    Examples:
    - "Decon_20190411_Cell1_dist1" -> Cell 1, distal dendrite 1
    - "Decon_20170706_Cell1_prox12" -> Cell 1, proximal dendrites 1&2

    Parameters
    ----------
    subfolder_name : str
        Name of subfolder containing image stack

    Returns
    -------
    Dict[str, str]
        Dictionary containing container name and description
    """
    # Extract cell number and location - case insensitive
    subfolder_name_lower = subfolder_name.lower()
    cell_match = re.search(r"cell(\d+)", subfolder_name_lower)
    cell_number = cell_match.group(1) if cell_match else "unknown"

    # Extract location (proximal/distal/medium) and dendrite identifier
    # Handle various patterns like dist1, dist12, dist34, distA, distC, proxBC, etc.
    if "dist" in subfolder_name_lower:
        location_match = re.search(r"dist([0-9a-z]+)", subfolder_name_lower)
        location = "Distal"
        dendrite_num = location_match.group(1) if location_match else "1"
    elif "prox" in subfolder_name_lower:
        location_match = re.search(r"prox([0-9a-z]+)", subfolder_name_lower)
        location = "Proximal"
        dendrite_num = location_match.group(1) if location_match else "1"
    elif "medium" in subfolder_name_lower:
        location_match = re.search(r"medium([0-9a-z]+)", subfolder_name_lower)
        location = "Medium"
        dendrite_num = location_match.group(1) if location_match else "1"
    elif "Prox" in subfolder_name:
        location_match = re.search(r"Prox_([0-9a-z]+)", subfolder_name)
        location = "Proximal"
        dendrite_num = location_match.group(1) if location_match else "1"
    else:
        raise ValueError(f"Could not determine location from subfolder name: {subfolder_name}")

    # Use your suggested naming scheme: Images{Proximal|Distal}{dendrite_nums}{?Real}{MeasurementNumber}
    container_name = f"Images{location}{dendrite_num}"

    # Handle "real" variants
    if "real" in subfolder_name_lower:
        container_name += "Real"

    # For Figure 7, extract measurement number from patterns like "prox2-001", "prox2-002"
    measurement_match = re.search(r"-(\d{3})$", subfolder_name)
    if measurement_match:
        measurement_num = measurement_match.group(1)
        container_name += measurement_num

    # Add distance information based on location
    if location == "Proximal":
        distance_info = "proximal (~40 μm from soma)"
    elif location == "Medium":
        distance_info = "medium (~60 μm from soma)"
    else:  # Distal
        distance_info = "distal (>80 μm from soma)"

    description = (
        f"Image stack of {location.lower()} dendrite {dendrite_num} from SPN cell {cell_number} "
        f"for spine density analysis. Location: {distance_info}. "
        f"Acquired with 0.15 μm pixels, 0.3 μm z-steps using two-photon microscopy. "
        f"Images deconvolved in AutoQuant X3.0.4 and analyzed using NeuronStudio."
    )

    return {
        "container_name": container_name,
        "description": description,
    }


def convert_spine_density_session_to_nwbfile(
    session_folder_path: Path,
    condition: str,
    figure_config: Dict[str, Any],
    session_id_parameters: Dict[str, Any],
    verbose: bool = False,
) -> NWBFile:
    """
    Convert a single session of spine density data to NWB format.

    This shared function handles the common conversion logic for all spine density
    experiments across different figures, with figure-specific variations handled through
    the configuration dictionary.

    Parameters
    ----------
    session_folder_path : Path
        Path to the session folder containing image stack subfolders
    condition : str
        Experimental condition (e.g., "LID off-state dSPN", "control iSPN", etc.)
    figure_config : Dict[str, Any]
        Figure-specific configuration containing:
        - metadata_key: str, key for session_specific_metadata.yaml
    session_id_parameters : Dict[str, Any]
        Session ID parameters using revised schema (excluding timestamp which is determined here):
        - fig: str, figure number (e.g., "F2", "F4")
        - meas_comp: str, measurement + compartment (e.g., "SpineDens")
        - cell_type: str, cell type (e.g., "dSPN", "iSPN")
        - state: str, experimental state (e.g., "OffState", "OnState")
        - pharm: str, pharmacology condition (e.g., "none", "M1RaThp")
        - geno: str, genotype (e.g., "WT", "CDGIKO")
    verbose : bool, default=False
        Enable verbose output showing detailed processing information

    Returns
    -------
    NWBFile
        NWB file with the converted data

    Notes
    -----
    This function implements the shared spine density conversion logic while allowing
    figure-specific customization through the configuration parameters.
    """
    # Extract configuration
    metadata_key = figure_config["metadata_key"]

    # Parse session information
    session_info = parse_session_info(session_folder_path)

    # Create canonical session ID with explicit parameters
    timestamp = session_info["session_start_time"].strftime("%Y%m%d%H%M%S")

    # Use session ID parameters passed by the calling script (revised schema)
    session_id = generate_canonical_session_id(
        fig=session_id_parameters["fig"],
        meas_comp=session_id_parameters["meas_comp"],
        cell_type=session_id_parameters["cell_type"],
        state=session_id_parameters["state"],
        pharm=session_id_parameters["pharm"],
        geno=session_id_parameters["geno"],
        timestamp=timestamp,
    )

    # Add session_id to session_info
    session_info["session_id"] = session_id

    # Load general and session-specific metadata from YAML files
    general_metadata_path = Path(__file__).parent.parent.parent / "general_metadata.yaml"
    general_metadata = load_dict_from_file(general_metadata_path)

    session_metadata_path = Path(__file__).parent.parent.parent / "session_specific_metadata.yaml"
    session_metadata_template = load_dict_from_file(session_metadata_path)
    script_template = session_metadata_template[metadata_key]

    # Get condition formatting for description
    condition_human_readable = format_condition[condition]["human_readable"]

    # Build pharmacology description using centralized mapping
    pharmacology_base = general_metadata["NWBFile"]["pharmacology"]
    pharm_token = session_id_parameters["pharm"]

    # Add pharmacology-specific text if applicable
    if pharm_token != "none" and pharm_token in PHARMACOLOGY_ADDITIONS:
        pharmacology_text = pharmacology_base + " " + PHARMACOLOGY_ADDITIONS[pharm_token]
    else:
        pharmacology_text = pharmacology_base

    # Create session-specific metadata from template with runtime substitutions
    session_specific_metadata = {
        "NWBFile": {
            "session_description": script_template["NWBFile"]["session_description"].format(
                condition=condition_human_readable,
                animal_id=session_info["animal_id"],
                date_str=session_info["date_str"],
            ),
            "identifier": str(uuid.uuid4()),
            "session_start_time": session_info["session_start_time"],
            "session_id": session_info["session_id"],
            "pharmacology": pharmacology_text,
            "keywords": script_template["NWBFile"]["keywords"],
        },
        "Subject": {
            "subject_id": f"SubjectRecordedAt{timestamp}",
            "description": script_template["Subject"]["description"].format(
                session_id=session_info["session_id"],
                date_str=session_info["date_str"],
                animal_id=session_info["animal_id"],
            ),
            "genotype": script_template["Subject"]["genotype"],
        },
    }

    # Deep merge with general metadata
    merged_metadata = dict_deep_update(general_metadata, session_specific_metadata)

    # Create NWB file using neuroconv helper function
    nwbfile = make_nwbfile_from_metadata(merged_metadata)

    # Process each image stack using TiffImageStackInterface
    subfolders = [f for f in session_folder_path.iterdir() if f.is_dir()]

    for subfolder in subfolders:
        # Parse container information
        container_info = parse_container_info(subfolder.name)

        # Create TiffImageStackInterface (handles file filtering, XML parsing, and device creation)
        interface = TiffImageStackInterface(subfolder=subfolder, container_info=container_info, verbose=verbose)

        # Add to NWB file (automatically creates microscope device and metadata)
        interface.add_to_nwbfile(nwbfile=nwbfile)

    if verbose:
        print(f"Created NWB file with {len(subfolders)} image stacks and session ID: {session_id}")

    return nwbfile
