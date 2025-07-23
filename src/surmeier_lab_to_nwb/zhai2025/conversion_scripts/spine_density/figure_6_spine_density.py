"""
Figure 6 Spine Density Conversion Script - Zhai et al. 2025
========================================================

This script converts two-photon spine density imaging data from Figure 6 of
Zhai et al. 2025 into NWB (Neurodata Without Borders) format. The data examines
dendritic spine density changes in a mouse model of Parkinson's disease
and levodopa-induced dyskinesia.

For detailed experimental context, protocols, data structure, and analysis methods,
see: /src/surmeier_lab_to_nwb/zhai2025/conversion_notes_folder/figure_6_conversion_notes.md
"""

import re
import uuid
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional
from zoneinfo import ZoneInfo

import xmltodict
from neuroconv.datainterfaces import ImageInterface
from neuroconv.tools import configure_and_write_nwbfile
from neuroconv.utils import dict_deep_update, load_dict_from_file
from pynwb import NWBFile
from pynwb.device import Device
from pynwb.file import Subject


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

    if verbose:
        print(f"    Parsing XML metadata from: {xml_file.name}")

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

    if verbose:
        print(f"    Extracted metadata keys: {list(metadata.keys())}")

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

    return xml_files[0]  # Return the single XML file


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

    if verbose:
        print(f"    Creating metadata for {len(tiff_files_with_z)} individual images")

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

        if verbose:
            print(
                f"      Added metadata for: {z_info['image_name']} (Z={z_info['z_slice_number']}, "
                f"depth={z_info['z_depth_um']:.1f}μm)"
            )

    return images_metadata


def extract_date_from_tiff_filename(tiff_file: Path) -> datetime:
    """
    Extract date from TIFF filename.

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
    for subfolder in session_folder.iterdir():
        if subfolder.is_dir():
            tiff_files = [f for f in subfolder.iterdir() if f.suffix.lower() == ".tif"]
            if tiff_files:
                tiff_date = extract_date_from_tiff_filename(tiff_files[0])
                if tiff_date:
                    session_start_time = tiff_date
                    break
            else:
                raise ValueError(f"Could not extract date from TIFF file in folder: {subfolder.name}")

    # Extract animal ID from folder name
    # Figure 6 uses iSPNs, adjust animal ID extraction for different naming pattern
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


def parse_container_info(subfolder_name: str, session_id: str) -> Dict[str, str]:
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

    # Extract location (proximal/distal) and dendrite identifier - case insensitive
    if "dist" in subfolder_name_lower:
        location_match = re.search(r"dist([a-z0-9]+)", subfolder_name_lower)
        location = "Distal"
        dendrite_num = location_match.group(1) if location_match else "1"
    elif "prox" in subfolder_name_lower:
        location_match = re.search(r"prox([a-z0-9]+)", subfolder_name_lower)
        location = "Proximal"
        dendrite_num = location_match.group(1) if location_match else "1"
    else:
        raise ValueError(f"Could not determine location from subfolder name: {subfolder_name}")

    container_name = f"Images{location}{dendrite_num}"

    # Handle "real" designation if present
    if "real" in subfolder_name:
        container_name += "Real"

    # Add distance information based on location
    distance_info = "proximal (~40 μm from soma)" if location == "Proximal" else "distal (>80 μm from soma)"

    # Note: Figure 6 data is from iSPNs (indirect pathway), not dSPNs
    description = (
        f"Image stack of {location.lower()} dendrite {dendrite_num} from iSPN cell {cell_number} "
        f"for spine density analysis. Location: {distance_info}. "
        f"Acquired with 0.15 μm pixels, 0.3 μm z-steps using two-photon microscopy. "
        f"Images deconvolved in AutoQuant X3.0.4 and analyzed using NeuronStudio."
    )

    return {
        "container_name": container_name,
        "description": description,
    }


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


def convert_data_to_nwb(session_folder_path: Path, condition: str, verbose: bool = False) -> NWBFile:
    """
    Convert spine density data to NWB format for Figure 6.

    Parameters
    ----------
    session_folder_path : Path
        Path to the top level folders in the conditions for spine density figure 6
    condition : str
        The experimental condition (e.g., 'control', 'M1R antagonist')
    verbose : bool, default=False
        Enable verbose output showing detailed processing information

    Returns
    -------
    NWBFile
        The populated NWB file object
    """
    # Parse session information
    session_info = parse_session_info(session_folder_path)
    if verbose:
        print(f"Session date: {session_info['date_str']}, Animal: {session_info['animal_id']}")

    # Load metadata from YAML file
    metadata_file_path = Path(__file__).parent.parent.parent / "metadata.yaml"
    paper_metadata = load_dict_from_file(metadata_file_path)

    if verbose:
        print(f"Loaded paper metadata from: {metadata_file_path}")
        print(f"Experiment description: {paper_metadata['NWBFile']['experiment_description'][:100]}...")

    # Create BIDS-style base session ID with detailed timestamp when available
    session_start_time = session_info["session_start_time"]
    if hasattr(session_start_time, "hour"):
        timestamp = session_start_time.strftime("%Y%m%d_%H%M%S")
    else:
        timestamp = session_start_time.strftime("%Y%m%d")

    base_session_id = f"figure6_SpineDensity_{condition.replace(' ', '_').replace('-', '_')}_{timestamp}"
    script_specific_id = f"Sub{session_info['animal_id']}"
    session_id = f"{base_session_id}_{script_specific_id}"

    # Create session-specific metadata for Figure 6
    session_specific_metadata = {
        "NWBFile": {
            "session_description": (
                f"Dendritic spine density assessment in indirect pathway spiny projection neurons (iSPNs) "
                f"for condition {condition}. Two-photon laser scanning microscopy was used to acquire "
                f"Z-stack images of dendritic segments at two locations: proximal (~40 μm from soma) "
                f"and distal (>80 μm from soma). Acquisition parameters: 0.15 μm pixels, 0.3 μm z-steps, "
                f"60x objective (NA=1.0), optical zoom 5.2x, 10 μs dwell time. Images were deconvolved "
                f"using AutoQuant X3.0.4 (MediaCybernetics) and semi-automated spine counting was performed "
                f"using 3D reconstructions in NeuronStudio (CNIC, Mount Sinai). This data is part of "
                f"Figure 6, which investigates the role of M1R signaling in spine density."
            ),
            "identifier": str(uuid.uuid4()),
            "session_start_time": session_info["session_start_time"],
            "session_id": session_id,
            "keywords": [
                "spine density",
                "dendritic spines",
                "two-photon microscopy",
                "iSPNs",
                "M1R",
                "muscarinic receptor",
            ],
        },
        "Subject": {
            "subject_id": f"dSPN_mouse_{session_info['animal_id']}",
            "description": (
                f"Adult mouse with unilateral 6-OHDA lesion (>95% dopamine depletion) modeling Parkinson's disease. "
                f"Received dyskinesiogenic levodopa treatment for spine density analysis. "
                f"Animal {session_info['animal_id']} recorded on {session_info['date_str']}. "
                f"M1R antagonist treatment: {'THP (3 mg/kg, i.p.)' if 'antagonist' in condition else 'Saline control'}."
            ),
            "genotype": "C57BL/6J wild-type",
        },
    }

    # Merge paper metadata with session-specific metadata
    merged_metadata = dict_deep_update(paper_metadata, session_specific_metadata)

    # Create NWB file with explicit arguments
    nwbfile = NWBFile(
        session_description=merged_metadata["NWBFile"]["session_description"],
        identifier=merged_metadata["NWBFile"]["identifier"],
        session_start_time=merged_metadata["NWBFile"]["session_start_time"],
        experimenter=merged_metadata["NWBFile"]["experimenter"],
        lab=merged_metadata["NWBFile"]["lab"],
        institution=merged_metadata["NWBFile"]["institution"],
        experiment_description=merged_metadata["NWBFile"]["experiment_description"],
        session_id=merged_metadata["NWBFile"]["session_id"],
        keywords=merged_metadata["NWBFile"]["keywords"],
    )

    # Create subject using merged metadata
    subject = Subject(
        subject_id=merged_metadata["Subject"]["subject_id"],
        species=merged_metadata["Subject"]["species"],
        strain=merged_metadata["Subject"]["strain"],
        description=merged_metadata["Subject"]["description"],
        genotype=merged_metadata["Subject"]["genotype"],
        sex=merged_metadata["Subject"]["sex"],
        age=merged_metadata["Subject"]["age"],
    )
    nwbfile.subject = subject

    # Create microscope device (only once, using first available XML file)
    microscope_device = None
    subfolders = [f for f in session_folder_path.iterdir() if f.is_dir()]

    if verbose:
        print(f"Found {len(subfolders)} image stacks in session")

    # Process each image stack
    for subfolder in subfolders:
        if verbose:
            print(f"  Processing: {subfolder.name}")

        # Parse container information
        container_info = parse_container_info(subfolder.name, session_id)

        if verbose:
            print(f"    Container name: {container_info['container_name']}")

        # Get all TIFF files in the folder, excluding projection images
        all_tiff_files = [f for f in subfolder.iterdir() if f.suffix.lower() == ".tif"]
        # Filter out projection images (_xy.tif, _xz.tif, _zy.tif, project_xy.tif, *_xy.tif, etc.) and keep only Z-stack images
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
        assert len(tiff_files) > 0, f"No Z-stack TIFF files found in {subfolder.name}"

        if verbose:
            print(f"    Found {len(tiff_files)} TIFF files in {subfolder.name}")

        # Find and parse XML metadata file
        xml_file = find_xml_metadata_file(subfolder)
        if not xml_file:
            raise FileNotFoundError(f"No XML metadata file found in {subfolder.name}")

        if verbose:
            print(f"    Found XML metadata file: {xml_file.name}")

        xml_metadata = parse_xml_metadata(xml_file, verbose=verbose)
        if not xml_metadata:
            raise ValueError(f"Failed to parse XML metadata from {xml_file.name}")

        if verbose:
            print(
                f"    Extracted metadata: pixel_size={xml_metadata.get('pixel_size_x', 'N/A')} μm, "
                f"objective={xml_metadata.get('objective_lens', 'N/A')}, "
                f"NA={xml_metadata.get('numerical_aperture', 'N/A')}"
            )

        # Create microscope device (only once, using first available metadata)
        if microscope_device is None:
            microscope_device = create_microscope_device(nwbfile, xml_metadata)
            if verbose:
                print(f"    Created microscope device: {microscope_device.name}")

        # Create ImageInterface with sorted TIFF files
        interface = ImageInterface(
            file_paths=tiff_files, images_container_metadata_key=container_info["container_name"]
        )

        # Get base metadata from interface
        metadata = interface.get_metadata()

        # Create custom metadata for individual images using interface's resolved file paths
        images_metadata = create_image_metadata(
            tiff_files=interface.file_paths, container_info=container_info, xml_metadata=xml_metadata, verbose=verbose
        )

        # Update the metadata with custom image information
        metadata["Images"][container_info["container_name"]].update(images_metadata)

        if verbose:
            print(f"    Adding interface to NWB file with {len(tiff_files)} individual images")

        # Add to NWB file using the interface
        interface.add_to_nwbfile(nwbfile=nwbfile, metadata=metadata)

        if verbose:
            print(f"    Successfully added image stack: {container_info['container_name']}")

    if verbose:
        print(f"Conversion completed for session: {session_info['animal_id']}")

    return nwbfile


if __name__ == "__main__":
    import logging

    from tqdm import tqdm

    # Control verbose output
    verbose = False  # Set to True for detailed output
    stub_test = True  # Set to True to process only first 2 files per condition for testing

    logging.getLogger("tifffile").setLevel(logging.ERROR)

    # Define the base path to the data - Figure 6 spine density (two-photon method, iSPNs)
    base_path = Path("/home/heberto/development/surmeier-lab-to-nwb/link_to_raw_data/Figure 6/Spine density/")

    # Create nwb_files directory at root level
    root_dir = Path(__file__).parent.parent.parent.parent.parent.parent  # Go up to repo root
    nwb_files_dir = root_dir / "nwb_files" / "spine_density" / "figure_6"
    nwb_files_dir.mkdir(parents=True, exist_ok=True)

    # Figure 6 conditions use iSPNs (indirect pathway)
    conditions = ["control", "M1R antagonist"]

    for condition in conditions:
        condition_path = base_path / condition
        assert condition_path.exists(), f"Base path does not exist: {condition_path}"

        if verbose:
            print(f"Processing spine density data for: {condition=}")

        # Get all session folders
        session_folders = [f for f in condition_path.iterdir() if f.is_dir()]
        session_folders.sort()

        # Apply stub_test filtering if enabled
        if stub_test:
            session_folders = session_folders[:2]
            if verbose:
                print(f"stub_test enabled: processing only first {len(session_folders)} session folders")

        if verbose:
            print(f"Found {len(session_folders)} session folders")

        # Use tqdm for progress bar when verbose is disabled
        session_iterator = tqdm(
            session_folders, desc=f"Converting Figure6 SpineDensity {condition}", disable=verbose, unit=" session"
        )

        for session_folder_path in session_iterator:
            if verbose:
                print(f"\nProcessing session: {session_folder_path.name}")

            # Convert data to NWB format
            nwbfile = convert_data_to_nwb(
                session_folder_path=session_folder_path,
                condition=condition,
                verbose=verbose,
            )

            # Create output filename
            condition_safe = condition.replace(" ", "_").replace("(", "").replace(")", "")
            nwbfile_path = nwb_files_dir / f"figure6_spine_density_{condition_safe}_{session_folder_path.name}.nwb"

            # Write NWB file
            configure_and_write_nwbfile(nwbfile, nwbfile_path=nwbfile_path)
            if verbose:
                print(f"Successfully saved: {nwbfile_path.name}")
