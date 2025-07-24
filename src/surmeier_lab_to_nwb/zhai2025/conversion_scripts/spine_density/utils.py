"""
Utility functions for spine density conversion scripts - Zhai et al. 2025
========================================================================

This module contains shared utility functions used across all spine density
conversion scripts for parsing XML metadata files and other common operations.
"""

import re
from pathlib import Path
from typing import Any, Dict, List, Optional

import xmltodict
from pynwb import NWBFile
from pynwb.device import Device


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
