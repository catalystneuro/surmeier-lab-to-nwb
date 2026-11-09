"""
Utility functions for spine density conversion scripts - Zhai et al. 2025
========================================================================

This module contains shared utility functions used across all spine density
conversion scripts for parsing XML metadata files and other common operations.
"""

import re
import uuid
from collections import Counter
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
    surgery_and_pharmacology,
)
from surmeier_lab_to_nwb.zhai2025.devices import get_or_create_bruker_ultima_model
from surmeier_lab_to_nwb.zhai2025.subjects import find_subject


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

    # The stack's description names the cell and the site; the lab's file names give the wrong cell for four cells
    description = (
        f"Z-slice {z_slice_num} at depth {z_depth_um:.1f} μm. Acquired with two-photon microscopy "
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

    # Link to the shared Bruker Ultima DeviceModel rather than using the deprecated
    # Device.manufacturer attribute (pynwb 3.1+ pattern).
    bruker_model = get_or_create_bruker_ultima_model(nwbfile)
    microscope_device = Device(name="TwoPhotonMicroscope", description=device_description, model=bruker_model)

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
        super().__init__(file_paths=tiff_files, metadata_key=container_info["container_name"])

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


# A cell folder: the day and a letter that numbers the cell. The day is MMDD (0714b), MMDDYYYY (04122019a), YYYYMMDD
# (20221004a) or, once, MMDD before the letter and the year after it (0411a2019).
_RE_CELL_FOLDER = re.compile(r"^(\d{4}|\d{8})([a-d])(\d{4})?$")
_LETTER_TO_CELL = {"a": 1, "b": 2, "c": 3, "d": 4}
# A stack folder ends with its site: prox, dist or medium, the dendrites the lab numbered (none, 1, 12, BC), and "real"
# before or after it for a second spot of the site; raw Prairie View stacks add the sequence number (-001).
_RE_STACK_SITE = re.compile(r"cell(\d+)[_-](real)?(prox|dist|medium)([0-9a-z]*?)(_?real)?(?:-(\d{3}))?$")
# Cell folders whose stack names give another cell. The folder letter agrees with the cells the lab's Data Connections
# sheet lists for spine density on that day, so it is used.
_STACK_NAME_CELL_MISMATCHES = {
    "LID off-state dSPN/1111b",  # stacks say Cell1, the sheet lists cells 2 and 3
    "PD dSPN/1216b",  # stacks say Cell1, the sheet lists cells 2, 3 and 4
    "LID off-state iSPN/0523b",  # stacks say Cell3, the sheet lists cells 1 and 2; one raw name says Cell2
    "control/08122020b",  # stacks say Cell3, the sheet lists cells 1 and 2
}
_LOCATION = {"prox": "Proximal", "dist": "Distal", "medium": "Medium"}
_DISTANCE = {
    "Proximal": "proximal (~40 μm from soma)",
    "Medium": "medium (~60 μm from soma)",
    "Distal": "distal (>80 μm from soma)",
}


def cell_number_of(cell_folder: Path) -> int:
    """The cell number of a cell folder, from its letter; the stack names inside must agree unless listed above."""
    match = _RE_CELL_FOLDER.match(cell_folder.name)
    if match is None:
        raise ValueError(f"Unrecognized spine density cell folder: {cell_folder}")
    cell_number = _LETTER_TO_CELL[match.group(2)]
    for stack in (p for p in cell_folder.iterdir() if p.is_dir()):
        site = _RE_STACK_SITE.search(stack.name.lower())
        if site is None:
            raise ValueError(f"Unrecognized spine density stack folder: {stack}")
        if int(site.group(1)) != cell_number and f"{cell_folder.parent.name}/{cell_folder.name}" not in (
            _STACK_NAME_CELL_MISMATCHES
        ):
            raise ValueError(f"{stack} names cell {site.group(1)} but is in the folder of cell {cell_number}")
    return cell_number


def cell_date(cell_folder: Path) -> str:
    """The recording day (YYYYMMDD) of a cell folder.

    The folder name gives the month and day, and the year when it has one. Otherwise the year is the one most of the
    cell's stack and image file names give for that month and day: a few names carry typos (20150523 among 20160523,
    20160602 for the 0502 folders, 20140426 for 04262019a).
    """
    match = _RE_CELL_FOLDER.match(cell_folder.name)
    if match is None:
        raise ValueError(f"Unrecognized spine density cell folder: {cell_folder}")
    digits, _, year_after = match.groups()
    if year_after:
        return year_after + digits
    if len(digits) == 8:
        return digits if digits.startswith("20") else digits[4:] + digits[:4]
    votes: Counter = Counter()
    for stack in (p for p in cell_folder.iterdir() if p.is_dir()):
        names = [stack.name] + [f.name for f in stack.iterdir() if f.suffix.lower() == ".tif"]
        for found in (date for name in names for date in re.findall(r"(?<!\d)(\d{8})(?!\d)", name)):
            # Names give YYYYMMDD or, in a few cells, MMDDYYYY
            date = found if found[:2] in ("19", "20", "21") else found[4:] + found[:4]
            if date[4:] == digits:
                votes[date] += 1
    ranked = votes.most_common(2)
    if not ranked or (len(ranked) == 2 and ranked[0][1] == ranked[1][1]):
        raise ValueError(f"No single year for {cell_folder} in its file names: {dict(votes)}")
    return ranked[0][0]


def group_cells_by_mouse_day(condition_folder: Path) -> Dict[str, List[Path]]:
    """Group the cell folders of one condition by recording day (YYYYMMDD); the lab recorded one mouse per day."""
    cells_per_day: Dict[str, List[Path]] = {}
    for cell_folder in sorted(p for p in condition_folder.iterdir() if p.is_dir()):
        cell_number_of(cell_folder)
        cells_per_day.setdefault(cell_date(cell_folder), []).append(cell_folder)
    return cells_per_day


def parse_container_info(subfolder_name: str, cell_number: int) -> Dict[str, str]:
    """
    Name and describe one image stack, e.g. "Decon_20170706_Cell1_prox12" -> ImagesCell1Proximal12.

    The name is Images, the cell, the location and the dendrites as the lab numbered them; a second spot of a site
    (the lab's "real") is told apart by the caller.

    Parameters
    ----------
    subfolder_name : str
        Name of subfolder containing image stack
    cell_number : int
        The cell number, from the cell folder's letter

    Returns
    -------
    Dict[str, str]
        Dictionary containing container name, site, whether the lab marked it real, and description
    """
    site = _RE_STACK_SITE.search(subfolder_name.lower())
    if site is None:
        raise ValueError(f"Could not determine location from subfolder name: {subfolder_name}")
    _, real_before, location_token, dendrites, real_after, _ = site.groups()
    location = _LOCATION[location_token]
    dendrites = dendrites.upper() if dendrites.isalpha() else dendrites
    container_name = f"ImagesCell{cell_number}{location}{dendrites}"

    description = (
        f"Image stack of {location.lower()} dendrite {dendrites or '(not numbered)'} from SPN cell {cell_number} "
        f"for spine density analysis. Location: {_DISTANCE[location]}. "
        f"Acquired with 0.15 μm pixels, 0.3 μm z-steps using two-photon microscopy. "
        f"Images deconvolved in AutoQuant X3.0.4 and analyzed using NeuronStudio."
    )

    return {
        "container_name": container_name,
        "real": bool(real_before or real_after),
        "description": description,
    }


def convert_spine_density_mouse_day_to_nwbfile(
    cell_folders: List[Path],
    condition: str,
    figure_config: Dict[str, Any],
    session_id_parameters: Dict[str, Any],
    verbose: bool = False,
) -> NWBFile:
    """
    Convert the spine density stacks of all cells of one mouse on one day to one NWB file.

    This shared function handles the common conversion logic for all spine density
    experiments across different figures, with figure-specific variations handled through
    the configuration dictionary.

    Parameters
    ----------
    cell_folders : List[Path]
        The cell folders of one recording day and condition, each holding image stack subfolders
    condition : str
        Experimental condition (e.g., "LID off-state dSPN", "control iSPN", etc.)
    figure_config : Dict[str, Any]
        Figure-specific configuration containing:
        - metadata_key: str, key for session_specific_metadata.yaml
        - animal: str, the mouse line as the Data Connections sheet writes it (e.g., "tD-Tomato")
        - category: str, the group as the Data Connections sheet writes it (e.g., "CONT", "LID OFF")
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

    # The deconvolved stacks carry no acquisition time, so the session starts at midnight of the recording day
    dates = {cell_date(cell_folder) for cell_folder in cell_folders}
    if len(dates) != 1:
        raise ValueError(f"The cell folders {cell_folders} are from different days")
    (date_token,) = dates
    session_start_time = datetime.strptime(date_token, "%Y%m%d").replace(tzinfo=ZoneInfo("America/Chicago"))
    date_str = session_start_time.strftime("%Y-%m-%d")
    cells = {cell_number_of(cell_folder): cell_folder for cell_folder in cell_folders}
    if len(cells) != len(cell_folders):
        raise ValueError(f"Two cell folders of one day have the same cell number: {cell_folders}")

    session_id = generate_canonical_session_id(
        cell_type=session_id_parameters["cell_type"],
        state=session_id_parameters["state"],
        pharm=session_id_parameters["pharm"],
        geno=session_id_parameters["geno"],
        date=date_token,
    )

    # Load general and session-specific metadata from YAML files
    general_metadata_path = Path(__file__).parent.parent.parent / "general_metadata.yaml"
    general_metadata = load_dict_from_file(general_metadata_path)

    session_metadata_path = Path(__file__).parent.parent.parent / "session_specific_metadata.yaml"
    session_metadata_template = load_dict_from_file(session_metadata_path)
    script_template = session_metadata_template[metadata_key]

    # Get condition formatting for description
    condition_human_readable = format_condition[condition]["human_readable"]
    cell_numbers = sorted(cells)
    cells_text = ", ".join(map(str, cell_numbers[:-1])) + " and " if len(cell_numbers) > 1 else ""
    cells_text += str(cell_numbers[-1])

    subject = find_subject(date_token, figure_config["animal"], figure_config["category"])

    # Surgery and pharmacology follow the group: controls had no levodopa, unlesioned controls no lesion either
    state = session_id_parameters["state"]
    surgery, pharmacology_base = surgery_and_pharmacology(state, general_metadata)
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
                cells=cells_text,
                date_str=date_str,
            ),
            "experiment_description": script_template["NWBFile"]["experiment_description"],
            "identifier": str(uuid.uuid4()),
            "session_start_time": session_start_time,
            "session_id": session_id,
            "surgery": surgery,
            "pharmacology": pharmacology_text,
            "keywords": script_template["NWBFile"]["keywords"],
        },
        "Subject": {
            # The subjects table's ID for this mouse, shared with its patch-clamp file
            "subject_id": subject.subject_id,
            # The control groups have their own description; the default one is for lesioned, levodopa-treated mice
            "description": script_template["Subject"]
            .get("description_by_state", {})
            .get(state, script_template["Subject"]["description"])
            .format(session_id=session_id, date_str=date_str),
            "strain": subject.strain,
            "genotype": subject.genotype,
        },
    }

    # Deep merge with general metadata
    merged_metadata = dict_deep_update(general_metadata, session_specific_metadata)
    if subject.strain is None:  # the paper gives no background for this line
        del merged_metadata["Subject"]["strain"]

    # Create NWB file using neuroconv helper function
    nwbfile = make_nwbfile_from_metadata(merged_metadata)

    # Name every stack; a site the lab imaged at two spots, the second marked "real", gets Spot1 and Spot2, and the
    # one cell with two stacks of a site (Figure 7, 2016-08-16 cell 1, prox2-001 and prox2-002) gets Stack1 and Stack2.
    stacks = []
    for cell_number in cell_numbers:
        for subfolder in sorted(p for p in cells[cell_number].iterdir() if p.is_dir()):
            stacks.append((subfolder, parse_container_info(subfolder.name, cell_number)))
    base_names = [info["container_name"] for _, info in stacks]
    for (subfolder, container_info), base_name in zip(stacks, base_names):
        same_site = [info for (_, info), name in zip(stacks, base_names) if name == base_name]
        if len(same_site) == 2 and {info["real"] for info in same_site} == {False, True}:
            spot = 2 if container_info["real"] else 1
            container_info["container_name"] = f"{base_name}Spot{spot}"
            container_info[
                "description"
            ] += f" The lab imaged this site at two spots of the dendrite; this is spot {spot}" + (
                ", which the lab's files mark 'real'." if spot == 2 else "."
            )
        elif len(same_site) > 1:
            number = [s for (s, _), name in zip(stacks, base_names) if name == base_name].index(subfolder) + 1
            container_info["container_name"] = f"{base_name}Stack{number}"
        elif container_info["real"]:
            container_info["description"] += " The lab's files mark this stack 'real'."

        # Create TiffImageStackInterface (handles file filtering, XML parsing, and device creation)
        interface = TiffImageStackInterface(subfolder=subfolder, container_info=container_info, verbose=verbose)

        # Add to NWB file (automatically creates microscope device and metadata)
        interface.add_to_nwbfile(nwbfile=nwbfile)

    if verbose:
        print(f"Created NWB file with {len(stacks)} image stacks and session ID: {session_id}")

    return nwbfile
