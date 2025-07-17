import re
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, Optional
from zoneinfo import ZoneInfo

import xmltodict
from neuroconv.datainterfaces import ImageInterface
from neuroconv.tools import configure_and_write_nwbfile
from neuroconv.utils import dict_deep_update, load_dict_from_file
from pynwb import NWBFile
from pynwb.device import Device


def parse_xml_metadata(xml_file: Path) -> Dict[str, Any]:
    """
    Parse XML metadata file to extract acquisition parameters.

    Args:
        xml_file: Path to XML metadata file

    Returns:
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

    # Extract metadata using xmltodict approach
    return _extract_metadata_from_dict(xml_dict)


def _extract_metadata_from_dict(xml_dict: Dict[str, Any]) -> Dict[str, Any]:
    """Extract metadata from xmltodict parsed dictionary."""
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

    Args:
        subfolder: Path to image stack subfolder

    Returns:
        Path to XML metadata file or None if not found
    """
    # Each folder contains exactly one XML file - find it by suffix
    xml_files = [f for f in subfolder.iterdir() if f.suffix == ".xml"]

    if xml_files:
        return xml_files[0]  # Return the single XML file

    return None


def extract_date_from_tiff_filename(tiff_file: Path) -> datetime:
    """
    Extract date from TIFF filename.

    Example filenames:
    - "20_ZSeries-20170707_Cell2_prox12-001_Cycle00001_Ch1_#.ome_Z026.tif"
    - "20_ZSeries-20190411_Cell1_dist1-001_Cycle00001_Ch1_#.ome_Z01.tif"

    Args:
        tiff_file: Path to TIFF file

    Returns:
        datetime object with extracted date
    """
    filename = tiff_file.name

    # Look for date pattern YYYYMMDD in filename
    date_match = re.search(r"(\d{8})", filename)
    if date_match:
        date_str = date_match.group(1)
        year = int(date_str[:4])
        month = int(date_str[4:6])
        day = int(date_str[6:8])
    else:
        raise ValueError(f"Could not extract date from filename: {filename}")

    # Illinois is in Central Time Zone
    central_tz = ZoneInfo("America/Chicago")

    return datetime(year, month, day, 0, 0, 0, tzinfo=central_tz)


def parse_session_info(session_folder: Path) -> Dict[str, Any]:
    """
        Parse session information, preferring TIFF filename dates over folder names.
    x
    """

    # Look for TIFF files in subfolders
    for subfolder in session_folder.iterdir():
        if subfolder.is_dir():
            tiff_files = [f for f in subfolder.iterdir() if f.suffix.lower() == ".tif"]
            tiff_date = extract_date_from_tiff_filename(tiff_files[0])
            if tiff_date:
                session_start_time = tiff_date
            else:
                raise ValueError(f"Could not extract date from TIFF file in folder: {subfolder.name}")

    # Extract animal ID from folder name
    # TODO: this is a guess in the moment, need to confirm
    folder_name = session_folder.name
    if "2019" in folder_name:
        animal_id = folder_name[8:] if len(folder_name) > 8 else "unknown"
    else:
        animal_id = folder_name[4:] if len(folder_name) > 4 else "unknown"

    return {
        "session_start_time": session_start_time,
        "animal_id": animal_id,
        "date_str": f"{session_start_time.year}-{session_start_time.month:02d}-{session_start_time.day:02d}",
        "session_id": f"{session_start_time.year}{session_start_time.month:02d}{session_start_time.day:02d}_{animal_id}",
    }


def parse_container_info(subfolder_name: str) -> Dict[str, str]:
    """
    Parse container information from subfolder names.

    Examples:
    - "Decon_20190411_Cell1_dist1" -> Cell 1, distal dendrite 1
    - "Decon_20170706_Cell1_prox12" -> Cell 1, proximal dendrites 1&2

    """
    # Extract cell number and location
    cell_match = re.search(r"Cell(\d+)", subfolder_name)
    cell_number = cell_match.group(1) if cell_match else "unknown"

    # Extract location (proximal/distal) and dendrite number
    if "dist" in subfolder_name:
        location_match = re.search(r"dist(\d+)", subfolder_name)
        location = "Distal"
        dendrite_num = location_match.group(1) if location_match else "1"
    elif "prox" in subfolder_name:
        location_match = re.search(r"prox(\d+)", subfolder_name)
        location = "Proximal"
        dendrite_num = location_match.group(1) if location_match else "1"
    elif "Prox" in subfolder_name:
        location_match = re.search(r"Prox_(\d+)", subfolder_name)
        location = "Proximal"
        dendrite_num = location_match.group(1) if location_match else "1"
    else:
        raise ValueError(f"Could not determine location from subfolder name: {subfolder_name}")

    container_name = f"ImagesCell{cell_number}{location}Dendrite{dendrite_num}"

    # TODO: figure out what real means here
    # Most likely the experiment was done twice and this is "correct"
    if "real" in subfolder_name:
        container_name += "Real"

    # Add distance information based on location
    distance_info = "proximal (~40 μm from soma)" if location == "Proximal" else "distal (>80 μm from soma)"

    description = (
        f"Image stack of {location.lower()} dendrite {dendrite_num} from dSPN cell {cell_number} "
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

    Args:
        nwbfile: NWB file to add device to
        xml_metadata: Metadata extracted from XML file

    Returns:
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
    Convert spine density data to NWB format.

    Parameters
    ----------
    session_folder_path : Path
        Path to the top level folders in the conditions for spiny density figure 2
    condition : str
        The experimental condition (e.g., 'control dSPN', 'LID on-state dSPN')
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
    metadata_file_path = Path(__file__).parent / "metadata.yaml"
    paper_metadata = load_dict_from_file(metadata_file_path)

    if verbose:
        print(f"Loaded paper metadata from: {metadata_file_path}")
        print(f"Experiment description: {paper_metadata['NWBFile']['experiment_description'][:100]}...")

    # Create session-specific metadata
    session_specific_metadata = {
        "NWBFile": {
            "session_description": (
                f"Dendritic spine density assessment in direct pathway spiny projection neurons (dSPNs) "
                f"for condition {condition}. Two-photon laser scanning microscopy was used to acquire "
                f"Z-stack images of dendritic segments at two locations: proximal (~40 μm from soma) "
                f"and distal (>80 μm from soma). Acquisition parameters: 0.15 μm pixels, 0.3 μm z-steps, "
                f"60x objective (NA=1.0), optical zoom 5.2x, 10 μs dwell time. Images were deconvolved "
                f"using AutoQuant X3.0.4 (MediaCybernetics) and semi-automated spine counting was performed "
                f"using 3D reconstructions in NeuronStudio (CNIC, Mount Sinai). On average, 2 proximal "
                f"and 2 distal dendrites were imaged and analyzed per neuron."
            ),
            "identifier": f"zhai2025_fig2_spine_density_{session_info['session_id']}_{condition.replace(' ', '_')}",
            "session_start_time": session_info["session_start_time"],
            "session_id": session_info["session_id"],
            "keywords": ["spine density", "dendritic spines", "two-photon microscopy"],
        }
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

    # Create microscope device using metadata from first available XML file
    microscope_device = None
    subfolders = [f for f in session_folder_path.iterdir() if f.is_dir()]

    if verbose:
        print(f"Found {len(subfolders)} image stacks in session")

    # Process each image stack
    for subfolder in subfolders:
        if verbose:
            print(f"  Processing: {subfolder.name}")

        # Parse container information
        container_info = parse_container_info(subfolder.name)

        if verbose:
            print(f"    Container name: {container_info['container_name']}")

        # Get all TIFF files in the folder
        tiff_files = sorted([f for f in subfolder.iterdir() if f.suffix.lower() == ".tif"])
        assert len(tiff_files) > 0, f"No TIFF files found in {subfolder.name}"

        if verbose:
            print(f"    Found {len(tiff_files)} TIFF files in {subfolder.name}")

        # Find and parse XML metadata file
        xml_file = find_xml_metadata_file(subfolder)
        if not xml_file:
            raise FileNotFoundError(f"No XML metadata file found in {subfolder.name}")

        if verbose:
            print(f"    Found XML metadata file: {xml_file.name}")

        xml_metadata = parse_xml_metadata(xml_file)
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

        interface = ImageInterface(file_paths=tiff_files)

        # Create comprehensive description with all acquisition parameters
        detailed_description = (
            f"{container_info['description']} "
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

        metadata = interface.get_metadata()
        metadata["Images"]["description"] = detailed_description

        if verbose:
            print(f"    Adding interface to NWB file with container: {container_info['container_name']}")

        interface.add_to_nwbfile(nwbfile=nwbfile, metadata=metadata, container_name=container_info["container_name"])

        if verbose:
            print(f"    Successfully added image stack: {container_info['container_name']}")

    if verbose:
        print(f"Conversion completed for session: {session_info['session_id']}")

    return nwbfile


if __name__ == "__main__":
    import logging

    from tqdm import tqdm

    # Control verbose output
    verbose = False  # Set to True for detailed output

    logging.getLogger("tifffile").setLevel(logging.ERROR)

    # Define the base path to the data
    base_path = Path("/home/heberto/development/surmeier-lab-to-nwb/link_to_raw_data/Figure 2_SF1A/Spine density/")

    # Create nwb_files directory at root level
    root_dir = Path(__file__).parent.parent.parent.parent  # Go up to repo root
    nwb_files_dir = root_dir / "nwb_files" / "figure_2_spine_density"
    nwb_files_dir.mkdir(parents=True, exist_ok=True)

    conditions = ["PD dSPN", "LID off-state dSPN", "control dSPN", "LID on-state dSPN"]
    for condition in conditions:
        condition_path = base_path / condition
        assert condition_path.exists(), f"Base path does not exist: {base_path}"

        print(f"Processing spine density data for: {condition=}")

        # Get all session folders
        session_folders = [f for f in condition_path.iterdir() if f.is_dir()]
        session_folders.sort()

        print(f"Found {len(session_folders)} session folders")

        # Use tqdm for progress bar when verbose is disabled
        session_iterator = (
            tqdm(session_folders, desc=f"Processing {condition}", disable=verbose) if not verbose else session_folders
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
            nwbfile_path = nwb_files_dir / f"figure2_spine_density_{condition_safe}_{session_folder_path.name}.nwb"

            # Write NWB file
            configure_and_write_nwbfile(nwbfile, nwbfile_path=nwbfile_path)
            print(f"Successfully saved: {nwbfile_path.name}")
