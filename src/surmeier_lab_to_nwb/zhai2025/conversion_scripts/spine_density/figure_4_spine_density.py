"""
Figure 4 Spine Density Conversion Script - Zhai et al. 2025
========================================================

This script converts two-photon spine density imaging data from Figure 4 of
Zhai et al. 2025 into NWB (Neurodata Without Borders) format. The data examines
dendritic spine density changes in indirect pathway SPNs (iSPNs) in a mouse model
of Parkinson's disease and levodopa-induced dyskinesia.

For detailed experimental context, protocols, data structure, and analysis methods,
see: /src/surmeier_lab_to_nwb/zhai2025/conversion_notes_folder/figure_4_conversion_notes.md
"""

import re
import uuid
from datetime import datetime
from pathlib import Path
from typing import Any, Dict
from zoneinfo import ZoneInfo

from neuroconv.datainterfaces import ImageInterface
from neuroconv.tools import configure_and_write_nwbfile
from neuroconv.utils import dict_deep_update, load_dict_from_file
from pynwb import NWBFile
from pynwb.file import Subject

from surmeier_lab_to_nwb.zhai2025.conversion_scripts.spine_density.utils import (
    create_image_metadata,
    create_microscope_device,
    find_xml_metadata_file,
    parse_xml_metadata,
)


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
    # Figure 4 uses iSPNs, adjust animal ID extraction for different naming pattern
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

    # Note: Figure 4 data is from iSPNs (indirect pathway), not dSPNs
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


def convert_data_to_nwb(session_folder_path: Path, condition: str, verbose: bool = False) -> NWBFile:
    """
    Convert spine density data to NWB format for Figure 4.

    Parameters
    ----------
    session_folder_path : Path
        Path to the top level folders in the conditions for spine density figure 4
    condition : str
        The experimental condition (e.g., 'control iSPN', 'LID on-state iSPN')
    verbose : bool, default=False
        Enable verbose output showing detailed processing information

    Returns
    -------
    NWBFile
        The populated NWB file object
    """
    # Parse session information
    session_info = parse_session_info(session_folder_path)

    # Create session ID following new pattern
    clean_condition = condition.replace(" ", "_").replace("-", "_")
    base_session_id = (
        f"figure4_SpineDensity_{clean_condition}_{session_info['session_start_time'].strftime('%Y%m%d_%H%M%S')}"
    )
    script_specific_id = f"Animal{session_info['animal_id']}"
    session_id = f"{base_session_id}_{script_specific_id}"

    # Add session_id to session_info
    session_info["session_id"] = session_id

    # Load general and session-specific metadata from YAML files
    general_metadata_path = Path(__file__).parent.parent.parent / "general_metadata.yaml"
    general_metadata = load_dict_from_file(general_metadata_path)

    session_metadata_path = Path(__file__).parent.parent.parent / "session_specific_metadata.yaml"
    session_metadata_template = load_dict_from_file(session_metadata_path)
    script_template = session_metadata_template["figure_4_spine_density"]

    # Create session-specific metadata from template with runtime substitutions
    session_specific_metadata = {
        "NWBFile": {
            "session_description": script_template["NWBFile"]["session_description"].format(
                condition=condition, animal_id=session_info["animal_id"], date_str=session_info["date_str"]
            ),
            "identifier": str(uuid.uuid4()),
            "session_start_time": session_info["session_start_time"],
            "experiment_description": script_template["NWBFile"]["experiment_description"],
            "session_id": session_info["session_id"],
            "keywords": script_template["NWBFile"]["keywords"],
        },
        "Subject": {
            "subject_id": f"iSPN_mouse_{session_info['session_id']}",
            "description": script_template["Subject"]["description"].format(
                session_id=session_info["session_id"], date_str=session_info["date_str"]
            ),
            "genotype": script_template["Subject"]["genotype"],
        },
    }

    # Deep merge with general metadata
    merged_metadata = dict_deep_update(general_metadata, session_specific_metadata)

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
        surgery=merged_metadata["NWBFile"]["surgery"],
        pharmacology=merged_metadata["NWBFile"]["pharmacology"],
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

    # Create microscope device using metadata from first available XML file
    microscope_device = None
    subfolders = [f for f in session_folder_path.iterdir() if f.is_dir()]

    # Process each image stack
    for subfolder in subfolders:

        # Parse container information
        container_info = parse_container_info(subfolder.name, session_info["session_id"])

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

        # Find and parse XML metadata file
        xml_file = find_xml_metadata_file(subfolder)
        if not xml_file:
            raise FileNotFoundError(f"No XML metadata file found in {subfolder.name}")

        xml_metadata = parse_xml_metadata(xml_file, verbose=verbose)
        if not xml_metadata:
            raise ValueError(f"Failed to parse XML metadata from {xml_file.name}")

        # Create microscope device (only once, using first available metadata)
        if microscope_device is None:
            microscope_device = create_microscope_device(nwbfile, xml_metadata)

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

        # Add to NWB file using the interface
        interface.add_to_nwbfile(nwbfile=nwbfile, metadata=metadata)

    return nwbfile


if __name__ == "__main__":
    import logging

    from tqdm import tqdm

    # Control verbose output
    verbose = False  # Set to True for detailed output
    stub_test = True  # Set to True to process only first 2 files per condition for testing

    logging.getLogger("tifffile").setLevel(logging.ERROR)

    # Define the base path to the data - Figure 4 spine density (two-photon method, iSPNs)
    base_path = Path("/home/heberto/development/surmeier-lab-to-nwb/link_to_raw_data/Figure 4_SF1B_SF5/Spine density/")

    # Create nwb_files directory at root level
    root_dir = Path(__file__).parent.parent.parent.parent.parent.parent  # Go up to repo root
    nwb_files_dir = root_dir / "nwb_files" / "spine_density" / "figure_4"
    nwb_files_dir.mkdir(parents=True, exist_ok=True)

    # Figure 4 conditions use iSPNs (indirect pathway)
    conditions = ["control iSPN", "LID off-state iSPN", "LID on-state iSPN", "PD iSPN"]

    for condition in conditions:
        condition_path = base_path / condition
        assert condition_path.exists(), f"Base path does not exist: {condition_path}"

        # Get all session folders
        session_folders = [f for f in condition_path.iterdir() if f.is_dir()]
        session_folders.sort()

        # Apply stub_test filtering if enabled
        if stub_test:
            session_folders = session_folders[:2]

        # Use tqdm for progress bar when verbose is disabled
        session_iterator = tqdm(
            session_folders, desc=f"Converting Figure4 SpineDensity {condition}", disable=verbose, unit=" session"
        )

        for session_folder_path in session_iterator:

            # Convert data to NWB format
            nwbfile = convert_data_to_nwb(
                session_folder_path=session_folder_path,
                condition=condition,
                verbose=verbose,
            )

            # Create output filename
            condition_safe = condition.replace(" ", "_").replace("-", "_").replace("(", "").replace(")", "")
            nwbfile_path = nwb_files_dir / f"figure4_spine_density_{condition_safe}_{session_folder_path.name}.nwb"

            # Write NWB file
            configure_and_write_nwbfile(nwbfile, nwbfile_path=nwbfile_path)
