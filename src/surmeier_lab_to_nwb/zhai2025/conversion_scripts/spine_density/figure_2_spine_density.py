"""
Figure 2 Spine Density Conversion Script - Zhai et al. 2025
========================================================

This script converts two-photon spine density imaging data from Figure 2 of
Zhai et al. 2025 into NWB (Neurodata Without Borders) format. The data examines
dendritic spine density changes in direct pathway SPNs (dSPNs) in a mouse model
of Parkinson's disease and levodopa-induced dyskinesia.

For detailed experimental context, protocols, data structure, and analysis methods,
see: /src/surmeier_lab_to_nwb/zhai2025/conversion_notes_folder/figure_2_conversion_notes.md
"""

import re
import uuid
from datetime import datetime
from pathlib import Path
from typing import Any, Dict
from zoneinfo import ZoneInfo

from neuroconv.tools import configure_and_write_nwbfile
from neuroconv.utils import dict_deep_update, load_dict_from_file
from pynwb import NWBFile
from pynwb.file import Subject

from surmeier_lab_to_nwb.zhai2025.conversion_scripts.spine_density.utils import (
    TiffImageStackInterface,
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

    container_name = f"Images{location}{dendrite_num}"

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

    # Create session ID with ++ separators (no dashes or underscores)
    condition_to_clean = {
        "control dSPN": "ControlDSPN",
        "LID off-state dSPN": "LIDOffStateDSPN",
        "LID on-state dSPN": "LIDOnStateDSPN",
        "PD dSPN": "PDDSN",
    }

    timestamp = session_info["session_start_time"].strftime("%Y%m%d%H%M%S")
    clean_condition = condition_to_clean.get(condition, condition.replace(" ", "").replace("-", ""))
    base_session_id = f"Figure2++SpineDensity++{clean_condition}++Timestamp++{timestamp}"
    script_specific_id = f"Animal++{session_info['animal_id']}"
    session_id = f"{base_session_id}++{script_specific_id}"

    # Add session_id to session_info
    session_info["session_id"] = session_id

    # Load general and session-specific metadata from YAML files
    general_metadata_path = Path(__file__).parent.parent.parent / "general_metadata.yaml"
    general_metadata = load_dict_from_file(general_metadata_path)

    session_metadata_path = Path(__file__).parent.parent.parent / "session_specific_metadata.yaml"
    session_metadata_template = load_dict_from_file(session_metadata_path)
    script_template = session_metadata_template["figure_2_spine_density"]

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
            "subject_id": f"dSPN_mouse_{session_info['session_id']}",
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
        slices=merged_metadata["NWBFile"]["slices"],
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

    # Process each image stack using TiffImageStackInterface
    subfolders = [f for f in session_folder_path.iterdir() if f.is_dir()]

    for subfolder in subfolders:
        # Parse container information
        container_info = parse_container_info(subfolder.name)

        # Create TiffImageStackInterface (handles file filtering, XML parsing, and device creation)
        interface = TiffImageStackInterface(subfolder=subfolder, container_info=container_info, verbose=verbose)

        # Add to NWB file (automatically creates microscope device and metadata)
        interface.add_to_nwbfile(nwbfile=nwbfile)

    return nwbfile


if __name__ == "__main__":
    import argparse
    import logging

    from tqdm import tqdm

    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Convert Figure 2 spine density data to NWB format")
    parser.add_argument(
        "--stub-test",
        type=bool,
        default=True,
        help="Process only first 2 files per condition for testing (default: True). Use --stub-test=False for full processing.",
    )

    args = parser.parse_args()
    stub_test = args.stub_test

    # Control verbose output
    verbose = False  # Set to True for detailed output

    logging.getLogger("tifffile").setLevel(logging.ERROR)

    # Define the base path to the data
    base_path = Path("/home/heberto/development/surmeier-lab-to-nwb/link_to_raw_data/Figure 2_SF1A/Spine density/")

    # Create nwb_files directory at root level
    root_dir = Path(__file__).parent.parent.parent.parent.parent.parent  # Go up to repo root
    nwb_files_dir = root_dir / "nwb_files" / "spine_density" / "figure_2"
    nwb_files_dir.mkdir(parents=True, exist_ok=True)

    conditions = ["PD dSPN", "LID off-state dSPN", "control dSPN", "LID on-state dSPN"]
    for condition in conditions:
        condition_path = base_path / condition
        assert condition_path.exists(), f"Base path does not exist: {base_path}"

        # Get all session folders
        session_folders = [f for f in condition_path.iterdir() if f.is_dir()]
        session_folders.sort()

        # Apply stub_test filtering if enabled
        if stub_test:
            session_folders = session_folders[:2]

        # Use tqdm for progress bar  when verbose is disabled
        session_iterator = tqdm(session_folders, desc=f"Converting Figure2 SpineDensity {condition}", unit=" session")

        for session_folder_path in session_iterator:

            # Convert data to NWB format
            nwbfile = convert_data_to_nwb(
                session_folder_path=session_folder_path,
                condition=condition,
                verbose=verbose,
            )

            # Create output filename using session_id from nwbfile
            nwbfile_path = nwb_files_dir / f"{nwbfile.session_id}.nwb"

            # Write NWB file
            configure_and_write_nwbfile(nwbfile, nwbfile_path=nwbfile_path)
