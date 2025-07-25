"""
Figure 1 Somatic Excitability Conversion Script - Zhai et al. 2025
==============================================================

This script converts somatic excitability data from Figure 1 of Zhai et al. 2025
into NWB (Neurodata Without Borders) format. The data contains whole-cell patch-clamp
recordings from direct pathway SPNs (dSPNs) examining somatic excitability changes
in a mouse model of Parkinson's disease and levodopa-induced dyskinesia.

For detailed experimental context, protocols, data structure, and analysis methods,
see: /src/surmeier_lab_to_nwb/zhai2025/conversion_notes_folder/figure_1_conversion_notes.md
"""

import re
from pathlib import Path
from typing import Any, Dict

from neuroconv.tools import configure_and_write_nwbfile
from pynwb import NWBFile

from surmeier_lab_to_nwb.zhai2025.conversion_scripts.conversion_utils import (
    str_to_bool,
)
from surmeier_lab_to_nwb.zhai2025.conversion_scripts.somatic_excitability.somatic_excitability_utils import (
    convert_somatic_excitability_session_to_nwbfile,
)
from surmeier_lab_to_nwb.zhai2025.interfaces import (
    PROTOCOL_STEP_TO_CURRENT,
)


def parse_session_info_from_folder_name(recording_folder: Path) -> Dict[str, Any]:
    """
    Parse essential recording information from somatic excitability recording folder names.
    Session start time comes from XML metadata.

    Expected folder name format: cell[N]-[XXX] (e.g., cell1-001, cell2-015)

    The session folder structure is: MMDDYYYY_N/cell[N]-[XXX]/
    - MMDDYYYY_N: Date and cell number (e.g., 02022017_1)
    - cell[N]-[XXX]: Recording folder with cell number and protocol step

    Protocol steps:
    - 001-006: Hyperpolarizing currents (-120 to -20 pA)
    - 007-021: Depolarizing currents (+20 to +300 pA)
    - 022-026: Extended depolarizing currents (+320 to +400 pA, not all cells)

    Parameters
    ----------
    recording_folder : Path
        Path to the recording folder

    Returns
    -------
    Dict[str, Any]
        Dictionary containing recording information including date as datetime object
    """
    folder_name = recording_folder.name
    session_folder = recording_folder.parent

    # Parse session folder (parent) for date and cell info
    # Format: MMDDYYYY_N
    session_folder_name = session_folder.name
    if "_" in session_folder_name:
        _, cell_number = session_folder_name.split("_")
    else:
        raise ValueError(f"Could not parse session folder name: {session_folder_name}")

    # Note: Date will be extracted from XML session start time, not from folder name

    # Parse recording folder name for protocol step
    # Format: cell[N]-[XXX] or cell[N]_SCH-[XXX] or cell[N]_sul-[XXX]
    # Use regex to handle optional "_SCH" or "_sul" suffix
    pattern = r"cell(\d+)(?:_SCH|_sul)?-(\d+)"
    match = re.match(pattern, folder_name)

    if not match:
        raise ValueError(f"Could not parse recording folder name: {folder_name}")

    cell_number_from_recording = match.group(1)
    protocol_step = match.group(2)

    # Verify cell number matches session
    if cell_number_from_recording != cell_number:
        raise ValueError(f"Cell number mismatch: session {cell_number} vs recording {cell_number_from_recording}")

    # Get current value for this protocol step using shared constant
    if protocol_step not in PROTOCOL_STEP_TO_CURRENT:
        raise ValueError(
            f"Protocol step {protocol_step} not found in PROTOCOL_STEP_TO_CURRENT mapping for {folder_name}"
        )

    current_pA = PROTOCOL_STEP_TO_CURRENT[protocol_step]

    # Format current with consistent notation
    current_formatted = f"{current_pA:+04d}pA"  # e.g., "+020pA", "-120pA"

    # Create step order (sequential numbering for analysis)
    step_order = int(protocol_step)

    return {
        "cell_number": cell_number,
        "protocol_step": protocol_step,
        "current_pA": current_pA,
        "current_formatted": current_formatted,
        "step_order": step_order,
        "recording_folder_name": folder_name,
    }


def convert_session_to_nwbfile(session_folder_path: Path, condition: str) -> NWBFile:
    """
    Convert a single session of Figure 1 somatic excitability data to NWB format.

    This is a wrapper function that calls the shared conversion function with
    Figure 1-specific configuration and session ID parameters.

    Parameters
    ----------
    session_folder_path : Path
        Path to the session folder containing cell recordings
    condition : str
        Experimental condition (e.g., "LID off-state", "LID on-state", "LID on-state with SCH")

    Returns
    -------
    NWBFile
        NWB file with the converted data
    """
    # Configuration for Figure 1 somatic excitability experiments
    figure_1_config = {
        "parse_function": parse_session_info_from_folder_name,
        "metadata_key": "figure_1_somatic_excitability",
    }

    # Local mapping for Figure 1 conditions to revised schema tokens
    figure_1_mappings = {
        "LID off-state": {"state": "OffState", "pharm": "none"},
        "LID on-state": {"state": "OnState", "pharm": "none"},
        "LID on-state with SCH": {"state": "OnState", "pharm": "D1RaSch"},
    }

    state = figure_1_mappings[condition]["state"]
    pharmacology = figure_1_mappings[condition]["pharm"]

    # Build session ID parameters using revised schema
    session_id_parameters = {
        "fig": "F1",
        "meas_comp": "SomExc",  # Somatic excitability
        "cell_type": "dSPN",  # Direct pathway SPN
        "state": state,
        "pharm": pharmacology,
        "geno": "WT",  # Wild-type for all Figure 1
    }

    return convert_somatic_excitability_session_to_nwbfile(
        session_folder_path=session_folder_path,
        condition=condition,
        figure_config=figure_1_config,
        session_id_parameters=session_id_parameters,
    )


if __name__ == "__main__":
    import argparse
    import logging
    import warnings

    from tqdm import tqdm

    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Convert Figure 1 somatic excitability data to NWB format")

    parser.add_argument(
        "--stub-test",
        type=str_to_bool,
        default=True,
        help="Process only first 2 files per condition for testing (default: True). Use --stub-test=False for full processing.",
    )

    args = parser.parse_args()
    stub_test = args.stub_test

    # Suppress tifffile warnings
    logging.getLogger("tifffile").setLevel(logging.ERROR)

    # Suppress specific warnings
    warnings.filterwarnings("ignore", message="invalid value encountered in divide")
    warnings.filterwarnings("ignore", message=".*no datetime before year 1.*")

    # Set the base path to your data
    base_path = Path("./link_to_raw_data/Figure 1/Somatic excitability")

    # Create nwb_files directory at root level
    root_dir = Path(__file__).parent.parent.parent.parent.parent.parent  # Go up to repo root
    nwb_files_dir = root_dir / "nwb_files" / "somatic_excitability" / "figure_1"
    nwb_files_dir.mkdir(parents=True, exist_ok=True)

    # Figure 1 somatic excitability conditions
    conditions = ["LID off-state", "LID on-state", "LID on-state with SCH"]

    for condition in conditions:
        condition_path = base_path / condition

        if not condition_path.exists():
            raise FileNotFoundError(f"Expected condition path does not exist: {condition_path}")

        # Get all session folders (each session = one cell)
        session_folders = [f for f in condition_path.iterdir() if f.is_dir()]
        session_folders.sort()

        # Apply stub_test filtering if enabled
        if stub_test:
            session_folders = session_folders[:2]

        # Use tqdm for progress bar
        session_iterator = tqdm(
            session_folders, desc=f"Converting Figure1 SomaticExcitability {condition}", unit=" session"
        )

        for session_folder in session_iterator:

            # Convert session data to NWB format with time alignment
            nwbfile = convert_session_to_nwbfile(
                session_folder_path=session_folder,
                condition=condition,
            )

            # Create output filename using session_id from nwbfile
            nwbfile_path = nwb_files_dir / f"{nwbfile.session_id}.nwb"

            # Write NWB file
            configure_and_write_nwbfile(nwbfile, nwbfile_path=nwbfile_path)
