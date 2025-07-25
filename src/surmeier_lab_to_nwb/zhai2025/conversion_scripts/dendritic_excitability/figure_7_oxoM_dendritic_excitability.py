"""
Figure 7 OxoM Dendritic Excitability Conversion Script - Zhai et al. 2025
=====================================================================

This script converts OxoM-induced dendritic excitability data from Figure 7 of
Zhai et al. 2025 into NWB (Neurodata Without Borders) format. The data combines
patch-clamp recordings with simultaneous two-photon imaging to examine the effects
of muscarinic agonist oxotremorine-M on dendritic excitability.

For detailed experimental context, protocols, data structure, and analysis methods,
see: /src/surmeier_lab_to_nwb/zhai2025/conversion_notes_folder/figure_7_conversion_notes.md
"""

import re
import warnings
from pathlib import Path
from typing import Any, Dict

from neuroconv.tools import configure_and_write_nwbfile
from pynwb import NWBFile

from surmeier_lab_to_nwb.zhai2025.conversion_scripts.conversion_utils import (
    str_to_bool,
)
from surmeier_lab_to_nwb.zhai2025.conversion_scripts.dendritic_excitability.dendritic_excitability_utils import (
    convert_dendritic_excitability_session_to_nwbfile,
)


def parse_session_info_from_folder_name(recording_folder: Path) -> Dict[str, Any]:
    """
    Parse essential recording information from Figure 7E oxotremorine-M dendritic excitability recording folder names.
    Session start time comes from XML metadata.

    Expected folder name format: [MMDDYYYY]_Cell[N]_[location][location_number]_trio-[trial_number]
    (e.g., 03012020_Cell1_dist1_trio-001, 03012020_Cell1_prox1_trio-004)

    The session folder structure is: [MMDD][letter][digit]/[MMDDYYYY]_Cell[N]_[location][location_number]_trio-[trial]/
    - [MMDD][letter][digit]: Date and animal identifier (e.g., 0301a, 0301b)
    - [MMDDYYYY]_Cell[N]_[location][location_number]_trio-[trial]: Recording with dendritic location and trial info

    Protocol: Three 2 nA current injections, 2 ms each, at ~50 Hz for dendritic calcium imaging
    oxoM Protocol Convention:
    - trio-001, -002, -003: BEFORE oxotremorine-M application (baseline measurements)
    - trio-004, -005, -006: AFTER oxotremorine-M application (drug response)
    - trio-007, -008, -009: AFTER oxotremorine-M application (additional recordings if present)

    Parameters
    ----------
    recording_folder : Path
        Path to the recording folder

    Returns
    -------
    Dict[str, Any]
        Dictionary containing recording information
    """
    folder_name = recording_folder.name
    session_folder = recording_folder.parent

    # Parse session folder (parent) for date and animal info
    # Format: [MMDD][letter][digit] (e.g., 0301a, 0301b, 0302a)
    session_folder_name = session_folder.name

    # Extract date and animal identifier
    # Pattern: 4 digits (MMDD) followed by a letter and optionally a digit
    session_pattern = r"(\d{4})([a-z])(\d*)"
    session_match = re.match(session_pattern, session_folder_name)

    if not session_match:
        raise ValueError(f"Could not parse session folder name: {session_folder_name}")

    animal_letter = session_match.group(2)  # letter (a, b, c, etc.)
    animal_digit = session_match.group(3) or ""  # optional digit
    animal_id = f"{animal_letter}{animal_digit}" if animal_digit else animal_letter

    # Note: Date will be extracted from XML session start time, not from folder name

    # Parse recording folder name for dendritic recording info
    # Format: [MMDDYYYY]_Cell[N]_[location][location_number]_trio-[trial_number]
    # Example: 03012020_Cell1_dist1_trio-001, 03012020_Cell1_prox1_trio-004
    pattern = r"(\d{8})_Cell(\d+)_([a-z]+)(\d+)_trio-(\d+)"
    match = re.match(pattern, folder_name)

    if not match:
        raise ValueError(f"Could not parse recording folder name: {folder_name}")

    date_from_folder = match.group(1)  # MMDDYYYY (for reference, but we'll use XML time)
    cell_number = match.group(2)
    location_type = match.group(3)  # "dist" or "prox"
    location_number = match.group(4)
    trial_number = match.group(5)

    # Determine recording phase based on trial number
    trial_num = int(trial_number)
    if trial_num <= 3:
        recording_phase = "baseline"
        phase_description = "before oxotremorine-M application"
    else:
        recording_phase = "post_oxoM"
        phase_description = "after oxotremorine-M application"

    # Create location identifier
    location_id = f"{location_type}{location_number}"

    # Determine approximate distance from soma based on location type
    if location_type == "prox":
        approximate_distance_um = 40  # proximal dendrites ~40 μm from soma
    elif location_type == "dist":
        approximate_distance_um = 90  # distal dendrites ~90 μm from soma
    else:
        approximate_distance_um = None

    # Create full location description
    location_full = "Proximal" if location_type == "prox" else "Distal"
    location_description = f"{location_full} dendrite {location_number}"

    # Map to standard schema used by shared function
    return {
        "cell_number": cell_number,
        "animal_id": animal_id,
        "location": location_type,  # Standard field name for shared function
        "location_number": location_number,
        "location_id": location_id,
        "location_full": location_full,
        "location_description": location_description,
        "trial_number": trial_number,
        "recording_phase": recording_phase,
        "phase_description": phase_description,
        "approximate_distance_um": approximate_distance_um,
        "recording_folder_name": folder_name,
        "experiment_type": "trio",  # Standard field for shared function
        "base_line_experiment_type": "",  # Standard field for shared function
        "variant": "",  # Standard field for shared function
    }


def convert_session_to_nwbfile(session_folder_path: Path, condition: str, verbose: bool = False) -> NWBFile:
    """
    Convert a single session of Figure 7E oxotremorine-M dendritic excitability data to NWB format with time alignment.

    This is a wrapper function that calls the shared conversion function with
    Figure 7 oxoM-specific configuration and session ID parameters.

    Parameters
    ----------
    session_folder_path : Path
        Path to the session folder containing recordings (both baseline and post-oxoM)
    condition : str
        Experimental condition (e.g., "WT oxoM treatment", "CDGI KO oxoM treatment")
    verbose : bool, default=False
        Enable verbose output

    Returns
    -------
    NWBFile
        NWB file with all recordings from the session (baseline and post-oxoM)

    Notes
    -----
    Each session includes recordings before and after oxotremorine-M application.
    Trial numbers 1-3 are baseline, 4+ are post-drug.
    """
    # Configuration for Figure 7 oxoM dendritic excitability experiments
    figure_7_oxom_config = {
        "parse_function": parse_session_info_from_folder_name,
        "metadata_key": "figure_7_oxoM_dendritic_excitability",
    }

    # Map condition to appropriate schema tokens
    condition_mappings = {
        "WT oxoM treatment": {"geno": "WT"},
        "CDGI KO oxoM treatment": {"geno": "CDGIKO"},
    }

    if condition not in condition_mappings:
        raise ValueError(f"Unknown condition: {condition}")

    geno_token = condition_mappings[condition]["geno"]

    # Build session ID parameters using revised schema
    session_id_parameters = {
        "fig": "F7",
        "meas_comp": "DendExc",  # Dendritic excitability
        "cell_type": "iSPN",  # Indirect pathway SPN
        "state": "OffState",  # OFF state for oxoM experiments
        "pharm": "M1RaOxoM",  # Muscarinic agonist oxotremorine-M
        "geno": geno_token,
    }

    return convert_dendritic_excitability_session_to_nwbfile(
        session_folder_path=session_folder_path,
        condition=condition,
        figure_config=figure_7_oxom_config,
        session_id_parameters=session_id_parameters,
        verbose=verbose,
    )


if __name__ == "__main__":
    import argparse
    import logging

    from tqdm import tqdm

    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Convert Figure 7E oxoM dendritic excitability data to NWB format")
    parser.add_argument(
        "--stub-test",
        type=str_to_bool,
        default=True,
        help="Process only first 2 files per genotype for testing (default: True). Use --stub-test=False for full processing.",
    )

    args = parser.parse_args()
    stub_test = args.stub_test

    verbose = False

    # Suppress tifffile warnings
    logging.getLogger("tifffile").setLevel(logging.ERROR)

    # Suppress specific warnings
    warnings.filterwarnings("ignore", message="invalid value encountered in divide")
    warnings.filterwarnings("ignore", message=".*no datetime before year 1.*")

    # Set the base path to your data
    base_path = Path("./link_to_raw_data/Figure 7/oxoM on DE")

    # Create nwb_files directory at root level
    root_dir = Path(__file__).parent.parent.parent.parent.parent.parent  # Go up to repo root
    nwb_files_dir = root_dir / "nwb_files" / "dendritic_excitability" / "figure_7_oxoM"
    nwb_files_dir.mkdir(parents=True, exist_ok=True)

    # Figure 7E oxoM folder-to-condition mapping
    folder_to_condition = {
        "WT": "WT oxoM treatment",
        "KO": "CDGI KO oxoM treatment",
    }

    for folder_name, condition in folder_to_condition.items():
        folder_path = base_path / folder_name

        if not folder_path.exists():
            raise FileNotFoundError(f"Expected folder path does not exist: {folder_path}")

        # Get all session folders (e.g., 0301a, 0301b, 0302a)
        # Each session folder contains recordings for both baseline and post-oxoM phases
        session_folders = [f for f in folder_path.iterdir() if f.is_dir()]
        session_folders.sort()

        # Apply stub_test filtering if enabled
        if stub_test:
            session_folders = session_folders[:2]

        # Use tqdm for progress bar
        session_iterator = tqdm(
            session_folders,
            desc=f"Converting Figure7E OxoM DendriticExcitability {folder_name}",
            unit=" session",
        )

        for session_folder in session_iterator:

            # Convert all recordings from this session to NWB format
            nwbfile = convert_session_to_nwbfile(
                session_folder_path=session_folder,
                condition=condition,
                verbose=verbose,
            )

            # Create output filename using session_id from nwbfile
            nwbfile_path = nwb_files_dir / f"{nwbfile.session_id}.nwb"

            # Write NWB file
            configure_and_write_nwbfile(nwbfile, nwbfile_path=nwbfile_path)
