"""
Figure 6 Dendritic Excitability Conversion Script - Zhai et al. 2025
================================================================

This script converts dendritic excitability data from Figure 6 of Zhai et al. 2025
into NWB (Neurodata Without Borders) format. The data combines patch-clamp recordings
with simultaneous two-photon imaging to examine dendritic excitability changes
in a mouse model of Parkinson's disease and levodopa-induced dyskinesia.

For detailed experimental context, protocols, data structure, and analysis methods,
see: /src/surmeier_lab_to_nwb/zhai2025/conversion_notes_folder/figure_6_conversion_notes.md
"""

import re
import warnings
from datetime import datetime
from pathlib import Path
from typing import Any, Dict

from neuroconv.tools import configure_and_write_nwbfile
from pynwb import NWBFile

from surmeier_lab_to_nwb.zhai2025.conversion_scripts.conversion_utils import (
    format_condition,
    str_to_bool,
)
from surmeier_lab_to_nwb.zhai2025.conversion_scripts.dendritic_excitability.dendritic_excitability_utils import (
    convert_dendritic_excitability_session_to_nwbfile,
)


def _normalize_figure6_folder_name(folder_name: str) -> str:
    """
    Normalize Figure 6 dendritic excitability folder names to handle dataset-specific anomalies.

    This method is very specific to Figure 6 datasets because some sessions have inverted
    naming patterns that deviate from the standard format used in other figures.

    Standard format: [date]_Cell[N]_[location][N]_trio-[trial]
    Examples: 02182020_Cell2_dist1_trio-001, 08122020_Cell1_prox2_trio-003

    Anomalous datasets with inverted structure:
    - 0425a: 04252023_cell1_trio_dist1-001 (should be 04252023_cell1_dist1_trio-001)
    - 0425b: 04252023_cell2_trio_prox2-003 (should be 04252023_cell2_prox2_trio-003)
    - 0430a: 04302023_cell1_trio_dist-001 (should be 04302023_cell1_dist1_trio-001)
    - 0430b: 04302023_cell2_trio_dist2-004 (should be 04302023_cell2_dist2_trio-004)

    The anomalous pattern has "trio_[location]" instead of "[location]_trio" and sometimes
    missing location numbers (e.g., "trio_dist" instead of "dist1_trio").

    Parameters
    ----------
    folder_name : str
        Original folder name from the dataset

    Returns
    -------
    str
        Normalized folder name following the standard format
    """
    # Handle trio_location[N] -> location[N]_trio pattern inversions
    # Use regex to properly capture location numbers and reorder them
    folder_name = re.sub(r"_trio_(dist|prox)(\d+)", r"_\1\2_trio", folder_name)

    # Handle missing location numbers in the inverted pattern
    # trio_dist-XXX -> dist1_trio-XXX (assumes location number 1 when missing)
    folder_name = re.sub(r"_trio_dist-", "_dist1_trio-", folder_name)
    folder_name = re.sub(r"_trio_prox-", "_prox1_trio-", folder_name)

    return folder_name


def parse_session_info_from_folder_name(recording_folder: Path) -> Dict[str, Any]:
    """
    Parse essential recording information from Figure 6 dendritic excitability recording folder names.
    Only extracts location and trial information - session start time comes from XML metadata.

    Expected folder name format: [date]_Cell[cell_number]_[location][location_number][variant]_[experiment_type]-[trial_number]
    Examples:
    - 02182020_Cell2_dist1_trio-001 (YYYYMMDD format)
    - 02182020_Cell3_prox2_trio-003 (YYYYMMDD format with different cell/location)

    Date formats:
    - YYYYMMDD: Year, month, day (e.g., 02182020 = February 18, 2020)
    - MMDDYYYY: Month, day, year (e.g., 02182016 = February 18, 2016)
    Format is auto-detected based on first digit (if '2', assumes YYYYMMDD format)

    Experiment types:
    - trio: Three current injections protocol (main experimental condition)
    - bsl: Baseline recording protocol (control/reference measurements)

    Variant:
    - real: Optional suffix indicating a re-recording or alternate version of the same location
    - (empty): Standard recording without variant

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

    # Normalize folder name to handle Figure 6 dataset-specific anomalies
    folder_name = _normalize_figure6_folder_name(folder_name)

    # Parse using regex pattern for dendritic recordings
    # Handle variations like "dist1real", "realprox1" and different experiment types like "trio" or "bsl"
    # Some folders may not have experiment type, make it optional (with or without underscore)
    # Handle case variations like "Cell" vs "cell"
    if folder_name.startswith("00"):
        folder_name = folder_name[
            1:
        ]  # There is one case,  LID on-state/0718a/007182016_cell1_dist1_trio-001, where there is an extra leading zero

    if "_sul" in folder_name:
        folder_name = folder_name.replace("_sul", "")  # Remove "_sul" suffix if present

    if "proxt" in folder_name:
        folder_name = folder_name.replace("proxt", "prox")  # LID on-state with sul/0508c has this typo

    # Parse dendritic excitability recording folder names
    # Format: [date]_Cell[cell_number]_[location][location_number][variant]_[experiment_type][_SCH]-[trial_number]
    pattern = (
        r"(?P<date>\d{8})"  # Date (8 digits): YYYYMMDD or MMDDYYYY
        r"_[Cc]ell(?P<cell_number>\d+)"  # Cell number: Cell1, cell2, etc.
        r"_(?:(?P<variant_prefix>real))?"  # Optional variant prefix: "real" before location
        r"(?P<location>dist|prox|)"  # Location: "dist" (distal) or "prox" (proximal)
        r"(?P<location_number>\d+)"  # Location number: 1, 2, etc.
        r"(?:(?P<variant_suffix>real))?"  # Optional variant suffix: "real" after location
        r"(?:_(?P<experiment_type>trio|bsl))?"  # Optional experiment type: "trio" or "bsl"
        r"(?:_(?P<sch_suffix>SCH))?"  # Optional SCH suffix for SCH condition
        r"-(?P<trial_number>\d+)"  # Trial number: 001, 002, etc.
    )
    match = re.match(pattern, folder_name)

    if not match:
        raise ValueError(f"Could not parse recording folder name: {folder_name}")

    date_str = match.group("date")
    cell_number = match.group("cell_number")
    location = match.group("location")
    location_number = match.group("location_number")
    variant_prefix = match.group("variant_prefix") or ""  # "real" before location
    variant_suffix = match.group("variant_suffix") or ""  # "real" after location
    variant = variant_prefix or variant_suffix  # Use whichever is present
    experiment_type = match.group("experiment_type") or ""  # "trio", "bsl", or empty
    trial_number = match.group("trial_number")

    # Handle different date formats based on first digit
    if date_str[0] == "2":
        # YYYYMMDD format (e.g., 20160523)
        year = int(date_str[:4])
        month = int(date_str[4:6])
        day = int(date_str[6:8])
    else:
        # MMDDYYYY format (e.g., 05232016)
        month = int(date_str[:2])
        day = int(date_str[2:4])
        year = int(date_str[4:8])

    # Validate date components
    if month < 1 or month > 12:
        raise ValueError(f"Invalid month '{month}' in date '{date_str}' from folder '{folder_name}'")
    if day < 1 or day > 31:
        raise ValueError(f"Invalid day '{day}' in date '{date_str}' from folder '{folder_name}'")

    # Determine location description
    location_full = "Proximal" if location == "prox" else "Distal"
    variant_suffix = f" {variant.capitalize()}" if variant else ""
    location_description = f"{location_full} dendrite {location_number}{variant_suffix}"
    base_line_experiment_type = "Baseline" if experiment_type == "bsl" else ""

    # Create datetime object for the date
    session_date = datetime(year, month, day)

    return {
        "cell_number": cell_number,
        "location": location,
        "location_number": location_number,
        "location_full": location_full,
        "location_description": location_description,
        "experiment_type": experiment_type,
        "trial_number": trial_number,
        "variant": variant,
        "base_line_experiment_type": base_line_experiment_type,
        "date": session_date,
    }


def convert_session_to_nwbfile(session_folder_path: Path, condition: str, verbose: bool = False) -> NWBFile:
    """
    Convert all Figure 6 dendritic excitability recordings from a session folder to a single NWB format file.

    This is a wrapper function that calls the shared conversion function with
    Figure 6-specific configuration and session ID parameters.

    Parameters
    ----------
    session_folder_path : Path
        Path to the session folder (corresponds to a single animal/subject)
    condition : str
        Experimental condition (e.g., "control", "M1R antagonist")
    verbose : bool, default=False
        Enable verbose output showing detailed processing information

    Returns
    -------
    NWBFile
        NWB file with all the converted data from the session

    Notes
    -----
    Each session folder corresponds to a single animal. The structure is:
    session_folder/recording_folder/files
    where recording_folder represents different recordings from the same animal.
    """
    # Configuration for Figure 6 dendritic excitability experiments
    figure_6_config = {
        "parse_function": parse_session_info_from_folder_name,
        "metadata_key": "figure_6_dendritic_excitability",
    }

    # Local mapping for Figure 6 conditions to revised schema tokens
    figure_6_mappings = {
        "control": {"state": "OffState", "pharm": "none"},
        "M1R antagonist": {"state": "OffState", "pharm": "M1RaThp"},
    }

    state = figure_6_mappings[condition]["state"]
    pharmacology = figure_6_mappings[condition]["pharm"]

    # Build session ID parameters using revised schema
    session_id_parameters = {
        "fig": "F6",
        "meas_comp": "DendExc",  # Dendritic excitability
        "cell_type": "iSPN",  # Indirect pathway SPN
        "state": state,
        "pharm": pharmacology,
        "geno": "WT",  # Wild-type for all Figure 6
    }

    return convert_dendritic_excitability_session_to_nwbfile(
        session_folder_path=session_folder_path,
        condition=condition,
        figure_config=figure_6_config,
        session_id_parameters=session_id_parameters,
        verbose=verbose,
    )


if __name__ == "__main__":
    import argparse
    import logging
    import warnings

    from tqdm import tqdm

    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Convert Figure 6 dendritic excitability data to NWB format")
    parser.add_argument(
        "--stub-test",
        type=str_to_bool,
        default=True,
        help="Process only first 2 files per condition for testing (default: True). Use --stub-test=False for full processing.",
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
    base_path = Path("./link_to_raw_data/Figure 6/Dendritic excitability")

    # Create nwb_files directory at root level
    root_dir = Path(__file__).parent.parent.parent.parent.parent.parent  # Go up to repo root
    nwb_files_dir = root_dir / "nwb_files" / "dendritic_excitability" / "figure_6"
    nwb_files_dir.mkdir(parents=True, exist_ok=True)

    # Figure 6 dendritic conditions
    conditions = ["control", "M1R antagonist"]

    for condition in conditions:
        condition_path = base_path / condition

        if not condition_path.exists():
            raise FileNotFoundError(f"Expected condition path does not exist: {condition_path}")

        # Get all session folders (e.g., 0217a, 0218b, etc.)
        # Note: Each session folder corresponds to a single animal
        session_folders = [f for f in condition_path.iterdir() if f.is_dir()]
        session_folders.sort()

        # Apply stub_test filtering if enabled
        if stub_test:
            session_folders = session_folders[:2]

        # Use tqdm for progress bar when verbose is disabled
        session_iterator = tqdm(
            session_folders,
            desc=f"Converting Figure6 DendriticExcitability {format_condition[condition]['human_readable']}",
            disable=verbose,
            unit=" session",
        )

        for session_folder in session_iterator:

            # Convert all recordings from this session to NWB format
            nwbfile = convert_session_to_nwbfile(
                session_folder_path=session_folder,
                condition=condition,
                verbose=verbose,
            )

            # Create output filename
            condition_safe = condition.replace(" ", "_").replace("-", "_").replace("(", "").replace(")", "")
            nwbfile_path = nwb_files_dir / f"{nwbfile.session_id}.nwb"

            # Write NWB file
            configure_and_write_nwbfile(nwbfile, nwbfile_path=nwbfile_path)
