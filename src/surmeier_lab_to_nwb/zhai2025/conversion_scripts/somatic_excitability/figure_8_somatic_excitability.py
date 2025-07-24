"""
Figure 8 Somatic Excitability Conversion Script - Zhai et al. 2025
==============================================================

This script converts somatic excitability data from Figure 8 of Zhai et al. 2025
into NWB (Neurodata Without Borders) format. The data contains whole-cell patch-clamp
recordings examining somatic excitability changes in a mouse model of Parkinson's
disease and levodopa-induced dyskinesia.

For detailed experimental context, protocols, data structure, and analysis methods,
see: /src/surmeier_lab_to_nwb/zhai2025/conversion_notes_folder/figure_8_conversion_notes.md
"""

import re
import uuid
from pathlib import Path
from typing import Any, Dict

from neuroconv.tools import configure_and_write_nwbfile
from neuroconv.utils import dict_deep_update, load_dict_from_file
from pynwb import NWBFile
from pynwb.file import Subject

from surmeier_lab_to_nwb.zhai2025.conversion_scripts.somatic_excitability.utils import (
    build_somatic_icephys_table_structure,
)
from surmeier_lab_to_nwb.zhai2025.interfaces import (
    PROTOCOL_STEP_TO_CURRENT,
    PrairieViewCurrentClampInterface,
)


def parse_session_info_from_folder_name(recording_folder: Path) -> Dict[str, Any]:
    """
    Parse essential recording information from M1R CRISPR somatic excitability recording folder names.
    Session start time comes from XML metadata.

    Expected folder name format: cell[N]-[XXX] (e.g., cell1-001, cell2-015)

    The session folder structure is: YYYYMMDD[X]/cell[N]-[XXX]/
    - YYYYMMDD[X]: Date and animal identifier (e.g., 20221004b, 20221012a)
    - cell[N]-[XXX]: Recording folder with cell number and protocol step

    Protocol steps:
    - 001-006: Hyperpolarizing currents (-120 to -20 pA)
    - 007-021: Depolarizing currents (+20 to +300 pA)

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

    # Parse session folder (parent) for date and animal info
    # Format: YYYYMMDD[X] (e.g., 20221004b, 20221012a)
    session_folder_name = session_folder.name
    # Extract animal identifier (letter at the end)
    if len(session_folder_name) >= 9 and session_folder_name[-1].isalpha():
        animal_id = session_folder_name[-1]
        date_part = session_folder_name[:-1]
    else:
        raise ValueError(f"Could not parse session folder name: {session_folder_name}")

    # Parse recording folder name for protocol step
    # Format: cell[N]-[XXX]
    pattern = r"cell(\d+)-(\d+)"
    match = re.match(pattern, folder_name)

    if not match:
        raise ValueError(f"Could not parse recording folder name: {folder_name}")

    cell_number = match.group(1)
    protocol_step = match.group(2)

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
        "animal_id": animal_id,
        "date_part": date_part,
        "protocol_step": protocol_step,
        "current_pA": current_pA,
        "current_formatted": current_formatted,
        "step_order": step_order,
        "recording_folder_name": folder_name,
    }


def convert_session_to_nwbfile(session_folder_path: Path, condition: str) -> NWBFile:
    """
    Convert a single session of Figure 8 M1R CRISPR somatic excitability data to NWB format with time alignment.

    Parameters
    ----------
    session_folder_path : Path
        Path to the session folder containing cell recordings
    condition : str
        Experimental condition (e.g., "M1R CRISPR", "interleaved control")

    Returns
    -------
    NWBFile
        NWB file with the converted data

    Notes
    -----
    This function implements temporal alignment by extracting precise timestamps from XML files
    and calculating t_start offsets for each recording relative to the earliest recording.
    """

    # Find all recording folders (current steps) for this cell
    recording_folders = [f for f in session_folder_path.iterdir() if f.is_dir()]
    recording_folders.sort()

    if not recording_folders:
        raise ValueError(f"No recording folders found in session folder: {session_folder_path}")

    # Parse session information from first recording folder (all should have same session info)
    first_recording_info = parse_session_info_from_folder_name(recording_folders[0])
    session_info = {
        "cell_number": first_recording_info["cell_number"],
        "animal_id": first_recording_info["animal_id"],
        "date_part": first_recording_info["date_part"],
    }

    # Calculate recording IDs, session start times, and create interface mappings
    session_start_times = []  # (timestamp, recording_folder, recording_id)
    recording_id_to_info = {}
    recording_id_to_folder = {}
    t_starts = {}  # t_starts[recording_id] = t_start_offset

    for recording_folder in recording_folders:
        # Parse recording information using unified function
        recording_info = parse_session_info_from_folder_name(recording_folder)

        # Create unique recording ID
        recording_id = f"Cell{recording_info['cell_number']}{recording_info['current_formatted']}"

        # Find XML file for this recording
        xml_file = recording_folder / f"{recording_folder.name}_Cycle00001_VoltageRecording_001.xml"
        if not xml_file.exists():
            raise FileNotFoundError(f"Expected XML file does not exist: {xml_file}")

        # Get session start time from XML metadata
        session_start_time = PrairieViewCurrentClampInterface.get_session_start_time_from_file(xml_file)
        if session_start_time is None:
            raise ValueError(f"Could not extract session start time from {xml_file}")

        # Store mappings
        recording_id_to_info[recording_id] = recording_info
        recording_id_to_folder[recording_id] = recording_folder
        session_start_times.append((session_start_time, recording_folder, recording_id))

    if not session_start_times:
        raise ValueError(f"No valid recordings found in session folder: {session_folder_path}")

    # Find the earliest session start time across all recordings
    earliest_time = min(session_start_times, key=lambda x: x[0])[0]

    # Calculate t_start offsets for temporal alignment
    for start_time, _, recording_id in session_start_times:
        # Calculate offset relative to overall session start time
        t_start_offset = (start_time - earliest_time).total_seconds()
        t_starts[recording_id] = t_start_offset

    # Use earliest time as session start time for NWB file
    session_start_time = earliest_time

    # Create session ID following pattern from figure_1_somatic_excitability.py
    condition_to_camel_case = {
        "LID off-state": "LIDOffState",
        "LID on-state": "LIDOnState",
        "LID on-state with SCH": "LIDOnStateWithSchD1Antagonist",
    }

    timestamp = session_start_time.strftime("%Y%m%d%H%M%S")
    clean_condition = condition_to_camel_case.get(condition, condition.replace(" ", "").replace("-", ""))
    base_session_id = f"Figure8++SomaticExcitability++{clean_condition}++{timestamp}"
    script_specific_id = f"Cell++{session_info['cell_number']}++Animal++{session_info['animal_id']}"
    session_id = f"{base_session_id}++{script_specific_id}"

    # Load general and session-specific metadata from YAML files
    general_metadata_path = Path(__file__).parent.parent.parent / "general_metadata.yaml"
    general_metadata = load_dict_from_file(general_metadata_path)

    session_metadata_path = Path(__file__).parent.parent.parent / "session_specific_metadata.yaml"
    session_metadata_template = load_dict_from_file(session_metadata_path)
    script_template = session_metadata_template["figure_8_somatic_excitability"]

    # Determine cell type based on condition
    if condition == "M1R CRISPR":
        cell_type = "iSPN"
    else:  # interleaved control
        cell_type = "iSPN"

    # Create session-specific metadata from template with runtime substitutions
    session_specific_metadata = {
        "NWBFile": {
            "session_description": script_template["NWBFile"]["session_description"].format(
                condition=condition, cell_number=session_info["cell_number"]
            ),
            "identifier": str(uuid.uuid4()),
            "session_start_time": session_start_time,
            "experiment_description": script_template["NWBFile"]["experiment_description"].format(
                cell_type=cell_type, condition=condition, num_current_steps=len(recording_folders)
            ),
            "session_id": session_id,
            "surgery": general_metadata["NWBFile"]["surgery"] + " " + script_template["NWBFile"]["surgery_addition"],
            "keywords": script_template["NWBFile"]["keywords"],
        },
        "Subject": {
            "subject_id": f"{cell_type}_M1R_CRISPR_mouse_{session_id}",
            "description": script_template["Subject"]["description"].format(cell_number=session_info["cell_number"]),
            "genotype": script_template["Subject"]["genotype"],
        },
    }

    # Deep merge with general metadata
    metadata = dict_deep_update(general_metadata, session_specific_metadata)

    # Create NWB file with merged metadata
    nwbfile = NWBFile(
        session_description=metadata["NWBFile"]["session_description"],
        identifier=metadata["NWBFile"]["identifier"],
        session_start_time=metadata["NWBFile"]["session_start_time"],
        experimenter=metadata["NWBFile"]["experimenter"],
        lab=metadata["NWBFile"]["lab"],
        institution=metadata["NWBFile"]["institution"],
        experiment_description=metadata["NWBFile"]["experiment_description"],
        session_id=metadata["NWBFile"]["session_id"],
        surgery=metadata["NWBFile"]["surgery"],
        pharmacology=metadata["NWBFile"]["pharmacology"],
        slices=metadata["NWBFile"]["slices"],
        keywords=metadata["NWBFile"]["keywords"],
    )

    # Create subject using merged metadata
    subject = Subject(
        subject_id=metadata["Subject"]["subject_id"],
        species=metadata["Subject"]["species"],
        strain=metadata["Subject"]["strain"],
        description=metadata["Subject"]["description"],
        genotype=metadata["Subject"]["genotype"],
        sex=metadata["Subject"]["sex"],
        age=metadata["Subject"]["age"],
    )
    nwbfile.subject = subject

    # Add custom columns to intracellular recording table for M1R CRISPR somatic experiment annotations
    intracellular_recording_table = nwbfile.get_intracellular_recordings()
    intracellular_recording_table.add_column(
        name="stimulus_current_pA",
        description="The current stimulus applied during the recording in picoamps",
    )
    intracellular_recording_table.add_column(
        name="protocol_step", description="Protocol step number (e.g., '001', '002', etc.)"
    )
    intracellular_recording_table.add_column(
        name="recording_id", description="Full recording identifier containing step and current information"
    )
    intracellular_recording_table.add_column(
        name="animal_id", description="Animal identifier for tracking across experimental sessions"
    )

    # Data structures for tracking icephys table indices
    recording_indices = []  # Store all intracellular recording indices

    # Process each recording using the calculated recording IDs and temporal alignment
    for recording_id, recording_folder in recording_id_to_folder.items():
        recording_info = recording_id_to_info[recording_id]
        xml_file = recording_folder / f"{recording_folder.name}_Cycle00001_VoltageRecording_001.xml"

        # Create interface for this recording
        icephys_metadata_key = f"PrairieView{recording_id}"
        interface = PrairieViewCurrentClampInterface(file_path=xml_file, icephys_metadata_key=icephys_metadata_key)

        # Apply temporal alignment offset
        interface.set_aligned_starting_time(t_starts[recording_id])

        # Get and update interface metadata
        interface_metadata = interface.get_metadata()

        # Update electrode description for M1R CRISPR somatic recording (consistent name for all recordings)
        electrode_name = "IntracellularElectrode"
        interface_metadata["Icephys"]["IntracellularElectrodes"][icephys_metadata_key].update(
            {
                "name": electrode_name,
                "description": (
                    f"Whole-cell patch clamp electrode recording from {cell_type} soma in the dorsolateral striatum - "
                    f"{condition} - Cell {session_info['cell_number']} from animal {session_info['animal_id']} - "
                    f"F-I protocol with {len(recording_folders)} current steps"
                ),
                "cell_id": f"Cell{session_info['cell_number']}",
                "location": "soma - dorsolateral striatum",
                "slice": general_metadata["NWBFile"]["slices"],
            }
        )

        # Update current clamp series metadata
        series_name = f"CurrentClampSeries{recording_info['protocol_step']}"
        interface_metadata["Icephys"]["CurrentClampSeries"][icephys_metadata_key].update(
            {
                "name": series_name,
                "description": (
                    f"Current clamp recording from {cell_type} - {condition} - "
                    f"Cell {session_info['cell_number']} from animal {session_info['animal_id']} - "
                    f"{recording_info['current_formatted']} current injection - "
                    f"F-I protocol step {recording_info['protocol_step']}"
                ),
            }
        )

        # Add intracellular data to NWB file
        interface.add_to_nwbfile(nwbfile=nwbfile, metadata=interface_metadata)

        # Add intracellular recording to icephys table with custom annotations
        current_clamp_series = nwbfile.acquisition[series_name]

        # Add intracellular recording entry with enhanced metadata
        recording_index = nwbfile.add_intracellular_recording(
            electrode=current_clamp_series.electrode,
            response=current_clamp_series,
            stimulus_current_pA=recording_info["current_pA"],
            protocol_step=recording_info["protocol_step"],
            recording_id=recording_id,
            animal_id=session_info["animal_id"],
        )

        # Track recording index for table building
        recording_indices.append(recording_index)

    # Build icephys table hierarchical structure using shared utility function
    build_somatic_icephys_table_structure(
        nwbfile=nwbfile,
        recording_indices=recording_indices,
        session_info=session_info,
        condition=condition,
        stimulus_type="F-I_protocol_M1R_CRISPR_somatic_excitability",
        include_animal_letter=True,
        animal_id_key="animal_id",
    )

    return nwbfile


if __name__ == "__main__":
    import argparse
    import logging
    import warnings

    from tqdm import tqdm

    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Convert Figure 8 somatic excitability data to NWB format")
    parser.add_argument(
        "--stub-test",
        type=bool,
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
    base_path = Path("./link_to_raw_data/Figure 8/M1R CRISPR SE")

    # Create nwb_files directory at root level
    root_dir = Path(__file__).parent.parent.parent.parent.parent.parent  # Go up to repo root
    nwb_files_dir = root_dir / "nwb_files" / "somatic_excitability" / "figure_8"
    nwb_files_dir.mkdir(parents=True, exist_ok=True)

    # Figure 8 M1R CRISPR somatic excitability conditions
    conditions = ["M1R CRISPR", "interleaved control"]

    for condition in conditions:
        condition_path = base_path / condition

        if not condition_path.exists():
            raise FileNotFoundError(f"Expected condition path does not exist: {condition_path}")

        # Get all session folders (each session = one cell from one animal on one date)
        session_folders = [f for f in condition_path.iterdir() if f.is_dir()]
        session_folders.sort()

        # Apply stub_test filtering if enabled
        if stub_test:
            session_folders = session_folders[:2]

        # Use tqdm for progress bar
        session_iterator = tqdm(
            session_folders,
            desc=f"Converting Figure8 SomaticExcitability {condition}",
            unit=" session",
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
