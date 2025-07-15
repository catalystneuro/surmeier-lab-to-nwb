import re
from pathlib import Path
from typing import Any, Dict

from neuroconv.tools import configure_and_write_nwbfile
from neuroconv.utils import dict_deep_update, load_dict_from_file
from pynwb import NWBFile
from pynwb.file import Subject

from surmeier_lab_to_nwb.zhai2025.intracellular_interfaces import (
    PROTOCOL_STEP_TO_CURRENT,
    PrairieViewCurrentClampInterface,
)


def parse_session_info_from_folder_name(recording_folder: Path) -> Dict[str, Any]:
    """
    Parse essential recording information from CDGI knockout somatic excitability recording folder names.
    Session start time comes from XML metadata.

    Expected folder name format: cell[N]-[XXX] (e.g., cell1-001, cell2-015)

    The session folder structure is: [MMDD][letter]/cell[N]-[XXX]/
    - [MMDD][letter]: Date and animal identifier (e.g., 0525a, 0525b)
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

    # Parse session folder (parent) for date and animal info
    # Format: [MMDD][letter] (e.g., 0525a, 0816b)
    session_folder_name = session_folder.name

    # Extract date and animal identifier
    # Pattern: 4 digits (MMDD) followed by a letter
    session_pattern = r"(\d{4})([a-z])"
    session_match = re.match(session_pattern, session_folder_name)

    if not session_match:
        raise ValueError(f"Could not parse session folder name: {session_folder_name}")

    animal_letter = session_match.group(2)  # letter (a, b, c, etc.)

    # Note: Date will be extracted from XML session start time, not from folder name

    # Parse recording folder name for protocol step
    # Format: cell[N]-[XXX] (no suffixes for Figure 7)
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
        "animal_letter": animal_letter,
        "protocol_step": protocol_step,
        "current_pA": current_pA,
        "current_formatted": current_formatted,
        "step_order": step_order,
        "recording_folder_name": folder_name,
    }


def convert_session_to_nwbfile(session_folder_path: Path, condition: str, verbose: bool = False) -> NWBFile:
    """
    Convert a single session of Figure 7 CDGI knockout somatic excitability data to NWB format with time alignment.

    Parameters
    ----------
    session_folder_path : Path
        Path to the session folder containing cell recordings
    condition : str
        Experimental condition (e.g., "KO off-state", "KO on-state")
    verbose : bool, default=False
        Enable verbose output showing detailed processing information

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
        "animal_letter": first_recording_info["animal_letter"],
    }

    print(
        f"Processing session folder: {session_folder_path.name} (Animal {session_info['animal_letter']}, Cell {session_info['cell_number']})"
    )
    print(f"  Found {len(recording_folders)} current step recordings")

    # Calculate recording IDs, session start times, and create interface mappings
    session_start_times = []  # (timestamp, recording_folder, recording_id)
    recording_id_to_info = {}
    recording_id_to_folder = {}
    t_starts = {}  # t_starts[recording_id] = t_start_offset

    if verbose:
        print(f"  Validating session start times and calculating recording IDs...")

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

        if verbose:
            print(f"    Recording {recording_folder.name}: timestamp = {session_start_time}")

    if not session_start_times:
        raise ValueError(f"No valid recordings found in session folder: {session_folder_path}")

    # Find the earliest session start time across all recordings
    earliest_time = min(session_start_times, key=lambda x: x[0])[0]
    earliest_folder = next(folder for start_time, folder, _ in session_start_times if start_time == earliest_time)

    print(f"  Overall session start time: {earliest_time}")
    print(f"    Earliest time source: recording {earliest_folder.name}")

    # Calculate t_start offsets for temporal alignment
    for start_time, folder, recording_id in session_start_times:
        # Calculate offset relative to overall session start time
        t_start_offset = (start_time - earliest_time).total_seconds()
        t_starts[recording_id] = t_start_offset

        if verbose:
            print(f"    Recording {folder.name} ({recording_id}) temporal alignment:")
            print(f"      t_start offset = {t_start_offset:.3f} seconds")

    # Use earliest time as session start time for NWB file
    session_start_time = earliest_time

    # Extract date from actual session start time and update session info
    session_date_str = session_start_time.strftime("%Y-%m-%d")
    session_id = (
        f"{session_start_time.strftime('%Y%m%d')}_{session_info['animal_letter']}_Cell{session_info['cell_number']}"
    )
    session_info.update(
        {
            "date_str": session_date_str,
            "session_id": session_id,
        }
    )

    print(f"Session date: {session_info['date_str']}")

    # Load metadata from YAML file
    metadata_file_path = Path(__file__).parent / "metadata.yaml"
    paper_metadata = load_dict_from_file(metadata_file_path)

    # Create session-specific metadata using precise session start time from XML

    session_specific_metadata = {
        "NWBFile": {
            "session_description": (
                f"Somatic excitability assessment in CDGI knockout iSPNs for condition '{condition}'. "
                f"Whole-cell patch clamp recording in current clamp mode with current injection steps "
                f"from -120 pA to +300 pA (500 ms duration each). Animal {session_info['animal_letter']}, "
                f"Cell {session_info['cell_number']} recorded on {session_info['date_str']}."
            ),
            "identifier": f"zhai2025_fig7_somatic_{session_info['session_id']}_{condition.replace(' ', '_')}",
            "session_start_time": session_start_time,
            "experiment_description": (
                f"CDGI knockout effects on iSPN somatic excitability during condition '{condition}'. "
                f"This experiment is part of Figure 7 from Zhai et al. 2025, investigating the role of "
                f"CalDAG-GEFI in iSPN excitability changes during LID. F-I protocol with {len(recording_folders)} current steps."
            ),
            "session_id": session_info["session_id"],
            "keywords": [
                "somatic excitability",
                "F-I relationship",
                "rheobase",
                "CDGI knockout",
                "CalDAG-GEFI",
            ],
        }
    }

    # Deep merge with paper metadata
    metadata = dict_deep_update(paper_metadata, session_specific_metadata)

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
        keywords=metadata["NWBFile"]["keywords"],
    )

    # Create subject metadata for CDGI knockout experiments (Figure 7)
    subject = Subject(
        subject_id=f"CDGI_KO_mouse_{session_info['session_id']}",
        species="Mus musculus",
        strain="C57Bl/6",
        description=(
            f"CDGI knockout mouse with unilateral 6-OHDA lesion in the medial forebrain bundle. "
            f"iSPNs identified by lack of Drd1-Tdtomato expression (negative selection). "
            f"Animal {session_info['animal_letter']}, Cell {session_info['cell_number']} recorded on {session_info['date_str']}."
        ),
        genotype="CalDAG-GEFI knockout on Drd1-Tdtomato bacterial artificial chromosome (BAC) transgenic background",
        sex="M",
        age="P49-P84",  # ISO format for 7-12 weeks (postnatal days)
    )
    nwbfile.subject = subject

    # Add custom columns to intracellular recording table for CDGI knockout experiment annotations
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
        name="animal_letter", description="Animal identifier letter for this experimental session"
    )

    # Data structures for tracking icephys table indices
    recording_indices = []  # Store all intracellular recording indices
    recording_to_metadata = {}  # Map recording index to metadata for table building
    sequential_recording_indices = []  # Store sequential recording indices

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

        # Update electrode description for CDGI knockout iSPN somatic recording (consistent name for all recordings)
        electrode_name = "IntracellularElectrode"
        interface_metadata["Icephys"]["IntracellularElectrodes"][icephys_metadata_key].update(
            {
                "name": electrode_name,
                "description": (
                    f"Whole-cell patch clamp electrode recording from CDGI knockout iSPN soma in the dorsolateral striatum - "
                    f"{condition} - Animal {session_info['animal_letter']}, Cell {session_info['cell_number']} - "
                    f"F-I protocol with {len(recording_folders)} current steps"
                ),
                "cell_id": session_info["cell_number"],
                "animal_id": session_info["animal_letter"],
                "location": "soma",
            }
        )

        # Update current clamp series metadata
        series_name = f"CurrentClamp{recording_info['protocol_step']}Series"
        interface_metadata["Icephys"]["CurrentClampSeries"][icephys_metadata_key].update(
            {
                "name": series_name,
                "description": (
                    f"Current clamp recording from CDGI knockout iSPN - {condition} - "
                    f"Animal {session_info['animal_letter']}, Cell {session_info['cell_number']} - "
                    f"{recording_info['current_formatted']} current injection - "
                    f"F-I protocol step {recording_info['protocol_step']}"
                ),
            }
        )

        if verbose:
            print(f"  Processing recording for folder: {recording_folder.name}")
            print(f"    Recording ID: {recording_id}")
            print(f"    Current: {recording_info['current_pA']} pA (step {recording_info['protocol_step']})")
            print(f"    Temporal alignment offset: {t_starts[recording_id]:.3f} seconds")

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
            animal_letter=session_info["animal_letter"],
        )

        # Track recording index and metadata for table building
        recording_indices.append(recording_index)
        recording_to_metadata[recording_index] = {
            "recording_id": recording_id,
            "recording_info": recording_info,
            "series_name": series_name,
        }

        if verbose:
            print(f"    Successfully processed recording: {recording_folder.name}")

    print(f"Successfully processed all recordings from session: {session_folder_path.name}")

    # Build icephys table hierarchical structure following PyNWB best practices
    if verbose:
        print(f"  Building icephys table structure for {len(recording_indices)} recordings...")

    # Step 1: Build simultaneous recordings (each current step is its own simultaneous group)
    simultaneous_recording_indices = []
    for recording_index in recording_indices:
        simultaneous_index = nwbfile.add_icephys_simultaneous_recording(
            recordings=[recording_index]  # Each current step is its own simultaneous group
        )
        simultaneous_recording_indices.append(simultaneous_index)

    # Step 2: Build sequential recording (group all current steps for this cell)
    sequential_index = nwbfile.add_icephys_sequential_recording(
        simultaneous_recordings=simultaneous_recording_indices,
        stimulus_type="F-I_protocol_somatic_excitability_CDGI_knockout",
    )
    sequential_recording_indices.append(sequential_index)

    # Step 3: Build repetitions table (for this single cell, it's just one repetition)
    repetitions_table = nwbfile.get_icephys_repetitions()
    repetitions_table.add_column(name="cell_number", description="Cell number identifier for this recording session")
    repetitions_table.add_column(
        name="animal_letter", description="Animal identifier letter for this experimental session"
    )

    repetition_index = nwbfile.add_icephys_repetition(
        sequential_recordings=sequential_recording_indices,
        cell_number=int(session_info["cell_number"]),
        animal_letter=session_info["animal_letter"],
    )

    # Step 4: Build experimental conditions table (group by CDGI knockout condition)
    experimental_conditions_table = nwbfile.get_icephys_experimental_conditions()
    experimental_conditions_table.add_column(
        name="condition", description="Experimental condition for CDGI knockout LID study"
    )

    experimental_condition_index = nwbfile.add_icephys_experimental_condition(
        repetitions=[repetition_index], condition=condition
    )

    if verbose:
        print(f"    Added experimental condition '{condition}' with 1 repetition")
        print(f"  Successfully built icephys table hierarchy:")
        print(f"    - {len(recording_indices)} intracellular recordings")
        print(f"    - {len(simultaneous_recording_indices)} simultaneous recordings")
        print(f"    - {len(sequential_recording_indices)} sequential recordings")
        print(f"    - 1 repetition (Animal {session_info['animal_letter']}, Cell {session_info['cell_number']})")
        print(f"    - 1 experimental condition ('{condition}')")

    return nwbfile


if __name__ == "__main__":
    import logging
    import warnings

    # Control verbose output from here
    VERBOSE = False  # Set to True for detailed output

    # Suppress tifffile warnings
    logging.getLogger("tifffile").setLevel(logging.ERROR)

    # Suppress specific warnings
    warnings.filterwarnings("ignore", message="invalid value encountered in divide")
    warnings.filterwarnings("ignore", message=".*no datetime before year 1.*")

    # Set the base path to your data
    base_path = Path("./link_to_raw_data/Figure 7/KO SE on vs off")

    # Create nwb_files directory at root level
    root_dir = Path(__file__).parent.parent.parent.parent  # Go up to repo root
    nwb_files_dir = root_dir / "nwb_files" / "figure_7_somatic_excitability"
    nwb_files_dir.mkdir(parents=True, exist_ok=True)

    # Figure 7 CDGI knockout somatic excitability conditions
    conditions = ["KO off-state", "KO on-state"]

    for condition in conditions:
        condition_path = base_path / condition

        if not condition_path.exists():
            raise FileNotFoundError(f"Expected condition path does not exist: {condition_path}")

        print(f"Processing CDGI knockout somatic excitability data for: {condition}")

        # Get all session folders (each session = one animal/cell)
        session_folders = [f for f in condition_path.iterdir() if f.is_dir()]
        session_folders.sort()

        print(f"Found {len(session_folders)} session folders")

        for session_folder in session_folders:
            print(f"\nProcessing session: {session_folder.name}")

            # Convert session data to NWB format with time alignment
            nwbfile = convert_session_to_nwbfile(
                session_folder_path=session_folder,
                condition=condition,
                verbose=VERBOSE,
            )

            # Create output filename
            condition_safe = condition.replace(" ", "_").replace("-", "_")
            nwbfile_path = nwb_files_dir / f"figure7_somatic_excitability_{condition_safe}_{session_folder.name}.nwb"

            # Write NWB file
            configure_and_write_nwbfile(nwbfile, nwbfile_path=nwbfile_path)
            print(f"Successfully saved: {nwbfile_path.name}")
