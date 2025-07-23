import re
import uuid
from pathlib import Path
from typing import Any, Dict

from neuroconv.tools import configure_and_write_nwbfile
from neuroconv.utils import dict_deep_update, load_dict_from_file
from pynwb import NWBFile
from pynwb.file import Subject

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


def convert_session_to_nwbfile(session_folder_path: Path, condition: str, verbose: bool = False) -> NWBFile:
    """
    Convert a single session of Figure 8 M1R CRISPR somatic excitability data to NWB format with time alignment.

    Parameters
    ----------
    session_folder_path : Path
        Path to the session folder containing cell recordings
    condition : str
        Experimental condition (e.g., "M1R CRISPR", "interleaved control")
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
        "animal_id": first_recording_info["animal_id"],
        "date_part": first_recording_info["date_part"],
    }

    if verbose:
        print(
            f"Processing session folder: {session_folder_path.name} (Cell {session_info['cell_number']}, Animal {session_info['animal_id']})"
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

    if verbose:
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

    # Create session ID following new pattern
    base_session_id = f"figure8_SomaticExcitability_{condition.replace(' ', '_').replace('-', '_')}_{session_start_time.strftime('%Y%m%d_%H%M%S')}"
    script_specific_id = f"Cell{session_info['cell_number']}_Animal{session_info['animal_id']}"
    session_id = f"{base_session_id}_{script_specific_id}"

    session_info.update(
        {
            "date_str": session_date_str,
            "session_id": session_id,
        }
    )

    if verbose:
        print(f"Session date: {session_info['date_str']}")

    # Load metadata from YAML file
    metadata_file_path = Path(__file__).parent / "metadata.yaml"
    paper_metadata = load_dict_from_file(metadata_file_path)

    # Determine cell type and description based on condition
    if condition == "M1R CRISPR":
        cell_type = "iSPN"
        cell_description = (
            f"iSPNs with M1R CRISPR deletion identified by lack of Drd1-Tdtomato expression (negative selection) "
            f"and confirmed M1R protein loss through CRISPR-Cas9 gene editing. "
            f"Cell {session_info['cell_number']} from animal {session_info['animal_id']} recorded on {session_info['date_str']}."
        )
        genetic_manipulation = "M1R deletion via CRISPR-Cas9 (AAV-Cas9 + AAV-gRNA-FusionRed)"
    else:  # interleaved control
        cell_type = "iSPN"
        cell_description = (
            f"Control iSPNs identified by lack of Drd1-Tdtomato expression (negative selection). "
            f"Control group received gRNA-FusionRed only (no Cas9). "
            f"Cell {session_info['cell_number']} from animal {session_info['animal_id']} recorded on {session_info['date_str']}."
        )
        genetic_manipulation = "Control (AAV-gRNA-FusionRed only, no M1R deletion)"

    # Create session-specific metadata using precise session start time from XML
    session_specific_metadata = {
        "NWBFile": {
            "session_description": (
                f"M1R CRISPR somatic excitability assessment in {cell_type}s "
                f"for condition '{condition}'. Whole-cell patch clamp recording in current clamp mode "
                f"with current injection steps from -120 pA to +300 pA (500 ms duration each). "
                f"Cell {session_info['cell_number']} from animal {session_info['animal_id']} recorded on {session_info['date_str']}."
            ),
            "identifier": str(uuid.uuid4()),
            "session_start_time": session_start_time,
            "experiment_description": (
                f"M1R CRISPR effects on {cell_type} somatic excitability during condition '{condition}'. "
                f"This experiment is part of Figure 8 from Zhai et al. 2025, investigating how M1R deletion "
                f"affects {cell_type} excitability and dyskinetic behaviors using CRISPR-Cas9 gene editing. "
                f"F-I protocol with {len(recording_folders)} current steps."
            ),
            "session_id": session_info["session_id"],
            "keywords": [
                "somatic excitability",
                "F-I relationship",
                "rheobase",
                "M1R CRISPR",
                "CRISPR-Cas9",
                "gene editing",
            ],
        },
        "Subject": {
            "subject_id": f"{cell_type}_M1R_CRISPR_mouse_{session_info['session_id']}",
            "description": cell_description,
            "genotype": f"Drd1-Tdtomato bacterial artificial chromosome (BAC) transgenic; {genetic_manipulation}",
        },
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
                "cell_id": session_info["cell_number"],
                "location": "soma - dorsolateral striatum",
                "slice": "280 μm sagittal brain slice from dorsolateral striatum (Paper Methods: 'Sagittal sections (280 μm thick) were cut using a Leica VT1200 vibratome')",
                "seal": "Gigaohm seal (whole-cell configuration) (Paper Methods: patch clamp methodology, whole-cell configuration implied)",
                "resistance": "3-5 MΩ (borosilicate glass pipette) (Protocol: Ex_vivo_mouse_brain_patch_clamp_recordings: 'Pipette resistance must be of 3 to 5 megaohms')",
                "filtering": "2 kHz low-pass filter (Paper Methods: 'signals were filtered at 2 kHz and digitized at 10 kHz')",
                "initial_access_resistance": "<20 MΩ (typical for whole-cell recordings) (Standard electrophysiology practice for healthy whole-cell recordings)",
            }
        )

        # Update current clamp series metadata
        series_name = f"CurrentClamp{recording_info['protocol_step']}Series"
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
            animal_id=session_info["animal_id"],
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

    if verbose:
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
        stimulus_type="F-I_protocol_M1R_CRISPR_somatic_excitability",
    )
    sequential_recording_indices.append(sequential_index)

    # Step 3: Build repetitions table (for this single cell, it's just one repetition)
    repetitions_table = nwbfile.get_icephys_repetitions()
    repetitions_table.add_column(name="cell_number", description="Cell number identifier for this recording session")
    repetitions_table.add_column(
        name="animal_id", description="Animal identifier for tracking across experimental sessions"
    )

    repetition_index = nwbfile.add_icephys_repetition(
        sequential_recordings=sequential_recording_indices,
        cell_number=int(session_info["cell_number"]),
        animal_id=session_info["animal_id"],
    )

    # Step 4: Build experimental conditions table (group by M1R CRISPR condition)
    experimental_conditions_table = nwbfile.get_icephys_experimental_conditions()
    experimental_conditions_table.add_column(
        name="condition", description="M1R CRISPR experimental condition (M1R deletion vs. control)"
    )

    nwbfile.add_icephys_experimental_condition(repetitions=[repetition_index], condition=condition)

    if verbose:
        print(f"    Added experimental condition '{condition}' with 1 repetition")
        print(f"  Successfully built icephys table hierarchy:")
        print(f"    - {len(recording_indices)} intracellular recordings")
        print(f"    - {len(simultaneous_recording_indices)} simultaneous recordings")
        print(f"    - {len(sequential_recording_indices)} sequential recordings")
        print(f"    - 1 repetition (Cell {session_info['cell_number']}, Animal {session_info['animal_id']})")
        print(f"    - 1 experimental condition ('{condition}')")

    return nwbfile


if __name__ == "__main__":
    import logging
    import warnings

    from tqdm import tqdm

    # Control verbose output from here
    verbose = False  # Set to True for detailed output
    stub_test = True  # Set to True to process only first 2 files per condition for testing

    # Suppress tifffile warnings
    logging.getLogger("tifffile").setLevel(logging.ERROR)

    # Suppress specific warnings
    warnings.filterwarnings("ignore", message="invalid value encountered in divide")
    warnings.filterwarnings("ignore", message=".*no datetime before year 1.*")

    # Set the base path to your data
    base_path = Path("./link_to_raw_data/Figure 8/M1R CRISPR SE")

    # Create nwb_files directory at root level
    root_dir = Path(__file__).parent.parent.parent.parent  # Go up to repo root
    nwb_files_dir = root_dir / "nwb_files" / "figure_8" / "somatic_excitability"
    nwb_files_dir.mkdir(parents=True, exist_ok=True)

    # Figure 8 M1R CRISPR somatic excitability conditions
    conditions = ["M1R CRISPR", "interleaved control"]

    for condition in conditions:
        condition_path = base_path / condition

        if not condition_path.exists():
            raise FileNotFoundError(f"Expected condition path does not exist: {condition_path}")

        if verbose:
            print(f"Processing M1R CRISPR somatic excitability data for: {condition}")

        # Get all session folders (each session = one cell from one animal on one date)
        session_folders = [f for f in condition_path.iterdir() if f.is_dir()]
        session_folders.sort()

        # Apply stub_test filtering if enabled
        if stub_test:
            session_folders = session_folders[:2]
            if verbose:
                print(f"stub_test enabled: processing only first {len(session_folders)} session folders")

        if verbose:
            print(f"Found {len(session_folders)} session folders")

        # Use tqdm for progress bar when verbose is disabled
        session_iterator = tqdm(
            session_folders,
            desc=f"Converting {condition} from figure_8_somatic_experiments to NWB",
        )

        for session_folder in session_iterator:
            if verbose:
                print(f"\nProcessing session: {session_folder.name}")

            # Convert session data to NWB format with time alignment
            nwbfile = convert_session_to_nwbfile(
                session_folder_path=session_folder,
                condition=condition,
                verbose=verbose,
            )

            # Create output filename
            condition_safe = condition.replace(" ", "_").replace("-", "_").replace("(", "").replace(")", "")
            nwbfile_path = nwb_files_dir / f"figure8_somatic_excitability_{condition_safe}_{session_folder.name}.nwb"

            # Write NWB file
            configure_and_write_nwbfile(nwbfile, nwbfile_path=nwbfile_path)
            if verbose:
                print(f"Successfully saved: {nwbfile_path.name}")
