"""
Utility functions for somatic excitability experiments.

This module contains shared utility functions used across multiple somatic excitability
conversion scripts to reduce code duplication and ensure consistency.
"""

from pathlib import Path
from typing import Any, Dict

from neuroconv.tools.nwb_helpers import make_nwbfile_from_metadata
from neuroconv.utils import dict_deep_update, load_dict_from_file
from pynwb import NWBFile

from surmeier_lab_to_nwb.zhai2025.conversion_scripts.conversion_utils import (
    GENOTYPE_DESCRIPTION_MAPPING,
    PHARMACOLOGY_ADDITIONS,
    format_condition,
    generate_canonical_session_id,
)
from surmeier_lab_to_nwb.zhai2025.interfaces import PrairieViewCurrentClampInterface


def build_somatic_icephys_table_structure(
    nwbfile: NWBFile,
    recording_indices: list[int],
    condition: str,
) -> int:
    """
    Build icephys table hierarchical structure following PyNWB best practices for somatic excitability experiments.

    This function creates the complete icephys table hierarchy:
    1. Simultaneous recordings (each current step is its own simultaneous group)
    2. Sequential recordings (group all current steps for this session)
    3. Repetitions table (single repetition per session)
    4. Experimental conditions table (group by LID condition)

    Parameters
    ----------
    nwbfile : NWBFile
        The NWB file to add the icephys table structure to
    recording_indices : list[int]
        List of intracellular recording indices from nwbfile.add_intracellular_recording calls
    condition : str
        Experimental condition (e.g., "LID off-state", "LID on-state", etc.)

    Returns
    -------
    int
        The experimental condition index from the final table level

    Notes
    -----
    This function assumes each current injection step is recorded as a separate simultaneous group,
    which is typical for F-I relationship protocols in somatic excitability experiments.
    """

    # Step 1: Build simultaneous recordings (each current step is its own simultaneous group)
    simultaneous_recording_indices = []
    for recording_index in recording_indices:
        simultaneous_index = nwbfile.add_icephys_simultaneous_recording(
            recordings=[recording_index]  # Each current step is its own simultaneous group
        )
        simultaneous_recording_indices.append(simultaneous_index)

    # Step 2: Build sequential recording (group all current steps for this session)
    sequential_index = nwbfile.add_icephys_sequential_recording(
        simultaneous_recordings=simultaneous_recording_indices,
        stimulus_type="F-I_protocol_somatic_excitability",
    )
    sequential_recording_indices = [sequential_index]

    # Step 3: Build repetitions table (single repetition per session)
    repetition_index = nwbfile.add_icephys_repetition(sequential_recordings=sequential_recording_indices)

    # Step 4: Build experimental conditions table (group by experimental condition)
    experimental_conditions_table = nwbfile.get_icephys_experimental_conditions()
    experimental_conditions_table.add_column(
        name="condition", description="Experimental condition for L-DOPA induced dyskinesia study"
    )

    experimental_condition_index = nwbfile.add_icephys_experimental_condition(
        repetitions=[repetition_index], condition=condition
    )

    return experimental_condition_index


def convert_somatic_excitability_session_to_nwbfile(
    session_folder_path: Path,
    condition: str,
    figure_config: Dict[str, Any],
    session_id_parameters: Dict[str, Any],
) -> NWBFile:
    """
    Convert a single session of somatic excitability data to NWB format with time alignment.

    This shared function handles the common conversion logic for all somatic excitability
    experiments across different figures, with figure-specific variations handled through
    the configuration dictionary.

    Parameters
    ----------
    session_folder_path : Path
        Path to the session folder containing cell recordings
    condition : str
        Experimental condition (e.g., "LID off-state", "control", "KO on-state", etc.)
    figure_config : Dict[str, Any]
        Figure-specific configuration containing:
        - parse_function: Callable to parse folder names
        - metadata_key: str, key for session_specific_metadata.yaml
    session_id_parameters : Dict[str, Any]
        Session ID parameters using revised schema (excluding timestamp which is determined here):
        - fig: str, figure number (e.g., "F1", "F3")
        - meas_comp: str, measurement + compartment (e.g., "SomExc", "DendExc")
        - cell_type: str, cell type (e.g., "dSPN", "iSPN")
        - state: str, experimental state (e.g., "OffState", "OnState")
        - pharm: str, pharmacology condition (e.g., "none", "D1RaSch")
        - geno: str, genotype (e.g., "WT", "CDGIKO")

    Returns
    -------
    NWBFile
        NWB file with the converted data

    Notes
    -----
    This function implements temporal alignment by extracting precise timestamps from XML files
    and calculating t_start offsets for each recording relative to the earliest recording.
    """
    # Extract configuration
    parse_function = figure_config["parse_function"]
    metadata_key = figure_config["metadata_key"]

    # Get cell type from session ID parameters
    cell_type = session_id_parameters["cell_type"]

    # Find all recording folders (current steps) for this cell
    recording_folders = [f for f in session_folder_path.iterdir() if f.is_dir()]
    recording_folders.sort()

    if not recording_folders:
        raise ValueError(f"No recording folders found in session folder: {session_folder_path}")

    # Parse session information from first recording folder (unused but kept for consistency)

    # Calculate recording IDs, session start times, and create interface mappings
    session_start_times = []  # (timestamp, recording_folder, recording_id)
    recording_id_to_info = {}
    recording_id_to_folder = {}
    t_starts = {}  # t_starts[recording_id] = t_start_offset

    for recording_folder in recording_folders:
        # Parse recording information using figure-specific function
        recording_info = parse_function(recording_folder)

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

    # Create canonical session ID with explicit parameters
    timestamp = session_start_time.strftime("%Y%m%d%H%M%S")

    # Use session ID parameters passed by the calling script (revised schema)
    session_id = generate_canonical_session_id(
        fig=session_id_parameters["fig"],
        meas_comp=session_id_parameters["meas_comp"],
        cell_type=session_id_parameters["cell_type"],
        state=session_id_parameters["state"],
        pharm=session_id_parameters["pharm"],
        geno=session_id_parameters["geno"],
        timestamp=timestamp,
    )

    # Load general and session-specific metadata from YAML files
    general_metadata_path = Path(__file__).parent.parent.parent / "general_metadata.yaml"
    general_metadata = load_dict_from_file(general_metadata_path)

    session_metadata_path = Path(__file__).parent.parent.parent / "session_specific_metadata.yaml"
    session_metadata_template = load_dict_from_file(session_metadata_path)
    script_template = session_metadata_template[metadata_key]

    # Get condition formatting for underscore version
    condition_underscore = format_condition[condition]["underscore"]

    # Build pharmacology description using centralized mapping
    pharmacology_base = general_metadata["NWBFile"]["pharmacology"]
    pharm_token = session_id_parameters["pharm"]

    # Add pharmacology-specific text if applicable
    if pharm_token != "none" and pharm_token in PHARMACOLOGY_ADDITIONS:
        pharmacology_text = pharmacology_base + " " + PHARMACOLOGY_ADDITIONS[pharm_token]
    else:
        pharmacology_text = pharmacology_base

    # Handle special metadata fields (like surgery for Figure 8)
    surgery_text = general_metadata["NWBFile"]["surgery"]
    if "surgery_addition" in script_template["NWBFile"]:
        surgery_text += " " + script_template["NWBFile"]["surgery_addition"]

    # Create session-specific metadata from template with runtime substitutions
    condition_human_readable = format_condition[condition]["human_readable"]
    session_specific_metadata = {
        "NWBFile": {
            "session_description": script_template["NWBFile"]["session_description"].format(
                condition=condition_human_readable
            ),
            "session_start_time": session_start_time,
            "session_id": session_id,
            "pharmacology": pharmacology_text,
            "keywords": script_template["NWBFile"]["keywords"],
        },
        "Subject": {
            "subject_id": f"SubjectRecordedAt{timestamp}",
            "description": script_template["Subject"]["description"],
            "genotype": script_template["Subject"]["genotype"],
        },
    }

    # Add surgery field if it was modified
    if surgery_text != general_metadata["NWBFile"]["surgery"]:
        session_specific_metadata["NWBFile"]["surgery"] = surgery_text

    # Merge general metadata with session-specific metadata
    metadata = dict_deep_update(general_metadata, session_specific_metadata)

    # Create NWB file using neuroconv helper function
    nwbfile = make_nwbfile_from_metadata(metadata)

    # Add custom columns to intracellular recording table for somatic experiment annotations
    intracellular_recording_table = nwbfile.get_intracellular_recordings()
    intracellular_recording_table.add_column(
        name="stimulus_current_pA",
        description="The current stimulus applied during the recording in picoamps",
    )
    intracellular_recording_table.add_column(
        name="protocol_step", description="Protocol step number (e.g., '001', '002', etc.)"
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

        # Build description based on figure-specific context using centralized mapping
        genotype = session_id_parameters["geno"]
        genotype_description = GENOTYPE_DESCRIPTION_MAPPING.get(genotype, "")

        # Update electrode description for somatic recording (consistent name for all recordings)
        electrode_name = "IntracellularElectrode"
        interface_metadata["Icephys"]["IntracellularElectrodes"][icephys_metadata_key].update(
            {
                "name": electrode_name,
                "description": (
                    f"Whole-cell patch clamp electrode recording from {genotype_description}{cell_type} soma in the dorsolateral striatum - "
                    f"{condition_human_readable} - F-I protocol with {len(recording_folders)} current steps"
                ),
                "cell_id": f"CellRecordedAt{timestamp}",
                "location": "soma - dorsolateral striatum",
                "slice": general_metadata["NWBFile"]["slices"],
            }
        )

        # Update current clamp series metadata
        series_name = f"CurrentClampSeries{recording_info['protocol_step']}"
        series_description_prefix = f"{genotype_description}{cell_type}" if genotype != "WT" else cell_type
        interface_metadata["Icephys"]["CurrentClampSeries"][icephys_metadata_key].update(
            {
                "name": series_name,
                "description": (
                    f"Current clamp recording from {series_description_prefix} - {condition_human_readable} - "
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
        )

        # Track recording index for table building
        recording_indices.append(recording_index)

    # Build icephys table hierarchical structure using shared utility function
    build_somatic_icephys_table_structure(
        nwbfile=nwbfile,
        recording_indices=recording_indices,
        condition=condition_underscore,
    )

    return nwbfile
