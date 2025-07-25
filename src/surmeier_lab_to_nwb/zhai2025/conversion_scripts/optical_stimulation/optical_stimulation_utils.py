"""
Utility functions for optical stimulation experiments.

This module contains shared utility functions used across multiple optical stimulation
conversion scripts to reduce code duplication and ensure consistency.
"""

import re
import uuid
from pathlib import Path
from typing import Any, Dict, List

from neuroconv.tools.nwb_helpers import make_nwbfile_from_metadata
from neuroconv.utils import dict_deep_update, load_dict_from_file
from pynwb import NWBFile

from surmeier_lab_to_nwb.zhai2025.conversion_scripts.conversion_utils import (
    format_condition,
    generate_canonical_session_id,
)
from surmeier_lab_to_nwb.zhai2025.interfaces import (
    PrairieViewOptogeneticsInterface,
    PrairieViewVoltageClampInterface,
)
from surmeier_lab_to_nwb.zhai2025.utils import time_series_duration


def parse_session_info_from_folder_name(recording_folder: Path) -> Dict[str, Any]:
    """
    Parse essential recording information from Sr²⁺-oEPSC recording folder names.
    Session start time comes from XML metadata.

    Expected folder name format: cell{N}_LED{X}-{YYY} (e.g., cell1_LED14-001, cell2_LED12-006)

    The session folder structure is: MMDDYYYY{letter}/cell{N}_LED{X}-{YYY}/
    - MMDDYYYY{letter}: Date and animal identifier (e.g., 07052023a, 07062023b)
    - cell{N}_LED{X}-{YYY}: Recording folder with cell number, LED intensity, and sweep number

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
    # Format: MMDDYYYY{letter} (e.g., 07052023a, 07062023b)
    session_folder_name = session_folder.name

    # Extract date and animal identifier
    # Pattern: 8 digits (MMDDYYYY) followed by a letter
    session_pattern = r"(\d{8})([a-z])"
    session_match = re.match(session_pattern, session_folder_name)

    if not session_match:
        raise ValueError(f"Could not parse session folder name: {session_folder_name}")

    session_letter = session_match.group(
        2
    )  # letter (a, b, c, etc.) - identifies different sessions within the same day

    # Note: Date will be extracted from XML session start time, not from folder name

    # Parse recording folder name for Sr²⁺-oEPSC protocol info
    # Format variations (covering both Figure 2 and Figure 4 patterns):
    # - cell{N}_LED{X}-{YYY} (standard)
    # - LED{X}-{YYY} (without cell prefix)
    # - cell{N}_LED{X}{letter}-{YYY} (LED with letters)
    # - cell{N}_LED{X.Y}-{YYY} (LED with decimals)
    # - cell{N}_LED{X}_{additional}-{YYY} (LED with additional parameters)
    # - cell{N}_LEDint{X}-{YYY} (LED with 'int' prefix - Figure 4)
    # - cell{N}LED{X}-{YYY} (missing underscore - Figure 4)
    pattern = r"(?:cell(\d+)_?)?LED(?:int)?([0-9.]+[a-z]*(?:_\w+)*)-(\d+)"
    match = re.match(pattern, folder_name)

    if not match:
        raise ValueError(f"Could not parse recording folder name: {folder_name}")

    cell_number = match.group(1) if match.group(1) else "1"  # Default to cell 1 if not specified
    led_intensity = match.group(2)
    sweep_number = match.group(3)

    # Create recording ID for unique identification
    recording_id = f"Cell{cell_number}_LED{led_intensity}_Sweep{sweep_number}"

    return {
        "cell_number": cell_number,
        "session_letter": session_letter,
        "led_intensity": led_intensity,
        "sweep_number": sweep_number,
        "recording_id": recording_id,
        "recording_folder_name": folder_name,
    }


def convert_optical_stimulation_session_to_nwbfile(
    session_folder_path: Path,
    condition: str,
    figure_config: Dict[str, Any],
    session_id_parameters: Dict[str, Any],
    verbose: bool = False,
) -> NWBFile:
    """
    Convert a single session of Sr²⁺-oEPSC optical stimulation data to NWB format.

    This shared function handles the common conversion logic for all optical stimulation
    experiments across different figures, with figure-specific variations handled through
    the configuration dictionary.

    Parameters
    ----------
    session_folder_path : Path
        Path to the session folder containing sweep recordings
    condition : str
        Experimental condition (e.g., "LID off-state", "LID on-state")
    figure_config : Dict[str, Any]
        Figure-specific configuration containing:
        - metadata_key: str, key for session_specific_metadata.yaml
        - cell_type: str, cell type for descriptions ("dSPN" or "iSPN")
        - subject_id_prefix: str, prefix for subject ID
    session_id_parameters : Dict[str, Any]
        Session ID parameters using revised schema (excluding timestamp which is determined here):
        - fig: str, figure number (e.g., "F2", "F4")
        - meas_comp: str, measurement + compartment (e.g., "oEPSC")
        - cell_type: str, cell type (e.g., "dSPN", "iSPN")
        - state: str, experimental state (e.g., "OffState", "OnState")
        - pharm: str, pharmacology condition (e.g., "none")
        - geno: str, genotype (e.g., "WT")
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
    # Extract configuration
    metadata_key = figure_config["metadata_key"]
    cell_type = figure_config["cell_type"]
    subject_id_prefix = figure_config["subject_id_prefix"]

    # Find all recording folders (sweeps) for this session
    recording_folders = [f for f in session_folder_path.iterdir() if f.is_dir()]
    recording_folders.sort()

    if not recording_folders:
        raise ValueError(f"No recording folders found in session folder: {session_folder_path}")

    # Parse session information from first recording folder (all should have same session info)
    first_recording_info = parse_session_info_from_folder_name(recording_folders[0])
    session_info = {
        "session_letter": first_recording_info["session_letter"],
    }

    # Calculate recording IDs, session start times, and create interface mappings
    session_start_times = []  # (timestamp, recording_folder, recording_id)
    recording_id_to_info = {}
    recording_id_to_folder = {}
    t_starts = {}  # t_starts[recording_id] = t_start_offset

    for recording_folder in recording_folders:
        # Parse recording information using unified function
        recording_info = parse_session_info_from_folder_name(recording_folder)

        recording_id = recording_info["recording_id"]

        # Find XML file for this recording
        xml_file = recording_folder / f"{recording_folder.name}_Cycle00001_VoltageRecording_001.xml"
        assert xml_file.exists(), f"Expected XML file does not exist: {xml_file}"

        # Get session start time from XML metadata
        session_start_time = PrairieViewVoltageClampInterface.get_session_start_time_from_file(xml_file)
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
    earliest_folder = next(folder for start_time, folder, _ in session_start_times if start_time == earliest_time)

    # Calculate t_start offsets for temporal alignment
    for start_time, folder, recording_id in session_start_times:
        # Calculate offset relative to overall session start time
        t_start_offset = (start_time - earliest_time).total_seconds()
        t_starts[recording_id] = t_start_offset

    # Use earliest time as session start time for NWB file
    session_start_time = earliest_time

    # Extract date from actual session start time and update session info
    session_date_str = session_start_time.strftime("%Y-%m-%d")
    timestamp = session_start_time.strftime("%Y%m%d%H%M%S")

    # Create canonical session ID with explicit parameters
    condition_human_readable = format_condition[condition]["human_readable"]

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

    session_info.update(
        {
            "date_str": session_date_str,
            "session_id": session_id,
        }
    )

    # Load general and session-specific metadata from YAML files
    general_metadata_path = Path(__file__).parent.parent.parent / "general_metadata.yaml"
    general_metadata = load_dict_from_file(general_metadata_path)

    session_metadata_path = Path(__file__).parent.parent.parent / "session_specific_metadata.yaml"
    session_metadata_template = load_dict_from_file(session_metadata_path)
    script_template = session_metadata_template[metadata_key]

    # Create session-specific metadata from template with runtime substitutions
    session_specific_metadata = {
        "NWBFile": {
            "session_description": script_template["NWBFile"]["session_description"].format(
                condition=condition_human_readable,
                session_letter=session_info["session_letter"],
                num_sweeps=len(recording_folders),
            ),
            "identifier": str(uuid.uuid4()),
            "session_start_time": session_start_time,
            "session_id": session_info["session_id"],
            "keywords": script_template["NWBFile"]["keywords"],
        },
        "Subject": {
            "subject_id": f"{subject_id_prefix}_{session_info['session_id']}",
            "description": script_template["Subject"]["description"].format(
                session_letter=session_info["session_letter"], date_str=session_info["date_str"]
            ),
            "genotype": script_template["Subject"]["genotype"],
        },
    }

    # Add surgery information for Figure 2 if specified
    if "surgery_addition" in script_template["NWBFile"]:
        session_specific_metadata["NWBFile"]["surgery"] = (
            general_metadata["NWBFile"]["surgery"] + " " + script_template["NWBFile"]["surgery_addition"]
        )

    # Deep merge with general metadata
    metadata = dict_deep_update(general_metadata, session_specific_metadata)

    # Create NWB file using neuroconv helper function
    nwbfile = make_nwbfile_from_metadata(metadata)

    # Add custom columns to intracellular recording table for Sr²⁺-oEPSC experiment annotations
    intracellular_recording_table = nwbfile.get_intracellular_recordings()
    intracellular_recording_table.add_column(
        name="led_intensity",
        description="LED intensity setting for optogenetic stimulation",
    )
    intracellular_recording_table.add_column(
        name="sweep_number", description="Sweep number within the experimental session"
    )
    intracellular_recording_table.add_column(
        name="cell_number", description="Cell number identifier for this recording session"
    )
    intracellular_recording_table.add_column(
        name="session_letter", description="Session identifier letter for this experimental day"
    )

    # Data structures for tracking icephys table indices
    recording_indices = []  # Store all intracellular recording indices
    recording_to_metadata = {}  # Map recording index to metadata for table building

    # Process each recording using the calculated recording IDs and temporal alignment
    # Sort recordings by start time to ensure chronological processing
    sorted_recording_items = sorted(
        recording_id_to_folder.items(), key=lambda item: t_starts[item[0]]  # Sort by t_starts value
    )

    for recording_index, (recording_id, recording_folder) in enumerate(sorted_recording_items):
        recording_info = recording_id_to_info[recording_id]

        # Calculate next recording start time (None for last recording)
        next_recording_start_time = None
        if recording_index < len(sorted_recording_items) - 1:  # Not the last recording
            next_recording_id = sorted_recording_items[recording_index + 1][0]  # Get recording_id from next item
            next_recording_start_time = t_starts[next_recording_id]

        # Find XML files
        xml_file = recording_folder / f"{recording_folder.name}_Cycle00001_VoltageRecording_001.xml"
        assert xml_file.exists(), f"Expected XML file does not exist: {xml_file}"

        # Create interface for this recording
        icephys_metadata_key = f"VoltageClamp{recording_id}"
        interface = PrairieViewVoltageClampInterface(file_path=xml_file, icephys_metadata_key=icephys_metadata_key)

        # Apply temporal alignment offset
        interface.set_aligned_starting_time(t_starts[recording_id])

        # Get and update interface metadata
        interface_metadata = interface.get_metadata()

        # Update electrode description for voltage clamp recording (consistent name for all recordings)
        electrode_name = "IntracellularElectrode"
        interface_metadata["Icephys"]["IntracellularElectrodes"][icephys_metadata_key].update(
            {
                "name": electrode_name,
                "description": (
                    f"Whole-cell patch clamp electrode recording from {cell_type} in the dorsolateral striatum - "
                    f"{condition} - Session {session_info['session_letter']} - "
                    f"Sr²⁺-oEPSC protocol with LED {recording_info['led_intensity']} intensity"
                ),
                "cell_id": f"Cell{recording_info['cell_number']}",
                "session_id": session_info["session_letter"],
                "location": "soma - dorsolateral striatum",
                "slice": general_metadata["NWBFile"]["slices"],
            }
        )

        # Update voltage clamp series metadata
        voltage_clamp_series_name = f"VoltageClampSeriesSweep{recording_info['sweep_number']}"
        interface_metadata["Icephys"]["VoltageClampSeries"][icephys_metadata_key].update(
            {
                "name": voltage_clamp_series_name,
                "description": (
                    f"Voltage clamp recording from {cell_type} - {condition} - "
                    f"Session {session_info['session_letter']}, Cell {recording_info['cell_number']} - "
                    f"LED intensity {recording_info['led_intensity']}, Sweep {recording_info['sweep_number']} - "
                    f"Sr²⁺-oEPSC protocol at -70 mV holding potential"
                ),
            }
        )

        # Add intracellular data to NWB file
        interface.add_to_nwbfile(nwbfile=nwbfile, metadata=interface_metadata)

        # Add optogenetic stimulus for this sweep
        voltage_output_files = list(recording_folder.glob("*_VoltageOutput_*.xml"))
        voltage_output_xml_path = voltage_output_files[0] if voltage_output_files else None
        assert (
            voltage_output_xml_path
        ), f"Expected voltage output XML file does not exist in {recording_folder}: {voltage_output_xml_path}"

        optogenetics_interface = PrairieViewOptogeneticsInterface(
            voltage_output_xml_path=voltage_output_xml_path,
        )

        metadata = optogenetics_interface.get_metadata()

        # Configure specific virus and injection metadata for Zhai2025 experiment
        metadata["OptogeneticExperimentMetadata"]["virus"].update(
            {
                "name": "AAV5-hSyn-hChR2-H134R-EYFP",
                "construct_name": "AAV5-hSyn-hChR2(H134R)-EYFP",
                "description": (
                    "Adeno-associated virus serotype 5 expressing channelrhodopsin-2 (H134R variant) "
                    "fused to enhanced yellow fluorescent protein under the human synapsin promoter. "
                    "Used for optogenetic activation of corticostriatal terminals. "
                    "Addgene #26973."
                ),
                "manufacturer": "Addgene",
                "titer_in_vg_per_ml": 1e13,  # Typical AAV titer for this construct
            }
        )

        metadata["OptogeneticExperimentMetadata"]["virus_injection"].update(
            {
                "name": "M1_cortex_injection",
                "description": (
                    "Stereotaxic injection of AAV5-hSyn-hChR2(H134R)-EYFP into M1 motor cortex "
                    "ipsilateral to 6-OHDA lesion for corticostriatal terminal labeling. "
                    "Coordinates: AP +1.15mm, ML -1.60mm, DV -1.55mm relative to Bregma. "
                    "Volume: 0.15 µL. Expression time: 4 weeks."
                ),
                "location": "M1 motor cortex",
                "hemisphere": "ipsilateral to 6-OHDA lesion",
                "reference": "Bregma at the cortical surface",
                "ap_in_mm": 1.15,  # mm anterior to Bregma
                "ml_in_mm": -1.60,  # mm medial to Bregma (negative = left)
                "dv_in_mm": -1.55,  # mm ventral to cortical surface (negative = deeper)
                "pitch_in_deg": 0.0,  # Vertical injection
                "yaw_in_deg": 0.0,  # No lateral angle
                "roll_in_deg": 0.0,  # No roll angle
                "volume_in_uL": 0.15,  # 0.15 µL injection volume
                "injection_date": None,  # Would need session-specific dates
            }
        )

        # Add optogenetic stimulus with 20ms delay (pulse occurs 20ms into the sweep)
        pulse_delay_seconds = 0.020

        # Get recording duration from voltage clamp series for proper epoch timing
        voltage_clamp_series = nwbfile.acquisition[voltage_clamp_series_name]
        recording_duration = time_series_duration(voltage_clamp_series)

        optogenetics_interface.add_to_nwbfile(
            nwbfile=nwbfile,
            metadata=metadata,
            starting_time=t_starts[recording_id] + pulse_delay_seconds,
            recording_duration=recording_duration,
            next_recording_start_time=next_recording_start_time,
        )

        # Add intracellular recording to icephys table with custom annotations

        # Add intracellular recording entry with enhanced metadata
        recording_index = nwbfile.add_intracellular_recording(
            electrode=voltage_clamp_series.electrode,
            response=voltage_clamp_series,
            led_intensity=recording_info["led_intensity"],
            sweep_number=recording_info["sweep_number"],
            cell_number=recording_info["cell_number"],
            session_letter=session_info["session_letter"],
        )

        # Track recording index and metadata for table building
        recording_indices.append(recording_index)
        recording_to_metadata[recording_index] = {
            "recording_id": recording_id,
            "recording_info": recording_info,
        }

    # Build icephys table hierarchical structure using shared utility function
    build_optical_icephys_table_structure(
        nwbfile=nwbfile,
        recording_indices=recording_indices,
        session_info=session_info,
        condition=condition,
        stimulus_type="Sr2+_oEPSC_optogenetic_protocol",
    )

    if verbose:
        print(f"Created NWB file with {len(recording_folders)} sweeps and session ID: {session_id}")

    return nwbfile


def build_optical_icephys_table_structure(
    nwbfile: NWBFile,
    recording_indices: List[int],
    session_info: Dict[str, Any],
    condition: str,
    stimulus_type: str = "Sr2+_oEPSC_optogenetic_protocol",
) -> int:
    """
    Build icephys table hierarchical structure following PyNWB best practices for optical stimulation experiments.

    This function creates the complete icephys table hierarchy:
    1. Simultaneous recordings (each sweep is its own simultaneous group)
    2. Sequential recordings (group all sweeps for this session)
    3. Repetitions table (for this session, it's one repetition)
    4. Experimental conditions table (group by LID condition)

    Parameters
    ----------
    nwbfile : NWBFile
        The NWB file to add the icephys table structure to
    recording_indices : List[int]
        List of intracellular recording indices from nwbfile.add_intracellular_recording calls
    session_info : Dict[str, Any]
        Session information dictionary containing at least 'session_letter'
    condition : str
        Experimental condition (e.g., "LID off-state", "LID on-state")
    stimulus_type : str, default="Sr2+_oEPSC_optogenetic_protocol"
        Type of stimulus protocol

    Returns
    -------
    int
        The experimental condition index from the final table level

    Notes
    -----
    This function assumes each sweep is recorded as a separate simultaneous group,
    which is typical for Sr²⁺-oEPSC protocols in optical stimulation experiments.
    """

    # Step 1: Build simultaneous recordings (each sweep is its own simultaneous group)
    simultaneous_recording_indices = []
    for recording_index in recording_indices:
        simultaneous_index = nwbfile.add_icephys_simultaneous_recording(
            recordings=[recording_index]  # Each sweep is its own simultaneous group
        )
        simultaneous_recording_indices.append(simultaneous_index)

    # Step 2: Build sequential recording (group all sweeps for this session)
    sequential_index = nwbfile.add_icephys_sequential_recording(
        simultaneous_recordings=simultaneous_recording_indices,
        stimulus_type=stimulus_type,
    )

    # Step 3: Build repetitions table (for this session, it's one repetition)
    repetitions_table = nwbfile.get_icephys_repetitions()
    repetitions_table.add_column(
        name="session_letter", description="Session identifier letter for this experimental day"
    )

    repetition_index = nwbfile.add_icephys_repetition(
        sequential_recordings=[sequential_index],
        session_letter=session_info["session_letter"],
    )

    # Step 4: Build experimental conditions table (group by LID condition)
    experimental_conditions_table = nwbfile.get_icephys_experimental_conditions()
    experimental_conditions_table.add_column(
        name="condition", description="Experimental condition for Sr²⁺-oEPSC study"
    )

    experimental_condition_index = nwbfile.add_icephys_experimental_condition(
        repetitions=[repetition_index], condition=condition
    )

    return experimental_condition_index
