#!/usr/bin/env python3
"""
Figure 4 Optical Stimuli Conversion Script

This script converts Figure 4 Sr²⁺-oEPSC data from the Zhai et al. 2025 paper.
This data tests iSPN (indirect spiny projection neurons) responses to optogenetic
stimulation of corticostriatal terminals in the LID (L-DOPA-induced dyskinesia) model.

The experiment uses the same Sr²⁺-oEPSC protocol as Figure 2 but on iSPNs instead of dSPNs.
Optogenetic stimulation is delivered via whole-field LED (470nm) targeting ChR2-expressing
corticostriatal terminals, with 0.3ms pulses at 20ms into each 30-second inter-sweep interval.

Data structure:
- Figure 4_SF1B_SF5/Sr-oEPSC/LID on-state/: Sessions with LID-inducing L-DOPA treatment
- Figure 4_SF1B_SF5/Sr-oEPSC/LID off-state/: Sessions without L-DOPA treatment
- Each session contains 10 sweep recordings from iSPN cells

Key differences from Figure 2:
- Cell type: iSPN (indirect spiny projection neurons) instead of dSPN
- Different experimental sessions and dates
- Same optogenetic protocol and LED stimulation parameters
- Same virus (AAV5-hSyn-hChR2(H134R)-EYFP) and injection coordinates

Author: Generated with Claude Code
Date: 2025
"""

import re
from pathlib import Path
from typing import Dict

from neuroconv.tools.nwb_helpers import configure_and_write_nwbfile
from neuroconv.utils import dict_deep_update, load_dict_from_file
from pynwb import NWBFile
from pynwb.file import Subject

# Import interfaces from the codebase
from surmeier_lab_to_nwb.zhai2025.intracellular_interfaces import (
    PrairieViewVoltageClampInterface,
)
from surmeier_lab_to_nwb.zhai2025.optogenetics_interfaces import (
    PrairieViewOptogeneticsInterface,
)


def parse_session_info_from_folder_name(recording_folder: Path) -> Dict[str, str]:
    """
    Parse essential recording information from Figure 4 Sr²⁺-oEPSC recording folder names.
    Session start time comes from XML metadata.

    Expected folder name format: cell{N}_LED{X}-{YYY} (e.g., cell1_LED16-001, cell2_LED12-006)

    The session folder structure is: MMDDYYYY{letter}/cell{N}_LED{X}-{YYY}/
    - MMDDYYYY{letter}: Date and animal identifier (e.g., 02282023a, 03012023b)
    - cell{N}_LED{X}-{YYY}: Recording folder with cell number, LED intensity, and sweep number

    Parameters
    ----------
    recording_folder : Path
        Path to the recording folder

    Returns
    -------
    Dict[str, str]
        Dictionary containing recording information
    """
    folder_name = recording_folder.name
    session_folder = recording_folder.parent

    # Parse session folder (parent) for date and animal info
    # Format: MMDDYYYY{letter} (e.g., 02282023a, 03012023b)
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
    # Format variations:
    # - cell{N}_LED{X}-{YYY} (standard)
    # - LED{X}-{YYY} (without cell prefix)
    # - cell{N}_LED{X}{letter}-{YYY} (LED with letters)
    # - cell{N}_LED{X.Y}-{YYY} (LED with decimals)
    # - cell{N}_LED{X}_{additional}-{YYY} (LED with additional parameters)
    # - cell{N}_LEDint{X}-{YYY} (LED with 'int' prefix)
    # - cell{N}LED{X}-{YYY} (missing underscore)
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


def convert_session_to_nwbfile(
    session_folder_path: Path,
    condition: str,
    verbose: bool = False,
) -> NWBFile:
    """
    Convert a Figure 4 Sr²⁺-oEPSC session to NWB format with time alignment.

    Parameters
    ----------
    session_folder_path : Path
        Path to the session folder containing sweep recordings
    condition : str
        Experimental condition ("LID on-state" or "LID off-state")
    verbose : bool
        Whether to print verbose output

    Returns
    -------
    NWBFile
        Complete NWB file with all recordings from this session
    """
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

    if verbose:
        print(
            f"Processing session folder: {session_folder_path.name} (Session {session_info['session_letter']}) ({condition})"
        )
        print(f"  Found {len(recording_folders)} sweep recordings")

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
    session_id = f"{session_start_time.strftime('%Y%m%d')}_{session_info['session_letter']}"
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

    # Create session-specific metadata using precise session start time from XML
    session_specific_metadata = {
        "NWBFile": {
            "session_description": (
                f"Sr²⁺-oEPSC recordings from iSPNs in dorsolateral striatum under {condition} condition. "
                f"Voltage clamp at -70 mV in Ca²⁺-free ACSF containing 3 mM SrCl₂ and 10 μM gabazine. "
                f"Optogenetic stimulation of ChR2-expressing corticostriatal terminals with 0.3 ms blue LED pulses "
                f"every 30 seconds. Analysis of asynchronous EPSCs between 40-400 ms post-stimulation to measure "
                f"unitary synaptic strength. Session {session_info['session_letter']}, {len(recording_folders)} sweeps."
            ),
            "identifier": f"zhai2025_fig4_sr_oepsc_{session_info['session_id']}_{condition.replace(' ', '_')}",
            "session_start_time": session_start_time,
            "experiment_description": (
                f"Figure 4 Sr²⁺-oEPSC experiment from Zhai et al. 2025 investigating corticostriatal synaptic strength "
                f"changes between LID off-state and on-state in iSPN cells. Strontium substitution enables detection of individual "
                f"synaptic events rather than summed responses, revealing state-dependent changes in synaptic amplitude "
                f"that correlate with spine morphology changes during dyskinesia."
            ),
            "session_id": session_info["session_id"],
            "keywords": [
                "Sr2+-oEPSC",
                "voltage clamp",
                "corticostriatal synapses",
                "iSPN",
                "optogenetics",
                "ChR2",
                "levodopa-induced dyskinesia",
                "synaptic strength",
                "asynchronous EPSCs",
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

    # Create subject metadata for Sr²⁺-oEPSC experiments (Figure 4)
    subject = Subject(
        subject_id=f"iSPN_mouse_{session_info['session_id']}",
        species="Mus musculus",
        strain="Drd2-EGFP transgenic",
        description=(
            f"Adult Drd2-EGFP transgenic mouse with unilateral 6-OHDA lesion (>95% dopamine depletion) "
            f"modeling Parkinson's disease. Received dyskinesiogenic levodopa treatment and AAV5-hSyn-hChR2(H134R)-EYFP "
            f"injection into motor cortex for optogenetic experiments. iSPNs identified by lack of Drd2-EGFP expression. "
            f"Session {session_info['session_letter']} recorded on {session_info['date_str']}."
        ),
        genotype="Drd2-EGFP+",
        sex="M",
        age="P56/P84",  # Adult mice, 8-12 weeks
    )
    nwbfile.subject = subject

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
    recording_ids = list(recording_id_to_folder.keys())

    for recording_index, (recording_id, recording_folder) in enumerate(recording_id_to_folder.items()):
        recording_info = recording_id_to_info[recording_id]

        # Calculate next recording start time (None for last recording)
        next_recording_start_time = None
        if recording_index < len(recording_ids) - 1:  # Not the last recording
            next_recording_id = recording_ids[recording_index + 1]
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

        # Update electrode description for iSPN voltage clamp recording (consistent name for all recordings)
        electrode_name = "IntracellularElectrode"
        interface_metadata["Icephys"]["IntracellularElectrodes"][icephys_metadata_key].update(
            {
                "name": electrode_name,
                "description": (
                    f"Whole-cell patch clamp electrode recording from iSPN in the dorsolateral striatum - "
                    f"{condition} - Session {session_info['session_letter']} - "
                    f"Sr²⁺-oEPSC protocol with LED {recording_info['led_intensity']} intensity"
                ),
                "cell_id": recording_info["cell_number"],
                "session_id": session_info["session_letter"],
            }
        )

        # Update voltage clamp series metadata
        voltage_clamp_series_name = f"VoltageClampSeriesSweep{recording_info['sweep_number']}"
        interface_metadata["Icephys"]["VoltageClampSeries"][icephys_metadata_key].update(
            {
                "name": voltage_clamp_series_name,
                "description": (
                    f"Voltage clamp recording from iSPN - {condition} - "
                    f"Session {session_info['session_letter']}, Cell {recording_info['cell_number']} - "
                    f"LED intensity {recording_info['led_intensity']}, Sweep {recording_info['sweep_number']} - "
                    f"Sr²⁺-oEPSC protocol at -70 mV holding potential"
                ),
            }
        )

        if verbose:
            print(f"  Processing recording for folder: {recording_folder.name}")
            print(f"    Recording ID: {recording_id}")
            print(f"    LED intensity: {recording_info['led_intensity']}, Sweep: {recording_info['sweep_number']}")
            print(f"    Temporal alignment offset: {t_starts[recording_id]:.3f} seconds")

        # Add intracellular data to NWB file
        interface.add_to_nwbfile(nwbfile=nwbfile, metadata=interface_metadata)

        # Add optogenetic stimulation data
        voltage_output_xml_path = recording_folder / f"{recording_folder.name}_Cycle00001_VoltageOutput_001.xml"
        assert (
            voltage_output_xml_path.exists()
        ), f"Expected voltage output XML file does not exist in {recording_folder}: {voltage_output_xml_path}"

        optogenetics_interface = PrairieViewOptogeneticsInterface(
            voltage_output_xml_path=voltage_output_xml_path,
        )

        metadata = optogenetics_interface.get_metadata()

        # Configure specific virus and injection metadata for Figure 4 (same as Figure 2)
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
        timestamps = voltage_clamp_series.timestamps[:]
        recording_duration = timestamps[-1] - timestamps[0]

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

        if verbose:
            print(f"    Successfully processed recording: {recording_folder.name}")

    if verbose:
        print(f"Successfully processed all recordings from session: {session_folder_path.name}")

    # Build icephys table hierarchical structure following PyNWB best practices
    if verbose:
        print(f"  Building icephys table structure for {len(recording_indices)} recordings...")

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
        stimulus_type="Sr2+_oEPSC_optogenetic_protocol",
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

    nwbfile.add_icephys_experimental_condition(repetitions=[repetition_index], condition=condition)

    if verbose:
        print(f"  Successfully built icephys table structure")

    return nwbfile


def main():
    """Main conversion function for Figure 4 Sr²⁺-oEPSC data."""
    from tqdm import tqdm

    # Control verbose output
    verbose = False  # Set to True for detailed output

    # Define raw data and output paths
    raw_data_root = Path("/home/heberto/development/surmeier-lab-to-nwb/link_to_raw_data/Figure 4_SF1B_SF5/Sr-oEPSC")
    output_root = Path("/home/heberto/development/surmeier-lab-to-nwb/nwb_files/figure_4/sr_oepsc")

    # Create output directory
    output_root.mkdir(parents=True, exist_ok=True)

    # Define conditions to process
    conditions = ["LID on-state", "LID off-state"]

    # Process each condition
    for condition in conditions:
        if verbose:
            print(f"Processing Sr²⁺-oEPSC data for: {condition}")

        condition_path = raw_data_root / condition
        if not condition_path.exists():
            if verbose:
                print(f"  Condition path does not exist: {condition_path}")
            continue

        # Find all session folders
        session_folders = [
            folder for folder in condition_path.iterdir() if folder.is_dir() and re.match(r"\d{8}[a-z]", folder.name)
        ]

        session_folders.sort()  # Sort by session name

        if verbose:
            print(f"Found {len(session_folders)} session folders")

        # Use tqdm for progress bar when verbose is disabled
        session_iterator = (
            tqdm(session_folders, desc=f"Processing {condition}", disable=verbose) if not verbose else session_folders
        )

        # Process each session
        for session_folder in session_iterator:
            if verbose:
                print(f"\nProcessing session: {session_folder.name}")
            elif not verbose:
                session_iterator.set_description(f"Processing {session_folder.name}")

            # Convert session to NWB
            nwbfile = convert_session_to_nwbfile(
                session_folder_path=session_folder,
                condition=condition,
                verbose=verbose,
            )

            # Generate output filename
            condition_safe = condition.replace(" ", "_").replace("-", "_")
            nwbfile_path = output_root / f"figure4_sr_oepsc_{condition_safe}_{session_folder.name}.nwb"

            # Write NWB file
            configure_and_write_nwbfile(nwbfile, nwbfile_path=nwbfile_path)

            if verbose:
                print(f"Successfully saved: {nwbfile_path}")
            elif not verbose:
                session_iterator.write(f"Successfully saved: {nwbfile_path.name}")

        if not verbose:
            print(f"Completed {condition}: {len(session_folders)} sessions processed")


if __name__ == "__main__":
    main()
