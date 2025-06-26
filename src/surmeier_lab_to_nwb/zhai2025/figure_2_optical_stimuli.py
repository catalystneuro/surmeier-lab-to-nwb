# -*- coding: utf-8 -*-
import re
from datetime import datetime
from pathlib import Path
from typing import Any, Dict
from zoneinfo import ZoneInfo

from neuroconv.tools import configure_and_write_nwbfile
from pynwb import NWBFile
from pynwb.file import Subject

from surmeier_lab_to_nwb.zhai2025.intracellular_interfaces import (
    PrairieViewVoltageClampInterface,
)
from surmeier_lab_to_nwb.zhai2025.optogenetics_interfaces import (
    PrairieViewOptogeneticsInterface,
)


def extract_date_from_folder_name(folder_name: str) -> datetime:
    """
    Extract experimental session date from Sr-oEPSC folder naming convention.

    The folder names follow the pattern MMDDYYYY{letter}, where the letter suffix
    distinguishes multiple sessions on the same day (e.g., '07052023a', '08162023b').
    This function parses the date component and returns a timezone-aware datetime
    object set to Central Time Zone (Chicago), which is the timezone of the
    Northwestern University lab where the experiments were conducted.

    Parameters
    ----------
    folder_name : str
        The folder name containing the embedded date in MMDDYYYY format.
        Examples: '07052023a', '08162023b', '07062023c'

    Returns
    -------
    datetime
        A timezone-aware datetime object representing the experimental session date,
        set to midnight Central Time Zone.

    Raises
    ------
    ValueError
        If the folder name does not contain a recognizable 8-digit date pattern.
    """
    # Look for date pattern MMDDYYYY in folder name
    date_match = re.search(r"(\d{8})", folder_name)
    if date_match:
        date_str = date_match.group(1)
        month = int(date_str[:2])
        day = int(date_str[2:4])
        year = int(date_str[4:8])
    else:
        raise ValueError(f"Could not extract date from folder name: {folder_name}")

    # Illinois is in Central Time Zone
    central_tz = ZoneInfo("America/Chicago")

    return datetime(year, month, day, 0, 0, 0, tzinfo=central_tz)


def parse_session_info(session_folder: Path) -> Dict[str, Any]:
    """
    Extract comprehensive session metadata from Sr-oEPSC experimental folder structure.

    This function analyzes the session folder name to extract temporal and organizational
    metadata needed for NWB file creation. The session folder naming convention follows
    the pattern MMDDYYYY{letter}, where the letter suffix allows multiple experimental
    sessions to be conducted on the same calendar date.

    Parameters
    ----------
    session_folder : Path
        Path object pointing to the experimental session folder.
        Expected folder name format: MMDDYYYY{letter} (e.g., '07052023a', '08162023b')

    Returns
    -------
    Dict[str, Any]
        A dictionary containing parsed session metadata with the following keys:

        - 'session_start_time' (datetime): Timezone-aware datetime object
        - 'date_str' (str): Human-readable date string in YYYY-MM-DD format
        - 'session_id' (str): Unique session identifier combining date and suffix
        - 'session_suffix' (str): Letter suffix distinguishing sessions on same date
    """
    folder_name = session_folder.name
    session_start_time = extract_date_from_folder_name(folder_name)

    # Extract session identifier from folder name (e.g., 'a', 'b', 'c')
    session_suffix = folder_name[-1] if folder_name[-1].isalpha() else ""

    return {
        "session_start_time": session_start_time,
        "date_str": f"{session_start_time.year}-{session_start_time.month:02d}-{session_start_time.day:02d}",
        "session_id": f"{session_start_time.year}{session_start_time.month:02d}{session_start_time.day:02d}{session_suffix}",
        "session_suffix": session_suffix,
    }


def parse_sweep_info(sweep_folder_name: str) -> Dict[str, Any]:
    """
    Extract detailed sweep metadata from Sr-oEPSC experimental folder naming.

    The sweep folder naming convention encodes important experimental parameters
    including the recorded cell identifier, LED stimulation parameters, and sweep
    number. This function parses these components to provide structured
    metadata for NWB file organization and identification.

    Parameters
    ----------
    sweep_folder_name : str
        The name of the sweep folder following the pattern:
        'cell{N}_LED{X}-{YYY}' or variations like 'cell{N}_LED{X}.{Z}-{YYY}'

        Examples:
        - 'cell1_LED14-001': Cell 1, LED stimulation 14, sweep 001
        - 'cell2_LED12-006': Cell 2, LED stimulation 12, sweep 006
        - 'cell1_LED9.1-001': Cell 1, LED stimulation 9.1, sweep 001

    Returns
    -------
    Dict[str, Any]
        A dictionary containing parsed sweep components:

        - 'cell_number' (str): Numeric identifier of the recorded cell
        - 'led_number' (str): LED stimulation parameter identifier
        - 'sweep_number' (str): Sequential sweep number within the series
        - 'sweep_name' (str): Complete original folder name for reference

    Notes
    -----
    The LED number may contain decimal points (e.g., '9.1') and represents
    different stimulation protocols or intensities used in the experiment.
    Each sweep contains a single voltage clamp recording with
    LED stimulation pulses delivered every 30 seconds during the sweep.
    """
    # Extract cell number
    cell_match = re.search(r"cell(\d+)", sweep_folder_name)
    cell_number = cell_match.group(1) if cell_match else "unknown"

    # Extract LED stimulation protocol number
    led_match = re.search(r"LED(\d+)", sweep_folder_name)
    led_number = led_match.group(1) if led_match else "unknown"

    # Extract sweep number
    sweep_match = re.search(r"-(\d+)$", sweep_folder_name)
    sweep_number = sweep_match.group(1) if sweep_match else "unknown"

    return {
        "cell_number": cell_number,
        "led_number": led_number,
        "sweep_number": sweep_number,
        "sweep_name": sweep_folder_name,
    }


def convert_data_to_nwb(session_folder_path: Path, condition: str) -> NWBFile:
    """
    Comprehensive conversion of Sr-oEPSC experimental session to NWB format.

    This function orchestrates the complete conversion of a single experimental
    session containing multiple Sr-oEPSC (strontium-substituted optogenetically
    evoked postsynaptic current) sweeps. It processes both voltage clamp
    electrophysiology recordings and associated optogenetic stimulation metadata,
    organizing them into a standardized NWB file structure suitable for
    computational analysis and data sharing.

    Parameters
    ----------
    session_folder_path : Path
        Path to the experimental session folder containing multiple sweep subfolders.
        Each sweep folder should contain Prairie View recording files including:
        - VoltageRecording XML and CSV files (electrophysiology data)
        - VoltageOutput XML files (stimulus protocol definitions)
        - Environment and reference files

    condition : str
        Experimental condition identifier describing the animal model state.
        Examples: "LID on-state", "LID off-state"
        This affects metadata and helps categorize the experimental context.

    Returns
    -------
    NWBFile
        A complete NWB file object containing:

        - Session metadata (timestamps, experimenter info, lab details)
        - Multiple VoltageClampSeries objects (one per sweep)
        - OptogeneticSeries objects with stimulus timing and parameters
        - Device information for recording and stimulation equipment
        - Electrode and stimulus site specifications

    Notes
    -----
    The conversion process includes:

    1. Session-level metadata extraction from folder naming conventions
    2. Sweep-by-sweep processing of electrophysiology recordings
    3. Stimulus parameter extraction from Prairie View XML configurations
    4. Creation of NWB-compliant data structures and relationships
    5. Comprehensive error handling for missing or malformed files

    Each sweep within the experimental session is processed independently, allowing for
    robust handling of incomplete datasets while preserving all available
    experimental data and metadata.
    """
    # Parse session information
    session_info = parse_session_info(session_folder_path)
    print(f"Session date: {session_info['date_str']}, Session: {session_info['session_id']}")

    # Create Subject information
    subject = Subject(
        subject_id=f"mouse_{session_info['session_id']}",
        description=(
            f"Adult Drd1-Tdtomato transgenic mouse used for Sr2+-oEPSC experiments under {condition} condition. "
            f"Mouse received unilateral 6-OHDA lesion (>95% dopamine depletion) in substantia nigra pars compacta "
            f"to model Parkinson's disease. For dyskinesia experiments, mice received dyskinesiogenic levodopa "
            f"treatment (6-12 mg/kg + benserazide) every other day for ≥5 sessions. AAV5-hSyn-hChR2(H134R)-EYFP "
            f"was stereotaxically injected into motor cortex ipsilateral to lesion for optogenetic stimulation "
            f"of corticostriatal terminals."
        ),
        species="Mus musculus",
        strain="Drd1-Tdtomato transgenic",
        sex="U",  # Unknown - not specified in paper
        age="P56D/",  # Adult mice, exact age not specified
        genotype="Drd1-Tdtomato+",
    )

    # Create NWB file with appropriate metadata
    nwbfile = NWBFile(
        session_description=(
            f"Sr2+-oEPSC (strontium-substituted optogenetically evoked postsynaptic current) recordings "
            f"in direct pathway spiny projection neurons (dSPNs) from dyskinetic mice under {condition} condition. "
            f"Voltage clamp recordings at -70 mV holding potential in Ca2+-free ACSF containing 3 mM SrCl2 "
            f"and 10 μM gabazine. ChR2-expressing corticostriatal terminals stimulated with blue LED pulses "
            f"(470 nm, 0.3 ms duration) to evoke asynchronous excitatory postsynaptic currents. "
            f"Recorded from Drd1-Tdtomato mice with unilateral 6-OHDA lesions and established LID."
        ),
        identifier=f"zhai2025_fig2_sr_oepsc_{condition.replace(' ', '_').replace('-', '_')}_{session_info['session_id']}",
        session_start_time=session_info["session_start_time"],
        experimenter=[
            "Zhai, Shenyu",
            "Cui, Qiaoling",
            "Wokosin, David",
            "Sun, Linqing",
            "Tkatch, Tatiana",
            "Crittenden, Jill R.",
            "Graybiel, Ann M.",
            "Surmeier, D. James",
        ],
        lab="Surmeier Lab",
        institution="Northwestern University",
        experiment_description=(
            f"Investigation of corticostriatal synaptic strength in levodopa-induced dyskinesia (LID) "
            f"using Sr2+-oEPSC measurements under {condition} condition. AAV5-hSyn-hChR2(H134R)-EYFP "
            f"was injected into motor cortex ipsilateral to 6-OHDA lesion to express channelrhodopsin-2 "
            f"in corticostriatal terminals. Mice received dyskinesiogenic levodopa treatment (6-12 mg/kg "
            f"+ 12 mg/kg benserazide) and were sacrificed either 30 min post-dose (on-state) or 24-48h "
            f"post-dose (off-state). Whole-cell voltage clamp recordings from dSPNs in Ca2+-free ACSF "
            f"with Sr2+ substitution to measure asynchronous release from optogenetically activated "
            f"corticostriatal synapses. Part of Figure 2 from Zhai et al. 2025 studying state-dependent "
            f"modulation of synaptic connectivity in Parkinson's disease model."
        ),
        session_id=session_info["session_id"],
        subject=subject,
        notes=(
            "EXPERIMENTAL PROTOCOLS: "
            "1. Animal Model: Drd1-Tdtomato mice with unilateral 6-OHDA lesions (>95% DA depletion). "
            "2. Dyskinesia Induction: Dyskinesiogenic levodopa doses (6 mg/kg first 2 sessions, "
            "12 mg/kg later sessions + 12 mg/kg benserazide) every other day for ≥5 sessions. "
            "3. State Definition: On-state = 30 min post-levodopa; Off-state = 24-48h post-levodopa. "
            "4. ChR2 Expression: AAV5-hSyn-hChR2(H134R)-EYFP injected into motor cortex. "
            "5. Recording Solution: Ca2+-free ACSF + 3 mM SrCl2 + 10 μM gabazine, 25 min incubation. "
            "6. Pipette Solution: 120 mM CsMeSO3, 5 mM NaCl, 0.25 mM EGTA, 10 mM HEPES, "
            "4 mM Mg-ATP, 0.3 mM Na-GTP, 10 mM TEA, 5 mM QX-314 (pH 7.25, 280-290 mOsm). "
            "7. Recording Parameters: Voltage clamp at -70 mV, pipette resistance 3-4 MΩ, "
            "stimulation every 30s, analysis window 40-400ms post-stimulus for asynchronous events. "
            "Each sweep represents one voltage clamp recording session."
        ),
    )

    # Get all sweep folders
    sweep_folders = [f for f in session_folder_path.iterdir() if f.is_dir() and "cell" in f.name.lower()]
    sweep_folders.sort()

    print(f"Found {len(sweep_folders)} sweeps in experimental session")

    # Process each sweep
    for sweep_folder in sweep_folders:
        print(f"  Processing: {sweep_folder.name}")

        # Parse sweep information
        sweep_info = parse_sweep_info(sweep_folder.name)

        # Find XML recording file
        xml_files = list(sweep_folder.glob("*_VoltageRecording_*.xml"))
        if not xml_files:
            raise FileNotFoundError(
                f"Required VoltageRecording XML file not found in sweep folder: {sweep_folder.name}"
            )

        xml_recording_file = xml_files[0]

        # Create voltage clamp interface
        voltage_clamp_interface = PrairieViewVoltageClampInterface(
            file_path=xml_recording_file,
            icephys_metadata_key=f"VoltageClamp_Cell{sweep_info['cell_number']}_LED{sweep_info['led_number']}_Sweep{sweep_info['sweep_number']}",
        )

        # Add voltage clamp data to NWB file
        voltage_clamp_interface.add_to_nwbfile(nwbfile=nwbfile)

        # Find VoltageOutput XML file and create optogenetics interface
        voltage_output_files = list(sweep_folder.glob("*_VoltageOutput_*.xml"))
        if voltage_output_files:
            optogenetics_interface = PrairieViewOptogeneticsInterface(
                voltage_output_xml_path=voltage_output_files[0],
                sweep_info=sweep_info,
                optogenetics_metadata_key=f"Optogenetics_Cell{sweep_info['cell_number']}_LED{sweep_info['led_number']}_Sweep{sweep_info['sweep_number']}",
            )

            # Add optogenetic stimulus to NWB file
            optogenetics_interface.add_to_nwbfile(nwbfile=nwbfile)
            print(
                f"    Found LED stimulus: {optogenetics_interface.stimulus_params['pulse_width_ms']}ms pulse at {optogenetics_interface.stimulus_params['pulse_amplitude_v']}V"
            )
        else:
            print(f"    Warning: No VoltageOutput XML file found in {sweep_folder.name}")

        print(
            f"    Successfully added sweep: Cell {sweep_info['cell_number']}, LED {sweep_info['led_number']}, Sweep {sweep_info['sweep_number']}"
        )

    return nwbfile


if __name__ == "__main__":
    import logging

    logging.getLogger("tifffile").setLevel(logging.ERROR)

    # Define the base path to the data
    base_path = Path(
        "/media/heberto/One Touch/Surmeier-CN-data-share/consolidated_data/LID_paper_Zhai_2025/Raw data for Figs/Figure 2_SF1A/Sr-oEPSC/"
    )

    conditions = ["LID on-state", "LID off-state"]

    for condition in conditions:
        condition_path = base_path / condition
        if not condition_path.exists():
            print(f"Warning: Condition path does not exist: {condition_path}")
            continue

        print(f"Processing Sr-oEPSC data for: {condition}")

        # Get all session folders
        session_folders = [f for f in condition_path.iterdir() if f.is_dir()]
        session_folders.sort()

        print(f"Found {len(session_folders)} session folders")

        for session_folder in session_folders:
            print(f"\nProcessing session: {session_folder.name}")

            # Convert data to NWB format
            nwbfile = convert_data_to_nwb(
                session_folder_path=session_folder,
                condition=condition,
            )

            # Generate output filename
            safe_condition = condition.replace(" ", "_").replace("-", "_")
            nwbfile_path = Path(f"sr_oepsc_{safe_condition}_{session_folder.name}.nwb")

            configure_and_write_nwbfile(nwbfile, nwbfile_path=nwbfile_path)
            print(f"Successfully saved: {nwbfile_path}")
