# -*- coding: utf-8 -*-
"""
Figure 2 Optical Stimuli Conversion Script - Zhai et al. 2025
============================================================

This script converts Sr²⁺-oEPSC (strontium-substituted optogenetically evoked
postsynaptic current) experimental data from Figure 2 of Zhai et al. 2025 into
NWB (Neurodata Without Borders) format.

EXPERIMENTAL CONTEXT:
===================
The data comes from experiments investigating corticostriatal synaptic strength
in a Parkinson's disease model with levodopa-induced dyskinesia (LID). Mice
received:
- Unilateral 6-OHDA lesions (>95% dopamine depletion)
- Dyskinesiogenic levodopa treatment
- AAV5-hSyn-hChR2(H134R)-EYFP injection into motor cortex

Two experimental conditions were compared:
- LID on-state: Sacrificed 30 min post-levodopa dose
- LID off-state: Sacrificed 24-48h post-levodopa dose

DATA STRUCTURE:
==============
The raw data follows a hierarchical organization:

Base Path/
├── LID on-state/
│   ├── 07052023a/          # Session A on July 5, 2023
│   │   ├── cell1_LED14-001/    # Cell 1, LED intensity 14, sweep 001
│   │   │   ├── *_VoltageRecording_*.xml/csv  # Electrophysiology data
│   │   │   ├── *_VoltageOutput_*.xml         # Optogenetic parameters
│   │   │   └── *_Environment.xml             # Recording configuration
│   │   ├── cell1_LED14-002/    # Next sweep (30 seconds later)
│   │   └── ...
│   ├── 07062023a/          # Session A on July 6, 2023
│   └── ...
└── LID off-state/
    ├── 07132023a/          # Different session dates
    └── ...

TIMING VERIFICATION:
Actual timestamps from XML files confirm 30-second intervals:
2023-07-05T15:19:12 → cell1_LED14-001
2023-07-05T15:19:42 → cell1_LED14-002 (+30 seconds)
2023-07-05T15:20:12 → cell1_LED14-003 (+30 seconds)
...

OPTOGENETIC PARAMETERS:
=====================
From VoltageOutput XML metadata:
- Pulse duration: 0.3 milliseconds (matches paper specification)
- LED voltage: 5V (corresponds to different intensity settings)
- Pulse timing: 20ms delay + 0.3ms pulse + recording window
- Repetitions: Single pulse per sweep
- Inter-sweep interval: 30 seconds

CONVERSION OUTPUT:
================
NWB files are saved to: nwb_files/figure_2_optical_stimuli/
Filename format: sr_oepsc_{condition}_{session_folder}.nwb

Each NWB file contains:
- Session metadata with precise timing
- Subject information with experimental condition
- Multiple VoltageClampSeries (one per sweep)
- OptogeneticSeries with stimulus parameters
- Device specifications for recording and stimulation equipment

This conversion enables standardized analysis of the hierarchical data structure
(Session → Cell → Sweep → Individual EPSCs) required for Figure 2's statistical
comparisons between LID on-state and off-state conditions.
"""
import re
from datetime import datetime
from pathlib import Path
from typing import Any, Dict
from zoneinfo import ZoneInfo

from neuroconv.tools import configure_and_write_nwbfile
from neuroconv.utils import dict_deep_update, load_dict_from_file
from pynwb import NWBFile
from pynwb.file import Subject

from surmeier_lab_to_nwb.zhai2025.intracellular_interfaces import (
    PrairieViewVoltageClampInterface,
)
from surmeier_lab_to_nwb.zhai2025.optogenetics_interfaces import (
    PrairieViewOptogeneticsInterface,
)


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
    Convert Sr²⁺-oEPSC experimental session to NWB format for Figure 2 analysis.

    This function converts a complete experimental session from the Zhai et al. 2025
    Figure 2 optogenetic experiments. Each session represents a single day of experiments
    with consistent conditions, containing multiple cells and sweeps for statistical analysis.

    DATA ORGANIZATION:
    ==================

    Sessions (Date-based folders):
    - Format: MMDDYYYY{letter} (e.g., 07052023a, 07062023b)
    - Represents: Single experimental day with consistent conditions
    - Contains: Multiple cells recorded under same protocol
    - Conditions: Either "LID on-state" or "LID off-state"

    Sweeps (Individual recordings):
    - Format: cell{N}_LED{X}-{YYY} (e.g., cell1_LED14-001)
    - Represents: One complete optogenetic stimulation + recording episode
    - Timing: Exact 30-second intervals between sweeps (verified from timestamps)
    - Contains: Single 0.3ms blue LED pulse + Sr²⁺-oEPSC recording

    OPTOGENETIC STIMULUS PARAMETERS (from XML metadata):
    ====================================================
    - Pulse duration: 0.3 milliseconds (matches paper exactly)
    - LED voltage: 5V drive (corresponds to different LED intensities)
    - Timing: 20ms baseline + 0.3ms pulse + ~380ms recording window
    - Frequency: Every 30 seconds between sweeps
    - Target: ChR2-expressing corticostriatal terminals

    EXPERIMENTAL TIMELINE PER SWEEP:
    ================================
    0-20 ms:     Baseline recording (LED off, 0V)
    20-20.3 ms:  Light pulse (LED on, 5V) → ChR2 activation
    20.3-400 ms: Recording window for asynchronous EPSCs
    400+ ms:     End sweep, wait 30 seconds for next sweep

    Parameters
    ----------
    session_folder_path : Path
        Path to experimental session folder (date-based naming).
        Contains multiple sweep subfolders with Prairie View files:
        - VoltageRecording XML/CSV: Electrophysiology data
        - VoltageOutput XML: Optogenetic stimulus parameters
        - Configuration files: Recording setup metadata

    condition : str
        Experimental condition: "LID on-state" or "LID off-state"
        - On-state: 30 min post-levodopa injection
        - Off-state: 24-48h post-levodopa injection

    Returns
    -------
    NWBFile
        Complete NWB file containing:
        - Session metadata with precise timing from folder names
        - VoltageClampSeries for each sweep (Sr²⁺-oEPSC recordings)
        - OptogeneticSeries with 0.3ms pulse parameters
        - Device specifications for recording and LED stimulation
        - Subject information with experimental condition

    Notes
    -----
    CONVERSION METHODOLOGY:
    1. Parse session start time from folder name (MMDDYYYY format)
    2. Load general metadata from YAML (lab, institution, experimenters)
    3. Create NWBFile with session-specific metadata outside processing loop
    4. Process each sweep independently:
       - Create VoltageClampInterface for electrophysiology data
       - Create OptogeneticsInterface for LED stimulus parameters
       - Add both to existing NWBFile
    5. Maintain hierarchical organization: Session → Cell → Sweep → EPSCs

    STATISTICAL ANALYSIS STRUCTURE:
    - Individual EPSC events detected within each sweep
    - Multiple sweeps per cell provide repeated measurements
    - Multiple cells per session provide biological replicates
    - Multiple sessions per condition provide experimental replicates

    This hierarchical organization enables proper statistical analysis while
    maintaining experimental traceability for Figure 2's comparative analysis.
    """
    # Parse session start time from folder name (MMDDYYYY pattern)
    folder_name = session_folder_path.name
    # Extract first 8 digits from folder name (MMDDYYYY format)
    date_str = folder_name[:8]
    # Parse MMDDYYYY format
    session_start_time = datetime.strptime(date_str, "%m%d%Y")
    # Add Central Time Zone (Illinois)
    central_tz = ZoneInfo("America/Chicago")
    session_start_time_with_tz = session_start_time.replace(tzinfo=central_tz)

    session_id = session_folder_path.name
    print(f"Processing session folder: {session_folder_path.name}")

    # Load metadata from YAML file
    metadata_file_path = Path(__file__).parent / "metadata.yaml"
    paper_metadata = load_dict_from_file(metadata_file_path)

    # Create Subject information
    subject = Subject(
        subject_id=f"mouse_{session_id}",
        description=(
            f"Adult Drd1-Tdtomato transgenic mouse with unilateral 6-OHDA lesion (>95% dopamine depletion) "
            f"modeling Parkinson's disease. Received dyskinesiogenic levodopa treatment and AAV5-hSyn-hChR2(H134R)-EYFP "
            f"injection into motor cortex for optogenetic experiments. Experimental condition: {condition}."
        ),
        species="Mus musculus",
        strain="Drd1-Tdtomato transgenic",
        sex="U",  # Unknown - not specified in paper
        age="P56D/",  # Adult mice, exact age not specified
        genotype="Drd1-Tdtomato+",
    )

    # Session description for this specific recording
    session_description = (
        f"Sr²⁺-oEPSC recordings from dorsal striatal projection neurons (dSPNs) under {condition} condition. "
        f"Voltage clamp at -70 mV in Ca²⁺-free ACSF with 3 mM SrCl₂ and 10 μM gabazine. "
        f"Optogenetic stimulation: 0.3ms blue LED pulses (5V drive) activating ChR2-expressing "
        f"corticostriatal terminals every 30 seconds. Recording window: 20ms baseline + 0.3ms pulse + "
        f"~380ms for asynchronous EPSC detection. Multiple sweeps per cell for statistical analysis."
    )

    # Update with session-specific metadata
    session_specific_metadata = {
        "NWBFile": {
            "session_description": session_description,
            "identifier": f"zhai2025_fig2_sr_oepsc_{condition.replace(' ', '_').replace('-', '_')}_{session_id}",
            "session_id": session_id,
            "session_start_time": session_start_time_with_tz,
            "keywords": ["Sr2+-oEPSC", "corticostriatal terminals", "dSPNs", "ChR2", "blue LED stimulation"],
            "notes": (
                "DETAILED EXPERIMENTAL PROTOCOLS: "
                "Animal preparation: Dyskinesiogenic levodopa (6-12 mg/kg + 12 mg/kg benserazide) every other day for ≥5 sessions. "
                "Recording solution: Ca2+-free ACSF + 3 mM SrCl2 + 10 μM gabazine, 25 min incubation. "
                "Pipette solution: 120 mM CsMeSO3, 5 mM NaCl, 0.25 mM EGTA, 10 mM HEPES, 4 mM Mg-ATP, "
                "0.3 mM Na-GTP, 10 mM TEA, 5 mM QX-314 (pH 7.25, 280-290 mOsm). "
                "Recording parameters: -70 mV holding potential, 3-4 MΩ pipette resistance, 30s stimulation interval, "
                "40-400ms analysis window for asynchronous events."
            ),
        }
    }

    # Deep update the metadata with session-specific information
    metadata = dict_deep_update(paper_metadata, session_specific_metadata)

    # Create NWBFile directly
    nwbfile = NWBFile(
        session_description=metadata["NWBFile"]["session_description"],
        identifier=metadata["NWBFile"]["identifier"],
        session_start_time=metadata["NWBFile"]["session_start_time"],
        experimenter=metadata["NWBFile"]["experimenter"],
        lab=metadata["NWBFile"]["lab"],
        institution=metadata["NWBFile"]["institution"],
        experiment_description=metadata["NWBFile"]["experiment_description"],
        session_id=metadata["NWBFile"]["session_id"],
        subject=subject,
        keywords=metadata["NWBFile"]["keywords"],
        notes=metadata["NWBFile"]["notes"],
    )

    # Get all sweep folders
    sweep_folders = [f for f in session_folder_path.iterdir() if f.is_dir()]

    print(f"Found {len(sweep_folders)} sweeps in experimental session")

    if not sweep_folders:
        raise ValueError(f"No sweep folders found in session: {session_folder_path} for condition: {condition}")

    # Process each sweep
    for sweep_index, sweep_folder in enumerate(sweep_folders):
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

        # Create voltage clamp interface with consistent sweep numbering
        sweep_number_padded = f"{sweep_index + 1:03d}"  # 1-indexed, zero-padded to 3 digits
        voltage_clamp_interface = PrairieViewVoltageClampInterface(
            file_path=xml_recording_file,
            icephys_metadata_key=f"VoltageClampSweep{sweep_number_padded}",
        )

        # Add voltage clamp data to NWB file
        voltage_clamp_interface.add_to_nwbfile(nwbfile=nwbfile)

        # Find VoltageOutput XML file and create optogenetics interface
        voltage_output_files = list(sweep_folder.glob("*_VoltageOutput_*.xml"))
        if not voltage_output_files:
            raise FileNotFoundError(f"Required VoltageOutput XML file not found in sweep folder: {sweep_folder.name}")

        # Create optogenetics interface with consistent sweep numbering
        optogenetics_interface = PrairieViewOptogeneticsInterface(
            voltage_output_xml_path=voltage_output_files[0],
            sweep_info=sweep_info,
            optogenetics_metadata_key=f"OptogeneticSweep{sweep_number_padded}",
        )

        # Add optogenetic stimulus to NWB file
        optogenetics_interface.add_to_nwbfile(nwbfile=nwbfile)

        print(
            f"    Found LED stimulus: {optogenetics_interface.stimulus_params['pulse_width_ms']}ms pulse at {optogenetics_interface.stimulus_params['pulse_amplitude_v']}V control voltage (optical power unknown)"
        )

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
            raise FileNotFoundError(f"Base path does not exist: {condition_path}")

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

            # Generate output filename and create output directory structure
            safe_condition = condition.replace(" ", "_").replace("-", "_")

            # Create nwb_files directory at root level
            root_dir = Path(__file__).parent.parent.parent.parent  # Go up to repo root
            nwb_files_dir = root_dir / "nwb_files" / "figure_2_optical_stimuli"
            nwb_files_dir.mkdir(parents=True, exist_ok=True)

            nwbfile_path = nwb_files_dir / f"sr_oepsc_{safe_condition}_{session_folder.name}.nwb"

            configure_and_write_nwbfile(nwbfile, nwbfile_path=nwbfile_path)
            print(f"Successfully saved: {nwbfile_path}")
