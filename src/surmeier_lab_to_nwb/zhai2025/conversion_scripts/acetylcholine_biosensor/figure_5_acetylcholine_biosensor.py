# -*- coding: utf-8 -*-
"""
Figure 5 Acetylcholine Biosensor Conversion Script - Zhai et al. 2025
====================================================================

This script converts GRABACh3.0 acetylcholine biosensor imaging data from Figure 5
of Zhai et al. 2025 into NWB (Neurodata Without Borders) format.

The data examines acetylcholine release dynamics using genetically encoded biosensors
in a mouse model of Parkinson's disease and levodopa-induced dyskinesia, revealing
state-dependent modulation of cholinergic signaling.

For detailed experimental context, protocols, data structure, and analysis methods,
see: /src/surmeier_lab_to_nwb/zhai2025/conversion_notes_folder/figure_5_conversion_notes.md
"""

import re
import uuid
from datetime import datetime
from pathlib import Path
from typing import Any
from zoneinfo import ZoneInfo

from neuroconv.converters import BrukerTiffSinglePlaneConverter
from neuroconv.tools import configure_and_write_nwbfile
from neuroconv.utils import dict_deep_update, load_dict_from_file
from pynwb import NWBFile
from pynwb.device import Device
from pynwb.epoch import TimeIntervals
from pynwb.file import Subject

from surmeier_lab_to_nwb.zhai2025.interfaces import (
    PrairieViewFluorescenceInterface,
)


def parse_session_info_from_folder_name(session_folder: Path) -> dict[str, Any]:
    """
    Parse session information from BOT session folder names.

    Expected format: BOT_[date]_[slice_info]_[treatment]_[stimulation]-[session_num]
    Examples:
    - BOT_04162024_slice2ROI1_50nMDA_burst-001
    - BOT_05242024_slice1A_ctr_single-001
    - BOT_04052024_slice2_ACh-001 (calibration)

    Parameters
    ----------
    session_folder : Path
        Path to the session folder

    Returns
    -------
    Dict[str, Any]
        Dictionary containing session information
    """
    session_name = session_folder.name

    # Parse session name: BOT_[date]_[slice_info]_[treatment]_[stimulation]-[session_num]
    pattern = r"BOT_(\d{8})_(.+?)_([^_]+)_([^-]+)-(\d+)"
    match = re.match(pattern, session_name)

    if not match:
        # Handle calibration sessions (ACh, TTX) - no stimulation protocol
        calibration_pattern = r"BOT_(\d{8})_(.+?)_([^-]+)-(\d+)"
        calibration_match = re.match(calibration_pattern, session_name)

        if calibration_match:
            date_str, slice_info, treatment, session_num = calibration_match.groups()

            # Check if this is a real calibration session (ACh, TTX) or missing stimulation type
            if treatment in ["ACh", "TTX"] or "ACh" in treatment or "TTX" in treatment:
                stimulation = "calibration"
            else:
                # Experimental session missing stimulation type - assume single pulse
                stimulation = "single"
                print(f"Warning: Session {session_name} missing stimulation type, assuming single pulse")
        else:
            raise ValueError(f"Could not parse session folder name: {session_name}")
    else:
        date_str, slice_info, treatment, stimulation, session_num = match.groups()

        # Check if stimulation field contains calibration markers
        if stimulation in ["ACh", "TTX"] or "ACh" in stimulation or "TTX" in stimulation:
            # This is actually a calibration session, adjust treatment and stimulation
            treatment = f"{treatment}_{stimulation}_calibration"
            stimulation = "calibration"

    # Parse date (MMDDYYYY format)
    month = int(date_str[0:2])
    day = int(date_str[2:4])
    year = int(date_str[4:8])
    session_date = datetime(year, month, day)

    # Map treatment abbreviations to full names
    treatment_mapping = {
        "ctr": "control",
        "50nMDA": "50nM_dopamine",
        "quin": "quinpirole",
        "sul": "sulpiride",
        "ACh": "acetylcholine_calibration",
        "TTX": "TTX_calibration",
    }

    treatment_full = treatment_mapping.get(treatment, treatment)

    # Map stimulation types (including typo correction)
    stimulation_mapping = {
        "single": "single_pulse",
        "sinlge": "single_pulse",  # Fix typo in raw data folders
        "singl": "single_pulse",  # Fix truncated version in raw data folders
        "burst": "burst_stimulation",
        "20Hz": "burst_stimulation",  # Frequency-based naming for burst stimulation
        "calibration": "calibration",
    }

    stimulation_full = stimulation_mapping.get(stimulation, stimulation)

    return {
        "session_name": session_name,
        "session_date": session_date,
        "date_str": session_date.strftime("%Y-%m-%d"),
        "slice_info": slice_info,
        "treatment": treatment_full,
        "stimulation": stimulation_full,
        "session_number": session_num,
        "is_calibration": treatment in ["ACh", "TTX"] or "ACh" in treatment or "TTX" in treatment,
    }


def convert_session_to_nwbfile(session_folder: Path, condition: str, verbose: bool = False) -> NWBFile:
    """
    Convert a single session of Figure 5 GRABACh3.0 data to NWB format.

    Parameters
    ----------
    session_folder : Path
        Path to the session folder
    condition : str
        Experimental condition ("UL control", "PD", "LID off")
    verbose : bool, default=False
        Enable verbose output

    Returns
    -------
    NWBFile
        NWB file with the converted data
    """
    # Define condition descriptions for metadata
    condition_descriptions = {
        "UL control": "Unlesioned control mouse with healthy striatum and intact dopaminergic system",
        "PD": "6-OHDA lesioned mouse (>95 per cent dopamine depletion) modeling Parkinson's disease",
        "LID off": "Dyskinetic mouse in off-state (24-48h post-levodopa) with established levodopa-induced dyskinesia",
    }

    # Parse session information
    session_info = parse_session_info_from_folder_name(session_folder)

    if verbose:
        print(f"Processing session: {session_info['session_name']}")
        print(f"  Treatment: {session_info['treatment']}")
        print(f"  Stimulation: {session_info['stimulation']}")

    # Find required files for BOT interface
    bot_csv_file = None
    xml_metadata_file = None

    for file in session_folder.iterdir():
        if file.name.endswith("-botData.csv"):
            bot_csv_file = file
        elif file.name.endswith(".xml") and "VoltageRecording" not in file.name and "VoltageOutput" not in file.name:
            xml_metadata_file = file

    if not bot_csv_file:
        raise FileNotFoundError(f"No botData.csv file found in {session_folder}")

    if not xml_metadata_file:
        raise FileNotFoundError(f"No XML metadata file found in {session_folder}")

    # Load metadata from YAML file
    metadata_file_path = Path(__file__).parent.parent.parent / "metadata.yaml"
    paper_metadata = load_dict_from_file(metadata_file_path)

    # Create session-specific metadata with Chicago timezone
    central_tz = ZoneInfo("America/Chicago")
    session_start_time = datetime.combine(session_info["session_date"], datetime.min.time()).replace(tzinfo=central_tz)

    # Create BIDS-style base session ID with detailed timestamp when available
    if hasattr(session_start_time, "hour"):
        timestamp = session_start_time.strftime("%Y%m%d_%H%M%S")
    else:
        timestamp = session_start_time.strftime("%Y%m%d")

    base_session_id = f"figure5_AcetylcholineGRAB_{condition.replace(' ', '_').replace('-', '_')}_{timestamp}_Sub{session_info['slice_info']}"
    session_id = f"{base_session_id}_Session{session_info['session_number']}"

    session_specific_metadata = {
        "NWBFile": {
            "session_description": (
                f"GRABACh3.0 acetylcholine sensor imaging session from striatal slice - {condition} condition. "
                f"Treatment: {session_info['treatment']}, Stimulation: {session_info['stimulation']}, "
                f"Session {session_info['session_number']} from {session_info['slice_info']}. Two-photon microscopy "
                f"at 920nm excitation measuring ACh release dynamics with electrical stimulation."
            ),
            "identifier": str(uuid.uuid4()),
            "session_start_time": session_start_time,
            "experiment_description": (
                f"Figure 5 GRABACh3.0 experiment from Zhai et al. 2025 investigating acetylcholine release "
                f"dynamics in Parkinson's disease and levodopa-induced dyskinesia. Single session recording "
                f"with {session_info['treatment']} treatment and {session_info['stimulation']} protocol."
            ),
            "session_id": session_id,
            "keywords": [
                "GRABACh3.0",
                "acetylcholine",
                "two-photon microscopy",
                "cholinergic interneurons",
                "levodopa-induced dyskinesia",
                "pharmacology",
                "electrical stimulation",
            ],
        },
        "Subject": {
            "subject_id": f"grabatch_mouse_{session_info['slice_info']}",
            "description": (
                f"{condition_descriptions[condition]}. Striatal injection of AAV-GRABACh3.0 for "
                f"acetylcholine sensor expression in cholinergic interneurons. Session recorded on "
                f"{session_info['date_str']} from {session_info['slice_info']}."
            ),
            "genotype": "Wild-type with AAV-GRABACh3.0",
        },
    }

    # Deep merge with paper metadata
    metadata = dict_deep_update(paper_metadata, session_specific_metadata)

    # Create NWB file
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

    # Add stimulation device and electrode information
    stimulation_device = Device(
        name="DeviceConcentricBipolarElectrode",
        description=(
            "Concentric bipolar stimulating electrode (CBAPD75, FHC) for electrical "
            "stimulation of striatal tissue. Placed 200 μm ventral to GRABACh3.0 "
            "imaging region for acetylcholine release experiments."
        ),
    )
    nwbfile.add_device(stimulation_device)

    # Add imaging device for BOT interface
    imaging_device = Device(
        name="default",  # Name expected by BrukerTiffSinglePlaneConverter
        description=(
            "Bruker two-photon microscope for acetylcholine GRAB biosensor imaging. "
            "Used for brightness over time (BOT) measurements of acetylcholine release dynamics."
        ),
    )
    nwbfile.add_device(imaging_device)

    # Add stimulation protocol information
    general_stimulation_description = (
        f"Electrical stimulation using concentric bipolar electrode (CBAPD75, FHC) "
        f"placed 200 μm ventral to imaging region. {session_info['stimulation']} protocol "
        f"as described in methods. Stimulation delivered at t=3s after baseline start. "
        f"XML metadata shows LED stimulator with different parameters but paper "
        f"specifications take precedence."
    )

    if session_info["stimulation"] == "single_pulse":
        session_stimulation_description = "Single electrical pulse: 1 ms duration, 0.3 mA amplitude"
    elif session_info["stimulation"] == "burst_stimulation":
        session_stimulation_description = "Burst stimulation: 20 pulses at 20 Hz, 1 ms duration, 0.3 mA amplitude each"
    else:
        session_stimulation_description = "Calibration protocol - no electrical stimulation"

    # Add stimulus information as session notes
    description = f"{general_stimulation_description} {session_stimulation_description}"

    # Add stimulation epochs if not calibration session
    if not session_info["is_calibration"]:
        # Create stimulus table with custom columns
        stimulus_table = TimeIntervals(
            name="stimulus_table",
            description=description,
        )

        # Add custom columns for stimulation parameters
        stimulus_table.add_column(name="stimulus_type", description="Type of stimulation protocol")
        stimulus_table.add_column(name="amplitude", description="Stimulation amplitude")
        stimulus_table.add_column(name="pulse_width", description="Pulse width duration")
        stimulus_table.add_column(name="frequency", description="Stimulation frequency")
        stimulus_table.add_column(name="electrode", description="Electrode model and type")
        stimulus_table.add_column(name="notes", description="Additional stimulation details")

        # Add stimulus timing (typically at 3s after baseline, lasting brief duration)
        stimulation_start_time = 3.0  # seconds after session start

        if session_info["stimulation"] == "single_pulse":
            stimulation_duration = 0.001  # 1 ms
        elif session_info["stimulation"] == "burst_stimulation":
            stimulation_duration = 1.0  # 20 pulses at 20 Hz = ~1 second
        else:
            raise ValueError(f"Unknown stimulation type for experimental session: {session_info['stimulation']}")

        stimulus_table.add_interval(
            start_time=stimulation_start_time,
            stop_time=stimulation_start_time + stimulation_duration,
            stimulus_type=session_info["stimulation"],
            amplitude="0.3 mA",
            pulse_width="1 ms" if "pulse" in session_info["stimulation"] else "1 ms per pulse",
            frequency="single" if session_info["stimulation"] == "single_pulse" else "20 Hz",
            electrode="CBAPD75 concentric bipolar",
            notes=session_stimulation_description,
        )

        # Add to stimulus namespace
        nwbfile.add_stimulus(stimulus_table)

    # Create Bruker TIFF converter for raw imaging data
    bruker_converter = BrukerTiffSinglePlaneConverter(folder_path=session_folder)

    # Get metadata for Bruker converter and add raw imaging data to NWB file
    bruker_metadata = bruker_converter.get_metadata()
    bruker_converter.add_to_nwbfile(nwbfile=nwbfile, metadata=bruker_metadata)

    # Create acetylcholine fluorescence interface for BOT data
    acetylcholine_fluorescence_interface = PrairieViewFluorescenceInterface(
        bot_csv_data_file_path=bot_csv_file, xml_metadata_file_path=xml_metadata_file
    )

    # Add fluorescence data to NWB file
    acetylcholine_fluorescence_interface.add_to_nwbfile(nwbfile=nwbfile)

    if verbose:
        print(f"  Successfully processed session")

    return nwbfile


if __name__ == "__main__":
    import logging
    import warnings

    from tqdm import tqdm

    # Control verbose output
    verbose = False  # Set to True for detailed output
    stub_test = True  # Set to True to process only first 2 files per condition for testing

    # Suppress warnings
    logging.getLogger("tifffile").setLevel(logging.ERROR)
    warnings.filterwarnings("ignore", message="invalid value encountered in divide")
    warnings.filterwarnings("ignore", message=".*no datetime before year 1.*")

    # Define the base path to the data
    base_path = Path(
        "/media/heberto/One Touch/Surmeier-CN-data-share/consolidated_data/LID_paper_Zhai_2025/Raw data for Figs/Figure 5_SF2"
    )
    if not base_path.exists():
        raise FileNotFoundError(f"Base path does not exist: {base_path}")

    # Create nwb_files directory at root level
    root_dir = Path(__file__).parent.parent.parent.parent.parent.parent  # Go up to repo root
    nwb_files_dir = root_dir / "nwb_files" / "acetylcholine_biosensor" / "figure_5"
    nwb_files_dir.mkdir(parents=True, exist_ok=True)

    # Figure 5 acetylcholine conditions
    conditions = ["UL control", "PD", "LID off"]

    for condition in conditions:
        condition_path = base_path / condition

        if not condition_path.exists():
            raise FileNotFoundError(f"Condition path does not exist: {condition_path}")

        if verbose:
            print(f"Processing acetylcholine GRABACh data for: {condition}")

        # Get all session folders for this condition
        # Structure: condition/parent_folder/BOT_session_folders
        all_sessions = []
        for parent_folder in condition_path.iterdir():
            if parent_folder.is_dir():
                # Get all BOT session folders in this parent folder
                bot_session_folders = [f for f in parent_folder.iterdir() if f.is_dir() and f.name.startswith("BOT_")]
                bot_session_folders.sort()

                for session_folder_path in bot_session_folders:
                    all_sessions.append(
                        {
                            "session_folder_path": session_folder_path,
                            "parent_folder": parent_folder,
                        }
                    )

        if verbose:
            print(f"Found {len(all_sessions)} session folders")

        # Apply stub_test filtering if enabled
        if stub_test:
            all_sessions = all_sessions[:2]
            if verbose:
                print(f"stub_test enabled: processing only first {len(all_sessions)} session folders")

        # Use tqdm for progress bar when verbose is disabled
        session_iterator = tqdm(
            all_sessions, desc=f"Converting Figure5 AcetylcholineBiosensor {condition}", unit=" session"
        )

        for session_info in session_iterator:
            session_folder_path = session_info["session_folder_path"]
            parent_folder = session_info["parent_folder"]

            if verbose:
                print(f"\nProcessing session: {session_folder_path.name}")

            # Convert session to NWB format
            nwbfile = convert_session_to_nwbfile(
                session_folder=session_folder_path,
                condition=condition,
                verbose=verbose,
            )

            # Create output filename
            condition_safe = condition.replace(" ", "_").replace("-", "_")
            nwbfile_path = (
                nwb_files_dir
                / f"figure5_biosensor_{condition_safe}_{parent_folder.name}_{session_folder_path.name}.nwb"
            )

            configure_and_write_nwbfile(nwbfile, nwbfile_path=nwbfile_path)

            if verbose:
                print(f"Successfully saved: {nwbfile_path.name}")
