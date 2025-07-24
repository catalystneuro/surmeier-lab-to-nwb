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
from neuroconv.tools.nwb_helpers import make_nwbfile_from_metadata
from neuroconv.utils import dict_deep_update, load_dict_from_file
from pynwb.device import Device
from pynwb.epoch import TimeIntervals

from surmeier_lab_to_nwb.zhai2025.interfaces import (
    PrairieViewFluorescenceInterface,
)


def parse_slice_session_info(slice_folder: Path) -> dict[str, Any]:
    """
    Parse session information from slice folder names (e.g., 04022024slice1ROI1).

    Parameters
    ----------
    slice_folder : Path
        Path to the slice folder containing BOT recordings

    Returns
    -------
    Dict[str, Any]
        Dictionary containing slice session information
    """
    folder_name = slice_folder.name

    # Parse slice folder name: [date][slice][ROI] (e.g., 04022024slice1ROI1)
    pattern = r"(\d{8})(slice\d+)(ROI\d+)?"
    match = re.match(pattern, folder_name)

    if not match:
        # Try alternative pattern without ROI
        pattern = r"(\d{8})(slice\d+[A-Z]?)"
        match = re.match(pattern, folder_name)

    if not match:
        raise ValueError(f"Could not parse slice folder name: {folder_name}")

    date_str = match.group(1)
    slice_info = match.group(2)
    roi_info = match.group(3) if match.lastindex >= 3 else None

    # Parse date (MMDDYYYY format)
    month = int(date_str[0:2])
    day = int(date_str[2:4])
    year = int(date_str[4:8])
    session_date = datetime(year, month, day)

    full_slice_info = slice_info + (roi_info or "")

    return {
        "slice_folder_name": folder_name,
        "session_date": session_date,
        "date_str": session_date.strftime("%Y-%m-%d"),
        "slice_info": full_slice_info,
        "subject_id": f"mouse_{date_str}",  # Use date as animal proxy
    }


def parse_trial_info_from_bot_folder(bot_folder: Path) -> dict[str, Any]:
    """
    Parse trial information from individual BOT recording folder names.

    Expected format: BOT_[date]_[slice_info]_[treatment]_[stimulation]-[trial_num]
    Examples:
    - BOT_04162024_slice2ROI1_50nMDA_burst-001
    - BOT_05242024_slice1A_ctr_single-001
    - BOT_04052024_slice2_ACh-001 (calibration)

    Parameters
    ----------
    bot_folder : Path
        Path to the BOT recording folder

    Returns
    -------
    Dict[str, Any]
        Dictionary containing trial information
    """
    trial_name = bot_folder.name

    # Parse trial name: BOT_[date]_[slice_info]_[treatment]_[stimulation]-[trial_num]
    pattern = r"BOT_(\d{8})_(.+?)_([^_]+)_([^-]+)-(\d+)"
    match = re.match(pattern, trial_name)

    if not match:
        # Handle calibration sessions (ACh, TTX) - no stimulation protocol
        calibration_pattern = r"BOT_(\d{8})_(.+?)_([^-]+)-(\d+)"
        calibration_match = re.match(calibration_pattern, trial_name)

        if calibration_match:
            date_str, slice_info, treatment, trial_num = calibration_match.groups()

            # Check if this is a real calibration session (ACh, TTX) or missing stimulation type
            if treatment in ["ACh", "TTX"] or "ACh" in treatment or "TTX" in treatment:
                stimulation = "calibration"
            else:
                # Experimental trial missing stimulation type - assume single pulse
                stimulation = "single"
                print(f"Warning: Trial {trial_name} missing stimulation type, assuming single pulse")
        else:
            raise ValueError(f"Could not parse BOT folder name: {trial_name}")
    else:
        date_str, slice_info, treatment, stimulation, trial_num = match.groups()

        # Check if stimulation field contains calibration markers
        if stimulation in ["ACh", "TTX"] or "ACh" in stimulation or "TTX" in stimulation:
            # This is actually a calibration trial, adjust treatment and stimulation
            treatment = f"{treatment}_{stimulation}_calibration"
            stimulation = "calibration"

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
        "trial_name": trial_name,
        "treatment": treatment_full,
        "stimulation": stimulation_full,
        "trial_number": trial_num,
        "is_calibration": treatment in ["ACh", "TTX"] or "ACh" in treatment or "TTX" in treatment,
    }


def convert_slice_session_to_nwbfile(slice_folder: Path, condition: str, session_times_json: dict = None) -> NWBFile:
    """
    Convert a single slice session (containing multiple BOT trials) to NWB format.

    Parameters
    ----------
    slice_folder : Path
        Path to the slice folder (e.g., 04022024slice1ROI1)
    condition : str
        Experimental condition ("UL control", "PD", "LID off")
    session_times_json : dict, optional
        Session timing data from JSON analysis

    Returns
    -------
    NWBFile
        NWB file with the converted data
    """
    # Parse slice session information
    session_info = parse_slice_session_info(slice_folder)

    # Get all BOT trial folders in this slice session
    bot_trial_folders = [f for f in slice_folder.iterdir() if f.is_dir() and f.name.startswith("BOT_")]

    if not bot_trial_folders:
        raise FileNotFoundError(f"No BOT trial folders found in {slice_folder}")

    # Sort trials for consistent ordering
    bot_trial_folders.sort()

    # Get session start time from the first trial's XML if available
    first_trial_folder = bot_trial_folders[0]
    xml_files = [
        f
        for f in first_trial_folder.iterdir()
        if f.name.endswith(".xml") and "VoltageRecording" not in f.name and "VoltageOutput" not in f.name
    ]

    if xml_files:
        # Use interface to get session start time
        temp_csv = [f for f in first_trial_folder.iterdir() if f.name.endswith("-botData.csv")][0]
        temp_interface = PrairieViewFluorescenceInterface(
            bot_csv_data_file_path=temp_csv, xml_metadata_file_path=xml_files[0]
        )
        actual_session_start_time = temp_interface.get_session_start_time()
    else:
        actual_session_start_time = None

    # Load general and session-specific metadata from YAML files
    general_metadata_path = Path(__file__).parent.parent.parent / "general_metadata.yaml"
    general_metadata = load_dict_from_file(general_metadata_path)

    session_metadata_path = Path(__file__).parent.parent.parent / "session_specific_metadata.yaml"
    session_metadata_template = load_dict_from_file(session_metadata_path)
    script_template = session_metadata_template["figure_5_acetylcholine_biosensor"]

    # Create session-specific metadata with Chicago timezone
    central_tz = ZoneInfo("America/Chicago")

    if actual_session_start_time:
        session_start_time = actual_session_start_time.replace(tzinfo=central_tz)
        timestamp = session_start_time.strftime("%Y%m%d%H%M%S")
    else:
        session_start_time = datetime.combine(session_info["session_date"], datetime.min.time()).replace(
            tzinfo=central_tz
        )
        timestamp = session_start_time.strftime("%Y%m%d")

    # Create session ID following standard pattern used in other conversion scripts
    condition_to_camel_case = {"UL control": "ULControl", "PD": "PD", "LID off": "LIDOff"}
    clean_condition = condition_to_camel_case.get(condition, condition.replace(" ", "").replace("-", ""))

    # Create session ID with ++ separators (no dashes or underscores)
    base_session_id = f"Figure5++AcetylcholineGRAB++{clean_condition}++{timestamp}"
    script_specific_id = f"Sub++{session_info['slice_folder_name']}"
    session_id = f"{base_session_id}++{script_specific_id}"

    # Handle surgery and pharmacology from template
    surgery_text = general_metadata["NWBFile"]["surgery"] + " " + script_template["NWBFile"]["surgery_addition"]
    pharmacology_text = (
        general_metadata["NWBFile"]["pharmacology"]
        + " All treatments applied during this session as described in trials table."
    )

    # Create session-specific metadata from template with runtime substitutions
    session_specific_metadata = {
        "NWBFile": {
            "session_description": f"Acetylcholine GRAB biosensor session for {condition} condition in slice {session_info['slice_info']}. Contains multiple trials with different treatments and stimulation protocols.",
            "identifier": str(uuid.uuid4()),
            "session_start_time": session_start_time,
            "session_id": session_id,
            "surgery": surgery_text,
            "pharmacology": pharmacology_text,
            "keywords": script_template["NWBFile"]["keywords"],
        },
        "Subject": {
            "subject_id": session_info["subject_id"],
            "description": script_template["Subject"]["description"].format(
                session_id=session_id, date_str=session_info["date_str"]
            ),
            "genotype": script_template["Subject"]["genotype"],
        },
    }

    # Deep merge with general metadata
    metadata = dict_deep_update(general_metadata, session_specific_metadata)

    # Create NWB file using neuroconv helper function
    nwbfile = make_nwbfile_from_metadata(metadata)

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

    # Create trials table for all BOT recordings in this session
    trials_table = TimeIntervals(
        name="trials", description="Individual BOT recording trials with different treatments and stimulation protocols"
    )

    # Add custom columns for trial metadata
    trials_table.add_column(name="trial_name", description="Name of the BOT recording trial")
    trials_table.add_column(name="treatment", description="Pharmacological treatment applied")
    trials_table.add_column(name="stimulation", description="Stimulation protocol used")
    trials_table.add_column(name="stimulus_type", description="Type of stimulation (single/burst/calibration)")
    trials_table.add_column(name="amplitude", description="Stimulation amplitude if applicable")
    trials_table.add_column(name="electrode", description="Electrode model and type")
    trials_table.add_column(name="is_calibration", description="Whether this is a calibration trial")
    trials_table.add_column(name="original_trial_number", description="Original trial number from BOT filename")

    # Process each BOT trial and add fluorescence data
    cumulative_time = 0.0
    trial_duration = 8.0  # Typical trial length in seconds

    for trial_index, bot_trial_folder in enumerate(bot_trial_folders):
        # Parse trial information
        trial_info = parse_trial_info_from_bot_folder(bot_trial_folder)

        # Find BOT CSV and XML files for this trial
        bot_csv_file = None
        xml_metadata_file = None

        for file in bot_trial_folder.iterdir():
            if file.name.endswith("-botData.csv"):
                bot_csv_file = file
            elif (
                file.name.endswith(".xml") and "VoltageRecording" not in file.name and "VoltageOutput" not in file.name
            ):
                xml_metadata_file = file

        if not bot_csv_file or not xml_metadata_file:
            print(f"Warning: Missing files in {bot_trial_folder}, skipping trial")
            continue

        # Create naming mappings for this trial (used for both fluorescence and raw data)
        treatment_to_camel_case = {
            "control": "Control",
            "50nM_dopamine": "50nMDopamine",
            "quinpirole": "Quinpirole",
            "sulpiride": "Sulpiride",
            "acetylcholine_calibration": "AcetylcholineCalibration",
            "TTX_calibration": "TTXCalibration",
        }

        stimulation_to_camel_case = {
            "single_pulse": "SinglePulse",
            "burst_stimulation": "BurstStimulation",
            "calibration": "Calibration",
        }

        # Get camelCase versions
        clean_condition = condition_to_camel_case.get(condition, condition.replace(" ", "").replace("-", ""))
        clean_treatment = treatment_to_camel_case.get(trial_info["treatment"], trial_info["treatment"].replace("_", ""))
        clean_stimulation = stimulation_to_camel_case.get(
            trial_info["stimulation"], trial_info["stimulation"].replace("_", "")
        )

        # Create base name with condition and stimulus info, avoiding duplications
        # For calibration trials, the treatment already contains "Calibration", so don't add stimulation
        if trial_info["is_calibration"]:
            base_name = f"{clean_condition}{clean_treatment}"
        else:
            base_name = f"{clean_condition}{clean_treatment}{clean_stimulation}"

        # Create acetylcholine fluorescence interface for this trial
        acetylcholine_interface = PrairieViewFluorescenceInterface(
            bot_csv_data_file_path=bot_csv_file, xml_metadata_file_path=xml_metadata_file
        )

        # Get precise trial start time and shift timestamps
        trial_start_time = acetylcholine_interface.get_session_start_time()
        time_shift = (trial_start_time - session_start_time.replace(tzinfo=None)).total_seconds()
        trial_start_shifted = cumulative_time + time_shift

        # Add trial to trials table (include actual trial number from filename in metadata)
        stimulus_desc = "None" if trial_info["is_calibration"] else "0.3 mA"
        electrode_desc = "None" if trial_info["is_calibration"] else "CBAPD75 concentric bipolar"
        actual_trial_number = int(trial_info["trial_number"])

        trials_table.add_interval(
            start_time=trial_start_shifted,
            stop_time=trial_start_shifted + trial_duration,
            trial_name=trial_info["trial_name"],
            treatment=trial_info["treatment"],
            stimulation=trial_info["stimulation"],
            stimulus_type=trial_info["stimulation"],
            amplitude=stimulus_desc,
            electrode=electrode_desc,
            is_calibration=trial_info["is_calibration"],
            original_trial_number=actual_trial_number,  # Store actual trial number from filename
        )

        # Add fluorescence data for this trial with trial-specific naming
        # Use sequential trial_index for unique naming within session, actual trial number stored in trials table
        fluorescence_suffix = f"{base_name}Trial{trial_index:03d}"
        acetylcholine_interface.add_to_nwbfile(nwbfile=nwbfile, trial_id=fluorescence_suffix)

        # Add raw imaging data for this trial
        # Add raw imaging data with custom naming
        bruker_converter = BrukerTiffSinglePlaneConverter(folder_path=bot_trial_folder)
        bruker_metadata = bruker_converter.get_metadata()

        # Customize series name in metadata before adding to NWB
        if "Ophys" in bruker_metadata and "TwoPhotonSeries" in bruker_metadata["Ophys"]:
            # Update the series names in metadata (it's a list of series metadata)
            for series_metadata in bruker_metadata["Ophys"]["TwoPhotonSeries"]:
                # Extract channel info and create better channel names
                channel_suffix = series_metadata["name"].replace("TwoPhotonSeries", "")
                if channel_suffix == "Ch2":
                    channel_name = "GRABCh"  # GRAB acetylcholine channel
                elif channel_suffix == "Dodt":
                    channel_name = "DoDT"  # Differential optical detection transmission
                else:
                    channel_name = channel_suffix.replace("Ch", "Channel")

                # Create series name: TwoPhotonSeries + condition + treatment + stimulation + channel + trial
                series_metadata["name"] = f"TwoPhotonSeries{base_name}{channel_name}Trial{trial_index:03d}"

        bruker_converter.add_to_nwbfile(nwbfile=nwbfile, metadata=bruker_metadata)

        # Update cumulative time for next trial
        cumulative_time = trial_start_shifted + trial_duration + 1.0  # 1s gap between trials

    # Add trials table to NWB file
    nwbfile.trials = trials_table

    return nwbfile


if __name__ == "__main__":
    import argparse
    import logging
    import warnings

    from tqdm import tqdm

    # Parse command line arguments
    def str_to_bool(v):
        if isinstance(v, bool):
            return v
        if v.lower() in ("yes", "true", "t", "y", "1"):
            return True
        elif v.lower() in ("no", "false", "f", "n", "0"):
            return False
        else:
            raise argparse.ArgumentTypeError("Boolean value expected.")

    parser = argparse.ArgumentParser(description="Convert Figure 5 acetylcholine biosensor data to NWB format")
    parser.add_argument(
        "--stub-test",
        type=str_to_bool,
        default=True,
        help="Process only first 2 files per condition for testing (default: True). Use --stub-test=False for full processing.",
    )

    args = parser.parse_args()
    stub_test = args.stub_test

    # Control verbose output
    verbose = False  # Set to True for detailed output

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

        # Get all slice session folders for this condition
        # Structure: condition/slice_folder (e.g., 04022024slice1ROI1)
        slice_session_folders = [f for f in condition_path.iterdir() if f.is_dir()]
        slice_session_folders.sort()

        # Apply stub_test filtering if enabled
        if stub_test:
            slice_session_folders = slice_session_folders[:2]

        # Use tqdm for progress bar
        session_iterator = tqdm(
            slice_session_folders, desc=f"Converting Figure5 AcetylcholineBiosensor {condition}", unit=" session"
        )

        for slice_folder in session_iterator:
            # Convert slice session to NWB format
            nwbfile = convert_slice_session_to_nwbfile(
                slice_folder=slice_folder,
                condition=condition,
            )

            # Create output filename using session_id
            nwbfile_path = nwb_files_dir / f"{nwbfile.session_id}.nwb"

            configure_and_write_nwbfile(nwbfile, nwbfile_path=nwbfile_path)
