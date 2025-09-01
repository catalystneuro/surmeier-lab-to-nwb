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
from datetime import datetime
from pathlib import Path
from typing import Any
from zoneinfo import ZoneInfo

from neuroconv.converters import BrukerTiffSinglePlaneConverter
from neuroconv.tools import configure_and_write_nwbfile
from neuroconv.tools.nwb_helpers import make_nwbfile_from_metadata
from neuroconv.utils import dict_deep_update, load_dict_from_file
from pynwb import NWBFile
from pynwb.device import Device
from pynwb.epoch import TimeIntervals

from surmeier_lab_to_nwb.zhai2025.conversion_scripts.conversion_utils import (
    generate_canonical_session_id,
    str_to_bool,
)
from surmeier_lab_to_nwb.zhai2025.interfaces import (
    PrairieViewFluorescenceInterface,
)
from surmeier_lab_to_nwb.zhai2025.interfaces.ophys_interfaces import (
    BrukerReferenceImagesInterface,
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
        "subject_id": f"SubjectRecordedAt{date_str.replace('-', '')}",  # Use SubjectRecordedAt pattern
    }


def parse_basic_trial_info_from_folder(bot_folder: Path) -> dict[str, str]:
    """
    Extract basic trial information from BOT folder name.

    Parameters
    ----------
    bot_folder : Path
        Path to the BOT recording folder

    Returns
    -------
    Dict[str, str]
        Dictionary containing treatment and repetition_number
    """
    trial_name = bot_folder.name

    # REGEX PATTERN BREAKDOWN:
    # BOT_\d{8}_[^_]+_(.+)-(\d+)
    #
    # BOT_           - Literal "BOT_" prefix
    # \d{8}          - Exactly 8 digits (date: MMDDYYYY)
    # _              - Literal underscore separator
    # [^_]+          - One or more non-underscore chars (slice info: slice1A, slice2ROI1, etc.)
    # _              - Literal underscore separator
    # (.+)           - Group 1: Treatment info (everything before dash, may include underscores)
    #                  Examples: "ctr", "50nMDA_burst", "ACh wash in", "sul_sinlge"
    # -              - Literal dash separator
    # (\d+)          - Group 2: Repetition number (001, 002, etc.)
    #
    # Examples:
    # BOT_04022024_slice2ROI1_ctr_single-001 → Groups: ('ctr_single', '001')
    # BOT_05242024_slice1A_ACh-001           → Groups: ('ACh', '001')
    # BOT_04022024_slice2ROI1_ACh wash in-001 → Groups: ('ACh wash in', '001')
    #
    pattern = r"BOT_\d{8}_[^_]+_(.+)-(\d+)"
    match = re.match(pattern, trial_name)

    if match:
        treatment_raw, trial_num = match.groups()

        # Clean up treatment name for special cases
        if treatment_raw == "ACh wash in":
            treatment = "ACh"
        elif "_" in treatment_raw:
            # Handle cases like "50nMDA_burst" or "ctr_single" - take everything before last underscore
            treatment = (
                "_".join(treatment_raw.split("_")[:-1])
                if treatment_raw.split("_")[-1] in ["single", "burst", "sinlge", "singl", "20Hz"]
                else treatment_raw
            )
        else:
            treatment = treatment_raw

        return {"treatment": treatment, "repetition_number": trial_num}
    else:
        raise ValueError(f"Could not parse BOT folder name: {trial_name}")


def get_stimulation_type_from_xml(bot_folder: Path) -> str:
    """
    Determine stimulation type from VoltageOutput XML file.

    Parameters
    ----------
    bot_folder : Path
        Path to the BOT recording folder

    Returns
    -------
    str
        Stimulation type: "calibration", "single_pulse", or "burst_stimulation"
    """
    voltage_xml = bot_folder / f"{bot_folder.name}_Cycle00001_VoltageOutput_001.xml"

    if not voltage_xml.exists():
        return "calibration"  # No electrical stimulation - drug-only calibration

    try:
        import xml.etree.ElementTree as ET

        tree = ET.parse(voltage_xml)
        root = tree.getroot()

        # Find the LED waveform (the actual stimulation channel)
        for waveform in root.findall("Waveform"):
            name_elem = waveform.find("Name")
            enabled_elem = waveform.find("Enabled")

            if (
                name_elem is not None
                and name_elem.text == "LED"
                and enabled_elem is not None
                and enabled_elem.text == "true"
            ):

                # Get pulse count from the waveform component
                pulse_train = waveform.find("WaveformComponent_PulseTrain")
                if pulse_train is not None:
                    pulse_count_elem = pulse_train.find("PulseCount")
                    if pulse_count_elem is not None:
                        pulse_count = int(pulse_count_elem.text)

                        if pulse_count == 1:
                            return "single_pulse"
                        elif pulse_count > 1:
                            return "burst_stimulation"

        # Fallback if we can't find LED waveform
        raise ValueError(f"Could not determine stimulation type from XML: {voltage_xml}")

    except Exception as e:
        raise ValueError(f"Error parsing VoltageOutput XML {voltage_xml}: {e}")


def parse_trial_info_from_bot_folder(bot_folder: Path) -> dict[str, Any]:
    """
    Parse complete trial information from BOT folder using both folder name and XML files.

    Parameters
    ----------
    bot_folder : Path
        Path to the BOT recording folder

    Returns
    -------
    Dict[str, Any]
        Dictionary containing complete trial information
    """
    # Get basic info from folder name
    basic_info = parse_basic_trial_info_from_folder(bot_folder)

    # Get stimulation type from actual XML files
    stimulation_type = get_stimulation_type_from_xml(bot_folder)

    # Map treatment abbreviations to full names
    treatment_mapping = {
        "ctr": "control",
        "50nMDA": "50nM_dopamine",
        "quin": "quinpirole",
        "sul": "sulpiride",
        "ACh": "acetylcholine_calibration",
        "TTX": "TTX_calibration",
    }

    treatment_full = treatment_mapping.get(basic_info["treatment"], basic_info["treatment"])

    return {
        "trial_name": bot_folder.name,
        "treatment": treatment_full,
        "stimulation": stimulation_type,
        "repetition_number": basic_info["repetition_number"],
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

    # First, instantiate all interfaces and collect their metadata with session times
    trial_interfaces_data = []

    for bot_trial_folder in bot_trial_folders:
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
            raise FileNotFoundError(f"Missing files in {bot_trial_folder}")

        # Parse trial information
        trial_info = parse_trial_info_from_bot_folder(bot_trial_folder)

        # Create interface to get session start time
        acetylcholine_interface = PrairieViewFluorescenceInterface(
            bot_csv_data_file_path=bot_csv_file, xml_metadata_file_path=xml_metadata_file
        )

        # Get session start time for this trial
        trial_start_time = acetylcholine_interface.get_session_start_time()

        # Store all the data we need for later processing
        trial_interfaces_data.append(
            {
                "bot_trial_folder": bot_trial_folder,
                "trial_info": trial_info,
                "acetylcholine_interface": acetylcholine_interface,
                "trial_start_time": trial_start_time,
            }
        )

    # Sort all trials by their actual session start time
    trial_interfaces_data.sort(key=lambda x: x["trial_start_time"])

    # Get the earliest session start time for the overall session
    earliest_trial_start_time = trial_interfaces_data[0]["trial_start_time"]
    session_start_time = earliest_trial_start_time.replace(tzinfo=central_tz)

    # Load general and session-specific metadata from YAML files
    general_metadata_path = Path(__file__).parent.parent.parent / "general_metadata.yaml"
    general_metadata = load_dict_from_file(general_metadata_path)

    session_metadata_path = Path(__file__).parent.parent.parent / "session_specific_metadata.yaml"
    session_metadata_template = load_dict_from_file(session_metadata_path)
    script_template = session_metadata_template["figure_5_acetylcholine_biosensor"]

    # Create session-specific metadata with Chicago timezone
    central_tz = ZoneInfo("America/Chicago")

    # Map Figure 5 specific conditions to explicit parameters
    condition_mapping = {"UL control": "control", "PD": "6-OHDA", "LID off": "off-state"}
    standardized_condition = condition_mapping.get(condition, condition)

    # Map condition to explicit state
    if standardized_condition == "control":
        state = "CTRL"
    elif standardized_condition == "6-OHDA":
        state = "PD"
    elif standardized_condition == "off-state":
        state = "OFF"
    else:
        raise ValueError(f"Unknown condition: {standardized_condition}")

    # Create canonical session ID with explicit parameters
    session_id = generate_canonical_session_id(
        fig="F5",
        meas_comp="AChFP",  # GRAB-ACh fiber-photometry
        cell_type="pan",  # Non cell-specific bulk signal
        state=state,
        pharm="none",  # No acute pharmacology
        geno="WT",  # Wild-type
        timestamp=session_start_time.strftime("%Y%m%d%H%M%S"),
    )

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
    trials_table.add_column(name="amplitude", description="Stimulation amplitude if applicable")
    trials_table.add_column(name="electrode", description="Electrode model and type")
    trials_table.add_column(name="original_trial_number", description="Original trial number from BOT filename")

    # Process each BOT trial and add fluorescence data (now in temporal order)
    cumulative_time = 0.0
    trial_duration = 8.0  # Typical trial length in seconds

    for trial_index, trial_data in enumerate(trial_interfaces_data):
        # Extract pre-computed data
        bot_trial_folder = trial_data["bot_trial_folder"]
        trial_info = trial_data["trial_info"]
        acetylcholine_interface = trial_data["acetylcholine_interface"]

        # Create naming mappings for this trial (used for both fluorescence and raw data)
        condition_to_camel_case = {
            "UL control": "ULControl",
            "PD": "PD",
            "LID off": "LIDOff",
        }

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
        if trial_info["stimulation"] == "calibration":
            base_name = f"{clean_condition}{clean_treatment}"
        else:
            base_name = f"{clean_condition}{clean_treatment}{clean_stimulation}"

        # Get precise trial start time and shift timestamps (interface already created)
        trial_start_time = trial_data["trial_start_time"]
        time_shift = (trial_start_time - session_start_time.replace(tzinfo=None)).total_seconds()
        trial_start_shifted = cumulative_time + time_shift

        # Add trial to trials table (include actual trial number from filename in metadata)
        stimulus_desc = "None" if trial_info["stimulation"] == "calibration" else "0.3 mA"
        electrode_desc = "None" if trial_info["stimulation"] == "calibration" else "CBAPD75 concentric bipolar"
        actual_trial_number = int(trial_info["repetition_number"])

        trials_table.add_interval(
            start_time=trial_start_shifted,
            stop_time=trial_start_shifted + trial_duration,
            trial_name=trial_info["trial_name"],
            treatment=trial_info["treatment"],
            stimulation=trial_info["stimulation"],
            amplitude=stimulus_desc,
            electrode=electrode_desc,
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

        # Add reference images if they exist for this trial
        references_folder = bot_trial_folder / "References"
        # Create container name specific to this trial
        ref_container_name = f"BackgroundReferences{base_name}Trial{trial_index:03d}"
        reference_interface = BrukerReferenceImagesInterface(
            references_folder_path=references_folder, container_name=ref_container_name
        )
        reference_interface.add_to_nwbfile(nwbfile=nwbfile)

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
