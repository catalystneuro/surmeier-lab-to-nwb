"""
Figure 7 OxoM Dendritic Excitability Conversion Script - Zhai et al. 2025
=====================================================================

This script converts OxoM-induced dendritic excitability data from Figure 7 of
Zhai et al. 2025 into NWB (Neurodata Without Borders) format. The data combines
patch-clamp recordings with simultaneous two-photon imaging to examine the effects
of muscarinic agonist oxotremorine-M on dendritic excitability.

For detailed experimental context, protocols, data structure, and analysis methods,
see: /src/surmeier_lab_to_nwb/zhai2025/conversion_notes_folder/figure_7_conversion_notes.md
"""

import re
import uuid
import warnings
from pathlib import Path
from typing import Any, Dict

from neuroconv.tools import configure_and_write_nwbfile
from neuroconv.utils import dict_deep_update, load_dict_from_file
from pynwb import NWBFile
from pynwb.file import Subject

from surmeier_lab_to_nwb.zhai2025.interfaces import (
    DendriticTrialsInterface,
    PrairieViewCurrentClampInterface,
    PrairieViewLineScanInterface,
)


def parse_session_info_from_folder_name(recording_folder: Path) -> Dict[str, Any]:
    """
    Parse essential recording information from Figure 7E oxotremorine-M dendritic excitability recording folder names.
    Session start time comes from XML metadata.

    Expected folder name format: [MMDDYYYY]_Cell[N]_[location][location_number]_trio-[trial_number]
    (e.g., 03012020_Cell1_dist1_trio-001, 03012020_Cell1_prox1_trio-004)

    The session folder structure is: [MMDD][letter][digit]/[MMDDYYYY]_Cell[N]_[location][location_number]_trio-[trial]/
    - [MMDD][letter][digit]: Date and animal identifier (e.g., 0301a, 0301b)
    - [MMDDYYYY]_Cell[N]_[location][location_number]_trio-[trial]: Recording with dendritic location and trial info

    Protocol: Three 2 nA current injections, 2 ms each, at ~50 Hz for dendritic calcium imaging
    oxoM Protocol Convention:
    - trio-001, -002, -003: BEFORE oxotremorine-M application (baseline measurements)
    - trio-004, -005, -006: AFTER oxotremorine-M application (drug response)
    - trio-007, -008, -009: AFTER oxotremorine-M application (additional recordings if present)

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
    # Format: [MMDD][letter][digit] (e.g., 0301a, 0301b, 0302a)
    session_folder_name = session_folder.name

    # Extract date and animal identifier
    # Pattern: 4 digits (MMDD) followed by a letter and optionally a digit
    session_pattern = r"(\d{4})([a-z])(\d*)"
    session_match = re.match(session_pattern, session_folder_name)

    if not session_match:
        raise ValueError(f"Could not parse session folder name: {session_folder_name}")

    animal_letter = session_match.group(2)  # letter (a, b, c, etc.)
    animal_digit = session_match.group(3) or ""  # optional digit
    animal_id = f"{animal_letter}{animal_digit}" if animal_digit else animal_letter

    # Note: Date will be extracted from XML session start time, not from folder name

    # Parse recording folder name for dendritic recording info
    # Format: [MMDDYYYY]_Cell[N]_[location][location_number]_trio-[trial_number]
    # Example: 03012020_Cell1_dist1_trio-001, 03012020_Cell1_prox1_trio-004
    pattern = r"(\d{8})_Cell(\d+)_([a-z]+)(\d+)_trio-(\d+)"
    match = re.match(pattern, folder_name)

    if not match:
        raise ValueError(f"Could not parse recording folder name: {folder_name}")

    date_from_folder = match.group(1)  # MMDDYYYY (for reference, but we'll use XML time)
    cell_number = match.group(2)
    location_type = match.group(3)  # "dist" or "prox"
    location_number = match.group(4)
    trial_number = match.group(5)

    # Determine recording phase based on trial number
    trial_num = int(trial_number)
    if trial_num <= 3:
        recording_phase = "baseline"
        phase_description = "before oxotremorine-M application"
    else:
        recording_phase = "post_oxoM"
        phase_description = "after oxotremorine-M application"

    # Create location identifier
    location_id = f"{location_type}{location_number}"

    # Determine approximate distance from soma based on location type
    if location_type == "prox":
        approximate_distance_um = 40  # proximal dendrites ~40 μm from soma
    elif location_type == "dist":
        approximate_distance_um = 90  # distal dendrites ~90 μm from soma
    else:
        approximate_distance_um = None

    # Create full location description
    location_full = "Proximal" if location_type == "prox" else "Distal"
    location_description = f"{location_full} dendrite {location_number}"

    return {
        "cell_number": cell_number,
        "animal_id": animal_id,
        "location_type": location_type,
        "location_number": location_number,
        "location_id": location_id,
        "location_full": location_full,
        "location_description": location_description,
        "trial_number": trial_number,
        "recording_phase": recording_phase,
        "phase_description": phase_description,
        "approximate_distance_um": approximate_distance_um,
        "recording_folder_name": folder_name,
    }


def convert_session_to_nwbfile(session_folder_path: Path, genotype: str, verbose: bool = False) -> NWBFile:
    """
    Convert a single session of Figure 7E oxotremorine-M dendritic excitability data to NWB format with time alignment.

    Parameters
    ----------
    session_folder_path : Path
        Path to the session folder containing dendritic recordings
    genotype : str
        Mouse genotype ("WT" for wildtype or "KO" for CDGI knockout)
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

    The oxoM protocol follows a paired before/after design:
    - Baseline recordings (trio-001 to -003): Measure dendritic excitability in normal ACSF
    - Drug application: Add oxotremorine-M to test M1 receptor responsiveness
    - Post-drug recordings (trio-004 to -006/009): Measure dendritic excitability after oxo-M
    """

    # Define genotype-specific variables for metadata
    genotype_full = "CalDAG-GEFI knockout" if genotype == "KO" else "wildtype control"
    genotype_bg = (
        "CalDAG-GEFI knockout on Drd1-Tdtomato bacterial artificial chromosome (BAC) transgenic background"
        if genotype == "KO"
        else "Drd1-Tdtomato bacterial artificial chromosome (BAC) transgenic"
    )

    # Find all recording folders for this session
    recording_folders = [f for f in session_folder_path.iterdir() if f.is_dir()]
    recording_folders.sort()

    if not recording_folders:
        raise ValueError(f"No recording folders found in session folder: {session_folder_path}")

    # Parse session information from first recording folder (all should have same session info)
    first_recording_info = parse_session_info_from_folder_name(recording_folders[0])
    session_info = {
        "cell_number": first_recording_info["cell_number"],
        "animal_id": first_recording_info["animal_id"],
    }

    # Group recordings by location and phase for analysis
    recordings_by_location = {}
    baseline_recordings = []
    post_oxoM_recordings = []

    # Calculate recording IDs, session start times, and create interface mappings
    ophys_session_start_times = []  # (ophys_time, recording_folder, recording_id)
    intracellular_session_start_times = []  # (intracellular_time, recording_folder, recording_id)
    recording_id_to_info = {}
    recording_id_to_folder = {}
    recording_id_to_location_id = {}
    t_starts = {}  # t_starts[recording_id][interface] = t_start_offset

    for recording_folder in recording_folders:
        # Parse recording information using unified function
        recording_info = parse_session_info_from_folder_name(recording_folder)

        # Create unique recording ID that includes location and trial info
        recording_id = (
            f"Cell{recording_info['cell_number']}_"
            f"{recording_info['location_id']}_"
            f"trio{recording_info['trial_number']}"
        )

        # Find XML files for this recording
        main_xml_file = recording_folder / f"{recording_folder.name}.xml"
        if not main_xml_file.exists():
            raise FileNotFoundError(f"Expected main XML file does not exist: {main_xml_file}")

        electrophysiology_xml_file = recording_folder / f"{recording_folder.name}_Cycle00001_VoltageRecording_001.xml"
        if not electrophysiology_xml_file.exists():
            raise FileNotFoundError(f"Expected electrophysiology XML file does not exist: {electrophysiology_xml_file}")

        # Get session start times from both sources
        ophys_session_start_time = PrairieViewLineScanInterface.get_session_start_time_from_file(main_xml_file)
        if ophys_session_start_time is None:
            raise ValueError(f"Could not extract ophys session start time from {main_xml_file}")

        intracellular_session_start_time = PrairieViewCurrentClampInterface.get_session_start_time_from_file(
            electrophysiology_xml_file
        )
        if intracellular_session_start_time is None:
            raise ValueError(f"Could not extract intracellular session start time from {electrophysiology_xml_file}")

        # Store mappings
        location_id = f"Cell{recording_info['cell_number']}{recording_info['location_full']}Dendrite{recording_info['location_number']}"
        recording_id_to_info[recording_id] = recording_info
        recording_id_to_folder[recording_id] = recording_folder
        recording_id_to_location_id[recording_id] = location_id

        # Store session start times separately by interface type
        ophys_session_start_times.append((ophys_session_start_time, recording_folder, recording_id))
        intracellular_session_start_times.append((intracellular_session_start_time, recording_folder, recording_id))

        # Group by location and phase for organization
        location_key = recording_info["location_id"]
        if location_key not in recordings_by_location:
            recordings_by_location[location_key] = {"baseline": [], "post_oxoM": []}
        recordings_by_location[location_key][recording_info["recording_phase"]].append(recording_id)

        if recording_info["recording_phase"] == "baseline":
            baseline_recordings.append(recording_id)
        else:
            post_oxoM_recordings.append(recording_id)

    if not ophys_session_start_times or not intracellular_session_start_times:
        raise ValueError(f"No valid recordings found in session folder: {session_folder_path}")

    # Find the earliest session start time from each interface type
    earliest_ophys_time = min(ophys_session_start_times, key=lambda x: x[0])[0]
    earliest_intracellular_time = min(intracellular_session_start_times, key=lambda x: x[0])[0]

    # Overall session start time is the earliest across all interfaces
    session_start_time = min(earliest_ophys_time, earliest_intracellular_time)

    # Determine which interface had the earliest time
    if session_start_time == earliest_ophys_time:
        earliest_folder = next(
            folder for start_time, folder, _ in ophys_session_start_times if start_time == session_start_time
        )
        earliest_interface = "line_scan_ophys"
    else:
        earliest_folder = next(
            folder for start_time, folder, _ in intracellular_session_start_times if start_time == session_start_time
        )
        earliest_interface = "intracellular_electrophysiology"

    # Calculate t_start offsets for temporal alignment with interface-specific timing
    for ophys_time, folder, recording_id in ophys_session_start_times:
        intracellular_time = next(time for time, _, rid in intracellular_session_start_times if rid == recording_id)

        # Calculate offsets relative to overall session start time
        ophys_t_start = (ophys_time - session_start_time).total_seconds()
        intracellular_t_start = (intracellular_time - session_start_time).total_seconds()

        # Initialize t_starts for this recording_id with interface-specific timing
        t_starts[recording_id] = {
            "intracellular": intracellular_t_start,
            "line_scan_structural_channel": ophys_t_start,  # Ch1/Alexa568 line scan uses ophys timing
            "line_scan_calcium_channel": ophys_t_start,  # Ch2/Fluo4 line scan uses ophys timing
        }

    # Extract date from actual session start time and update session info
    session_date_str = session_start_time.strftime("%Y-%m-%d")

    # Create session ID following new pattern
    base_session_id = f"figure7_DendriticExcitability_oxoM_{genotype}_{session_start_time.strftime('%Y%m%d_%H%M%S')}"
    script_specific_id = f"Animal{session_info['animal_id']}_Cell{session_info['cell_number']}"
    session_id = f"{base_session_id}_{script_specific_id}"

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
    script_template = session_metadata_template["figure_7_oxoM_dendritic_excitability"]

    # Handle conditional pharmacology for oxotremorine-M treatment
    pharmacology_addition = ""
    if "pharmacology_conditions" in script_template["NWBFile"]:
        if "oxotremorine_M" in script_template["NWBFile"]["pharmacology_conditions"]:
            pharmacology_addition = " " + script_template["NWBFile"]["pharmacology_conditions"]["oxotremorine_M"]

    # Create session-specific metadata from template with runtime substitutions
    genotype_description = "CDGI knockout" if genotype == "KO" else "wildtype"

    session_specific_metadata = {
        "NWBFile": {
            "session_description": script_template["NWBFile"]["session_description"].format(
                genotype_description=genotype_description,
                animal_id=session_info["animal_id"],
                cell_number=session_info["cell_number"],
                date_str=session_info["date_str"],
            ),
            "identifier": str(uuid.uuid4()),
            "session_start_time": session_start_time,
            "session_id": session_info["session_id"],
            "pharmacology": general_metadata["NWBFile"]["pharmacology"] + pharmacology_addition,
            "keywords": script_template["NWBFile"]["keywords"],
        },
        "Subject": {
            "subject_id": f"CDGI_{genotype}_oxoM_mouse_{session_info['animal_id']}",
            "description": script_template["Subject"]["description"].format(
                genotype_full=genotype_full,
                animal_id=session_info["animal_id"],
                cell_number=session_info["cell_number"],
                date_str=session_info["date_str"],
            ),
            "genotype": script_template["Subject"]["genotype"] if genotype == "KO" else "Wild-type CDGI",
        },
    }

    # Merge general metadata with session-specific metadata
    metadata = dict_deep_update(general_metadata, session_specific_metadata)

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
        surgery=metadata["NWBFile"]["surgery"],
        pharmacology=metadata["NWBFile"]["pharmacology"],
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

    # Add custom columns to intracellular recording table for dendritic experiment annotations
    intracellular_recording_table = nwbfile.get_intracellular_recordings()
    intracellular_recording_table.add_column(
        name="recording_phase",
        description="Phase of oxotremorine-M protocol (baseline or post_oxoM)",
    )
    intracellular_recording_table.add_column(
        name="location_type", description="Dendritic location type (prox=proximal, dist=distal)"
    )
    intracellular_recording_table.add_column(
        name="location_id", description="Specific dendritic location identifier (e.g., dist1, prox2)"
    )
    intracellular_recording_table.add_column(
        name="approximate_distance_um", description="Approximate distance from soma in micrometers"
    )
    intracellular_recording_table.add_column(name="trial_number", description="Trial number within location and phase")
    intracellular_recording_table.add_column(
        name="animal_id", description="Animal identifier for this experimental session"
    )
    intracellular_recording_table.add_column(
        name="genotype_comparison", description="Genotype for M1 receptor responsiveness study (WT or KO)"
    )

    # Data structures for tracking icephys table indices
    recording_indices = []  # Store all intracellular recording indices
    recording_to_metadata = {}  # Map recording index to metadata for table building
    location_to_recording_indices = {}  # Group recordings by location for repetitions table
    sequential_recording_indices = []  # Store sequential recording indices

    # Process each recording using the calculated recording IDs and temporal alignment
    for recording_id, recording_folder in recording_id_to_folder.items():
        recording_info = recording_id_to_info[recording_id]
        location_id = recording_id_to_location_id[recording_id]
        main_xml_file = recording_folder / f"{recording_folder.name}.xml"
        electrophysiology_xml_file = recording_folder / f"{recording_folder.name}_Cycle00001_VoltageRecording_001.xml"

        # Create interfaces for the two known channels
        structural_ophys_key = f"PrairieViewLineScan{recording_id}Alexa568"
        calcium_ophys_key = f"PrairieViewLineScan{recording_id}Fluo4"

        structural_interface = PrairieViewLineScanInterface(
            file_path=main_xml_file,
            channel_name="Ch1",
            ophys_metadata_key=structural_ophys_key,
        )

        calcium_interface = PrairieViewLineScanInterface(
            file_path=main_xml_file,
            channel_name="Ch2",
            ophys_metadata_key=calcium_ophys_key,
        )

        # Apply temporal alignment offsets using precise mapping with descriptive interface names
        structural_interface.set_aligned_starting_time(t_starts[recording_id]["line_scan_structural_channel"])
        calcium_interface.set_aligned_starting_time(t_starts[recording_id]["line_scan_calcium_channel"])

        # Create intracellular recording interface
        icephys_metadata_key = f"PrairieView{recording_id}"
        intracellular_interface = PrairieViewCurrentClampInterface(
            file_path=electrophysiology_xml_file,
            icephys_metadata_key=icephys_metadata_key,
        )

        # Apply temporal alignment offset
        intracellular_interface.set_aligned_starting_time(t_starts[recording_id]["intracellular"])

        # Get and update intracellular metadata
        intracellular_metadata = intracellular_interface.get_metadata()

        # Update electrode description for CDGI dendritic recording
        electrode_name = f"IntracellularElectrode{location_id}"
        intracellular_metadata["Icephys"]["IntracellularElectrodes"][icephys_metadata_key].update(
            {
                "name": electrode_name,
                "description": (
                    f"Whole-cell patch clamp electrode recording from {genotype_description} iSPN dendrite "
                    f"in the dorsolateral striatum - Animal {session_info['animal_id']}, "
                    f"Cell {session_info['cell_number']} - Location {recording_info['location_id']} "
                    f"(~{recording_info['approximate_distance_um']}μm from soma) - oxotremorine-M protocol"
                ),
                "cell_id": session_info["cell_number"],
                "animal_id": session_info["animal_id"],
                "location": f"dendrite (~{recording_info['approximate_distance_um']}μm from soma)",
                "slice": "280 μm sagittal brain slice from dorsolateral striatum (Paper Methods: 'Sagittal sections (280 μm thick) were cut using a Leica VT1200 vibratome')",
                "seal": "Gigaohm seal (whole-cell configuration) (Paper Methods: patch clamp methodology, whole-cell configuration implied)",
                "resistance": "3-5 MΩ (borosilicate glass pipette) (Protocol: Ex_vivo_mouse_brain_patch_clamp_recordings: 'Pipette resistance must be of 3 to 5 megaohms')",
                "filtering": "2 kHz low-pass filter (Paper Methods: 'signals were filtered at 2 kHz and digitized at 10 kHz')",
                "initial_access_resistance": "<20 MΩ (typical for whole-cell recordings) (Standard electrophysiology practice for healthy whole-cell recordings)",
            }
        )

        # Update current clamp series metadata with dendritic information
        series_name = f"CurrentClamp{recording_info['location_id']}_trio{recording_info['trial_number']}Series"
        intracellular_metadata["Icephys"]["CurrentClampSeries"][icephys_metadata_key].update(
            {
                "name": series_name,
                "description": (
                    f"Current clamp recording from {genotype_description} iSPN dendrite - "
                    f"Animal {session_info['animal_id']}, Cell {session_info['cell_number']} - "
                    f"Location {recording_info['location_id']} (~{recording_info['approximate_distance_um']}μm) - "
                    f"{recording_info['phase_description']} - trio {recording_info['trial_number']}"
                ),
            }
        )

        # Add intracellular data to NWB file
        intracellular_interface.add_to_nwbfile(nwbfile=nwbfile, metadata=intracellular_metadata)

        # Add intracellular recording to icephys table with custom annotations
        current_clamp_series = nwbfile.acquisition[series_name]

        # Add intracellular recording entry with enhanced metadata
        recording_index = nwbfile.add_intracellular_recording(
            electrode=current_clamp_series.electrode,
            response=current_clamp_series,
            recording_phase=recording_info["recording_phase"],
            location_type=recording_info["location_type"],
            location_id=recording_info["location_id"],
            approximate_distance_um=recording_info["approximate_distance_um"],
            trial_number=int(recording_info["trial_number"]),
            animal_id=session_info["animal_id"],
            genotype_comparison=genotype,
        )

        # Track recording index and metadata for table building
        recording_indices.append(recording_index)
        recording_to_metadata[recording_index] = {
            "recording_id": recording_id,
            "location_id": location_id,
            "recording_info": recording_info,
            "series_name": series_name,
        }

        # Group recordings by location for repetitions table
        if location_id not in location_to_recording_indices:
            location_to_recording_indices[location_id] = []
        location_to_recording_indices[location_id].append(recording_index)

        # Process structural channel (Ch1/Alexa568)
        structural_metadata = structural_interface.get_metadata()
        # Apply fluorophore-specific metadata based on experimental knowledge
        structural_metadata["Devices"][structural_ophys_key]["name"] = "BrukerUltima"
        structural_metadata["Devices"][structural_ophys_key][
            "description"
        ] = "Bruker Ultima two-photon microscope for line scan imaging. 810 nm excitation laser (Chameleon Ultra II, Coherent). Signals filtered at 2 kHz and digitized at 10 kHz."
        structural_metadata["Ophys"]["ImagingPlanes"][structural_ophys_key][
            "name"
        ] = f"ImagingPlane{recording_id}Alexa568"
        structural_metadata["Ophys"]["ImagingPlanes"][structural_ophys_key][
            "description"
        ] = f"Line scan imaging plane for {location_id} using Alexa Fluor 568 structural dye. Line scan parameters: 64 pixels per line, 10 μs dwell time, ~640 μs per line."
        structural_metadata["Ophys"]["ImagingPlanes"][structural_ophys_key]["indicator"] = "Alexa Fluor 568"
        structural_metadata["Ophys"]["PlaneSegmentation"][structural_ophys_key][
            "name"
        ] = f"PlaneSegmentation{recording_id}Alexa568"
        structural_metadata["Ophys"]["PlaneSegmentation"][structural_ophys_key][
            "description"
        ] = f"Line scan ROI segmentation for {location_id} structural imaging. Detected by Hamamatsu R3982 side-on PMT (580-620 nm)."
        structural_metadata["Ophys"]["RoiResponseSeries"][structural_ophys_key][
            "name"
        ] = f"RoiResponseSeriesAlexa568{recording_id}"
        structural_metadata["Ophys"]["RoiResponseSeries"][structural_ophys_key][
            "description"
        ] = f"Structural reference fluorescence from Alexa Fluor 568 hydrazide (50 μM) - {location_id}. Ca2+-insensitive dye to visualize dendrites."
        structural_metadata["Acquisition"]["SourceImages"][structural_ophys_key][
            "name"
        ] = f"ImageAlexa568{recording_id}"
        structural_metadata["Acquisition"]["SourceImages"][structural_ophys_key][
            "description"
        ] = f"Source image for Alexa Fluor 568 structural reference - {location_id}. Field of view with scan line overlay."
        structural_metadata["TimeSeries"][structural_ophys_key]["name"] = f"TimeSeriesLineScanRawAlexa568{recording_id}"
        structural_metadata["TimeSeries"][structural_ophys_key][
            "description"
        ] = f"Line scan raw data for Alexa Fluor 568 structural reference - {location_id}. Typical acquisition: 2500 lines (time points)."

        # Add structural data to NWB file
        structural_interface.add_to_nwbfile(nwbfile=nwbfile, metadata=structural_metadata)

        # Process calcium channel (Ch2/Fluo4)
        calcium_metadata = calcium_interface.get_metadata()
        # Apply fluorophore-specific metadata based on experimental knowledge
        calcium_metadata["Devices"][calcium_ophys_key]["name"] = "BrukerUltima"
        calcium_metadata["Devices"][calcium_ophys_key][
            "description"
        ] = "Bruker Ultima two-photon microscope for line scan imaging. 810 nm excitation laser (Chameleon Ultra II, Coherent). Signals filtered at 2 kHz and digitized at 10 kHz."
        ophys_metadata = calcium_metadata["Ophys"]
        ophys_metadata["ImagingPlanes"][calcium_ophys_key]["name"] = f"ImagingPlane{recording_id}Fluo4"
        ophys_metadata["ImagingPlanes"][calcium_ophys_key][
            "description"
        ] = f"Line scan imaging plane for {location_id} using Fluo-4 calcium indicator. Line scan parameters: 64 pixels per line, 10 μs dwell time, ~640 μs per line. Temporal resolution: ~1.6 seconds for 2500 lines."
        ophys_metadata["ImagingPlanes"][calcium_ophys_key]["indicator"] = "Fluo-4"
        ophys_metadata["PlaneSegmentation"][calcium_ophys_key]["name"] = f"PlaneSegmentation{recording_id}Fluo4"
        ophys_metadata["PlaneSegmentation"][calcium_ophys_key][
            "description"
        ] = f"Line scan ROI segmentation for {location_id} calcium imaging. Detected by Hamamatsu H7422P-40 GaAsP PMT (490-560 nm)."
        ophys_metadata["RoiResponseSeries"][calcium_ophys_key]["name"] = f"RoiResponseSeriesFluo4{recording_id}"
        ophys_metadata["RoiResponseSeries"][calcium_ophys_key][
            "description"
        ] = f"Calcium fluorescence from Fluo-4 (100 μM) - {location_id}. Ca2+-sensitive dye for measuring back-propagating action potential-evoked calcium transients. Magnitude serves as surrogate estimate of dendritic depolarization extent."
        calcium_metadata["Acquisition"]["SourceImages"][calcium_ophys_key]["name"] = f"ImageFluo4{recording_id}"
        calcium_metadata["Acquisition"]["SourceImages"][calcium_ophys_key][
            "description"
        ] = f"Source image for Fluo-4 calcium indicator - {location_id}. Field of view with scan line overlay."
        calcium_metadata["TimeSeries"][calcium_ophys_key]["name"] = f"TimeSeriesLineScanRawFluo4{recording_id}"
        calcium_metadata["TimeSeries"][calcium_ophys_key][
            "description"
        ] = f"Line scan raw data for Fluo-4 calcium indicator - {location_id}. Typical acquisition: 2500 lines (time points). Kymograph structure: (C, T, X) where C=channels, T=time/lines, X=pixels along scan line."

        # Add calcium data to NWB file
        calcium_interface.add_to_nwbfile(nwbfile=nwbfile, metadata=calcium_metadata)

    # Build icephys table hierarchical structure following PyNWB best practices

    # Step 1: Build simultaneous recordings (each location/trial is its own simultaneous group)
    simultaneous_recording_indices = []
    for recording_index in recording_indices:
        simultaneous_index = nwbfile.add_icephys_simultaneous_recording(
            recordings=[recording_index]  # Each dendritic recording is its own simultaneous group
        )
        simultaneous_recording_indices.append(simultaneous_index)

    # Step 2: Build sequential recording (group all dendritic recordings for this cell)
    sequential_index = nwbfile.add_icephys_sequential_recording(
        simultaneous_recordings=simultaneous_recording_indices,
        stimulus_type="dendritic_excitability_oxotremorine_M_protocol_CDGI",
    )
    sequential_recording_indices.append(sequential_index)

    # Step 3: Build repetitions table (for this single cell, it's just one repetition)
    repetitions_table = nwbfile.get_icephys_repetitions()
    repetitions_table.add_column(name="cell_number", description="Cell number identifier for this recording session")
    repetitions_table.add_column(name="animal_id", description="Animal identifier for this experimental session")

    repetition_index = nwbfile.add_icephys_repetition(
        sequential_recordings=sequential_recording_indices,
        cell_number=int(session_info["cell_number"]),
        animal_id=session_info["animal_id"],
    )

    # Step 4: Build experimental conditions table (group by genotype)
    experimental_conditions_table = nwbfile.get_icephys_experimental_conditions()
    experimental_conditions_table.add_column(name="genotype", description="Mouse genotype for CDGI knockout study")

    experimental_condition_index = nwbfile.add_icephys_experimental_condition(
        repetitions=[repetition_index], genotype=genotype
    )

    # Use utility function to add trials table with proper chronological ordering
    # Add trials table using interface
    trials_interface = DendriticTrialsInterface(
        recording_indices=recording_indices, recording_to_metadata=recording_to_metadata, t_starts=t_starts
    )
    trials_interface.add_to_nwbfile(nwbfile, verbose=verbose)

    return nwbfile


if __name__ == "__main__":
    import logging

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
    base_path = Path("./link_to_raw_data/Figure 7/oxoM on DE")

    # Create nwb_files directory at root level
    root_dir = Path(__file__).parent.parent.parent.parent.parent.parent  # Go up to repo root
    nwb_files_dir = root_dir / "nwb_files" / "dendritic_excitability" / "figure_7"
    nwb_files_dir.mkdir(parents=True, exist_ok=True)

    # Figure 7 oxotremorine-M dendritic excitability genotypes
    genotypes = ["WT", "KO"]

    for genotype in genotypes:
        genotype_path = base_path / genotype

        if not genotype_path.exists():
            continue

        # Get all session folders (each session = one animal/cell)
        session_folders = [f for f in genotype_path.iterdir() if f.is_dir()]
        session_folders.sort()

        # Apply stub_test filtering if enabled
        if stub_test:
            session_folders = session_folders[:2]

        if not session_folders:
            continue

        # Process each session folder with progress bar
        session_iterator = (
            tqdm(
                session_folders,
                desc=f"Converting Figure7 OxoMDendriticExcitability {genotype}",
                unit=" session",
                position=0,
                leave=True,
                disable=not verbose,
            )
            if not verbose
            else session_folders
        )

        for session_folder in session_iterator:

            # Convert session data to NWB format with time alignment
            nwbfile = convert_session_to_nwbfile(
                session_folder_path=session_folder,
                genotype=genotype,
                verbose=verbose,
            )

            # Create output filename
            nwbfile_path = nwb_files_dir / f"figure7E_oxoM_dendritic_excitability_{genotype}_{session_folder.name}.nwb"

            # Write NWB file
            configure_and_write_nwbfile(nwbfile, nwbfile_path=nwbfile_path)
