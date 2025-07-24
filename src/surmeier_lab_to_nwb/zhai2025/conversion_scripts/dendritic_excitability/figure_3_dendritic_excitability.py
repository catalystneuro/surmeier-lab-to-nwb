"""
Figure 3 Dendritic Excitability Conversion Script - Zhai et al. 2025
================================================================

This script converts dendritic excitability data from Figure 3 of Zhai et al. 2025
into NWB (Neurodata Without Borders) format. The data combines patch-clamp recordings
from indirect pathway SPNs (iSPNs) with simultaneous two-photon imaging to examine
dendritic excitability changes in a mouse model of Parkinson's disease and
levodopa-induced dyskinesia.

For detailed experimental context, protocols, data structure, and analysis methods,
see: /src/surmeier_lab_to_nwb/zhai2025/conversion_notes_folder/figure_3_conversion_notes.md
"""

import re
import warnings
from datetime import datetime
from pathlib import Path
from typing import Any, Dict

from neuroconv.tools import configure_and_write_nwbfile
from neuroconv.tools.nwb_helpers import make_nwbfile_from_metadata
from neuroconv.utils import dict_deep_update, load_dict_from_file
from pynwb import NWBFile

from surmeier_lab_to_nwb.zhai2025.conversion_scripts.conversion_utils import (
    FOLDER_TO_PAPER_CONDITION,
    format_condition,
    str_to_bool,
)
from surmeier_lab_to_nwb.zhai2025.conversion_scripts.dendritic_excitability.dendritic_excitability_utils import (
    build_dendritic_icephys_table_structure,
)
from surmeier_lab_to_nwb.zhai2025.interfaces import (
    DendriticTrialsInterface,
    PrairieViewCurrentClampInterface,
    PrairieViewLineScanInterface,
)
from surmeier_lab_to_nwb.zhai2025.interfaces.ophys_interfaces import (
    BrukerReferenceImagesInterface,
)


def parse_session_info_from_folder_name(recording_folder: Path) -> Dict[str, Any]:
    """
    Parse essential recording information from dendritic excitability recording folder names.
    Only extracts location and trial information - session start time comes from XML metadata.

    Expected folder name format: [date]_Cell[cell_number]_[location][location_number][variant]_[experiment_type]-[trial_number]
    Examples:
    - 05232016_Cell1_dist1_trio-001 (MMDDYYYY format)
    - 20160523_Cell1_dist1real_trio-001 (YYYYMMDD format with variant)

    Date formats:
    - MMDDYYYY: Month, day, year (e.g., 05232016 = May 23, 2016)
    - YYYYMMDD: Year, month, day (e.g., 20160523 = May 23, 2016)
    Format is auto-detected based on first digit (if '2', assumes YYYYMMDD format)

    Experiment types:
    - trio: Three current injections protocol (main experimental condition)
    - bsl: Baseline recording protocol (control/reference measurements)

    Variant:
    - real: Optional suffix indicating a re-recording or alternate version of the same location
    - (empty): Standard recording without variant

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

    # Parse using regex pattern for dendritic recordings
    # Handle variations like "dist1real", "realprox1" and different experiment types like "trio" or "bsl"
    # Some folders may not have experiment type, make it optional (with or without underscore)
    # Handle case variations like "Cell" vs "cell"
    if folder_name.startswith("00"):
        folder_name = folder_name[
            1:
        ]  # There is one case,  LID on-state/0718a/007182016_cell1_dist1_trio-001, where there is an extra leading zero

    if "_sul" in folder_name:
        folder_name = folder_name.replace("_sul", "")  # Remove "_sul" suffix if present

    if "proxt" in folder_name:
        folder_name = folder_name.replace("proxt", "prox")  # LID on-state with sul/0508c has this typo

    # Parse dendritic excitability recording folder names
    # Format: [date]_Cell[cell_number]_[location][location_number][variant]_[experiment_type]-[trial_number]
    pattern = (
        r"(?P<date>\d{8})"  # Date (8 digits): YYYYMMDD or MMDDYYYY
        r"_[Cc]ell(?P<cell_number>\d+)"  # Cell number: Cell1, cell2, etc.
        r"_(?:(?P<variant_prefix>real))?"  # Optional variant prefix: "real" before location
        r"(?P<location>dist|prox|)"  # Location: "dist" (distal) or "prox" (proximal)
        r"(?P<location_number>\d+)"  # Location number: 1, 2, etc.
        r"(?:(?P<variant_suffix>real))?"  # Optional variant suffix: "real" after location
        r"(?:_(?P<experiment_type>trio|bsl))?"  # Optional experiment type: "trio" or "bsl"
        r"-(?P<trial_number>\d+)"  # Trial number: 001, 002, etc.
    )
    match = re.match(pattern, folder_name)

    if not match:
        raise ValueError(f"Could not parse recording folder name: {folder_name}")

    date_str = match.group("date")
    cell_number = match.group("cell_number")
    location = match.group("location")
    location_number = match.group("location_number")
    variant_prefix = match.group("variant_prefix") or ""  # "real" before location
    variant_suffix = match.group("variant_suffix") or ""  # "real" after location
    variant = variant_prefix or variant_suffix  # Use whichever is present
    experiment_type = match.group("experiment_type") or ""  # "trio", "bsl", or empty
    trial_number = match.group("trial_number")

    # Handle different date formats based on first digit
    if date_str[0] == "2":
        # YYYYMMDD format (e.g., 20160523)
        year = int(date_str[:4])
        month = int(date_str[4:6])
        day = int(date_str[6:8])
    else:
        # MMDDYYYY format (e.g., 05232016)
        month = int(date_str[:2])
        day = int(date_str[2:4])
        year = int(date_str[4:8])

    # Validate date components
    if month < 1 or month > 12:
        raise ValueError(f"Invalid month '{month}' in date '{date_str}' from folder '{folder_name}'")
    if day < 1 or day > 31:
        raise ValueError(f"Invalid day '{day}' in date '{date_str}' from folder '{folder_name}'")

    # Determine location description
    location_full = "Proximal" if location == "prox" else "Distal"
    variant_suffix = f" {variant.capitalize()}" if variant else ""
    location_description = f"{location_full} dendrite {location_number}{variant_suffix}"
    base_line_experiment_type = "Baseline" if experiment_type == "bsl" else ""

    # Create datetime object for the date
    session_date = datetime(year, month, day)

    return {
        "cell_number": cell_number,
        "location": location,
        "location_number": location_number,
        "location_full": location_full,
        "location_description": location_description,
        "experiment_type": experiment_type,
        "trial_number": trial_number,
        "variant": variant,
        "base_line_experiment_type": base_line_experiment_type,
        "date": session_date,
    }


def convert_session_to_nwbfile(session_folder_path: Path, condition: str, verbose: bool = False) -> NWBFile:
    """
    Convert all dendritic excitability recordings from a session folder to a single NWB format file.

    Parameters
    ----------
    session_folder_path : Path
        Path to the session folder (corresponds to a single subject)
    condition : str
        Experimental condition (e.g., "LID off-state", "LID on-state", "LID on-state with sul")
    verbose : bool, default=False
        Enable verbose output showing detailed processing information

    Returns
    -------
    NWBFile
        NWB file with all the converted data from the session

    Notes
    -----
    Each session folder corresponds to a single subject. All recordings within the session
    are combined into a single NWBFile.
    """

    # Get all recording folders within the session folder
    recording_folders = [f for f in session_folder_path.iterdir() if f.is_dir()]
    recording_folders.sort()

    if not recording_folders:
        raise ValueError(f"No recording folders found in session folder: {session_folder_path}")

    # Calculate recording IDs, session start times, and create interface mappings
    ophys_session_start_times = []  # (ophys_time, recording_folder, recording_id)
    intracellular_session_start_times = []  # (intracellular_time, recording_folder, recording_id)
    recording_id_to_location_id = {}
    recording_id_to_folder = {}
    t_starts = {}  # t_starts[recording_id][interface] = t_start_offset

    for recording_folder in recording_folders:
        # Find main experiment XML file (ophys)
        main_xml_file = recording_folder / f"{recording_folder.name}.xml"
        if not main_xml_file.exists():
            raise FileNotFoundError(f"Expected main XML file does not exist: {main_xml_file}")

        # Find electrophysiology XML file (intracellular)
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

        # Get unique identifiers for recording to name objects
        recording_info = parse_session_info_from_folder_name(recording_folder)
        repetition_id = f"{recording_info['base_line_experiment_type']}Trial{recording_info['trial_number']}{recording_info['variant']}"
        location_id = f"{recording_info['location_full']}Dendrite{recording_info['location_number']}"
        recording_id = f"{location_id}{repetition_id}"

        # Store mappings
        recording_id_to_location_id[recording_id] = location_id
        recording_id_to_folder[recording_id] = recording_folder

        # Store session start times separately by interface type
        ophys_session_start_times.append((ophys_session_start_time, recording_folder, recording_id))
        intracellular_session_start_times.append((intracellular_session_start_time, recording_folder, recording_id))

    if not ophys_session_start_times or not intracellular_session_start_times:
        raise ValueError(f"No valid recordings found in session folder: {session_folder_path}")

    # Find the earliest session start time from each interface type
    earliest_ophys_time = min(ophys_session_start_times, key=lambda x: x[0])[0]
    earliest_intracellular_time = min(intracellular_session_start_times, key=lambda x: x[0])[0]

    # Overall session start time is the earliest across all interfaces
    session_start_time = min(earliest_ophys_time, earliest_intracellular_time)

    # Calculate t_start offsets for temporal alignment with interface-specific timing
    for ophys_time, folder, recording_id in ophys_session_start_times:
        intracellular_time = next(time for time, _, rid in intracellular_session_start_times if rid == recording_id)

        # Calculate offsets relative to overall session start time
        ophys_t_start = (ophys_time - session_start_time).total_seconds()
        intracellular_t_start = (intracellular_time - session_start_time).total_seconds()

        # Initialize t_starts for this recording_id with interface-specific timing
        # Use descriptive names for clarity about which channel/interface this refers to
        t_starts[recording_id] = {
            "intracellular": intracellular_t_start,
            "line_scan_structural_channel": ophys_t_start,  # Ch1/Alexa568 line scan uses ophys timing
            "line_scan_calcium_channel": ophys_t_start,  # Ch2/Fluo4 line scan uses ophys timing
        }

    # Load general and session-specific metadata from YAML files
    general_metadata_path = Path(__file__).parent.parent.parent / "general_metadata.yaml"
    general_metadata = load_dict_from_file(general_metadata_path)

    session_metadata_path = Path(__file__).parent.parent.parent / "session_specific_metadata.yaml"
    session_metadata_template = load_dict_from_file(session_metadata_path)
    script_template = session_metadata_template["figure_3_dendritic_excitability"]

    # Create session-specific metadata using session start time from XML

    # Get first recording info for session description, all of them are equal
    first_recording_folder = next(iter(recording_id_to_folder.values()))
    first_recording_info = parse_session_info_from_folder_name(first_recording_folder)

    # Create session ID following pattern from somatic excitability scripts
    cell_type = "iSPN"  # Indirect pathway SPN for dendritic excitability experiments
    timestamp = session_start_time.strftime("%Y%m%d%H%M%S")
    condition_camel_case = format_condition[condition]["CamelCase"]  # Use centralized mapping
    base_session_id = f"Figure3++DendriticExcitability++{condition_camel_case}++{timestamp}"
    script_specific_id = f"{cell_type}++Cell++{first_recording_info['cell_number']}"
    session_id = f"{base_session_id}++{script_specific_id}"

    # Handle pharmacology conditions dynamically
    pharmacology_text = general_metadata["NWBFile"]["pharmacology"]
    if "sul" in condition.lower() and "pharmacology_conditions" in script_template:
        if "sulpiride" in script_template["pharmacology_conditions"]:
            pharmacology_text += " " + script_template["pharmacology_conditions"]["sulpiride"]

    # Create session-specific metadata from template with runtime substitutions
    condition_underscore = format_condition[condition]["underscore"]
    condition_human_readable = format_condition[condition]["human_readable"]
    session_specific_metadata = {
        "NWBFile": {
            "session_description": script_template["NWBFile"]["session_description"].format(
                condition=condition_human_readable,
                cell_number=first_recording_info["cell_number"],
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

    # Deep merge with general metadata
    metadata = dict_deep_update(general_metadata, session_specific_metadata)

    # Create NWB file using neuroconv helper function
    nwbfile = make_nwbfile_from_metadata(metadata)

    # Add custom columns to intracellular recording table for dendritic experiment annotations
    intracellular_recording_table = nwbfile.get_intracellular_recordings()
    intracellular_recording_table.add_column(
        name="stimulus_protocol", description="Current injection protocol used for dendritic excitability testing"
    )
    intracellular_recording_table.add_column(
        name="dendrite_distance_um", description="Approximate distance from soma in micrometers"
    )
    intracellular_recording_table.add_column(
        name="dendrite_type", description="Type of dendritic location: Distal or Proximal"
    )
    intracellular_recording_table.add_column(
        name="dendrite_number", description="Number identifier for the specific dendritic location (1, 2, etc.)"
    )
    intracellular_recording_table.add_column(
        name="trial_number", description="Trial number for this dendritic location (1, 2, 3)"
    )
    intracellular_recording_table.add_column(
        name="recording_id", description="Full recording identifier containing location and trial information"
    )

    # Data structures for tracking icephys table indices
    recording_indices = []  # Store all intracellular recording indices
    recording_to_metadata = {}  # Map recording index to metadata for table building
    location_to_recording_indices = {}  # Group recordings by location for repetitions table
    sequential_recording_indices = []  # Store sequential recording indices

    # Process each recording using the calculated recording IDs
    for recording_id, recording_folder in recording_id_to_folder.items():
        location_id = recording_id_to_location_id[recording_id]
        recording_info = parse_session_info_from_folder_name(recording_folder)
        main_xml_file = recording_folder / f"{recording_folder.name}.xml"

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

        # Find electrophysiology XML file (exact name from Figure 3 notes)
        electrophysiology_xml_file = recording_folder / f"{recording_folder.name}_Cycle00001_VoltageRecording_001.xml"

        if not electrophysiology_xml_file.exists():
            raise FileNotFoundError(f"Expected electrophysiology XML file does not exist: {electrophysiology_xml_file}")

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

        # Update electrode description for iSPN dendritic recording
        # One electrode per location (not for trial)
        electrode_name = f"IntracellularElectrode{location_id}"
        intracellular_metadata["Icephys"]["IntracellularElectrodes"][icephys_metadata_key].update(
            {
                "name": electrode_name,
                "description": (
                    f"Recording from {cell_type} {recording_info['location_description']} - {condition_human_readable} - "
                    f"Cell {recording_info['cell_number']} - Trial {recording_info['trial_number']}"
                ),
                "cell_id": f"CellRecordedAt{timestamp}",
                "location": recording_info["location_description"],
                "slice": general_metadata["NWBFile"]["slices"],
            }
        )

        # Update current clamp series metadata
        series_name = f"CurrentClampSeries{recording_id}"
        intracellular_metadata["Icephys"]["CurrentClampSeries"][icephys_metadata_key].update(
            {
                "name": series_name,
                "description": (
                    f"Current clamp recording from {cell_type} {recording_info['location_description']} - "
                    f"{condition_human_readable} - Cell {recording_info['cell_number']} - Trial {recording_info['trial_number']} - "
                    f"Three 2 nA current injections, 2 ms each, at 50 Hz"
                ),
            }
        )

        # Add intracellular data to NWB file
        intracellular_interface.add_to_nwbfile(nwbfile=nwbfile, metadata=intracellular_metadata)

        # Add intracellular recording to icephys table with custom annotations
        current_clamp_series = nwbfile.acquisition[series_name]

        # Calculate dendrite distance based on location type (approximate values from literature)
        dendrite_distance_um = 90 if recording_info["location"] == "dist" else 40  # distal ~90μm, proximal ~40μm

        # Add intracellular recording entry with essential metadata annotations
        recording_index = nwbfile.add_intracellular_recording(
            electrode=current_clamp_series.electrode,
            response=current_clamp_series,
            stimulus_protocol="3x2nA_2ms_50Hz_dendritic_excitability",
            dendrite_distance_um=dendrite_distance_um,
            dendrite_type=recording_info["location_full"],  # "Distal" or "Proximal"
            dendrite_number=int(recording_info["location_number"]),
            trial_number=int(recording_info["trial_number"]),
            recording_id=recording_id,
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
        ] = "Bruker two-photon microscope for line scan imaging"
        structural_metadata["Ophys"]["ImagingPlanes"][structural_ophys_key]["name"] = f"ImagingPlane{location_id}"
        structural_metadata["Ophys"]["ImagingPlanes"][structural_ophys_key][
            "description"
        ] = f"Line scan imaging plane for {location_id} using Alexa Fluor 568 structural dye"
        structural_metadata["Ophys"]["ImagingPlanes"][structural_ophys_key]["indicator"] = "Alexa Fluor 568"
        structural_metadata["Ophys"]["PlaneSegmentation"][structural_ophys_key][
            "name"
        ] = f"PlaneSegmentation{recording_id}"
        structural_metadata["Ophys"]["PlaneSegmentation"][structural_ophys_key][
            "description"
        ] = f"Line scan ROI segmentation for {location_id} structural imaging"
        structural_metadata["Ophys"]["RoiResponseSeries"][structural_ophys_key][
            "name"
        ] = f"RoiResponseSeriesAlexa568{recording_id}"
        structural_metadata["Ophys"]["RoiResponseSeries"][structural_ophys_key][
            "description"
        ] = f"Structural reference fluorescence from Alexa Fluor 568 hydrazide (50 μM) - {location_id}"
        structural_metadata["Acquisition"]["SourceImages"][structural_ophys_key][
            "name"
        ] = f"ImageAlexa568{recording_id}"
        structural_metadata["Acquisition"]["SourceImages"][structural_ophys_key][
            "description"
        ] = f"Source image for Alexa Fluor 568 structural reference - {location_id}"
        structural_metadata["TimeSeries"][structural_ophys_key]["name"] = f"TimeSeriesLineScanRawAlexa568{recording_id}"
        structural_metadata["TimeSeries"][structural_ophys_key][
            "description"
        ] = f"Line scan raw data for Alexa Fluor 568 structural reference - {location_id}"

        # Add structural data to NWB file
        structural_interface.add_to_nwbfile(nwbfile=nwbfile, metadata=structural_metadata)

        # Process calcium channel (Ch2/Fluo4)
        calcium_metadata = calcium_interface.get_metadata()
        # Apply fluorophore-specific metadata based on experimental knowledge
        calcium_metadata["Devices"][calcium_ophys_key]["name"] = "BrukerUltima"
        calcium_metadata["Devices"][calcium_ophys_key][
            "description"
        ] = "Bruker two-photon microscope for line scan imaging"
        ophys_metadata = calcium_metadata["Ophys"]
        ophys_metadata["ImagingPlanes"][calcium_ophys_key]["name"] = f"ImagingPlane{location_id}"
        ophys_metadata["ImagingPlanes"][calcium_ophys_key][
            "description"
        ] = f"Line scan imaging plane for {location_id} using Fluo-4 calcium indicator"
        ophys_metadata["ImagingPlanes"][calcium_ophys_key]["indicator"] = "Fluo-4"
        ophys_metadata["PlaneSegmentation"][calcium_ophys_key]["name"] = f"PlaneSegmentation{recording_id}"
        ophys_metadata["PlaneSegmentation"][calcium_ophys_key][
            "description"
        ] = f"Line scan ROI segmentation for {location_id} calcium imaging"
        ophys_metadata["RoiResponseSeries"][calcium_ophys_key]["name"] = f"RoiResponseSeriesFluo4{recording_id}"
        ophys_metadata["RoiResponseSeries"][calcium_ophys_key][
            "description"
        ] = f"Calcium fluorescence from Fluo-4 (100 μM) - {location_id}"
        calcium_metadata["Acquisition"]["SourceImages"][calcium_ophys_key]["name"] = f"ImageFluo4{recording_id}"
        calcium_metadata["Acquisition"]["SourceImages"][calcium_ophys_key][
            "description"
        ] = f"Source image for Fluo-4 calcium indicator - {location_id}"
        calcium_metadata["TimeSeries"][calcium_ophys_key]["name"] = f"TimeSeriesLineScanRawFluo4{recording_id}"
        calcium_metadata["TimeSeries"][calcium_ophys_key][
            "description"
        ] = f"Line scan raw data for Fluo-4 calcium indicator - {location_id}"

        # Add calcium data to NWB file
        calcium_interface.add_to_nwbfile(nwbfile=nwbfile, metadata=calcium_metadata)

        # Add reference images for this recording
        references_folder = recording_folder / "References"
        ref_container_name = f"ImagesBackground{recording_id}"
        reference_interface = BrukerReferenceImagesInterface(
            references_folder_path=references_folder, container_name=ref_container_name
        )
        reference_interface.add_to_nwbfile(nwbfile=nwbfile)

    # Build icephys table hierarchical structure using shared utility function
    build_dendritic_icephys_table_structure(
        nwbfile=nwbfile,
        recording_indices=recording_indices,
        recording_to_metadata=recording_to_metadata,
        session_info={"cell_number": first_recording_info["cell_number"]},
        condition=condition_underscore,
        stimulus_type="dendritic_excitability_current_injection",
        verbose=verbose,
    )

    # Add trials table using interface
    trials_interface = DendriticTrialsInterface(
        recording_indices=recording_indices, recording_to_metadata=recording_to_metadata, t_starts=t_starts
    )
    trials_interface.add_to_nwbfile(nwbfile, verbose=verbose)

    return nwbfile


if __name__ == "__main__":
    import argparse
    import logging

    from tqdm import tqdm

    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Convert Figure 3 dendritic excitability data to NWB format")
    parser.add_argument(
        "--stub-test",
        type=str_to_bool,
        default=True,
        help="Process only first 2 files per condition for testing (default: True). Use --stub-test=False for full processing.",
    )

    args = parser.parse_args()
    stub_test = args.stub_test

    # Set to True to enable verbose output
    verbose = False

    # Suppress tifffile warnings
    logging.getLogger("tifffile").setLevel(logging.ERROR)

    # Suppress specific warnings
    warnings.filterwarnings("ignore", message="invalid value encountered in divide")
    warnings.filterwarnings("ignore", message=".*no datetime before year 1.*")

    # Set the base path to your data
    base_path = Path("./link_to_raw_data/Figure 3/Dendritic excitability")

    # Create nwb_files directory at root level
    root_dir = Path(__file__).parent.parent.parent.parent.parent.parent  # Go up to repo root
    nwb_files_dir = root_dir / "nwb_files" / "dendritic_excitability" / "figure_3"
    nwb_files_dir.mkdir(parents=True, exist_ok=True)

    # Figure 3 dendritic conditions
    conditions = ["LID off-state", "LID on-state", "LID on-state with sul"]

    for condition in conditions:
        condition_path = base_path / condition

        if not condition_path.exists():
            raise FileNotFoundError(f"Expected condition path does not exist: {condition_path}")

        # Normalize folder name to paper condition name
        paper_condition = FOLDER_TO_PAPER_CONDITION.get(condition, condition)

        # Get all session folders (e.g., 0523a, 0523b, etc.)
        # Note: Each session folder corresponds to a single subject
        session_folders = [f for f in condition_path.iterdir() if f.is_dir()]
        session_folders.sort()

        if stub_test:
            session_folders = session_folders[:2]  # Take only the first two sessions for stub test

        # Use tqdm for progress bar when verbose is disabled
        session_iterator = tqdm(
            session_folders,
            desc=f"Converting Figure3 DendriticExcitability {format_condition[paper_condition]['human_readable']}",
            unit=" session",
        )

        for session_folder in session_iterator:

            # Convert all recordings from this session to NWB format
            # All recordings from the session are combined into a single NWBFile
            nwbfile = convert_session_to_nwbfile(
                session_folder_path=session_folder,
                condition=paper_condition,
                verbose=verbose,
            )

            # Create output filename
            condition_safe = condition.replace(" ", "_").replace("-", "_").replace("(", "").replace(")", "")
            nwbfile_path = nwb_files_dir / f"{nwbfile.session_id}.nwb"

            # Write NWB file
            configure_and_write_nwbfile(nwbfile, nwbfile_path=nwbfile_path)
