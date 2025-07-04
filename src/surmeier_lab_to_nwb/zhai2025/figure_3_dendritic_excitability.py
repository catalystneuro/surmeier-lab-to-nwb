import re
import warnings
from datetime import datetime
from pathlib import Path
from typing import Any, Dict

from neuroconv.tools import configure_and_write_nwbfile
from neuroconv.utils import dict_deep_update, load_dict_from_file
from pynwb import NWBFile
from pynwb.file import Subject

from surmeier_lab_to_nwb.zhai2025.intracellular_interfaces import (
    PrairieViewCurrentClampInterface,
)
from surmeier_lab_to_nwb.zhai2025.ophys_interfaces import (
    PrairieViewLineScanInterface,
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


def convert_session_to_nwbfile(session_folder_path: Path, condition: str) -> NWBFile:
    """
    Convert all dendritic excitability recordings from a session folder to a single NWB format file.

    Parameters
    ----------
    session_folder_path : Path
        Path to the session folder (corresponds to a single subject)
    condition : str
        Experimental condition (e.g., "LID off-state", "LID on-state", "LID on-state with sul")

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

    print(f"Processing session folder: {session_folder_path.name} (corresponds to one subject)")
    print(f"  Found {len(recording_folders)} recordings")

    # Get session start times for all recordings to find the earliest one
    session_start_times = []
    recording_infos = []
    for recording_folder in recording_folders:
        # Find main experiment XML file
        main_xml_file = recording_folder / f"{recording_folder.name}.xml"
        if not main_xml_file.exists():
            raise FileNotFoundError(f"Expected main XML file does not exist: {main_xml_file}")

        # Get session start time using static method
        session_start_time = PrairieViewLineScanInterface.get_session_start_time_from_file(main_xml_file)
        recording_info = parse_session_info_from_folder_name(recording_folder)

        session_start_times.append((session_start_time, recording_folder))
        recording_infos.append((recording_folder, recording_info))

    if not session_start_times:
        raise ValueError(f"No valid recordings found in session folder: {session_folder_path}")

    # Find the earliest session start time
    session_start_time, earliest_folder = min(session_start_times, key=lambda x: x[0])
    print(f"  Session start time: {session_start_time} (from {earliest_folder.name})")

    # Create a dictionary of recording info for easy lookup
    recording_info_dict = dict(recording_infos)

    # Load metadata from YAML file
    metadata_file_path = Path(__file__).parent / "metadata.yaml"
    paper_metadata = load_dict_from_file(metadata_file_path)

    # Create session-specific metadata using session start time from XML
    session_date_str = session_start_time.strftime("%Y-%m-%d")

    # Get first recording info for session description, all of them are equal
    first_recording_info = recording_infos[0][1]

    session_specific_metadata = {
        "NWBFile": {
            "session_description": (
                f"Dendritic excitability assessment in indirect pathway spiny projection neurons (iSPNs) "
                f"for condition '{condition}'. Combined patch clamp electrophysiology and two-photon "
                f"laser scanning microscopy. Brief current steps (three 2 nA injections, 2 ms each, at 50 Hz) "
                f"with Ca2+ imaging. Cell {first_recording_info['cell_number']} in date {first_recording_info['date'].strftime('%Y-%m-%d')}."
            ),
            "identifier": f"zhai2025_fig3_dendritic_{session_folder_path.name}_{condition.replace(' ', '_')}",
            "session_start_time": session_start_time,
            "experiment_description": (
                f"Dendritic excitability changes in iSPNs during condition '{condition}'. "
                f"This experiment is part of Figure 3 from Zhai et al. 2025, investigating how LID affects "
                f"iSPN excitability and the role of D2 receptor signaling. Multiple recordings from one subject."
            ),
            "session_id": f"{session_folder_path.name}_{condition.replace(' ', '_')}",
            "keywords": [
                "intracellular electrophysiology",
                "patch clamp",
                "two-photon microscopy",
                "calcium imaging",
                "dendritic excitability",
                "line scan",
                "current injection",
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

    # Create subject metadata for iSPN experiments
    subject = Subject(
        subject_id=f"iSPN_mouse_{session_folder_path.name}",
        species="Mus musculus",
        strain="C57Bl/6",
        description=(
            f"Experimental mouse with unilateral 6-OHDA lesion in the medial forebrain bundle. "
            f"iSPNs identified by lack of Drd1-Tdtomato expression (negative selection). "
            f"Session {session_folder_path.name} recorded on {session_date_str}."
        ),
        genotype="Drd1-Tdtomato bacterial artificial chromosome (BAC) transgenic",
        sex="M",
        age="P49-P84",  # ISO format for 7-12 weeks (postnatal days)
    )
    nwbfile.subject = subject

    # Process each recording
    for recording_folder in recording_folders:
        recording_info = recording_info_dict[recording_folder]
        main_xml_file = recording_folder / f"{recording_folder.name}.xml"

        # Create line scan interface for this recording
        line_scan_interface = PrairieViewLineScanInterface(xml_metadata_file_path=main_xml_file)

        repetition_id = f"{recording_info['base_line_experiment_type']}Trial{recording_info['trial_number']}{recording_info['variant']}"
        location_id = f"{recording_info['location_full']}Dendrite{recording_info['location_number']}"
        recording_id = f"{location_id}{repetition_id}"

        print(f"  Processing recording for folder: {recording_folder.name}")
        print(f"    Recording ID: {recording_id}")

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

        # Get and update intracellular metadata
        intracellular_metadata = intracellular_interface.get_metadata()

        # Update electrode description for iSPN dendritic recording
        # One electrode per location (not for trial)
        electrode_name = f"IntracellularElectrode{location_id}"
        intracellular_metadata["Icephys"]["IntracellularElectrodes"][icephys_metadata_key].update(
            {
                "name": electrode_name,
                "description": (
                    f"Recording from iSPN {recording_info['location_description']} - {condition} - "
                    f"Cell {recording_info['cell_number']} - Trial {recording_info['trial_number']}"
                ),
                "cell_id": recording_info["cell_number"],
                "location": recording_info["location_description"],
            }
        )

        # Update current clamp series metadata
        series_name = f"CurrentClampSeries{recording_id}"
        intracellular_metadata["Icephys"]["CurrentClampSeries"][icephys_metadata_key].update(
            {
                "name": series_name,
                "description": (
                    f"Current clamp recording from iSPN {recording_info['location_description']} - "
                    f"{condition} - Cell {recording_info['cell_number']} - Trial {recording_info['trial_number']} - "
                    f"Three 2 nA current injections, 2 ms each, at 50 Hz"
                ),
            }
        )

        # Add intracellular data to NWB file
        intracellular_interface.add_to_nwbfile(nwbfile=nwbfile, metadata=intracellular_metadata)

        # Get line scan metadata and update with recording-specific information
        line_scan_metadata = line_scan_interface.get_metadata()

        line_scan_metadata["recording_id"] = recording_id
        line_scan_metadata["location_id"] = location_id

        # Add line scan data to NWB file with updated metadata
        line_scan_interface.add_to_nwbfile(nwbfile=nwbfile, metadata=line_scan_metadata)

        print(f"    Added line scan imaging data")
        print(f"    Successfully processed recording: {recording_folder.name}")

    print(f"Successfully processed all recordings from session: {session_folder_path.name}")

    return nwbfile


if __name__ == "__main__":
    import logging

    # Suppress tifffile warnings
    logging.getLogger("tifffile").setLevel(logging.ERROR)

    # Suppress specific warnings
    warnings.filterwarnings("ignore", message="invalid value encountered in divide")
    warnings.filterwarnings("ignore", message=".*no datetime before year 1.*")

    # Set the base path to your data
    base_path = Path("./link_to_raw_data/Figure 3/Dendritic excitability")

    # Create nwb_files directory at root level
    root_dir = Path(__file__).parent.parent.parent.parent  # Go up to repo root
    nwb_files_dir = root_dir / "nwb_files" / "figure_3_dendritic_excitability"
    nwb_files_dir.mkdir(parents=True, exist_ok=True)

    # Figure 3 dendritic conditions
    conditions = ["LID off-state", "LID on-state", "LID on-state with sul"]

    for condition in conditions:
        condition_path = base_path / condition

        if not condition_path.exists():
            raise FileNotFoundError(f"Expected condition path does not exist: {condition_path}")

        print(f"Processing dendritic excitability data for: {condition}")

        # Get all session folders (e.g., 0523a, 0523b, etc.)
        # Note: Each session folder corresponds to a single subject
        session_folders = [f for f in condition_path.iterdir() if f.is_dir()]
        session_folders.sort()

        print(f"Found {len(session_folders)} session folders")

        for session_folder in session_folders:
            print(f"\nProcessing session folder: {session_folder.name}")

            # Convert all recordings from this session to NWB format
            # All recordings from the session are combined into a single NWBFile
            nwbfile = convert_session_to_nwbfile(
                session_folder_path=session_folder,
                condition=condition,
            )

            # Create output filename
            condition_safe = condition.replace(" ", "_").replace("(", "").replace(")", "")
            nwbfile_path = nwb_files_dir / f"figure3_dendritic_excitability_{condition_safe}_{session_folder.name}.nwb"

            # Write NWB file
            configure_and_write_nwbfile(nwbfile, nwbfile_path=nwbfile_path)
            print(f"Successfully saved: {nwbfile_path.name}")
