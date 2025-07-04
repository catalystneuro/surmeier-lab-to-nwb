from datetime import datetime
from pathlib import Path
from typing import Any, Dict
from zoneinfo import ZoneInfo

from neuroconv.tools import configure_and_write_nwbfile
from neuroconv.utils import dict_deep_update, load_dict_from_file
from pynwb import NWBFile
from pynwb.file import Subject

from surmeier_lab_to_nwb.zhai2025.intracellular_interfaces import (
    PrairieViewCurrentClampInterface,
)


def parse_session_info(session_folder: Path) -> Dict[str, Any]:
    """
    Parse session information from folder structure.

    Expected folder name format: MMDDYYYY_N (e.g., 05232016_1)

    Parameters
    ----------
    session_folder : Path
        Path to the session folder

    Returns
    -------
    Dict[str, Any]
        Dictionary containing session information
    """
    folder_name = session_folder.name

    # Parse date and cell number from folder name
    # Format: MMDDYYYY_N
    if "_" in folder_name:
        date_str, cell_number = folder_name.split("_")
    else:
        raise ValueError(f"Could not parse session folder name: {folder_name}")

    # Parse date (MMDDYYYY format)
    if len(date_str) == 8:
        month = int(date_str[:2])
        day = int(date_str[2:4])
        year = int(date_str[4:8])
    else:
        raise ValueError(f"Could not parse date from folder name: {folder_name}")

    # Illinois is in Central Time Zone
    central_tz = ZoneInfo("America/Chicago")
    session_start_time = datetime(year, month, day, 0, 0, 0, tzinfo=central_tz)

    return {
        "session_start_time": session_start_time,
        "cell_number": cell_number,
        "date_str": f"{year}-{month:02d}-{day:02d}",
        "session_id": f"{year}{month:02d}{day:02d}_{cell_number}",
    }


def convert_session_to_nwbfile(session_folder_path: Path, condition: str) -> NWBFile:
    """
    Convert a single session of Figure 3 somatic excitability data to NWB format.

    Parameters
    ----------
    session_folder_path : Path
        Path to the session folder containing cell recordings
    condition : str
        Experimental condition (e.g., "LID off-state", "LID on-state", "LID on-state with sul (iSPN)")

    Returns
    -------
    NWBFile
        NWB file with the converted data
    """

    # Parse session information
    session_info = parse_session_info(session_folder_path)
    print(f"Session date: {session_info['date_str']}, Cell: {session_info['cell_number']}")

    # Load metadata from YAML file
    metadata_file_path = Path(__file__).parent / "metadata.yaml"
    paper_metadata = load_dict_from_file(metadata_file_path)

    # Create session-specific metadata
    session_specific_metadata = {
        "NWBFile": {
            "session_description": (
                f"Somatic excitability assessment in indirect pathway spiny projection neurons (iSPNs) "
                f"for condition '{condition}'. Whole-cell patch clamp recording in current clamp mode "
                f"with current injection steps from -120 pA to +300 pA (500 ms duration each). "
                f"Cell {session_info['cell_number']} recorded on {session_info['date_str']}."
            ),
            "identifier": f"zhai2025_fig3_somatic_{session_info['session_id']}_{condition.replace(' ', '_')}",
            "session_start_time": session_info["session_start_time"],
            "experiment_description": (
                f"Somatic excitability changes in iSPNs during condition '{condition}'. "
                f"This experiment is part of Figure 3 from Zhai et al. 2025, investigating how LID affects "
                f"iSPN excitability and the role of D2 receptor signaling."
            ),
            "session_id": session_info["session_id"],
            "keywords": [
                "somatic excitability",
                "F-I relationship" "rheobase",
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
        subject_id=f"iSPN_mouse_{session_info['session_id']}",
        species="Mus musculus",
        strain="C57Bl/6",
        description=(
            f"Experimental mouse with unilateral 6-OHDA lesion in the medial forebrain bundle. "
            f"iSPNs identified by lack of Drd1-Tdtomato expression (negative selection). "
            f"Cell {session_info['cell_number']} recorded on {session_info['date_str']}."
        ),
        genotype="Drd1-Tdtomato bacterial artificial chromosome (BAC) transgenic",
        sex="M",
        age="P49-P84",  # ISO format for 7-12 weeks (postnatal days)
    )
    nwbfile.subject = subject

    # Add columns to intracellular recording tables
    intracellular_recording_table = nwbfile.get_intracellular_recordings()
    intracellular_recording_table.add_column(
        name="stimulus_current_pA",
        description="The current stimulus applied during the recording in picoamps",
    )

    # Input current values for each recording
    protocol_step_to_current = {
        "001": -120,  # pA
        "002": -100,  # pA
        "003": -80,  # pA
        "004": -60,  # pA
        "005": -40,  # pA
        "006": -20,  # pA
        "007": 20,  # pA
        "008": 40,  # pA
        "009": 60,  # pA
        "010": 80,  # pA
        "011": 100,  # pA
        "012": 120,  # pA
        "013": 140,  # pA
        "014": 160,  # pA
        "015": 180,  # pA
        "016": 200,  # pA
        "017": 220,  # pA
        "018": 240,  # pA
        "019": 260,  # pA
        "020": 280,  # pA
        "021": 300,  # pA
        "022": 320,  # pA  # some cells may go further than 300 pA
        "023": 340,  # pA
        "024": 360,  # pA
        "025": 380,  # pA
        "026": 400,  # pA
    }

    # Find all recording folders (current steps) for this cell
    recording_folders = [f for f in session_folder_path.iterdir() if f.is_dir()]
    recording_folders.sort()

    print(f"Found {len(recording_folders)} current step recordings in session")

    # Store recording indices for building intracellular tables
    recording_indices = []

    # Process each recording (current step)
    for recording_folder in recording_folders:
        print(f"  Processing: {recording_folder.name}")

        # Extract protocol step from recording folder name (e.g., "cell1-001" -> "001")
        if "-" in recording_folder.name:
            protocol_step = recording_folder.name.split("-")[1]
        else:
            raise ValueError(f"Could not extract protocol step from recording folder name: {recording_folder.name}")

        # Get current value for this protocol step
        if protocol_step not in protocol_step_to_current:
            raise ValueError(
                f"Protocol step {protocol_step} not found in protocol_step_to_current mapping for {recording_folder.name}"
            )

        current_pA = protocol_step_to_current[protocol_step]

        # Format current with consistent 3-digit format and sign
        current_formatted = f"{current_pA:+04d}pA"  # e.g., "+020pA", "-120pA"

        # Find XML file for this recording (exact name from Figure 3 notes)
        xml_file = recording_folder / f"{recording_folder.name}_Cycle00001_VoltageRecording_001.xml"

        if not xml_file.exists():
            raise FileNotFoundError(f"Expected XML file does not exist: {xml_file}")

        # Create interface for this recording
        interface = PrairieViewCurrentClampInterface(
            file_path=xml_file, icephys_metadata_key=f"PrairieView_{recording_folder.name}"
        )

        # Get metadata and update electrode description
        metadata = interface.get_metadata()

        # Update electrode metadata for iSPN
        electrode_name = "IntracellularElectrode"
        metadata["Icephys"]["IntracellularElectrodes"][electrode_name].update(
            {
                "description": f"Recording from iSPN in the dorsolateral striatum - {condition} - Cell {session_info['cell_number']} - {current_formatted} step",
                "cell_id": session_info["cell_number"],
            }
        )

        # Update CurrentClampSeries metadata with current annotation
        series_key = f"PrairieView_{recording_folder.name}"
        series_name = f"CurrentClampSeries{current_formatted}"

        metadata["Icephys"]["CurrentClampSeries"][series_key].update(
            {
                "name": series_name,
                "description": f"Current clamp recording from iSPN - {condition} - Cell {session_info['cell_number']} - {current_formatted} current injection",
            }
        )

        # Add to NWB file and get recording index
        interface.add_to_nwbfile(nwbfile=nwbfile, metadata=metadata)

        # Find the current clamp series that was just added
        current_clamp_series = nwbfile.acquisition[series_name]

        # Add intracellular recording entry
        recording_index = nwbfile.add_intracellular_recording(
            electrode=current_clamp_series.electrode,
            response=current_clamp_series,
            stimulus_current_pA=current_pA,
        )

        recording_indices.append(recording_index)
        print(f"Successfully added recording: {recording_folder.name} ({current_formatted})")

    # Build intracellular recording structure
    # Add simultaneous recordings (one per current step)
    simultaneous_recording_indices = []
    for recording_index in recording_indices:
        simultaneous_index = nwbfile.add_icephys_simultaneous_recording(recordings=[recording_index])
        simultaneous_recording_indices.append(simultaneous_index)

    # Add sequential recording for this cell (grouping all current steps)
    nwbfile.add_icephys_sequential_recording(
        simultaneous_recordings=simultaneous_recording_indices,
        stimulus_type="F-I protocol",
    )

    print(f"Added intracellular recording structure for cell {session_info['cell_number']}")

    return nwbfile


if __name__ == "__main__":
    # Set the base path to your data
    base_path = Path("./link_to_raw_data/Figure 3/Somatic excitability")

    # Create nwb_files directory at root level
    root_dir = Path(__file__).parent.parent.parent.parent  # Go up to repo root
    nwb_files_dir = root_dir / "nwb_files" / "figure_3_somatic_excitability"
    nwb_files_dir.mkdir(parents=True, exist_ok=True)

    # Figure 3 conditions
    conditions = ["LID off-state", "LID on-state", "LID on-state with sul (iSPN)"]

    for condition in conditions:
        condition_path = base_path / condition

        if not condition_path.exists():
            raise FileNotFoundError(f"Expected condition path does not exist: {condition_path}")

        print(f"Processing somatic excitability data for: {condition}")

        # Get all session folders (each session = one cell)
        session_folders = [f for f in condition_path.iterdir() if f.is_dir()]
        session_folders.sort()

        print(f"Found {len(session_folders)} session folders")

        for session_folder in session_folders:
            print(f"\nProcessing session: {session_folder.name}")

            # Convert session data to NWB format
            nwbfile = convert_session_to_nwbfile(
                session_folder_path=session_folder,
                condition=condition,
            )

            # Create output filename
            condition_safe = condition.replace(" ", "_").replace("(", "").replace(")", "")
            nwbfile_path = nwb_files_dir / f"figure3_somatic_excitability_{condition_safe}_{session_folder.name}.nwb"

            # Write NWB file
            configure_and_write_nwbfile(nwbfile, nwbfile_path=nwbfile_path)
            print(f"Successfully saved: {nwbfile_path}")
