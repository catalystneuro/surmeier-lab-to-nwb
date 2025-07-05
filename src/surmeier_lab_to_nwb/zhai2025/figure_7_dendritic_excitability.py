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
    Parse essential recording information from Figure 7 dendritic excitability recording folder names.
    Only extracts location and trial information - session start time comes from XML metadata.

    Expected folder name format: [date]_Cell[cell_number]_[location][location_number][variant]_[experiment_type]-[trial_number]
    Examples:
    - 05262016_Cell1_dist1_trio-001 (MMDDYYYY format)
    - 05262016_Cell2_prox2_trio-003 (MMDDYYYY format with different cell/location)

    Date formats:
    - YYYYMMDD: Year, month, day (e.g., 20160523)
    - MMDDYYYY: Month, day, year (e.g., 05232016)
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
    # Format: [date]_Cell[cell_number]_[location][location_number][variant]_[experiment_type][_SCH]-[trial_number]
    pattern = (
        r"(?P<date>\d{8})"  # Date (8 digits): YYYYMMDD or MMDDYYYY
        r"_[Cc]ell(?P<cell_number>\d+)"  # Cell number: Cell1, cell2, etc.
        r"_(?:(?P<variant_prefix>real))?"  # Optional variant prefix: "real" before location
        r"(?P<location>dist|prox|)"  # Location: "dist" (distal) or "prox" (proximal)
        r"(?P<location_number>\d+)"  # Location number: 1, 2, etc.
        r"(?:(?P<variant_suffix>real))?"  # Optional variant suffix: "real" after location
        r"(?:_(?P<experiment_type>[Tt]rio|bsl))?"  # Optional experiment type: "trio"/"Trio" or "bsl"
        r"(?:_(?P<sch_suffix>SCH))?"  # Optional SCH suffix for SCH condition
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
    Convert all Figure 7 dendritic excitability recordings from a session folder to a single NWB format file.

    Parameters
    ----------
    session_folder_path : Path
        Path to the session folder (corresponds to a single animal/subject)
    condition : str
        Experimental condition (e.g., "KO off-state", "KO on-state")
    verbose : bool, default=False
        Enable verbose output showing detailed processing information

    Returns
    -------
    NWBFile
        NWB file with all the converted data from the session

    Notes
    -----
    Each session folder corresponds to a single animal. The structure is:
    session_folder/recording_folder/files
    where recording_folder represents different recordings from the same animal.
    """

    print(f"Processing session folder: {session_folder_path.name} (corresponds to one animal)")

    # Get all recording folders within the session folder
    all_recording_folders = [f for f in session_folder_path.iterdir() if f.is_dir()]
    all_recording_folders.sort()

    if not all_recording_folders:
        raise ValueError(f"No recording folders found in session folder: {session_folder_path}")

    print(f"  Total recordings found: {len(all_recording_folders)}")

    # Calculate recording IDs, session start times, and create interface mappings
    ophys_session_start_times = []  # (ophys_time, recording_folder, recording_id)
    intracellular_session_start_times = []  # (intracellular_time, recording_folder, recording_id)
    recording_id_to_location_id = {}
    recording_id_to_folder = {}
    t_starts = {}  # t_starts[recording_id][interface] = t_start_offset

    if verbose:
        print(f"  Validating session start times and calculating recording IDs...")

    for recording_folder in all_recording_folders:
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

        # Compare session start times
        time_diff = abs((ophys_session_start_time - intracellular_session_start_time).total_seconds())
        print(f"    Times for {recording_folder.name}:")
        print(f"      Ophys time: {ophys_session_start_time}")
        print(f"      Intracellular time: {intracellular_session_start_time}")
        print(f"      Difference: {time_diff:.1f} seconds")

        # Get unique identifiers for recording to name objects
        recording_info = parse_session_info_from_folder_name(recording_folder)
        repetition_id = f"{recording_info['base_line_experiment_type']}Trial{recording_info['trial_number']}{recording_info['variant']}"
        location_id = f"Cell{recording_info['cell_number']}{recording_info['location_full']}Dendrite{recording_info['location_number']}"
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

    print(f"  Overall session start time: {session_start_time}")
    print(f"    Earliest time source: {earliest_interface} interface from recording {earliest_folder.name}")
    print(f"    Earliest line scan ophys time: {earliest_ophys_time}")
    print(f"    Earliest intracellular electrophysiology time: {earliest_intracellular_time}")

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

        if verbose:
            print(f"    Recording {folder.name} ({recording_id}) temporal alignment:")
            print(f"      Line scan interfaces (structural Ch1 + calcium Ch2) t_start = {ophys_t_start:.3f} seconds")
            print(f"      Intracellular electrophysiology interface t_start = {intracellular_t_start:.3f} seconds")

    # Get first recording info for session description
    first_recording_folder = next(iter(recording_id_to_folder.values()))
    first_recording_info = parse_session_info_from_folder_name(first_recording_folder)

    # Load metadata from YAML file
    metadata_file_path = Path(__file__).parent / "metadata.yaml"
    paper_metadata = load_dict_from_file(metadata_file_path)

    # Create session-specific metadata using session start time from XML
    session_date_str = session_start_time.strftime("%Y-%m-%d")

    session_specific_metadata = {
        "NWBFile": {
            "session_description": (
                f"Figure 7 dendritic excitability assessment in CDGI knockout mice for condition '{condition}'. "
                f"Combined patch clamp electrophysiology and two-photon laser scanning microscopy (2PLSM). "
                f"Brief current steps (three 2 nA injections, 2 ms each, at 50 Hz) with simultaneous Ca2+ imaging "
                f"to assess back-propagating action potential invasion. CDGI knockout disrupts M1 muscarinic signaling "
                f"pathway to test role in dendritic excitability adaptations. Recorded on {session_date_str}. "
                f"Total recordings: {len(all_recording_folders)}."
            ),
            "identifier": f"zhai2025_fig7_dendritic_{session_folder_path.name}_{condition.replace(' ', '_')}",
            "session_start_time": session_start_time,
            "experiment_description": (
                f"Figure 7 dendritic excitability changes in CDGI knockout mice during condition '{condition}'. "
                f"This experiment is part of Figure 7 from Zhai et al. 2025, investigating the role of M1 muscarinic "
                f"signaling in LID-induced dendritic adaptations using CDGI knockout mice. CDGI (Calmodulin-dependent "
                f"protein kinase kinase) is essential for M1R-mediated adenylyl cyclase inhibition. Methodology: "
                f"combination of patch clamp electrophysiology and 2PLSM. Cells filled with Ca2+-sensitive dye Fluo-4 "
                f"(100 μM) and Ca2+-insensitive dye Alexa Fluor 568 (50 μM). Somatically delivered current steps evoke "
                f"spikes that back-propagate into dendrites. Ca2+ signals at proximal (~40 μm) and distal (~90 μm) "
                f"locations serve as surrogate estimate of dendritic depolarization extent. CDGI KO mice: conditional "
                f"knockout targeting striatal spiny projection neurons."
            ),
            "session_id": f"{session_folder_path.name}_{condition.replace(' ', '_')}",
            "keywords": [
                "calcium imaging",
                "dendritic excitability",
                "current injection",
                "back-propagating action potentials",
                "CDGI knockout",
                "M1 muscarinic receptor",
                "levodopa-induced dyskinesia",
                "Fluo-4",
                "Alexa Fluor 568",
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

    # Create subject metadata for CDGI knockout experiments (Figure 7)
    subject = Subject(
        subject_id=f"CDGI_KO_mouse_{session_folder_path.name}",
        species="Mus musculus",
        strain="C57Bl/6",
        description=(
            f"CDGI conditional knockout mouse with unilateral 6-OHDA lesion in the medial forebrain bundle (MFB). "
            f"CDGI (Calmodulin-dependent protein kinase kinase) knockout disrupts M1 muscarinic receptor signaling "
            f"pathway in striatal spiny projection neurons. Animal {session_folder_path.name} recorded on {session_date_str}. "
            f"Lesion assessment: drug-free forelimb-use asymmetry test (cylinder test). LID induction: dyskinesiogenic "
            f"doses of levodopa (6 mg/kg first two sessions, 12 mg/kg later sessions, supplemented with 12 mg/kg "
            f"benserazide) every other day for at least five sessions. Recording condition: {condition}."
        ),
        genotype="CDGI conditional knockout (Camk2g-flox/flox; Dlx5/6-Cre)",
        sex="M",
        age="P49-P84",  # ISO format for 7-12 weeks (postnatal days)
    )
    nwbfile.subject = subject

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
    intracellular_recording_table.add_column(
        name="cdgi_genotype", description="CDGI knockout genotype status for this recording"
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

        if verbose:
            print(f"  Processing recording for folder: {recording_folder.name}")
            print(f"    Recording ID: {recording_id}")
            print(f"    Location ID: {location_id}")
            print(
                f"    Line scan structural channel (Ch1/Alexa568) temporal alignment offset: {t_starts[recording_id]['line_scan_structural_channel']:.3f} seconds"
            )
            print(
                f"    Line scan calcium channel (Ch2/Fluo4) temporal alignment offset: {t_starts[recording_id]['line_scan_calcium_channel']:.3f} seconds"
            )
            print(
                f"    Intracellular electrophysiology temporal alignment offset: {t_starts[recording_id]['intracellular']:.3f} seconds"
            )

        # Find electrophysiology XML file (exact name from Figure 7 notes)
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

        # Update electrode description for CDGI knockout dendritic recording (Figure 7)
        # One electrode per cell and location combination
        electrode_name = f"IntracellularElectrode{location_id}"
        intracellular_metadata["Icephys"]["IntracellularElectrodes"][icephys_metadata_key].update(
            {
                "name": electrode_name,
                "description": (
                    f"Recording from CDGI knockout mouse {recording_info['location_description']} - {condition} - "
                    f"Cell {recording_info['cell_number']} - Trial {recording_info['trial_number']} - "
                    f"Brief current steps (three 2 nA injections, 2 ms each, at 50 Hz) with simultaneous "
                    f"two-photon line scan imaging of calcium transients. "
                    f"CDGI genotype: Conditional knockout (Camk2g-flox/flox; Dlx5/6-Cre)"
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
                    f"Current clamp recording from CDGI knockout mouse {recording_info['location_description']} - "
                    f"{condition} - Cell {recording_info['cell_number']} - Trial {recording_info['trial_number']} - "
                    f"Three 2 nA current injections, 2 ms each, at 50 Hz. Stimulus protocol: "
                    f"PulseCount=3, PulseWidth=2ms, PulseSpacing=18ms (50Hz), FirstPulseDelay=900ms. "
                    f"CDGI knockout disrupts M1R-mediated adenylyl cyclase inhibition pathway."
                ),
            }
        )

        # Add intracellular data to NWB file
        intracellular_interface.add_to_nwbfile(nwbfile=nwbfile, metadata=intracellular_metadata)

        # Add intracellular recording to icephys table with custom annotations
        current_clamp_series = nwbfile.acquisition[series_name]

        # Calculate dendrite distance based on location type (approximate values from literature)
        dendrite_distance_um = 90 if recording_info["location"] == "dist" else 40  # distal ~90μm, proximal ~40μm

        # Determine CDGI genotype status
        cdgi_genotype = "CDGI KO (Camk2g-flox/flox; Dlx5/6-Cre)"

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
            cdgi_genotype=cdgi_genotype,
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
        structural_metadata["Ophys"]["ImagingPlanes"][structural_ophys_key]["name"] = f"ImagingPlane{location_id}"
        structural_metadata["Ophys"]["ImagingPlanes"][structural_ophys_key][
            "description"
        ] = f"Line scan imaging plane for {location_id} using Alexa Fluor 568 structural dye. Line scan parameters: 64 pixels per line, 10 μs dwell time, ~640 μs per line."
        structural_metadata["Ophys"]["ImagingPlanes"][structural_ophys_key]["indicator"] = "Alexa Fluor 568"
        structural_metadata["Ophys"]["PlaneSegmentation"][structural_ophys_key][
            "name"
        ] = f"PlaneSegmentation{recording_id}"
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
        ophys_metadata["ImagingPlanes"][calcium_ophys_key]["name"] = f"ImagingPlane{location_id}"
        ophys_metadata["ImagingPlanes"][calcium_ophys_key][
            "description"
        ] = f"Line scan imaging plane for {location_id} using Fluo-4 calcium indicator. Line scan parameters: 64 pixels per line, 10 μs dwell time, ~640 μs per line. Temporal resolution: ~1.6 seconds for 2500 lines."
        ophys_metadata["ImagingPlanes"][calcium_ophys_key]["indicator"] = "Fluo-4"
        ophys_metadata["PlaneSegmentation"][calcium_ophys_key]["name"] = f"PlaneSegmentation{recording_id}"
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

        if verbose:
            print(f"    Added line scan imaging data")
            print(f"    Successfully processed recording: {recording_folder.name}")

    print(f"Successfully processed all recordings from session: {session_folder_path.name}")

    # Build icephys table hierarchical structure following PyNWB best practices
    if verbose:
        print(f"  Building icephys table structure for {len(recording_indices)} recordings...")

    # Step 1: Build simultaneous recordings (each trial is its own simultaneous group)
    simultaneous_recording_indices = []
    for recording_index in recording_indices:
        metadata = recording_to_metadata[recording_index]
        simultaneous_index = nwbfile.add_icephys_simultaneous_recording(
            recordings=[recording_index]  # Each trial is its own simultaneous group
        )
        simultaneous_recording_indices.append(simultaneous_index)

    # Step 2: Build sequential recordings (each trial is its own sequence)
    for simultaneous_index in simultaneous_recording_indices:
        sequential_index = nwbfile.add_icephys_sequential_recording(
            simultaneous_recordings=[simultaneous_index],  # Each trial is its own sequence as requested
            stimulus_type="dendritic_excitability_current_injection",
        )
        sequential_recording_indices.append(sequential_index)

    # Step 3: Build repetitions table (group trials by dendritic location)
    repetitions_table = nwbfile.get_icephys_repetitions()
    repetitions_table.add_column(
        name="dendrite_distance_um", description="Approximate distance from soma in micrometers for this location"
    )
    repetitions_table.add_column(name="dendrite_type", description="Type of dendritic location: Distal or Proximal")
    repetitions_table.add_column(
        name="dendrite_number", description="Number identifier for the specific dendritic location"
    )
    repetitions_table.add_column(name="cdgi_genotype", description="CDGI knockout genotype status for this repetition")

    repetition_indices = []
    for location_id, location_recording_indices in location_to_recording_indices.items():
        # Get metadata from first recording at this location for location info
        first_recording_index = location_recording_indices[0]
        first_metadata = recording_to_metadata[first_recording_index]
        recording_info = first_metadata["recording_info"]

        # Get corresponding sequential recording indices for this location
        location_sequential_indices = []
        for recording_index in location_recording_indices:
            # Find the sequential index that corresponds to this recording
            seq_index = recording_indices.index(recording_index)
            location_sequential_indices.append(sequential_recording_indices[seq_index])

        dendrite_distance_um = 90 if recording_info["location"] == "dist" else 40
        cdgi_genotype = "CDGI KO (Camk2g-flox/flox; Dlx5/6-Cre)"

        repetition_index = nwbfile.add_icephys_repetition(
            sequential_recordings=location_sequential_indices,
            dendrite_distance_um=dendrite_distance_um,
            dendrite_type=recording_info["location_full"],  # "Distal" or "Proximal"
            dendrite_number=int(recording_info["location_number"]),
            cdgi_genotype=cdgi_genotype,
        )
        repetition_indices.append(repetition_index)

    # Step 4: Build experimental conditions table (group all repetitions by CDGI knockout condition)
    # Note: Cell information is per-session data, stored at NWBFile level, not per-recording
    experimental_conditions_table = nwbfile.get_icephys_experimental_conditions()
    experimental_conditions_table.add_column(
        name="condition", description="Experimental condition for CDGI knockout study in LID"
    )
    experimental_conditions_table.add_column(name="cdgi_genotype", description="CDGI knockout genotype details")

    cdgi_genotype = "CDGI KO (Camk2g-flox/flox; Dlx5/6-Cre)"
    experimental_condition_index = nwbfile.add_icephys_experimental_condition(
        repetitions=repetition_indices, condition=condition, cdgi_genotype=cdgi_genotype
    )

    if verbose:
        print(f"    Added experimental condition '{condition}' with {len(repetition_indices)} repetitions")
        print(f"  Successfully built icephys table hierarchy:")
        print(f"    - {len(recording_indices)} intracellular recordings")
        print(f"    - {len(simultaneous_recording_indices)} simultaneous recordings")
        print(f"    - {len(sequential_recording_indices)} sequential recordings")
        print(f"    - {len(repetition_indices)} repetitions (grouped by location)")
        print(f"    - 1 experimental condition ('{condition}')")

    return nwbfile


if __name__ == "__main__":
    import logging

    from tqdm import tqdm

    # Control verbose output from here
    VERBOSE = True  # Set to True for detailed output

    # Suppress tifffile warnings
    logging.getLogger("tifffile").setLevel(logging.ERROR)

    # Suppress specific warnings
    warnings.filterwarnings("ignore", message="invalid value encountered in divide")
    warnings.filterwarnings("ignore", message=".*no datetime before year 1.*")

    # Set the base path to your data
    base_path = Path("./link_to_raw_data/Figure 7/KO DE on vs off")

    # Create nwb_files directory at root level
    root_dir = Path(__file__).parent.parent.parent.parent  # Go up to repo root
    nwb_files_dir = root_dir / "nwb_files" / "figure_7_dendritic_excitability"
    nwb_files_dir.mkdir(parents=True, exist_ok=True)

    # Figure 7 dendritic conditions
    conditions = ["KO off-state", "KO on-state"]

    for condition in conditions:
        condition_path = base_path / condition

        if not condition_path.exists():
            raise FileNotFoundError(f"Expected condition path does not exist: {condition_path}")

        print(f"Processing Figure 7 dendritic excitability data for: {condition}")

        # Get all session folders (e.g., 0525a, 0525b, etc.)
        # Note: Each session folder corresponds to a single animal
        session_folders = [f for f in condition_path.iterdir() if f.is_dir()]
        session_folders.sort()

        print(f"Found {len(session_folders)} session folders")

        # Process each session folder with progress bar
        with tqdm(session_folders, desc=f"Processing {condition}", unit="session", position=0, leave=True) as pbar:
            for session_folder in pbar:
                pbar.write(f"\nProcessing session folder: {session_folder.name}")

                # Convert all recordings from this session to NWB format
                nwbfile = convert_session_to_nwbfile(
                    session_folder_path=session_folder,
                    condition=condition,
                    verbose=VERBOSE,
                )

                # Create output filename
                condition_safe = condition.replace(" ", "_").replace("(", "").replace(")", "")
                nwbfile_path = (
                    nwb_files_dir / f"figure7_dendritic_excitability_{condition_safe}_{session_folder.name}.nwb"
                )

                # Write NWB file
                configure_and_write_nwbfile(nwbfile, nwbfile_path=nwbfile_path)
                pbar.write(f"Successfully saved: {nwbfile_path.name}")
