"""
Utility functions for dendritic excitability experiments.

This module contains shared utility functions used across multiple dendritic excitability
conversion scripts to reduce code duplication and ensure consistency.
"""

from pathlib import Path
from typing import Any, Dict, List

from neuroconv.tools.nwb_helpers import make_nwbfile_from_metadata
from neuroconv.utils import dict_deep_update, load_dict_from_file
from pynwb import NWBFile

from surmeier_lab_to_nwb.zhai2025.conversion_scripts.conversion_utils import (
    GENOTYPE_DESCRIPTION_MAPPING,
    PHARMACOLOGY_ADDITIONS,
    format_condition,
    generate_canonical_session_id,
)
from surmeier_lab_to_nwb.zhai2025.interfaces import (
    PrairieViewCurrentClampInterface,
    PrairieViewLineScanInterface,
)


def build_dendritic_icephys_table_structure(
    nwbfile: NWBFile,
    recording_indices: List[int],
    recording_to_metadata: Dict[int, Dict[str, Any]],
    condition: str,
    stimulus_type: str = "dendritic_excitability_current_injection",
    verbose: bool = False,
    **extra_condition_kwargs: Any,
) -> int:
    """
    Build icephys table hierarchical structure for dendritic excitability experiments.

    This function creates the complete icephys table hierarchy:
    1. Simultaneous recordings (each dendritic trial is its own simultaneous group)
    2. Sequential recordings (each trial is its own sequence as per existing implementation)
    3. Repetitions table (group trials by dendritic location)
    4. Experimental conditions table (group by experimental condition)

    Parameters
    ----------
    nwbfile : NWBFile
        The NWB file to add the icephys table structure to
    recording_indices : List[int]
        List of intracellular recording indices from nwbfile.add_intracellular_recording calls
    recording_to_metadata : Dict[int, Dict[str, Any]]
        Dictionary mapping recording index to metadata (recording_id, recording_info, series_name)
    session_info : Dict[str, Any]
        Session information dictionary (currently unused but maintained for API compatibility)
    condition : str
        Experimental condition (e.g., "LID off-state", "LID on-state", etc.)
    stimulus_type : str, default="dendritic_excitability_current_injection"
        Type of stimulus protocol
    verbose : bool, default=False
        Enable verbose output
    **extra_condition_kwargs : Any
        Additional keyword arguments to pass to add_icephys_experimental_condition

    Returns
    -------
    int
        The experimental condition index from the final table level

    Notes
    -----
    This function matches the existing complex dendritic excitability icephys table structure
    where each trial is its own sequence and repetitions are grouped by dendritic location.
    """

    if verbose:
        print(f"  Building icephys table structure for {len(recording_indices)} recordings...")

    # Step 1: Build simultaneous recordings (each trial is its own simultaneous group)
    simultaneous_recording_indices = []
    for recording_index in recording_indices:
        simultaneous_index = nwbfile.add_icephys_simultaneous_recording(
            recordings=[recording_index]  # Each trial is its own simultaneous group
        )
        simultaneous_recording_indices.append(simultaneous_index)

    # Step 2: Build sequential recordings (each trial is its own sequence)
    sequential_recording_indices = []
    for simultaneous_index in simultaneous_recording_indices:
        sequential_index = nwbfile.add_icephys_sequential_recording(
            simultaneous_recordings=[simultaneous_index],  # Each trial is its own sequence as requested
            stimulus_type=stimulus_type,
        )
        sequential_recording_indices.append(sequential_index)

    # Step 3: Group recordings by dendritic location for repetitions table
    location_to_recording_indices = {}  # Group recordings by location for repetitions table
    for recording_index in recording_indices:
        recording_metadata = recording_to_metadata[recording_index]
        recording_info = recording_metadata["recording_info"]

        # Create location identifier from location and number
        # Handle different field naming conventions
        if "location" in recording_info and "location_number" in recording_info:
            # Standard format (e.g., figure_1_dendritic_excitability.py)
            location_id = f"{recording_info['location']}_{recording_info['location_number']}"
        elif "location_type" in recording_info and "location_number" in recording_info:
            # OxoM format (e.g., figure_7_oxoM_dendritic_excitability.py)
            location_id = f"{recording_info['location_type']}_{recording_info['location_number']}"
        elif "location_id" in recording_info:
            # Already has location_id
            location_id = recording_info["location_id"]
        else:
            raise ValueError(f"Cannot determine location_id from recording_info: {recording_info.keys()}")

        if location_id not in location_to_recording_indices:
            location_to_recording_indices[location_id] = []
        location_to_recording_indices[location_id].append(recording_index)

    # Step 4: Build repetitions table (group trials by dendritic location)
    repetitions_table = nwbfile.get_icephys_repetitions()
    repetitions_table.add_column(
        name="dendrite_distance_um", description="Approximate distance from soma in micrometers for this location"
    )
    repetitions_table.add_column(name="dendrite_type", description="Type of dendritic location: Distal or Proximal")
    repetitions_table.add_column(
        name="dendrite_number", description="Number identifier for the specific dendritic location"
    )

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

        # Handle different field naming conventions for distance calculation
        if "location" in recording_info:
            # Standard format
            dendrite_distance_um = 90 if recording_info["location"] == "dist" else 40
            dendrite_type = recording_info.get(
                "location_full", "Distal" if recording_info["location"] == "dist" else "Proximal"
            )
        elif "location_type" in recording_info:
            # OxoM format
            dendrite_distance_um = 90 if recording_info["location_type"] == "dist" else 40
            dendrite_type = "Distal" if recording_info["location_type"] == "dist" else "Proximal"
        elif "approximate_distance_um" in recording_info:
            # Direct distance specification
            dendrite_distance_um = recording_info["approximate_distance_um"]
            # Infer type from distance
            dendrite_type = "Distal" if dendrite_distance_um > 60 else "Proximal"
        else:
            # Default values
            dendrite_distance_um = 40  # Default to proximal
            dendrite_type = "Proximal"

        # Get dendrite number
        dendrite_number = int(recording_info.get("location_number", 1))

        repetition_index = nwbfile.add_icephys_repetition(
            sequential_recordings=location_sequential_indices,
            dendrite_distance_um=dendrite_distance_um,
            dendrite_type=dendrite_type,
            dendrite_number=dendrite_number,
        )
        repetition_indices.append(repetition_index)

    # Step 5: Build experimental conditions table (group all repetitions by condition)
    experimental_conditions_table = nwbfile.get_icephys_experimental_conditions()
    experimental_conditions_table.add_column(
        name="condition", description="Experimental condition for L-DOPA induced dyskinesia study"
    )

    # Add any extra columns requested by specific scripts
    for key in extra_condition_kwargs:
        if key not in ["repetitions", "condition"]:
            if key == "m1r_treatment":
                experimental_conditions_table.add_column(
                    name=key, description="M1R treatment type (antagonist or control)"
                )
            elif key == "cdgi_genotype":
                experimental_conditions_table.add_column(
                    name=key, description="CDGI genotype for CDGI knockout experiments"
                )

    experimental_condition_index = nwbfile.add_icephys_experimental_condition(
        repetitions=repetition_indices, condition=condition, **extra_condition_kwargs
    )

    if verbose:
        print(f"    Completed icephys table structure with {len(repetition_indices)} repetitions")

    return experimental_condition_index


def convert_dendritic_excitability_session_to_nwbfile(
    session_folder_path: Path,
    condition: str,
    figure_config: Dict[str, Any],
    session_id_parameters: Dict[str, Any],
    verbose: bool = False,
) -> NWBFile:
    """
    Convert a single session of dendritic excitability data to NWB format with time alignment.

    This shared function handles the common conversion logic for all dendritic excitability
    experiments across different figures, with figure-specific variations handled through
    the configuration dictionary.

    Parameters
    ----------
    session_folder_path : Path
        Path to the session folder containing cell recordings
    condition : str
        Experimental condition (e.g., "LID off-state", "control", "KO on-state", etc.)
    figure_config : Dict[str, Any]
        Figure-specific configuration containing:
        - parse_function: Callable to parse folder names
        - metadata_key: str, key for session_specific_metadata.yaml
    session_id_parameters : Dict[str, Any]
        Session ID parameters using revised schema (excluding timestamp which is determined here):
        - fig: str, figure number (e.g., "F1", "F3")
        - meas_comp: str, measurement + compartment (e.g., "DendExc")
        - cell_type: str, cell type (e.g., "dSPN", "iSPN")
        - state: str, experimental state (e.g., "OffState", "OnState")
        - pharm: str, pharmacology condition (e.g., "none", "D1RaSch")
        - geno: str, genotype (e.g., "WT", "CDGIKO")
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
    """
    # Extract configuration
    parse_function = figure_config["parse_function"]
    metadata_key = figure_config["metadata_key"]

    # Get cell type from session ID parameters
    cell_type = session_id_parameters["cell_type"]

    # Get all recording folders within the session folder
    all_recording_folders = [f for f in session_folder_path.iterdir() if f.is_dir()]
    all_recording_folders.sort()

    if not all_recording_folders:
        raise ValueError(f"No recording folders found in session folder: {session_folder_path}")

    # Calculate recording IDs, session start times, and create interface mappings
    ophys_session_start_times = []  # (ophys_time, recording_folder, recording_id)
    intracellular_session_start_times = []  # (intracellular_time, recording_folder, recording_id)
    recording_id_to_location_id = {}
    recording_id_to_folder = {}
    t_starts = {}  # t_starts[recording_id][interface] = t_start_offset

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

        # Get unique identifiers for recording to name objects
        recording_info = parse_function(recording_folder)
        repetition_id = f"{recording_info['base_line_experiment_type']}Trial{recording_info['trial_number']}{recording_info['variant']}"
        location_id = f"{recording_info['location_full']}Dendrite{recording_info['location_number']}"
        recording_id = f"{location_id}{repetition_id}"

        # Store mappings
        recording_id_to_location_id[recording_id] = location_id
        recording_id_to_folder[recording_id] = recording_folder
        ophys_session_start_times.append((ophys_session_start_time, recording_folder, recording_id))
        intracellular_session_start_times.append((intracellular_session_start_time, recording_folder, recording_id))

    if not ophys_session_start_times or not intracellular_session_start_times:
        raise ValueError(f"No valid recordings found in session folder: {session_folder_path}")

    # Find the earliest session start time across all recordings (both ophys and intracellular)
    earliest_ophys_time = min(ophys_session_start_times, key=lambda x: x[0])[0]
    earliest_intracellular_time = min(intracellular_session_start_times, key=lambda x: x[0])[0]
    overall_earliest_time = min(earliest_ophys_time, earliest_intracellular_time)

    # Calculate t_start offsets for temporal alignment using specific interface names
    for ophys_time, _, recording_id in ophys_session_start_times:
        # Calculate offsets relative to overall session start time with improved precision
        ophys_offset = (ophys_time - overall_earliest_time).total_seconds()
        t_starts[recording_id] = {
            "line_scan_structural_channel": ophys_offset,
            "line_scan_calcium_channel": ophys_offset,  # Both line scan channels start together
        }

    for intracellular_time, _, recording_id in intracellular_session_start_times:
        intracellular_offset = (intracellular_time - overall_earliest_time).total_seconds()
        # Ensure the key exists before updating
        if recording_id not in t_starts:
            t_starts[recording_id] = {}
        t_starts[recording_id]["intracellular"] = intracellular_offset

    # Use overall earliest time as session start time for NWB file
    session_start_time = overall_earliest_time

    # Create canonical session ID with explicit parameters
    timestamp = session_start_time.strftime("%Y%m%d%H%M%S")

    # Use session ID parameters passed by the calling script (revised schema)
    session_id = generate_canonical_session_id(
        fig=session_id_parameters["fig"],
        meas_comp=session_id_parameters["meas_comp"],
        cell_type=session_id_parameters["cell_type"],
        state=session_id_parameters["state"],
        pharm=session_id_parameters["pharm"],
        geno=session_id_parameters["geno"],
        timestamp=timestamp,
    )

    # Load general and session-specific metadata from YAML files
    general_metadata_path = Path(__file__).parent.parent.parent / "general_metadata.yaml"
    general_metadata = load_dict_from_file(general_metadata_path)

    session_metadata_path = Path(__file__).parent.parent.parent / "session_specific_metadata.yaml"
    session_metadata_template = load_dict_from_file(session_metadata_path)
    script_template = session_metadata_template[metadata_key]

    # Get condition formatting for underscore version
    condition_underscore = format_condition[condition]["underscore"]

    # Build pharmacology description using centralized mapping
    pharmacology_base = general_metadata["NWBFile"]["pharmacology"]
    pharm_token = session_id_parameters["pharm"]

    # Add pharmacology-specific text if applicable
    if pharm_token != "none" and pharm_token in PHARMACOLOGY_ADDITIONS:
        pharmacology_text = pharmacology_base + " " + PHARMACOLOGY_ADDITIONS[pharm_token]
    else:
        pharmacology_text = pharmacology_base

    # Create session-specific metadata from template with runtime substitutions
    condition_human_readable = format_condition[condition]["human_readable"]
    session_specific_metadata = {
        "NWBFile": {
            "session_description": script_template["NWBFile"]["session_description"].format(
                condition=condition_human_readable
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

    # Process each recording using the calculated recording IDs
    for recording_id, recording_folder in recording_id_to_folder.items():
        location_id = recording_id_to_location_id[recording_id]
        recording_info = parse_function(recording_folder)
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

        # Find electrophysiology XML file (exact name from Figure 1 notes)
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

        # Build description based on figure-specific context using centralized mapping
        genotype = session_id_parameters["geno"]
        genotype_description = GENOTYPE_DESCRIPTION_MAPPING.get(genotype, "")

        # Update electrode description for dendritic recording
        # One electrode per cell and location combination
        electrode_name = f"IntracellularElectrode{location_id}"
        intracellular_metadata["Icephys"]["IntracellularElectrodes"][icephys_metadata_key].update(
            {
                "name": electrode_name,
                "description": (
                    f"Recording from {genotype_description}{cell_type} {recording_info['location_description']} - {condition_human_readable} - "
                    f"Trial {recording_info['trial_number']} - "
                    f"Brief current steps (three 2 nA injections, 2 ms each, at 50 Hz) with simultaneous "
                    f"two-photon line scan imaging of calcium transients"
                ),
                "cell_id": f"CellRecordedAt{timestamp}",
                "location": recording_info["location_description"],
                "slice": general_metadata["NWBFile"]["slices"],
            }
        )

        # Update current clamp series metadata
        series_name = f"CurrentClampSeries{recording_id}"
        series_description_prefix = f"{genotype_description}{cell_type}" if genotype != "WT" else cell_type
        intracellular_metadata["Icephys"]["CurrentClampSeries"][icephys_metadata_key].update(
            {
                "name": series_name,
                "description": (
                    f"Current clamp recording from {series_description_prefix} {recording_info['location_description']} - "
                    f"{condition_human_readable} - Trial {recording_info['trial_number']} - "
                    f"Three 2 nA current injections, 2 ms each, at 50 Hz. Stimulus protocol: "
                    f"PulseCount=3, PulseWidth=2ms, PulseSpacing=18ms (50Hz), FirstPulseDelay=900ms"
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
        ] = "Bruker Ultima two-photon microscope for line scan imaging. 810 nm excitation laser (Chameleon Ultra II, Coherent). Signals filtered at 2 kHz and digitized at 10 kHz."
        structural_metadata["Ophys"]["ImagingPlanes"][structural_ophys_key]["name"] = f"ImagingPlane{location_id}"
        structural_metadata["Ophys"]["ImagingPlanes"][structural_ophys_key][
            "description"
        ] = f"Line scan imaging plane for {location_id} using Alexa Fluor 568 structural dye. Line scan parameters: 64 pixels per line, 10 μs dwell time, ~640 μs per line."
        structural_metadata["Ophys"]["ImagingPlanes"][structural_ophys_key]["indicator"] = "Alexa Fluor 568"

        # Update PlaneSegmentation metadata for structural channel
        structural_metadata["Ophys"]["PlaneSegmentation"][structural_ophys_key][
            "name"
        ] = f"PlaneSegmentation{recording_id}"
        structural_metadata["Ophys"]["PlaneSegmentation"][structural_ophys_key][
            "description"
        ] = f"Line scan ROI segmentation for {location_id} structural imaging. Detected by Hamamatsu R3982 side-on PMT (580-620 nm)."

        # Update RoiResponseSeries metadata for structural channel
        structural_metadata["Ophys"]["RoiResponseSeries"][structural_ophys_key][
            "name"
        ] = f"RoiResponseSeriesAlexa568{recording_id}"
        structural_metadata["Ophys"]["RoiResponseSeries"][structural_ophys_key][
            "description"
        ] = f"Structural reference fluorescence from Alexa Fluor 568 hydrazide (50 μM) - {location_id}. Ca2+-insensitive dye to visualize dendrites."

        # Update SourceImages metadata for structural channel
        structural_metadata["Acquisition"]["SourceImages"][structural_ophys_key][
            "name"
        ] = f"ImageAlexa568{recording_id}"
        structural_metadata["Acquisition"]["SourceImages"][structural_ophys_key][
            "description"
        ] = f"Source image for Alexa Fluor 568 structural reference - {location_id}. Field of view with scan line overlay."

        # Update TimeSeries metadata for structural channel
        structural_metadata["TimeSeries"][structural_ophys_key]["name"] = f"TimeSeriesLineScanRawAlexa568{recording_id}"
        structural_metadata["TimeSeries"][structural_ophys_key][
            "description"
        ] = f"Line scan raw data for Alexa Fluor 568 structural reference - {location_id}. Typical acquisition: 2500 lines (time points)."

        # Update TwoPhotonSeries metadata for structural channel
        structural_metadata["Ophys"]["TwoPhotonSeries"][structural_ophys_key][
            "name"
        ] = f"TwoPhotonSeries{recording_id}Alexa568"
        structural_metadata["Ophys"]["TwoPhotonSeries"][structural_ophys_key][
            "description"
        ] = f"Line scan imaging data from {location_id} using Alexa Fluor 568 structural dye (Ch1). High-speed line scans across dendritic segments during current injection protocol. Acquisition: {recording_info['trial_number']} of dendritic excitability protocol."

        # Process calcium channel (Ch2/Fluo-4)
        calcium_metadata = calcium_interface.get_metadata()
        # Apply fluorophore-specific metadata
        calcium_metadata["Devices"][calcium_ophys_key]["name"] = "BrukerUltima"
        calcium_metadata["Devices"][calcium_ophys_key][
            "description"
        ] = "Bruker Ultima two-photon microscope for line scan imaging. 810 nm excitation laser (Chameleon Ultra II, Coherent). Signals filtered at 2 kHz and digitized at 10 kHz."
        calcium_metadata["Ophys"]["ImagingPlanes"][calcium_ophys_key]["name"] = f"ImagingPlane{location_id}"
        calcium_metadata["Ophys"]["ImagingPlanes"][calcium_ophys_key][
            "description"
        ] = f"Line scan imaging plane for {location_id} using Fluo-4 calcium indicator. Line scan parameters: 64 pixels per line, 10 μs dwell time, ~640 μs per line."
        calcium_metadata["Ophys"]["ImagingPlanes"][calcium_ophys_key]["indicator"] = "Fluo-4"

        # Update PlaneSegmentation metadata for calcium channel
        calcium_metadata["Ophys"]["PlaneSegmentation"][calcium_ophys_key]["name"] = f"PlaneSegmentation{recording_id}"
        calcium_metadata["Ophys"]["PlaneSegmentation"][calcium_ophys_key][
            "description"
        ] = f"Line scan ROI segmentation for {location_id} calcium imaging. Detected by Hamamatsu H7422P-40 GaAsP PMT (490-560 nm)."

        # Update RoiResponseSeries metadata for calcium channel
        calcium_metadata["Ophys"]["RoiResponseSeries"][calcium_ophys_key][
            "name"
        ] = f"RoiResponseSeriesFluo4{recording_id}"
        calcium_metadata["Ophys"]["RoiResponseSeries"][calcium_ophys_key][
            "description"
        ] = f"Calcium fluorescence from Fluo-4 (100 μM) - {location_id}. Ca2+-sensitive dye for measuring back-propagating action potential-evoked calcium transients. Magnitude serves as surrogate estimate of dendritic depolarization extent."

        # Update SourceImages metadata for calcium channel
        calcium_metadata["Acquisition"]["SourceImages"][calcium_ophys_key]["name"] = f"ImageFluo4{recording_id}"
        calcium_metadata["Acquisition"]["SourceImages"][calcium_ophys_key][
            "description"
        ] = f"Source image for Fluo-4 calcium indicator - {location_id}. Field of view with scan line overlay."

        # Update TimeSeries metadata for calcium channel
        calcium_metadata["TimeSeries"][calcium_ophys_key]["name"] = f"TimeSeriesLineScanRawFluo4{recording_id}"
        calcium_metadata["TimeSeries"][calcium_ophys_key][
            "description"
        ] = f"Line scan raw data for Fluo-4 calcium indicator - {location_id}. Typical acquisition: 2500 lines (time points). Kymograph structure: (C, T, X) where C=channels, T=time/lines, X=pixels along scan line."

        # Update TwoPhotonSeries metadata for calcium channel
        calcium_metadata["Ophys"]["TwoPhotonSeries"][calcium_ophys_key]["name"] = f"TwoPhotonSeries{recording_id}Fluo4"
        calcium_metadata["Ophys"]["TwoPhotonSeries"][calcium_ophys_key][
            "description"
        ] = f"Line scan imaging data from {location_id} using Fluo-4 calcium indicator (Ch2). High-speed line scans measuring calcium transients during current injection protocol. Acquisition: {recording_info['trial_number']} of dendritic excitability protocol."

        # Add both ophys channels to NWB file
        structural_interface.add_to_nwbfile(nwbfile=nwbfile, metadata=structural_metadata)
        calcium_interface.add_to_nwbfile(nwbfile=nwbfile, metadata=calcium_metadata)

    # Build icephys table hierarchical structure using shared utility function
    build_dendritic_icephys_table_structure(
        nwbfile=nwbfile,
        recording_indices=recording_indices,
        recording_to_metadata=recording_to_metadata,
        condition=condition_underscore,
        verbose=verbose,
    )

    if verbose:
        print(f"Created NWB file with {len(recording_indices)} recordings and session ID: {session_id}")

    return nwbfile
