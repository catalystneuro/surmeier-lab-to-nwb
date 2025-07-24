"""
Utility functions for dendritic excitability experiments.

This module contains shared utility functions used across multiple dendritic excitability
conversion scripts to reduce code duplication and ensure consistency.
"""

from typing import Any, Dict, List

from pynwb import NWBFile


def build_dendritic_icephys_table_structure(
    nwbfile: NWBFile,
    recording_indices: List[int],
    recording_to_metadata: Dict[int, Dict[str, Any]],
    session_info: Dict[str, Any],
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
