"""
Utility functions for dendritic excitability experiments.

This module contains shared utility functions used across multiple dendritic excitability
conversion scripts to reduce code duplication and ensure consistency.
"""

from typing import Any, Dict, List

from pynwb import NWBFile

from surmeier_lab_to_nwb.zhai2025.interfaces.trials_interface import (
    DendriticTrialsInterface,
)


def build_icephys_table_structure(
    nwbfile: NWBFile,
    recording_indices: List[int],
    recording_to_metadata: Dict[int, Dict[str, Any]],
    t_starts: Dict[str, Dict[str, float]],
    session_info: Dict[str, Any],
    condition: str,
    stimulus_type: str = "3x2nA_2ms_50Hz_dendritic_excitability",
    include_animal_id: bool = False,
    animal_id_key: str = "animal_id",
    verbose: bool = False,
) -> int:
    """
    Build icephys table hierarchical structure and trials table for dendritic excitability experiments.

    This function creates the complete icephys table hierarchy and a trials table:
    1. Simultaneous recordings (each dendritic location/trial is its own simultaneous group)
    2. Sequential recordings (group all trials for this cell)
    3. Repetitions table (for this single cell/animal, it's just one repetition)
    4. Experimental conditions table (group by experimental condition)
    5. Trials table (detailed trial-by-trial information for dendritic protocol)

    Parameters
    ----------
    nwbfile : NWBFile
        The NWB file to add the icephys table structure to
    recording_indices : List[int]
        List of intracellular recording indices from nwbfile.add_intracellular_recording calls
    recording_to_metadata : Dict[int, Dict[str, Any]]
        Dictionary mapping recording index to metadata (recording_id, recording_info, series_name)
    t_starts : Dict[str, Dict[str, float]]
        Dictionary mapping recording_id to timing information with keys like 'intracellular'
    session_info : Dict[str, Any]
        Session information dictionary containing cell_number and optionally animal identifiers
    condition : str
        Experimental condition (e.g., "LID off-state", "LID on-state", etc.)
    stimulus_type : str, default="3x2nA_2ms_50Hz_dendritic_excitability"
        Type of stimulus protocol
    include_animal_id : bool, default=False
        Whether to include animal identifier column in repetitions table
    animal_id_key : str, default="animal_id"
        The key in session_info for the animal identifier ("animal_id" or "animal_letter")
    verbose : bool, default=False
        Enable verbose output

    Returns
    -------
    int
        The experimental condition index from the final table level

    Notes
    -----
    This function assumes each dendritic recording trial is recorded as a separate simultaneous group,
    which is typical for dendritic excitability protocols with multiple locations and trials.
    """

    if verbose:
        print(f"  Building icephys table structure for {len(recording_indices)} recordings...")

    # Step 1: Build simultaneous recordings (each dendritic trial is its own simultaneous group)
    simultaneous_recording_indices = []
    for recording_index in recording_indices:
        simultaneous_index = nwbfile.add_icephys_simultaneous_recording(
            recordings=[recording_index]  # Each dendritic trial is its own simultaneous group
        )
        simultaneous_recording_indices.append(simultaneous_index)

    # Step 2: Build sequential recording (group all trials for this cell)
    sequential_index = nwbfile.add_icephys_sequential_recording(
        simultaneous_recordings=simultaneous_recording_indices,
        stimulus_type=stimulus_type,
    )
    sequential_recording_indices = [sequential_index]

    # Step 3: Build repetitions table (for this single cell, it's just one repetition)
    repetitions_table = nwbfile.get_icephys_repetitions()
    repetitions_table.add_column(name="cell_number", description="Cell number identifier for this recording session")

    # Add animal identifier column if requested (used by some figures)
    if include_animal_id:
        animal_column_name = animal_id_key
        animal_description = (
            "Animal identifier for tracking across experimental sessions"
            if animal_id_key == "animal_id"
            else "Animal identifier letter for this experimental session"
        )
        repetitions_table.add_column(name=animal_column_name, description=animal_description)

    # Prepare repetition arguments
    repetition_kwargs = {
        "sequential_recordings": sequential_recording_indices,
        "cell_number": int(session_info["cell_number"]),
    }

    # Add animal identifier if required
    if include_animal_id:
        repetition_kwargs[animal_id_key] = session_info[animal_id_key]

    repetition_index = nwbfile.add_icephys_repetition(**repetition_kwargs)

    # Step 4: Build experimental conditions table (group by experimental condition)
    experimental_conditions_table = nwbfile.get_icephys_experimental_conditions()
    experimental_conditions_table.add_column(
        name="condition", description="Experimental condition for L-DOPA induced dyskinesia study"
    )

    experimental_condition_index = nwbfile.add_icephys_experimental_condition(
        repetitions=[repetition_index], condition=condition
    )

    # Step 5: Add trials table using DendriticTrialsInterface
    if verbose:
        print(f"  Building trials table for {len(recording_indices)} trials...")

    trials_interface = DendriticTrialsInterface(
        recording_indices=recording_indices, recording_to_metadata=recording_to_metadata, t_starts=t_starts
    )
    trials_interface.add_to_nwbfile(nwbfile, verbose=verbose)

    if verbose:
        print(f"    Completed icephys table structure and trials table")

    return experimental_condition_index
