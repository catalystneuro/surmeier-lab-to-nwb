"""
Utility functions for somatic excitability experiments.

This module contains shared utility functions used across multiple somatic excitability
conversion scripts to reduce code duplication and ensure consistency.
"""

from typing import Any, Dict, List

from pynwb import NWBFile


def build_somatic_icephys_table_structure(
    nwbfile: NWBFile,
    recording_indices: List[int],
    session_info: Dict[str, Any],
    condition: str,
    stimulus_type: str,
    include_animal_letter: bool = False,
    animal_id_key: str = "animal_letter",
) -> int:
    """
    Build icephys table hierarchical structure following PyNWB best practices for somatic excitability experiments.

    This function creates the complete icephys table hierarchy:
    1. Simultaneous recordings (each current step is its own simultaneous group)
    2. Sequential recordings (group all current steps for this cell)
    3. Repetitions table (for this single cell/animal, it's just one repetition)
    4. Experimental conditions table (group by LID condition)

    Parameters
    ----------
    nwbfile : NWBFile
        The NWB file to add the icephys table structure to
    recording_indices : List[int]
        List of intracellular recording indices from nwbfile.add_intracellular_recording calls
    session_info : Dict[str, Any]
        Session information dictionary containing at least 'cell_number' and optionally 'animal_letter'
    condition : str
        Experimental condition (e.g., "LID off-state", "LID on-state", etc.)
    stimulus_type : str
        Type of stimulus protocol (e.g., "F-I_protocol_somatic_excitability")
    include_animal_letter : bool, default=False
        Whether to include animal identifier column in repetitions table
    animal_id_key : str, default="animal_letter"
        The key in session_info for the animal identifier ("animal_letter" or "animal_id")

    Returns
    -------
    int
        The experimental condition index from the final table level

    Notes
    -----
    This function assumes each current injection step is recorded as a separate simultaneous group,
    which is typical for F-I relationship protocols in somatic excitability experiments.
    """

    # Step 1: Build simultaneous recordings (each current step is its own simultaneous group)
    simultaneous_recording_indices = []
    for recording_index in recording_indices:
        simultaneous_index = nwbfile.add_icephys_simultaneous_recording(
            recordings=[recording_index]  # Each current step is its own simultaneous group
        )
        simultaneous_recording_indices.append(simultaneous_index)

    # Step 2: Build sequential recording (group all current steps for this cell)
    sequential_index = nwbfile.add_icephys_sequential_recording(
        simultaneous_recordings=simultaneous_recording_indices,
        stimulus_type=stimulus_type,
    )
    sequential_recording_indices = [sequential_index]

    # Step 3: Build repetitions table (for this single cell, it's just one repetition)
    repetitions_table = nwbfile.get_icephys_repetitions()
    repetitions_table.add_column(name="cell_number", description="Cell number identifier for this recording session")

    # Add animal identifier column if requested (used by some figures)
    if include_animal_letter:
        animal_column_name = animal_id_key
        animal_description = (
            "Animal identifier letter for this experimental session"
            if animal_id_key == "animal_letter"
            else "Animal identifier for tracking across experimental sessions"
        )
        repetitions_table.add_column(name=animal_column_name, description=animal_description)

    # Prepare repetition arguments
    repetition_kwargs = {
        "sequential_recordings": sequential_recording_indices,
        "cell_number": int(session_info["cell_number"]),
    }

    # Add animal identifier if required
    if include_animal_letter:
        repetition_kwargs[animal_id_key] = session_info[animal_id_key]

    repetition_index = nwbfile.add_icephys_repetition(**repetition_kwargs)

    # Step 4: Build experimental conditions table (group by LID condition)
    experimental_conditions_table = nwbfile.get_icephys_experimental_conditions()
    experimental_conditions_table.add_column(
        name="condition", description="Experimental condition for L-DOPA induced dyskinesia study"
    )

    experimental_condition_index = nwbfile.add_icephys_experimental_condition(
        repetitions=[repetition_index], condition=condition
    )

    return experimental_condition_index
