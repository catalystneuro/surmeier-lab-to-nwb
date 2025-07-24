"""
Utility functions for somatic excitability experiments.

This module contains shared utility functions used across multiple somatic excitability
conversion scripts to reduce code duplication and ensure consistency.
"""

from pynwb import NWBFile


def build_somatic_icephys_table_structure(
    nwbfile: NWBFile,
    recording_indices: list[int],
    condition: str,
) -> int:
    """
    Build icephys table hierarchical structure following PyNWB best practices for somatic excitability experiments.

    This function creates the complete icephys table hierarchy:
    1. Simultaneous recordings (each current step is its own simultaneous group)
    2. Sequential recordings (group all current steps for this session)
    3. Repetitions table (single repetition per session)
    4. Experimental conditions table (group by LID condition)

    Parameters
    ----------
    nwbfile : NWBFile
        The NWB file to add the icephys table structure to
    recording_indices : list[int]
        List of intracellular recording indices from nwbfile.add_intracellular_recording calls
    condition : str
        Experimental condition (e.g., "LID off-state", "LID on-state", etc.)

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

    # Step 2: Build sequential recording (group all current steps for this session)
    sequential_index = nwbfile.add_icephys_sequential_recording(
        simultaneous_recordings=simultaneous_recording_indices,
        stimulus_type="F-I_protocol_somatic_excitability",
    )
    sequential_recording_indices = [sequential_index]

    # Step 3: Build repetitions table (single repetition per session)
    repetition_index = nwbfile.add_icephys_repetition(sequential_recordings=sequential_recording_indices)

    # Step 4: Build experimental conditions table (group by experimental condition)
    experimental_conditions_table = nwbfile.get_icephys_experimental_conditions()
    experimental_conditions_table.add_column(
        name="condition", description="Experimental condition for L-DOPA induced dyskinesia study"
    )

    experimental_condition_index = nwbfile.add_icephys_experimental_condition(
        repetitions=[repetition_index], condition=condition
    )

    return experimental_condition_index
