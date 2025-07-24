"""
Utility functions for optical stimulation experiments.

This module contains shared utility functions used across multiple optical stimulation
conversion scripts to reduce code duplication and ensure consistency.
"""

from typing import Any, Dict, List

from pynwb import NWBFile


def build_optical_icephys_table_structure(
    nwbfile: NWBFile,
    recording_indices: List[int],
    session_info: Dict[str, Any],
    condition: str,
    stimulus_type: str = "Sr2+_oEPSC_optogenetic_protocol",
) -> int:
    """
    Build icephys table hierarchical structure following PyNWB best practices for optical stimulation experiments.

    This function creates the complete icephys table hierarchy:
    1. Simultaneous recordings (each sweep is its own simultaneous group)
    2. Sequential recordings (group all sweeps for this session)
    3. Repetitions table (for this session, it's one repetition)
    4. Experimental conditions table (group by LID condition)

    Parameters
    ----------
    nwbfile : NWBFile
        The NWB file to add the icephys table structure to
    recording_indices : List[int]
        List of intracellular recording indices from nwbfile.add_intracellular_recording calls
    session_info : Dict[str, Any]
        Session information dictionary containing at least 'session_letter'
    condition : str
        Experimental condition (e.g., "LID off-state", "LID on-state")
    stimulus_type : str, default="Sr2+_oEPSC_optogenetic_protocol"
        Type of stimulus protocol

    Returns
    -------
    int
        The experimental condition index from the final table level

    Notes
    -----
    This function assumes each sweep is recorded as a separate simultaneous group,
    which is typical for Sr²⁺-oEPSC protocols in optical stimulation experiments.
    """

    # Step 1: Build simultaneous recordings (each sweep is its own simultaneous group)
    simultaneous_recording_indices = []
    for recording_index in recording_indices:
        simultaneous_index = nwbfile.add_icephys_simultaneous_recording(
            recordings=[recording_index]  # Each sweep is its own simultaneous group
        )
        simultaneous_recording_indices.append(simultaneous_index)

    # Step 2: Build sequential recording (group all sweeps for this session)
    sequential_index = nwbfile.add_icephys_sequential_recording(
        simultaneous_recordings=simultaneous_recording_indices,
        stimulus_type=stimulus_type,
    )

    # Step 3: Build repetitions table (for this session, it's one repetition)
    repetitions_table = nwbfile.get_icephys_repetitions()
    repetitions_table.add_column(
        name="session_letter", description="Session identifier letter for this experimental day"
    )

    repetition_index = nwbfile.add_icephys_repetition(
        sequential_recordings=[sequential_index],
        session_letter=session_info["session_letter"],
    )

    # Step 4: Build experimental conditions table (group by LID condition)
    experimental_conditions_table = nwbfile.get_icephys_experimental_conditions()
    experimental_conditions_table.add_column(
        name="condition", description="Experimental condition for Sr²⁺-oEPSC study"
    )

    experimental_condition_index = nwbfile.add_icephys_experimental_condition(
        repetitions=[repetition_index], condition=condition
    )

    return experimental_condition_index
