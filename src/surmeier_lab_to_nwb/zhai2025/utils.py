"""Utility functions for NWB conversion."""

from pynwb.base import TimeSeries
from pynwb.file import NWBFile


def add_dendritic_trials_table(
    nwbfile: NWBFile, recording_indices: list, recording_to_metadata: dict, t_starts: dict, verbose: bool = False
) -> None:
    """
    Add trials table for dendritic excitability experiments to NWB file.

    This function adds a standardized trials table that captures the experimental
    protocol for dendritic excitability experiments: 3x 2nA current injections,
    2ms each, at 50Hz frequency, with simultaneous two-photon line scan imaging.

    Parameters
    ----------
    nwbfile : NWBFile
        NWB file to add trials table to
    recording_indices : list
        List of recording indices from intracellular recordings
    recording_to_metadata : dict
        Dictionary mapping recording index to metadata
    t_starts : dict
        Dictionary mapping recording_id to timing information
    verbose : bool, default=False
        Enable verbose output
    """
    if verbose:
        print(f"  Building trials table for {len(recording_indices)} trials...")

    # Add trials table columns specific to dendritic excitability
    nwbfile.add_trial_column(
        name="stimulus_current_nA",
        description="Current amplitude of stimulus injections in nanoamps (2 nA for all trials)",
    )
    nwbfile.add_trial_column(
        name="stimulus_pulse_width_ms",
        description="Duration of each current pulse in milliseconds (2 ms for all trials)",
    )
    nwbfile.add_trial_column(
        name="stimulus_pulse_count", description="Number of current pulses in the stimulus train (3 for all trials)"
    )
    nwbfile.add_trial_column(
        name="stimulus_frequency_Hz", description="Frequency of current pulse train in Hz (50 Hz for all trials)"
    )
    nwbfile.add_trial_column(
        name="dendrite_location",
        description="Dendritic location of recording (proximal ~40μm or distal ~90μm from soma)",
    )
    nwbfile.add_trial_column(
        name="dendrite_distance_um", description="Approximate distance of recording site from soma in micrometers"
    )
    nwbfile.add_trial_column(name="trial_number", description="Trial number for this dendritic location (1, 2, 3)")
    nwbfile.add_trial_column(
        name="recording_id", description="Full recording identifier containing location and trial information"
    )
    nwbfile.add_trial_column(name="cell_number", description="Cell number identifier")

    # Collect trial data and sort by start_time to ensure chronological order
    trial_data_list = []

    for recording_index in recording_indices:
        recording_metadata = recording_to_metadata[recording_index]
        recording_info = recording_metadata["recording_info"]

        # Calculate trial timing from t_starts (already aligned to session start)
        recording_id = recording_metadata["recording_id"]
        trial_start_time = t_starts[recording_id]["intracellular"]

        # Estimate trial duration from intracellular recording
        # Typical duration is ~3-5 seconds for dendritic recordings
        trial_duration = 4.0  # seconds (approximate based on protocol)
        trial_stop_time = trial_start_time + trial_duration

        # Calculate dendrite distance based on location type
        dendrite_distance_um = 90 if recording_info["location"] == "dist" else 40

        # Store trial data for sorting
        trial_data = {
            "start_time": trial_start_time,
            "stop_time": trial_stop_time,
            "stimulus_current_nA": 2.0,
            "stimulus_pulse_width_ms": 2.0,
            "stimulus_pulse_count": 3,
            "stimulus_frequency_Hz": 50.0,
            "dendrite_location": recording_info["location_full"],
            "dendrite_distance_um": dendrite_distance_um,
            "trial_number": int(recording_info["trial_number"]),
            "recording_id": recording_id,
            "cell_number": recording_info["cell_number"],
        }
        trial_data_list.append(trial_data)

    # Sort trials by start_time to ensure ascending chronological order
    trial_data_list.sort(key=lambda x: x["start_time"])

    # Add sorted trials to trials table
    for trial_data in trial_data_list:
        nwbfile.add_trial(**trial_data)

    if verbose:
        print(f"    Added {len(recording_indices)} trials to trials table")
        print(f"    Trial structure: 3 trials each at proximal and distal dendritic locations")
        print(f"    Stimulus protocol: 3x 2nA current injections, 2ms each, at 50Hz")


def time_series_duration(time_series: TimeSeries) -> float:
    """
    Calculate the duration of a TimeSeries using its timing information.

    If the TimeSeries has explicit timestamps, uses those to calculate duration.
    If it has starting_time and rate, calculates duration from data shape and rate.

    Parameters
    ----------
    time_series : TimeSeries
        The TimeSeries object to calculate duration for

    Returns
    -------
    float
        Duration in seconds
    """
    if hasattr(time_series, "timestamps") and time_series.timestamps is not None:
        # Use explicit timestamps
        timestamps = time_series.timestamps[:]
        return float(timestamps[-1] - timestamps[0])
    elif hasattr(time_series, "starting_time") and hasattr(time_series, "rate"):
        # Use starting_time and rate
        data_length = len(time_series.data)
        return (data_length - 1) / time_series.rate
    else:
        raise ValueError("TimeSeries must have either timestamps or starting_time+rate")
