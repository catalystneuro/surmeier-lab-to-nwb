"""Trials interface for dendritic excitability experiments."""

from neuroconv.basedatainterface import BaseDataInterface
from pynwb.file import NWBFile


class DendriticTrialsInterface(BaseDataInterface):
    """Interface for adding trials table to dendritic excitability experiments.

    This interface adds a standardized trials table that captures the experimental
    protocol for dendritic excitability experiments: 3x 2nA current injections,
    2ms each, at 50Hz frequency, with simultaneous two-photon line scan imaging.
    """

    def __init__(self, recording_indices: list, recording_to_metadata: dict, t_starts: dict):
        """
        Initialize the trials interface.

        Parameters
        ----------
        recording_indices : list
            List of recording indices from intracellular recordings
        recording_to_metadata : dict
            Dictionary mapping recording index to metadata
        t_starts : dict
            Dictionary mapping recording_id to timing information
        """
        self.recording_indices = recording_indices
        self.recording_to_metadata = recording_to_metadata
        self.t_starts = t_starts

    def add_to_nwbfile(self, nwbfile: NWBFile, verbose: bool = False) -> None:
        """
        Add trials table for dendritic excitability experiments to NWB file.

        Parameters
        ----------
        nwbfile : NWBFile
            NWB file to add trials table to
        verbose : bool, default=False
            Enable verbose output
        """
        if verbose:
            print(f"  Building trials table for {len(self.recording_indices)} trials...")

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

        for recording_index in self.recording_indices:
            recording_metadata = self.recording_to_metadata[recording_index]
            recording_info = recording_metadata["recording_info"]

            # Calculate trial timing from t_starts (already aligned to session start)
            recording_id = recording_metadata["recording_id"]
            trial_start_time = self.t_starts[recording_id]["intracellular"]

            # Estimate trial duration from intracellular recording
            # Typical duration is ~3-5 seconds for dendritic recordings
            trial_duration = 4.0  # seconds (approximate based on protocol)
            trial_stop_time = trial_start_time + trial_duration

            # Calculate dendrite distance based on location type
            # Handle different key names used across different scripts
            location_key = "location" if "location" in recording_info else "location_type"
            dendrite_distance_um = 90 if recording_info[location_key] == "dist" else 40

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
            print(f"    Added {len(self.recording_indices)} trials to trials table")
            print(f"    Trial structure: 3 trials each at proximal and distal dendritic locations")
            print(f"    Stimulus protocol: 3x 2nA current injections, 2ms each, at 50Hz")
