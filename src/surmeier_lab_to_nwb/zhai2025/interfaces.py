"""Primary classes for converting Prairie View data."""

from datetime import datetime
from pathlib import Path
from typing import Literal, Optional

import pandas as pd
import xmltodict
from neuroconv.basedatainterface import BaseDataInterface
from neuroconv.utils import DeepDict
from pynwb.file import NWBFile
from pynwb.icephys import CurrentClampSeries


class PrairieViewIntracellularRecordingInterface(BaseDataInterface):
    """Interface for Prairie View intracellular recording data."""

    keywords = ("intracellular", "recording", "prairie view")

    def __init__(
        self,
        xml_file_path: str | Path,
        recording_csv_file_path: Optional[str | Path] = None,
        cell_id: str = "CellID",
        trial_id: str = "",
        position: Literal["", "Proximal", "Distal"] = "",
        condition: Literal["", "LIDOffState", "LIDOnState", "LIDOnStateWithSCH"] = "",
    ):
        """
        Initialize the interface with paths to the XML and CSV files.

        Parameters
        ----------
        xml_file_path : Union[str, Path]
            Path to the XML file containing recording metadata
        recording_csv_file_path : Optional[Union[str, Path]], default=None
            Path to the CSV file containing recording data. If None, it will be inferred from the XML file.
        """
        self.xml_file_path = Path(xml_file_path)
        self.cell_id = cell_id
        self.trial_id = trial_id
        self.position = position
        self.condition = condition

        # Load XML data
        with open(self.xml_file_path, "r") as xml_file:
            xml_content = xml_file.read()
            self.xml_recording_dict = xmltodict.parse(xml_content)

        # Extract basic metadata to determine CSV file path
        vrec_session_entry = self.xml_recording_dict["VRecSessionEntry"]

        # Determine recording CSV file path if not provided
        if recording_csv_file_path is None:
            recording_data_file_path_stem = vrec_session_entry["DataFile"]
            self.recording_csv_file_path = self.xml_file_path.parent / f"{recording_data_file_path_stem}.csv"
        else:
            self.recording_csv_file_path = Path(recording_csv_file_path)

        if not self.recording_csv_file_path.is_file():
            raise FileNotFoundError(f"Recording data file not found: {self.recording_csv_file_path}")

    def get_metadata(self) -> DeepDict:
        """
        Extract metadata from the XML file.

        Returns
        -------
        DeepDict
            Dictionary containing metadata
        """
        metadata = super().get_metadata()

        # Extract basic metadata
        vrec_session_entry = self.xml_recording_dict["VRecSessionEntry"]
        experiment_metadata = vrec_session_entry["Experiment"]
        acquisition_time = experiment_metadata["AcquisitionTime"]
        sampling_rate = float(experiment_metadata["Rate"])

        # Extract signal information
        signal_list = experiment_metadata["SignalList"]
        assert len(signal_list) == 1, "Only one signal is expected in the XML file"
        signal_list = signal_list["VRecSignal"]
        enabled_signals = [signal for signal in signal_list if signal["Enabled"] == "true"]
        expected_signal = enabled_signals[0]

        # Extract unit metadata
        unit_metadata = expected_signal["Unit"]
        unit = unit_metadata["UnitName"]
        device = unit_metadata["PatchclampDevice"]

        # Extract session start time
        session_start_time = None
        if "DateTime" in vrec_session_entry:
            datetime_str = vrec_session_entry["DateTime"]
            try:
                # Parse ISO format datetime string
                session_start_time = datetime.fromisoformat(datetime_str.replace("Z", "+00:00"))
            except (ValueError, AttributeError):
                # If parsing fails, leave as None
                pass

        # Add recording-specific metadata
        metadata["Ecephys"] = {
            "Device": [{"name": device, "description": f"Prairie View intracellular recording device: {device}"}],
            "IntracelllularSeries": [
                {
                    "name": "PrairieViewIntracellularRecording",
                    "description": f"Intracellular recording from {device}",
                    "unit": unit,
                    "rate": sampling_rate,
                }
            ],
        }

        # Add session start time if available
        if session_start_time is not None:
            metadata["NWBFile"] = {"session_start_time": session_start_time}

        return metadata

    def add_to_nwbfile(self, nwbfile: NWBFile, metadata: Optional[dict] = None, electrode=None):
        """
        Add the recording data to the NWB file.

        Parameters
        ----------
        nwbfile : NWBFile
            NWB file to add the data to
        metadata : dict
            Metadata dictionary
        electrode : optional
            Electrode to use for the recording. If None, a new electrode will be created.
        """
        # Extract metadata from XML
        vrec_session_entry = self.xml_recording_dict["VRecSessionEntry"]
        experiment_metadata = vrec_session_entry["Experiment"]
        signal_list = experiment_metadata["SignalList"]["VRecSignal"]
        enabled_signals = [signal for signal in signal_list if signal["Enabled"] == "true"]
        expected_signal = enabled_signals[0]
        gain = float(expected_signal["Gain"])
        unit_metadata = expected_signal["Unit"]
        multiplier = float(unit_metadata["Multiplier"])
        divisor = float(unit_metadata["Divisor"])
        unit = unit_metadata["UnitName"]
        # device_name = unit_metadata["PatchclampDevice"]

        # Load recording data
        recording_data_df = pd.read_csv(
            self.recording_csv_file_path,
            header=0,
        )
        timestamps = recording_data_df["Time(ms)"] / 1_000.0  # Convert to seconds
        voltage_mV = recording_data_df[" Primary"] * multiplier / divisor
        voltage_volts = voltage_mV.values / 1_000.0  # Convert mV to volts

        # Create or get device
        device_name = "MultiClamp 700B"
        device_description = f"{device_name} amplifier (Axon Instruments/Molecular Devices) with signals filtered at 2 kHz and digitized at 10 kHz"
        if device_name in nwbfile.devices:
            device = nwbfile.devices[device_name]
        else:
            device = nwbfile.create_device(name=device_name, description=device_description)

        # Create or use provided electrode
        electrode_name = f"IntracellularElectrodeCondition{self.condition}Cell{self.cell_id}{self.position}"
        available_electrodes = nwbfile.icephys_electrodes
        if electrode_name in available_electrodes:
            electrode = available_electrodes[electrode_name]
        else:
            electrode = nwbfile.create_icephys_electrode(
                name=electrode_name,
                description=f"Prairie View electrode for {device_name}",
                device=device,
            )

        # Create current clamp series
        name = f"CurrentClampSeriesCondition{self.condition}Cell{self.cell_id}{self.position}Trial{self.trial_id}"
        current_clamp_series = CurrentClampSeries(
            name=name,
            description=f"Intracellular recording from {device_name}",
            data=voltage_volts,
            timestamps=timestamps.values,
            electrode=electrode,
        )

        # Add current clamp series to acquisition
        nwbfile.add_acquisition(current_clamp_series)
