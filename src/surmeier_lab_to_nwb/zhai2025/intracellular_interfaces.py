from datetime import datetime
from pathlib import Path
from typing import Dict, Optional

import xmltodict
from neuroconv.basedatainterface import BaseDataInterface
from neuroconv.utils import DeepDict
from pynwb.file import NWBFile
from pynwb.icephys import CurrentClampSeries, VoltageClampSeries

# Protocol step to current mapping for F-I (frequency-current) relationship experiments
# Mapping of protocol step numbers to current injection values in picoamps (pA).
#
# This protocol is used for somatic excitability experiments in Figures 1 and 3
# of Zhai et al. 2025. The protocol involves stepping current injections from
# -120 pA to +300 pA in 20 pA increments for 500 ms duration each.
#
# Key experimental details:
# - Step duration: 500 ms
# - Current range: -120 pA to +300 pA
# - Step size: 20 pA
# - Total steps: 21 (some cells may have additional steps beyond +300 pA)
# - Used for F-I relationship analysis and rheobase measurements
#
# Protocol steps:
# - 001-006: Hyperpolarizing currents (-120 to -20 pA)
# - 007-021: Depolarizing currents (+20 to +300 pA)
# - 022-026: Extended depolarizing currents (+320 to +400 pA, not all cells)
#
# References:
# - Figure 1: dSPN somatic excitability in LID model
# - Figure 3: iSPN somatic excitability and D2 receptor signaling
PROTOCOL_STEP_TO_CURRENT: Dict[str, int] = {
    "001": -120,  # pA - Strong hyperpolarization
    "002": -100,  # pA
    "003": -80,  # pA
    "004": -60,  # pA
    "005": -40,  # pA
    "006": -20,  # pA - Weak hyperpolarization
    "007": 20,  # pA - Weak depolarization (often near rheobase)
    "008": 40,  # pA
    "009": 60,  # pA
    "010": 80,  # pA
    "011": 100,  # pA
    "012": 120,  # pA
    "013": 140,  # pA
    "014": 160,  # pA
    "015": 180,  # pA
    "016": 200,  # pA
    "017": 220,  # pA
    "018": 240,  # pA
    "019": 260,  # pA
    "020": 280,  # pA
    "021": 300,  # pA - Standard maximum current
    "022": 320,  # pA - Extended protocol (some cells only)
    "023": 340,  # pA
    "024": 360,  # pA
    "025": 380,  # pA
    "026": 400,  # pA - Extended maximum current
}


class PrairieViewPathClampBaseInterface(BaseDataInterface):
    """Interface for Prairie View intracellular data."""

    keywords = ("intracellular electrophysiology", "patch clamp", "prairie view")

    def __init__(
        self,
        file_path: str | Path,
        recording_csv_file_path: Optional[str | Path] = None,
        icephys_metadata_key: str = "PrairieView",
    ):
        """
        Initialize the interface with paths to the XML and CSV files.

        Parameters
        ----------
        file_path : str | Path
            Path to the XML file containing recording metadata
        recording_csv_file_path : Optional[str | Path], default=None
            Path to the CSV file containing recording data. If None, it will be inferred from the XML file.
            Use this if the xml file is not available.
        """

        self.file_path = Path(file_path)
        self.icephys_metadata_key = icephys_metadata_key
        self._t_start = 0.0  # Default starting time

        # Load XML data
        with open(self.file_path, "r") as xml_file:
            xml_content = xml_file.read()
            self.xml_recording_dict = xmltodict.parse(xml_content)

        # Extract basic metadata to determine CSV file path
        vrec_session_entry = self.xml_recording_dict["VRecSessionEntry"]

        # Determine recording CSV file path if not provided
        if recording_csv_file_path is None:
            recording_data_file_path_stem = vrec_session_entry["DataFile"]
            self.recording_csv_file_path = self.file_path.parent / f"{recording_data_file_path_stem}.csv"
        else:
            self.recording_csv_file_path = Path(recording_csv_file_path)

        if not self.recording_csv_file_path.is_file():
            raise FileNotFoundError(f"Recording data file not found: {self.recording_csv_file_path}")

        experiment_metadata = vrec_session_entry["Experiment"]

        # Extract signal information
        signal_list = experiment_metadata["SignalList"]
        assert len(signal_list) == 1, "At the moment a single enabled signal is supported"
        signal_list = signal_list["VRecSignal"]
        enabled_signals = [signal for signal in signal_list if signal["Enabled"] == "true"]
        self.signal_metadata = enabled_signals[0]

    def get_session_start_time(self) -> Optional[datetime]:
        """
        Extract session start time from the XML file.

        Returns
        -------
        Optional[datetime]
            Session start time parsed from XML DateTime field, or None if not available
        """
        datetime_str = self.xml_recording_dict["VRecSessionEntry"].get("DateTime", None)
        if datetime_str is not None:
            try:
                # Parse ISO format datetime string
                return datetime.fromisoformat(datetime_str.replace("Z", "+00:00"))
            except (ValueError, AttributeError):
                # If parsing fails, return None
                pass
        return None

    def get_metadata(self) -> DeepDict:
        """
        Extract metadata from the XML file.

        Returns
        -------
        DeepDict
            Dictionary containing metadata
        """
        metadata = super().get_metadata()

        # Add session start time if available
        session_start_time = self.get_session_start_time()
        if session_start_time is not None:
            metadata["NWBFile"]["session_start_time"] = session_start_time

        # Extract unit metadata
        device = self.signal_metadata["Unit"]["PatchclampDevice"]
        self.device_name = device
        device_metadata = {
            "name": self.device_name,
            "description": f"Prairie View intracellular recording device: {self.device_name}",
        }
        metadata["Devices"][self.icephys_metadata_key] = device_metadata
        icephys_metadata = metadata["Icephys"]

        self.electrode_name = "IntracellularElectrode"
        icephys_metadata["IntracellularElectrodes"][self.icephys_metadata_key] = {
            "name": self.electrode_name,
            "description": "Intracellular electrodes",
        }

        return metadata

    def _get_or_create_device(self, nwbfile: NWBFile, metadata: dict):
        """
        Get or create a device for the intracellular recording.

        Parameters
        ----------
        nwbfile : NWBFile
            NWB file to add the device to
        metadata : dict
            Metadata dictionary

        Returns
        -------
        Device
            The device object
        """
        device_metadata = metadata["Devices"][self.icephys_metadata_key]
        default_device_name = "DevicePrairieViewIntracelular"
        device_name = device_metadata.get("name", default_device_name)

        # Create or get device
        if device_name in nwbfile.devices:
            device = nwbfile.devices[device_name]
        else:
            device_description = device_metadata.get(
                "description", f"Prairie View intracellular recording device: {device_name}"
            )
            device = nwbfile.create_device(name=device_name, description=device_description)

        return device

    def _get_or_create_electrode(self, nwbfile: NWBFile, metadata: dict, device):
        """
        Get or create an electrode for the intracellular recording.

        Parameters
        ----------
        nwbfile : NWBFile
            NWB file to add the electrode to
        metadata : dict
            Metadata dictionary
        device : Device
            Device object for the electrode

        Returns
        -------
        IntracellularElectrode
            The electrode object
        """
        available_electrodes = nwbfile.icephys_electrodes
        electrode_metadata = metadata["Icephys"]["IntracellularElectrodes"][self.icephys_metadata_key]

        default_electrode_name = "IntracellularElectrode"
        electrode_name = electrode_metadata.get("name", default_electrode_name)

        if electrode_name in available_electrodes:
            electrode = available_electrodes[electrode_name]
        else:
            description = electrode_metadata.get(
                "description", f"Prairie View intracellular electrode for {device.name}"
            )
            electrode = nwbfile.create_icephys_electrode(
                name=electrode_name,
                description=description,
                device=device,
            )

        return electrode

    def add_to_nwbfile(self, nwbfile: NWBFile, metadata: Optional[dict] = None):
        """
        Add device and electrode to NWBFile and return them for child classes.

        Parameters
        ----------
        nwbfile : NWBFile
            NWB file to add the device and electrode to
        metadata : dict, optional
            Metadata dictionary

        Returns
        -------
        tuple
            (device, electrode) objects for use by child classes
        """
        metadata = metadata or self.get_metadata()
        device = self._get_or_create_device(nwbfile, metadata)
        electrode = self._get_or_create_electrode(nwbfile, metadata, device)
        return device, electrode

    def set_aligned_starting_time(self, aligned_starting_time: float) -> None:
        """
        Set the aligned starting time for all time series data.

        Parameters
        ----------
        aligned_starting_time : float
            Starting time in seconds relative to session start for temporal alignment
        """
        self._t_start = aligned_starting_time


class PrairieViewCurrentClampInterface(PrairieViewPathClampBaseInterface):
    """Interface for Prairie View current clamp data."""

    keywords = ("intracellular electrophysiology", "patch clamp", "prairie view", "current clamp")

    def __init__(
        self,
        file_path: str | Path,
        recording_csv_file_path: Optional[str | Path] = None,
        icephys_metadata_key: str = "PrairieView",
    ):
        """
        Initialize the interface with paths to the XML and CSV files.

        Parameters
        ----------
        file_path : str | Path
            Path to the XML file containing recording metadata
        recording_csv_file_path : Optional[str | Path], default=None
            Path to the CSV file containing recording data. If None, it will be inferred from the XML file.
            Use this if the xml file is not available.
        """
        super().__init__(
            file_path=file_path,
            recording_csv_file_path=recording_csv_file_path,
            icephys_metadata_key=icephys_metadata_key,
        )

    def get_metadata(self) -> DeepDict:
        """
        Extract metadata from the XML file.

        Returns
        -------
        DeepDict
            Dictionary containing metadata
        """
        metadata = super().get_metadata()

        icephys_metadata = metadata["Icephys"]
        icephys_metadata["CurrentClampSeries"][self.icephys_metadata_key] = {
            "name": "CurrentClampSeries",
            "description": "Prairie View current clamp series",
            "device": self.device_name,
            "electrode": self.electrode_name,
        }

        return metadata

    @staticmethod
    def get_session_start_time_from_file(file_path: str | Path) -> Optional[datetime]:
        """
        Extract session start time from XML file without creating interface instance.

        Parameters
        ----------
        file_path : str | Path
            Path to the XML file containing recording metadata

        Returns
        -------
        Optional[datetime]
            Session start time parsed from XML DateTime field, or None if not available
        """
        from pathlib import Path

        import xmltodict

        file_path = Path(file_path)

        # Load XML data
        with open(file_path, "r") as xml_file:
            xml_content = xml_file.read()
            xml_recording_dict = xmltodict.parse(xml_content)

        # Extract session start time
        datetime_str = xml_recording_dict["VRecSessionEntry"].get("DateTime", None)
        if datetime_str is not None:
            try:
                # Parse ISO format datetime string
                return datetime.fromisoformat(datetime_str.replace("Z", "+00:00"))
            except (ValueError, AttributeError):
                # If parsing fails, return None
                pass
        return None

    def add_to_nwbfile(self, nwbfile: NWBFile, metadata: Optional[dict] = None):
        """
        Add the recording data to the NWB file.

        Parameters
        ----------
        nwbfile : NWBFile
            NWB file to add the data to
        metadata : dict
            Metadata dictionary
        """
        device, electrode = super().add_to_nwbfile(nwbfile, metadata)
        metadata = metadata or self.get_metadata()

        patch_clamp_series_metadata = metadata["Icephys"]["CurrentClampSeries"][self.icephys_metadata_key]

        import pandas as pd

        gain = float(self.signal_metadata["Gain"])
        multiplier = float(self.signal_metadata["Unit"]["Multiplier"])
        divisor = float(self.signal_metadata["Unit"]["Divisor"])
        unit = self.signal_metadata["Unit"]["UnitName"]

        # Load recording data
        recording_data_df = pd.read_csv(
            self.recording_csv_file_path,
            header=0,
        )
        timestamps = recording_data_df["Time(ms)"] / 1_000.0  # Convert to seconds
        voltage_mV = recording_data_df[" Primary"] * multiplier / divisor
        voltage_volts = voltage_mV.values / 1_000.0  # Convert mV to volts

        # Apply aligned starting time offset
        if self._t_start != 0.0:
            timestamps = timestamps + self._t_start

        # Create current clamp series
        name = patch_clamp_series_metadata["name"]
        current_clamp_series = CurrentClampSeries(
            name=name,
            description=patch_clamp_series_metadata["description"],
            data=voltage_volts,
            timestamps=timestamps.values,
            electrode=electrode,
        )

        # Add current clamp series to acquisition
        nwbfile.add_acquisition(current_clamp_series)


class PrairieViewVoltageClampInterface(PrairieViewPathClampBaseInterface):
    """Interface for Prairie View voltage clamp data."""

    keywords = ("intracellular electrophysiology", "patch clamp", "prairie view")

    def __init__(
        self,
        file_path: str | Path,
        recording_csv_file_path: Optional[str | Path] = None,
        icephys_metadata_key: str = "PrairieView",
    ):
        """
        Initialize the interface with paths to the XML and CSV files.

        Parameters
        ----------
        file_path : str | Path
            Path to the XML file containing recording metadata
        recording_csv_file_path : Optional[str | Path], default=None
            Path to the CSV file containing recording data. If None, it will be inferred from the XML file.
            Use this if the xml file is not available.
        """
        super().__init__(
            file_path=file_path,
            recording_csv_file_path=recording_csv_file_path,
            icephys_metadata_key=icephys_metadata_key,
        )

    def get_metadata(self) -> DeepDict:
        """
        Extract metadata from the XML file.

        Returns
        -------
        DeepDict
            Dictionary containing metadata
        """
        metadata = super().get_metadata()

        icephys_metadata = metadata["Icephys"]
        icephys_metadata["VoltageClampSeries"][self.icephys_metadata_key] = {
            "name": self.icephys_metadata_key,
            "description": "Prairie View voltage clamp series",
            "device": self.device_name,
            "electrode": self.electrode_name,
        }

        return metadata

    def add_to_nwbfile(self, nwbfile: NWBFile, metadata: Optional[dict] = None, starting_time: Optional[float] = None):
        """
        Add the recording data to the NWB file.

        Parameters
        ----------
        nwbfile : NWBFile
            NWB file to add the data to
        metadata : dict
            Metadata dictionary
        starting_time : float, optional
            Start time of the recording relative to session start (in seconds).
            If None, uses timestamps from data. If provided, uses starting_time
            for absolute temporal synchronization within the experimental session.
        """
        device, electrode = super().add_to_nwbfile(nwbfile, metadata)
        metadata = metadata or self.get_metadata()

        patch_clamp_series_metadata = metadata["Icephys"]["VoltageClampSeries"][self.icephys_metadata_key]

        import pandas as pd

        gain = float(self.signal_metadata["Gain"])
        multiplier = float(self.signal_metadata["Unit"]["Multiplier"])
        divisor = float(self.signal_metadata["Unit"]["Divisor"])
        unit = self.signal_metadata["Unit"]["UnitName"]

        # Load recording data
        recording_data_df = pd.read_csv(
            self.recording_csv_file_path,
            header=0,
        )
        timestamps = recording_data_df["Time(ms)"] / 1_000.0  # Convert to seconds
        current_pA = recording_data_df[" Primary"] * multiplier / divisor
        nano_amperes_factor = 1e9
        current_amps = current_pA.values / nano_amperes_factor

        name = patch_clamp_series_metadata["name"]

        # Use absolute timing if provided, otherwise use relative timestamps
        if starting_time is not None:
            # Calculate sampling rate from timestamps
            sampling_rate = 1.0 / (timestamps.iloc[1] - timestamps.iloc[0])
            voltage_clamp_series = VoltageClampSeries(
                name=name,
                description=patch_clamp_series_metadata["description"],
                data=current_amps,
                starting_time=starting_time,
                rate=sampling_rate,
                electrode=electrode,
            )
        else:
            # Apply aligned starting time offset
            if self._t_start != 0.0:
                timestamps = timestamps + self._t_start

            voltage_clamp_series = VoltageClampSeries(
                name=name,
                description=patch_clamp_series_metadata["description"],
                data=current_amps,
                timestamps=timestamps.values,
                electrode=electrode,
            )

        # Add voltage clamp series to acquisition
        nwbfile.add_acquisition(voltage_clamp_series)
