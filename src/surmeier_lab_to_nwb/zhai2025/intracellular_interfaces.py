from datetime import datetime
from pathlib import Path
from typing import Optional

import xmltodict
from neuroconv.basedatainterface import BaseDataInterface
from neuroconv.utils import DeepDict
from pynwb.file import NWBFile
from pynwb.icephys import CurrentClampSeries, VoltageClampSeries


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

    def get_metadata(self) -> DeepDict:
        """
        Extract metadata from the XML file.

        Returns
        -------
        DeepDict
            Dictionary containing metadata
        """
        metadata = super().get_metadata()

        # Extract session start time
        session_start_time = None
        datetime_str = self.xml_recording_dict["VRecSessionEntry"].get("DateTime", None)
        if datetime_str is not None:
            try:
                # Parse ISO format datetime string
                session_start_time = datetime.fromisoformat(datetime_str.replace("Z", "+00:00"))
            except (ValueError, AttributeError):
                # If parsing fails, leave as None
                pass

        # Add session start time if available
        if session_start_time is not None:
            metadata["NWBFile"]["session_start_time"] = session_start_time

        # Extract unit metadata
        device = self.signal_metadata["Unit"]["PatchclampDevice"]
        self.device_name = device
        device_metadata = {
            "name": self.device_name,
            "description": f"Prairie View intracellular recording device: {self.device_name}",
        }
        metadata["Devices"][self.device_name] = device_metadata
        icephys_metadata = metadata["Icephys"]

        self.electrode_name = "IntracellularElectrode"
        icephys_metadata["IntracellularElectrodes"][self.electrode_name] = {
            "name": self.electrode_name,
            "description": "Intracellular electrodes",
        }

        return metadata


class PrairieViewCurrentClampInterface(PrairieViewPathClampBaseInterface):
    """Interface for Prairie View current clamp data."""

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
        icephys_metadata["CurrentClampSeries"][self.icephys_metadata_key] = {
            "name": "CurrentClampSeries",
            "description": "Prairie View current clamp series",
            "device": self.device_name,
            "electrode": self.electrode_name,
        }

        return metadata

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

        metadata = metadata or self.get_metadata()

        patch_clamp_series_metadata = metadata["Icephys"]["CurrentClampSeries"][self.icephys_metadata_key]
        device_tag = patch_clamp_series_metadata["device"]
        device_metadata = metadata["Devices"][device_tag]

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

        # Create or use provided electrode
        available_electrodes = nwbfile.icephys_electrodes
        electrode_tag = patch_clamp_series_metadata.get("electrode", "IntracellularElectrode")

        electrode_metadata = metadata["Icephys"]["IntracellularElectrodes"][electrode_tag]

        default_electrode_name = "IntracellularElectrode"
        electrode_name = electrode_metadata.get("name", default_electrode_name)

        if electrode_name in available_electrodes:
            electrode = available_electrodes[electrode_name]
        else:
            description = electrode_metadata.get(
                "description", f"Prairie View intracellular electrode for {device_name}"
            )
            electrode = nwbfile.create_icephys_electrode(
                name=electrode_name,
                description=description,
                device=device,
            )

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
            "name": "VoltageClampSeries",
            "description": "Prairie View voltage clamp series",
            "device": self.device_name,
            "electrode": self.electrode_name,
        }

        return metadata

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

        metadata = metadata or self.get_metadata()

        patch_clamp_series_metadata = metadata["Icephys"]["VoltageClampSeries"][self.icephys_metadata_key]
        device_tag = patch_clamp_series_metadata["device"]
        device_metadata = metadata["Devices"][device_tag]

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

        # Create or use provided electrode
        available_electrodes = nwbfile.icephys_electrodes
        electrode_tag = patch_clamp_series_metadata.get("electrode", "IntracellularElectrode")

        electrode_metadata = metadata["Icephys"]["IntracellularElectrodes"][electrode_tag]

        default_electrode_name = "IntracellularElectrode"
        electrode_name = electrode_metadata.get("name", default_electrode_name)

        if electrode_name in available_electrodes:
            electrode = available_electrodes[electrode_name]
        else:
            description = electrode_metadata.get(
                "description", f"Prairie View intracellular electrode for {device_name}"
            )
            electrode = nwbfile.create_icephys_electrode(
                name=electrode_name,
                description=description,
                device=device,
            )

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
        voltage_clamp_series = VoltageClampSeries(
            name=name,
            description=patch_clamp_series_metadata["description"],
            data=current_amps,
            timestamps=timestamps.values,
            electrode=electrode,
        )

        # Add voltage clamp series to acquisition
        nwbfile.add_acquisition(voltage_clamp_series)
