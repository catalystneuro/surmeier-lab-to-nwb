# -*- coding: utf-8 -*-
from pathlib import Path
from typing import Optional

import xmltodict
from neuroconv.basedatainterface import BaseDataInterface
from neuroconv.utils import DeepDict
from pynwb.file import NWBFile
from pynwb.ogen import OptogeneticSeries, OptogeneticStimulusSite


class PrairieViewOptogeneticsInterface(BaseDataInterface):
    """Interface for Prairie View optogenetic stimulation data from Sr2+-oEPSC experiments."""

    keywords = ("optogenetics", "led stimulation", "prairie view", "sr-oepsc")

    def __init__(
        self,
        voltage_output_xml_path: str | Path,
        optogenetics_metadata_key: str = "PrairieViewOptogenetics",
    ):
        """
        Initialize the optogenetics interface.

        Parameters
        ----------
        voltage_output_xml_path : str | Path
            Path to the Prairie View VoltageOutput XML configuration file.
            Expected format: '*_Cycle00001_VoltageOutput_001.xml'

        optogenetics_metadata_key : str, default="PrairieViewOptogenetics"
            Unique identifier for this optogenetics interface instance.
        """
        self.voltage_output_xml_path = Path(voltage_output_xml_path)
        self.optogenetics_metadata_key = optogenetics_metadata_key

        # Parse stimulus parameters from XML
        self.stimulus_params = self._parse_optogenetic_stimulus()

    def _parse_optogenetic_stimulus(self) -> dict:
        """
        Extract optogenetic LED stimulation parameters from Prairie View XML.

        Returns
        -------
        dict
            Dictionary containing stimulus parameters:
            - 'pulse_count' (int): Number of LED pulses in the stimulus train
            - 'pulse_width_ms' (float): Duration of each LED pulse in milliseconds
            - 'pulse_amplitude_v' (float): LED drive voltage amplitude in volts
            - 'first_pulse_delay_ms' (float): Delay before first pulse in milliseconds
            - 'rest_potential_v' (float): Baseline voltage when LED is off
            - 'ao_line' (str): Analog output line number controlling the LED
        """
        with open(self.voltage_output_xml_path, "r") as xml_file:
            xml_content = xml_file.read()
            voltage_output_dict = xmltodict.parse(xml_content)

        # Navigate to LED signal parameters - this structure has Experiment as root
        if "Experiment" in voltage_output_dict:
            waveforms = voltage_output_dict["Experiment"]["Waveform"]

            # Find LED waveform (assuming it's named "LED" and enabled)
            led_waveform = None
            if isinstance(waveforms, list):
                for waveform in waveforms:
                    if waveform.get("Name") == "LED" and waveform.get("Enabled") == "true":
                        led_waveform = waveform
                        break
            elif isinstance(waveforms, dict) and waveforms.get("Name") == "LED":
                led_waveform = waveforms

            if led_waveform is None:
                raise ValueError(f"Could not find enabled LED waveform in {self.voltage_output_xml_path}")

            # Extract pulse train parameters
            pulse_train = led_waveform["WaveformComponent_PulseTrain"]
        else:
            raise ValueError(f"Unexpected XML structure in {self.voltage_output_xml_path}")

        return {
            "pulse_count": int(pulse_train["PulseCount"]),
            "pulse_width_ms": float(pulse_train["PulseWidth"]),
            "pulse_amplitude_v": float(pulse_train["PulsePotentialStart"]),
            "first_pulse_delay_ms": float(pulse_train["FirstPulseDelay"]),
            "rest_potential_v": float(pulse_train["RestPotential"]),
            "ao_line": led_waveform["AOLine"],
        }

    def get_metadata(self) -> DeepDict:
        """
        Extract optogenetic metadata for NWB file.

        Returns
        -------
        DeepDict
            Dictionary containing optogenetic device and stimulus site metadata.
        """
        metadata = super().get_metadata()

        # Device metadata
        device_name = "BlueLED_470nm"
        metadata["Devices"][device_name] = {
            "name": device_name,
            "description": (
                "Blue LED (470 nm) for whole-field optogenetic stimulation controlled via Prairie View. "
                "Drives 5V pulses to activate ChR2-expressing corticostriatal terminals."
            ),
        }

        # Optogenetic stimulus site metadata
        site_name = "corticostriatal_terminals"
        metadata["Optogenetics"] = {
            "OptogeneticStimulusSites": {
                site_name: {
                    "name": site_name,
                    "device": device_name,
                    "description": (
                        "Optogenetic stimulation site in dorsal striatum targeting ChR2-expressing corticostriatal terminals. "
                        "ChR2 expression from motor cortex injection ipsilateral to 6-OHDA lesion. "
                        "LED illumination activates projections to evoke Sr2+-oEPSCs in recorded dSPNs."
                    ),
                    "excitation_lambda": 470.0,  # nm, blue LED
                    "location": "dorsal striatum (ipsilateral to 6-OHDA lesion)",
                }
            },
            "OptogeneticSeries": {
                self.optogenetics_metadata_key: {
                    "name": self.optogenetics_metadata_key,
                    "description": (
                        f"LED stimulus: {self.stimulus_params['pulse_width_ms']} ms blue LED pulse "
                        f"delivered at {self.stimulus_params['first_pulse_delay_ms']} ms delay. "
                        f"Control voltage: {self.stimulus_params['pulse_amplitude_v']} V. "
                        f"NOTE: Actual optical power unknown - only electrical control parameters available. "
                        f"Pulses occur every 30 seconds to evoke asynchronous EPSCs."
                    ),
                    "site": site_name,
                }
            },
        }

        return metadata

    def add_to_nwbfile(self, nwbfile: NWBFile, metadata: Optional[dict] = None, starting_time: Optional[float] = None):
        """
        Add optogenetic stimulation data and metadata to NWB file.

        Parameters
        ----------
        nwbfile : NWBFile
            The target NWB file object where optogenetic data will be stored.

        metadata : dict, optional
            Metadata dictionary. If not provided, will use get_metadata().

        starting_time : float, optional
            Start time of the optogenetic stimulus relative to session start (in seconds).
            This should include the sweep offset + pulse delay (typically 20ms) for
            proper temporal synchronization with the voltage clamp recording.
        """
        metadata = metadata or self.get_metadata()

        optogenetics_metadata = metadata["Optogenetics"]
        device_metadata = metadata["Devices"]

        # Create LED device
        device_name = "BlueLED_470nm"
        if device_name not in nwbfile.devices:
            device_info = device_metadata[device_name]
            led_device = nwbfile.create_device(name=device_info["name"], description=device_info["description"])
        else:
            led_device = nwbfile.devices[device_name]

        # Create optogenetic stimulus site
        site_name = "corticostriatal_terminals"
        if site_name not in nwbfile.ogen_sites:
            site_info = optogenetics_metadata["OptogeneticStimulusSites"][site_name]
            ogen_site = OptogeneticStimulusSite(
                name=site_info["name"],
                device=led_device,
                description=site_info["description"],
                excitation_lambda=site_info["excitation_lambda"],
                location=site_info["location"],
            )
            nwbfile.add_ogen_site(ogen_site)
        else:
            ogen_site = nwbfile.ogen_sites[site_name]

        # Create stimulus time series data
        pulse_width_s = self.stimulus_params["pulse_width_ms"] / 1000.0
        pulse_delay_s = self.stimulus_params["first_pulse_delay_ms"] / 1000.0

        # NOTE: Optical power not available in source data
        # The XML files only contain LED control voltage (5V), not actual optical power
        # For NWB compliance, we use a placeholder value and document the limitation
        optical_power_unknown = 0.001  # Placeholder: 1 mW (typical optogenetics range)

        # Create simple pulse stimulus data
        # Assuming typical recording duration of ~1 second
        sampling_rate = 10000  # 10 kHz
        duration_s = 1.0
        n_samples = int(duration_s * sampling_rate)

        stimulus_data = [0.0] * n_samples
        pulse_start_sample = int(pulse_delay_s * sampling_rate)
        pulse_end_sample = int((pulse_delay_s + pulse_width_s) * sampling_rate)

        # Set optical power during pulse period (placeholder value)
        # Actual optical power unknown - only control voltage ({control_voltage}V) available
        for i in range(pulse_start_sample, min(pulse_end_sample, n_samples)):
            stimulus_data[i] = optical_power_unknown

        # Create OptogeneticSeries with proper timing
        series_info = optogenetics_metadata["OptogeneticSeries"][self.optogenetics_metadata_key]

        # Use provided starting_time for temporal synchronization, or default to 0.0
        stimulus_starting_time = starting_time if starting_time is not None else 0.0

        stimulus_series = OptogeneticSeries(
            name=series_info["name"],
            description=series_info["description"],
            data=stimulus_data,
            site=ogen_site,
            starting_time=stimulus_starting_time,
            rate=float(sampling_rate),
        )

        # Add stimulus series to NWB file
        nwbfile.add_stimulus(stimulus_series)
