# -*- coding: utf-8 -*-
from pathlib import Path
from typing import Any, Dict, Optional

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
        sweep_info: Dict[str, Any],
        optogenetics_metadata_key: str = "PrairieViewOptogenetics",
    ):
        """
        Initialize the optogenetics interface.

        Parameters
        ----------
        voltage_output_xml_path : str | Path
            Path to the Prairie View VoltageOutput XML configuration file.
            Expected format: '*_Cycle00001_VoltageOutput_001.xml'

        sweep_info : Dict[str, Any]
            Sweep metadata dictionary containing:
            - 'cell_number': Identifier of the recorded cell
            - 'led_number': LED stimulation protocol identifier
            - 'sweep_number': Sequential sweep number
            - 'sweep_name': Complete sweep folder name

        optogenetics_metadata_key : str, default="PrairieViewOptogenetics"
            Unique identifier for this optogenetics interface instance.
        """
        self.voltage_output_xml_path = Path(voltage_output_xml_path)
        self.sweep_info = sweep_info
        self.optogenetics_metadata_key = optogenetics_metadata_key

        # Parse stimulus parameters from XML
        self.stimulus_params = self._parse_optogenetic_stimulus()

    def _parse_optogenetic_stimulus(self) -> Dict[str, Any]:
        """
        Extract optogenetic LED stimulation parameters from Prairie View XML.

        Returns
        -------
        Dict[str, Any]
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
                "Blue LED (470 nm wavelength) for whole-field optogenetic stimulation of "
                "AAV5-hSyn-hChR2(H134R)-EYFP expressing corticostriatal terminals. LED controlled "
                "via Prairie View analog output line 2 with 5V drive voltage for 0.3 ms pulses "
                "delivered every 30 seconds to evoke Sr2+-oEPSCs in dSPNs."
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
                        "Optogenetic stimulation site targeting AAV5-hSyn-hChR2(H134R)-EYFP expressing "
                        "corticostriatal terminals in dorsal striatum. ChR2 expression achieved via "
                        "stereotaxic injection into motor cortex ipsilateral to 6-OHDA lesion. "
                        "Whole-field LED illumination activates corticostriatal projections to evoke "
                        "Sr2+-oEPSCs in dSPNs under voltage clamp at -70 mV holding potential."
                    ),
                    "excitation_lambda": 470.0,  # nm, blue LED
                    "location": "dorsal striatum (ipsilateral to 6-OHDA lesion)",
                }
            },
            "OptogeneticSeries": {
                self.optogenetics_metadata_key: {
                    "name": f"OptogeneticStimulus_Cell{self.sweep_info['cell_number']}_LED{self.sweep_info['led_number']}_Sweep{self.sweep_info['sweep_number']}",
                    "description": (
                        f"Blue LED optogenetic stimulus (470 nm, {self.stimulus_params['pulse_width_ms']} ms, "
                        f"{self.stimulus_params['pulse_amplitude_v']} V) for Sr2+-oEPSC experiment. "
                        f"Cell {self.sweep_info['cell_number']}, LED protocol {self.sweep_info['led_number']}, "
                        f"Sweep {self.sweep_info['sweep_number']}. Stimulus delivered at "
                        f"{self.stimulus_params['first_pulse_delay_ms']} ms to activate ChR2-expressing "
                        f"corticostriatal terminals and evoke asynchronous EPSCs in Ca2+-free ACSF "
                        f"with 3 mM Sr2+ substitution. LED pulses occur every 30 seconds during the sweep."
                    ),
                    "site": site_name,
                }
            },
        }

        return metadata

    def add_to_nwbfile(self, nwbfile: NWBFile, metadata: Optional[dict] = None):
        """
        Add optogenetic stimulation data and metadata to NWB file.

        Parameters
        ----------
        nwbfile : NWBFile
            The target NWB file object where optogenetic data will be stored.

        metadata : dict, optional
            Metadata dictionary. If not provided, will use get_metadata().
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
        pulse_amplitude = self.stimulus_params["pulse_amplitude_v"]

        # Create simple pulse stimulus data
        # Assuming typical recording duration of ~1 second
        sampling_rate = 10000  # 10 kHz
        duration_s = 1.0
        n_samples = int(duration_s * sampling_rate)

        stimulus_data = [0.0] * n_samples
        pulse_start_sample = int(pulse_delay_s * sampling_rate)
        pulse_end_sample = int((pulse_delay_s + pulse_width_s) * sampling_rate)

        # Set pulse amplitude during pulse period
        for i in range(pulse_start_sample, min(pulse_end_sample, n_samples)):
            stimulus_data[i] = pulse_amplitude

        # Create OptogeneticSeries
        series_info = optogenetics_metadata["OptogeneticSeries"][self.optogenetics_metadata_key]
        stimulus_series = OptogeneticSeries(
            name=series_info["name"],
            description=series_info["description"],
            data=stimulus_data,
            site=ogen_site,
            starting_time=0.0,
            rate=float(sampling_rate),
        )

        # Add stimulus series to NWB file
        nwbfile.add_stimulus(stimulus_series)
