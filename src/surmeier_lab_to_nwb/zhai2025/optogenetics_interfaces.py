# -*- coding: utf-8 -*-
from pathlib import Path
from typing import Optional

import xmltodict
from ndx_optogenetics import (
    ExcitationSource,
    ExcitationSourceModel,
    OpticalFiber,
    OpticalFiberLocationsTable,
    OpticalFiberModel,
    OptogeneticEpochsTable,
    OptogeneticExperimentMetadata,
)
from neuroconv.basedatainterface import BaseDataInterface
from neuroconv.utils import DeepDict
from pynwb.file import NWBFile
from pynwb.ogen import OptogeneticStimulusSite


class PrairieViewOptogeneticsInterface(BaseDataInterface):
    """Interface for Prairie View optogenetic stimulation data from Sr2+-oEPSC experiments.

    This interface processes optogenetic stimulation data using the ndx-optogenetics extension,
    creating detailed OptogeneticEpochsTable entries with precise timing information for each
    sweep and stimulus pulse. The interface extracts LED stimulus parameters from Prairie View
    XML files and creates comprehensive metadata including excitation sources, optical fibers,
    and stimulus sites.

    Key features:
    - Extracts stimulus timing and parameters from Prairie View VoltageOutput XML files
    - Creates OptogeneticEpochsTable with one row per sweep/recording
    - Provides detailed device metadata using ndx-optogenetics structure
    - Maintains temporal synchronization with voltage clamp recordings
    """

    keywords = ("optogenetics", "led stimulation", "prairie view", "sr-oepsc")

    def __init__(
        self,
        voltage_output_xml_path: str | Path,
    ):
        """
        Initialize the optogenetics interface.

        Parameters
        ----------
        voltage_output_xml_path : str | Path
            Path to the Prairie View VoltageOutput XML configuration file.
            Expected format: '*_Cycle00001_VoltageOutput_001.xml'
        """
        self.voltage_output_xml_path = Path(voltage_output_xml_path)

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
        Extract optogenetic metadata for NWB file using ndx-optogenetics structure.

        Returns
        -------
        DeepDict
            Dictionary containing optogenetic device and stimulus site metadata.
        """
        metadata = super().get_metadata()

        # Device metadata - using ndx-optogenetics structure
        device_name = "BlueLED470nm"
        metadata["Devices"][device_name] = {
            "name": device_name,
            "description": (
                "Blue LED (470 nm) for whole-field optogenetic stimulation controlled via Prairie View. "
                "Drives 5V pulses to activate ChR2-expressing corticostriatal terminals."
            ),
        }

        # Optogenetic stimulus site metadata - maintaining compatibility with existing structure
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
        }

        # Add ndx-optogenetics metadata structure
        metadata["OptogeneticExperimentMetadata"] = {
            "stimulation_software": "Prairie View",
            "excitation_source_model": {
                "name": "LED470nmModel",
                "description": "Blue LED model for 470nm optogenetic stimulation",
                "illumination_type": "LED",
                "wavelength_range_in_nm": [455.0, 485.0],  # Typical blue LED FWHM range
            },
            "excitation_source": {
                "name": "LED470nmSource",
                "description": "Blue LED excitation source for ChR2 activation",
                "wavelength_in_nm": 470.0,
                "power_in_W": 0.001,  # Placeholder - actual power unknown
            },
            "optical_fiber_model": {
                "name": "WholeFieldFiberModel",
                "description": "Whole-field illumination setup for optogenetic stimulation",
                "fiber_name": "Whole-field LED",
                "manufacturer": "Prairie View",
                "numerical_aperture": 0.0,  # Not applicable for whole-field
                "core_diameter_in_um": 0.0,  # Not applicable for whole-field
            },
            "optical_fiber": {
                "name": "WholeFieldFiber",
                "description": "Whole-field illumination optical system",
            },
            "optical_fiber_locations": {
                "reference": "Dorsal striatum slice preparation - Sr2+-oEPSC protocol",
                "locations": [
                    {
                        "implanted_fiber_description": (
                            "Whole-field LED illumination targeting ChR2-expressing corticostriatal terminals. "
                            "Protocol: 0.3ms blue LED pulse at 20ms into 420ms sweep. "
                            "Detection window: 40-400ms post-stimulus. "
                            "Inter-sweep interval: 30 seconds for Sr2+ clearance."
                        ),
                        "location": "dorsal striatum",
                        "hemisphere": "ipsilateral to 6-OHDA lesion",
                        "ap_in_mm": 0.0,  # Not applicable for slice preparation
                        "ml_in_mm": 0.0,  # Not applicable for slice preparation
                        "dv_in_mm": 0.0,  # Not applicable for slice preparation
                    }
                ],
            },
        }

        return metadata

    def add_to_nwbfile(self, nwbfile: NWBFile, metadata: Optional[dict] = None, starting_time: Optional[float] = None):
        """
        Add optogenetic stimulation metadata and epoch data to NWB file using ndx-optogenetics structure.

        This method creates comprehensive optogenetic metadata including excitation sources,
        optical fibers, and stimulus sites. It also adds a row to the OptogeneticEpochsTable
        for this specific sweep/recording with detailed timing and stimulus parameters.

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
        ndx_metadata = metadata["OptogeneticExperimentMetadata"]

        # Create LED device (maintaining compatibility)
        device_name = "BlueLED470nm"
        if device_name not in nwbfile.devices:
            device_info = device_metadata[device_name]
            led_device = nwbfile.create_device(name=device_info["name"], description=device_info["description"])
        else:
            led_device = nwbfile.devices[device_name]

        # Create ndx-optogenetics devices
        # Excitation source model
        excitation_source_model_info = ndx_metadata["excitation_source_model"]
        excitation_source_model_name = excitation_source_model_info["name"]
        if excitation_source_model_name not in nwbfile.devices:
            excitation_source_model = ExcitationSourceModel(
                name=excitation_source_model_name,
                description=excitation_source_model_info["description"],
                illumination_type=excitation_source_model_info["illumination_type"],
                wavelength_range_in_nm=excitation_source_model_info["wavelength_range_in_nm"],
            )
            nwbfile.add_device(excitation_source_model)
        else:
            excitation_source_model = nwbfile.devices[excitation_source_model_name]

        # Excitation source
        excitation_source_info = ndx_metadata["excitation_source"]
        excitation_source_name = excitation_source_info["name"]
        if excitation_source_name not in nwbfile.devices:
            excitation_source = ExcitationSource(
                name=excitation_source_name,
                description=excitation_source_info["description"],
                wavelength_in_nm=excitation_source_info["wavelength_in_nm"],
                power_in_W=excitation_source_info["power_in_W"],
                model=excitation_source_model,
            )
            nwbfile.add_device(excitation_source)
        else:
            excitation_source = nwbfile.devices[excitation_source_name]

        # Optical fiber model
        optical_fiber_model_info = ndx_metadata["optical_fiber_model"]
        optical_fiber_model_name = optical_fiber_model_info["name"]
        if optical_fiber_model_name not in nwbfile.devices:
            optical_fiber_model = OpticalFiberModel(
                name=optical_fiber_model_name,
                description=optical_fiber_model_info["description"],
                fiber_name=optical_fiber_model_info["fiber_name"],
                manufacturer=optical_fiber_model_info["manufacturer"],
                numerical_aperture=optical_fiber_model_info["numerical_aperture"],
                core_diameter_in_um=optical_fiber_model_info["core_diameter_in_um"],
            )
            nwbfile.add_device(optical_fiber_model)
        else:
            optical_fiber_model = nwbfile.devices[optical_fiber_model_name]

        # Optical fiber
        optical_fiber_info = ndx_metadata["optical_fiber"]
        optical_fiber_name = optical_fiber_info["name"]
        if optical_fiber_name not in nwbfile.devices:
            optical_fiber = OpticalFiber(
                name=optical_fiber_name,
                description=optical_fiber_info["description"],
                model=optical_fiber_model,
            )
            nwbfile.add_device(optical_fiber)
        else:
            optical_fiber = nwbfile.devices[optical_fiber_name]

        # Optical fiber locations table
        fiber_locations_info = ndx_metadata["optical_fiber_locations"]
        optical_fiber_locations_table = OpticalFiberLocationsTable(
            name="optical_fiber_locations_table",
            description="Optical fiber locations for optogenetic stimulation",
            reference=fiber_locations_info["reference"],
        )

        # Add location data
        for location_data in fiber_locations_info["locations"]:
            optical_fiber_locations_table.add_row(
                implanted_fiber_description=location_data["implanted_fiber_description"],
                location=location_data["location"],
                hemisphere=location_data["hemisphere"],
                ap_in_mm=location_data["ap_in_mm"],
                ml_in_mm=location_data["ml_in_mm"],
                dv_in_mm=location_data["dv_in_mm"],
                excitation_source=excitation_source,
                optical_fiber=optical_fiber,
            )

        # Create OptogeneticExperimentMetadata
        if "optogenetic_experiment_metadata" not in nwbfile.lab_meta_data:
            optogenetic_experiment_metadata = OptogeneticExperimentMetadata(
                stimulation_software=ndx_metadata["stimulation_software"],
                optical_fiber_locations_table=optical_fiber_locations_table,
            )
            nwbfile.add_lab_meta_data(optogenetic_experiment_metadata)

        # Create or get OptogeneticEpochsTable
        if "optogenetic_epochs_table" not in nwbfile.intervals:
            optogenetic_epochs_table = OptogeneticEpochsTable(
                name="optogenetic_epochs_table",
                description="Optogenetic stimulation epochs with detailed parameters",
                target_tables={"optical_fiber_locations": optical_fiber_locations_table},
            )
            nwbfile.add_time_intervals(optogenetic_epochs_table)
        else:
            optogenetic_epochs_table = nwbfile.intervals["optogenetic_epochs_table"]

        # Add epoch data for this specific sweep/recording with detailed timing
        pulse_width_ms = self.stimulus_params["pulse_width_ms"]  # 0.3ms pulse
        pulse_delay_ms = self.stimulus_params["first_pulse_delay_ms"]  # 20ms delay

        # Use provided starting_time for temporal synchronization, or default to 0.0
        stimulus_starting_time = starting_time if starting_time is not None else 0.0

        # Each epoch represents one complete sweep with detailed timing structure
        # Based on Sr2+-oEPSC experimental protocol:
        # - Sweep starts at stimulus_starting_time - pulse_delay (to capture baseline)
        # - LED pulse occurs at 20ms into sweep (stimulus_starting_time)
        # - Critical detection window: 40-400ms post-stimulus (20.04-20.4s into sweep)
        # - Total sweep duration: ~0.42 seconds (420ms) to capture full analysis window

        sweep_start_time = stimulus_starting_time - (pulse_delay_ms / 1000.0)  # 20ms before pulse
        sweep_duration_s = 0.42  # 420ms total sweep duration (baseline + pulse + detection window)
        sweep_stop_time = sweep_start_time + sweep_duration_s

        optogenetic_epochs_table.add_row(
            start_time=sweep_start_time,
            stop_time=sweep_stop_time,
            stimulation_on=True,
            pulse_length_in_ms=pulse_width_ms,  # 0.3ms LED pulse
            period_in_ms=30000.0,  # 30 seconds between pulses (inter-sweep interval)
            number_pulses_per_pulse_train=1,  # Single pulse per sweep
            number_trains=1,  # Single train per sweep
            optical_fiber_locations=[0],  # Reference to first row in fiber locations table
        )

        # Create optogenetic stimulus site (maintaining compatibility)
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

        # NOTE: OptogeneticSeries removed - timing information now handled by OptogeneticEpochsTable
        # The epochs table provides more detailed and accurate stimulus timing information
        # including pulse parameters, inter-pulse intervals, and precise temporal alignment
