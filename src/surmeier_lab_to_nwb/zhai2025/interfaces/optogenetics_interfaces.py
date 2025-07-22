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
    OptogeneticVirus,
    OptogeneticViruses,
    OptogeneticVirusInjection,
    OptogeneticVirusInjections,
)
from neuroconv.basedatainterface import BaseDataInterface
from neuroconv.utils import DeepDict
from pynwb.file import NWBFile
from pynwb.ogen import OptogeneticStimulusSite


class PrairieViewOptogeneticsInterface(BaseDataInterface):
    """Interface for Prairie View optogenetic stimulation data from Sr2+-oEPSC experiments.

    This interface processes optogenetic stimulation data using the ndx-optogenetics extension,
    creating comprehensive metadata including virus information, injection coordinates, device
    specifications, and detailed timing information for each sweep phase. The interface extracts
    LED stimulus parameters from Prairie View XML files and creates a complete representation
    of the optogenetic experimental setup.

    Key features:
    - Extracts stimulus timing and parameters from Prairie View VoltageOutput XML files
    - Creates OptogeneticEpochsTable with four rows per sweep/recording:
      1. Baseline Period (0-20ms): Pre-stimulus baseline measurement
      2. Stimulation Period (20-20.3ms): 0.3ms LED pulse activation
      3. Detection Window (40-400ms post-stimulus): Critical EPSC analysis period
      4. Inter-Sweep Interval (420ms-30s): Recovery period between sweeps
    - Creates comprehensive virus metadata (AAV5-hSyn-hChR2(H134R)-EYFP, Addgene #26973)
    - Records stereotaxic injection coordinates (M1 motor cortex: AP +1.15, ML -1.60, DV -1.55mm)
    - Provides detailed device metadata (LED sources, stimulation sites)
    - Maintains temporal synchronization with voltage clamp recordings
    - Follows ndx-optogenetics best practices for structured optogenetic metadata
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
            "virus": {
                "name": "DefaultOptogeneticVirus",
                "construct_name": "Default optogenetic virus construct",
                "description": "Default optogenetic virus - should be overridden in script metadata",
                "manufacturer": "Unknown",
                "titer_in_vg_per_ml": 1e12,  # Default titer
            },
            "virus_injection": {
                "name": "DefaultVirusInjection",
                "description": "Default virus injection - should be overridden in script metadata",
                "location": "Unknown",
                "hemisphere": "Unknown",
                "reference": "Bregma at the cortical surface",
                "ap_in_mm": 0.0,
                "ml_in_mm": 0.0,
                "dv_in_mm": 0.0,
                "pitch_in_deg": 0.0,
                "yaw_in_deg": 0.0,
                "roll_in_deg": 0.0,
                "volume_in_uL": 0.0,
                "injection_date": None,
            },
            "optical_fiber_model": {
                "name": "PlaceholderFiberModel",
                "description": "Placeholder fiber model - no actual optical fiber used in this experiment",
                "fiber_name": "Whole-field LED (no fiber)",
                "manufacturer": "Not applicable",
                "numerical_aperture": 0.0,  # Not applicable for whole-field LED
                "core_diameter_in_um": 0.0,  # Not applicable for whole-field LED
            },
            "optical_fiber": {
                "name": "PlaceholderFiber",
                "description": "Placeholder fiber device - no actual optical fiber used, LED provides whole-field illumination",
            },
            "stimulation_location": {
                "reference": "Dorsal striatum slice preparation - Sr2+-oEPSC protocol",
                "location": "dorsal striatum",
                "hemisphere": "ipsilateral to 6-OHDA lesion",
                "description": (
                    "Whole-field LED illumination targeting ChR2-expressing corticostriatal terminals. "
                    "Protocol: 0.3ms blue LED pulse at 20ms into 420ms sweep. "
                    "Detection window: 40-400ms post-stimulus. "
                    "Inter-sweep interval: 30 seconds for Sr2+ clearance. "
                    "No optical fiber used - LED provides broad illumination of brain slice."
                ),
                "ap_in_mm": 0.0,  # Not applicable for slice preparation
                "ml_in_mm": 0.0,  # Not applicable for slice preparation
                "dv_in_mm": 0.0,  # Not applicable for slice preparation
            },
        }

        return metadata

    def add_to_nwbfile(
        self,
        nwbfile: NWBFile,
        metadata: Optional[dict] = None,
        starting_time: Optional[float] = None,
        recording_duration: Optional[float] = None,
        next_recording_start_time: Optional[float] = None,
    ):
        """
        Add optogenetic stimulation metadata and epoch data to NWB file using ndx-optogenetics structure.

        This method creates comprehensive optogenetic metadata including excitation sources
        and stimulus sites. It also adds a row to the OptogeneticEpochsTable
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

        recording_duration : float, optional
            Duration of the voltage clamp recording in seconds. If not provided,
            defaults to 0.42 seconds (420ms) which is typical for Sr2+-oEPSC recordings.

        next_recording_start_time : float, optional
            Start time of the next recording in the session (in seconds). If provided,
            the inter_sweep_interval stage will end at this time. If None, no
            inter_sweep_interval stage will be added (for the last recording).
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

        # Create placeholder optical fiber model (required by ndx-optogenetics API)
        # NOTE: No actual optical fiber is used - this is a placeholder to satisfy API requirements
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

        # Create placeholder optical fiber (required by ndx-optogenetics API)
        # NOTE: No actual optical fiber is used - this is a placeholder to satisfy API requirements
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

        # Create OpticalFiberLocationsTable for excitation source tracking
        # NOTE: This table is required by ndx-optogenetics but no actual optical fibers are used
        # The table tracks the excitation source location for whole-field LED illumination
        stimulation_location_info = ndx_metadata["stimulation_location"]
        optical_fiber_locations_table = OpticalFiberLocationsTable(
            name="optical_fiber_locations_table",
            description="Stimulation location tracking for whole-field LED illumination (no optical fibers used)",
            reference=stimulation_location_info["reference"],
        )

        # Add stimulation location data (excitation source and placeholder fiber)
        optical_fiber_locations_table.add_row(
            implanted_fiber_description=stimulation_location_info["description"],
            location=stimulation_location_info["location"],
            hemisphere=stimulation_location_info["hemisphere"],
            ap_in_mm=stimulation_location_info["ap_in_mm"],
            ml_in_mm=stimulation_location_info["ml_in_mm"],
            dv_in_mm=stimulation_location_info["dv_in_mm"],
            excitation_source=excitation_source,
            optical_fiber=optical_fiber,  # Placeholder fiber to satisfy API requirements
        )

        # Create OptogeneticVirus and OptogeneticVirusInjection objects
        virus_info = ndx_metadata["virus"]
        injection_info = ndx_metadata["virus_injection"]

        # Create OptogeneticVirus object with comprehensive description
        enhanced_virus_description = (
            f"{virus_info['description']} "
            f"Construct: {virus_info['construct_name']} from {virus_info['manufacturer']}. "
            f"Titer: {virus_info['titer_in_vg_per_ml']:.0e} vg/ml."
        )

        optogenetic_virus = OptogeneticVirus(
            name=virus_info["name"],
            construct_name=virus_info["construct_name"],
            description=enhanced_virus_description,
            manufacturer=virus_info["manufacturer"],
            titer_in_vg_per_ml=virus_info["titer_in_vg_per_ml"],
        )

        # Create OptogeneticVirusInjection object with comprehensive description
        enhanced_injection_description = (
            f"{injection_info['description']} "
            f"Coordinates: AP {injection_info['ap_in_mm']}mm, ML {injection_info['ml_in_mm']}mm, "
            f"DV {injection_info['dv_in_mm']}mm relative to {injection_info['reference']}. "
            f"Injection volume: {injection_info['volume_in_uL']}µL."
        )

        optogenetic_virus_injection = OptogeneticVirusInjection(
            name=injection_info["name"],
            description=enhanced_injection_description,
            location=injection_info["location"],
            hemisphere=injection_info["hemisphere"],
            reference=injection_info["reference"],
            ap_in_mm=injection_info["ap_in_mm"],
            ml_in_mm=injection_info["ml_in_mm"],
            dv_in_mm=injection_info["dv_in_mm"],
            pitch_in_deg=injection_info["pitch_in_deg"],
            yaw_in_deg=injection_info["yaw_in_deg"],
            roll_in_deg=injection_info["roll_in_deg"],
            volume_in_uL=injection_info["volume_in_uL"],
            injection_date=injection_info["injection_date"],
            virus=optogenetic_virus,
        )

        # Create OptogeneticViruses and OptogeneticVirusInjections containers
        optogenetic_viruses = OptogeneticViruses(optogenetic_virus=[optogenetic_virus])
        optogenetic_virus_injections = OptogeneticVirusInjections(
            optogenetic_virus_injections=[optogenetic_virus_injection]
        )

        # Create OptogeneticExperimentMetadata
        if "optogenetic_experiment_metadata" not in nwbfile.lab_meta_data:
            optogenetic_experiment_metadata = OptogeneticExperimentMetadata(
                stimulation_software=ndx_metadata["stimulation_software"],
                optical_fiber_locations_table=optical_fiber_locations_table,
                optogenetic_viruses=optogenetic_viruses,
                optogenetic_virus_injections=optogenetic_virus_injections,
            )
            nwbfile.add_lab_meta_data(optogenetic_experiment_metadata)

        # Create or get OptogeneticEpochsTable
        if "optogenetic_epochs_table" not in nwbfile.intervals:
            optogenetic_epochs_table = OptogeneticEpochsTable(
                name="optogenetic_epochs_table",
                description="Optogenetic stimulation epochs with detailed parameters",
                target_tables={"optical_fiber_locations": optical_fiber_locations_table},
            )

            # Add custom column for stage name
            optogenetic_epochs_table.add_column(
                name="stage_name",
                description="Name of the experimental stage: pre_stimulation, stimulation, detection, post_detection, or inter_sweep_interval",
            )

            nwbfile.add_time_intervals(optogenetic_epochs_table)
        else:
            optogenetic_epochs_table = nwbfile.intervals["optogenetic_epochs_table"]

        # Add epoch rows for this specific sweep/recording with detailed timing
        # Number of rows depends on whether next_recording_start_time is provided
        pulse_width_ms = self.stimulus_params["pulse_width_ms"]  # 0.3ms pulse
        pulse_delay_ms = self.stimulus_params["first_pulse_delay_ms"]  # 20ms delay

        # Use provided starting_time for temporal synchronization, or default to 0.0
        stimulus_starting_time = starting_time if starting_time is not None else 0.0

        # Use provided recording_duration, or default to 0.42 seconds (420ms)
        voltage_clamp_duration = recording_duration if recording_duration is not None else 0.42

        # Calculate absolute times for each phase based on Sr2+-oEPSC experimental protocol
        sweep_start_time = stimulus_starting_time - (pulse_delay_ms / 1000.0)  # 20ms before pulse
        stimulus_start_time = stimulus_starting_time  # LED pulse start
        stimulus_end_time = stimulus_starting_time + (pulse_width_ms / 1000.0)  # LED pulse end
        detection_start_time = stimulus_starting_time + 0.040  # 40ms post-stimulus
        detection_end_time = stimulus_starting_time + 0.400  # 400ms post-stimulus
        voltage_clamp_end_time = sweep_start_time + voltage_clamp_duration  # End of voltage clamp recording
        inter_sweep_start_time = voltage_clamp_end_time  # Start of inter-sweep interval

        # Ensure post_detection phase has duration if voltage clamp extends beyond detection
        if voltage_clamp_end_time <= detection_end_time:
            # If voltage clamp ends before or at detection end, extend it slightly
            voltage_clamp_end_time = detection_end_time + 0.020  # Add 20ms buffer
            inter_sweep_start_time = voltage_clamp_end_time

        # Phase 1: Baseline Period (0-20ms)
        optogenetic_epochs_table.add_row(
            start_time=sweep_start_time,
            stop_time=stimulus_start_time,
            stimulation_on=False,  # No stimulation during baseline
            pulse_length_in_ms=float("nan"),  # NaN when stimulation is off
            period_in_ms=float("nan"),  # NaN when stimulation is off
            number_pulses_per_pulse_train=-1,  # -1 when stimulation is off
            number_trains=-1,  # -1 when stimulation is off
            intertrain_interval_in_ms=float("nan"),  # NaN when stimulation is off
            power_in_mW=0.0,  # No power during baseline
            optical_fiber_locations=[0],  # References excitation source location
            stage_name="pre_stimulation",
        )

        # Phase 2: Stimulation Period (20-20.3ms)
        optogenetic_epochs_table.add_row(
            start_time=stimulus_start_time,
            stop_time=stimulus_end_time,
            stimulation_on=True,
            pulse_length_in_ms=pulse_width_ms,  # 0.3ms LED pulse
            period_in_ms=float("nan"),  # NaN because only single pulse per train
            number_pulses_per_pulse_train=1,  # Single pulse per train
            number_trains=1,  # Single train per epoch
            intertrain_interval_in_ms=float("nan"),  # Not applicable for single train
            power_in_mW=float("nan"),  # Power not measured/known
            optical_fiber_locations=[0],  # References excitation source location
            stage_name="stimulation",
        )

        # Phase 3: Detection Window (40-400ms post-stimulus)
        optogenetic_epochs_table.add_row(
            start_time=detection_start_time,
            stop_time=detection_end_time,
            stimulation_on=False,  # No stimulation during detection
            pulse_length_in_ms=float("nan"),  # NaN when stimulation is off
            period_in_ms=float("nan"),  # NaN when stimulation is off
            number_pulses_per_pulse_train=-1,  # -1 when stimulation is off
            number_trains=-1,  # -1 when stimulation is off
            intertrain_interval_in_ms=float("nan"),  # NaN when stimulation is off
            power_in_mW=0.0,  # No power during detection
            optical_fiber_locations=[0],  # References excitation source location
            stage_name="detection",
        )

        # Phase 4: Post-Detection Period (400ms - end of voltage clamp)
        optogenetic_epochs_table.add_row(
            start_time=detection_end_time,
            stop_time=voltage_clamp_end_time,
            stimulation_on=False,  # No stimulation during post-detection
            pulse_length_in_ms=float("nan"),  # NaN when stimulation is off
            period_in_ms=float("nan"),  # NaN when stimulation is off
            number_pulses_per_pulse_train=-1,  # -1 when stimulation is off
            number_trains=-1,  # -1 when stimulation is off
            intertrain_interval_in_ms=float("nan"),  # NaN when stimulation is off
            power_in_mW=0.0,  # No power during post-detection
            optical_fiber_locations=[0],  # References excitation source location
            stage_name="post_detection",
        )

        # Phase 5: Inter-Sweep Interval (end of voltage clamp - next sweep)
        # Only add this phase if next_recording_start_time is provided (not the last recording)
        if next_recording_start_time is not None:
            optogenetic_epochs_table.add_row(
                start_time=inter_sweep_start_time,
                stop_time=next_recording_start_time,  # Use actual next recording start time
                stimulation_on=False,  # No stimulation during inter-sweep
                pulse_length_in_ms=float("nan"),  # NaN when stimulation is off
                period_in_ms=float("nan"),  # NaN when stimulation is off
                number_pulses_per_pulse_train=-1,  # -1 when stimulation is off
                number_trains=-1,  # -1 when stimulation is off
                intertrain_interval_in_ms=float("nan"),  # NaN when stimulation is off
                power_in_mW=0.0,  # No power during inter-sweep
                optical_fiber_locations=[0],  # References excitation source location
                stage_name="inter_sweep_interval",
            )

        # Create optogenetic stimulus site (maintaining compatibility)
        # This is used to add location data as that is done in the fibers and this
        # experiment does not have fibers
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

        # Create OptogeneticSeries to satisfy NWB requirements that stimulus sites be referenced
        # While timing information is primarily in OptogeneticEpochsTable, we need a basic
        # OptogeneticSeries to link the stimulus site and provide stimulus trace data
        import numpy as np
        from pynwb.ogen import OptogeneticSeries

        # Create a simple stimulus trace based on the stimulus timing
        # Generate stimulus data: 0 during baseline, 1 during stimulation, 0 during detection
        stimulus_duration = stimulus_end_time - stimulus_start_time  # Should be 0.3ms
        total_duration = voltage_clamp_end_time - sweep_start_time

        # Create timeline with high temporal resolution (1 kHz sampling for stimulus precision)
        sampling_rate = 1000.0  # Hz
        num_samples = int(total_duration * sampling_rate)
        stimulus_data = np.zeros(num_samples)

        # Set stimulus data to 1 during stimulation period
        stimulus_start_sample = int((stimulus_start_time - sweep_start_time) * sampling_rate)
        stimulus_end_sample = int((stimulus_end_time - sweep_start_time) * sampling_rate)
        stimulus_data[stimulus_start_sample:stimulus_end_sample] = 1.0

        # Create OptogeneticSeries with sweep naming convention to match VoltageClampSeries
        # Count existing OptogeneticSeries to determine sweep number
        existing_optogenetic_series = [
            name for name in nwbfile.stimulus.keys() if name.startswith("OptogeneticSeriesSweep")
        ]
        sweep_number = len(existing_optogenetic_series) + 1
        optogenetic_series_name = f"OptogeneticSeriesSweep{sweep_number:03d}"

        optogenetic_series = OptogeneticSeries(
            name=optogenetic_series_name,
            description=f"LED stimulation trace - 0.3ms pulse during voltage clamp sweep",
            data=stimulus_data,
            site=ogen_site,
            starting_time=sweep_start_time,
            rate=sampling_rate,
        )
        nwbfile.add_stimulus(optogenetic_series)
