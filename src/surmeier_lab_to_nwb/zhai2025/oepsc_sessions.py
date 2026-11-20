"""Sr2+-oEPSC sessions (Figures 2 and 4): one NWB file per mouse and recording day.

Whole-cell voltage clamp of striatal spiny projection neurons while a blue LED flash activates ChR2-expressing motor
cortex axons, in Ca2+-free solution with Sr2+, so the evoked synaptic currents are asynchronous.
"""

from __future__ import annotations

import uuid
import xml.etree.ElementTree as ET
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

from ndx_ophys_devices import (
    Effector,
    ExcitationSource,
    ExcitationSourceModel,
    ViralVector,
    ViralVectorInjection,
)
from ndx_optogenetics import (
    OptogeneticEffectors,
    OptogeneticExperimentMetadata,
    OptogeneticPulsesTable,
    OptogeneticSitesTable,
    OptogeneticViruses,
    OptogeneticVirusInjections,
)
from neuroconv.tools.nwb_helpers import (
    configure_and_write_nwbfile,
    make_nwbfile_from_metadata,
)
from pynwb import NWBFile
from pynwb.device import Device

from surmeier_lab_to_nwb.zhai2025.metadata import (
    device_metadata,
    session_metadata,
    shared_metadata,
)
from surmeier_lab_to_nwb.zhai2025.mouse_day_merger import MouseDayBundle, SweepRecord
from surmeier_lab_to_nwb.zhai2025.patch_clamp_sessions import (
    _CONDITION_LABEL,
    TIMEZONE,
    _add_sweeps_through_neuroconv,
    _build_electrode,
    _build_session_id,
    _find_voltage_output_xml,
    _label_protocol_runs,
    _patch_clamp_subject_description,
)


@dataclass
class LightEvokedProtocol:
    """The commands of one Sr2+-oEPSC sweep, from its VoltageOutput XML."""

    step_delay_ms: float  # amplifier command step
    step_width_ms: float
    step_amplitude_volts: float
    flash_delay_ms: float  # LED pulse
    flash_width_ms: float


def _single_pulse(xml_path: Path, waveform_name: str) -> tuple[ET.Element, ET.Element]:
    """The enabled waveform of this name and its only component, which must be a single pulse."""
    root = ET.parse(xml_path).getroot()
    waveforms = [
        w for w in root.findall("Waveform") if w.findtext("Name") == waveform_name and w.findtext("Enabled") == "true"
    ]
    if len(waveforms) != 1:
        raise ValueError(f"Expected one enabled {waveform_name} waveform in {xml_path}")
    components = [c for c in waveforms[0] if c.tag.startswith("WaveformComponent")]
    if len(components) != 1 or components[0].tag != "WaveformComponent_PulseTrain":
        raise ValueError(f"Expected the {waveform_name} waveform in {xml_path} to be one pulse train")
    if components[0].findtext("PulseCount") != "1" or components[0].findtext("Repetitions") != "1":
        raise ValueError(f"Expected the {waveform_name} waveform in {xml_path} to be a single pulse")
    return waveforms[0], components[0]


def _read_light_evoked_protocol(sweep_path: Path) -> LightEvokedProtocol:
    """Read the amplifier step and the LED pulse of one Sr2+-oEPSC sweep."""
    xml_path = _find_voltage_output_xml(sweep_path)
    if xml_path is None:
        raise FileNotFoundError(f"No VoltageOutput XML in {sweep_path}")
    amplifier, step = _single_pulse(xml_path, "700B")
    if amplifier.findtext("Units") != "mV":
        raise ValueError(f"Expected the 700B command of {xml_path} in mV")
    # The analog output in volts times the amplifier's command sensitivity (mV per V) gives the command in mV.
    step_mV = float(step.findtext("PulsePotentialStart")) * float(amplifier.findtext("UnitScaleFactor"))
    _led, flash = _single_pulse(xml_path, "LED")
    return LightEvokedProtocol(
        step_delay_ms=float(step.findtext("FirstPulseDelay")),
        step_width_ms=float(step.findtext("PulseWidth")),
        step_amplitude_volts=step_mV * 1e-3,
        flash_delay_ms=float(flash.findtext("FirstPulseDelay")),
        flash_width_ms=float(flash.findtext("PulseWidth")),
    )


def _add_optogenetics_metadata(nwbfile: NWBFile, optogenetics: dict) -> OptogeneticSitesTable:
    """The LED, the ChR2 virus and its injection, from the paper's Methods; returns the sites table."""
    led_model = ExcitationSourceModel(**optogenetics["ExcitationSourceModel"])
    nwbfile.add_device_model(led_model)
    led = ExcitationSource(**optogenetics["ExcitationSource"], model=led_model)
    nwbfile.add_device(led)
    virus = ViralVector(**optogenetics["ViralVector"])
    injection = ViralVectorInjection(**optogenetics["ViralVectorInjection"], viral_vector=virus)
    effector = Effector(**optogenetics["Effector"], viral_vector_injection=injection)
    sites_table = OptogeneticSitesTable(**optogenetics["OptogeneticSitesTable"])
    sites_table.add_row(excitation_source=led, effector=effector)
    nwbfile.add_lab_meta_data(
        OptogeneticExperimentMetadata(
            stimulation_software=optogenetics["stimulation_software"],
            optogenetic_sites_table=sites_table,
            optogenetic_viruses=OptogeneticViruses(viral_vectors=[virus]),
            optogenetic_virus_injections=OptogeneticVirusInjections(viral_vector_injections=[injection]),
            optogenetic_effectors=OptogeneticEffectors(effectors=[effector]),
        )
    )
    return sites_table


def merge_oepsc_mouse_day(bundle: MouseDayBundle, output_path: Path) -> Path:
    """Build one Sr2+-oEPSC NWB file for a mouse-day (Figures 2 and 4) and write it to disk."""
    session = session_metadata("oepsc")
    subject_entry = bundle.subject_entry
    conditions = {_CONDITION_LABEL[sw.condition_subfolder] for sweeps in bundle.cells.values() for sw in sweeps}
    if len(conditions) != 1:
        raise ValueError(f"{subject_entry.subject_id} has Sr2+-oEPSC recordings under {sorted(conditions)}")
    (condition,) = conditions
    session_start_time = bundle.session_start_time.replace(tzinfo=TIMEZONE)
    analysis_window = session["analysis_window"]

    metadata = shared_metadata()
    shared = metadata["NWBFile"]
    metadata["NWBFile"].update(
        session_description=session["NWBFile"]["session_description"].format(
            figure=session["figure_by_line"][subject_entry.animal],
            n_cells=len(bundle.cell_indices()),
            cells=session["cells_by_line"][subject_entry.animal],
            condition=condition,
        ),
        experiment_description=session["NWBFile"]["experiment_description"],
        identifier=str(uuid.uuid4()),
        session_start_time=session_start_time,
        session_id=_build_session_id(bundle),
        surgery=session["NWBFile"]["surgery"].format(surgery=shared["surgery"]),
        pharmacology=session["NWBFile"]["pharmacology"].format(pharmacology=shared["pharmacology"]),
    )
    metadata["Subject"].update(
        subject_id=subject_entry.subject_id,
        strain=subject_entry.strain,
        genotype=subject_entry.genotype,
        description=_patch_clamp_subject_description(bundle, virus=session["Subject"]["virus"]),
    )
    # The subjects table gives no strain or genotype for some lines, which the file then leaves out
    metadata["Subject"] = {key: value for key, value in metadata["Subject"].items() if value is not None}
    nwbfile = make_nwbfile_from_metadata(metadata)
    amplifier = Device(**device_metadata("MultiClamp700B"))
    nwbfile.add_device(amplifier)
    optogenetics = session["Optogenetics"]
    sites_table = _add_optogenetics_metadata(nwbfile, optogenetics)

    pulses = optogenetics["OptogeneticPulsesTable"]
    pulses_table = OptogeneticPulsesTable(
        name=pulses["name"],
        description=pulses["description"].format(analysis_window=analysis_window),
        target_tables={"optogenetic_sites": sites_table},
    )
    for name, description in pulses["columns"].items():
        pulses_table.add_column(name=name, description=description)
    nwbfile.add_time_intervals(pulses_table)

    icephys = session["Icephys"]
    irt = nwbfile.get_intracellular_recordings()
    irt_values: dict[int, dict] = {}
    sweep_per_row: dict[int, SweepRecord] = {}
    irt_rows_per_cell: dict[int, list[int]] = defaultdict(list)
    for cell_index in bundle.cell_indices():
        electrode = _build_electrode(nwbfile, cell_index, amplifier)
        sweeps = bundle.sweeps_for(cell_index)
        series, stimulus = icephys["VoltageClampSeries"], icephys["VoltageClampStimulusSeries"]
        rows_and_starts = _add_sweeps_through_neuroconv(
            nwbfile,
            sweeps=sweeps,
            electrode=electrode,
            name=series["name"].format(cell=cell_index),
            stimulus_name=stimulus["name"].format(cell=cell_index),
            description=series["description"].format(
                cell=cell_index, analysis_window=analysis_window, n_sweeps=len(sweeps)
            ),
            stimulus_description=stimulus["description"].format(cell=cell_index),
            session_start_time=session_start_time,
        )
        for sw, (row_index, sweep_start_time) in zip(sweeps, rows_and_starts, strict=True):
            irt_values[row_index] = {
                "cell_id": cell_index,
                "condition": condition,
                "site": "Soma",
                "sweep_start_time": sweep_start_time,
            }
            sweep_per_row[row_index] = sw
            irt_rows_per_cell[cell_index].append(row_index)
            protocol = _read_light_evoked_protocol(sw.recording_path)
            flash_start = sweep_start_time + protocol.flash_delay_ms * 1e-3
            pulses_table.add_row(
                start_time=flash_start,
                stop_time=flash_start + protocol.flash_width_ms * 1e-3,
                power_in_mW=float("nan"),
                wavelength_in_nm=470.0,
                optogenetic_sites=[0],
                cell_id=cell_index,
                led_setting=sw.led_setting,
            )

    _label_protocol_runs(irt, irt_rows_per_cell.values(), sweep_per_row)
    for name, description in icephys["intracellular_recordings_columns"].items():
        irt.add_column(name=name, description=description, data=[irt_values[row][name] for row in range(len(irt))])

    # One sequential recording and one repetition per cell, all under the mouse-day's one condition.
    repetitions_table = nwbfile.get_icephys_repetitions()
    for name, description in icephys["repetitions_columns"].items():
        repetitions_table.add_column(name=name, description=description)
    repetition_indices = []
    for cell_index, irt_rows in irt_rows_per_cell.items():
        simultaneous = [nwbfile.add_icephys_simultaneous_recording(recordings=[row]) for row in irt_rows]
        sequential = nwbfile.add_icephys_sequential_recording(
            simultaneous_recordings=simultaneous, stimulus_type=icephys["stimulus_type"]
        )
        repetition_indices.append(
            nwbfile.add_icephys_repetition(
                sequential_recordings=[sequential],
                cell_id=cell_index,
                protocol=icephys["stimulus_type"],
                condition=condition,
            )
        )
    conditions_table = nwbfile.get_icephys_experimental_conditions()
    for name, description in icephys["experimental_conditions_columns"].items():
        conditions_table.add_column(name=name, description=description)
    nwbfile.add_icephys_experimental_condition(repetitions=repetition_indices, condition=condition)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    configure_and_write_nwbfile(nwbfile=nwbfile, nwbfile_path=output_path, backend="hdf5", metadata=shared_metadata())
    return output_path
