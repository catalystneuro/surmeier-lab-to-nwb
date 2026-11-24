"""Patch-clamp sessions (Figures 1, 3, 6, 7 and 9): one NWB file per mouse and recording day.

Each file holds every cell recorded from the mouse that day: one CurrentClampSeries per cell and measurement (somatic
current steps, dendritic trials) written through NeuroConv's Bruker interface, the two-photon line scans of the dendritic
trials, the intracellular recordings tables and a RecordingsIndexTable that gathers each sweep's rows. The Sr2+-oEPSC and
GRAB-ACh sessions, in `oepsc_sessions` and `grab_ach_sessions`, share the helpers here.
"""

from __future__ import annotations

import logging
import re
import uuid
import warnings
import xml.etree.ElementTree as ET
from collections import defaultdict
from dataclasses import dataclass
from datetime import datetime
from decimal import Decimal
from pathlib import Path
from zoneinfo import ZoneInfo

# Bruker Prairie View leaves UIC CreateTime/LastSavedTime as julianday=0 in the
# line-scan TIFFs. tifffile warns when parsing these unparsable fields, but the
# pixel data reads correctly and the real acquisition timestamp lives in the
# PVScan @date attribute of the master XML (which we read separately). Suppress
# the noise.
logging.getLogger("tifffile").setLevel(logging.ERROR)
warnings.filterwarnings("ignore", message=".*no datetime before year 1.*")
warnings.filterwarnings("ignore", message="invalid value encountered in divide", module="tifffile")

import numpy as np
import pandas as pd
from hdmf.common import DynamicTable, DynamicTableRegion, VectorData
from neuroconv.datainterfaces import BrukerVoltageRecordingInterface
from neuroconv.tools.nwb_helpers import (
    configure_and_write_nwbfile,
    make_nwbfile_from_metadata,
)
from pynwb import NWBFile
from pynwb.base import TimeSeries, TimeSeriesReference, TimeSeriesReferenceVectorData
from pynwb.device import Device
from pynwb.icephys import IntracellularElectrode

from surmeier_lab_to_nwb.zhai2025.behavior_sessions import (
    _CRISPR_INJECTION,
    _CRISPR_VIRUSES,
)
from surmeier_lab_to_nwb.zhai2025.devices import get_or_create_bruker_ultima_model
from surmeier_lab_to_nwb.zhai2025.interfaces.bruker_line_scan_interface import (
    BrukerLineScanInterface,
)
from surmeier_lab_to_nwb.zhai2025.metadata import (
    device_metadata,
    session_metadata,
    shared_metadata,
)
from surmeier_lab_to_nwb.zhai2025.mouse_day_merger import (
    MouseDayBundle,
    SweepRecord,
)

TIMEZONE = ZoneInfo("America/Chicago")  # Northwestern University


# --------- line scans ----------------------------------------------------------------------


_PATCH_CLAMP = session_metadata("patch_clamp")
_LINE_SCANS = _PATCH_CLAMP["LineScans"]
# The dye of each detector channel, which names the objects in the file, and its description
_DYE_NAME = {channel: text["name"] for channel, text in _LINE_SCANS["channels"].items()}
_DYE_DESCRIPTION = {channel: text["description"] for channel, text in _LINE_SCANS["channels"].items()}
_DENDRITIC_TRIAL = " " + _LINE_SCANS["dendritic_trial"]


def _spot_sentence(spot: str) -> str:
    """Description sentence for a site the lab imaged at two spots of the dendrite; empty for one spot."""
    if not spot:
        return ""
    return " " + _LINE_SCANS["spot"].format(spot=spot[-1], real=_LINE_SCANS["spot_real"] if spot == "Spot2" else ".")


def _add_line_scan(
    nwbfile: NWBFile,
    *,
    sweeps: list[SweepRecord],
    sweep_start_times: list[float],
    channel: str,
    name: str,
    roi_name: str,
    about: str,
    kymograph_note: str,
    roi_response_note: str,
    images: list[tuple[str, str]],
) -> tuple[TimeSeries, TimeSeries, dict[str, tuple[int, int]]]:
    """Write one channel of the line scans of sweeps that share a scan line, each sweep at its start time.

    The series are named with `name` after their type (`LineScan{name}`, `RoiResponseSeries{name}`), and the ROI is
    the row `roi_name` of the dye's `PlaneSegmentation{dye}`, which every line scan of that dye shares. `about` says
    what was scanned and the notes end the descriptions of the two series. `images` gives the name and description
    of each sweep's source image, in the order of `sweeps`. Returns the line scan, the fluorescence and the first
    sample and number of samples of each sweep in them, by folder name.
    """
    interface = BrukerLineScanInterface(
        folder_paths=[sweep.recording_path for sweep in sweeps],
        channel_name=channel,
        roi_id=roi_name,
        metadata_key=name,
    )
    for sweep, sweep_start_time in zip(sweeps, sweep_start_times):
        interface.temporal_alignment[sweep.recording_path.name].move_start_to(sweep_start_time)
    trial = interface.segmentation_extractor.get_trials()[0]
    start_x, start_y, stop_x, stop_y = trial.line
    x0, x1 = trial.profile_pixels
    dye = _DYE_NAME[channel]

    metadata = interface.get_metadata()
    metadata["Devices"]["bruker_device"]["name"] = device_metadata("BrukerUltima")["name"]  # added by _ensure_devices
    channel_text = _LINE_SCANS["channels"][channel]
    dye_description = channel_text["description"]
    ophys = metadata["Ophys"]
    plane = _LINE_SCANS["ImagingPlane"]
    ophys["ImagingPlanes"][name].update(
        name=plane["name"].format(dye=dye),
        description=plane["description"].format(dye_description=dye_description),
        excitation_lambda=plane["excitation_lambda"],
        indicator=channel_text["indicator"],
        location=plane["location"],
        optical_channel=[
            dict(
                name=_LINE_SCANS["OpticalChannel"]["name"].format(dye=dye),
                description=channel_text["detector"],
                emission_lambda=np.nan,
            )
        ],
    )
    segmentation = _LINE_SCANS["PlaneSegmentation"]
    ophys["SegmentationTables"] = {
        f"LineScan{dye}": dict(
            name=segmentation["name"].format(dye=dye),
            description=segmentation["description"].format(dye_description=dye_description),
        )
    }
    ophys["PlaneSegmentations"][name]["segmentation_table_metadata_key"] = f"LineScan{dye}"
    roi_response = _LINE_SCANS["RoiResponseSeries"]
    ophys["RoiResponses"][name]["raw"].update(
        name=roi_response["name"].format(name=name),
        description=roi_response["description"].format(about=about, x0=x0, x1=x1, note=roi_response_note),
    )
    kymograph = _LINE_SCANS["LineScanKymograph"]
    ophys["LineScanKymographs"][name].update(
        name=kymograph["name"].format(name=name),
        description=kymograph["description"].format(
            about=about, start_x=start_x, start_y=start_y, stop_x=stop_x, stop_y=stop_y, note=kymograph_note
        ),
    )
    ophys["LineScanSourceImages"][name] = {
        sweep.recording_path.name: dict(name=image_name, description=image_description)
        for sweep, (image_name, image_description) in zip(sweeps, images)
    }
    interface.add_to_nwbfile(nwbfile=nwbfile, metadata=metadata)
    return (
        nwbfile.acquisition[f"LineScan{name}"],
        nwbfile.processing["ophys"]["Fluorescence"][f"RoiResponseSeries{name}"],
        interface.get_trial_sample_ranges(),
    )


# --------- reading per-sweep stimulus protocol ---------------------------------------------


@dataclass
class StimulusProtocol:
    """Current-clamp command of one sweep, read from its VoltageOutput XML."""

    name: str  # e.g., "Step 5_-40 pA", "Trio bAP"
    amplitude_pA: float  # amplitude of the injected current in picoamperes: the step, or each dendritic-trial pulse
    pulses_ms: list[tuple[float, float]]  # (start, width) of each current pulse, in ms from the start of the sweep
    duration_ms: float  # length of the whole command, its components back to back


def _parse_voltage_output(xml_path: Path) -> StimulusProtocol:
    """Read the current-clamp command from a Bruker VoltageOutput XML.

    The enabled '700B' waveform (the amplifier command) is a sequence of pulse-train components played back to back:
    each waits its FirstPulseDelay, then gives PulseCount pulses of PulseWidth separated by PulseSpacing. Somatic
    sweeps have their step in the first component (200 ms in, 500 ms long). Dendritic sweeps have three 2 ms pulses
    in the second ("Theta burst"), 900 ms after a 100 ms first component, so at 1000, 1020 and 1040 ms, when the
    recorded spikes occur. The other components are 0 pA.
    """
    root = ET.parse(xml_path).getroot()
    waveforms = [w for w in root.findall("Waveform") if w.findtext("Enabled", default="false").lower() == "true"]
    if [w.findtext("Name") for w in waveforms] != ["700B"] or waveforms[0].findtext("Units") != "pA":
        raise ValueError(f"{xml_path} does not have the 700B current command as its only output")
    scale_factor = Decimal(waveforms[0].findtext("UnitScaleFactor"))
    amplitudes: set[float] = set()
    pulses_ms: list[tuple[float, float]] = []
    component_start = 0.0
    for component in waveforms[0]:
        if not component.tag.startswith("WaveformComponent"):
            continue
        values = {child.tag: float(child.text) for child in component if child.tag != "Name"}
        if (
            component.tag != "WaveformComponent_PulseTrain"
            or values["Repetitions"] != 1
            or values["PulsePotentialDelta"] != 0
            or values["RestPotential"] != 0
        ):
            raise ValueError(f"{xml_path}: component {component.findtext('Name')} is not one train of equal pulses")
        count, width, spacing = int(values["PulseCount"]), values["PulseWidth"], values["PulseSpacing"]
        first_start = component_start + values["FirstPulseDelay"]
        # In decimal, so 0.55 V at 400 pA/V is 220 pA and not 220.00000000000003
        amplitude = float(Decimal(component.findtext("PulsePotentialStart")) * scale_factor)
        if amplitude != 0:
            amplitudes.add(amplitude)
            pulses_ms.extend((first_start + number * (width + spacing), width) for number in range(count))
        component_start = first_start + count * width + (count - 1) * spacing
    if len(amplitudes) != 1:
        raise ValueError(f"{xml_path} injects current at {sorted(amplitudes)} pA, not at one amplitude")
    return StimulusProtocol(
        name=root.findtext("Name", default="").strip(),
        amplitude_pA=amplitudes.pop(),
        pulses_ms=pulses_ms,
        duration_ms=component_start,
    )


def _find_voltage_output_xml(sweep_path: Path) -> Path | None:
    for candidate in sweep_path.glob("*_VoltageOutput_001.xml"):
        return candidate
    return None


def _voltage_recording_csv(sweep: SweepRecord) -> Path:
    """The sweep's VoltageRecording CSV."""
    csv_paths = list(sweep.recording_path.glob("*_VoltageRecording_001.csv"))
    if len(csv_paths) != 1:
        raise FileNotFoundError(f"Expected one VoltageRecording CSV in {sweep.recording_path}")
    return csv_paths[0]


def _acquisition_stem(sweep: SweepRecord) -> str:
    """The name Prairie View gave the sweep's acquisition, which NeuroConv's Bruker interface labels a protocol run by."""
    return re.sub(r"_Cycle\d+_VoltageRecording_\d+$", "", _voltage_recording_csv(sweep).stem)


def _recorded_seconds(sweep: SweepRecord) -> float:
    """How long the sweep's voltage recording lasts, from the time column of its CSV."""
    time_ms = pd.read_csv(_voltage_recording_csv(sweep), usecols=[0]).iloc[:, 0].to_numpy()
    return len(time_ms) / round(1000.0 / (time_ms[1] - time_ms[0]))


def _voltage_recording_interface(sweeps: list[SweepRecord], **options) -> BrukerVoltageRecordingInterface:
    """NeuroConv's Bruker interface over the sweeps' amplifier recordings (the Primary signal).

    It refuses a sweep without its master XML beside it, whose start would come from its VoltageRecording DateTime,
    which can be seconds early.
    """
    with warnings.catch_warnings():
        warnings.filterwarnings("error", message=".*have no PVScan XML beside them")
        return BrukerVoltageRecordingInterface(
            file_paths=[_voltage_recording_csv(sweep) for sweep in sweeps], response_signal_name="Primary", **options
        )


def _recording_start(sweep: SweepRecord, session_start_time: datetime) -> float:
    """Start of the sweep's voltage recording in seconds from `session_start_time`: the master XML's Sequence time
    plus the recording's offset in it, as NeuroConv's Bruker interface reads it."""
    interface = _voltage_recording_interface([sweep])
    return (interface.get_metadata()["NWBFile"]["session_start_time"] - session_start_time).total_seconds()


def _sample_time(series: TimeSeries, index: int) -> float:
    if series.timestamps is not None:
        return float(series.timestamps[index])
    return series.starting_time + index / series.rate


def _add_sweeps_through_neuroconv(
    nwbfile: NWBFile,
    *,
    sweeps: list[SweepRecord],
    electrode: IntracellularElectrode,
    name: str,
    stimulus_name: str,
    description: str,
    stimulus_description: str,
    session_start_time: datetime,
) -> list[tuple[int, float]]:
    """Write the sweeps' recordings as one series, and their amplifier commands as its stimulus, with one
    intracellular recordings row per sweep, through NeuroConv's Bruker interface.

    The interface places each sweep at its start (the master XML's Sequence time plus the recording's offset in it)
    and draws each sweep's command from its VoltageOutput XML over the command's own span. It writes each row's
    protocol_run as the acquisition stem of the series' first sweep; the caller relabels the rows by protocol run.
    Returns each sweep's row and start in seconds from `session_start_time`, in the order of `sweeps`.
    """
    interface = _voltage_recording_interface(sweeps, stimulus_command="700B", metadata_key=name)
    metadata = interface.get_metadata()
    interface.temporal_alignment.shift_times(
        (metadata["NWBFile"]["session_start_time"] - session_start_time).total_seconds()
    )
    series_metadata = metadata["Icephys"]["PatchClampSeries"][name]
    # The cell's electrode is in the file already; the interface reuses it by name
    metadata["Icephys"]["IntracellularElectrodes"][series_metadata["electrode_metadata_key"]]["name"] = electrode.name
    series_metadata.update(name=name, description=description)
    metadata["Icephys"]["PatchClampStimulusSeries"][name].update(name=stimulus_name, description=stimulus_description)

    recordings = nwbfile.get_intracellular_recordings()
    first_row = len(recordings)
    interface.add_to_nwbfile(nwbfile, metadata=metadata)
    rows = range(first_row, len(recordings))
    if len(rows) != len(sweeps):
        raise ValueError(f"NeuroConv wrote {len(rows)} rows for the {len(sweeps)} sweeps of {name}")
    series = nwbfile.acquisition[name]
    responses = recordings.category_tables["responses"]["response"]
    session_start_date = session_start_time.replace(tzinfo=None)
    rows_and_starts = []
    for sweep, row in zip(sweeps, rows):
        start = _sample_time(series, responses.data[row].idx_start)
        # The PVScan date is the Sequence time truncated to the second, and the recording starts a few milliseconds
        # into the Sequence, which ties each row to its sweep
        if not 0 <= start - (sweep.pv_scan_date - session_start_date).total_seconds() < 1.1:
            raise ValueError(f"Row {row} of {name} starts at {start} s, outside the second of {sweep.recording_path}")
        rows_and_starts.append((row, start))
    return rows_and_starts


# --------- core merger logic ---------------------------------------------------


def _series_name(cell_index: int, modality: str, neurodata_type: str = "CurrentClampSeries") -> str:
    suffix = "Somatic" if modality == "soma_excitability" else "Dendritic"
    return f"{neurodata_type}Cell{cell_index}{suffix}"


_PUBLISHED_FIGURE = _PATCH_CLAMP["published_figure"]
_CELLS_OF_TYPE = {
    cell_type: (text["singular"], text["plural"]) for cell_type, text in _PATCH_CLAMP["cells_of_type"].items()
}


def _patch_clamp_session_description(bundle: MouseDayBundle, cell_type: str) -> str:
    """What was recorded from the mouse that day: figure, cells, conditions in recording order, and the two assays."""
    singular, plural = _CELLS_OF_TYPE[cell_type]
    n_cells = len(bundle.cell_indices())
    conditions: list[str] = []
    for sweep in sorted((sw for sweeps in bundle.cells.values() for sw in sweeps), key=lambda sw: sw.pv_scan_date):
        condition = _sweep_condition(sweep)
        if condition not in conditions:
            conditions.append(condition)
    drugs = {_STATE_AND_DRUG[condition][1] for condition in conditions}
    somatic = [cell for cell in bundle.cell_indices() if bundle.sweeps_for(cell, "soma_excitability")]
    dendritic = [cell for cell in bundle.cell_indices() if bundle.sweeps_for(cell, "dend_excitability")]

    text = _PATCH_CLAMP["NWBFile"]["session_description"]
    parts = [
        text["recordings"].format(
            figure=_PUBLISHED_FIGURE[bundle.figure_code],
            n_cells=n_cells,
            cells=singular if n_cells == 1 else plural,
            conditions_label="Condition" if len(conditions) == 1 else "Conditions",
            conditions="; ".join(conditions),
        )
    ]
    if somatic:
        parts.append(text["somatic"].format(n_cells=len(somatic), cells="cell" if len(somatic) == 1 else "cells"))
    if dendritic:
        parts.append(text["dendritic"].format(n_cells=len(dendritic), cells="cell" if len(dendritic) == 1 else "cells"))
    if "none" in drugs and drugs & {"SCH 23390", "sulpiride", "oxo-M"}:
        parts.append(text["added_drug"].format(drug=" and ".join(sorted(drugs - {"none"}))))
    parts.append(text["tables"])
    return " ".join(parts)


def _protocol_description(protocol_label: str) -> str:
    """Readable protocol for the repetitions table from the internal label ('soma_fi_curve', 'dend_prox1')."""
    protocols = _PATCH_CLAMP["Icephys"]["protocol"]
    if protocol_label == "soma_fi_curve":
        return protocols["soma_fi_curve"]
    site = protocol_label.removeprefix("dend_")
    dendrite_type = "Proximal" if site.startswith("prox") else "Distal"
    return protocols["dendritic"].format(site=f"{dendrite_type}{site[4:]}")


# The lab's condition folders named as in the paper's figure legends, each condition's state of the mouse and drug in
# the bath, and the text the drugs added during the recordings add to the pharmacology
_CONDITION_LABEL = _PATCH_CLAMP["condition_by_folder"]
_STATE_AND_DRUG = {
    condition: (text["state"], text["drug"]) for condition, text in _PATCH_CLAMP["state_and_drug"].items()
}
_ADDED_DRUG_PHARMACOLOGY = {
    drug: " " + text for drug, text in _PATCH_CLAMP["NWBFile"]["added_drug_pharmacology"].items()
}


def _sweep_condition(sweep: SweepRecord) -> str:
    """The sweep's condition, named as in the paper's figure legends."""
    condition = _CONDITION_LABEL[sweep.condition_subfolder]
    if sweep.suffix == "oxoM":
        condition += " with oxo-M"
    return condition


_NAIVE_SURGERY = _PATCH_CLAMP["NWBFile"]["naive_surgery"]
_OXO_M_PHARMACOLOGY = _PATCH_CLAMP["NWBFile"]["oxo_m_pharmacology"]
_M1R_ANTAGONIST_PHARMACOLOGY = {
    category: " " + text for category, text in _PATCH_CLAMP["NWBFile"]["m1r_antagonist_pharmacology"].items()
}
_CRISPR_REPORTER = " " + _PATCH_CLAMP["NWBFile"]["crispr_reporter"]
_LINE_DESCRIPTION = _PATCH_CLAMP["Subject"]["line_description"]
_GROUP_DESCRIPTION = _PATCH_CLAMP["Subject"]["group_description"]
_CONTROL_DESCRIPTION = _PATCH_CLAMP["Subject"]["control_description"]


def _patch_clamp_subject_description(bundle: MouseDayBundle, virus: str = "") -> str:
    """The mouse line and what was done to the mouse, for the patch-clamp and Sr2+-oEPSC files."""
    subject_entry = bundle.subject_entry
    if subject_entry.category == "CONT":
        group = _CONTROL_DESCRIPTION[bundle.figure_code]
    else:
        group = _GROUP_DESCRIPTION[subject_entry.category]
    return " ".join(
        sentence
        for sentence in [
            _LINE_DESCRIPTION[subject_entry.animal],
            group,
            virus,
            _PATCH_CLAMP["Subject"]["identification"],
        ]
        if sentence
    )


def _ensure_devices(nwbfile: NWBFile) -> tuple[Device, Device]:
    """Add Bruker microscope and Multiclamp amplifier devices."""
    bruker_model = get_or_create_bruker_ultima_model(nwbfile)
    bruker_name = device_metadata("BrukerUltima")["name"]
    if bruker_name in nwbfile.devices:
        bruker = nwbfile.devices[bruker_name]
    else:
        bruker = Device(**device_metadata("BrukerUltima"), model=bruker_model)
        nwbfile.add_device(bruker)

    amp_name = device_metadata("MultiClamp700B")["name"]
    if amp_name in nwbfile.devices:
        amp = nwbfile.devices[amp_name]
    else:
        amp = Device(**device_metadata("MultiClamp700B"))
        nwbfile.add_device(amp)
    return bruker, amp


def _build_electrode(nwbfile: NWBFile, cell_index: int, device: Device) -> IntracellularElectrode:
    """The cell's whole-cell patch pipette, at the cell body; the dendritic sites were only imaged."""
    electrode = _PATCH_CLAMP["Icephys"]["IntracellularElectrode"]
    name = electrode["name"].format(cell=cell_index)
    if name in nwbfile.icephys_electrodes:
        return nwbfile.icephys_electrodes[name]
    return nwbfile.create_icephys_electrode(
        name=name,
        description=electrode["description"].format(cell=cell_index),
        device=device,
        location=electrode["location"],
        cell_id=electrode["cell_id"].format(cell=cell_index),
    )


def _add_patch_clamp_series(
    nwbfile: NWBFile,
    sweeps: list[SweepRecord],
    cell_index: int,
    modality: str,
    electrode: IntracellularElectrode,
    session_start_time: datetime,
) -> list[tuple[SweepRecord, int, float, StimulusProtocol]]:
    """Write all sweeps of (cell, modality) as one CurrentClampSeries and its stimulus through NeuroConv's Bruker
    interface, which adds one intracellular recordings row per sweep.

    Returns, per sweep in recording order, the sweep, its row, its start in seconds from the session start and its
    stimulus protocol, parsed from the VoltageOutput XML beside it for the stimulus_current_pA and
    stimulus_protocol_name columns.
    """
    protocols: list[StimulusProtocol] = []
    # Durations of the sweeps whose recording stops before their current command ends.
    short_sweep_seconds: list[float] = []
    for sw in sweeps:
        vo_xml = _find_voltage_output_xml(sw.recording_path)
        if vo_xml is None:
            raise FileNotFoundError(f"No VoltageOutput file with the current command in {sw.recording_path}")
        protocol = _parse_voltage_output(vo_xml)
        recorded_seconds = _recorded_seconds(sw)
        if sum(protocol.pulses_ms[-1]) > recorded_seconds * 1e3:
            if modality != "dend_excitability":
                raise ValueError(f"{sw.recording_path} is recorded for less time than its current command lasts")
            short_sweep_seconds.append(recorded_seconds)
        protocols.append(protocol)

    series_name = _series_name(cell_index, modality)
    series_text = _PATCH_CLAMP["Icephys"]["CurrentClampSeries"]
    recording = series_text["recording"][modality].format(dendritic_trial=_DENDRITIC_TRIAL)
    short_sweeps = ""
    if short_sweep_seconds:
        seconds = " or ".join(f"{duration:g}" for duration in sorted(set(short_sweep_seconds)))
        short_sweeps = " " + series_text["short_sweeps"].format(n_short=len(short_sweep_seconds), seconds=seconds)
    numbering = "".join(f" {note}" for note in sorted({sw.numbering_note for sw in sweeps if sw.numbering_note}))
    rows_and_starts = _add_sweeps_through_neuroconv(
        nwbfile,
        sweeps=sweeps,
        electrode=electrode,
        name=series_name,
        stimulus_name=_series_name(cell_index, modality, "CurrentClampStimulusSeries"),
        description=series_text["description"].format(
            cell=cell_index, recording=recording, n_sweeps=len(sweeps), short_sweeps=short_sweeps, numbering=numbering
        ),
        stimulus_description=_PATCH_CLAMP["Icephys"]["CurrentClampStimulusSeries"]["description"].format(
            series=series_name, short_sweeps=short_sweeps
        ),
        session_start_time=session_start_time,
    )
    return [
        (sw, row, start, protocol)
        for sw, (row, start), protocol in zip(sweeps, rows_and_starts, protocols, strict=True)
    ]


def _build_session_id(bundle: MouseDayBundle) -> str:
    """Build the 5-token session_id from the mouse-day bundle."""
    # Cell-type from animal reporter
    cell_type_map = {
        "tD-Tomato": "dSPN",  # Drd1-tdTomato
        "eGFP": "iSPN",  # Drd2-eGFP
        "CDGIko X eGFP": "iSPN",
        "Adora2-Cre": "iSPN",
    }
    cell_type = cell_type_map.get(bundle.subject_entry.animal, "pan")
    # State from category
    state_map = {
        "LID OFF": "OffState",
        "LID ON": "OnState",
        "CONT": "UnlesionedControl",
        "6-OHDA": "LesionedControl",
        "M1R ant": "M1RAntag",
        "M1R CRSPR": "M1RCRISPR",
        "NAIVE": "Naive",
    }
    state = state_map.get(bundle.subject_entry.category, "Unknown")
    all_conditions = {s.condition_subfolder for sweeps in bundle.cells.values() for s in sweeps}
    # The Data Connections sheet files the vehicle controls of Figure 6 and the gRNA-only controls of Figure 8 as
    # CONT, but they are lesioned, levodopa-treated mice recorded in the off-state.
    if all_conditions & {"control", "interleaved control"}:
        state = "OffState"
    # Pharm from the most distinguishing condition in the bundle. If a paired-
    # recording cell exists, surface the drug used (SCH or sul).
    pharm = "none"
    if any("SCH" in c for c in all_conditions):
        pharm = "D1RaSch"
    elif any("sul" in c for c in all_conditions):
        pharm = "D2RaSul"
    elif any(s.suffix == "oxoM" for sweeps in bundle.cells.values() for s in sweeps):
        pharm = "M1RaOxoM"
    # Genotype default to WT; specialize for KO lines.
    geno_map = {
        "tD-Tomato": "WT",
        "eGFP": "WT",
        "CDGIko X eGFP": "CDGIKO",
        "Adora2-Cre": "M1RCRISPR",
    }
    geno = geno_map.get(bundle.subject_entry.animal, "WT")
    # The gRNA-only controls of Figure 8 are Adora2-Cre mice without the Cas9 virus, so their M1R is intact
    if bundle.subject_entry.category != "M1R CRSPR" and geno == "M1RCRISPR":
        geno = "WT"
    return f"{cell_type}++{state}++{pharm}++{geno}++{bundle.subject_entry.date_compact}"


def _label_protocol_runs(irt, rows_per_protocol_run, sweep_per_row: dict[int, SweepRecord]) -> None:
    """Give each row the protocol run it belongs to (one sequential recording), labelled as NeuroConv's Bruker
    interface labels a run of its own: the acquisition stem of its first sweep. The interface wrote one label for
    all the rows of a series, whose sweeps span several protocol runs."""
    protocol_runs = irt["protocol_run"].data
    labels = []
    for rows in rows_per_protocol_run:
        label = _acquisition_stem(sweep_per_row[rows[0]])
        labels.append(label)
        for row in rows:
            protocol_runs[row] = label
    if len(set(labels)) != len(labels):
        raise ValueError(f"Two protocol runs of one file start with the same acquisition: {sorted(labels)}")


def merge_one_mouse_day(bundle: MouseDayBundle, output_path: Path) -> Path:
    """Build one per-mouse-day merged NWB file and write it to disk."""
    session_start_time_naive = bundle.session_start_time
    session_start_time = session_start_time_naive.replace(tzinfo=TIMEZONE)

    metadata = shared_metadata()
    subject_entry = bundle.subject_entry

    session_id = _build_session_id(bundle)
    session_description = _patch_clamp_session_description(bundle, cell_type=session_id.split("++")[0])
    surgery = metadata["NWBFile"]["surgery"]
    pharmacology = metadata["NWBFile"]["pharmacology"]
    if subject_entry.category == "NAIVE":
        surgery = _NAIVE_SURGERY
        pharmacology = _OXO_M_PHARMACOLOGY
    elif bundle.figure_code == "F6":
        pharmacology = pharmacology.rstrip() + _M1R_ANTAGONIST_PHARMACOLOGY[subject_entry.category]
    elif bundle.figure_code == "F8":
        viruses = _CRISPR_VIRUSES[subject_entry.category]
        surgery = surgery.rstrip() + _CRISPR_INJECTION.format(viruses=viruses) + _CRISPR_REPORTER
    drugs = {_STATE_AND_DRUG[_sweep_condition(sw)][1] for sweeps in bundle.cells.values() for sw in sweeps}
    for drug in sorted(drugs & _ADDED_DRUG_PHARMACOLOGY.keys()):
        pharmacology = pharmacology.rstrip() + _ADDED_DRUG_PHARMACOLOGY[drug]
    metadata["NWBFile"].update(
        session_description=session_description,
        experiment_description=_PATCH_CLAMP["NWBFile"]["experiment_description"],
        identifier=str(uuid.uuid4()),
        session_start_time=session_start_time,
        session_id=session_id,
        surgery=surgery,
        pharmacology=pharmacology,
    )
    metadata["Subject"].update(
        subject_id=subject_entry.subject_id,
        strain=subject_entry.strain,
        genotype=subject_entry.genotype,
        description=_patch_clamp_subject_description(bundle),
    )
    # The subjects table gives no strain or genotype for some lines, which the file then leaves out
    metadata["Subject"] = {key: value for key, value in metadata["Subject"].items() if value is not None}
    nwbfile = make_nwbfile_from_metadata(metadata)

    bruker, amp = _ensure_devices(nwbfile)

    # icephys chain bookkeeping
    # Key: (cell_id, protocol_family, condition) -> list of IRT row indices.
    # The condition is included so paired-recording cells (baseline + post-drug
    # sweeps on the same cell) split into separate Sequential/Repetition/ExperimentalConditions
    # rows rather than collapsing into one.
    irt_rows_per_sequential: dict[tuple[int, str, str], list[int]] = defaultdict(list)

    # Per-IRT-row metadata for the RecordingsIndexTable built post-chain.
    # Maps irt_row_index -> dict(cell_id, condition, dendrite_type, dendrite_distance_um,
    #                          site, stimulus_current_pA, sweep_start_time_s).
    index_row_metadata: dict[int, dict] = {}

    # The scan line at a location sometimes moves or changes length between trials.
    # scan_lines_per_location[(cell_id, location_type, location_number, spot)] maps each distinct
    # line (its end points) to its number at that location, in order of first use.
    # spot is "" for a site imaged at one spot, "Spot1" or "Spot2" for a site imaged at two.
    scan_lines_per_location: dict[tuple[int, str, int, str], dict[tuple, int]] = defaultdict(dict)

    # The dendritic trials of each (cell_id, location_type, location_number, spot, scan_line_number), in recording
    # order, as (IRT row, sweep, sweep start time); each group becomes one line-scan series per channel.
    line_scan_trials: dict[tuple[int, str, int, str, int], list[tuple[int, SweepRecord, float]]] = defaultdict(list)
    # The number of each trial at its location, in recording order, by IRT row; it names the source images.
    trial_number_at_location: dict[int, int] = {}

    # The intracellular recordings rows come from NeuroConv's Bruker interface, one series at a time; the columns
    # this script adds take their values once every row exists, from irt_values[row][column].
    irt = nwbfile.get_intracellular_recordings()
    irt_column_descriptions = _PATCH_CLAMP["Icephys"]["intracellular_recordings_columns"]
    irt_values: dict[int, dict] = {}
    sweep_per_row: dict[int, SweepRecord] = {}

    # Process each cell.
    for cell_index in bundle.cell_indices():
        # One pipette per cell, at the cell body, for both the somatic and the dendritic sweeps.
        cell_electrode = _build_electrode(nwbfile, cell_index, amp)
        # Somatic block (one PatchClampSeries for all somatic sweeps of this cell).
        soma_sweeps = bundle.sweeps_for(cell_index, "soma_excitability")
        if soma_sweeps:
            for sw, row_index, sweep_start_time, protocol in _add_patch_clamp_series(
                nwbfile, soma_sweeps, cell_index, "soma_excitability", cell_electrode, session_start_time
            ):
                condition = _CONDITION_LABEL[sw.condition_subfolder]
                irt_values[row_index] = {
                    "cell_id": cell_index,
                    "dendrite_type": "Soma",
                    "dendrite_distance_um": 0.0,
                    "sweep_start_time": sweep_start_time,
                    "stimulus_current_pA": protocol.amplitude_pA,
                    "stimulus_protocol_name": protocol.name,
                    "condition": condition,
                    "state": _STATE_AND_DRUG[condition][0],
                    "drug": _STATE_AND_DRUG[condition][1],
                    "site": "Soma",
                }
                sweep_per_row[row_index] = sw
                irt_rows_per_sequential[(cell_index, "soma_fi_curve", condition)].append(row_index)
                # Per-row metadata for the RecordingsIndexTable (built post-chain).
                index_row_metadata[row_index] = {
                    "cell_id": cell_index,
                    "condition": condition,
                    "dendrite_type": "Soma",
                    "dendrite_distance_um": 0.0,
                    "site": "Soma",
                    "stimulus_current_pA": protocol.amplitude_pA,
                    "sweep_start_time_s": sweep_start_time,
                }

        # Dendritic block: one series for all dendritic sweeps of the cell (Pattern B), each row tagged with its site;
        # each site and condition becomes one SequentialRecordings row.
        dend_sweeps = bundle.sweeps_for(cell_index, "dend_excitability")
        if not dend_sweeps:
            continue
        # Sites the lab imaged at two spots on the dendrite: the second spot's folders say "real". A "real" site without
        # plain folders of the same site is the region the lab meant after naming a folder after the wrong one.
        sites_with_two_spots = {
            (sw.dendrite_location, sw.location_number) for sw in dend_sweeps if "real" in sw.site_label
        } & {(sw.dendrite_location, sw.location_number) for sw in dend_sweeps if "real" not in sw.site_label}
        for sw, row_index, sweep_start_time, protocol in _add_patch_clamp_series(
            nwbfile, dend_sweeps, cell_index, "dend_excitability", cell_electrode, session_start_time
        ):
            dendrite_type = "Proximal" if sw.dendrite_location == "prox" else "Distal"
            dendrite_distance = 40.0 if sw.dendrite_location == "prox" else 90.0
            spot = ""
            if (sw.dendrite_location, sw.location_number) in sites_with_two_spots:
                spot = "Spot2" if "real" in sw.site_label else "Spot1"
            site = f"{dendrite_type}{sw.location_number}{spot}"
            condition = _CONDITION_LABEL[sw.condition_subfolder]
            if sw.suffix == "oxoM":
                condition += " with oxo-M"

            irt_values[row_index] = {
                "cell_id": cell_index,
                "dendrite_type": dendrite_type,
                "dendrite_distance_um": dendrite_distance,
                "sweep_start_time": sweep_start_time,
                "stimulus_current_pA": protocol.amplitude_pA,
                "stimulus_protocol_name": protocol.name,
                "condition": condition,
                "state": _STATE_AND_DRUG[condition][0],
                "drug": _STATE_AND_DRUG[condition][1],
                "site": site,
            }
            sweep_per_row[row_index] = sw
            protocol_label = f"dend_{sw.dendrite_location}{sw.location_number}"
            irt_rows_per_sequential[(cell_index, protocol_label, condition)].append(row_index)
            index_row_metadata[row_index] = {
                "cell_id": cell_index,
                "condition": condition,
                "dendrite_type": dendrite_type,
                "dendrite_distance_um": dendrite_distance,
                "site": site,
                "stimulus_current_pA": protocol.amplitude_pA,
                "sweep_start_time_s": sweep_start_time,
            }

            # Both channels share the scan line. Each spot of a site imaged at two spots is a location of its own.
            location = (cell_index, sw.dendrite_location, sw.location_number, spot)
            scan_lines = scan_lines_per_location[location]
            scan_line_number = scan_lines.setdefault(
                BrukerLineScanInterface.get_scan_line(sw.recording_path), len(scan_lines) + 1
            )
            line_scan_trials[(*location, scan_line_number)].append((row_index, sw, sweep_start_time))
            trial_number_at_location[row_index] = sum(
                len(trials) for key, trials in line_scan_trials.items() if key[:4] == location
            )

    _label_protocol_runs(irt, irt_rows_per_sequential.values(), sweep_per_row)
    for name, description in irt_column_descriptions.items():
        irt.add_column(name=name, description=description, data=[irt_values[row][name] for row in range(len(irt))])

    # ---- Line scans ---------------------------------------------------------------
    # One line scan and one fluorescence series per (cell, location, scan line, channel), each trial at its own
    # start. Each trial's samples in them go to the RecordingsIndexTable:
    # irt_row_index -> channel_name -> (series, idx_start, count).
    line_scan_slices: dict[int, dict[str, tuple[TimeSeries, int, int]]] = defaultdict(dict)
    roi_response_slices: dict[int, dict[str, tuple[TimeSeries, int, int]]] = defaultdict(dict)
    for (cell_index, loc_type, loc_num, spot, scan_line_number), trials in line_scan_trials.items():
        location = (cell_index, loc_type, loc_num, spot)
        trials_at_location = sum(len(entries) for key, entries in line_scan_trials.items() if key[:4] == location)
        dendrite_type_pretty = "Proximal" if loc_type == "prox" else "Distal"
        spot_note = _LINE_SCANS["spot_note"].format(spot=spot[-1]) if spot else ""
        site_note = _DENDRITIC_TRIAL + _spot_sentence(spot)
        line_suffix = ""
        n_scan_lines = len(scan_lines_per_location[location])
        if n_scan_lines > 1:
            line_suffix = f"Line{scan_line_number}"
            site_note += " " + _LINE_SCANS["moved_line"].format(n_lines=n_scan_lines, line=scan_line_number)
        image_text = _LINE_SCANS["source_image"]
        for chan_name in ("Ch1", "Ch2"):
            dye = _DYE_NAME[chan_name]
            images = [
                (
                    image_text["name"].format(
                        dye=dye,
                        cell=cell_index,
                        region=dendrite_type_pretty,
                        site_number=loc_num,
                        spot=spot,
                        trial=trial_number_at_location[row_index],
                    ),
                    image_text["description"].format(
                        cell=cell_index,
                        region=dendrite_type_pretty,
                        site_number=loc_num,
                        spot_note=spot_note,
                        dye_description=_DYE_DESCRIPTION[chan_name],
                        trial=trial_number_at_location[row_index],
                        n_trials=trials_at_location,
                        row=row_index,
                        start_time=sweep_start,
                    ),
                )
                for row_index, _, sweep_start in trials
            ]

            def series_note(column: str) -> str:
                return " " + _LINE_SCANS["series_note"].format(column=column, site_note=site_note)

            kymograph, roi_response, sample_ranges = _add_line_scan(
                nwbfile,
                sweeps=[sw for _, sw, _ in trials],
                sweep_start_times=[sweep_start for _, _, sweep_start in trials],
                channel=chan_name,
                name=f"{dye}Cell{cell_index}{dendrite_type_pretty}{loc_num}{spot}{line_suffix}",
                roi_name=f"Cell{cell_index}{dendrite_type_pretty}{loc_num}{spot}{line_suffix}",
                about=_LINE_SCANS["about"].format(
                    cell=cell_index,
                    region=dendrite_type_pretty,
                    site_number=loc_num,
                    dye_description=_DYE_DESCRIPTION[chan_name],
                ),
                kymograph_note=series_note(f"line_scan_{dye.lower()}"),
                roi_response_note=series_note(f"roi_response_{dye.lower()}"),
                images=images,
            )
            for row_index, sw, _ in trials:
                idx_start, count = sample_ranges[sw.recording_path.name]
                line_scan_slices[row_index][chan_name] = (kymograph, idx_start, count)
                roi_response_slices[row_index][chan_name] = (roi_response, idx_start, count)

    # Background scans (`_bsl`): one line scan per site with the line moved off the dendrite and no stimulus.
    # They are not trials, so they stay out of the intracellular tables.
    for cell_index in bundle.cell_indices():
        for sw in bundle.sweeps_for(cell_index, "dend_background"):
            site = f"{'Proximal' if sw.dendrite_location == 'prox' else 'Distal'}{sw.location_number}"
            background = _LINE_SCANS["background"]
            background_note = " " + background["note"]
            for chan_name in ("Ch1", "Ch2"):
                dye = _DYE_NAME[chan_name]
                _add_line_scan(
                    nwbfile,
                    sweeps=[sw],
                    sweep_start_times=[_recording_start(sw, session_start_time)],
                    channel=chan_name,
                    name=f"Background{dye}Cell{cell_index}{site}",
                    roi_name=f"BackgroundCell{cell_index}{site}",
                    about=background["about"].format(
                        cell=cell_index, site=site, dye_description=_DYE_DESCRIPTION[chan_name]
                    ),
                    kymograph_note=background_note,
                    roi_response_note=background_note,
                    images=[
                        (
                            background["source_image"]["name"].format(dye=dye, cell=cell_index, site=site),
                            background["source_image"]["description"].format(
                                cell=cell_index, site=site, dye_description=_DYE_DESCRIPTION[chan_name]
                            ),
                        )
                    ],
                )

    # ---- icephys hierarchical chain --------------------------------------------
    # Build the chain on top of the IRT rows.
    # 1. SimultaneousRecordings: one row per IRT row (single-electrode rig).
    simultaneous_indices_per_irt: dict[int, int] = {}
    for irt_indices in irt_rows_per_sequential.values():
        for irt_index in irt_indices:
            sim_index = nwbfile.add_icephys_simultaneous_recording(recordings=[irt_index])
            simultaneous_indices_per_irt[irt_index] = sim_index

    # 2. SequentialRecordings: one per (cell, protocol family, condition).
    # Paired-recording cells with both baseline and post-drug sweeps will get TWO
    # sequential rows for their soma_fi_curve protocol (one per condition).
    sequential_indices_per_key: dict[tuple[int, str, str], int] = {}
    sequential_index_per_irt: dict[int, int] = {}
    for key, irt_indices in irt_rows_per_sequential.items():
        sim_indices = [simultaneous_indices_per_irt[i] for i in irt_indices]
        cell_id, protocol, condition = key
        stimulus_type = _PATCH_CLAMP["Icephys"]["stimulus_type"]["soma" if "soma" in protocol else "dendritic"]
        seq_index = nwbfile.add_icephys_sequential_recording(
            simultaneous_recordings=sim_indices,
            stimulus_type=stimulus_type,
        )
        sequential_indices_per_key[key] = seq_index
        for irt_index in irt_indices:
            sequential_index_per_irt[irt_index] = seq_index

    # 3. Repetitions: one per (cell, protocol family, condition).
    repetitions_table = nwbfile.get_icephys_repetitions()
    for name, description in _PATCH_CLAMP["Icephys"]["repetitions_columns"].items():
        repetitions_table.add_column(name=name, description=description)
    # Track which repetitions belong to which condition for the ExperimentalConditions table.
    repetitions_per_condition: dict[str, list[int]] = defaultdict(list)
    repetition_index_per_irt: dict[int, int] = {}
    for key, seq_index in sequential_indices_per_key.items():
        cell_id, protocol, condition = key
        rep_index = nwbfile.add_icephys_repetition(
            sequential_recordings=[seq_index],
            cell_id=cell_id,
            protocol=_protocol_description(protocol),
            condition=condition,
        )
        repetitions_per_condition[condition].append(rep_index)
        for irt_index in irt_rows_per_sequential[key]:
            repetition_index_per_irt[irt_index] = rep_index

    # 4. ExperimentalConditions: one row per unique condition across the mouse-day.
    # For non-paired mouse-days this collapses to a single row. For paired (e.g., ON+SCH)
    # the table has two rows, each pointing to the repetitions acquired under that condition.
    ec_table = nwbfile.get_icephys_experimental_conditions()
    if "condition" not in ec_table.colnames:
        ec_table.add_column(
            name="condition", description=_PATCH_CLAMP["Icephys"]["experimental_conditions_columns"]["condition"]
        )
    ec_index_per_condition: dict[str, int] = {}
    for ec_index, (condition, rep_indices) in enumerate(repetitions_per_condition.items()):
        nwbfile.add_icephys_experimental_condition(
            repetitions=rep_indices,
            condition=condition,
        )
        ec_index_per_condition[condition] = ec_index

    # ---- RecordingsIndexTable -------------------------------------------------
    # Denormalized lookup that pre-joins the icephys chain for analyst convenience.
    # One row per IRT row; carries all the cross-references and per-sweep metadata
    # in flat columns. The canonical chain is still authoritative; this is a query
    # helper.
    _build_recordings_index_table(
        nwbfile,
        irt_table=irt,
        sr_table=nwbfile.icephys_simultaneous_recordings,
        seq_table=nwbfile.icephys_sequential_recordings,
        reps_table=nwbfile.icephys_repetitions,
        ec_table=nwbfile.icephys_experimental_conditions,
        index_row_metadata=index_row_metadata,
        simultaneous_index_per_irt=simultaneous_indices_per_irt,
        sequential_index_per_irt=sequential_index_per_irt,
        repetition_index_per_irt=repetition_index_per_irt,
        ec_index_per_condition=ec_index_per_condition,
        line_scan_slices=line_scan_slices,
        roi_response_slices=roi_response_slices,
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    # NeuroConv's helper applies its default HDF5 chunking and compression
    configure_and_write_nwbfile(nwbfile=nwbfile, nwbfile_path=output_path, backend="hdf5", metadata=shared_metadata())
    return output_path


def _build_recordings_index_table(
    nwbfile,
    *,
    irt_table,
    sr_table,
    seq_table,
    reps_table,
    ec_table,
    index_row_metadata: dict[int, dict],
    simultaneous_index_per_irt: dict[int, int],
    sequential_index_per_irt: dict[int, int],
    repetition_index_per_irt: dict[int, int],
    ec_index_per_condition: dict[str, int],
    line_scan_slices: dict[int, dict[str, tuple[TimeSeries, int, int]]],
    roi_response_slices: dict[int, dict[str, tuple]] | None = None,
) -> None:
    """Build and attach the RecordingsIndexTable.

    The table denormalizes one row per IRT row with the cross-references from
    the icephys hierarchical chain plus paired line-scan slice references
    (raw kymographs AND ROI responses).
    For sweeps without line scans (somatic), line_scan_alexa568/fluo4 and
    roi_response_alexa568/fluo4 entries use idx_start=-1, count=-1 (null-reference).
    """
    roi_response_slices = roi_response_slices or {}
    n_rows = len(irt_table)
    if n_rows == 0:
        return

    # Find any line-scan TimeSeries to use as a placeholder for null entries.
    # TimeSeriesReferenceVectorData requires a valid TimeSeries field even for
    # null entries; we use idx_start=-1 to signal "no actual reference here".
    fallback_ts: TimeSeries | None = None
    for slices in line_scan_slices.values():
        for ts, _, _ in slices.values():
            fallback_ts = ts
            break
        if fallback_ts is not None:
            break

    # Build columns
    cell_id_col: list[int] = []
    condition_col: list[str] = []
    dendrite_type_col: list[str] = []
    dendrite_distance_col: list[float] = []
    site_col: list[str] = []
    stim_pA_col: list[float] = []
    sweep_start_col: list[float] = []
    irt_indices_col: list[int] = []
    sr_indices_col: list[int] = []
    seq_indices_col: list[int] = []
    rep_indices_col: list[int] = []
    ec_indices_col: list[int] = []
    line_scan_ch1_refs: list[TimeSeriesReference] = []
    line_scan_ch2_refs: list[TimeSeriesReference] = []
    roi_response_ch1_refs: list[TimeSeriesReference] = []
    roi_response_ch2_refs: list[TimeSeriesReference] = []

    # Find a fallback ROI response TimeSeries for null entries (same pattern as kymo fallback).
    fallback_rrs = None
    for slices in roi_response_slices.values():
        for ts, _, _ in slices.values():
            fallback_rrs = ts
            break
        if fallback_rrs is not None:
            break

    for irt_index in range(n_rows):
        meta = index_row_metadata.get(irt_index, {})
        cell_id_col.append(int(meta.get("cell_id", -1)))
        condition_col.append(str(meta.get("condition", "")))
        dendrite_type_col.append(str(meta.get("dendrite_type", "")))
        dendrite_distance_col.append(float(meta.get("dendrite_distance_um", 0.0)))
        site_col.append(str(meta.get("site", "")))
        stim_pA_col.append(float(meta.get("stimulus_current_pA", 0.0)))
        sweep_start_col.append(float(meta.get("sweep_start_time_s", 0.0)))
        irt_indices_col.append(irt_index)
        sr_indices_col.append(int(simultaneous_index_per_irt.get(irt_index, -1)))
        seq_indices_col.append(int(sequential_index_per_irt.get(irt_index, -1)))
        rep_indices_col.append(int(repetition_index_per_irt.get(irt_index, -1)))
        ec_indices_col.append(int(ec_index_per_condition.get(meta.get("condition", ""), -1)))

        slices = line_scan_slices.get(irt_index, {})
        for chan_name, target_list in (("Ch1", line_scan_ch1_refs), ("Ch2", line_scan_ch2_refs)):
            entry = slices.get(chan_name)
            if entry is not None:
                ts, idx_start, count = entry
                target_list.append(TimeSeriesReference(int(idx_start), int(count), ts))
            elif fallback_ts is not None:
                # Null reference: idx_start=-1, count=-1
                target_list.append(TimeSeriesReference(-1, -1, fallback_ts))
            else:
                target_list.append(None)  # type: ignore[arg-type]

        roi_slices = roi_response_slices.get(irt_index, {})
        for chan_name, target_list in (("Ch1", roi_response_ch1_refs), ("Ch2", roi_response_ch2_refs)):
            entry = roi_slices.get(chan_name)
            if entry is not None:
                ts, idx_start, count = entry
                target_list.append(TimeSeriesReference(int(idx_start), int(count), ts))
            elif fallback_rrs is not None:
                target_list.append(TimeSeriesReference(-1, -1, fallback_rrs))
            else:
                target_list.append(None)  # type: ignore[arg-type]

    # Build the DynamicTable. Need to construct columns before passing to
    # DynamicTable so DynamicTableRegion can reference the icephys tables.
    text = _PATCH_CLAMP["Icephys"]["RecordingsIndexTable"]
    described = text["columns"]
    columns: list = [
        VectorData(name="cell_id", description=described["cell_id"], data=cell_id_col),
        VectorData(name="condition", description=described["condition"], data=condition_col),
        VectorData(name="dendrite_type", description=described["dendrite_type"], data=dendrite_type_col),
        VectorData(
            name="dendrite_distance_um", description=described["dendrite_distance_um"], data=dendrite_distance_col
        ),
        VectorData(name="site", description=described["site"], data=site_col),
        VectorData(name="stimulus_current_pA", description=described["stimulus_current_pA"], data=stim_pA_col),
        VectorData(name="sweep_start_time_s", description=described["sweep_start_time_s"], data=sweep_start_col),
        DynamicTableRegion(name="irt_row", description=described["irt_row"], data=irt_indices_col, table=irt_table),
        DynamicTableRegion(
            name="simultaneous_row", description=described["simultaneous_row"], data=sr_indices_col, table=sr_table
        ),
        DynamicTableRegion(
            name="sequential_row", description=described["sequential_row"], data=seq_indices_col, table=seq_table
        ),
        DynamicTableRegion(
            name="repetition_row", description=described["repetition_row"], data=rep_indices_col, table=reps_table
        ),
        DynamicTableRegion(
            name="experimental_condition_row",
            description=described["experimental_condition_row"],
            data=ec_indices_col,
            table=ec_table,
        ),
    ]

    # Line-scan columns only if at least one sweep has line-scan slices.
    if fallback_ts is not None:
        for name, data in (("line_scan_alexa568", line_scan_ch1_refs), ("line_scan_fluo4", line_scan_ch2_refs)):
            columns.append(TimeSeriesReferenceVectorData(name=name, description=described[name], data=data))
    if fallback_rrs is not None:
        for name, data in (
            ("roi_response_alexa568", roi_response_ch1_refs),
            ("roi_response_fluo4", roi_response_ch2_refs),
        ):
            columns.append(TimeSeriesReferenceVectorData(name=name, description=described[name], data=data))

    table = DynamicTable(name="RecordingsIndexTable", description=text["description"], columns=columns)
    module = nwbfile.create_processing_module(name="icephys", description=text["processing_module_description"])
    module.add(table)
