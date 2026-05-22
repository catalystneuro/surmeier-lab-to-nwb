"""Prototype: build one per-mouse-day merged NWB file for Figure 1.

Demonstrates the per-mouse-day plan end-to-end on a specific test mouse-day
(Subject20170202LIDOFFtDTomato, 3 patched cells, somatic + partial dendritic).

What this prototype DOES:

- Subject from the Data Connections registry
- Devices: Bruker Ultima + Multiclamp 700B (linked to DeviceModel per devices.py)
- One IntracellularElectrode per (cell, compartment), all under one Subject
- Pattern B PatchClampSeries: one CurrentClampSeries per (cell, modality) with all sweeps
  concatenated end-to-end (no gap padding) and matching CurrentClampStimulusSeries
- IntracellularRecordingsTable with: electrode, stimulus, response (idx_start + count),
  plus custom columns dendrite_type, dendrite_distance_um, cell_id, sweep_start_time
- Full 5/5 hierarchical chain:
    ExperimentalConditions -> Repetitions -> SequentialRecordings -> SimultaneousRecordings -> IntracellularRecordings
- sweep_start_time column populated from each sweep's PVScan @date

What this prototype SKIPS for now (add later):
- Line-scan TwoPhotonSeries paired with dendritic sweeps
- The stimulus current values (just stores zeros for the stimulus series)
- Detailed per-figure metadata in surgery/pharmacology prose
"""

from __future__ import annotations

import logging
import re
import uuid
import warnings
import xml.etree.ElementTree as ET
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from zoneinfo import ZoneInfo

# Bruker Prairie View leaves UIC CreateTime/LastSavedTime as julianday=0 in the
# line-scan TIFFs. tifffile warns when parsing these unparsable fields, but the
# pixel data reads correctly and the real acquisition timestamp lives in the
# PVScan @date attribute of the master XML (which we read separately). Suppress
# the noise.
logging.getLogger("tifffile").setLevel(logging.ERROR)
warnings.filterwarnings("ignore", message=".*no datetime before year 1.*")

import numpy as np
import pandas as pd
from hdmf.common import DynamicTable, DynamicTableRegion, VectorData
from neuroconv.utils import load_dict_from_file
from pynwb import NWBHDF5IO, NWBFile
from pynwb.base import TimeSeries, TimeSeriesReference, TimeSeriesReferenceVectorData
from pynwb.device import Device
from pynwb.file import Subject
from pynwb.icephys import (
    CurrentClampSeries,
    CurrentClampStimulusSeries,
    IntracellularElectrode,
)

from surmeier_lab_to_nwb.zhai2025.devices import get_or_create_bruker_ultima_model
from surmeier_lab_to_nwb.zhai2025.mouse_day_merger import (
    MouseDayBundle,
)

REPO_ROOT = Path(__file__).resolve().parents[1]
TIMEZONE = ZoneInfo("America/Chicago")  # Northwestern University

GENERAL_METADATA_PATH = REPO_ROOT / "src" / "surmeier_lab_to_nwb" / "zhai2025" / "general_metadata.yaml"
_GENERAL_METADATA_CACHE: dict | None = None


def _general_metadata() -> dict:
    global _GENERAL_METADATA_CACHE
    if _GENERAL_METADATA_CACHE is None:
        _GENERAL_METADATA_CACHE = load_dict_from_file(GENERAL_METADATA_PATH)
    return _GENERAL_METADATA_CACHE


# --------- reading per-sweep line-scan kymographs -------------------------------------------


def _line_scan_kymographs(sweep_path: Path) -> dict[str, tuple[np.ndarray, float | None]]:
    """Read line-scan kymograph TIFFs for one dendritic sweep.

    Returns {channel_name: (data, rate_hz_or_None)}. Channels are typically:
        Ch1 = Alexa Fluor 568 reference (red) — R0 for ratiometric calculations
        Ch2 = Fluo-4 calcium indicator (green) — G in (G-G0)/G0R0
    """
    import tifffile

    out: dict[str, tuple[np.ndarray, float | None]] = {}
    for channel in ("Ch1", "Ch2"):
        candidates = list(sweep_path.glob(f"*Cycle00001_{channel}_*.ome.tif"))
        if not candidates:
            continue
        tif = candidates[0]
        try:
            arr = tifffile.imread(tif)
        except Exception:
            continue
        # Some recordings package both channels in one 3D TIFF (channels, lines, pixels).
        if arr.ndim == 3:
            channel_index = {"Ch1": 0, "Ch2": 1}.get(channel)
            if channel_index is None or channel_index >= arr.shape[0]:
                continue
            arr = arr[channel_index]
        if arr.ndim != 2:
            continue
        # Sample rate: from PrairieView's scanLinePeriod in the master XML.
        rate = _read_scan_line_rate(sweep_path)
        out[channel] = (arr, rate)
    return out


def _read_scan_line_rate(sweep_path: Path) -> float | None:
    """Extract the line-scan sample rate (lines/second) from the master XML.

    Returns 1/scanLinePeriod if available, otherwise None (in which case the
    caller should fall back to timestamps).
    """
    xml_path = sweep_path / f"{sweep_path.name}.xml"
    if not xml_path.exists():
        return None
    try:
        text = xml_path.read_text(errors="replace", encoding="utf-8")
    except OSError:
        return None
    match = re.search(r'\bkey="scanLinePeriod"[^/>]*\bvalue="([^"]+)"', text)
    if not match:
        return None
    try:
        period = float(match.group(1))
    except (TypeError, ValueError):
        return None
    return 1.0 / period if period > 0 else None


# --------- reading per-sweep stimulus protocol ---------------------------------------------


@dataclass
class StimulusProtocol:
    """Parsed current-step stimulus protocol from a VoltageOutput XML."""

    name: str  # e.g., "Step 5_-40 pA"
    amplitude_pA: float  # commanded step amplitude in picoamperes
    first_pulse_delay_ms: float  # delay before the step in milliseconds
    pulse_width_ms: float  # step duration in milliseconds


def _parse_voltage_output(xml_path: Path) -> StimulusProtocol | None:
    """Parse a Bruker VoltageOutput XML to extract the current-clamp protocol step.

    Returns a `StimulusProtocol` for the first enabled '700B' waveform component,
    or None if the file can't be parsed.
    """
    try:
        tree = ET.parse(xml_path)
    except (ET.ParseError, OSError):
        return None
    root = tree.getroot()
    name = root.findtext("Name", default="").strip()
    # Find the '700B' waveform (the patch amplifier output).
    for waveform in root.findall("Waveform"):
        if waveform.findtext("Name") != "700B":
            continue
        if waveform.findtext("Enabled", default="false").lower() != "true":
            continue
        scale_factor = float(waveform.findtext("UnitScaleFactor", default="1"))
        # The first PulseTrain component carries the step amplitude.
        pulse_train = waveform.find("WaveformComponent_PulseTrain")
        if pulse_train is None:
            continue
        try:
            potential_volts = float(pulse_train.findtext("PulsePotentialStart", default="0"))
            pulse_width_ms = float(pulse_train.findtext("PulseWidth", default="0"))
            first_pulse_delay_ms = float(pulse_train.findtext("FirstPulseDelay", default="0"))
        except (TypeError, ValueError):
            continue
        amplitude_pA = potential_volts * scale_factor
        return StimulusProtocol(
            name=name,
            amplitude_pA=amplitude_pA,
            first_pulse_delay_ms=first_pulse_delay_ms,
            pulse_width_ms=pulse_width_ms,
        )
    return None


def _find_voltage_output_xml(sweep_path: Path) -> Path | None:
    for candidate in sweep_path.glob("*_VoltageOutput_001.xml"):
        return candidate
    return None


def _synthesize_stimulus_waveform(
    protocol: StimulusProtocol,
    n_samples: int,
    rate_hz: float,
) -> np.ndarray:
    """Build the commanded current waveform (in amperes) from a parsed protocol.

    Layout: zeros for `first_pulse_delay_ms`, then amplitude for `pulse_width_ms`,
    then zeros to the end of the sweep. All converted from pA to amperes.
    """
    waveform = np.zeros(n_samples, dtype=np.float32)
    start_sample = int(round(protocol.first_pulse_delay_ms * 1e-3 * rate_hz))
    width_samples = int(round(protocol.pulse_width_ms * 1e-3 * rate_hz))
    end_sample = min(start_sample + width_samples, n_samples)
    if start_sample < n_samples:
        waveform[start_sample:end_sample] = protocol.amplitude_pA * 1e-12
    return waveform


# --------- reading per-sweep voltage data ---------------------------------------------------


def _read_signal_gain(sweep_path: Path) -> tuple[float, float]:
    """Parse the Primary signal's gain (multiplier/divisor) from the VoltageRecording XML.

    Returns (multiplier, divisor). Defaults to (1, 1) if not parseable.
    The conversion is `value_mV = raw * (multiplier / divisor)`.
    """
    xml_candidates = list(sweep_path.glob("*_VoltageRecording_001.xml"))
    if not xml_candidates:
        return 1.0, 1.0
    try:
        text = xml_candidates[0].read_text(errors="replace", encoding="utf-8")
    except OSError:
        return 1.0, 1.0
    # Match the FIRST Unit block (Primary signal). The pattern captures the
    # multiplier and divisor that follow the Unit tag.
    m = re.search(
        r"<Unit>\s*<UnitName>mV</UnitName>\s*<Multiplier>([^<]+)</Multiplier>\s*<Divisor>([^<]+)</Divisor>", text
    )
    if not m:
        return 1.0, 1.0
    try:
        return float(m.group(1)), float(m.group(2))
    except (TypeError, ValueError):
        return 1.0, 1.0


def _read_sweep_voltage(sweep_path: Path) -> tuple[np.ndarray, float]:
    """Read voltage trace and sample rate for one sweep.

    Applies the multiplier/divisor gain conversion from the VoltageRecording XML
    so the stored data is in true volts (Prairie View convention is to store the
    Primary channel scaled by 1/multiplier × divisor; we invert that and then
    convert mV -> volts for NWB storage).

    Returns (voltage_volts, sample_rate_hz).
    """
    csv_candidates = list(sweep_path.glob("*_VoltageRecording_001.csv"))
    if not csv_candidates:
        csv_candidates = list(sweep_path.glob("*.csv"))
    if not csv_candidates:
        raise FileNotFoundError(f"No voltage CSV in {sweep_path}")
    csv_path = csv_candidates[0]
    df = pd.read_csv(csv_path, header=0)
    time_ms = df.iloc[:, 0].values
    raw_primary = df.iloc[:, 1].values
    multiplier, divisor = _read_signal_gain(sweep_path)
    # Apply gain conversion to get true mV, then convert to volts for NWB.
    voltage_mV = raw_primary * (multiplier / divisor) if divisor != 0 else raw_primary
    voltage_v = voltage_mV / 1000.0
    if len(time_ms) < 2:
        sample_rate = 10000.0
    else:
        dt_s = (time_ms[1] - time_ms[0]) / 1000.0
        sample_rate = round(1.0 / dt_s)
    return voltage_v.astype(np.float32), float(sample_rate)


# --------- core merger logic ---------------------------------------------------


def _electrode_name(cell_index: int, compartment: str, location_number: int | None = None) -> str:
    """Stable name for an IntracellularElectrode."""
    if compartment == "soma":
        return f"Cell{cell_index}_Soma"
    if compartment == "prox":
        return f"Cell{cell_index}_Proximal{location_number}"
    if compartment == "dist":
        return f"Cell{cell_index}_Distal{location_number}"
    raise ValueError(f"Unknown compartment: {compartment}")


def _series_name(cell_index: int, modality: str) -> str:
    suffix = "Somatic" if modality == "soma_excitability" else "Dendritic"
    return f"CurrentClampSeries_Cell{cell_index}_{suffix}"


def _condition_label(condition_subfolder: str) -> str:
    """Map raw condition subfolder name to the paper's prose label."""
    return {
        "LID off-state": "LID off-state",
        "LID on-state": "LID on-state",
        "LID on-state with SCH": "LID on-state with SCH",
    }.get(condition_subfolder, condition_subfolder)


def _ensure_devices(nwbfile: NWBFile) -> tuple[Device, Device]:
    """Add Bruker microscope and Multiclamp amplifier devices."""
    bruker_model = get_or_create_bruker_ultima_model(nwbfile)
    bruker_name = "BrukerUltima"
    if bruker_name in nwbfile.devices:
        bruker = nwbfile.devices[bruker_name]
    else:
        bruker = Device(
            name=bruker_name,
            description="Bruker Ultima In Vitro multiphoton microscope (mouse-day merge prototype).",
            model=bruker_model,
        )
        nwbfile.add_device(bruker)

    amp_name = "MultiClamp700B"
    if amp_name in nwbfile.devices:
        amp = nwbfile.devices[amp_name]
    else:
        amp = Device(
            name=amp_name,
            description="Molecular Devices MultiClamp 700B patch-clamp amplifier.",
        )
        nwbfile.add_device(amp)
    return bruker, amp


def _build_electrode(
    nwbfile: NWBFile,
    cell_index: int,
    compartment: str,
    location_number: int | None,
    device: Device,
) -> IntracellularElectrode:
    name = _electrode_name(cell_index, compartment, location_number)
    if name in nwbfile.icephys_electrodes:
        return nwbfile.icephys_electrodes[name]
    if compartment == "soma":
        description = (
            f"Patch-clamp electrode at the soma of Cell {cell_index}, "
            "recorded as part of the mouse-day icephys merge."
        )
        location = "Caudoputamen"
    elif compartment == "prox":
        description = (
            f"Patch-clamp electrode at proximal dendrite location {location_number} "
            f"of Cell {cell_index} (~40 um from soma)."
        )
        location = "Caudoputamen"
    else:  # dist
        description = (
            f"Patch-clamp electrode at distal dendrite location {location_number} "
            f"of Cell {cell_index} (~90 um from soma)."
        )
        location = "Caudoputamen"
    electrode = nwbfile.create_icephys_electrode(
        name=name,
        description=description,
        device=device,
        location=location,
        cell_id=f"Cell{cell_index}",
    )
    return electrode


def _build_concatenated_series(
    nwbfile: NWBFile,
    bundle: MouseDayBundle,
    cell_index: int,
    modality: str,
    electrode: IntracellularElectrode,
    session_start_time,
) -> (
    tuple[
        CurrentClampSeries,
        CurrentClampStimulusSeries,
        list[tuple[int, int, float, StimulusProtocol | None]],
    ]
    | None
):
    """Concatenate all sweeps of (cell, modality) into one PatchClampSeries.

    Returns (response_series, stimulus_series, sweep_slices) where each slice carries
    (idx_start, count, sweep_start_time_seconds, stimulus_protocol_or_None). The stimulus
    protocol is parsed from the VoltageOutput XML alongside each sweep and used to
    populate the IRT's `stimulus_current_pA` column plus the synthesized stimulus waveform.
    """
    sweeps = bundle.sweeps_for(cell_index, modality)
    if not sweeps:
        return None

    voltage_chunks: list[np.ndarray] = []
    stim_chunks: list[np.ndarray] = []
    sweep_slices: list[tuple[int, int, float, StimulusProtocol | None]] = []
    sample_rate: float | None = None
    cumulative = 0
    for sw in sweeps:
        try:
            voltage, rate = _read_sweep_voltage(sw.recording_path)
        except FileNotFoundError:
            continue
        if sample_rate is None:
            sample_rate = rate
        elif abs(rate - sample_rate) > 1.0:
            print(f"WARNING: rate mismatch for {sw.recording_path.name} ({rate} vs {sample_rate}); skipping")
            continue
        # Parse the stimulus protocol from the VoltageOutput XML for this sweep.
        vo_xml = _find_voltage_output_xml(sw.recording_path)
        protocol = _parse_voltage_output(vo_xml) if vo_xml is not None else None
        if protocol is not None:
            stim_waveform = _synthesize_stimulus_waveform(
                protocol,
                n_samples=len(voltage),
                rate_hz=float(sample_rate),
            )
        else:
            stim_waveform = np.zeros_like(voltage)
        voltage_chunks.append(voltage)
        stim_chunks.append(stim_waveform)
        offset_s = (sw.pv_scan_date - session_start_time).total_seconds()
        sweep_slices.append((cumulative, len(voltage), offset_s, protocol))
        cumulative += len(voltage)

    if not voltage_chunks:
        return None

    data_concat = np.concatenate(voltage_chunks)
    stim_concat = np.concatenate(stim_chunks)
    series_name = _series_name(cell_index, modality)
    response = CurrentClampSeries(
        name=series_name,
        data=data_concat,
        electrode=electrode,
        gain=1.0,
        rate=float(sample_rate),
        starting_time=sweep_slices[0][2],
        description=(
            f"Concatenated current-clamp recordings for Cell {cell_index} {modality}. "
            f"{len(sweep_slices)} sweeps merged end-to-end (Pattern B). "
            "Per-sweep wall-clock offsets stored in IntracellularRecordingsTable.sweep_start_time."
        ),
    )
    nwbfile.add_acquisition(response)

    stim = CurrentClampStimulusSeries(
        name=f"{series_name}_Stimulus",
        data=stim_concat,
        electrode=electrode,
        gain=1.0,
        rate=float(sample_rate),
        starting_time=sweep_slices[0][2],
        description=(
            f"Commanded current waveform paired with {series_name}. Synthesized from each "
            "sweep's Bruker VoltageOutput XML (step amplitude × pulse width); the per-sweep "
            "amplitude is also exposed via the IntracellularRecordingsTable.stimulus_current_pA column."
        ),
    )
    nwbfile.add_stimulus(stim)

    return response, stim, sweep_slices


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
        "CONT": "Control",
        "6-OHDA": "LesionedControl",
        "M1R ant": "M1RAntag",
        "M1R CRSPR": "M1RCRISPR",
    }
    state = state_map.get(bundle.subject_entry.category, "Unknown")
    # Pharm from the most distinguishing condition in the bundle. If a paired-
    # recording cell exists, surface the drug used (SCH or sul).
    all_conditions = {s.condition_subfolder for sweeps in bundle.cells.values() for s in sweeps}
    pharm = "none"
    if any("SCH" in c for c in all_conditions):
        pharm = "D1RaSch"
    elif any("sul" in c for c in all_conditions):
        pharm = "D2RaSul"
    # Genotype default to WT; specialize for KO lines.
    geno_map = {
        "tD-Tomato": "WT",
        "eGFP": "WT",
        "CDGIko X eGFP": "CDGIKO",
        "Adora2-Cre": "M1RCRISPR",
    }
    geno = geno_map.get(bundle.subject_entry.animal, "WT")
    return f"{cell_type}++{state}++{pharm}++{geno}++{bundle.subject_entry.date_compact}"


def merge_one_mouse_day(bundle: MouseDayBundle, output_path: Path) -> Path:
    """Build one per-mouse-day merged NWB file and write it to disk."""
    session_start_time_naive = bundle.session_start_time
    session_start_time = session_start_time_naive.replace(tzinfo=TIMEZONE)

    general = _general_metadata()
    nwb_meta = general.get("NWBFile", {})
    subject_meta = general.get("Subject", {})
    subject_entry = bundle.subject_entry

    session_description = (
        f"Per-mouse-day icephys merge for {subject_entry.subject_id}: "
        f"{len(bundle.cell_indices())} patched cell(s), "
        f"category={subject_entry.category}, animal/reporter={subject_entry.animal}. "
        f"All cells from this mouse on this day collected in one NWB file with the icephys "
        f"hierarchical chain populated (IntracellularRecordings -> Simultaneous -> Sequential "
        f"-> Repetitions -> ExperimentalConditions)."
    )
    nwbfile = NWBFile(
        session_description=session_description,
        experiment_description=nwb_meta.get(
            "experiment_description",
            "Excitability and synaptic plasticity of striatal projection neurons in 6-OHDA-lesioned mice "
            "undergoing levodopa-induced dyskinesia.",
        ),
        identifier=str(uuid.uuid4()),
        session_start_time=session_start_time,
        session_id=_build_session_id(bundle),
        experimenter=nwb_meta.get("experimenter", ["Zhai, Shenyu"]),
        institution=nwb_meta.get("institution", "Northwestern University"),
        lab=nwb_meta.get("lab", "Surmeier Lab"),
        keywords=nwb_meta.get("keywords", []),
        related_publications=nwb_meta.get("related_publications", []),
        surgery=nwb_meta.get("surgery"),
        pharmacology=nwb_meta.get("pharmacology"),
        slices=nwb_meta.get("slices"),
    )
    nwbfile.subject = Subject(
        subject_id=subject_entry.subject_id,
        species=subject_meta.get("species", "Mus musculus"),
        strain=subject_meta.get("strain", "C57BL/6"),
        sex=subject_meta.get("sex", "M"),
        age=subject_meta.get("age", "P7W/P12W"),
        description=(
            f"Reporter line {subject_entry.animal}, category {subject_entry.category}, "
            "ipsilateral 6-OHDA lesion of the medial forebrain bundle. "
            f"See general_metadata.yaml surgery section for the lesion protocol."
        ),
    )

    _bruker, amp = _ensure_devices(nwbfile)

    # icephys chain bookkeeping
    # Key: (cell_id, protocol_family, condition_subfolder) -> list of IRT row indices.
    # The condition subfolder is included so paired-recording cells (baseline + post-drug
    # sweeps on the same cell) split into separate Sequential/Repetition/ExperimentalConditions
    # rows rather than collapsing into one.
    irt_rows_per_sequential: dict[tuple[int, str, str], list[int]] = defaultdict(list)

    # Per-IRT-row metadata for the RecordingsIndexTable built post-chain.
    # Maps irt_row_idx -> dict(cell_id, condition_subfolder, dendrite_type, dendrite_distance_um,
    #                          stimulus_current_pA, sweep_start_time_s).
    index_row_metadata: dict[int, dict] = {}

    # Buffer for per-sweep line-scan kymographs to consolidate later.
    # Key: (cell_id, location_type, location_number, channel_name)
    # Value: list of (irt_row_idx, kymo_data_2d, rate_hz_or_None, sweep_start_time_s)
    # After all dendritic IRT rows are added, we concatenate each buffer along axis 0
    # (the time/lines axis) to make one TimeSeries per (cell, location, channel) -- Pattern B
    # for line scans. Per-sweep slice info is preserved for the RecordingsIndexTable.
    line_scan_buffer: dict[tuple[int, str, int, str], list[tuple[int, np.ndarray, float | None, float]]] = defaultdict(
        list
    )

    # Add custom columns to IntracellularRecordingsTable BEFORE adding rows.
    # The `electrodes` aligned sub-table doesn't exist until the first add; using
    # the convenience add_intracellular_recording call below pre-creates it.
    irt = nwbfile.get_intracellular_recordings()
    irt.add_column(name="cell_id", description="Cell index within the mouse-day (1, 2, 3, 4).")
    irt.add_column(
        name="dendrite_type",
        description="Compartment of the recording electrode: Soma, Proximal, or Distal.",
    )
    irt.add_column(
        name="dendrite_distance_um",
        description="Approximate distance from soma in micrometers (0 for soma).",
    )
    irt.add_column(
        name="sweep_start_time",
        description=(
            "Per-sweep wall-clock offset (seconds) from session_start_time. "
            "Sourced from the Prairie View XML PVScan @date attribute at 1-second "
            "precision. Inter-sweep gaps (amplifier idle between sweeps) are not "
            "carried in the parent PatchClampSeries' uniform rate."
        ),
    )
    irt.add_column(
        name="stimulus_current_pA",
        description=(
            "Commanded current-step amplitude in picoamperes for this sweep. Parsed from "
            "the Bruker VoltageOutput XML alongside each sweep (PulsePotentialStart × UnitScaleFactor). "
            "0 if the protocol could not be parsed."
        ),
    )
    irt.add_column(
        name="stimulus_protocol_name",
        description="Stimulus protocol label from Bruker VoltageOutput XML (e.g., 'Step 5_-40 pA').",
    )
    irt.add_column(
        name="condition_subfolder",
        description=(
            "Raw condition subfolder this sweep came from (e.g., 'LID off-state', 'LID on-state', "
            "'LID on-state with SCH'). For paired-recording mouse-days, distinguishes baseline vs "
            "post-drug sweeps within the same file."
        ),
    )

    # Process each cell.
    for cell_index in bundle.cell_indices():
        # Somatic block (one PatchClampSeries for all somatic sweeps of this cell).
        soma_electrode = _build_electrode(nwbfile, cell_index, "soma", None, amp)
        soma_result = _build_concatenated_series(
            nwbfile,
            bundle,
            cell_index,
            "soma_excitability",
            soma_electrode,
            session_start_time_naive,
        )
        if soma_result is not None:
            response, stim, sweep_slices = soma_result
            sorted_soma = sorted(bundle.sweeps_for(cell_index, "soma_excitability"), key=lambda s: s.pv_scan_date)
            for sw, (idx_start, count, sweep_start_time, protocol) in zip(sorted_soma, sweep_slices):
                row_idx = nwbfile.add_intracellular_recording(
                    electrode=soma_electrode,
                    stimulus=stim,
                    stimulus_start_index=idx_start,
                    stimulus_index_count=count,
                    response=response,
                    response_start_index=idx_start,
                    response_index_count=count,
                    cell_id=cell_index,
                    dendrite_type="Soma",
                    dendrite_distance_um=0.0,
                    sweep_start_time=sweep_start_time,
                    stimulus_current_pA=protocol.amplitude_pA if protocol else 0.0,
                    stimulus_protocol_name=protocol.name if protocol else "",
                    condition_subfolder=sw.condition_subfolder,
                )
                irt_rows_per_sequential[(cell_index, "soma_fi_curve", sw.condition_subfolder)].append(row_idx)
                # Per-row metadata for the RecordingsIndexTable (built post-chain).
                index_row_metadata[row_idx] = {
                    "cell_id": cell_index,
                    "condition_subfolder": sw.condition_subfolder,
                    "dendrite_type": "Soma",
                    "dendrite_distance_um": 0.0,
                    "stimulus_current_pA": protocol.amplitude_pA if protocol else 0.0,
                    "sweep_start_time_s": float(sweep_start_time),
                }

        # Dendritic block: group sweeps by (compartment, location_number) so each
        # location becomes one SequentialRecordings row.
        dend_sweeps = bundle.sweeps_for(cell_index, "dend_excitability")
        if not dend_sweeps:
            continue
        # Build per-location series. To preserve Pattern B (one series per
        # (cell, modality)), concatenate ALL dendritic sweeps but tag the IRT row
        # with the per-recording compartment.
        dend_subbundle = MouseDayBundle(
            subject_entry=bundle.subject_entry,
            cells={cell_index: dend_sweeps},
        )
        # Need an electrode per (location_type, location_number); create them all up front.
        location_electrodes: dict[tuple[str, int], IntracellularElectrode] = {}
        for sw in dend_sweeps:
            key = (sw.dendrite_location, sw.location_number)
            if key not in location_electrodes:
                location_electrodes[key] = _build_electrode(
                    nwbfile, cell_index, sw.dendrite_location, sw.location_number, amp
                )

        # Single concatenated dendritic series for the whole cell.
        dend_result = _build_concatenated_series(
            nwbfile,
            dend_subbundle,
            cell_index,
            "dend_excitability",
            # The "electrode" arg on the series is required but mostly informational.
            # Use the first electrode encountered; the IRT row's electrode column is
            # the load-bearing one for per-sweep electrode identity.
            next(iter(location_electrodes.values())),
            session_start_time_naive,
        )
        if dend_result is None:
            continue
        dend_response, dend_stim, dend_sweep_slices = dend_result
        # Add IRT rows with correct per-sweep electrode (from location).
        # The sweep order in dend_sweep_slices matches the chronological order used
        # in _build_concatenated_series (which sorts by pv_scan_date).
        sorted_sweeps = sorted(dend_sweeps, key=lambda s: s.pv_scan_date)
        for sw, (idx_start, count, sweep_start_time, protocol) in zip(sorted_sweeps, dend_sweep_slices):
            electrode = location_electrodes[(sw.dendrite_location, sw.location_number)]
            dendrite_type = "Proximal" if sw.dendrite_location == "prox" else "Distal"
            dendrite_distance = 40.0 if sw.dendrite_location == "prox" else 90.0

            row_idx = nwbfile.add_intracellular_recording(
                electrode=electrode,
                stimulus=dend_stim,
                stimulus_start_index=idx_start,
                stimulus_index_count=count,
                response=dend_response,
                response_start_index=idx_start,
                response_index_count=count,
                cell_id=cell_index,
                dendrite_type=dendrite_type,
                dendrite_distance_um=dendrite_distance,
                sweep_start_time=sweep_start_time,
                stimulus_current_pA=protocol.amplitude_pA if protocol else 0.0,
                stimulus_protocol_name=protocol.name if protocol else "",
                condition_subfolder=sw.condition_subfolder,
            )
            protocol_label = f"dend_{sw.dendrite_location}{sw.location_number}"
            irt_rows_per_sequential[(cell_index, protocol_label, sw.condition_subfolder)].append(row_idx)
            index_row_metadata[row_idx] = {
                "cell_id": cell_index,
                "condition_subfolder": sw.condition_subfolder,
                "dendrite_type": dendrite_type,
                "dendrite_distance_um": dendrite_distance,
                "stimulus_current_pA": protocol.amplitude_pA if protocol else 0.0,
                "sweep_start_time_s": float(sweep_start_time),
            }

            # Buffer line-scan kymographs for later Pattern B consolidation.
            # Each sweep contributes its kymo to the (cell, location, channel) buffer;
            # we'll concatenate along axis 0 and write one TimeSeries per buffer.
            kymos = _line_scan_kymographs(sw.recording_path)
            for chan_name, (data, rate) in kymos.items():
                buffer_key = (cell_index, sw.dendrite_location, sw.location_number, chan_name)
                line_scan_buffer[buffer_key].append((row_idx, data, rate, float(sweep_start_time)))

    # ---- Pattern B for line scans -----------------------------------------------
    # Consolidate buffered per-sweep kymographs into one TimeSeries per
    # (cell, location, channel). Each contributing sweep gets a slice (idx_start, count)
    # into the consolidated array; we record these for the RecordingsIndexTable.
    # Per-sweep slice info: irt_row_idx -> channel_name -> (timeseries, idx_start, count)
    line_scan_slices: dict[int, dict[str, tuple[TimeSeries, int, int]]] = defaultdict(dict)
    for buffer_key, entries in line_scan_buffer.items():
        cell_index, loc_type, loc_num, chan_name = buffer_key
        # Sort entries chronologically so concat order matches wall-clock.
        entries.sort(key=lambda e: e[3])  # sort by sweep_start_time_s
        arrays = [e[1] for e in entries]
        # Sanity: all kymos for one (cell, location, channel) should share the
        # same pixel-along-line dimension. Drop any with mismatched widths.
        ref_width = arrays[0].shape[1] if arrays[0].ndim == 2 else None
        clean: list[tuple[int, np.ndarray, float | None, float]] = []
        for e in entries:
            arr = e[1]
            if arr.ndim != 2 or (ref_width is not None and arr.shape[1] != ref_width):
                continue
            clean.append(e)
        if not clean:
            continue
        concat = np.concatenate([e[1] for e in clean], axis=0)
        rate = clean[0][2] if clean[0][2] is not None else 1.0
        first_start = clean[0][3]
        dendrite_type_pretty = "Proximal" if loc_type == "prox" else "Distal"
        ts_name = f"LineScanCell{cell_index}_{dendrite_type_pretty}{loc_num}_{chan_name}"
        ts_description = (
            f"Bruker line-scan kymograph for Cell {cell_index} {dendrite_type_pretty} location "
            f"{loc_num}, channel {chan_name}. All trio sweeps at this location concatenated along "
            f"the time axis (rows = line repetitions, columns = pixels along the scan line). Per-sweep "
            f"slices are exposed via processing/recordings_index/RecordingsIndexTable (line_scan_{chan_name.lower()} "
            f"column). starting_time is the wall-clock of the first contributing sweep relative to session_start_time."
        )
        ts = TimeSeries(
            name=ts_name,
            description=ts_description,
            data=concat,
            unit="a.u.",
            starting_time=float(first_start),
            rate=float(rate),
        )
        nwbfile.add_acquisition(ts)
        cumulative = 0
        for irt_idx, data, _rate, _start in clean:
            n_lines = data.shape[0]
            line_scan_slices[irt_idx][chan_name] = (ts, cumulative, n_lines)
            cumulative += n_lines

    # ---- icephys hierarchical chain --------------------------------------------
    # Build the chain on top of the IRT rows.
    # 1. SimultaneousRecordings: one row per IRT row (single-electrode rig).
    simultaneous_indices_per_irt: dict[int, int] = {}
    for irt_indices in irt_rows_per_sequential.values():
        for irt_idx in irt_indices:
            sim_idx = nwbfile.add_icephys_simultaneous_recording(recordings=[irt_idx])
            simultaneous_indices_per_irt[irt_idx] = sim_idx

    # 2. SequentialRecordings: one per (cell, protocol family, condition).
    # Paired-recording cells with both baseline and post-drug sweeps will get TWO
    # sequential rows for their soma_fi_curve protocol (one per condition).
    sequential_indices_per_key: dict[tuple[int, str, str], int] = {}
    sequential_index_per_irt: dict[int, int] = {}
    for key, irt_indices in irt_rows_per_sequential.items():
        sim_indices = [simultaneous_indices_per_irt[i] for i in irt_indices]
        cell_id, protocol, condition_sub = key
        if "soma" in protocol:
            stimulus_type = "F-I_curve"
        else:
            stimulus_type = "dendritic_trio"
        seq_idx = nwbfile.add_icephys_sequential_recording(
            simultaneous_recordings=sim_indices,
            stimulus_type=stimulus_type,
        )
        sequential_indices_per_key[key] = seq_idx
        for irt_idx in irt_indices:
            sequential_index_per_irt[irt_idx] = seq_idx

    # 3. Repetitions: one per (cell, protocol family, condition).
    repetitions_table = nwbfile.get_icephys_repetitions()
    repetitions_table.add_column(
        name="cell_id",
        description="Cell index within the mouse-day this repetition belongs to.",
    )
    repetitions_table.add_column(
        name="protocol",
        description="Protocol family identifier: soma_fi_curve | dend_proxN | dend_distN.",
    )
    repetitions_table.add_column(
        name="condition_subfolder",
        description=(
            "Raw lab condition subfolder this repetition belongs to. For paired-recording "
            "cells, the same protocol can have multiple repetitions across conditions (baseline + post-drug)."
        ),
    )
    # Track which repetitions belong to which condition for the ExperimentalConditions table.
    repetitions_per_condition: dict[str, list[int]] = defaultdict(list)
    repetition_index_per_irt: dict[int, int] = {}
    for key, seq_idx in sequential_indices_per_key.items():
        cell_id, protocol, condition_sub = key
        rep_idx = nwbfile.add_icephys_repetition(
            sequential_recordings=[seq_idx],
            cell_id=cell_id,
            protocol=protocol,
            condition_subfolder=condition_sub,
        )
        repetitions_per_condition[condition_sub].append(rep_idx)
        for irt_idx in irt_rows_per_sequential[key]:
            repetition_index_per_irt[irt_idx] = rep_idx

    # 4. ExperimentalConditions: one row per unique condition across the mouse-day.
    # For non-paired mouse-days this collapses to a single row. For paired (e.g., ON+SCH)
    # the table has two rows, each pointing to the repetitions acquired under that condition.
    ec_table = nwbfile.get_icephys_experimental_conditions()
    if "condition" not in ec_table.colnames:
        ec_table.add_column(
            name="condition",
            description="Paper-prose label for the experimental condition.",
        )
    ec_index_per_condition: dict[str, int] = {}
    for ec_idx, (condition_sub, rep_indices) in enumerate(repetitions_per_condition.items()):
        nwbfile.add_icephys_experimental_condition(
            repetitions=rep_indices,
            condition=_condition_label(condition_sub),
        )
        ec_index_per_condition[condition_sub] = ec_idx

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
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with NWBHDF5IO(output_path, "w") as io:
        io.write(nwbfile)
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
) -> None:
    """Build and attach the RecordingsIndexTable.

    The table denormalizes one row per IRT row with the cross-references from
    the icephys hierarchical chain plus paired line-scan slice references.
    For sweeps without line scans (somatic), line_scan_ch1/ch2 entries use
    idx_start=-1, count=-1 (null-reference convention).
    """
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
    stim_pA_col: list[float] = []
    sweep_start_col: list[float] = []
    irt_indices_col: list[int] = []
    sr_indices_col: list[int] = []
    seq_indices_col: list[int] = []
    rep_indices_col: list[int] = []
    ec_indices_col: list[int] = []
    line_scan_ch1_refs: list[TimeSeriesReference] = []
    line_scan_ch2_refs: list[TimeSeriesReference] = []

    for irt_idx in range(n_rows):
        meta = index_row_metadata.get(irt_idx, {})
        cell_id_col.append(int(meta.get("cell_id", -1)))
        condition_col.append(str(meta.get("condition_subfolder", "")))
        dendrite_type_col.append(str(meta.get("dendrite_type", "")))
        dendrite_distance_col.append(float(meta.get("dendrite_distance_um", 0.0)))
        stim_pA_col.append(float(meta.get("stimulus_current_pA", 0.0)))
        sweep_start_col.append(float(meta.get("sweep_start_time_s", 0.0)))
        irt_indices_col.append(irt_idx)
        sr_indices_col.append(int(simultaneous_index_per_irt.get(irt_idx, -1)))
        seq_indices_col.append(int(sequential_index_per_irt.get(irt_idx, -1)))
        rep_indices_col.append(int(repetition_index_per_irt.get(irt_idx, -1)))
        condition_sub = meta.get("condition_subfolder", "")
        ec_indices_col.append(int(ec_index_per_condition.get(condition_sub, -1)))

        slices = line_scan_slices.get(irt_idx, {})
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

    # Build the DynamicTable. Need to construct columns before passing to
    # DynamicTable so DynamicTableRegion can reference the icephys tables.
    columns: list = [
        VectorData(name="cell_id", description="Cell index within the mouse-day.", data=cell_id_col),
        VectorData(
            name="condition_subfolder",
            description="Raw lab condition subfolder for this sweep (e.g., 'LID off-state').",
            data=condition_col,
        ),
        VectorData(
            name="dendrite_type",
            description="Compartment of the recording electrode: Soma, Proximal, or Distal.",
            data=dendrite_type_col,
        ),
        VectorData(
            name="dendrite_distance_um",
            description="Approximate distance from soma in micrometers (0 for Soma, 40 for Proximal, 90 for Distal).",
            data=dendrite_distance_col,
        ),
        VectorData(
            name="stimulus_current_pA",
            description="Commanded current-step amplitude in picoamperes for this sweep.",
            data=stim_pA_col,
        ),
        VectorData(
            name="sweep_start_time_s",
            description="Wall-clock start time of this sweep in seconds, relative to session_start_time.",
            data=sweep_start_col,
        ),
        DynamicTableRegion(
            name="irt_row",
            description="Reference to the corresponding row in IntracellularRecordingsTable.",
            data=irt_indices_col,
            table=irt_table,
        ),
        DynamicTableRegion(
            name="simultaneous_row",
            description="Reference to the corresponding row in icephys_simultaneous_recordings.",
            data=sr_indices_col,
            table=sr_table,
        ),
        DynamicTableRegion(
            name="sequential_row",
            description="Reference to the corresponding row in icephys_sequential_recordings.",
            data=seq_indices_col,
            table=seq_table,
        ),
        DynamicTableRegion(
            name="repetition_row",
            description="Reference to the corresponding row in icephys_repetitions.",
            data=rep_indices_col,
            table=reps_table,
        ),
        DynamicTableRegion(
            name="experimental_condition_row",
            description="Reference to the corresponding row in icephys_experimental_conditions.",
            data=ec_indices_col,
            table=ec_table,
        ),
    ]

    # Line-scan columns only if at least one sweep has line-scan slices.
    if fallback_ts is not None:
        columns.append(
            TimeSeriesReferenceVectorData(
                name="line_scan_ch1",
                description=(
                    "Reference to the paired Alexa Fluor 568 reference channel slice (R0 for ratiometric "
                    "Ca2+ calculations). idx_start=-1 for sweeps without line-scan imaging (somatic)."
                ),
                data=line_scan_ch1_refs,
            )
        )
        columns.append(
            TimeSeriesReferenceVectorData(
                name="line_scan_ch2",
                description=(
                    "Reference to the paired Fluo-4 calcium signal channel slice (G in (G-G0)/G0R0). "
                    "idx_start=-1 for sweeps without line-scan imaging (somatic)."
                ),
                data=line_scan_ch2_refs,
            )
        )

    table = DynamicTable(
        name="RecordingsIndexTable",
        description=(
            "Denormalized lookup for fast analyst queries. One row per icephys sweep, mirroring the order "
            "of IntracellularRecordingsTable. Carries per-sweep metadata (cell_id, condition_subfolder, "
            "dendrite_type, dendrite_distance_um, stimulus_current_pA, sweep_start_time_s) plus "
            "DynamicTableRegion pointers to the full icephys hierarchical chain "
            "(irt_row, simultaneous_row, sequential_row, repetition_row, experimental_condition_row) "
            "and TimeSeriesReference slices into the paired line-scan TimeSeries (line_scan_ch1, line_scan_ch2). "
            "This is a query helper; the canonical icephys chain is still authoritative."
        ),
        columns=columns,
    )

    module = nwbfile.create_processing_module(
        name="recordings_index",
        description=(
            "Cross-modality and cross-table lookup helpers for the per-mouse-day icephys data. "
            "Contains the RecordingsIndexTable that denormalizes the canonical icephys chain "
            "for one-stop pandas-style filtering."
        ),
    )
    module.add(table)


# Library module: `merge_one_mouse_day` is imported by the production runner
# (`scripts/full_per_mouse_day_merge.py`). No `main()` here — running the
# production runner is the right entry point for both single mouse-days and
# the whole dataset.
