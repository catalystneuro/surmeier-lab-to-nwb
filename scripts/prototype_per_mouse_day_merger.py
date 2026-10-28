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
from pynwb import NWBFile
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


# --------- reading per-sweep line-scan field-of-view + pixel-mask + ROI ---------------------


def _read_dendritic_artifacts(sweep_path: Path, channel: str) -> dict:
    """Read source image, pixel mask, and ROI response for one dendritic sweep + channel.

    Delegates to the project's existing `BrukerLineScanSegmentationExtractor` for
    line endpoints, pixel mask, ROI fluorescence profile, and frame shape — keeping
    a single source of truth for the Bruker XML / CSV parsing. The source TIFF is
    resolved from the same XML the extractor already loaded.

    Returns a dict with keys:
        source_image: ndarray (H, W) grayscale
        pixel_mask_image: ndarray (H, W) bool — True at line-scan pixel locations
        frame_shape: (H, W)
        roi_response: 1D ndarray (samples,) — fluorescence-along-line vs time
        rate: float Hz
        pixel_xs / pixel_ys: 1D int arrays of line pixel coords
    Raises if the sweep folder can't be parsed or any of these is missing.
    """
    import tifffile

    from surmeier_lab_to_nwb.zhai2025.extractors.bruker_line_scan_segmentation_extractor import (
        BrukerLineScanSegmentationExtractor,
    )

    extractor = BrukerLineScanSegmentationExtractor(folder_path=sweep_path, channel_name=channel)

    frame_h, frame_w = extractor.get_frame_shape()
    pixel_mask_image = np.squeeze(extractor._roi_masks.data).astype(bool)

    # Recover line endpoints from the extractor for the line-pixel array.
    sx, sy = extractor._line_start_x, extractor._line_start_y
    ex, ey = extractor._line_stop_x, extractor._line_stop_y
    if sy == ey:
        xs = np.arange(min(sx, ex), max(sx, ex) + 1)
        ys = np.full_like(xs, sy)
    elif sx == ex:
        ys = np.arange(min(sy, ey), max(sy, ey) + 1)
        xs = np.full_like(ys, sx)
    else:
        n_steps = max(abs(ex - sx), abs(ey - sy)) + 1
        xs = np.round(np.linspace(sx, ex, n_steps)).astype(int)
        ys = np.round(np.linspace(sy, ey, n_steps)).astype(int)
    xs = np.clip(xs, 0, frame_w - 1)
    ys = np.clip(ys, 0, frame_h - 1)

    # ROI response trace + sample rate from the extractor.
    if not extractor._roi_responses:
        raise ValueError(f"No ROI response for {channel} in {sweep_path}")
    traces = extractor._roi_responses[0].data
    roi_response = np.squeeze(np.asarray(traces, dtype=np.float32))
    if extractor._sampling_frequency:
        roi_rate = float(extractor._sampling_frequency)
    elif len(extractor._timestamps) > 1:
        roi_rate = 1.0 / float(np.mean(np.diff(extractor._timestamps)))
    else:
        raise ValueError(f"No sampling rate for the {channel} ROI response in {sweep_path}")

    # Source TIFF: resolve filename from the already-parsed XML on the extractor.
    sequence = extractor._xml_metadata["PVScan"]["Sequence"]
    if isinstance(sequence, list):
        sequence = next((s for s in sequence if int(s.get("@cycle", 0)) == extractor.cycle), sequence[0])
    frame_meta = sequence["Frame"]
    if isinstance(frame_meta, list):
        frame_meta = frame_meta[0]
    file_entries = frame_meta["File"]
    if not isinstance(file_entries, list):
        file_entries = [file_entries]
    source_filename = next(
        (e["@source"] for e in file_entries if e.get("@channelName") == channel and "@source" in e),
        None,
    )
    if source_filename is None:
        raise ValueError(f"No source image listed for {channel} in {sweep_path}")
    arr = tifffile.imread(sweep_path / source_filename)
    if arr.ndim == 3:
        arr = (
            np.mean(arr[:, :, :3], axis=2).astype(np.float32) if arr.shape[2] >= 3 else arr[:, :, 0].astype(np.float32)
        )
    source_image = np.asarray(arr, dtype=np.float32)

    return {
        "source_image": source_image,
        "pixel_mask_image": pixel_mask_image,
        "frame_shape": (frame_h, frame_w),
        "roi_response": roi_response,
        "rate": roi_rate,
        "pixel_xs": xs.astype(int),
        "pixel_ys": ys.astype(int),
    }


# --------- reading per-sweep line-scan kymographs -------------------------------------------


def _line_scan_kymographs(sweep_path: Path) -> dict[str, tuple[np.ndarray, float]]:
    """Read line-scan kymograph TIFFs for one dendritic sweep.

    Returns {channel_name: (data, rate_hz)}. Channels are typically:
        Ch1 = Alexa Fluor 568 reference (red) — R0 for ratiometric calculations
        Ch2 = Fluo-4 calcium indicator (green) — G in (G-G0)/G0R0
    """
    import tifffile

    out: dict[str, tuple[np.ndarray, float]] = {}
    for channel in ("Ch1", "Ch2"):
        candidates = list(sweep_path.glob(f"*Cycle00001_{channel}_*.ome.tif"))
        if not candidates:
            raise FileNotFoundError(f"No {channel} line-scan TIFF in {sweep_path}")
        tif = candidates[0]
        arr = tifffile.imread(tif)
        # Some recordings package both channels in one 3D TIFF (channels, lines, pixels).
        if arr.ndim == 3:
            channel_index = {"Ch1": 0, "Ch2": 1}[channel]
            if channel_index >= arr.shape[0]:
                raise ValueError(f"{tif} has {arr.shape[0]} channels, expected {channel}")
            arr = arr[channel_index]
        if arr.ndim != 2:
            raise ValueError(f"{tif} has {arr.ndim} dimensions, expected 2 (lines, pixels)")
        # Sample rate: from PrairieView's scanLinePeriod in the master XML.
        rate = _read_scan_line_rate(sweep_path)
        if rate is None:
            raise ValueError(f"No scanLinePeriod in the XML of {sweep_path}")
        out[channel] = (arr, rate)
    return out


# Ch1 and Ch2 are the microscope's detector numbers. Objects in the file are named after the dye
# instead; only Ch2 rises after the spikes, which confirms this assignment.
_DYE_NAME = {"Ch1": "Alexa568", "Ch2": "Fluo4"}
_DYE_DESCRIPTION = {
    "Ch1": "Alexa Fluor 568 (structural dye, red detector 580-620 nm)",
    "Ch2": "Fluo-4 (calcium indicator, green detector 490-560 nm)",
}


# What one dendritic trial is (lab protocol step 23 and the paper's methods); the lab calls it a "trio".
_DENDRITIC_TRIAL = (
    " In each trial, three 2 nA current injections of 2 ms at 50 Hz at the cell body make the cell fire spikes "
    "that travel back into the dendrites."
)


def _spot_sentence(spot: str) -> str:
    """Description sentence for a site the lab imaged at two spots of the dendrite; empty for one spot."""
    if not spot:
        return ""
    sentence = (
        " The lab imaged this site at two spots on the dendrite, moving the microscope in between; "
        f"this series holds the trials at spot {spot[-1]}"
    )
    return sentence + (", which the lab's files mark 'real'." if spot == "Spot2" else ".")


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


def _series_name(cell_index: int, modality: str, neurodata_type: str = "CurrentClampSeries") -> str:
    suffix = "Somatic" if modality == "soma_excitability" else "Dendritic"
    return f"{neurodata_type}Cell{cell_index}{suffix}"


def _protocol_description(protocol_label: str) -> str:
    """Readable protocol for the repetitions table from the internal label ('soma_fi_curve', 'dend_prox1')."""
    if protocol_label == "soma_fi_curve":
        return "F-I curve"
    site = protocol_label.removeprefix("dend_")
    dendrite_type = "Proximal" if site.startswith("prox") else "Distal"
    return f"back-propagating spikes at {dendrite_type}{site[4:]}"


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
            description="Bruker Ultima In Vitro multiphoton microscope.",
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
    if modality == "soma_excitability":
        recording = "while current steps were injected at the cell body to measure its firing (F-I curve)."
    else:
        recording = "during dendritic imaging." + _DENDRITIC_TRIAL
    response = CurrentClampSeries(
        name=series_name,
        data=data_concat,
        electrode=electrode,
        gain=1.0,
        rate=float(sample_rate),
        starting_time=sweep_slices[0][2],
        description=(
            f"Current-clamp recordings from cell {cell_index}, {recording} The {len(sweep_slices)} sweeps are "
            "stored back to back; each sweep's start time is in the sweep_start_time column of the intracellular "
            "recordings table."
        ),
    )
    nwbfile.add_acquisition(response)

    stim = CurrentClampStimulusSeries(
        name=_series_name(cell_index, modality, "CurrentClampStimulusSeries"),
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
    # Maps irt_row_index -> dict(cell_id, condition_subfolder, dendrite_type, dendrite_distance_um,
    #                          site, stimulus_current_pA, sweep_start_time_s).
    index_row_metadata: dict[int, dict] = {}

    # The scan line at a location sometimes moves or changes length between trials.
    # scan_lines_per_location[(cell_id, location_type, location_number, spot)] maps each distinct
    # line (its pixel coordinates) to its number at that location, in order of first use.
    # spot is "" for a site imaged at one spot, "Spot1" or "Spot2" for a site imaged at two.
    scan_lines_per_location: dict[tuple[int, str, int, str], dict[tuple, int]] = defaultdict(dict)

    # Buffer for per-sweep line-scan kymographs to consolidate later.
    # Key: (cell_id, location_type, location_number, spot, scan_line_number, channel_name)
    # Value: list of (irt_row_index, kymo_data_2d, rate_hz, sweep_start_time_s)
    # After all dendritic IRT rows are added, we concatenate each buffer along axis 0
    # (the time/lines axis) to make one TimeSeries per (cell, location, scan line, channel) -- Pattern B
    # for line scans. Per-sweep slice info is preserved for the RecordingsIndexTable.
    line_scan_buffer: dict[tuple[int, str, int, str, int, str], list[tuple[int, np.ndarray, float, float]]] = (
        defaultdict(list)
    )

    # Parallel buffers for the rich line-scan objects:
    # - ROI responses: Pattern B-consolidated per (cell, location, scan line, channel) — time-series data
    #   benefits from continuous concatenation across trios.
    # - Source images: kept per-trial (one image per sweep × channel) so cross-trial FOV
    #   verification is possible. Static reference data, cheap to keep separately.
    # - Pixel masks: per (cell, location, scan line, channel).
    roi_response_buffer: dict[tuple[int, str, int, str, int, str], list[tuple[int, np.ndarray, float, float]]] = (
        defaultdict(list)
    )
    # source_image_buffer[(cell, loc, num, spot, chan)] -> list of (irt_row_index, sweep_start, image_array)
    source_image_buffer: dict[tuple[int, str, int, str, str], list[tuple[int, float, np.ndarray]]] = defaultdict(list)
    pixel_mask_buffer: dict[tuple[int, str, int, str, int, str], dict] = {}  # mask_image, xs, ys, frame_shape

    # Add custom columns to IntracellularRecordingsTable BEFORE adding rows.
    # The `electrodes` aligned sub-table doesn't exist until the first add; using
    # the convenience add_intracellular_recording call below pre-creates it.
    irt = nwbfile.get_intracellular_recordings()
    irt.add_column(
        name="cell_id", description="Cell number within this file (the cells recorded from this mouse on this day)."
    )
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
    irt.add_column(
        name="site",
        description=(
            "Recording site, named like the line-scan series: 'Soma', or 'Proximal'/'Distal' with the site number "
            "(e.g., 'Proximal1'). Sites the lab imaged at two spots on the dendrite, moving the microscope in "
            "between, add 'Spot1' or 'Spot2' (e.g., 'Proximal1Spot2'); the lab's files mark the second spot 'real'."
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
                row_index = nwbfile.add_intracellular_recording(
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
                    site="Soma",
                )
                irt_rows_per_sequential[(cell_index, "soma_fi_curve", sw.condition_subfolder)].append(row_index)
                # Per-row metadata for the RecordingsIndexTable (built post-chain).
                index_row_metadata[row_index] = {
                    "cell_id": cell_index,
                    "condition_subfolder": sw.condition_subfolder,
                    "dendrite_type": "Soma",
                    "dendrite_distance_um": 0.0,
                    "site": "Soma",
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
        # Sites the lab imaged at two spots on the dendrite; the second spot's folders say "real".
        sites_with_two_spots = {
            (sw.dendrite_location, sw.location_number) for sw in dend_sweeps if "real" in sw.site_label
        }
        for sw, (idx_start, count, sweep_start_time, protocol) in zip(sorted_sweeps, dend_sweep_slices):
            electrode = location_electrodes[(sw.dendrite_location, sw.location_number)]
            dendrite_type = "Proximal" if sw.dendrite_location == "prox" else "Distal"
            dendrite_distance = 40.0 if sw.dendrite_location == "prox" else 90.0
            spot = ""
            if (sw.dendrite_location, sw.location_number) in sites_with_two_spots:
                spot = "Spot2" if "real" in sw.site_label else "Spot1"
            site = f"{dendrite_type}{sw.location_number}{spot}"

            row_index = nwbfile.add_intracellular_recording(
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
                site=site,
            )
            protocol_label = f"dend_{sw.dendrite_location}{sw.location_number}"
            irt_rows_per_sequential[(cell_index, protocol_label, sw.condition_subfolder)].append(row_index)
            index_row_metadata[row_index] = {
                "cell_id": cell_index,
                "condition_subfolder": sw.condition_subfolder,
                "dendrite_type": dendrite_type,
                "dendrite_distance_um": dendrite_distance,
                "site": site,
                "stimulus_current_pA": protocol.amplitude_pA if protocol else 0.0,
                "sweep_start_time_s": float(sweep_start_time),
            }

            # Buffer line-scan kymographs and rich artifacts (source image, pixel mask,
            # ROI response) for later Pattern B consolidation. Both channels share the
            # scan line, so its number at this location comes from the Ch1 pixel coordinates.
            # Each spot of a site imaged at two spots is a location of its own.
            kymos = _line_scan_kymographs(sw.recording_path)
            artifacts = {chan_name: _read_dendritic_artifacts(sw.recording_path, chan_name) for chan_name in kymos}
            location = (cell_index, sw.dendrite_location, sw.location_number, spot)
            scan_line = (tuple(artifacts["Ch1"]["pixel_xs"]), tuple(artifacts["Ch1"]["pixel_ys"]))
            scan_lines = scan_lines_per_location[location]
            scan_line_number = scan_lines.setdefault(scan_line, len(scan_lines) + 1)
            for chan_name, (data, rate) in kymos.items():
                buffer_key = (*location, scan_line_number, chan_name)
                line_scan_buffer[buffer_key].append((row_index, data, rate, float(sweep_start_time)))

            # Source images are kept per trial; pixel masks are one per scan line; ROI
            # responses get Pattern B-consolidated across the trials of each scan line.
            for chan_name, art in artifacts.items():
                source_image_buffer[(*location, chan_name)].append(
                    (row_index, float(sweep_start_time), art["source_image"])
                )
                buffer_key = (*location, scan_line_number, chan_name)
                if buffer_key not in pixel_mask_buffer:
                    pixel_mask_buffer[buffer_key] = {
                        "mask_image": art["pixel_mask_image"],
                        "xs": art["pixel_xs"],
                        "ys": art["pixel_ys"],
                        "frame_shape": art["frame_shape"],
                    }
                roi_response_buffer[buffer_key].append(
                    (row_index, art["roi_response"], art["rate"], float(sweep_start_time))
                )

    # ---- Pattern B for line scans -----------------------------------------------
    # Consolidate buffered per-sweep kymographs into one TimeSeries per
    # (cell, location, channel). Each contributing sweep gets a slice (idx_start, count)
    # into the consolidated array; we record these for the RecordingsIndexTable.
    # Per-sweep slice info: irt_row_index -> channel_name -> (timeseries, idx_start, count)
    line_scan_slices: dict[int, dict[str, tuple[TimeSeries, int, int]]] = defaultdict(dict)
    for buffer_key, entries in line_scan_buffer.items():
        cell_index, loc_type, loc_num, spot, scan_line_number, chan_name = buffer_key
        # Sort entries chronologically so concat order matches wall-clock.
        entries.sort(key=lambda e: e[3])  # sort by sweep_start_time_s
        concat = np.concatenate([e[1] for e in entries], axis=0)
        rate = entries[0][2]
        first_start = entries[0][3]
        dendrite_type_pretty = "Proximal" if loc_type == "prox" else "Distal"
        dye = _DYE_NAME[chan_name]
        ts_name = f"LineScan{dye}Cell{cell_index}{dendrite_type_pretty}{loc_num}{spot}"
        ts_description = (
            f"Bruker line-scan kymograph for Cell {cell_index} {dendrite_type_pretty} location "
            f"{loc_num}, {_DYE_DESCRIPTION[chan_name]}. All trials at this site, concatenated in recording order "
            f"along the time axis (rows = line repetitions, columns = pixels along the scan line); each trial's slice "
            f"is in the line_scan_{dye.lower()} column of processing/recordings_index/RecordingsIndexTable. "
            f"starting_time is the start of the first trial.{_DENDRITIC_TRIAL}"
        )
        ts_description += _spot_sentence(spot)
        n_scan_lines = len(scan_lines_per_location[(cell_index, loc_type, loc_num, spot)])
        if n_scan_lines > 1:
            ts_name = f"{ts_name}Line{scan_line_number}"
            ts_description += (
                f" The scan line at this location was moved or changed length between trials, so its trials "
                f"are split into {n_scan_lines} series; this one holds the trials recorded with line "
                f"{scan_line_number} ({concat.shape[1]} pixels)."
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
        for irt_index, data, _rate, _start in entries:
            n_lines = data.shape[0]
            line_scan_slices[irt_index][chan_name] = (ts, cumulative, n_lines)
            cumulative += n_lines

    # ---- Rich line-scan objects: source images + pixel masks + ROI responses --------
    # Source images go in `acquisition/ImageLineScanSource` (Images container).
    # Pixel masks go in `processing/ophys/ImageSegmentation` (PlaneSegmentation).
    # ROI responses (consolidated Pattern B) go in `processing/ophys/Fluorescence`.
    roi_response_slices: dict[int, dict[str, tuple]] = defaultdict(dict)
    if source_image_buffer or pixel_mask_buffer or roi_response_buffer:
        _build_line_scan_rich_objects(
            nwbfile,
            source_image_buffer=source_image_buffer,
            pixel_mask_buffer=pixel_mask_buffer,
            roi_response_buffer=roi_response_buffer,
            roi_response_slices=roi_response_slices,
            scan_lines_per_location=scan_lines_per_location,
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
        cell_id, protocol, condition_sub = key
        if "soma" in protocol:
            stimulus_type = "F-I curve"
        else:
            stimulus_type = "back-propagating spikes"
        seq_index = nwbfile.add_icephys_sequential_recording(
            simultaneous_recordings=sim_indices,
            stimulus_type=stimulus_type,
        )
        sequential_indices_per_key[key] = seq_index
        for irt_index in irt_indices:
            sequential_index_per_irt[irt_index] = seq_index

    # 3. Repetitions: one per (cell, protocol family, condition).
    repetitions_table = nwbfile.get_icephys_repetitions()
    repetitions_table.add_column(
        name="cell_id",
        description="Cell number within this file that this repetition belongs to.",
    )
    repetitions_table.add_column(
        name="protocol",
        description=(
            "Protocol of this repetition: 'F-I curve' (current steps at the cell body) or 'back-propagating spikes' "
            "followed by the dendritic site imaged (e.g., 'back-propagating spikes at Proximal1')."
        ),
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
    for key, seq_index in sequential_indices_per_key.items():
        cell_id, protocol, condition_sub = key
        rep_index = nwbfile.add_icephys_repetition(
            sequential_recordings=[seq_index],
            cell_id=cell_id,
            protocol=_protocol_description(protocol),
            condition_subfolder=condition_sub,
        )
        repetitions_per_condition[condition_sub].append(rep_index)
        for irt_index in irt_rows_per_sequential[key]:
            repetition_index_per_irt[irt_index] = rep_index

    # 4. ExperimentalConditions: one row per unique condition across the mouse-day.
    # For non-paired mouse-days this collapses to a single row. For paired (e.g., ON+SCH)
    # the table has two rows, each pointing to the repetitions acquired under that condition.
    ec_table = nwbfile.get_icephys_experimental_conditions()
    if "condition" not in ec_table.colnames:
        ec_table.add_column(
            name="condition",
            description="Experimental condition as named in the paper (e.g., 'LID off-state').",
        )
    ec_index_per_condition: dict[str, int] = {}
    for ec_index, (condition_sub, rep_indices) in enumerate(repetitions_per_condition.items()):
        nwbfile.add_icephys_experimental_condition(
            repetitions=rep_indices,
            condition=_condition_label(condition_sub),
        )
        ec_index_per_condition[condition_sub] = ec_index

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
    # Use neuroconv's helper to auto-apply default HDF5 chunking + gzip compression
    # via get_default_backend_configuration (saves ~50% disk vs plain contiguous writes).
    from neuroconv.tools.nwb_helpers import configure_and_write_nwbfile

    configure_and_write_nwbfile(nwbfile=nwbfile, nwbfile_path=output_path, backend="hdf5")
    return output_path


def _build_line_scan_rich_objects(
    nwbfile,
    *,
    source_image_buffer: dict,
    pixel_mask_buffer: dict,
    roi_response_buffer: dict,
    roi_response_slices: dict,
    scan_lines_per_location: dict,
) -> None:
    """Add the line-scan source images, pixel-mask PlaneSegmentations, and Pattern B
    consolidated RoiResponseSeries to the NWB file.

    Parameters
    ----------
    source_image_buffer
        Map (cell, location_type, location_num, spot, channel) -> 2D ndarray (the field-of-view
        image with the scan line overlay drawn by Prairie View).
    pixel_mask_buffer
        Map (cell, location_type, location_num, spot, scan_line_number, channel) -> dict with
        'mask_image' (HxW bool), 'xs', 'ys', 'frame_shape'.
    roi_response_buffer
        Map same-key as pixel_mask_buffer -> list of (irt_index, 1D roi_response, rate_hz,
        sweep_start_s) tuples across trials with that scan line.
    roi_response_slices
        Output dict, irt_index -> {channel_name: (timeseries, idx_start, count)}, populated
        by this function so the caller can wire it into the RecordingsIndexTable.
    scan_lines_per_location
        Map (cell, location_type, location_num, spot) -> {scan line pixel coordinates: scan_line_number}.
    """
    from pynwb.device import Device as _Device
    from pynwb.image import Image, Images
    from pynwb.ophys import (
        Fluorescence,
        ImageSegmentation,
        OpticalChannel,
        RoiResponseSeries,
    )

    # Source images -> nwbfile.acquisition['ImageLineScanSource']
    # One image per TRIAL (per cell, channel, location, trio sweep). Multiple images
    # at the same (cell, channel, location) capture potential FOV drift across trios.
    if source_image_buffer:
        if "ImageLineScanSource" not in nwbfile.acquisition:
            container = Images(
                name="ImageLineScanSource",
                description=(
                    "Source images for line scan experiments: field of view snapshots with "
                    "scan-line overlay, used to identify the recorded region. One image per "
                    "trial and dye."
                ),
            )
            nwbfile.add_acquisition(container)
        else:
            container = nwbfile.acquisition["ImageLineScanSource"]
        for key in sorted(source_image_buffer.keys()):
            cell, loc, num, spot, chan = key
            entries = source_image_buffer[key]
            entries.sort(key=lambda e: e[1])  # by sweep_start_time
            dendrite_type_pretty = "Proximal" if loc == "prox" else "Distal"
            spot_note = f", spot {spot[-1]} of 2" if spot else ""
            indicator = _DYE_NAME[chan]
            for trial_index, (row_index, sweep_start, arr) in enumerate(entries, start=1):
                name = f"Image{indicator}Cell{cell}{dendrite_type_pretty}{num}{spot}Trial{trial_index}"
                container.add_image(
                    Image(
                        name=name,
                        description=(
                            f"Source image (field of view with scan line overlay) for Cell {cell} "
                            f"{dendrite_type_pretty} location {num}{spot_note}, {_DYE_DESCRIPTION[chan]}, trial {trial_index} "
                            f"of {len(entries)} in recording order "
                            f"(row {row_index} of the intracellular recordings table, start time {sweep_start:.1f} s)."
                        ),
                        data=arr,
                    )
                )

    # Need an ImagingPlane for the PlaneSegmentations. Use a single shared one.
    if pixel_mask_buffer or roi_response_buffer:
        ophys_module = nwbfile.processing.get("ophys")
        if ophys_module is None:
            ophys_module = nwbfile.create_processing_module(
                name="ophys",
                description=(
                    "Optical physiology data: line-scan ROI responses + pixel-mask segmentations "
                    "for dendritic excitability recordings."
                ),
            )

        # Find or create a device + imaging plane
        device_name = "BrukerFluorescenceMicroscope"
        if device_name in nwbfile.devices:
            device = nwbfile.devices[device_name]
        else:
            device = _Device(name=device_name, description="Bruker Prairie View two-photon microscope.")
            nwbfile.add_device(device)

        imaging_plane_name = "ImagingPlaneLineScan"
        if imaging_plane_name in nwbfile.imaging_planes:
            imaging_plane = nwbfile.imaging_planes[imaging_plane_name]
        else:
            optical_channel = OpticalChannel(
                name="OpticalChannelLineScan",
                description="Composite optical channel for two-channel line-scan recording.",
                emission_lambda=float("nan"),
            )
            imaging_plane = nwbfile.create_imaging_plane(
                name=imaging_plane_name,
                description=("Two-photon line-scan imaging plane in the dorsolateral striatum (Caudoputamen)."),
                device=device,
                excitation_lambda=810.0,
                indicator="Fluo-4 (Ca2+) + Alexa Fluor 568 (reference)",
                location="Caudoputamen",
                optical_channel=optical_channel,
            )

    # PlaneSegmentation with pixel masks: one row per (cell, channel, location, scan line)
    plane_seg_row_for_key: dict[tuple, int] = {}
    if pixel_mask_buffer:
        # Find or create ImageSegmentation in ophys
        if "ImageSegmentation" in ophys_module.data_interfaces:
            img_seg = ophys_module["ImageSegmentation"]
        else:
            img_seg = ImageSegmentation(name="ImageSegmentation")
            ophys_module.add(img_seg)

        ps_name = "PlaneSegmentationLineScan"
        if ps_name in img_seg.plane_segmentations:
            plane_seg = img_seg.plane_segmentations[ps_name]
        else:
            plane_seg = img_seg.create_plane_segmentation(
                name=ps_name,
                description=(
                    "Line-scan ROIs: each row is the scan line for one (cell, channel, "
                    "dendrite location, scan line). The pixel_mask gives the (x, y, weight) coords of "
                    "the line on the source image."
                ),
                imaging_plane=imaging_plane,
            )
        for key in sorted(pixel_mask_buffer.keys()):
            cell, loc, num, spot, scan_line_number, chan = key
            info = pixel_mask_buffer[key]
            xs = info["xs"]
            ys = info["ys"]
            # pixel_mask is a list of (x, y, weight) tuples
            pixel_mask = [(int(x), int(y), 1.0) for x, y in zip(xs, ys)]
            plane_seg_row_for_key[key] = len(plane_seg.id)
            plane_seg.add_roi(pixel_mask=pixel_mask)

    # Fluorescence with consolidated RoiResponseSeries — Pattern B per (cell, location, scan line, channel)
    if roi_response_buffer:
        if "Fluorescence" in ophys_module.data_interfaces:
            fluo = ophys_module["Fluorescence"]
        else:
            fluo = Fluorescence(name="Fluorescence")
            ophys_module.add(fluo)

        for key in sorted(roi_response_buffer.keys()):
            cell, loc, num, spot, scan_line_number, chan = key
            entries = roi_response_buffer[key]
            entries.sort(key=lambda e: e[3])
            arrays = [e[1] for e in entries]
            concat = np.concatenate(arrays, axis=0)
            rate = entries[0][2]
            first_start = entries[0][3]

            dendrite_type_pretty = "Proximal" if loc == "prox" else "Distal"
            indicator = _DYE_NAME[chan]
            rrs_name = f"RoiResponseSeries{indicator}Cell{cell}{dendrite_type_pretty}{num}{spot}"
            rrs_description = (
                f"Line-scan fluorescence profile (average across pixels along the scan line) "
                f"for Cell {cell} {dendrite_type_pretty} location {num}, {_DYE_DESCRIPTION[chan]}. "
                f"All trials at this site, concatenated in recording order; each trial's slice is in the "
                f"roi_response_{indicator.lower()} column of processing/recordings_index/RecordingsIndexTable."
                f"{_DENDRITIC_TRIAL}"
            )
            rrs_description += _spot_sentence(spot)
            n_scan_lines = len(scan_lines_per_location[(cell, loc, num, spot)])
            if n_scan_lines > 1:
                rrs_name = f"{rrs_name}Line{scan_line_number}"
                rrs_description += (
                    f" The scan line at this location was moved or changed length between trials, so its "
                    f"trials are split into {n_scan_lines} series; this one holds the trials recorded with "
                    f"line {scan_line_number}."
                )

            # Build a DynamicTableRegion pointing at this row of PlaneSegmentationLineScan
            rois_region = plane_seg.create_roi_table_region(
                description=f"Line-scan ROI for Cell {cell} {dendrite_type_pretty} location {num}",
                region=[plane_seg_row_for_key[key]],
            )
            rrs = RoiResponseSeries(
                name=rrs_name,
                description=rrs_description,
                data=concat.astype(np.float32),
                rois=rois_region,
                unit="a.u.",
                starting_time=float(first_start),
                rate=float(rate),
            )
            fluo.add_roi_response_series(rrs)

            cumulative = 0
            for irt_index, data, _rate, _start in entries:
                n = data.shape[0]
                roi_response_slices[irt_index][chan] = (rrs, cumulative, n)
                cumulative += n


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
        condition_col.append(str(meta.get("condition_subfolder", "")))
        dendrite_type_col.append(str(meta.get("dendrite_type", "")))
        dendrite_distance_col.append(float(meta.get("dendrite_distance_um", 0.0)))
        site_col.append(str(meta.get("site", "")))
        stim_pA_col.append(float(meta.get("stimulus_current_pA", 0.0)))
        sweep_start_col.append(float(meta.get("sweep_start_time_s", 0.0)))
        irt_indices_col.append(irt_index)
        sr_indices_col.append(int(simultaneous_index_per_irt.get(irt_index, -1)))
        seq_indices_col.append(int(sequential_index_per_irt.get(irt_index, -1)))
        rep_indices_col.append(int(repetition_index_per_irt.get(irt_index, -1)))
        condition_sub = meta.get("condition_subfolder", "")
        ec_indices_col.append(int(ec_index_per_condition.get(condition_sub, -1)))

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
    columns: list = [
        VectorData(name="cell_id", description="Cell number within this file.", data=cell_id_col),
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
            name="site",
            description="Recording site (e.g., 'Soma', 'Proximal1', 'Proximal1Spot2'); see the same column in "
            "IntracellularRecordingsTable.",
            data=site_col,
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
            description="This sweep's row in the intracellular recordings table.",
            data=irt_indices_col,
            table=irt_table,
        ),
        DynamicTableRegion(
            name="simultaneous_row",
            description="This sweep's row in the simultaneous recordings table.",
            data=sr_indices_col,
            table=sr_table,
        ),
        DynamicTableRegion(
            name="sequential_row",
            description="This sweep's row in the sequential recordings table.",
            data=seq_indices_col,
            table=seq_table,
        ),
        DynamicTableRegion(
            name="repetition_row",
            description="This sweep's row in the repetitions table.",
            data=rep_indices_col,
            table=reps_table,
        ),
        DynamicTableRegion(
            name="experimental_condition_row",
            description="This sweep's row in the experimental conditions table.",
            data=ec_indices_col,
            table=ec_table,
        ),
    ]

    # Line-scan columns only if at least one sweep has line-scan slices.
    if fallback_ts is not None:
        columns.append(
            TimeSeriesReferenceVectorData(
                name="line_scan_alexa568",
                description=(
                    "This sweep's slice of the Alexa Fluor 568 line scan (structural dye). Somatic sweeps have no "
                    "line scan and hold idx_start=-1."
                ),
                data=line_scan_ch1_refs,
            )
        )
        columns.append(
            TimeSeriesReferenceVectorData(
                name="line_scan_fluo4",
                description=(
                    "This sweep's slice of the Fluo-4 line scan (calcium indicator). Somatic sweeps have no line scan "
                    "and hold idx_start=-1."
                ),
                data=line_scan_ch2_refs,
            )
        )
    if fallback_rrs is not None:
        columns.append(
            TimeSeriesReferenceVectorData(
                name="roi_response_alexa568",
                description=(
                    "This sweep's slice of the Alexa Fluor 568 fluorescence (mean along the scan line). Somatic "
                    "sweeps have no imaging and hold idx_start=-1."
                ),
                data=roi_response_ch1_refs,
            )
        )
        columns.append(
            TimeSeriesReferenceVectorData(
                name="roi_response_fluo4",
                description=(
                    "This sweep's slice of the Fluo-4 fluorescence (mean along the scan line). Somatic sweeps have "
                    "no imaging and hold idx_start=-1."
                ),
                data=roi_response_ch2_refs,
            )
        )

    table = DynamicTable(
        name="RecordingsIndexTable",
        description=(
            "One row per sweep, in the order of the intracellular recordings table, so a sweep can be looked up "
            "in one place. Each row holds the sweep's details (cell_id, condition_subfolder, dendrite_type, "
            "dendrite_distance_um, site, stimulus_current_pA, sweep_start_time_s), its rows in the intracellular "
            "recordings, simultaneous recordings, sequential recordings, repetitions and experimental conditions "
            "tables, and, for dendritic sweeps, its slices of the line-scan and fluorescence series. The same "
            "information is also held in those tables and series."
        ),
        columns=columns,
    )

    module = nwbfile.create_processing_module(
        name="recordings_index",
        description=(
            "RecordingsIndexTable, which gathers for each sweep its details, its rows in the intracellular tables "
            "and its slices of the line-scan and fluorescence series."
        ),
    )
    module.add(table)


# Library module: `merge_one_mouse_day` is imported by the production runner
# (`scripts/full_per_mouse_day_merge.py`). No `main()` here — running the
# production runner is the right entry point for both single mouse-days and
# the whole dataset.
