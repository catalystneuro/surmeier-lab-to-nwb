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
from datetime import datetime, timedelta
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
from neuroconv.utils import load_dict_from_file
from pynwb import NWBFile
from pynwb.base import TimeSeries, TimeSeriesReference, TimeSeriesReferenceVectorData
from pynwb.device import Device
from pynwb.file import Subject
from pynwb.icephys import (
    CurrentClampSeries,
    CurrentClampStimulusSeries,
    IntracellularElectrode,
    VoltageClampSeries,
    VoltageClampStimulusSeries,
)
from pynwb.image import Image

from surmeier_lab_to_nwb.zhai2025.behavior_sessions import (
    _CRISPR_INJECTION,
    _CRISPR_VIRUSES,
)
from surmeier_lab_to_nwb.zhai2025.devices import get_or_create_bruker_ultima_model
from surmeier_lab_to_nwb.zhai2025.interfaces.bruker_line_scan_interface import (
    BrukerLineScanInterface,
)
from surmeier_lab_to_nwb.zhai2025.mouse_day_merger import (
    MouseDayBundle,
    SweepRecord,
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


# --------- line scans ----------------------------------------------------------------------


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


_INDICATOR = {"Ch1": "Alexa Fluor 568", "Ch2": "Fluo-4"}
_DETECTOR = {"Ch1": "Red detector, 580-620 nm.", "Ch2": "Green detector, 490-560 nm."}


def _add_line_scan(
    nwbfile: NWBFile,
    *,
    sweeps: list[SweepRecord],
    sweep_start_times: list[float],
    channel: str,
    name: str,
    about: str,
    kymograph_note: str,
    roi_response_note: str,
    images: list[tuple[str, str]],
) -> tuple[TimeSeries, TimeSeries, dict[str, tuple[int, int]]]:
    """Write one channel of the line scans of sweeps that share a scan line, each sweep at its start time.

    The objects are named with `name` after their type (`LineScan{name}`, `RoiResponseSeries{name}`,
    `PlaneSegmentation{name}`), `about` says what was scanned and the notes end the descriptions of the two series.
    `images` gives the name and description of each sweep's source image, in the order of `sweeps`. Returns the
    line scan, the fluorescence and the first sample and number of samples of each sweep in them, by folder name.
    """
    interface = BrukerLineScanInterface(
        folder_paths=[sweep.recording_path for sweep in sweeps], channel_name=channel, metadata_key=name
    )
    for sweep, sweep_start_time in zip(sweeps, sweep_start_times):
        interface.temporal_alignment[sweep.recording_path.name].move_start_to(sweep_start_time)
    trial = interface.segmentation_extractor.get_trials()[0]
    start_x, start_y, stop_x, stop_y = trial.line
    x0, x1 = trial.profile_pixels
    dye = _DYE_NAME[channel]

    metadata = interface.get_metadata()
    metadata["Devices"]["bruker_device"]["name"] = "BrukerUltima"  # the microscope `_ensure_devices` added
    ophys = metadata["Ophys"]
    ophys["ImagingPlanes"][name].update(
        name=f"ImagingPlane{dye}",
        description=f"Two-photon line scans of dendrites in the dorsolateral striatum, {_DYE_DESCRIPTION[channel]}.",
        excitation_lambda=810.0,
        indicator=_INDICATOR[channel],
        location="Caudoputamen",
        optical_channel=[dict(name=f"OpticalChannel{dye}", description=_DETECTOR[channel], emission_lambda=np.nan)],
    )
    ophys["PlaneSegmentations"][name].update(
        name=f"PlaneSegmentation{name}",
        description=(
            f"The ROI of RoiResponseSeries{name}: pixels {x0} to {x1} along the scan line, counted from its start "
            "at 0, which Prairie View averaged into its line profile. Its coordinates are pixels of the source images."
        ),
    )
    ophys["RoiResponses"][name]["raw"].update(
        name=f"RoiResponseSeries{name}",
        description=(
            f"Fluorescence of {about}: the mean over pixels {x0} to {x1} along the scan line (the ROI) at each scan, "
            f"from Prairie View's LineProfileData file.{roi_response_note}"
        ),
    )
    ophys["LineScanKymographs"][name].update(
        name=f"LineScan{name}",
        description=(
            f"Bruker line scan of {about}: one row per scan of the line, one column per pixel along it, from pixel "
            f"({start_x}, {start_y}) to ({stop_x}, {stop_y}) of the source images in ImageLineScanSource."
            f"{kymograph_note}"
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


def _synthesize_stimulus_waveform(
    protocol: StimulusProtocol,
    n_samples: int,
    rate_hz: float,
) -> np.ndarray:
    """Build the commanded current waveform (in amperes) from a parsed protocol: its pulses on a zero baseline.

    The waveform covers the recorded sweep, or the whole command when the recording stops before it ends.
    """
    n_samples = max(n_samples, int(round(protocol.duration_ms * 1e-3 * rate_hz)))
    waveform = np.zeros(n_samples, dtype=np.float32)
    for start_ms, width_ms in protocol.pulses_ms:
        start_sample = int(round(start_ms * 1e-3 * rate_hz))
        end_sample = start_sample + int(round(width_ms * 1e-3 * rate_hz))
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


def _series_name(cell_index: int, modality: str, neurodata_type: str = "CurrentClampSeries") -> str:
    suffix = "Somatic" if modality == "soma_excitability" else "Dendritic"
    return f"{neurodata_type}Cell{cell_index}{suffix}"


# The lab's figure folders follow the preprint's numbering. The published paper (Science Advances) gave the behaviour of
# preprint Figures 7 and 8 figures of their own, so the M1R deletion recordings of the lab's Figure 8 are Figure 9.
_PUBLISHED_FIGURE = {"F1": "Figure 1", "F3": "Figure 3", "F6": "Figure 6", "F7": "Figure 7", "F8": "Figure 9"}
_CELLS_OF_TYPE = {
    "dSPN": ("direct-pathway spiny projection neuron (dSPN)", "direct-pathway spiny projection neurons (dSPNs)"),
    "iSPN": ("indirect-pathway spiny projection neuron (iSPN)", "indirect-pathway spiny projection neurons (iSPNs)"),
}


def _patch_clamp_session_description(bundle: MouseDayBundle, cell_type: str) -> str:
    """What was recorded from the mouse that day: figure, cells, conditions in recording order, and the two assays."""
    singular, plural = _CELLS_OF_TYPE[cell_type]
    n_cells = len(bundle.cell_indices())
    conditions: list[str] = []
    for sweep in sorted((sw for sweeps in bundle.cells.values() for sw in sweeps), key=lambda sw: sw.pv_scan_date):
        condition = _CONDITION_LABEL[sweep.condition_subfolder]
        if sweep.suffix == "oxoM":
            condition += " with oxo-M"
        if condition not in conditions:
            conditions.append(condition)
    somatic = [cell for cell in bundle.cell_indices() if bundle.sweeps_for(cell, "soma_excitability")]
    dendritic = [cell for cell in bundle.cell_indices() if bundle.sweeps_for(cell, "dend_excitability")]

    description = (
        f"Patch-clamp recordings ({_PUBLISHED_FIGURE[bundle.figure_code]}): whole-cell current clamp of {n_cells} "
        f"{singular if n_cells == 1 else plural} from one mouse. "
        f"{'Condition' if len(conditions) == 1 else 'Conditions'}: {'; '.join(conditions)}."
    )
    if somatic:
        description += (
            f" Somatic excitability of {len(somatic)} {'cell' if len(somatic) == 1 else 'cells'}: 500 ms current "
            "steps injected at the cell body."
        )
    if dendritic:
        description += (
            f" Dendritic excitability of {len(dendritic)} {'cell' if len(dendritic) == 1 else 'cells'}: spikes evoked "
            "by brief current injections at the cell body, imaged with two-photon line scans across proximal and "
            "distal dendrites filled with Alexa Fluor 568 and Fluo-4."
        )
    return (
        description + " The cell, condition and site of each sweep are in the intracellular recordings table and the "
        "RecordingsIndexTable."
    )


def _protocol_description(protocol_label: str) -> str:
    """Readable protocol for the repetitions table from the internal label ('soma_fi_curve', 'dend_prox1')."""
    if protocol_label == "soma_fi_curve":
        return "F-I curve"
    site = protocol_label.removeprefix("dend_")
    dendrite_type = "Proximal" if site.startswith("prox") else "Distal"
    return f"back-propagating spikes at {dendrite_type}{site[4:]}"


# The lab's condition folders, named as in the paper's figure legends. The lab named the sulpiride
# folders differently for somatic and dendritic recordings; both are one condition.
_CONDITION_LABEL = {
    "LID off-state": "LID off-state",
    "LID on-state": "LID on-state",
    "LID on-state with SCH": "LID on-state with SCH 23390",
    "LID on-state with sul": "LID on-state with sulpiride",
    "LID on-state with sul (iSPN)": "LID on-state with sulpiride",
    "KO off-state": "CDGI knockout, LID off-state",
    "KO on-state": "CDGI knockout, LID on-state",
    "control": "LID off-state, vehicle control",
    "M1R antagonist": "LID off-state with M1R antagonists",
    "interleaved control": "LID off-state, control (gRNA only)",
    "M1R CRISPR": "LID off-state, M1R deletion",
    "UL control": "unlesioned control",
    "PD": "6-OHDA-lesioned",
    "LID off": "LID off-state",
    "WT": "naive",  # Figure 7E; a trial after oxo-M adds " with oxo-M"
    "KO": "CDGI knockout, naive",
}

# Figure 7E: the paper calls the mice of the oxo-M experiment naive; the lab's Data Connections sheet labels their days
# "LID ON". We follow the paper and say so in the subject description.
_NAIVE_SURGERY = "None: naive mouse, with no 6-OHDA lesion, as the paper describes the mice of Figure 7D-E."
_OXO_M_PHARMACOLOGY = (
    "No levodopa. During recording, the muscarinic receptor agonist oxotremorine M (oxo-M, 10 µM) was applied in the "
    "bath: each cell was recorded at all its dendritic sites, then again at the same sites after oxo-M was added (the "
    "lab's ReadMe: trials 1 to 3 at a site are before oxo-M, trials 4 to 6 or 7 to 9 after). Each trial's condition "
    "says whether oxo-M was present."
)
# Figure 6: the M1R antagonist or saline, injected in the off-state; the antagonist mice's slices were kept in a second
# antagonist.
_M1R_ANTAGONIST_PHARMACOLOGY = {
    "M1R ant": " About 2 hours after the last levodopa dose, at the start of the off-state, the M1R antagonist "
    "trihexyphenidyl hydrochloride (3 mg/kg) was injected i.p.; brain slices were made 16 hours later and kept in the "
    "selective M1R antagonist VU 0255035 (5 µM) throughout the recordings.",
    "CONT": " About 2 hours after the last levodopa dose, at the start of the off-state, saline was injected i.p. (the "
    "control for the M1R antagonist trihexyphenidyl); brain slices were made 16 hours later.",
}
# Figure 8: the CRISPR or control virus, as in the behaviour files, plus the reporter virus of the ex vivo experiments.
_CRISPR_REPORTER = (
    " As in all ex vivo experiments that needed the CRISPR-expressing iSPNs identified, a third virus, "
    "AAV9-EF1a-DIO-eYFP (Addgene #27056), was mixed in (1:1:1)."
)

# The subject description of the patch-clamp and Sr2+-oEPSC files: the mouse line, then what was done to the mouse, in
# the paper's words. Keyed by the line and group as the lab's Data Connections sheet writes them.
_LINE_DESCRIPTION = {
    "tD-Tomato": (
        "Drd1-tdTomato BAC (bacterial artificial chromosome) transgenic mouse, whose direct-pathway spiny projection "
        "neurons (dSPNs) express tdTomato."
    ),
    "eGFP": (
        "Drd2-eGFP BAC (bacterial artificial chromosome) transgenic mouse, whose indirect-pathway spiny projection "
        "neurons (iSPNs) express eGFP."
    ),
    "CDGIko X eGFP": (
        "CalDAG-GEFI (CDGI) knockout mouse crossed with the Drd2-eGFP BAC (bacterial artificial chromosome) "
        "transgenic line: it lacks CDGI, and its indirect-pathway spiny projection neurons (iSPNs) express eGFP."
    ),
    "Adora2-Cre": (
        "Adora2a-Cre transgenic mouse, whose indirect-pathway spiny projection neurons (iSPNs) express Cre "
        "recombinase."
    ),
}
_LESION = "It received a unilateral 6-hydroxydopamine (6-OHDA) lesion of the medial forebrain bundle"
_LEVODOPA = "levodopa (with benserazide) every other day for at least five sessions to induce dyskinesia"
_LESION_AND_LEVODOPA = f"{_LESION}, then {_LEVODOPA}"
_M1R_INJECTION = (
    ". About 2 hours after the last dose, at the start of the off-state, it was injected intraperitoneally with "
    "{drug}; brain slices were made 16 hours later."
)
_CRISPR = (
    f"{_LESION} and, after the cylinder test, an injection of {{viruses}} into the dorsolateral striatum on the side "
    f"of the lesion. About a month later it was given {_LEVODOPA}, and it was studied in the off-state of "
    "levodopa-induced dyskinesia (LID)."
)
_GROUP_DESCRIPTION = {
    "LID OFF": _LESION_AND_LEVODOPA + "; brain slices were made 24 to 48 hours after the last dose, in the off-state "
    "of levodopa-induced dyskinesia (LID).",
    "LID ON": _LESION_AND_LEVODOPA + "; brain slices were made 30 min after the last dose, in the on-state of "
    "levodopa-induced dyskinesia (LID).",
    "M1R ant": _LESION_AND_LEVODOPA
    + _M1R_INJECTION.format(drug="the M1 muscarinic receptor (M1R) antagonist trihexyphenidyl (3 mg/kg)"),
    "M1R CRSPR": _CRISPR.format(
        viruses="adeno-associated viruses (AAVs) carrying Cre-dependent Cas9 and guide RNAs against Chrm1, the gene "
        "of the M1 muscarinic receptor (M1R), which delete M1R from iSPNs,"
    ),
    "NAIVE": (
        "It is naive: no 6-hydroxydopamine (6-OHDA) lesion and no levodopa, as the paper describes the mice of Figure "
        "7D-E. The lab's Data Connections sheet labels this day 'LID ON' (lesioned and levodopa-treated) instead; the "
        "recordings themselves do not say."
    ),
}
# The lab's "CONT" is a different control in each figure.
_CONTROL_DESCRIPTION = {
    "F6": _LESION_AND_LEVODOPA
    + _M1R_INJECTION.format(drug="saline (the control for the M1 muscarinic receptor antagonist)"),
    "F8": _CRISPR.format(viruses="only the guide RNA virus, without Cas9 (control),"),
}


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
            "The lab gave this mouse no number; its Data Connections sheet identifies it by date, line and condition.",
        ]
        if sentence
    )


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


def _build_electrode(nwbfile: NWBFile, cell_index: int, device: Device) -> IntracellularElectrode:
    """The cell's whole-cell patch pipette, at the cell body; the dendritic sites were only imaged."""
    name = f"ElectrodeCell{cell_index}"
    if name in nwbfile.icephys_electrodes:
        return nwbfile.icephys_electrodes[name]
    return nwbfile.create_icephys_electrode(
        name=name,
        description=f"Whole-cell patch pipette at the cell body of cell {cell_index}.",
        device=device,
        location="Caudoputamen",
        cell_id=f"Cell{cell_index}",
    )


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
        list[tuple[int, int, float, StimulusProtocol, int, int]],
    ]
    | None
):
    """Concatenate all sweeps of (cell, modality) into one PatchClampSeries.

    Returns (response_series, stimulus_series, sweep_slices) where each slice carries
    (idx_start, count, sweep_start_time_seconds, stimulus_protocol, stimulus_start, stimulus_count); the
    stimulus slice differs from the response slice only for a sweep recorded shorter than its command. The stimulus
    protocol is parsed from the VoltageOutput XML alongside each sweep and used to
    populate the IRT's `stimulus_current_pA` column plus the synthesized stimulus waveform.
    """
    sweeps = bundle.sweeps_for(cell_index, modality)
    if not sweeps:
        return None

    voltage_chunks: list[np.ndarray] = []
    stim_chunks: list[np.ndarray] = []
    sweep_slices: list[tuple[int, int, float, StimulusProtocol, int, int]] = []
    sample_rate: float | None = None
    cumulative = 0
    stimulus_cumulative = 0
    # Durations of the sweeps whose recording stops before their current command ends.
    short_sweep_seconds: list[float] = []
    for sw in sweeps:
        voltage, rate = _read_sweep_voltage(sw.recording_path)
        if sample_rate is None:
            sample_rate = rate
        elif abs(rate - sample_rate) > 1.0:
            raise ValueError(f"{sw.recording_path} is sampled at {rate} Hz, the other sweeps at {sample_rate} Hz")
        # Parse the stimulus protocol from the VoltageOutput XML for this sweep.
        vo_xml = _find_voltage_output_xml(sw.recording_path)
        if vo_xml is None:
            raise FileNotFoundError(f"No VoltageOutput file with the current command in {sw.recording_path}")
        protocol = _parse_voltage_output(vo_xml)
        if sum(protocol.pulses_ms[-1]) > len(voltage) / sample_rate * 1e3:
            if modality != "dend_excitability":
                raise ValueError(f"{sw.recording_path} is recorded for less time than its current command lasts")
            short_sweep_seconds.append(len(voltage) / sample_rate)
        stim_waveform = _synthesize_stimulus_waveform(
            protocol,
            n_samples=len(voltage),
            rate_hz=float(sample_rate),
        )
        voltage_chunks.append(voltage)
        stim_chunks.append(stim_waveform)
        offset_s = (sw.pv_scan_date - session_start_time).total_seconds()
        sweep_slices.append((cumulative, len(voltage), offset_s, protocol, stimulus_cumulative, len(stim_waveform)))
        cumulative += len(voltage)
        stimulus_cumulative += len(stim_waveform)

    if not voltage_chunks:
        return None

    data_concat = np.concatenate(voltage_chunks)
    stim_concat = np.concatenate(stim_chunks)
    series_name = _series_name(cell_index, modality)
    if modality == "soma_excitability":
        recording = "while current steps were injected at the cell body to measure its firing (F-I curve)."
    else:
        recording = "during dendritic imaging." + _DENDRITIC_TRIAL
    short_sweeps = ""
    if short_sweep_seconds:
        seconds = " or ".join(f"{duration:g}" for duration in sorted(set(short_sweep_seconds)))
        short_sweeps = (
            f" {len(short_sweep_seconds)} of the sweeps were recorded for only {seconds} s, which ends before the "
            "current injections, so their voltage holds no spikes; the injections were given, as the calcium signal "
            "in each of those trials' line scans shows. Their stimulus is the whole command, so it is longer than "
            "their recording."
        )
    numbering = "".join(f" {note}" for note in sorted({sw.numbering_note for sw in sweeps if sw.numbering_note}))
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
            "recordings table." + short_sweeps + numbering
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
            "sweep's Bruker VoltageOutput XML (every current pulse, with its amplitude, width and timing); the "
            "per-sweep amplitude is also exposed via the IntracellularRecordingsTable.stimulus_current_pA column."
            + short_sweeps
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


def merge_one_mouse_day(bundle: MouseDayBundle, output_path: Path) -> Path:
    """Build one per-mouse-day merged NWB file and write it to disk."""
    session_start_time_naive = bundle.session_start_time
    session_start_time = session_start_time_naive.replace(tzinfo=TIMEZONE)

    general = _general_metadata()
    nwb_meta = general.get("NWBFile", {})
    subject_meta = general.get("Subject", {})
    subject_entry = bundle.subject_entry

    session_id = _build_session_id(bundle)
    session_description = _patch_clamp_session_description(bundle, cell_type=session_id.split("++")[0])
    surgery = nwb_meta.get("surgery")
    pharmacology = nwb_meta.get("pharmacology")
    if subject_entry.category == "NAIVE":
        surgery = _NAIVE_SURGERY
        pharmacology = _OXO_M_PHARMACOLOGY
    elif bundle.figure_code == "F6":
        pharmacology = pharmacology.rstrip() + _M1R_ANTAGONIST_PHARMACOLOGY[subject_entry.category]
    elif bundle.figure_code == "F8":
        viruses = _CRISPR_VIRUSES[subject_entry.category]
        surgery = surgery.rstrip() + _CRISPR_INJECTION.format(viruses=viruses) + _CRISPR_REPORTER
    nwbfile = NWBFile(
        session_description=session_description,
        experiment_description=nwb_meta.get(
            "experiment_description",
            "Excitability and synaptic plasticity of striatal projection neurons in 6-OHDA-lesioned mice "
            "undergoing levodopa-induced dyskinesia.",
        ),
        identifier=str(uuid.uuid4()),
        session_start_time=session_start_time,
        session_id=session_id,
        experimenter=nwb_meta.get("experimenter", ["Zhai, Shenyu"]),
        institution=nwb_meta.get("institution", "Northwestern University"),
        lab=nwb_meta.get("lab", "Surmeier Lab"),
        keywords=nwb_meta.get("keywords", []),
        related_publications=nwb_meta.get("related_publications", []),
        surgery=surgery,
        pharmacology=pharmacology,
        slices=nwb_meta.get("slices"),
    )
    nwbfile.subject = Subject(
        subject_id=subject_entry.subject_id,
        species=subject_meta.get("species", "Mus musculus"),
        strain=subject_entry.strain,
        genotype=subject_entry.genotype,
        sex=subject_meta.get("sex", "M"),
        age=subject_meta.get("age", "P7W/P12W"),
        description=_patch_clamp_subject_description(bundle),
    )

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

    # Add custom columns to IntracellularRecordingsTable BEFORE adding rows.
    # The `electrodes` aligned sub-table doesn't exist until the first add; using
    # the convenience add_intracellular_recording call below pre-creates it.
    irt = nwbfile.get_intracellular_recordings()
    irt.add_column(
        name="cell_id", description="Cell number within this file (the cells recorded from this mouse on this day)."
    )
    irt.add_column(
        name="dendrite_type",
        description="Dendritic region imaged during the sweep (Proximal or Distal); Soma for somatic sweeps, "
        "which have no imaging. The recording electrode is always at the cell body.",
    )
    irt.add_column(
        name="dendrite_distance_um",
        description="Approximate distance of the imaged site from the cell body in micrometers: 40 for Proximal and "
        "90 for Distal, as set in the lab protocol; 0 for somatic sweeps.",
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
            "Commanded current amplitude in picoamperes for this sweep: the step for a somatic sweep, each of the "
            "three pulses for a dendritic one. Parsed from the Bruker VoltageOutput XML alongside each sweep "
            "(PulsePotentialStart × UnitScaleFactor)."
        ),
    )
    irt.add_column(
        name="stimulus_protocol_name",
        description="Stimulus protocol label from Bruker VoltageOutput XML (e.g., 'Step 5_-40 pA').",
    )
    irt.add_column(
        name="condition",
        description=(
            "Experimental condition, named as in the paper's figure legends (e.g., 'LID on-state', "
            "'LID on-state with sulpiride'). Where a cell was recorded before and after a drug, its sweeps "
            "differ here."
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
        # One pipette per cell, at the cell body, for both the somatic and the dendritic sweeps.
        cell_electrode = _build_electrode(nwbfile, cell_index, amp)
        # Somatic block (one PatchClampSeries for all somatic sweeps of this cell).
        soma_result = _build_concatenated_series(
            nwbfile,
            bundle,
            cell_index,
            "soma_excitability",
            cell_electrode,
            session_start_time_naive,
        )
        if soma_result is not None:
            response, stim, sweep_slices = soma_result
            sorted_soma = sorted(bundle.sweeps_for(cell_index, "soma_excitability"), key=lambda s: s.pv_scan_date)
            for sw, (idx_start, count, sweep_start_time, protocol, stimulus_start, stimulus_count) in zip(
                sorted_soma, sweep_slices
            ):
                condition = _CONDITION_LABEL[sw.condition_subfolder]
                row_index = nwbfile.add_intracellular_recording(
                    electrode=cell_electrode,
                    stimulus=stim,
                    stimulus_start_index=stimulus_start,
                    stimulus_index_count=stimulus_count,
                    response=response,
                    response_start_index=idx_start,
                    response_index_count=count,
                    cell_id=cell_index,
                    dendrite_type="Soma",
                    dendrite_distance_um=0.0,
                    sweep_start_time=sweep_start_time,
                    stimulus_current_pA=protocol.amplitude_pA,
                    stimulus_protocol_name=protocol.name,
                    condition=condition,
                    site="Soma",
                )
                irt_rows_per_sequential[(cell_index, "soma_fi_curve", condition)].append(row_index)
                # Per-row metadata for the RecordingsIndexTable (built post-chain).
                index_row_metadata[row_index] = {
                    "cell_id": cell_index,
                    "condition": condition,
                    "dendrite_type": "Soma",
                    "dendrite_distance_um": 0.0,
                    "site": "Soma",
                    "stimulus_current_pA": protocol.amplitude_pA,
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
            figure_code=bundle.figure_code,
            cells={cell_index: dend_sweeps},
        )
        # Single concatenated dendritic series for the whole cell.
        dend_result = _build_concatenated_series(
            nwbfile,
            dend_subbundle,
            cell_index,
            "dend_excitability",
            cell_electrode,
            session_start_time_naive,
        )
        if dend_result is None:
            continue
        dend_response, dend_stim, dend_sweep_slices = dend_result
        # Add IRT rows.
        # The sweep order in dend_sweep_slices matches the chronological order used
        # in _build_concatenated_series (which sorts by pv_scan_date).
        sorted_sweeps = sorted(dend_sweeps, key=lambda s: s.pv_scan_date)
        # Sites the lab imaged at two spots on the dendrite; the second spot's folders say "real".
        sites_with_two_spots = {
            (sw.dendrite_location, sw.location_number) for sw in dend_sweeps if "real" in sw.site_label
        }
        for sw, (idx_start, count, sweep_start_time, protocol, stimulus_start, stimulus_count) in zip(
            sorted_sweeps, dend_sweep_slices
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

            row_index = nwbfile.add_intracellular_recording(
                electrode=cell_electrode,
                stimulus=dend_stim,
                stimulus_start_index=stimulus_start,
                stimulus_index_count=stimulus_count,
                response=dend_response,
                response_start_index=idx_start,
                response_index_count=count,
                cell_id=cell_index,
                dendrite_type=dendrite_type,
                dendrite_distance_um=dendrite_distance,
                sweep_start_time=sweep_start_time,
                stimulus_current_pA=protocol.amplitude_pA,
                stimulus_protocol_name=protocol.name,
                condition=condition,
                site=site,
            )
            protocol_label = f"dend_{sw.dendrite_location}{sw.location_number}"
            irt_rows_per_sequential[(cell_index, protocol_label, condition)].append(row_index)
            index_row_metadata[row_index] = {
                "cell_id": cell_index,
                "condition": condition,
                "dendrite_type": dendrite_type,
                "dendrite_distance_um": dendrite_distance,
                "site": site,
                "stimulus_current_pA": protocol.amplitude_pA,
                "sweep_start_time_s": float(sweep_start_time),
            }

            # Both channels share the scan line. Each spot of a site imaged at two spots is a location of its own.
            location = (cell_index, sw.dendrite_location, sw.location_number, spot)
            scan_lines = scan_lines_per_location[location]
            scan_line_number = scan_lines.setdefault(
                BrukerLineScanInterface.get_scan_line(sw.recording_path), len(scan_lines) + 1
            )
            line_scan_trials[(*location, scan_line_number)].append((row_index, sw, float(sweep_start_time)))
            trial_number_at_location[row_index] = sum(
                len(trials) for key, trials in line_scan_trials.items() if key[:4] == location
            )

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
        spot_note = f", spot {spot[-1]} of 2" if spot else ""
        site_note = _DENDRITIC_TRIAL + _spot_sentence(spot)
        line_suffix = ""
        n_scan_lines = len(scan_lines_per_location[location])
        if n_scan_lines > 1:
            line_suffix = f"Line{scan_line_number}"
            site_note += (
                f" The scan line at this location was moved or changed length between trials, so its trials "
                f"are split into {n_scan_lines} series; this one holds the trials recorded with line "
                f"{scan_line_number}."
            )
        for chan_name in ("Ch1", "Ch2"):
            dye = _DYE_NAME[chan_name]
            images = [
                (
                    f"Image{dye}Cell{cell_index}{dendrite_type_pretty}{loc_num}{spot}"
                    f"Trial{trial_number_at_location[row_index]}",
                    f"Source image (snapshot of the field of view) for Cell {cell_index} {dendrite_type_pretty} "
                    f"location {loc_num}{spot_note}, {_DYE_DESCRIPTION[chan_name]}, trial "
                    f"{trial_number_at_location[row_index]} of {trials_at_location} in recording order (row "
                    f"{row_index} of the intracellular recordings table, start time {sweep_start:.1f} s).",
                )
                for row_index, _, sweep_start in trials
            ]
            series_note = (
                " All trials at this site, in recording order and each at its own start time; each trial's samples "
                "are in the {column} column of processing/recordings_index/RecordingsIndexTable." + site_note
            )
            kymograph, roi_response, sample_ranges = _add_line_scan(
                nwbfile,
                sweeps=[sw for _, sw, _ in trials],
                sweep_start_times=[sweep_start for _, _, sweep_start in trials],
                channel=chan_name,
                name=f"{dye}Cell{cell_index}{dendrite_type_pretty}{loc_num}{spot}{line_suffix}",
                about=f"Cell {cell_index} {dendrite_type_pretty} location {loc_num}, {_DYE_DESCRIPTION[chan_name]}",
                kymograph_note=series_note.format(column=f"line_scan_{dye.lower()}"),
                roi_response_note=series_note.format(column=f"roi_response_{dye.lower()}"),
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
            background_note = (
                " One line scan with the scan line moved off the dendrite and no stimulus, recorded right after the "
                "trials at this site. The paper subtracts background from both channels before analysis; the lab "
                "took such scans only for this cell, on 2017-05-25 and 2017-05-26."
            )
            for chan_name in ("Ch1", "Ch2"):
                dye = _DYE_NAME[chan_name]
                _add_line_scan(
                    nwbfile,
                    sweeps=[sw],
                    sweep_start_times=[(sw.pv_scan_date - session_start_time_naive).total_seconds()],
                    channel=chan_name,
                    name=f"Background{dye}Cell{cell_index}{site}",
                    about=f"the background scan for Cell {cell_index} {site}, {_DYE_DESCRIPTION[chan_name]}",
                    kymograph_note=background_note,
                    roi_response_note=background_note,
                    images=[
                        (
                            f"Image{dye}Cell{cell_index}{site}Background",
                            f"Source image (snapshot of the field of view) for the background scan of Cell "
                            f"{cell_index} {site}, {_DYE_DESCRIPTION[chan_name]}.",
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
        name="condition",
        description=(
            "Experimental condition of this repetition, named as in the paper's figure legends. A cell recorded "
            "before and after a drug has one repetition per condition for the same protocol."
        ),
    )
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
            name="condition",
            description="Experimental condition, named as in the paper's figure legends (e.g., 'LID off-state').",
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
    # Use neuroconv's helper to auto-apply default HDF5 chunking + gzip compression
    # via get_default_backend_configuration (saves ~50% disk vs plain contiguous writes).
    from neuroconv.tools.nwb_helpers import configure_and_write_nwbfile

    configure_and_write_nwbfile(nwbfile=nwbfile, nwbfile_path=output_path, backend="hdf5")
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
    columns: list = [
        VectorData(name="cell_id", description="Cell number within this file.", data=cell_id_col),
        VectorData(
            name="condition",
            description="Experimental condition of this sweep; see the same column in the intracellular recordings table.",
            data=condition_col,
        ),
        VectorData(
            name="dendrite_type",
            description="Dendritic region imaged during the sweep (Proximal or Distal); Soma for somatic sweeps.",
            data=dendrite_type_col,
        ),
        VectorData(
            name="dendrite_distance_um",
            description="Approximate distance of the imaged site from the cell body in micrometers (0 for Soma, 40 for "
            "Proximal, 90 for Distal).",
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
            description=(
                "Commanded current amplitude in picoamperes for this sweep: the step for a somatic sweep, each of "
                "the three pulses for a dendritic one."
            ),
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
            "in one place. Each row holds the sweep's details (cell_id, condition, dendrite_type, "
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


# --------- Sr2+-oEPSC (Figures 2 and 4) ---------------------------------------------


# Reporter line -> (figure, recorded cells), for the session description.
_OEPSC_CELLS = {
    "tD-Tomato": ("Figure 2", "direct-pathway spiny projection neurons (dSPNs)"),
    "eGFP": ("Figure 4", "indirect-pathway spiny projection neurons (iSPNs)"),
}

_OEPSC_ANALYSIS_WINDOW = (
    "The paper counted asynchronous events from 40 ms to 400 ms after each flash onset (Methods, Data acquisition "
    "and statistical analysis)."
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


def _read_sweep_current(sweep_path: Path) -> tuple[np.ndarray, float]:
    """Read the current of one voltage-clamp sweep in amperes, and its sample rate.

    The VoltageRecording XML gives the recorded signal's unit and gain: pA = value * multiplier / divisor.
    """
    xml_candidates = list(sweep_path.glob("*_VoltageRecording_001.xml"))
    csv_candidates = list(sweep_path.glob("*_VoltageRecording_001.csv"))
    if len(xml_candidates) != 1 or len(csv_candidates) != 1:
        raise FileNotFoundError(f"Expected one VoltageRecording XML and CSV in {sweep_path}")
    signals = [s for s in ET.parse(xml_candidates[0]).getroot().iter("VRecSignal") if s.findtext("Enabled") == "true"]
    if len(signals) != 1 or signals[0].findtext("Unit/UnitName") != "pA":
        raise ValueError(f"Expected one recorded signal, in pA, in {xml_candidates[0]}")
    multiplier = float(signals[0].findtext("Unit/Multiplier"))
    divisor = float(signals[0].findtext("Unit/Divisor"))
    frame = pd.read_csv(csv_candidates[0], header=0)
    time_ms = frame.iloc[:, 0].to_numpy()
    current_amperes = frame.iloc[:, 1].to_numpy() * multiplier / divisor * 1e-12
    sample_rate = round(1000.0 / (time_ms[1] - time_ms[0]))
    return current_amperes.astype(np.float32), float(sample_rate)


def _add_optogenetics_metadata(nwbfile: NWBFile) -> OptogeneticSitesTable:
    """The LED, the ChR2 virus and its injection, from the paper's Methods; returns the sites table."""
    led_model = ExcitationSourceModel(
        name="CoolLEDModel",
        manufacturer="CoolLED",
        model_number="p100-470",
        source_type="LED",
        excitation_mode="one-photon",
        description="CoolLED p100-470 blue LED (CoolLED, Andover, UK), as named in the paper's Methods.",
    )
    nwbfile.add_device_model(led_model)
    led = ExcitationSource(
        name="BlueLED",
        description=(
            "Blue LED (470 nm) giving whole-field illumination of the slice. Prairie View sent it a 5 V pulse for "
            "each flash; the paper does not report the light power."
        ),
        model=led_model,
    )
    nwbfile.add_device(led)
    virus = ViralVector(
        name="ChR2Virus",
        construct_name="AAV5-hSyn-hChR2(H134R)-eYFP",
        manufacturer="Addgene (#26973)",
        titer_in_vg_per_ml=float("nan"),
        description="Expresses channelrhodopsin-2 (H134R) fused to eYFP under the human synapsin promoter. The paper "
        "does not report the titer.",
    )
    injection = ViralVectorInjection(
        name="ChR2VirusInjection",
        description=(
            "0.15 µl injected into the M1 motor cortex on the side of the 6-OHDA lesion, at ML -1.60 mm, AP 1.15 mm "
            "and DV 1.55 mm relative to Bregma as given in the paper (DV stored as -1.55, ventral). Recordings were "
            "made about 3 to 4 weeks after the injection."
        ),
        location="Primary motor area",
        # In the paper negative ML is left (it injects "the left DLS" at ML -2.35), and the lesion protocol scores
        # right-paw use, so the lesion and this injection are on the left.
        hemisphere="left",
        reference="Bregma",
        ap_in_mm=1.15,
        ml_in_mm=-1.60,
        dv_in_mm=-1.55,
        volume_in_uL=0.15,
        viral_vector=virus,
    )
    effector = Effector(
        name="ChR2",
        label="hChR2(H134R)-eYFP",
        description="Channelrhodopsin-2 in motor cortex neurons; the LED activates their axons in the striatum.",
        viral_vector_injection=injection,
    )
    sites_table = OptogeneticSitesTable(
        name="optogenetic_sites_table",
        description="The one stimulation site: whole-field LED illumination of ChR2-expressing cortical axons in the "
        "slice.",
    )
    sites_table.add_row(excitation_source=led, effector=effector)
    nwbfile.add_lab_meta_data(
        OptogeneticExperimentMetadata(
            stimulation_software="Prairie View 5.3",
            optogenetic_sites_table=sites_table,
            optogenetic_viruses=OptogeneticViruses(viral_vectors=[virus]),
            optogenetic_virus_injections=OptogeneticVirusInjections(viral_vector_injections=[injection]),
            optogenetic_effectors=OptogeneticEffectors(effectors=[effector]),
        )
    )
    return sites_table


def merge_oepsc_mouse_day(bundle: MouseDayBundle, output_path: Path) -> Path:
    """Build one Sr2+-oEPSC NWB file for a mouse-day (Figures 2 and 4) and write it to disk."""
    session_start_time_naive = bundle.session_start_time
    general = _general_metadata()
    nwb_meta = general["NWBFile"]
    subject_meta = general["Subject"]
    subject_entry = bundle.subject_entry
    figure, cells = _OEPSC_CELLS[subject_entry.animal]
    conditions = {_CONDITION_LABEL[sw.condition_subfolder] for sweeps in bundle.cells.values() for sw in sweeps}
    if len(conditions) != 1:
        raise ValueError(f"{subject_entry.subject_id} has Sr2+-oEPSC recordings under {sorted(conditions)}")
    (condition,) = conditions

    nwbfile = NWBFile(
        session_description=(
            f"Sr2+-oEPSC recordings ({figure}): whole-cell voltage clamp of {len(bundle.cell_indices())} {cells} "
            f"from one mouse in the {condition}. In Ca2+-free solution with 3 mM SrCl2, a blue LED flash on motor "
            "cortex axons expressing ChR2 evokes asynchronous excitatory synaptic currents. One sweep per flash, "
            "most of them 30 s apart."
        ),
        experiment_description=(
            "Strength of cortical synapses onto striatal spiny projection neurons in 6-OHDA-lesioned mice in the "
            "off- and on-states of levodopa-induced dyskinesia (LID)."
        ),
        identifier=str(uuid.uuid4()),
        session_start_time=session_start_time_naive.replace(tzinfo=TIMEZONE),
        session_id=_build_session_id(bundle),
        experimenter=nwb_meta["experimenter"],
        institution=nwb_meta["institution"],
        lab=nwb_meta["lab"],
        keywords=nwb_meta["keywords"],
        related_publications=nwb_meta["related_publications"],
        surgery=(
            nwb_meta["surgery"] + "ChR2 virus (AAV5-hSyn-hChR2(H134R)-eYFP, 0.15 µl) injected into the M1 motor "
            "cortex on the side of the lesion, at ML -1.60, AP 1.15 and DV 1.55 mm relative to Bregma; recordings "
            "about 3 to 4 weeks later."
        ),
        pharmacology=(
            nwb_meta["pharmacology"] + "During recording the slice was superfused with Ca2+-free ACSF containing "
            "3 mM SrCl2 and 10 µM gabazine (to block GABA-A currents), starting 25 min before recording."
        ),
        slices=nwb_meta["slices"],
    )
    nwbfile.subject = Subject(
        subject_id=subject_entry.subject_id,
        species=subject_meta["species"],
        strain=subject_entry.strain,
        genotype=subject_entry.genotype,
        sex=subject_meta["sex"],
        age=subject_meta["age"],
        description=_patch_clamp_subject_description(
            bundle,
            virus="A channelrhodopsin-2 (ChR2) virus was injected into the motor cortex on the side of the lesion; see "
            "surgery.",
        ),
    )
    amplifier = Device(name="MultiClamp700B", description="Molecular Devices MultiClamp 700B patch-clamp amplifier.")
    nwbfile.add_device(amplifier)
    sites_table = _add_optogenetics_metadata(nwbfile)

    pulses_table = OptogeneticPulsesTable(
        name="optogenetic_pulses_table",
        description=(
            "One row per LED flash: the start and end of the 5 V pulse Prairie View sent to the LED, 20 ms into each "
            "sweep, at 470 nm. power_in_mW is NaN because the paper does not report the light power. "
            + _OEPSC_ANALYSIS_WINDOW
        ),
        target_tables={"optogenetic_sites": sites_table},
    )
    pulses_table.add_column(name="cell_id", description="Cell recorded during this flash.")
    pulses_table.add_column(
        name="led_setting",
        description=(
            "LED setting the lab wrote in the recording folder names (20 for 'cell1_LED20-001'); the files do not "
            "say what it measures."
        ),
    )
    nwbfile.add_time_intervals(pulses_table)

    irt = nwbfile.get_intracellular_recordings()
    irt.add_column(
        name="cell_id", description="Cell number within this file (the cells recorded from this mouse on this day)."
    )
    irt.add_column(name="condition", description="Experimental condition, named as in the paper's figure legends.")
    irt.add_column(name="site", description="Recording site: 'Soma', the cell body, where the pipette was.")
    irt.add_column(
        name="sweep_start_time",
        description=(
            "Start of the sweep in seconds from session_start_time, from the Prairie View XML at 1-second precision. "
            "The voltage-clamp series store the sweeps back to back, so their own time axis does not show the gaps "
            "between sweeps."
        ),
    )

    irt_rows_per_cell: dict[int, list[int]] = defaultdict(list)
    for cell_index in bundle.cell_indices():
        electrode = _build_electrode(nwbfile, cell_index, amplifier)
        sweeps = bundle.sweeps_for(cell_index)
        currents, commands, sweep_slices = [], [], []
        sample_rate = None
        cumulative = 0
        for sw in sweeps:
            current, rate = _read_sweep_current(sw.recording_path)
            if sample_rate is None:
                sample_rate = rate
            elif rate != sample_rate:
                raise ValueError(f"{sw.recording_path} is sampled at {rate} Hz, the other sweeps at {sample_rate} Hz")
            protocol = _read_light_evoked_protocol(sw.recording_path)
            command = np.zeros(len(current), dtype=np.float32)
            step_start = round(protocol.step_delay_ms * 1e-3 * rate)
            command[step_start : step_start + round(protocol.step_width_ms * 1e-3 * rate)] = (
                protocol.step_amplitude_volts
            )
            sweep_start_time = (sw.pv_scan_date - session_start_time_naive).total_seconds()
            currents.append(current)
            commands.append(command)
            sweep_slices.append((cumulative, len(current), sweep_start_time, protocol, sw))
            cumulative += len(current)

        response = VoltageClampSeries(
            name=f"VoltageClampSeriesCell{cell_index}",
            description=(
                f"Whole-cell voltage-clamp current from cell {cell_index}, held at -70 mV (paper Methods). Each 1 s "
                "sweep has one blue LED flash at 20 ms (see the optogenetic pulses table) and a -5 mV, 50 ms command "
                f"step at 920 ms. {_OEPSC_ANALYSIS_WINDOW} The {len(sweeps)} sweeps are stored back to back; each "
                "sweep's start time is in the sweep_start_time column of the intracellular recordings table."
            ),
            data=np.concatenate(currents),
            electrode=electrode,
            rate=sample_rate,
            starting_time=sweep_slices[0][2],
        )
        nwbfile.add_acquisition(response)
        stimulus = VoltageClampStimulusSeries(
            name=f"VoltageClampStimulusSeriesCell{cell_index}",
            description=(
                f"Command potential Prairie View sent to the amplifier during the sweeps of cell {cell_index}, read "
                "from each sweep's VoltageOutput file: 0, then a -5 mV step of 50 ms at 920 ms. The amplifier added "
                "it to the -70 mV holding potential (paper Methods), which the files do not record."
            ),
            data=np.concatenate(commands),
            electrode=electrode,
            rate=sample_rate,
            starting_time=sweep_slices[0][2],
        )
        nwbfile.add_stimulus(stimulus)

        for start_index, count, sweep_start_time, protocol, sw in sweep_slices:
            row_index = nwbfile.add_intracellular_recording(
                electrode=electrode,
                stimulus=stimulus,
                stimulus_start_index=start_index,
                stimulus_index_count=count,
                response=response,
                response_start_index=start_index,
                response_index_count=count,
                cell_id=cell_index,
                condition=condition,
                site="Soma",
                sweep_start_time=sweep_start_time,
            )
            irt_rows_per_cell[cell_index].append(row_index)
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

    # One sequential recording and one repetition per cell, all under the mouse-day's one condition.
    repetitions_table = nwbfile.get_icephys_repetitions()
    repetitions_table.add_column(name="cell_id", description="Cell this repetition belongs to.")
    repetitions_table.add_column(name="protocol", description="Protocol of this repetition.")
    repetitions_table.add_column(name="condition", description="Experimental condition of this repetition.")
    repetition_indices = []
    for cell_index, irt_rows in irt_rows_per_cell.items():
        simultaneous = [nwbfile.add_icephys_simultaneous_recording(recordings=[row]) for row in irt_rows]
        sequential = nwbfile.add_icephys_sequential_recording(
            simultaneous_recordings=simultaneous, stimulus_type="Sr2+-oEPSC"
        )
        repetition_indices.append(
            nwbfile.add_icephys_repetition(
                sequential_recordings=[sequential], cell_id=cell_index, protocol="Sr2+-oEPSC", condition=condition
            )
        )
    conditions_table = nwbfile.get_icephys_experimental_conditions()
    conditions_table.add_column(
        name="condition", description="Experimental condition, named as in the paper's figure legends."
    )
    nwbfile.add_icephys_experimental_condition(repetitions=repetition_indices, condition=condition)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    from neuroconv.tools.nwb_helpers import configure_and_write_nwbfile

    configure_and_write_nwbfile(nwbfile=nwbfile, nwbfile_path=output_path, backend="hdf5")
    return output_path


# --------- GRAB-ACh (Figure 5) ---------------------------------------------------------


# The bath drug, from the words before the stimulus in each recording's name, with the paper's concentrations.
_GRAB_ACH_DRUG = {
    "ctr": "control",
    "50nMDA": "50 nM dopamine",
    "quin": "10 µM quinpirole",
    "sul": "10 µM sulpiride",
    "TTX": "10 µM tetrodotoxin",
    "sul_TTX": "10 µM tetrodotoxin",  # 2024-04-12, after sulpiride, like every tetrodotoxin trial
    "ACh": "100 µM acetylcholine",
    "ACh wash in": "100 µM acetylcholine",
}
# Stimulus words in the recording names, misspellings included; the stimulus is read from the VoltageOutput file.
_GRAB_ACH_STIMULUS_WORDS = {"single", "sinlge", "singl", "burst", "20Hz"}

# Prairie View channel number -> channel name in the file, and what the channel records.
_GRAB_ACH_CHANNELS = {"2": "GRABACh3", "3": "Dodt"}
_GRAB_ACH_CHANNEL_DESCRIPTION = {
    "GRABACh3": (
        "GRAB-ACh3.0 fluorescence (acetylcholine sensor), Hamamatsu H7422P-40 GaAsP detector with a 490 to 560 nm "
        "filter (Prairie View channel 2)"
    ),
    "Dodt": (
        "laser light transmitted through the slice, recorded by the Dodt contrast detector (Prairie View channel 3); "
        "it shows the tissue, not fluorescence"
    ),
}

_GRAB_ACH_INJECTION = (
    "GRAB-ACh3.0 virus (AAV5-hSyn-ACh3.0(ACh4.3), WZ Biosciences YL001003-AV5, 0.6 µl) injected into the "
    "dorsolateral striatum{side} at ML -2.2, AP 0.8 and DV 3.3 mm relative to Bregma; slices were made about 3 weeks "
    "later."
)

_GRAB_ACH_BATH = (
    "During imaging the slice was perfused with ACSF at 32 to 34 °C containing 50 µM sodium metabisulfite to slow "
    "dopamine oxidation. Drugs were applied in the bath in this order, at least 5 min each: control, 50 nM dopamine "
    "(lesioned mice only), 10 µM quinpirole (D2 receptor agonist) and 10 µM sulpiride (D2 receptor antagonist). "
    "10 µM tetrodotoxin and 100 µM acetylcholine chloride were applied last, to find the sensor's lowest and highest "
    "brightness. The paper does not say whether a drug was washed out before the next. Each trial's drug is in the "
    "trials table."
)

_GRAB_ACH_SUBJECT = {
    "CONT": "unlesioned control",
    "6-OHDA": "unilateral 6-OHDA lesion, no levodopa",
    "LID OFF": "unilateral 6-OHDA lesion and levodopa treatment, recorded in the off-state",
}


@dataclass
class SensorTrial:
    """One GRAB-ACh recording, read from its Prairie View files. Times are seconds from the session start."""

    sweep: SweepRecord
    start_time: float
    frame_times: np.ndarray
    frame_period: float
    frames: dict[str, np.ndarray]  # channel -> (frames, x, y)
    brightness: dict[str, np.ndarray]  # channel -> Prairie View's brightness of the whole frame, one value per frame
    references: dict[str, np.ndarray]  # channel -> Prairie View's 16-bit reference image, (x, y)
    zoom: float
    microns_per_pixel: float
    drug: str
    stimulation: str
    stimulus_start_time: float  # NaN without stimulation


def _read_grab_ach_stimulus(voltage_output_path: Path) -> tuple[int, float, float, float]:
    """Pulse count, width (ms), period (ms) and first-pulse delay (ms) of the stimulus in a VoltageOutput file.

    The stimulator is triggered from the output Prairie View calls LED, with 5 V pulses; burst recordings also carry a
    0 V pulse train there, which sends nothing.
    """
    root = ET.parse(voltage_output_path).getroot()
    waveforms = [w for w in root.findall("Waveform") if w.findtext("Enabled") == "true"]
    if [w.findtext("Name") for w in waveforms] != ["LED"]:
        raise ValueError(f"Expected one enabled output, LED, in {voltage_output_path}")
    trains = [
        c
        for c in waveforms[0]
        if c.tag.startswith("WaveformComponent") and float(c.findtext("PulsePotentialStart", "0")) != 0
    ]
    if len(trains) != 1 or trains[0].tag != "WaveformComponent_PulseTrain" or trains[0].findtext("Repetitions") != "1":
        raise ValueError(f"Expected one pulse train on the LED output of {voltage_output_path}")
    train = trains[0]
    width = float(train.findtext("PulseWidth"))
    return (
        int(train.findtext("PulseCount")),
        width,
        width + float(train.findtext("PulseSpacing")),
        float(train.findtext("FirstPulseDelay")),
    )


def _read_grab_ach_trial(sweep: SweepRecord, session_start_time: datetime) -> SensorTrial:
    """Read the frames, brightness, reference images and stimulus of one GRAB-ACh recording."""
    import tifffile

    leaf = sweep.recording_path
    xml_path = leaf / f"{leaf.name}.xml"
    root = ET.parse(xml_path).getroot()
    state = {e.get("key"): e for e in root.find("PVStateShard").findall("PVStateValue")}
    width = int(state["pixelsPerLine"].get("value"))
    height = int(state["linesPerFrame"].get("value"))
    microns = {v.get("index"): float(v.get("value")) for v in state["micronsPerPixel"].findall("IndexedValue")}
    if microns["XAxis"] != microns["YAxis"]:
        raise ValueError(f"Pixels are not square in {xml_path}")
    sequences = root.findall("Sequence")
    if len(sequences) != 1 or sequences[0].get("type") != "BrightnessOverTime":
        raise ValueError(f"Expected one brightness-over-time sequence in {xml_path}")
    sequence = sequences[0]

    # The sequence start, to the microsecond; the PVScan date holds the same moment truncated to the second.
    hours, minutes, seconds = sequence.get("time").split(":")
    started = datetime.combine(sweep.pv_scan_date.date(), datetime.min.time()) + timedelta(
        hours=int(hours), minutes=int(minutes), seconds=float(seconds)
    )
    if not timedelta(0) <= started - sweep.pv_scan_date < timedelta(seconds=1):
        raise ValueError(f"The sequence time of {xml_path} does not match its date")
    start_time = (started - session_start_time).total_seconds()

    frame_elements = sequence.findall("Frame")
    files = {"2": [], "3": []}
    for frame in frame_elements:
        names = {f.get("channel"): (f.get("channelName"), f.get("filename")) for f in frame.findall("File")}
        if names.keys() != files.keys() or names["2"][0] != "Ch2" or names["3"][0] != "Dodt":
            raise ValueError(f"Expected channels Ch2 and Dodt in every frame of {xml_path}")
        for channel, (_name, filename) in names.items():
            files[channel].append(filename)
    relative_times = np.array([float(frame.get("absoluteTime")) for frame in frame_elements])
    frames = {}
    for channel, filenames in files.items():
        # is_ome=False reads only this file's frame; the first file's OME header would pull in the whole series.
        stack = np.stack([tifffile.imread(leaf / filename, is_ome=False) for filename in filenames])
        if stack.shape[1:] != (height, width):
            raise ValueError(f"Frames of {leaf} are {stack.shape[1:]}, the XML says {(height, width)}")
        frames[_GRAB_ACH_CHANNELS[channel]] = stack.transpose(0, 2, 1)

    # Prairie View's brightness-over-time regions: each the whole frame, one per channel, in botData.csv.
    bots = sequence.find("PVBOTs")
    regions = bots.findall("Region")
    if sorted(r.get("channel") for r in regions) != ["2", "3"] or any(
        (r.get("x"), r.get("y"), r.get("width"), r.get("height")) != ("0", "0", str(width), str(height))
        for r in regions
    ):
        raise ValueError(f"Expected one whole-frame brightness region per channel in {xml_path}")
    table = pd.read_csv(leaf / bots.get("botData"))
    if len(table) != len(relative_times) or not np.allclose(table["Timestamp"], relative_times, atol=1e-6):
        raise ValueError(f"The brightness samples of {leaf} do not match its frames")
    brightness = {
        _GRAB_ACH_CHANNELS[r.get("channel")]: table[f"Region {number}"].to_numpy(dtype=np.float32)
        for number, r in enumerate(regions, start=1)
    }

    references = {
        name: tifffile.imread(leaf / "References" / f"{leaf.name}-Cycle00001-Ch{channel}-16bit-Reference.tif").T
        for channel, name in _GRAB_ACH_CHANNELS.items()
    }

    words = sweep.suffix.rsplit("_", 1)
    drug_words = words[0] if len(words) == 2 and words[1] in _GRAB_ACH_STIMULUS_WORDS else sweep.suffix
    if drug_words not in _GRAB_ACH_DRUG:
        raise ValueError(f"Unknown drug in the name of {leaf}")

    voltage_output = sequence.find("VoltageOutput")
    if voltage_output is None:
        stimulation, stimulus_start_time = "none", float("nan")
    else:
        count, pulse_width, period, delay = _read_grab_ach_stimulus(leaf / voltage_output.get("filename"))
        if count == 1:
            stimulation = f"single {pulse_width:g} ms pulse"
        else:
            stimulation = f"{count} pulses of {pulse_width:g} ms at {1000 / period:g} Hz"
        stimulus_start_time = start_time + float(voltage_output.get("absoluteTime")) + delay * 1e-3

    return SensorTrial(
        sweep=sweep,
        start_time=start_time,
        frame_times=start_time + relative_times,
        frame_period=float(state["framePeriod"].get("value")),
        frames=frames,
        brightness=brightness,
        references=references,
        zoom=float(state["opticalZoom"].get("value")),
        microns_per_pixel=microns["XAxis"],
        drug=_GRAB_ACH_DRUG[drug_words],
        stimulation=stimulation,
        stimulus_start_time=stimulus_start_time,
    )


def merge_grab_ach_mouse_day(bundle: MouseDayBundle, output_path: Path) -> Path:
    """Build one GRAB-ACh NWB file for a mouse-day (Figure 5) and write it to disk."""
    from pynwb.image import Images
    from pynwb.ophys import (
        Fluorescence,
        ImageSegmentation,
        OpticalChannel,
        RoiResponseSeries,
        TwoPhotonSeries,
    )

    session_start_time_naive = bundle.session_start_time
    general = _general_metadata()
    nwb_meta = general["NWBFile"]
    subject_meta = general["Subject"]
    subject_entry = bundle.subject_entry
    conditions = {_CONDITION_LABEL[sw.condition_subfolder] for sweeps in bundle.cells.values() for sw in sweeps}
    if len(conditions) != 1:
        raise ValueError(f"{subject_entry.subject_id} has GRAB-ACh recordings under {sorted(conditions)}")
    (condition,) = conditions
    regions = bundle.cell_indices()
    slices = sorted({slice_number for slice_number, _roi in regions})

    if subject_entry.category == "CONT":
        surgery = "No 6-OHDA lesion (unlesioned control). " + _GRAB_ACH_INJECTION.format(side="")
        pharmacology = _GRAB_ACH_BATH
    else:
        surgery = nwb_meta["surgery"] + _GRAB_ACH_INJECTION.format(side=" on the side of the lesion")
        if subject_entry.category == "6-OHDA":
            # The general text ends with the levodopa treatment, which these mice did not receive.
            pharmacology = nwb_meta["pharmacology"][: nwb_meta["pharmacology"].index("L-DOPA priming")] + _GRAB_ACH_BATH
        else:
            pharmacology = (
                nwb_meta["pharmacology"]
                + "Slices were made 24 to 48 hours after the last levodopa dose (off-state). "
                + _GRAB_ACH_BATH
            )

    nwbfile = NWBFile(
        session_description=(
            f"GRAB-ACh3.0 imaging (Figure 5): two-photon movies of acetylcholine release in the dorsolateral striatum, "
            f"evoked by electrical stimulation, in {len(regions)} regions of {len(slices)} slices from one "
            f"{condition} mouse. Each region was filmed in a series of trials while the drug in the bath changed; the "
            "trials table gives each trial's region, drug, stimulus and frames."
        ),
        experiment_description=(
            "Acetylcholine release in the striatum of unlesioned, 6-OHDA-lesioned and LID off-state mice, and its "
            "suppression by dopamine D2 receptors."
        ),
        identifier=str(uuid.uuid4()),
        session_start_time=session_start_time_naive.replace(tzinfo=TIMEZONE),
        session_id=_build_session_id(bundle),
        experimenter=nwb_meta["experimenter"],
        institution=nwb_meta["institution"],
        lab=nwb_meta["lab"],
        keywords=nwb_meta["keywords"],
        related_publications=nwb_meta["related_publications"],
        surgery=surgery,
        pharmacology=pharmacology,
        slices=nwb_meta["slices"],
    )
    nwbfile.subject = Subject(
        subject_id=subject_entry.subject_id,
        species=subject_meta["species"],
        strain=subject_entry.strain,
        genotype=subject_entry.genotype,
        sex=subject_meta["sex"],
        age=subject_meta["age"],
        description=(
            f"Wild-type C57BL/6J mouse, {_GRAB_ACH_SUBJECT[subject_entry.category]}; GRAB-ACh3.0 virus injected into "
            "the dorsolateral striatum, see surgery. The lab gave these mice no number; its Data Connections sheet "
            "identifies each by date, line and condition."
        ),
    )
    microscope = Device(
        name="BrukerUltima",
        description="Bruker Ultima In Vitro multiphoton microscope.",
        model=get_or_create_bruker_ultima_model(nwbfile),
    )
    nwbfile.add_device(microscope)
    nwbfile.add_device(
        Device(
            name="StimulatingElectrode",
            description=(
                "Concentric bipolar stimulating electrode (CBAPD75, FHC) about 200 µm ventral to the imaged region, "
                "giving 1 ms, 0.3 mA pulses (paper Methods). Prairie View timed each pulse with a 5 V trigger."
            ),
        )
    )
    ophys = nwbfile.create_processing_module(
        name="ophys", description="Prairie View's brightness of the whole frame for each region and channel."
    )
    segmentation = ImageSegmentation(name="ImageSegmentation")
    ophys.add(segmentation)
    fluorescence = Fluorescence(name="Fluorescence")
    ophys.add(fluorescence)
    references = Images(
        name="ImagesPrairieViewReferences",
        description=(
            "Reference images Prairie View saved with each recording, 16-bit, one per channel and trial. They show the "
            "field of view but are not frames of the trial, and consecutive trials often share the same image. "
            "Prairie View's 8-bit screen captures of the same images are left out."
        ),
    )
    nwbfile.add_acquisition(references)

    trial_rows = []  # keyword arguments of add_trial, per trial
    trial_references = []  # (index in trial_rows, number in recording order at the region, trials there, images)
    for slice_number, roi_number in regions:
        region_name = f"Slice{slice_number}Roi{roi_number}"
        trials = sorted(
            (
                _read_grab_ach_trial(sw, session_start_time_naive)
                for sw in bundle.sweeps_for((slice_number, roi_number))
            ),
            key=lambda trial: trial.start_time,
        )
        microns_per_pixel = {trial.microns_per_pixel for trial in trials}
        if len(microns_per_pixel) != 1:
            raise ValueError(f"The pixel size changes between the trials of slice {slice_number} ROI {roi_number}")
        imaging_plane = nwbfile.create_imaging_plane(
            name=f"ImagingPlane{region_name}",
            description=(
                f"Slice {slice_number}, ROI {roi_number}: a field of the dorsolateral striatum, with the stimulating "
                "electrode about 200 µm ventral to it (paper Methods)."
            ),
            optical_channel=[
                OpticalChannel(name=name, description=_GRAB_ACH_CHANNEL_DESCRIPTION[name], emission_lambda=float("nan"))
                for name in _GRAB_ACH_CHANNELS.values()
            ],
            device=microscope,
            excitation_lambda=920.0,
            indicator="GRAB-ACh3.0",
            location="Caudoputamen",
            grid_spacing=[microns_per_pixel.pop()] * 2,
            grid_spacing_unit="micrometers",
        )

        zooms = sorted({trial.zoom for trial in trials})
        for zoom in zooms:
            # Trials keep their number in recording order at the region across zooms.
            numbered = [(number, trial) for number, trial in enumerate(trials, start=1) if trial.zoom == zoom]
            group = [trial for _number, trial in numbered]
            if len({trial.frames["GRABACh3"].shape[1:] for trial in group}) != 1:
                raise ValueError(f"The frame size changes at zoom {zoom:g} in slice {slice_number} ROI {roi_number}")
            name = region_name + (f"Zoom{zoom:g}" if len(zooms) > 1 else "")
            frame_shape = group[0].frames["GRABACh3"].shape[1:]
            rate = 1 / group[0].frame_period
            layout = (
                f"the {len(group)} trials at this region joined in recording order, {frame_shape[0]} × {frame_shape[1]} "
                f"pixels at {rate:.2f} frames per second. Each frame's time is in timestamps, so the gaps between "
                "trials are kept."
            )
            if len(zooms) > 1:
                layout += (
                    f" This region was filmed at {len(zooms)} zooms, so its trials are split into {len(zooms)} movies; "
                    f"this one holds the trials at zoom {zoom:g}."
                )
            timestamps = np.concatenate([trial.frame_times for trial in group])
            movies, traces = {}, {}
            plane_segmentation = segmentation.create_plane_segmentation(
                name=f"PlaneSegmentation{name}",
                description=(
                    f"Prairie View's brightness-over-time region for slice {slice_number} ROI {roi_number}: the whole "
                    "frame, the same in every trial."
                ),
                imaging_plane=imaging_plane,
            )
            plane_segmentation.add_roi(image_mask=np.ones(frame_shape, dtype=bool))
            for channel in _GRAB_ACH_CHANNELS.values():
                what = _GRAB_ACH_CHANNEL_DESCRIPTION[channel]
                movies[channel] = TwoPhotonSeries(
                    name=f"TwoPhotonSeries{channel}{name}",
                    description=(
                        f"Movie of slice {slice_number} ROI {roi_number}, {what}: {layout} Each trial's frames are in "
                        f"the two_photon_{channel.lower()} column of the trials table."
                    ),
                    data=np.concatenate([trial.frames[channel] for trial in group]),
                    imaging_plane=imaging_plane,
                    unit="a.u.",
                    timestamps=timestamps if channel == "GRABACh3" else movies["GRABACh3"],
                )
                nwbfile.add_acquisition(movies[channel])
                traces[channel] = RoiResponseSeries(
                    name=f"RoiResponseSeries{channel}{name}",
                    description=(
                        f"Brightness of the whole frame of slice {slice_number} ROI {roi_number}, {what}, computed by "
                        f"Prairie View during recording (its botData file), one value per frame of "
                        f"TwoPhotonSeries{channel}{name}. It differs slightly from the mean of the saved frame. The "
                        "paper measured a 50 µm circle in the middle of the field instead (Methods)."
                    ),
                    data=np.concatenate([trial.brightness[channel] for trial in group]),
                    rois=plane_segmentation.create_roi_table_region(description="The whole frame.", region=[0]),
                    unit="a.u.",
                    timestamps=movies["GRABACh3"],
                )
                fluorescence.add_roi_response_series(traces[channel])

            first_frame = 0
            for number, trial in numbered:
                count = len(trial.frame_times)
                trial_rows.append(
                    dict(
                        start_time=float(trial.frame_times[0]),
                        stop_time=float(trial.frame_times[-1]) + trial.frame_period,
                        region=region_name,
                        drug=trial.drug,
                        stimulation=trial.stimulation,
                        stimulus_start_time=trial.stimulus_start_time,
                        **{
                            f"{kind}_{channel.lower()}": TimeSeriesReference(first_frame, count, series[channel])
                            for kind, series in (("two_photon", movies), ("roi_response", traces))
                            for channel in _GRAB_ACH_CHANNELS.values()
                        },
                    )
                )
                trial_references.append((len(trial_rows) - 1, number, len(trials), trial.references))
                first_frame += count

    nwbfile.add_trial_column(
        name="region",
        description=(
            "Region filmed, as in the series names: slice and ROI within the slice (Slice1Roi2). The lab named the only "
            "ROI of a slice 'A' in the May and June recordings; it is ROI 1 here."
        ),
    )
    nwbfile.add_trial_column(
        name="drug",
        description=(
            "Drug in the bath, from the recording's name, at the concentration the paper gives. 'control' is ACSF "
            "alone. Tetrodotoxin and acetylcholine trials find the sensor's lowest and highest brightness."
        ),
    )
    nwbfile.add_trial_column(
        name="stimulation",
        description=(
            "Electrical stimulation in this trial, from the recording's Prairie View VoltageOutput file. The paper "
            "gives 0.3 mA per pulse."
        ),
    )
    nwbfile.add_trial_column(
        name="stimulus_start_time",
        description=(
            "Onset of the first pulse in seconds from session_start_time, NaN without stimulation: the start of "
            "Prairie View's VoltageOutput plus its first-pulse delay (10 s, or 3 s on the May and June days)."
        ),
    )
    for kind in ("two_photon", "roi_response"):
        for channel in _GRAB_ACH_CHANNELS.values():
            series_type = "TwoPhotonSeries" if kind == "two_photon" else "RoiResponseSeries"
            nwbfile.add_trial_column(
                name=f"{kind}_{channel.lower()}",
                description=f"The trial's frames in the {series_type}{channel} series of its region.",
                col_cls=TimeSeriesReferenceVectorData,
            )

    # The trials table lists the trials of all regions in time order.
    order = sorted(range(len(trial_rows)), key=lambda index: trial_rows[index]["start_time"])
    for index in order:
        nwbfile.add_trial(**trial_rows[index])
    row_of = {index: row for row, index in enumerate(order)}
    for index, number, trials_in_region, images in trial_references:
        row = trial_rows[index]
        for channel, image in images.items():
            references.add_image(
                Image(
                    name=f"Image{channel}{row['region']}Trial{number}",
                    description=(
                        f"Prairie View's reference image, {_GRAB_ACH_CHANNEL_DESCRIPTION[channel]}, for "
                        f"{row['region']} trial {number} of {trials_in_region} in recording order (row "
                        f"{row_of[index]} of the trials table, start time {row['start_time']:.1f} s)."
                    ),
                    data=image,
                )
            )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    from neuroconv.tools.nwb_helpers import configure_and_write_nwbfile

    configure_and_write_nwbfile(nwbfile=nwbfile, nwbfile_path=output_path, backend="hdf5")
    return output_path


# Library module: `merge_one_mouse_day` is imported by the production runner
# (`scripts/full_per_mouse_day_merge.py`). No `main()` here — running the
# production runner is the right entry point for both single mouse-days and
# the whole dataset.
