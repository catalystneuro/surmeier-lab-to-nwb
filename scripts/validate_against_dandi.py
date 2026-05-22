"""Cross-validate: run the same F-I curve analysis on:
  (1) the OLD DANDI archive files (per-modality, per-condition session structure), AND
  (2) our NEW local merged per-mouse-day files

For F1 dSPN somatic excitability (Fig 1F). If both give the same SCH-reversal
anomaly, the effect is in the data itself, not our merger. If they differ
substantially, we have a bug.
"""

from __future__ import annotations

import os
from collections import defaultdict
from pathlib import Path

import h5py
import numpy as np
import remfile
from dandi.dandiapi import DandiAPIClient
from pynwb import NWBHDF5IO

REPO_ROOT = Path(__file__).resolve().parents[1]


def _load_env() -> None:
    env_path = REPO_ROOT / ".env"
    if env_path.exists():
        for line in env_path.read_text().splitlines():
            line = line.strip()
            if not line or line.startswith("#") or "=" not in line:
                continue
            k, v = line.split("=", 1)
            os.environ.setdefault(k.strip(), v.strip())


_load_env()
DANDISET_ID = "001538"
STEP_START_MS = 200.0
STEP_END_MS = 700.0
DVDT_THRESHOLD_MVPMS = 20.0


def count_aps(voltage_v: np.ndarray, rate_hz: float) -> int:
    voltage_mv = voltage_v * 1000.0
    sample_period_ms = 1000.0 / rate_hz
    dvdt = np.diff(voltage_mv) / sample_period_ms
    above = dvdt > DVDT_THRESHOLD_MVPMS
    crossings = np.where((~above[:-1]) & above[1:])[0]
    refractory_samples = int(2.0 / sample_period_ms)
    last_crossing = -refractory_samples - 1
    n = 0
    for c in crossings:
        if c - last_crossing >= refractory_samples:
            n += 1
            last_crossing = c
    return n


# ---- DANDI side ------------------------------------------------------------


def extract_dandi_fi():
    """Pull F1 dSPN somatic files from DANDI, compute per-cell AP counts per step."""
    token = os.environ.get("DANDI_API_TOKEN")
    if not token:
        raise RuntimeError("DANDI_API_TOKEN not set")

    # Per-cell counts keyed by (condition, cell_session_id) -> {step: ap_count}
    out: dict[str, dict[str, dict[float, int]]] = defaultdict(lambda: defaultdict(dict))
    with DandiAPIClient(token=token) as client:
        ds = client.get_dandiset(DANDISET_ID, "draft")
        assets = list(ds.get_assets())
    f1_soma_assets = [a for a in assets if "F1++SomExc++dSPN" in a.path]
    print(f"DANDI: found {len(f1_soma_assets)} F1 dSPN somatic assets")
    for asset in f1_soma_assets:
        # Parse condition from session_id
        # Path looks like sub-X/sub-X_ses-F1++SomExc++dSPN++<state>++<pharm>++WT++<ts>_icephys.nwb
        ses = asset.path.split("ses-", 1)[1].split("_", 1)[0]
        parts = ses.split("++")
        if len(parts) < 5:
            continue
        state, pharm = parts[3], parts[4]
        if state == "OffState" and pharm == "none":
            condition = "LID off-state"
        elif state == "OnState" and pharm == "none":
            condition = "LID on-state"
        elif state == "OnState" and pharm == "D1RaSch":
            condition = "LID on-state with SCH"
        else:
            continue
        s3_url = asset.get_content_url(follow_redirects=1, strip_query=True)
        try:
            with remfile.File(s3_url) as f, h5py.File(f, "r") as h5:
                # Find current-clamp series in acquisition
                acq = h5["acquisition"]
                for name in acq:
                    obj = acq[name]
                    if obj.attrs.get("neurodata_type", b"").decode() != "CurrentClampSeries":
                        continue
                    # Find the matching stimulus to get step amplitude
                    stim_pA = None
                    stim_grp = h5.get(f"stimulus/presentation/{name}_Stimulus") or h5.get(
                        f"stimulus/presentation/{name}"
                    )
                    if stim_grp is not None and "data" in stim_grp:
                        stim_data = stim_grp["data"][...]
                        # Use the peak of the stimulus during the step window as the amplitude
                        rate = float(obj.attrs.get("rate", 10000.0))
                        step_s = int(STEP_START_MS * 1e-3 * rate)
                        step_e = min(int(STEP_END_MS * 1e-3 * rate), len(stim_data))
                        if step_s < step_e:
                            # Stimulus is in amperes; convert to pA
                            step_window = stim_data[step_s:step_e]
                            peak_a = step_window[np.argmax(np.abs(step_window))]
                            stim_pA = round(peak_a * 1e12)
                    if stim_pA is None:
                        # Fall back: try to parse from series name (some old files encoded it)
                        continue
                    rate = float(obj.attrs.get("rate", 10000.0))
                    starting_time = float(obj.attrs.get("starting_time", 0.0))
                    voltage = obj["data"][...]
                    step_s = int(STEP_START_MS * 1e-3 * rate)
                    step_e = min(int(STEP_END_MS * 1e-3 * rate), len(voltage))
                    if step_s >= step_e:
                        continue
                    # Apply gain conversion from data attributes
                    conversion = float(obj["data"].attrs.get("conversion", 1.0))
                    voltage_v = voltage[step_s:step_e].astype(np.float64) * conversion
                    n_aps = count_aps(voltage_v, rate)
                    # cell_session_id distinguishes per-cell in DANDI
                    out[condition][ses][float(stim_pA)] = n_aps
        except Exception as e:
            print(f"  ERROR opening {asset.path}: {e}")
    return out


# ---- Local merged side ---------------------------------------------------


def extract_local_fi():
    """Same as DANDI but reads from local merged F1 files."""
    out: dict[str, dict[tuple[str, int], dict[float, list[int]]]] = defaultdict(
        lambda: defaultdict(lambda: defaultdict(list))
    )
    for nwb_path in sorted((REPO_ROOT / "nwb_files" / "merged" / "F1").glob("*.nwb")):
        with NWBHDF5IO(nwb_path, "r") as io:
            nwb = io.read()
            df = nwb.intracellular_recordings.to_dataframe()
            for _, row in df.iterrows():
                if row[("intracellular_recordings", "dendrite_type")] != "Soma":
                    continue
                cond = row[("intracellular_recordings", "condition_subfolder")]
                if cond not in ("LID off-state", "LID on-state", "LID on-state with SCH"):
                    continue
                stim_pA = float(row[("intracellular_recordings", "stimulus_current_pA")])
                resp_ref = row[("responses", "response")]
                ts = resp_ref.timeseries
                idx_start = int(resp_ref.idx_start)
                count = int(resp_ref.count)
                voltage = ts.data[idx_start : idx_start + count]
                rate = float(ts.rate)
                step_s = int(STEP_START_MS * 1e-3 * rate)
                step_e = min(int(STEP_END_MS * 1e-3 * rate), len(voltage))
                if step_s >= step_e:
                    continue
                window = np.asarray(voltage[step_s:step_e])
                n_aps = count_aps(window, rate)
                cell_id = int(row[("intracellular_recordings", "cell_id")])
                cell_key = (nwb.subject.subject_id, cell_id)
                out[cond][cell_key][stim_pA].append(n_aps)
    return out


def main():
    print("=== Extracting from DANDI (old per-modality files) ===")
    dandi = extract_dandi_fi()
    print(f"DANDI cells per condition: {[(c, len(cells)) for c, cells in dandi.items()]}")

    print()
    print("=== Extracting from local merged files (per-mouse-day) ===")
    local = extract_local_fi()
    print(f"Local cells per condition: {[(c, len(cells)) for c, cells in local.items()]}")

    print()
    print("=== F-I comparison at key current steps ===")
    steps_to_show = [100, 140, 180, 220, 260, 300]
    print(f"{'Condition':30s}  {'Source':12s}  " + "  ".join(f"{s:>5d}pA" for s in steps_to_show))

    for cond in ["LID off-state", "LID on-state", "LID on-state with SCH"]:
        # DANDI side: per-step mean across cells
        dandi_per_step = defaultdict(list)
        for ses, by_step in dandi.get(cond, {}).items():
            for step, n in by_step.items():
                dandi_per_step[step].append(n)
        # Local side: collapse per-cell first
        local_per_step = defaultdict(list)
        for cell_key, by_step in local.get(cond, {}).items():
            for step, counts in by_step.items():
                local_per_step[step].append(int(np.mean(counts)))
        dandi_means = [
            np.mean(dandi_per_step.get(float(s), [])) if dandi_per_step.get(float(s)) else float("nan")
            for s in steps_to_show
        ]
        local_means = [
            np.mean(local_per_step.get(float(s), [])) if local_per_step.get(float(s)) else float("nan")
            for s in steps_to_show
        ]
        dandi_n = max(len(dandi_per_step.get(float(s), [])) for s in steps_to_show)
        local_n = max(len(local_per_step.get(float(s), [])) for s in steps_to_show)
        d_str = "  ".join(f"{m:6.1f}" if not np.isnan(m) else "    -" for m in dandi_means)
        l_str = "  ".join(f"{m:6.1f}" if not np.isnan(m) else "    -" for m in local_means)
        print(f"{cond:30s}  DANDI (n={dandi_n:<3d})  {d_str}")
        print(f"{cond:30s}  LOCAL (n={local_n:<3d})  {l_str}")
        print()


if __name__ == "__main__":
    main()
