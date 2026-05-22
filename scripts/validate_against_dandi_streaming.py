"""F-I cross-validation: DANDI (streaming, full cohort) vs local merged files.

Uses the streaming pattern from the existing reproduction notebook:
- DandiAPIClient(token=...).authenticate(token=...)
- asset.get_content_url(follow_redirects=1, strip_query=False)
- remfile.File(s3_url) -> h5py.File -> NWBHDF5IO

No downloads. Just remote h5py reads via HTTP range requests.
"""

from __future__ import annotations

import os
from collections import defaultdict
from pathlib import Path

import h5py
import numpy as np
import remfile
from dandi.dandiapi import DandiAPIClient
from dotenv import load_dotenv
from pynwb import NWBHDF5IO

REPO_ROOT = Path(__file__).resolve().parents[1]
load_dotenv(REPO_ROOT / ".env")

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
    last = -refractory_samples - 1
    n = 0
    for c in crossings:
        if c - last >= refractory_samples:
            n += 1
            last = c
    return n


def cond_from_session_id(ses: str) -> str | None:
    if "OffState++none" in ses:
        return "LID off-state"
    if "OnState++none" in ses:
        return "LID on-state"
    if "OnState++D1RaSch" in ses:
        return "LID on-state with SCH"
    return None


def stream_dandi_f1_soma_fi() -> dict[str, dict[str, dict[float, list[int]]]]:
    """Stream every F1 dSPN somatic file. Returns {cond: {cell_session_id: {step_pA: [ap_counts]}}}."""
    token = os.environ["DANDI_API_TOKEN"]
    client = DandiAPIClient(token=token)
    client.authenticate(token=token)
    ds = client.get_dandiset(DANDISET_ID, "draft")
    assets = [a for a in ds.get_assets() if "F1++SomExc++dSPN" in a.path]
    print(f"DANDI: streaming {len(assets)} F1 dSPN somatic files")

    out: dict[str, dict[str, dict[float, list[int]]]] = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
    for i, asset in enumerate(assets):
        ses = asset.path.split("ses-", 1)[1].split("_", 1)[0]
        cond = cond_from_session_id(ses)
        if cond is None:
            continue
        try:
            s3_url = asset.get_content_url(follow_redirects=1, strip_query=False)
            with remfile.File(s3_url) as f:
                with h5py.File(f, mode="r") as h5:
                    with NWBHDF5IO(file=h5, load_namespaces=True) as io:
                        nwb = io.read()
                        df = nwb.intracellular_recordings.to_dataframe()
                        for _, row in df.iterrows():
                            stim_pA = round(float(row[("intracellular_recordings", "stimulus_current_pA")]))
                            resp_ref = row[("responses", "response")]
                            ts = resp_ref.timeseries
                            idx_start = int(resp_ref.idx_start) if resp_ref.idx_start is not None else 0
                            count = (
                                int(resp_ref.count)
                                if resp_ref.count is not None and resp_ref.count > 0
                                else len(ts.data)
                            )
                            voltage = (
                                ts.data[idx_start : idx_start + count]
                                if (idx_start > 0 or count < len(ts.data))
                                else ts.data[:]
                            )
                            rate = float(ts.rate)
                            step_s = int(STEP_START_MS * 1e-3 * rate)
                            step_e = min(int(STEP_END_MS * 1e-3 * rate), len(voltage))
                            if step_s >= step_e:
                                continue
                            window = np.asarray(voltage[step_s:step_e])
                            n_aps = count_aps(window, rate)
                            out[cond][ses][float(stim_pA)].append(n_aps)
            if (i + 1) % 10 == 0:
                print(f"  streamed {i+1}/{len(assets)}")
        except Exception as e:
            print(f"  FAIL {asset.path[:90]}: {e}")
    return out


def local_f1_soma_fi() -> dict[str, dict[tuple, dict[float, list[int]]]]:
    out: dict[str, dict[tuple, dict[float, list[int]]]] = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
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
                stim_pA = round(float(row[("intracellular_recordings", "stimulus_current_pA")]))
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
                out[cond][cell_key][float(stim_pA)].append(n_aps)
    return out


def per_step_means(data, collapse_sweeps_first: bool):
    """data: {cell_key: {step: [counts]}} -> {step: list of per-cell means}."""
    per_step: dict[float, list[float]] = defaultdict(list)
    for cell, by_step in data.items():
        for step, counts in by_step.items():
            if collapse_sweeps_first:
                per_step[step].append(float(np.mean(counts)))
            else:
                per_step[step].extend(counts)
    return per_step


def main():
    print("=== Streaming DANDI F1 dSPN somatic ===")
    dandi = stream_dandi_f1_soma_fi()

    print()
    print("=== Reading local merged F1 ===")
    local = local_f1_soma_fi()

    target_steps = [-100, -60, -20, 20, 60, 100, 140, 180, 220, 260, 300]
    print()
    print(f"{'Condition':30s}  {'Source':15s}  " + "  ".join(f"{s:>5d}pA" for s in target_steps))
    print("-" * 130)
    for cond in ["LID off-state", "LID on-state", "LID on-state with SCH"]:
        d_per_step = per_step_means(dandi.get(cond, {}), collapse_sweeps_first=True)
        l_per_step = per_step_means(local.get(cond, {}), collapse_sweeps_first=True)
        n_d = len(dandi.get(cond, {}))
        n_l = len(local.get(cond, {}))

        def fmt(per_step):
            vals = []
            for s in target_steps:
                v = per_step.get(float(s), [])
                vals.append(f"{np.mean(v):6.2f}" if v else "    -")
            return "  ".join(vals)

        print(f"{cond:30s}  DANDI (n={n_d:<3d})  {fmt(d_per_step)}")
        print(f"{cond:30s}  LOCAL (n={n_l:<3d})  {fmt(l_per_step)}")
        # Compute per-step difference (LOCAL - DANDI)
        diffs = []
        for s in target_steps:
            dv = d_per_step.get(float(s), [])
            lv = l_per_step.get(float(s), [])
            if dv and lv:
                diffs.append(np.mean(lv) - np.mean(dv))
            else:
                diffs.append(float("nan"))
        diff_str = "  ".join(f"{d:+6.2f}" if not np.isnan(d) else "    -" for d in diffs)
        print(f"{cond:30s}  {'DIFF (L-D)':15s}  {diff_str}")
        print()


if __name__ == "__main__":
    main()
