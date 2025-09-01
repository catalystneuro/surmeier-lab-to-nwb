#!/usr/bin/env python3
"""
Figure 2G–H oEPSC Reproduction (Streaming from DANDI)

Streams F2 dSPN Sr2+-oEPSC NWB files directly from DANDI, detects events using a
MAD-based method with near-event merging, and generates cumulative probability and
subject-level (per-file mean) box plots.

Requirements:
- DANDI_API_TOKEN set in environment (dataset is draft)
- pip: dandi, pynwb, remfile, h5py, numpy, pandas, matplotlib, scipy, python-dotenv

Outputs:
- analysis_outputs/figure_2gh_from_dandi/figure_2gh_cumulative_probability.png
- analysis_outputs/figure_2gh_from_dandi/figure_2gh_amplitude_boxplot.png
- analysis_outputs/figure_2gh_from_dandi/figure_2gh_reproduction_report.md
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Dict, List, Tuple

import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import remfile
from dandi.dandiapi import DandiAPIClient
from dotenv import load_dotenv
from pynwb import NWBHDF5IO
from scipy import stats


def _get_ts_slice_in_si_units(ts, idx_start: int, count: int) -> Tuple[np.ndarray, np.ndarray]:
    data = np.asarray(ts.data[idx_start : idx_start + count], dtype=float)
    if ts.timestamps is not None:
        t = np.asarray(ts.timestamps[idx_start : idx_start + count], dtype=float)
    else:
        rate = float(ts.rate)
        start = float(ts.starting_time) + idx_start / rate
        t = start + np.arange(count, dtype=float) / rate
    return data, t


def _group_contiguous_indices(indices: np.ndarray) -> List[np.ndarray]:
    if indices.size == 0:
        return []
    splits = np.where(np.diff(indices) > 1)[0] + 1
    return np.split(indices, splits)


def merge_nearby_events(
    event_times_ms: np.ndarray,
    event_amplitudes_pA: np.ndarray,
    event_polarity: np.ndarray,
    merge_distance_ms: float = 5.0,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    if event_times_ms.size == 0:
        return event_times_ms, event_amplitudes_pA, event_polarity
    order = np.argsort(event_times_ms)
    t = np.asarray(event_times_ms)[order]
    a = np.asarray(event_amplitudes_pA)[order]
    p = np.asarray(event_polarity, dtype=int)[order]
    gaps = np.diff(t)
    splits = np.where(gaps > merge_distance_ms)[0] + 1
    groups = np.split(np.arange(t.size), splits)
    kept_idx = []
    for grp in groups:
        j = grp[np.argmax(a[grp])]  # keep max amplitude within window
        kept_idx.append(int(j))
    kept_idx = np.asarray(kept_idx, dtype=int)
    return t[kept_idx], a[kept_idx], p[kept_idx]


def detect_events_mad(
    data_pA: np.ndarray,
    timestamps_ms: np.ndarray,
    detection_start_ms: float,
    detection_stop_ms: float,
    detection_window_shift_ms: float = 100.0,
    threshold_k: float = 5.0,
    min_inter_event_ms: float = 5.0,
) -> Dict[str, np.ndarray]:
    det_start = detection_start_ms + detection_window_shift_ms
    det_stop = detection_stop_ms
    mask = (timestamps_ms >= det_start) & (timestamps_ms <= det_stop)
    if not np.any(mask):
        return {"times_ms": np.array([]), "amps_pA": np.array([]), "polarity": np.array([])}
    x = data_pA[mask]
    t = timestamps_ms[mask]
    med = float(np.median(x))
    mad = float(np.median(np.abs(x - med)))
    mad_std = mad * 1.4826
    pos_thr = med + threshold_k * mad_std
    neg_thr = med - threshold_k * mad_std
    pos_idx = np.flatnonzero(x > pos_thr)
    neg_idx = np.flatnonzero(x < neg_thr)
    times, amps, pol = [], [], []
    for grp in _group_contiguous_indices(pos_idx):
        j = grp[np.argmax(x[grp])]
        times.append(t[j])
        amps.append(x[j] - med)
        pol.append(1)
    for grp in _group_contiguous_indices(neg_idx):
        j = grp[np.argmin(x[grp])]
        times.append(t[j])
        amps.append(med - x[j])
        pol.append(-1)
    if not times:
        return {"times_ms": np.array([]), "amps_pA": np.array([]), "polarity": np.array([])}
    order = np.argsort(times)
    t = np.asarray(times)[order]
    a = np.asarray(amps)[order]
    p = np.asarray(pol, dtype=int)[order]
    if t.size and min_inter_event_ms > 0:
        t, a, p = merge_nearby_events(t, a, p, merge_distance_ms=min_inter_event_ms)
    return {"times_ms": t, "amps_pA": a, "polarity": p}


def main() -> None:
    load_dotenv()
    token = os.getenv("DANDI_API_TOKEN")
    if not token:
        raise SystemExit("DANDI_API_TOKEN not set. Please export it or add to .env")

    client = DandiAPIClient(token=token)
    client.authenticate(token=token)
    dandiset = client.get_dandiset("001538", "draft")
    assets = list(dandiset.get_assets())

    def get_session_id(asset_path: str) -> str:
        bottom = asset_path.split("/")[1]
        ses = bottom.split("_")[1]
        return ses.split("-")[1]

    def is_f2_oepsc_dspn(session_id: str) -> bool:
        parts = session_id.split("++")
        return len(parts) >= 3 and parts[0] == "F2" and parts[1] == "oEPSC" and parts[2] == "dSPN"

    oepsc_assets = [a for a in assets if is_f2_oepsc_dspn(get_session_id(a.path))]
    if not oepsc_assets:
        raise SystemExit("No F2 oEPSC dSPN assets found in DANDI.")

    print(f"Found {len(oepsc_assets)} F2 oEPSC assets")

    rows: List[Dict] = []
    for a in oepsc_assets:
        try:
            s3_url = a.get_content_url(follow_redirects=1, strip_query=False)
            rf = remfile.File(s3_url)
            f = h5py.File(rf, mode="r")
            io = NWBHDF5IO(file=f, load_namespaces=True)
            nwb = io.read()
        except Exception as e:
            print(f"Skipping {a.path} due to open/read error: {e}")
            continue

        # Intervals table for timing (relative offsets per sweep)
        opto = nwb.intervals["optogenetic_epochs_table"].to_dataframe()
        pre = opto[opto.stage_name == "pre_stimulation"].iloc[0]
        det = opto[opto.stage_name == "detection"].iloc[0]
        stim = opto[opto.stage_name == "stimulation"].iloc[0]
        det0_start = float(det.start_time) - float(pre.start_time)
        det0_stop = float(det.stop_time) - float(pre.start_time)
        stim0_start = float(stim.start_time) - float(pre.start_time)

        ice = nwb.get_intracellular_recordings()
        for i in range(len(ice.id)):
            row = ice[i]
            resp_ref = row[("responses", "response")].iloc[0]
            ts = resp_ref.timeseries
            data_A, t_s = _get_ts_slice_in_si_units(ts, resp_ref.idx_start, resp_ref.count)
            data_pA = data_A * 1e12
            t_ms = (t_s - t_s[0]) * 1000.0
            detres = detect_events_mad(
                data_pA=data_pA,
                timestamps_ms=t_ms,
                detection_start_ms=det0_start * 1000.0,
                detection_stop_ms=det0_stop * 1000.0,
                detection_window_shift_ms=100.0,
                threshold_k=5.0,
                min_inter_event_ms=5.0,
            )
            condition = "LID on-state" if "OnState" in a.path else "LID off-state"
            for t_ev, amp, pol in zip(detres["times_ms"], detres["amps_pA"], detres["polarity"]):
                rows.append(
                    dict(
                        nwb_file=a.path.split("/")[-1],
                        condition=condition,
                        oepsc_amplitude_pA=float(amp),
                        event_time_ms=float(t_ev),
                        event_polarity=int(pol),
                        stim_start_ms=stim0_start * 1000.0,
                        detection_start_ms=det0_start * 1000.0,
                        detection_stop_ms=det0_stop * 1000.0,
                    )
                )

    df = pd.DataFrame(rows)
    if df.empty:
        raise SystemExit("No events detected from DANDI assets.")
    print(f"Detected {len(df)} events across {df['nwb_file'].nunique()} files")

    out_dir = Path("analysis_outputs/figure_2gh_from_dandi")
    out_dir.mkdir(parents=True, exist_ok=True)

    # Cumulative probability (all events)
    fig, ax = plt.subplots(1, 1, figsize=(6, 4))
    for label, color in [("LID off-state", "black"), ("LID on-state", "gray")]:
        vals = df[df.condition == label]["oepsc_amplitude_pA"].values
        if vals.size == 0:
            continue
        s = np.sort(vals)
        cp = np.arange(1, s.size + 1) / s.size * 100
        ax.plot(s, cp, color=color, lw=1.5, label=label.split()[1])
    ax.set_xlabel("oEPSC amplitude (pA)")
    ax.set_ylabel("Cumulative Probability (%)")
    ax.set_xlim(0, 80)
    ax.set_ylim(0, 100)
    ax.legend(frameon=False, loc="upper left")
    ax.set_title("Figure 2G: dSPN oEPSC Cumulative Probability (streamed)")
    cdf_path = out_dir / "figure_2gh_cumulative_probability.png"
    fig.savefig(cdf_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)

    # Subject-level (per-file mean) box plot
    per_file = df.groupby(["condition", "nwb_file"])["oepsc_amplitude_pA"].mean().reset_index()
    order = ["LID off-state", "LID on-state"]
    labels = ["off-state", "on-state"]
    plot_data = [per_file[per_file.condition == c]["oepsc_amplitude_pA"].values for c in order]
    fig, ax = plt.subplots(1, 1, figsize=(4, 5))
    ax.boxplot(
        plot_data,
        labels=labels,
        patch_artist=True,
        boxprops=dict(facecolor="white", color="black", linewidth=1),
        whiskerprops=dict(color="black", linewidth=1),
        capprops=dict(color="black", linewidth=1),
        medianprops=dict(color="black", linewidth=1.5),
        flierprops=dict(marker="o", markerfacecolor="gray", markersize=3, markeredgecolor="black", alpha=0.7),
    )
    for i, data in enumerate(plot_data):
        x = np.random.normal(i + 1, 0.04, len(data))
        ax.scatter(x, data, color="gray", alpha=0.8, s=20, zorder=3)
    ax.set_ylabel("oEPSC amplitude (pA)")
    ax.set_ylim(0, 25)
    ax.set_yticks(np.arange(0, 26, 5))
    ax.set_title("Figure 2H: dSPN oEPSC Amplitude (streamed)")
    off = plot_data[0]
    on = plot_data[1]
    if off.size > 1 and on.size > 1:
        stat, p = stats.ttest_ind(off, on, equal_var=False)
        y = 23
        ax.text(
            1.5,
            y,
            "***" if p < 1e-3 else "**" if p < 1e-2 else "*" if p < 0.05 else "n.s.",
            ha="center",
            va="bottom",
            fontsize=12,
        )
    box_path = out_dir / "figure_2gh_amplitude_boxplot.png"
    fig.savefig(box_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)

    # Minimal report
    md = out_dir / "figure_2gh_reproduction_report.md"
    with md.open("w") as f:
        f.write("# Figure 2G–H oEPSC Reproduction (Streaming)\n\n")
        f.write(f"Detected events: {len(df)} across {df['nwb_file'].nunique()} files\n\n")
        for c in order:
            vals = df[df.condition == c]["oepsc_amplitude_pA"].values
            if vals.size:
                f.write(f"- {c}: n={vals.size}, mean={vals.mean():.2f} pA, median={np.median(vals):.2f} pA\n")
        f.write("\n")
        f.write(f"- CDF: {cdf_path}\n")
        f.write(f"- Box: {box_path}\n")

    print("Saved:", cdf_path)
    print("Saved:", box_path)
    print("Saved:", md)


if __name__ == "__main__":
    main()
