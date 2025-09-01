#!/usr/bin/env python3
"""
Figure 1E Somatic Excitability Reproduction (Streaming from DANDI)

This script reproduces Figure 1E F–I analysis by streaming NWB files directly
from DANDI (no local downloads). It mirrors the robust logic from the local
reproduction script but pulls data from the dandiset.

Requirements:
- Environment variable DANDI_API_TOKEN (dataset is currently draft)
- Packages: dandi, pynwb, remfile, h5py, python-dotenv, numpy, pandas, matplotlib

Outputs:
- analysis_outputs/figure_1e_from_dandi/figure_1e_fi_curves.png
- analysis_outputs/figure_1e_from_dandi/figure_1e_rheobase_comparison.png
- analysis_outputs/figure_1e_from_dandi/figure_1e_reproduction_report.md
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


def count_action_potentials(voltage_trace_mV: np.ndarray, timestamps_s: np.ndarray, threshold_mV: float = 0.0) -> int:
    """Threshold-crossing spike count within a 500 ms stimulus window.

    Matches approach in the local reproduction script; uses 200–700 ms from sweep start.
    """
    if voltage_trace_mV.size == 0:
        return 0
    if timestamps_s[-1] - timestamps_s[0] < 0.8:
        return 0
    start_t = timestamps_s[0] + 0.2
    end_t = timestamps_s[0] + 0.7
    i0 = np.searchsorted(timestamps_s, start_t)
    i1 = np.searchsorted(timestamps_s, end_t)
    if i0 >= i1 or i1 > voltage_trace_mV.size:
        i0 = voltage_trace_mV.size // 4
        i1 = 3 * voltage_trace_mV.size // 4
    x = voltage_trace_mV[i0:i1]
    if x.size == 0:
        return 0
    spikes = 0
    i = 0
    while i < x.size:
        if x[i] > threshold_mV:
            spikes += 1
            while i < x.size and x[i] > threshold_mV:
                i += 1
        else:
            i += 1
    return spikes


def get_session_id(asset_path: str) -> str:
    bottom = asset_path.split("/")[1]
    ses = bottom.split("_")[1]
    return ses.split("-")[1]


def get_fig(session_id: str) -> str:
    return session_id.split("++")[0]


def get_meas(session_id: str) -> str:
    return session_id.split("++")[1]


def is_f1_somexc(session_id: str) -> bool:
    return get_fig(session_id) == "F1" and get_meas(session_id) == "SomExc"


def calculate_rheobase(current_steps: List[float], spike_counts: List[int]) -> float:
    for current, spikes in zip(current_steps, spike_counts):
        if spikes >= 1:
            return current
    return np.nan


def create_fi_curves(df: pd.DataFrame) -> Tuple[plt.Figure, plt.Axes]:
    fig, ax = plt.subplots(1, 1, figsize=(6, 4))
    condition_styles = {
        "LID off-state": {"color": "black", "marker": "o", "linestyle": "-", "label": "off-state"},
        "LID on-state": {"color": "black", "marker": "s", "linestyle": "-", "label": "on-state"},
        "LID on-state with SCH": {
            "color": "gray",
            "marker": "^",
            "linestyle": "--",
            "label": "on-state+D1R\nantagonist",
        },
    }

    def safe_mean_sem(values: np.ndarray) -> Tuple[float, float]:
        if len(values) == 0:
            return np.nan, 0.0
        if len(values) == 1:
            return float(values[0]), 0.0
        return float(np.mean(values)), float(np.std(values, ddof=1) / np.sqrt(len(values)))

    for condition in ["LID off-state", "LID on-state", "LID on-state with SCH"]:
        if condition not in df["condition"].unique():
            continue
        condition_data = df[df["condition"] == condition]
        summary_rows = []
        for current, group in condition_data.groupby("current_pA"):
            m, s = safe_mean_sem(group["spike_count"].values)
            summary_rows.append({"current_pA": current, "mean_spikes": m, "sem_spikes": s})
        summary = pd.DataFrame(summary_rows).sort_values("current_pA")
        summary = summary[(summary["current_pA"] >= 0) & (summary["current_pA"] <= 300)]
        style = condition_styles[condition]
        ax.errorbar(
            summary["current_pA"],
            summary["mean_spikes"],
            yerr=summary["sem_spikes"],
            marker=style["marker"],
            color=style["color"],
            linestyle=style["linestyle"],
            linewidth=1.5,
            markersize=4,
            capsize=3,
            capthick=1,
            label=style["label"],
            markerfacecolor="white" if condition == "LID off-state" else style["color"],
            markeredgecolor=style["color"],
            markeredgewidth=1,
        )

    ax.set_xlabel("current (pA)", fontsize=12)
    ax.set_ylabel("number of APs", fontsize=12)
    ax.set_xlim(0, 300)
    ax.set_ylim(0, 18)
    ax.set_xticks([0, 100, 200, 300])
    ax.set_yticks([0, 5, 10, 15])
    ax.legend(loc="upper left", frameon=False, fontsize=10)
    ax.set_title("Figure 1E: dSPN Somatic Excitability (F-I Curves)", fontsize=14, pad=20)
    return fig, ax


def create_rheobase_comparison(cell_stats: pd.DataFrame) -> Tuple[plt.Figure, plt.Axes]:
    fig, ax = plt.subplots(1, 1, figsize=(4, 5))
    order = ["LID off-state", "LID on-state", "LID on-state with SCH"]
    labels = ["off-state", "on-state", "on+SCH"]
    plot_data = [cell_stats[cell_stats["condition"] == c]["rheobase_pA"].dropna().values for c in order]
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
    ax.set_ylabel("rheobase (pA)", fontsize=12)
    ax.set_ylim(0, 300)
    ax.set_title("Figure 1E: Rheobase Comparison", fontsize=14, pad=20)
    return fig, ax


def generate_markdown_report(cell_stats: pd.DataFrame, output_dir: Path) -> str:
    from datetime import datetime

    ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    valid = cell_stats.dropna(subset=["rheobase_pA"]).copy()
    md = [
        f"# Figure 1E Somatic Excitability Analysis Report\n",
        f"**Generated:** {ts}\n",
        f"**Data Source:** DANDI streaming (NWB icephys)\n",
        f"**Analysis:** F-I curves and rheobase measurements for dSPN somatic excitability\n\n",
        "## Executive Summary\n\n",
        f"- **Total Cells Analyzed:** {len(valid)}\n",
        f"- **Conditions:** {valid['condition'].nunique()}\n",
        f"- **Total Recordings:** {int(valid['n_recordings'].sum())}\n\n",
        "## Results by Condition\n\n",
    ]
    for cond in ["LID off-state", "LID on-state", "LID on-state with SCH"]:
        sub = valid[valid["condition"] == cond]
        if sub.empty:
            continue
        r = sub["rheobase_pA"].dropna()
        md += [
            f"### {cond}\n\n",
            f"- **Sample Size:** {len(sub)} cells\n",
            f"- **Rheobase:** {r.mean():.1f} ± {r.sem():.1f} pA (mean ± SEM)\n",
            f"- **Range:** {r.min():.1f} - {r.max():.1f} pA\n\n",
        ]
    md += [
        "## Figures\n\n",
        "### F-I Curves\n",
        "![F-I Curves](figure_1e_fi_curves.png)\n\n",
        "### Rheobase Comparison\n",
        "![Rheobase Comparison](figure_1e_rheobase_comparison.png)\n\n",
        "---\n\n",
        "*Report generated by Figure 1E Somatic Excitability Streaming Script*\n",
    ]
    report_path = output_dir / "figure_1e_reproduction_report.md"
    report_path.write_text("".join(md))
    return str(report_path)


def main() -> None:
    load_dotenv()
    token = os.getenv("DANDI_API_TOKEN")
    if not token:
        raise SystemExit("DANDI_API_TOKEN not set. Please export it or add to .env")

    client = DandiAPIClient(token=token)
    client.authenticate(token=token)
    dandiset = client.get_dandiset("001538", "draft")
    assets = list(dandiset.get_assets())

    som_assets = [a for a in assets if is_f1_somexc(get_session_id(a.path))]
    if not som_assets:
        raise SystemExit("No F1 SomExc assets found in the dandiset.")

    print(f"Found {len(som_assets)} F1 SomExc assets")

    # Collect recording-level rows across all assets
    all_rows: List[Dict] = []

    for a in som_assets:
        try:
            s3_url = a.get_content_url(follow_redirects=1, strip_query=False)
            rf = remfile.File(s3_url)
            f = h5py.File(rf, mode="r")
            io = NWBHDF5IO(file=f, load_namespaces=True)
            nwb = io.read()
        except Exception as e:
            print(f"Skipping {a.path} due to open/read error: {e}")
            continue

        try:
            # Prefer using the intracellular_recordings table with custom columns, if present
            if hasattr(nwb, "intracellular_recordings"):
                rec_df = nwb.intracellular_recordings.to_dataframe()
            else:
                rec_df = nwb.get_intracellular_recordings().to_dataframe()
        except Exception as e:
            print(f"No intracellular_recordings table for {a.path}: {e}")
            continue

        steps: List[float] = []
        spikes: List[int] = []

        # Iterate over rows; robustly obtain protocol_step and current value
        for _, row in rec_df.iterrows():
            # Protocol step
            protocol_step = row.get(("intracellular_recordings", "protocol_step"), None)
            # Current injected in pA (custom column may exist)
            current_pA = row.get(("intracellular_recordings", "stimulus_current_pA"), None)

            # Get the voltage response TimeSeries by two strategies
            ts = None

            # Strategy A: if acquisition series named by protocol_step exists
            if protocol_step is not None:
                series_name = (
                    f"CurrentClampSeries{int(protocol_step):03d}"
                    if not str(protocol_step).startswith("CurrentClampSeries")
                    else str(protocol_step)
                )
                if series_name in nwb.acquisition:
                    ts = nwb.acquisition[series_name]

            # Strategy B: follow the reference in the table
            if ts is None:
                try:
                    response_ref = row[("responses", "response")]
                    if hasattr(response_ref, "iloc"):
                        response_ref = response_ref.iloc[0]
                    ts = response_ref.timeseries
                except Exception:
                    ts = None

            if ts is None:
                continue

            # Timestamps robustly: either explicit or derived from starting_time+rate
            if ts.timestamps is not None:
                t_s = np.asarray(ts.timestamps, dtype=float)
            else:
                rate = float(ts.rate)
                t_s = float(ts.starting_time) + np.arange(ts.data.shape[0], dtype=float) / rate

            # Voltage in mV
            v_mV = np.asarray(ts.data, dtype=float) * 1000.0

            # If current not provided as custom column, estimate from stimulus reference
            if current_pA is None or pd.isna(current_pA):
                try:
                    stim_ref = row[("stimuli", "stimulus")]
                    if hasattr(stim_ref, "iloc"):
                        stim_ref = stim_ref.iloc[0]
                    stim_ts = stim_ref.timeseries
                    if stim_ts is not None:
                        iA = np.asarray(stim_ts.data, dtype=float)
                        current_pA = float(np.median(iA * 1e12))
                except Exception:
                    current_pA = np.nan

            steps.append(float(current_pA) if current_pA is not None else np.nan)
            spikes.append(count_action_potentials(v_mV, t_s, threshold_mV=0.0))

        # Session/condition metadata from path
        session_id = get_session_id(a.path)
        parts = session_id.split("++")
        condition = "unknown"
        if len(parts) >= 5:
            state, pharm = parts[3], parts[4]
            if state == "OffState" and pharm == "none":
                condition = "LID off-state"
            elif state == "OnState" and pharm == "none":
                condition = "LID on-state"
            elif state == "OnState" and pharm == "D1RaSch":
                condition = "LID on-state with SCH"

        try:
            subject_id = nwb.subject.subject_id if nwb.subject is not None else session_id
        except Exception:
            subject_id = session_id

        rec_df = pd.DataFrame(dict(current_pA=steps, spike_count=spikes)).dropna().sort_values("current_pA")
        for _, r in rec_df.iterrows():
            all_rows.append(
                dict(
                    session_id=session_id,
                    subject_id=subject_id,
                    condition=condition,
                    current_pA=float(r.current_pA),
                    spike_count=int(r.spike_count),
                    nwb_file=a.path.split("/")[-1],
                )
            )
        print(f"Processed {a.path.split('/')[-1]} with {len(rec_df)} sweeps")

    # Build DF and outputs
    out_dir = Path("analysis_outputs/figure_1e_from_dandi")
    out_dir.mkdir(parents=True, exist_ok=True)
    df_all = pd.DataFrame(all_rows)
    if df_all.empty:
        raise SystemExit("No somatic excitability sweeps parsed from DANDI.")

    # Create F–I curves across conditions
    fig1, _ = create_fi_curves(df_all)
    fig1_path = out_dir / "figure_1e_fi_curves.png"
    fig1.savefig(fig1_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig1)
    print(f"Saved: {fig1_path}")

    # Cell-level stats (per session_id)
    cell_rows = []
    for (session_id, condition), g in df_all.groupby(["session_id", "condition"]):
        g = g.sort_values("current_pA")
        rheo = calculate_rheobase(g["current_pA"].tolist(), g["spike_count"].tolist())
        cell_rows.append(
            dict(
                session_id=session_id,
                condition=condition,
                subject_id=g["subject_id"].iloc[0],
                rheobase_pA=rheo,
                n_recordings=len(g),
            )
        )
    cell_stats = pd.DataFrame(cell_rows)

    # Rheobase comparison plot
    fig2, _ = create_rheobase_comparison(cell_stats)
    fig2_path = out_dir / "figure_1e_rheobase_comparison.png"
    fig2.savefig(fig2_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig2)
    print(f"Saved: {fig2_path}")

    # Report
    report_path = generate_markdown_report(cell_stats, out_dir)
    print(f"Saved: {report_path}")


if __name__ == "__main__":
    main()
