#!/usr/bin/env python3
"""
Figure 1I Dendritic Excitability Reproduction (Streaming from DANDI)

Streams F1 dSPN Dendritic Excitability NWB files from DANDI, computes proximal
and distal ΔG/R0 AUCs from Fluorescence ROIResponseSeries, and plots the distal/proximal
ratio per condition.

Requirements:
- DANDI_API_TOKEN set in environment (dataset is draft)
- pip: dandi, pynwb, remfile, h5py, numpy, pandas, matplotlib, python-dotenv

Outputs:
- analysis_outputs/figure_1i_from_dandi/figure_1i_dendritic_excitability.png
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Dict, List

import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import remfile
from dandi.dandiapi import DandiAPIClient
from dotenv import load_dotenv
from pynwb import NWBHDF5IO


def get_session_id(asset_path: str) -> str:
    bottom = asset_path.split("/")[1]
    ses = bottom.split("_")[1]
    return ses.split("-")[1]


def is_f1_dendexc_dspn(session_id: str) -> bool:
    parts = session_id.split("++")
    return len(parts) >= 3 and parts[0] == "F1" and parts[1] == "DendExc" and parts[2] == "dSPN"


def calc_delta_g_over_r0_auc(fluo4_series, alexa568_series) -> float:
    try:
        fluo4 = np.asarray(fluo4_series.data, dtype=float)
        alexa = np.asarray(alexa568_series.data, dtype=float)
        if fluo4.size != alexa.size or fluo4.size == 0:
            return np.nan
        # Background subtraction (constants from lab analysis)
        bg_pmt1, bg_pmt2 = 55.0, 161.0
        alexa_corr = alexa - bg_pmt1
        fluo4_corr = fluo4 - bg_pmt2
        # 6-point rolling average
        alexa_s = pd.Series(alexa_corr).rolling(6, min_periods=1).mean().values
        fluo4_s = pd.Series(fluo4_corr).rolling(6, min_periods=1).mean().values
        alexa_s = np.where(alexa_s <= 0, np.nan, alexa_s)
        # G/R0 baseline from 0.5 to 0.95 s
        if fluo4_series.timestamps is not None:
            ts = np.asarray(fluo4_series.timestamps, dtype=float)
            dt = ts[1] - ts[0] if ts.size > 1 else 0.00064
        else:
            dt = 0.00064
        i0 = min(int(0.5 / dt), fluo4_s.size - 100)
        i1 = min(int(0.95 / dt), fluo4_s.size - 50)
        j1 = min(int(1.5 / dt), fluo4_s.size - 1)
        if i0 >= i1 or i1 >= j1:
            return np.nan
        g0 = np.nanmean(fluo4_s[i0:i1])
        r0 = np.nanmean(alexa_s[i0:i1])
        if not np.isfinite(g0) or not np.isfinite(r0) or r0 <= 0:
            return np.nan
        delta_g_over_r0 = (fluo4_s - g0) / r0
        resp = delta_g_over_r0[i1:j1]
        resp = resp[np.isfinite(resp)]
        if resp.size == 0:
            return np.nan
        return float(np.sum(resp) * dt)
    except Exception:
        return np.nan


def main() -> None:
    load_dotenv()
    token = os.getenv("DANDI_API_TOKEN")
    if not token:
        raise SystemExit("DANDI_API_TOKEN not set. Please export it or add to .env")

    client = DandiAPIClient(token=token)
    client.authenticate(token=token)
    dandiset = client.get_dandiset("001538", "draft")
    assets = list(dandiset.get_assets())
    dend_assets = [a for a in assets if is_f1_dendexc_dspn(get_session_id(a.path))]
    if not dend_assets:
        raise SystemExit("No F1 DendExc dSPN assets found.")
    print(f"Found {len(dend_assets)} F1 DendExc assets")
    # Optional cap to avoid long runtimes; set MAX_ASSETS env to override
    try:
        max_assets = int(os.getenv("MAX_ASSETS", "20"))
    except Exception:
        max_assets = 20

    rows: List[Dict] = []
    for a in dend_assets[:max_assets]:
        try:
            s3 = a.get_content_url(follow_redirects=1, strip_query=False)
            rf = remfile.File(s3)
            f = h5py.File(rf, "r")
            io = NWBHDF5IO(file=f, load_namespaces=True)
            nwb = io.read()
        except Exception as e:
            print(f"Skipping {a.path} due to open/read error: {e}")
            continue

        # Determine condition from session id
        sid = nwb.session_id
        parts = sid.split("++")
        condition = "unknown"
        if len(parts) >= 5:
            state, pharm = parts[3], parts[4]
            if state == "OffState" and pharm == "none":
                condition = "LID off-state"
            elif state == "OnState" and pharm == "none":
                condition = "LID on-state"
            elif state == "OnState" and pharm == "D1RaSch":
                condition = "LID on-state with SCH"

        if "ophys" not in nwb.processing or "Fluorescence" not in nwb.processing["ophys"].data_interfaces:
            print(f"No Fluorescence interface in {a.path}")
            continue
        fl = nwb.processing["ophys"].data_interfaces["Fluorescence"]
        rrs = fl.roi_response_series
        # Build location-trial suffixes from Fluo4 names
        suffixes = [
            name.replace("RoiResponseSeriesFluo4", "")
            for name in rrs.keys()
            if name.startswith("RoiResponseSeriesFluo4")
        ]
        # Aggregate proximal and distal AUCs
        prox_aucs, dist_aucs = [], []
        for suf in suffixes:
            fluo4_name = f"RoiResponseSeriesFluo4{suf}"
            alexa_name = f"RoiResponseSeriesAlexa568{suf}"
            if fluo4_name in rrs and alexa_name in rrs:
                auc = calc_delta_g_over_r0_auc(rrs[fluo4_name], rrs[alexa_name])
                target = dist_aucs if "Distal" in suf else prox_aucs if "Proximal" in suf else None
                if target is not None and np.isfinite(auc):
                    target.append(auc)
        if prox_aucs and dist_aucs:
            idx = float(np.nanmean(dist_aucs) / np.nanmean(prox_aucs)) if np.nanmean(prox_aucs) != 0 else np.nan
            rows.append(dict(nwb_file=a.path.split("/")[-1], condition=condition, dend_ex_idx=idx))

    df = pd.DataFrame(rows)
    if df.empty:
        raise SystemExit("No dendritic excitability indices computed.")
    print(f"Computed indices for {len(df)} files")

    # Plot per-condition box/points
    out = Path("analysis_outputs/figure_1i_from_dandi")
    out.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(1, 1, figsize=(4, 5))
    order = ["LID off-state", "LID on-state", "LID on-state with SCH"]
    labels = ["off", "on", "on+SCH"]
    pdata = [df[df.condition == c]["dend_ex_idx"].values for c in order]
    ax.boxplot(
        pdata,
        labels=labels,
        patch_artist=True,
        boxprops=dict(facecolor="white", color="black", linewidth=1),
        whiskerprops=dict(color="black", linewidth=1),
        capprops=dict(color="black", linewidth=1),
        medianprops=dict(color="black", linewidth=1.5),
        flierprops=dict(marker="o", markerfacecolor="gray", markersize=3, markeredgecolor="black", alpha=0.7),
    )
    for i, data in enumerate(pdata):
        x = np.random.normal(i + 1, 0.04, len(data))
        ax.scatter(x, data, color="gray", alpha=0.8, s=20, zorder=3)
    ax.set_ylabel("Dendritic excitability index (distal/proximal)")
    ax.set_ylim(0, max(2.0, np.nanmax([d.max() if d.size else 0 for d in pdata]) * 1.2))
    ax.set_title("Figure 1I: Dendritic Excitability (streamed)")
    fig_path = out / "figure_1i_dendritic_excitability.png"
    fig.savefig(fig_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print("Saved:", fig_path)


if __name__ == "__main__":
    main()
