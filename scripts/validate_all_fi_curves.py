"""Reproduce F-I curves from all somatic-excitability panels.

For each figure with somatic excitability data, extract per-cell AP counts per
current step and plot mean+SEM F-I curves. Compare qualitatively against the
paper's published findings.
"""

from __future__ import annotations

from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from pynwb import NWBHDF5IO

REPO_ROOT = Path(__file__).resolve().parents[1]
MERGED_ROOT = REPO_ROOT / "nwb_files" / "merged"
OUTPUT_PNG = MERGED_ROOT / "validation_all_fi_curves.png"

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


def extract_panel(figure: str, condition_filter) -> dict[float, dict[tuple[str, int], list[int]]]:
    """For one figure, extract AP counts per (step, cell) for somatic sweeps matching condition_filter.

    condition_filter is a callable: condition_subfolder_string -> bool.
    """
    per_step_per_cell: dict[float, dict[tuple[str, int], list[int]]] = defaultdict(lambda: defaultdict(list))
    fig_dir = MERGED_ROOT / figure
    if not fig_dir.exists():
        return per_step_per_cell
    for nwb_path in sorted(fig_dir.glob("*.nwb")):
        with NWBHDF5IO(nwb_path, "r") as io:
            nwb = io.read()
            df = nwb.intracellular_recordings.to_dataframe()
            for _, row in df.iterrows():
                if row[("intracellular_recordings", "dendrite_type")] != "Soma":
                    continue
                cond = row[("intracellular_recordings", "condition_subfolder")]
                if not condition_filter(cond):
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
                per_step_per_cell[stim_pA][cell_key].append(n_aps)
    return per_step_per_cell


def collapse_to_per_cell(per_step_per_cell: dict[float, dict[tuple[str, int], list[int]]]) -> dict[float, list[int]]:
    per_cell: dict[float, list[int]] = defaultdict(list)
    for step, by_cell in per_step_per_cell.items():
        for cell_key, counts in by_cell.items():
            per_cell[step].append(int(np.mean(counts)))
    return per_cell


# Define panels: (figure, panel_label, [(condition_label, color, predicate), ...])
PANELS = [
    (
        "F1",
        "Fig 1F dSPN somatic",
        [
            ("LID off-state", "black", lambda c: c == "LID off-state"),
            ("LID on-state", "red", lambda c: c == "LID on-state"),
            ("LID on + SCH", "blue", lambda c: c == "LID on-state with SCH"),
        ],
    ),
    (
        "F3",
        "Fig 3D iSPN somatic",
        [
            ("LID off-state", "black", lambda c: c == "LID off-state"),
            ("LID on-state", "red", lambda c: c == "LID on-state"),
            ("LID on + sul", "blue", lambda c: "with sul" in c),
        ],
    ),
    (
        "F6",
        "Fig 6C iSPN somatic (M1R antagonist)",
        [("control", "black", lambda c: c == "control"), ("M1R antagonist", "red", lambda c: c == "M1R antagonist")],
    ),
    (
        "F7",
        "Fig 7C CDGI-KO iSPN somatic",
        [
            ("KO off-state", "black", lambda c: c == "KO off-state"),
            ("KO on-state", "red", lambda c: c == "KO on-state"),
        ],
    ),
    (
        "F8",
        "Fig 8D M1R CRISPR iSPN somatic",
        [
            ("interleaved control", "black", lambda c: c == "interleaved control"),
            ("M1R CRISPR", "red", lambda c: c == "M1R CRISPR"),
        ],
    ),
]


def plot_panel(ax, fig: str, title: str, conditions: list):
    for cond_label, color, predicate in conditions:
        raw = extract_panel(fig, predicate)
        per_cell = collapse_to_per_cell(raw)
        if not per_cell:
            continue
        steps = sorted(per_cell.keys())
        means = []
        sems = []
        ns = []
        for s in steps:
            counts = np.asarray(per_cell[s], dtype=float)
            means.append(np.mean(counts))
            sems.append(np.std(counts, ddof=1) / np.sqrt(len(counts)) if len(counts) > 1 else 0.0)
            ns.append(len(counts))
        max_n = max(ns) if ns else 0
        ax.errorbar(
            steps,
            means,
            yerr=sems,
            fmt="o-",
            color=color,
            label=f"{cond_label} (n={max_n})",
            capsize=2,
            markersize=4,
            linewidth=1.2,
            alpha=0.85,
        )
    ax.set_title(title, fontsize=10)
    ax.set_xlabel("Current step (pA)", fontsize=9)
    ax.set_ylabel("# APs / 500 ms", fontsize=9)
    ax.legend(fontsize=8, loc="upper left")
    ax.grid(True, alpha=0.3)


def main():
    fig, axes = plt.subplots(2, 3, figsize=(15, 9))
    flat = axes.flatten()
    for i, (fig_code, title, conditions) in enumerate(PANELS):
        plot_panel(flat[i], fig_code, title, conditions)
    flat[-1].axis("off")
    fig.suptitle(
        "F-I curves from merged per-mouse-day NWB files (qualitative paper reproduction)", fontsize=12, y=0.998
    )
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(OUTPUT_PNG, dpi=110)
    plt.close(fig)
    print(f"Saved: {OUTPUT_PNG}")

    # Numerical summary at 300 pA per condition.
    print()
    print(f"{'Panel':50s}  {'condition':25s}  mean APs @300 pA (n)")
    print("-" * 110)
    for fig_code, title, conditions in PANELS:
        for cond_label, _, predicate in conditions:
            raw = extract_panel(fig_code, predicate)
            per_cell = collapse_to_per_cell(raw)
            counts_300 = per_cell.get(300.0, [])
            if counts_300:
                m = np.mean(counts_300)
                print(f"  {title:48s}  {cond_label:25s}  {m:6.1f}  (n={len(counts_300)})")


if __name__ == "__main__":
    main()
