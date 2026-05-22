"""Reproduce Figure 7C (CDGI KO somatic excitability) from merged NWB files.

The paper's Figure 7C plots current-response curves (number of action potentials
per 500-ms current step, across step amplitudes from ~-120 to ~300 pA) for two
conditions: off-state (n=13 cells / 5 mice) and on-state (n=8 cells / 4 mice).

This script:
  1. Walks the merged F7 NWB files
  2. For each somatic sweep, extracts the voltage trace via IRT idx_start + count
  3. Detects action potentials in the step window (200-700 ms) via dV/dt threshold (20 mV/ms)
  4. Groups cells by condition_subfolder
  5. Plots the mean+SEM F-I curve per condition
  6. Saves the plot for visual comparison with the paper's Fig 7C

If the resulting plot shows higher firing rates for off-state cells than on-state cells
(the paper's main finding), the merger is correctly preserving the experimental structure.
"""

from __future__ import annotations

from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from pynwb import NWBHDF5IO

REPO_ROOT = Path(__file__).resolve().parents[1]
MERGED_F7 = REPO_ROOT / "nwb_files" / "merged" / "F7"
OUTPUT_PNG = REPO_ROOT / "nwb_files" / "merged" / "validation_fig7c.png"

# Step window from the paper (within-sweep): "[200, 700]ms" -> 200 ms to 700 ms.
STEP_START_MS = 200.0
STEP_END_MS = 700.0
# AP threshold: dV/dt > 20 mV/ms (field-standard, mentioned in the paper).
DVDT_THRESHOLD_MVPMS = 20.0


def count_aps_in_window(voltage_v: np.ndarray, rate_hz: float) -> int:
    """Count action potentials in voltage trace via dV/dt threshold crossing.

    voltage_v: 1-D array of membrane voltage in volts.
    rate_hz: sample rate in Hz.
    """
    voltage_mv = voltage_v * 1000.0  # V -> mV
    sample_period_ms = 1000.0 / rate_hz
    # First derivative in mV/ms.
    dvdt = np.diff(voltage_mv) / sample_period_ms
    # Find upward threshold crossings.
    above = dvdt > DVDT_THRESHOLD_MVPMS
    crossings = np.where((~above[:-1]) & above[1:])[0]
    # Enforce a refractory of ~2 ms to avoid double-counting one AP.
    refractory_samples = int(2.0 / sample_period_ms)
    last_crossing = -refractory_samples - 1
    n_aps = 0
    for c in crossings:
        if c - last_crossing >= refractory_samples:
            n_aps += 1
            last_crossing = c
    return n_aps


def extract_fi_curves() -> dict[str, dict[float, list[int]]]:
    """Walk F7 merged files. Returns {condition: {step_pA: [ap_counts_per_cell]}}."""
    # Cell-level AP count per current step. Average sweeps per cell first.
    # Structure: {condition: {step_pA: {cell_key: [ap_counts]}}}
    raw: dict[str, dict[float, dict[tuple[str, int], list[int]]]] = defaultdict(
        lambda: defaultdict(lambda: defaultdict(list))
    )
    for nwb_path in sorted(MERGED_F7.glob("*.nwb")):
        with NWBHDF5IO(nwb_path, "r") as io:
            nwb = io.read()
            df = nwb.intracellular_recordings.to_dataframe()
            for _, row in df.iterrows():
                dtype = row[("intracellular_recordings", "dendrite_type")]
                if dtype != "Soma":
                    continue
                cond = row[("intracellular_recordings", "condition_subfolder")]
                cell_id = int(row[("intracellular_recordings", "cell_id")])
                stim_pA = float(row[("intracellular_recordings", "stimulus_current_pA")])
                # Pull the voltage trace for this sweep via the response reference.
                resp_ref = row[("responses", "response")]
                ts = resp_ref.timeseries
                idx_start = int(resp_ref.idx_start)
                count = int(resp_ref.count)
                voltage = ts.data[idx_start : idx_start + count]
                rate = float(ts.rate)
                # Within-sweep step window (200-700 ms) -> sample indices.
                step_start_idx = int(STEP_START_MS * 1e-3 * rate)
                step_end_idx = int(STEP_END_MS * 1e-3 * rate)
                if step_end_idx > len(voltage):
                    step_end_idx = len(voltage)
                if step_start_idx >= step_end_idx:
                    continue
                window = voltage[step_start_idx:step_end_idx]
                n_aps = count_aps_in_window(np.asarray(window), rate)
                cell_key = (nwb.subject.subject_id, cell_id)
                raw[cond][stim_pA][cell_key].append(n_aps)
    # Collapse multi-sweep-per-cell to per-cell mean.
    per_cell: dict[str, dict[float, list[int]]] = defaultdict(lambda: defaultdict(list))
    for cond, steps in raw.items():
        for step_pA, by_cell in steps.items():
            for cell_key, sweep_counts in by_cell.items():
                per_cell[cond][step_pA].append(int(np.mean(sweep_counts)))
    return per_cell


def plot_fi_curves(per_cell: dict[str, dict[float, list[int]]], out_path: Path) -> None:
    """Plot mean+SEM F-I curves per condition."""
    fig, ax = plt.subplots(figsize=(6, 5))
    colors = {"KO off-state": "tab:blue", "KO on-state": "tab:red"}
    for condition, steps in per_cell.items():
        sorted_steps = sorted(steps.keys())
        means = []
        sems = []
        n_cells_per_step = []
        for s in sorted_steps:
            counts = np.asarray(steps[s], dtype=float)
            means.append(np.mean(counts))
            sems.append(np.std(counts, ddof=1) / np.sqrt(len(counts)) if len(counts) > 1 else 0.0)
            n_cells_per_step.append(len(counts))
        means = np.asarray(means)
        sems = np.asarray(sems)
        steps_arr = np.asarray(sorted_steps)
        max_n = max(n_cells_per_step)
        label = f"{condition} (n_max={max_n} cells)"
        c = colors.get(condition, "black")
        ax.errorbar(steps_arr, means, yerr=sems, fmt="o-", color=c, label=label, capsize=3)
    ax.set_xlabel("Current step (pA)")
    ax.set_ylabel("Number of APs (within 500-ms step)")
    ax.set_title("Figure 7C reproduction — CDGI KO somatic excitability\n(from merged per-mouse-day NWB files)")
    ax.legend(loc="upper left")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_path, dpi=120)
    plt.close(fig)
    print(f"Saved: {out_path}")


def main():
    per_cell = extract_fi_curves()
    print(f"Conditions found: {sorted(per_cell.keys())}")
    for cond, steps in per_cell.items():
        n_cells_total = max((len(v) for v in steps.values()), default=0)
        n_steps = len(steps)
        print(f"  {cond}: {n_cells_total} cells (max across steps), {n_steps} distinct current steps")
    plot_fi_curves(per_cell, OUTPUT_PNG)
    # Print numerical summary for spot-check.
    print()
    print("Mean APs per current step (paper Fig 7C visual check):")
    print(f"{'Step (pA)':>10s}  {'off-state':>18s}  {'on-state':>18s}")
    all_steps = sorted({s for steps in per_cell.values() for s in steps.keys()})
    for step in all_steps:
        off_counts = per_cell.get("KO off-state", {}).get(step, [])
        on_counts = per_cell.get("KO on-state", {}).get(step, [])
        off_str = f"{np.mean(off_counts):.1f} (n={len(off_counts)})" if off_counts else "-"
        on_str = f"{np.mean(on_counts):.1f} (n={len(on_counts)})" if on_counts else "-"
        print(f"{step:>10.0f}  {off_str:>18s}  {on_str:>18s}")


if __name__ == "__main__":
    main()
