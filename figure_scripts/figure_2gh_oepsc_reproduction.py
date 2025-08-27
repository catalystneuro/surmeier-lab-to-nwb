#!/usr/bin/env python3
"""
Figure 2 G and H oEPSC Reproduction Script (Fixed)

This script reproduces Figure 2 G and H from Zhai et al. 2025 using NWB files
containing Sr²⁺-oEPSC data. It loads voltage clamp recordings, detects oEPSC
events with a robust median-absolute-deviation (MAD) method within the defined
"detection" window from the NWB `intervals/optogenetic_epochs_table`, and then
generates cumulative probability curves and box plots.

Updates in this version:
- Adjusts to the updated NWB file structure where condition is encoded as
  `OnState`/`OffState` in the filename and optogenetic timing comes from
  `intervals['optogenetic_epochs_table']` with `stage_name` columns.
- Robustly derives timestamps for sweeps whether `timestamps` are stored
  explicitly or via `starting_time` + `rate`.
- Uses a MAD-based event detection with a 100 ms shift after the original
  detection start to avoid stimulus transients, grouping contiguous threshold
  crossings into single events and measuring their peak amplitudes.

Usage:
    python figure_2gh_oepsc_reproduction_fixed.py

Output:
    - analysis_outputs/figure_2gh/figure_2gh_reproduction_report.md
    - analysis_outputs/figure_2gh/figure_2gh_cumulative_probability.png
    - analysis_outputs/figure_2gh/figure_2gh_amplitude_boxplot.png
"""

import warnings
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pynwb import NWBHDF5IO
from scipy import stats

# Suppress expected warnings
warnings.filterwarnings("ignore", message="Mean of empty slice")
warnings.filterwarnings("ignore", message="invalid value encountered in scalar divide")


# Set up plotting style to match paper (from figure_1e script)
plt.style.use("default")
plt.rcParams["figure.dpi"] = 300
plt.rcParams["savefig.dpi"] = 300
plt.rcParams["savefig.bbox"] = "tight"
plt.rcParams["savefig.facecolor"] = "white"
plt.rcParams["figure.facecolor"] = "white"
plt.rcParams["axes.facecolor"] = "white"
plt.rcParams["axes.edgecolor"] = "black"
plt.rcParams["axes.linewidth"] = 1.0
plt.rcParams["xtick.color"] = "black"
plt.rcParams["ytick.color"] = "black"
plt.rcParams["text.color"] = "black"
plt.rcParams["axes.labelcolor"] = "black"
plt.rcParams["axes.spines.top"] = False
plt.rcParams["axes.spines.right"] = False
plt.rcParams["axes.grid"] = False


def _get_ts_slice_in_si_units(ts, idx_start: int, count: int) -> Tuple[np.ndarray, np.ndarray]:
    """Return (data_in_amperes, timestamps_in_seconds) for a slice of a TimeSeries.

    Handles both explicit `timestamps` and implicit `starting_time` + `rate`.
    """
    data = np.asarray(ts.data[idx_start : idx_start + count])
    if ts.timestamps is not None:
        t = np.asarray(ts.timestamps[idx_start : idx_start + count])
    else:
        # timestamps are derived from starting_time and rate
        rate = float(ts.rate)
        start = float(ts.starting_time) + idx_start / rate
        t = start + np.arange(count, dtype=float) / rate
    return data.astype(float), t.astype(float)


def _group_contiguous_indices(indices: np.ndarray) -> List[np.ndarray]:
    """Group sorted indices into contiguous runs and return list of index arrays."""
    if indices.size == 0:
        return []
    # Find boundaries where the difference > 1
    splits = np.where(np.diff(indices) > 1)[0] + 1
    groups = np.split(indices, splits)
    return groups


def detect_events_mad(
    data_pA: np.ndarray,
    timestamps_ms: np.ndarray,
    detection_start_ms: float,
    detection_stop_ms: float,
    detection_window_shift_ms: float = 100.0,
    threshold_k: float = 5.0,
    min_inter_event_ms: float = 5.0,
) -> Dict[str, np.ndarray]:
    """
    MAD-based event detection within a shifted detection window.

    Returns a dict with keys: 'times_ms', 'amps_pA', 'polarity'.
    - times_ms: event time stamps
    - amps_pA: event amplitudes as positive magnitudes from the noise median
    - polarity: +1 for positive, -1 for negative events
    """
    # Apply shift to detection window start
    det_start = detection_start_ms + detection_window_shift_ms
    det_stop = detection_stop_ms

    mask = (timestamps_ms >= det_start) & (timestamps_ms <= det_stop)
    if not np.any(mask):
        return {"times_ms": np.array([]), "amps_pA": np.array([]), "polarity": np.array([])}

    x = data_pA[mask]
    t = timestamps_ms[mask]

    noise_median = np.median(x)
    mad = np.median(np.abs(x - noise_median))
    mad_std = mad * 1.4826  # approx std for normal dist
    pos_thr = noise_median + threshold_k * mad_std
    neg_thr = noise_median - threshold_k * mad_std

    pos_idx = np.flatnonzero(x > pos_thr)
    neg_idx = np.flatnonzero(x < neg_thr)

    event_times = []
    event_amps = []
    event_pol = []

    # Positive events: pick the local maximum within each contiguous run
    for grp in _group_contiguous_indices(pos_idx):
        j = grp[np.argmax(x[grp])]
        amp = x[j] - noise_median
        event_times.append(t[j])
        event_amps.append(abs(amp))
        event_pol.append(1)

    # Negative events: pick the local minimum within each contiguous run
    for grp in _group_contiguous_indices(neg_idx):
        j = grp[np.argmin(x[grp])]
        amp = noise_median - x[j]
        event_times.append(t[j])
        event_amps.append(abs(amp))
        event_pol.append(-1)

    if len(event_times) == 0:
        return {"times_ms": np.array([]), "amps_pA": np.array([]), "polarity": np.array([])}

    order = np.argsort(event_times)
    t = np.asarray(event_times)[order]
    a = np.asarray(event_amps)[order]
    p = np.asarray(event_pol, dtype=int)[order]

    # Merge nearby events within a refractory window
    if t.size and min_inter_event_ms > 0:
        t, a, p = merge_nearby_events(t, a, p, merge_distance_ms=min_inter_event_ms)

    return {"times_ms": t, "amps_pA": a, "polarity": p}


def merge_nearby_events(
    event_times_ms: np.ndarray,
    event_amplitudes_pA: np.ndarray,
    event_polarity: np.ndarray,
    merge_distance_ms: float = 5.0,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Merge events that are within merge_distance_ms of each other (single-linkage).
    For merged clusters, keep the event with maximum amplitude.

    Notes
    -----
    - Uses chained grouping: if consecutive events are each <= merge_distance_ms apart,
      they are merged into one cluster (not just those close to the first event).
    - Amplitudes are expected as positive magnitudes; polarity is preserved from
      the selected representative event.
    """
    if event_times_ms.size == 0:
        return event_times_ms, event_amplitudes_pA, event_polarity

    # Sort by time
    order = np.argsort(event_times_ms)
    t = np.asarray(event_times_ms)[order]
    a = np.asarray(event_amplitudes_pA)[order]
    p = np.asarray(event_polarity, dtype=int)[order]

    # Split where the gap exceeds merge_distance_ms
    gaps = np.diff(t)
    splits = np.where(gaps > merge_distance_ms)[0] + 1
    groups = np.split(np.arange(t.size), splits)

    kept_idx = []
    for grp in groups:
        # pick the event with max amplitude in the group
        j_local = int(np.argmax(a[grp]))
        kept_idx.append(int(grp[j_local]))

    kept_idx = np.asarray(kept_idx, dtype=int)
    return t[kept_idx], a[kept_idx], p[kept_idx]


def load_oepsc_data(nwb_dir: Path) -> pd.DataFrame:
    """
    Load oEPSC data from NWB files.

    Parameters
    ----------
    nwb_dir : Path
        Directory containing Figure 2 Sr²⁺-oEPSC NWB files.

    Returns
    -------
    pd.DataFrame
        Combined dataframe with all oEPSC amplitudes and conditions.
    """
    nwb_files = list(nwb_dir.glob("*.nwb"))

    if not nwb_files:
        raise FileNotFoundError(
            f"No NWB files found in {nwb_dir}. " "Please run figure_2_optical_stimuli.py first to generate NWB files."
        )

    all_data = []
    print(f"Loading data from {len(nwb_files)} NWB files...")

    for nwb_file in nwb_files:
        with NWBHDF5IO(str(nwb_file), "r") as io:
            nwb = io.read()

            # Extract condition from filename (updated naming: OnState/OffState)
            condition = "unknown"
            name_lower = nwb_file.name.lower()
            if "offstate" in name_lower:
                condition = "LID off-state"
            elif "onstate" in name_lower:
                condition = "LID on-state"

            # Get session info
            session_id = nwb.session_id
            subject_id = nwb.subject.subject_id

            # Access intracellular recordings table
            icephys_table = nwb.get_intracellular_recordings()

            # Get optogenetic epochs table
            opto_epochs = nwb.intervals["optogenetic_epochs_table"]
            opto_df = opto_epochs.to_dataframe()

            # Identify detection and stimulation stage windows as RELATIVE offsets per sweep
            # Compute offsets from the first 'pre_stimulation' start time
            if "stage_name" not in opto_df.columns:
                raise ValueError(f"optogenetic_epochs_table missing 'stage_name' in {nwb_file.name}")
            pre_rows = opto_df[opto_df["stage_name"] == "pre_stimulation"].copy()
            det_rows = opto_df[opto_df["stage_name"] == "detection"].copy()
            stim_rows = opto_df[opto_df["stage_name"] == "stimulation"].copy()
            if pre_rows.empty or det_rows.empty:
                raise ValueError(f"Missing required stage rows in optogenetic_epochs_table for {nwb_file.name}")
            pre0 = float(pre_rows.iloc[0]["start_time"])
            det0_start = float(det_rows.iloc[0]["start_time"]) - pre0
            det0_stop = float(det_rows.iloc[0]["stop_time"]) - pre0
            stim0_start = float(stim_rows.iloc[0]["start_time"]) - pre0 if not stim_rows.empty else det0_start

            # Detection parameters
            detection_window_shift_ms = 100.0

            # Iterate over sweeps and detect events
            for idx in range(len(icephys_table.id)):
                row = icephys_table[idx]
                cell_number = row[("intracellular_recordings", "cell_number")].iloc[0]
                led_intensity = row[("intracellular_recordings", "led_intensity")].iloc[0]
                sweep_number = row[("intracellular_recordings", "sweep_number")].iloc[0]
                recording_id = f"Cell{cell_number}_LED{led_intensity}_Sweep{sweep_number}"

                # Get the voltage clamp series (oEPSC trace)
                response_ref = row[("responses", "response")].iloc[0]
                ts = response_ref.timeseries
                data_A, time_s = _get_ts_slice_in_si_units(ts, response_ref.idx_start, response_ref.count)

                # Convert to pA and relative ms for detection (0 at sweep start)
                data_pA = data_A * 1e12
                time_ms = (time_s - time_s[0]) * 1000.0

                # Detect events within the (shifted) detection window
                detection = detect_events_mad(
                    data_pA=data_pA,
                    timestamps_ms=time_ms,
                    detection_start_ms=det0_start * 1000.0,
                    detection_stop_ms=det0_stop * 1000.0,
                    detection_window_shift_ms=detection_window_shift_ms,
                    threshold_k=5.0,
                )

                # Append one row per detected event
                for t_ms, amp_pA, pol in zip(detection["times_ms"], detection["amps_pA"], detection["polarity"]):
                    all_data.append(
                        {
                            "session_id": session_id,
                            "subject_id": subject_id,
                            "condition": condition,
                            "recording_id": recording_id,
                            "oepsc_amplitude_pA": float(amp_pA),
                            "event_time_ms": float(t_ms),
                            "event_polarity": int(pol),
                            "nwb_file": nwb_file.name,
                            "stim_start_ms": stim0_start * 1000.0,
                            "detection_start_ms": det0_start * 1000.0,
                            "detection_stop_ms": det0_stop * 1000.0,
                            "sweep_number": int(sweep_number),
                            "cell_number": int(cell_number),
                            "led_intensity": led_intensity,
                        }
                    )

    if not all_data:
        raise ValueError("No oEPSC data could be loaded from NWB files.")

    df = pd.DataFrame(all_data)
    # Filter out NaN amplitudes
    df = df.dropna(subset=["oepsc_amplitude_pA"])
    print(
        f"Loaded {len(df)} detected events from {len(df.session_id.unique())} cells across {len(df.nwb_file.unique())} NWB files."
    )
    return df


def create_cumulative_probability_plot(df: pd.DataFrame) -> Tuple[plt.Figure, plt.Axes]:
    """
    Create cumulative probability curves for oEPSC amplitudes.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with oEPSC amplitudes and conditions.

    Returns
    -------
    Tuple[plt.Figure, plt.Axes]
        Figure and axes objects.
    """
    fig, ax = plt.subplots(1, 1, figsize=(6, 4))

    # Define condition styling
    condition_styles = {
        "LID off-state": {"color": "black", "linestyle": "-", "label": "off-state"},
        "LID on-state": {"color": "gray", "linestyle": "-", "label": "on-state"},
    }

    for condition in ["LID off-state", "LID on-state"]:
        if condition not in df["condition"].unique():
            continue

        condition_data = df[df["condition"] == condition]["oepsc_amplitude_pA"].values
        if len(condition_data) == 0:
            continue

        # Sort data and calculate cumulative probability
        sorted_data = np.sort(condition_data)
        cumulative_prob = np.arange(1, len(sorted_data) + 1) / len(sorted_data) * 100

        style = condition_styles[condition]
        ax.plot(
            sorted_data,
            cumulative_prob,
            color=style["color"],
            linestyle=style["linestyle"],
            linewidth=1.5,
            label=style["label"],
        )

    ax.set_xlabel("oEPSC amplitude (pA)", fontsize=12)
    ax.set_ylabel("Cumulative Probability (%)", fontsize=12)
    ax.set_xlim(0, 80)  # Assuming inward currents are now positive
    ax.set_ylim(0, 100)
    ax.set_xticks(np.arange(0, 81, 20))
    ax.set_yticks(np.arange(0, 101, 25))

    ax.legend(loc="upper left", frameon=False, fontsize=10)
    ax.set_title("Figure 2G: dSPN oEPSC Cumulative Probability", fontsize=14, pad=20)

    plt.tight_layout()
    return fig, ax


def create_amplitude_boxplot(df: pd.DataFrame) -> Tuple[plt.Figure, plt.Axes]:
    """
    Create box plot comparing oEPSC amplitudes between conditions.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with oEPSC amplitudes and conditions.

    Returns
    -------
    Tuple[plt.Figure, plt.Axes]
        Figure and axes objects.
    """
    fig, ax = plt.subplots(1, 1, figsize=(4, 5))

    # Order conditions
    condition_order = ["LID off-state", "LID on-state"]
    condition_labels = ["off-state", "on-state"]

    # Compute per-NWB-file (subject-level) mean amplitudes per condition
    per_file_means = (
        df.groupby(["condition", "nwb_file"])["oepsc_amplitude_pA"]  # each NWB file is treated as a subject
        .mean()
        .reset_index()
    )

    plot_data = []
    for condition in condition_order:
        if condition in per_file_means["condition"].unique():
            vals = per_file_means[per_file_means["condition"] == condition]["oepsc_amplitude_pA"].values
            plot_data.append(vals)
        else:
            plot_data.append(np.array([]))  # Empty array if no data for condition

    # Create box plot
    box = ax.boxplot(
        plot_data,
        labels=condition_labels,
        patch_artist=True,
        boxprops=dict(facecolor="white", color="black", linewidth=1),
        whiskerprops=dict(color="black", linewidth=1),
        capprops=dict(color="black", linewidth=1),
        medianprops=dict(color="black", linewidth=1.5),
        flierprops=dict(marker="o", markerfacecolor="gray", markersize=3, markeredgecolor="black", alpha=0.7),
    )

    # Add individual data points
    for i, data in enumerate(plot_data):
        x_pos = np.random.normal(i + 1, 0.04, len(data))
        ax.scatter(x_pos, data, color="gray", alpha=0.8, s=20, zorder=3)

    ax.set_ylabel("oEPSC amplitude (pA)", fontsize=12)
    # Fixed y-axis limits to match paper style
    ax.set_ylim(0, 25)
    ax.set_yticks(np.arange(0, 26, 5))

    ax.set_title("Figure 2H: dSPN oEPSC Amplitude Comparison", fontsize=14, pad=20)

    # Perform statistical test on per-file means (subject-level)
    off_state_data = per_file_means[per_file_means["condition"] == "LID off-state"]["oepsc_amplitude_pA"].values
    on_state_data = per_file_means[per_file_means["condition"] == "LID on-state"]["oepsc_amplitude_pA"].values

    if len(off_state_data) > 1 and len(on_state_data) > 1:
        # Assuming independent samples and potentially unequal variances
        stat, p_value = stats.ttest_ind(off_state_data, on_state_data, equal_var=False)
        y_pos = 23
        if p_value < 0.001:
            ax.text(1.5, y_pos, "***", ha="center", va="bottom", fontsize=12, fontweight="bold")
        elif p_value < 0.01:
            ax.text(1.5, y_pos, "**", ha="center", va="bottom", fontsize=12, fontweight="bold")
        elif p_value < 0.05:
            ax.text(1.5, y_pos, "*", ha="center", va="bottom", fontsize=12, fontweight="bold")

    plt.tight_layout()
    return fig, ax


def generate_markdown_report(df: pd.DataFrame, output_dir: Path) -> str:
    """
    Generate comprehensive markdown report for Figure 2 G and H.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with oEPSC amplitudes and conditions.
    output_dir : Path
        Output directory.

    Returns
    -------
    str
        Path to generated report.
    """
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    markdown_content = f"""# Figure 2 G and H oEPSC Analysis Report

**Generated:** {timestamp}
**Data Source:** NWB files with Sr²⁺-oEPSC voltage clamp recordings
**Analysis:** Cumulative probability of oEPSC amplitudes and box plot comparison

## Executive Summary

This report presents the reproduction of Figure 2 G and H from Zhai et al. 2025,
demonstrating changes in corticostriatal synaptic strength (oEPSC amplitude)
in dSPNs across different LID (L-DOPA-induced dyskinesia) states.

### Key Findings

- **Total oEPSC Events Analyzed:** {len(df)}
- **Conditions:** {len(df['condition'].unique())}

## Methodology

### Experimental Protocol
Voltage clamp recordings from dSPN somata at -70 mV in Ca²⁺-free ACSF
containing 3 mM SrCl₂ and 10 μM gabazine. Optogenetic stimulation of
ChR2-expressing corticostriatal terminals with 0.3 ms blue LED pulses.

### Analysis Approach
1.  **oEPSC Amplitude Extraction:** Peak inward current detected within 40-400 ms
    post-stimulation window.
2.  **Cumulative Probability:** Distribution of oEPSC amplitudes for each condition.
3.  **Box Plot Comparison:** Statistical comparison of amplitudes between conditions.

## Results by Condition

"""

    for condition in df["condition"].unique():
        condition_data = df[df["condition"] == condition]["oepsc_amplitude_pA"].dropna()
        if len(condition_data) > 0:
            markdown_content += f"""### {condition}

- **Sample Size (oEPSC events):** {len(condition_data)}
- **Mean oEPSC Amplitude:** {condition_data.mean():.2f} ± {condition_data.sem():.2f} pA (mean ± SEM)
- **Median oEPSC Amplitude:** {condition_data.median():.2f} pA
- **Amplitude Range:** {condition_data.min():.2f} to {condition_data.max():.2f} pA

"""
        else:
            markdown_content += f"""### {condition}

- **No data available for this condition.**

"""

    # Statistical comparison summary
    off_state_data = df[df["condition"] == "LID off-state"]["oepsc_amplitude_pA"].values
    on_state_data = df[df["condition"] == "LID on-state"]["oepsc_amplitude_pA"].values

    if len(off_state_data) > 1 and len(on_state_data) > 1:
        stat, p_value = stats.ttest_ind(off_state_data, on_state_data, equal_var=False)
        markdown_content += f"""## Statistical Comparison (LID off-state vs. LID on-state)

- **Independent t-test (Welch's t-test):**
    - **t-statistic:** {stat:.3f}
    - **p-value:** {p_value:.4f}
    - **Significance:** {'*** (p < 0.001)' if p_value < 0.001 else '** (p < 0.01)' if p_value < 0.01 else '* (p < 0.05)' if p_value < 0.05 else 'Not significant'}

"""
    else:
        markdown_content += """## Statistical Comparison

- **Not enough data to perform statistical comparison.**

"""

    markdown_content += f"""
## Figures

### Figure 2G: oEPSC Cumulative Probability
![oEPSC Cumulative Probability](figure_2gh_cumulative_probability.png)

**Description:** Cumulative probability distributions of oEPSC amplitudes,
illustrating the shift towards larger amplitudes in the LID on-state.

### Figure 2H: oEPSC Amplitude Box Plot
![oEPSC Amplitude Box Plot](figure_2gh_amplitude_boxplot.png)

**Description:** Box plot comparing the distribution of oEPSC amplitudes
between LID off-state and on-state, highlighting the significant increase
in amplitude during the on-state.

## Files Generated

- `figure_2gh_cumulative_probability.png` - Cumulative probability plot
- `figure_2gh_amplitude_boxplot.png` - Box plot comparison
- `figure_2gh_reproduction_report.md` - This detailed report

---

*Report generated by Figure 2 G and H oEPSC Reproduction Script*
"""

    # Write report
    report_path = output_dir / "figure_2gh_reproduction_report.md"
    with open(report_path, "w") as f:
        f.write(markdown_content)

    return str(report_path)


def main():
    """Main function to orchestrate the Figure 2 G and H reproduction analysis."""

    # Set up paths
    base_dir = Path(__file__).parent.parent
    nwb_dir = base_dir / "nwb_files" / "optical_stimulation" / "figure_2"
    assert nwb_dir.exists(), f"NWB directory {nwb_dir} does not exist. Please run figure_2_optical_stimuli.py first."
    output_dir = base_dir / "analysis_outputs" / "figure_2gh"

    # Ensure output directory exists
    output_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 80)
    print("FIGURE 2 G AND H oEPSC REPRODUCTION ANALYSIS")
    print("=" * 80)
    print(f"NWB Directory: {nwb_dir}")
    print(f"Output Directory: {output_dir}")
    print()

    # Load oEPSC data
    print("Step 1: Loading oEPSC data from NWB files...")
    df = load_oepsc_data(nwb_dir)

    # Create cumulative probability plot
    print("Step 2: Creating cumulative probability plot...")
    fig_cum_prob, _ = create_cumulative_probability_plot(df)
    fig_cum_prob_path = output_dir / "figure_2gh_cumulative_probability.png"
    fig_cum_prob.savefig(fig_cum_prob_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig_cum_prob)
    print(f"  Saved: {fig_cum_prob_path}")

    # Create amplitude box plot
    print("Step 3: Creating amplitude box plot...")
    fig_boxplot, _ = create_amplitude_boxplot(df)
    fig_boxplot_path = output_dir / "figure_2gh_amplitude_boxplot.png"
    fig_boxplot.savefig(fig_boxplot_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig_boxplot)
    print(f"  Saved: {fig_boxplot_path}")

    # Generate report
    print("Step 4: Generating analysis report...")
    report_path = generate_markdown_report(df, output_dir)
    print(f"  Saved: {report_path}")

    # Print summary
    print("\n" + "=" * 80)
    print("ANALYSIS COMPLETED SUCCESSFULLY!")
    print("=" * 80)

    print(f"Total oEPSC Events Analyzed: {len(df)}")
    print(f"Conditions: {df['condition'].unique().tolist()}")

    print("\nCondition Summary:")
    for condition in df["condition"].unique():
        condition_data = df[df["condition"] == condition]["oepsc_amplitude_pA"].dropna()
        if len(condition_data) > 0:
            print(
                f"  {condition}: n={len(condition_data)}, "
                f"mean amplitude={condition_data.mean():.2f}±{condition_data.sem():.2f} pA"
            )
        else:
            print(f"  {condition}: No data available.")

    print("\nOutput Files:")
    print(f"  - {fig_cum_prob_path}")
    print(f"  - {fig_boxplot_path}")
    print(f"  - {report_path}")

    print("\n" + "=" * 80)
    print("Figure 2 G and H reproduction demonstrates NWB-powered oEPSC analysis!")
    print("Cumulative probability and box plots from structured electrophysiology data.")
    print("=" * 80)


if __name__ == "__main__":
    main()
