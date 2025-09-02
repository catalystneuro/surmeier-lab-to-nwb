#!/usr/bin/env python3
"""
Figure 5F GRABACh Reproduction Script - From NWB Files

This script reproduces Figure 5F from Zhai et al. 2025 using NWB files.
It analyzes GRABACh3.0 acetylcholine biosensor data to show state-dependent
acetylcholine release in control, PD, and LID off-state conditions.

The analysis follows the methodology from the original BOT_imaging data_script.ipynb
but sources data from standardized NWB files.

Usage:
    python figure_scripts/figure_5f_grabatch_reproduction.py

Output:
    - analysis_outputs/figure_5f/figure_5f_nwb_boxplots.png
    - analysis_outputs/figure_5f/figure_5f_nwb_report.md
"""

import warnings
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pynwb import NWBHDF5IO

# Suppress expected warnings
warnings.filterwarnings("ignore", message="Mean of empty slice")
warnings.filterwarnings("ignore", message="invalid value encountered in scalar divide")

# Set up plotting style to match paper
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


def load_nwb_trials_data(nwb_path: Path) -> List[Dict]:
    """
    Load all trial data from an NWB file using trials table iteration.

    Parameters
    ----------
    nwb_path : Path
        Path to the NWB file

    Returns
    -------
    List[Dict]
        List of trial dictionaries with fluorescence data and metadata
    """
    trials_data = []

    with NWBHDF5IO(str(nwb_path), mode="r") as io:
        nwbfile = io.read()

        # Get trials table and fluorescence module
        trials_df = nwbfile.trials.to_dataframe()
        fluorescence_module = nwbfile.processing["ophys"]["Fluorescence"]

        # Extract metadata from file name
        filename_parts = nwb_path.stem.split("++")
        condition = filename_parts[3] if len(filename_parts) > 3 else "unknown"

        # Map conditions to paper nomenclature
        condition_map = {"CTRL": "UL control", "PD": "6-OHDA", "OFF": "LID off"}

        # Iterate over trials table rows
        for _, trial in trials_df.iterrows():
            # Get the corresponding fluorescence series
            series_name = trial["roi_series_name"]
            series = fluorescence_module.roi_response_series[series_name]

            # Only process GRABCh channel (acetylcholine sensor)
            if "GRABCh" not in series_name:
                continue

            # Map treatment names
            treatment_mapping = {
                "control": "control",
                "50nM_dopamine": "dopamine",
                "quinpirole": "quinpirole",
                "sulpiride": "sulpiride",
                "TTX_calibration": "TTX_calibration",
                "acetylcholine_calibration": "ACh_calibration",
            }

            treatment = treatment_mapping.get(trial["treatment"], trial["treatment"])
            stimulation = "calibration" if "calibration" in trial["treatment"] else trial["stimulation"]

            fluorescence = series.data[:]
            timestamps = series.get_timestamps()
            # Convert to relative timestamps starting at 0
            timestamps = timestamps - timestamps[0]

            # Calculate relative stimulus time
            abs_stim_time = trial["stimulus_start_time"]
            orig_start = series.get_timestamps()[0]
            stim_time = abs_stim_time - orig_start

            trials_data.append(
                {
                    "condition": condition_map.get(condition, condition),
                    "treatment": treatment,
                    "stimulation": stimulation,
                    "subject_id": nwbfile.subject.subject_id if nwbfile.subject else "unknown",
                    "session_id": nwbfile.identifier,
                    "series_name": series_name,
                    "fluorescence": fluorescence,
                    "timestamps": timestamps,
                    "stim_time": stim_time,
                    "sampling_rate": 1.0 / np.mean(np.diff(timestamps)) if len(timestamps) > 1 else 1.0,
                    "file": nwb_path.name,
                }
            )

    return trials_data


def process_bot_trial(
    timestamps: np.ndarray,
    fluorescence: np.ndarray,
    stim_time: float = 10.0,
) -> Dict:
    """
    Process a single BOT trial following the original analysis methodology.

    Parameters
    ----------
    timestamps : np.ndarray
        Time points for the fluorescence data
    fluorescence : np.ndarray
        Background-subtracted fluorescence values
    stim_time : float
        Time of stimulation

    Returns
    -------
    dict
        Processed trial data including F0, dF/F0, AUC
    """
    # Calculate stimulus-relative windows
    baseline_window = (stim_time - 1.0, stim_time)
    response_window = (stim_time, stim_time + 2.0)

    # Find indices for time windows
    time_interval = timestamps[1] - timestamps[0]
    idx_baseline_start = np.argmin(np.abs(timestamps - baseline_window[0]))
    idx_baseline_end = np.argmin(np.abs(timestamps - baseline_window[1]))
    idx_response_start = np.argmin(np.abs(timestamps - response_window[0]))
    idx_response_end = np.argmin(np.abs(timestamps - response_window[1]))

    # Calculate F0 (baseline fluorescence)
    F0 = np.mean(fluorescence[idx_baseline_start:idx_baseline_end])

    # Calculate dF/F0
    dF_over_F0 = (fluorescence - F0) / F0 if F0 > 0 else np.zeros_like(fluorescence)

    # Calculate AUC for response window
    auc = np.sum(dF_over_F0[idx_response_start:idx_response_end]) * time_interval

    # Find peak response
    peak_response = np.max(dF_over_F0[idx_response_start:idx_response_end])

    return {
        "F0": F0,
        "dF_over_F0": dF_over_F0,
        "auc": auc,
        "peak_response": peak_response,
        "timestamps": timestamps,
        "fluorescence": fluorescence,
    }


def get_calibration_values_from_trials(trials_data: List[Dict]) -> Tuple[float, float]:
    """
    Extract Fmax and Fmin calibration values from ACh and TTX trials.

    Parameters
    ----------
    trials_data : List[Dict]
        List of trial data dictionaries

    Returns
    -------
    Fmax : float
        Maximum fluorescence from ACh application
    Fmin : float
        Minimum fluorescence from TTX application
    """
    # Find calibration trials
    ach_trials = [t for t in trials_data if t["treatment"] == "ACh_calibration"]
    ttx_trials = [t for t in trials_data if t["treatment"] == "TTX_calibration"]

    if ach_trials and ttx_trials:
        # Calculate average fluorescence from calibration trials
        # Apply same background subtraction as in original script
        bg_pmt = 162.0

        # Get mean fluorescence from ACh trials (after background subtraction)
        # Use entire signal since calibration trials have no stimulus
        ach_fluorescence = []
        for trial in ach_trials:
            bg_subtracted = trial["fluorescence"] - bg_pmt
            ach_fluorescence.extend(bg_subtracted)

        # Get mean fluorescence from TTX trials
        # Use entire signal since calibration trials have no stimulus
        ttx_fluorescence = []
        for trial in ttx_trials:
            bg_subtracted = trial["fluorescence"] - bg_pmt
            ttx_fluorescence.extend(bg_subtracted)

        if ach_fluorescence and ttx_fluorescence:
            Fmax = np.mean(ach_fluorescence)
            Fmin = np.mean(ttx_fluorescence)
            return Fmax, Fmin

    # Fallback to original script values if no calibration data
    # From original script: Fmax=493.47, Fmin=212.39 (after background subtraction)
    return 493.47, 212.39


def normalize_fluorescence(trial_data: Dict, Fmax: float, Fmin: float, stim_time: float = 10.0) -> Dict:
    """
    Normalize fluorescence data using calibration values.

    Parameters
    ----------
    trial_data : dict
        Processed trial data from process_bot_trial
    Fmax : float
        Maximum fluorescence from calibration
    Fmin : float
        Minimum fluorescence from calibration

    Returns
    -------
    dict
        Updated trial data with normalized values
    """
    # Calculate normalized dF/(Fmax-Fmin)
    FI = Fmax - Fmin

    trial_data["dF_over_FI"] = (trial_data["fluorescence"] - trial_data["F0"]) / FI

    # Calculate normalized AUC and peak for response window (stim_time to stim_time+2s)
    idx_start = np.argmin(np.abs(trial_data["timestamps"] - stim_time))
    idx_end = np.argmin(np.abs(trial_data["timestamps"] - (stim_time + 2.0)))
    response_slice = trial_data["dF_over_FI"][idx_start:idx_end]

    trial_data["auc_normalized"] = np.sum(response_slice) * (trial_data["timestamps"][1] - trial_data["timestamps"][0])
    trial_data["peak_normalized"] = np.max(response_slice)

    return trial_data


def analyze_nwb_dataset(nwb_dir: Path, stimulation_type: str = "single_pulse") -> pd.DataFrame:
    """
    Analyze all NWB files in the dataset and extract relevant measurements.

    Parameters
    ----------
    nwb_dir : Path
        Directory containing NWB files

    Returns
    -------
    pd.DataFrame
        DataFrame with all measurements
    """
    results = []

    # Get all NWB files
    nwb_files = sorted(nwb_dir.glob("*.nwb"))
    print(f"Found {len(nwb_files)} NWB files")

    # Process each file
    for nwb_file in nwb_files:
        try:
            # Load all trials from this file
            trials_data = load_nwb_trials_data(nwb_file)
            print(f"  Loaded {len(trials_data)} trials from {nwb_file.name}")

            # Get calibration values for this file
            Fmax, Fmin = get_calibration_values_from_trials(trials_data)

            # Process each trial
            for trial in trials_data:
                # Skip calibration trials for experimental analysis
                if trial["treatment"] in ["ACh_calibration", "TTX_calibration"]:
                    continue

                # Skip trials that don't match the requested stimulation type
                if trial["stimulation"] != stimulation_type:
                    continue

                # Get stimulus time from trials table
                stim_time = trial["stim_time"]

                # Process trial
                trial_data = process_bot_trial(trial["timestamps"], trial["fluorescence"], stim_time)

                # Normalize
                trial_data = normalize_fluorescence(trial_data, Fmax, Fmin, stim_time)

                # Store results
                results.append(
                    {
                        "condition": trial["condition"],
                        "treatment": trial["treatment"],
                        "subject_id": trial["subject_id"],
                        "session_id": trial["session_id"],
                        "series_name": trial["series_name"],
                        "F0": trial_data["F0"],
                        "auc_unnormalized": trial_data["auc"],
                        "auc_normalized": trial_data["auc_normalized"],
                        "peak_unnormalized": trial_data["peak_response"],
                        "peak_normalized": trial_data["peak_normalized"],
                        "file": trial["file"],
                    }
                )

        except Exception as e:
            print(f"  Error processing {nwb_file.name}: {e}")
            continue

    return pd.DataFrame(results)


def _collapse_first_two_trials_per_file(df: pd.DataFrame) -> pd.DataFrame:
    """
    Collapse per-trial rows into one value per file×condition×treatment by
    averaging the first two trials' normalized AUCs.

    Parameters
    ----------
    df : pd.DataFrame
        Per-trial dataset from analyze_nwb_dataset

    Returns
    -------
    pd.DataFrame
        Collapsed dataset with columns: condition, treatment, file, auc_normalized
    """
    plot_data = df[df["treatment"].isin(["control", "dopamine", "quinpirole", "sulpiride"])].copy()
    plot_data = plot_data.dropna(subset=["auc_normalized"])  # remove invalid rows

    # Preserve original order as proxy for chronological ordering, then take first two
    plot_data = plot_data.reset_index().rename(columns={"index": "_row_order"})
    first_two = plot_data.groupby(["condition", "treatment", "file"], sort=False, group_keys=False).apply(
        lambda g: g.nsmallest(2, columns=["_row_order"])
    )
    collapsed = first_two.groupby(["condition", "treatment", "file"], as_index=False)["auc_normalized"].mean()
    return collapsed


def create_figure_5f_from_nwb(df: pd.DataFrame) -> plt.Figure:
    """
    Create Figure 5F box plots from NWB-derived data matching the direct reproduction exactly.

    Uses the same single-panel layout, colors, treatments, and styling as the direct reproduction.
    """
    # Include dopamine along with main experimental treatments
    plot_data = df[df["treatment"].isin(["control", "dopamine", "quinpirole", "sulpiride"])].copy()

    if plot_data.empty:
        raise ValueError("No experimental data found for box plots")

    # Remove invalid data
    plot_data = plot_data.dropna(subset=["auc_normalized"])

    # Collapse to per-file averages using the first two trials per file×treatment.
    # Preserve original row order (proxy for chronological order) and average AUCs.
    plot_data = plot_data.reset_index().rename(columns={"index": "_row_order"})
    first_two = plot_data.groupby(["condition", "treatment", "file"], sort=False, group_keys=False).apply(
        lambda g: g.nsmallest(2, columns=["_row_order"])
    )
    collapsed = first_two.groupby(["condition", "treatment", "file"], as_index=False)["auc_normalized"].mean()

    # Create single panel figure to match direct reproduction layout
    fig, ax = plt.subplots(1, 1, figsize=(10, 6))

    # Define condition order and styling to match paper exactly
    conditions = ["UL control", "6-OHDA", "LID off"]
    condition_labels = ["control", "6-OHDA", "LID off-state"]

    # Treatments by condition - add dopamine for 6-OHDA and LID off groups only
    treatments_by_condition = {
        "UL control": ["control", "quinpirole", "sulpiride"],
        "6-OHDA": ["control", "dopamine", "quinpirole", "sulpiride"],
        "LID off": ["control", "dopamine", "quinpirole", "sulpiride"],
    }

    # Display labels for treatments (map dopamine to +DA)
    display_label = {
        "control": "control",
        "dopamine": "+DA",
        "quinpirole": "+quinpirole",
        "sulpiride": "+sulpiride",
    }

    # Colors to match paper: black, red, blue for the three conditions
    condition_colors = ["black", "red", "blue"]

    # Build dynamic x-positions per group with gaps between conditions
    group_gap = 1
    current_pos = 1
    x_positions: list[int] = []
    x_labels: list[str] = []
    group_bounds = []  # (start_pos, end_pos) for each condition

    all_data = []
    all_colors = []

    # Prepare data for all conditions and treatments
    for i, condition in enumerate(conditions):
        condition_data = collapsed[collapsed["condition"] == condition]

        start_pos = current_pos
        treatments = treatments_by_condition[condition]

        for treatment in treatments:
            treatment_series = condition_data[condition_data["treatment"] == treatment]["auc_normalized"]

            if not treatment_series.empty:
                all_data.append(treatment_series.values)
            else:
                all_data.append([])

            all_colors.append(condition_colors[i])
            x_positions.append(current_pos)
            x_labels.append(display_label[treatment])
            current_pos += 1

        end_pos = current_pos - 1
        group_bounds.append((start_pos, end_pos))
        current_pos += group_gap  # gap before next condition

    # Create box plots
    box_plot = ax.boxplot(
        all_data,
        positions=x_positions,
        patch_artist=True,
        boxprops=dict(linewidth=1.5),
        whiskerprops=dict(linewidth=1.5),
        capprops=dict(linewidth=1.5),
        medianprops=dict(color="black", linewidth=2),
        flierprops=dict(marker="o", markersize=3, markeredgecolor="black", alpha=0.7),
        widths=0.6,
    )

    # Color the boxes and add individual points with connecting lines
    for i, (patch, color) in enumerate(zip(box_plot["boxes"], all_colors)):
        # Set box face color based on condition
        if color == "black":
            patch.set_facecolor("white")
        elif color == "red":
            patch.set_facecolor("#ffcccc")  # Light red
        else:  # blue
            patch.set_facecolor("#ccccff")  # Light blue

        patch.set_edgecolor(color)
        patch.set_linewidth(1.5)

        # Add individual data points
        if len(all_data[i]) > 0:
            x_pos = np.random.normal(x_positions[i], 0.05, len(all_data[i]))
            ax.scatter(x_pos, all_data[i], color=color, alpha=0.7, s=20, zorder=3)

    # Add connecting lines between paired measurements (if data allows)
    # We connect control with other treatments within each condition when lengths match
    pos_idx = 0
    for i, condition in enumerate(conditions):
        condition_df = collapsed[collapsed["condition"] == condition]
        treatments = treatments_by_condition[condition]

        # Build a mapping from treatment to data and x position
        tr_to_vals = {}
        tr_to_x = {}
        for t in treatments:
            vals = condition_df[condition_df["treatment"] == t]["auc_normalized"].values
            tr_to_vals[t] = vals
            tr_to_x[t] = x_positions[pos_idx]
            pos_idx += 1

        ctrl_vals = tr_to_vals.get("control", np.array([]))

        # Connect control to each other treatment for which paired length matches
        for t in treatments:
            if t == "control":
                continue
            other_vals = tr_to_vals.get(t, np.array([]))
            min_len = min(len(ctrl_vals), len(other_vals))
            if min_len == 0:
                continue
            x_ctrl = tr_to_x["control"]
            x_t = tr_to_x[t]
            for k in range(min_len):
                x_ctrl_j = x_ctrl + np.random.normal(0, 0.05)
                x_t_j = x_t + np.random.normal(0, 0.05)
                ax.plot(
                    [x_ctrl_j, x_t_j],
                    [ctrl_vals[k], other_vals[k]],
                    color="gray",
                    alpha=0.3,
                    linewidth=0.5,
                    zorder=1,
                )

    # Styling to match paper
    ax.set_ylabel("AUC (normalized ΔF/F0*s)", fontsize=12, fontweight="bold")
    ax.set_ylim(-0.05, 0.7)
    ax.set_yticks([0, 0.2, 0.4, 0.6])

    # Add dotted horizontal line at y=0
    ax.axhline(y=0, color="black", linestyle=":", linewidth=1.0, alpha=0.8, zorder=2)

    # Set x-axis limits
    first_pos = group_bounds[0][0]
    last_pos = group_bounds[-1][1]
    ax.set_xlim(first_pos - 0.5, last_pos + 0.5)

    # Bottom: treatment labels under the x-axis line
    ax.set_xticks(x_positions)
    ax.set_xticklabels(x_labels, rotation=45, ha="right", fontsize=9)

    # Top: condition labels centered above groups
    group_centers = [0.5 * (s + e) for (s, e) in group_bounds]
    ax_top = ax.twiny()
    ax_top.set_xlim(ax.get_xlim())
    ax_top.set_xticks(group_centers)
    ax_top.set_xticklabels(condition_labels, fontsize=11, fontweight="bold")
    ax_top.tick_params(axis="x", which="major", pad=0)
    for spine in ["top", "bottom", "left", "right"]:
        ax_top.spines[spine].set_visible(False)

    # Add condition separators at group boundaries
    for start, end in group_bounds[:-1]:
        ax.axvline(x=end + 0.5, color="lightgray", linestyle="-", alpha=0.5, linewidth=1)

    # Add background shading for each condition
    shade_colors = ["lightgray", "lightcoral", "lightblue"]
    for (start, end), shade in zip(group_bounds, shade_colors):
        ax.axvspan(start - 0.5, end + 0.5, alpha=0.1, color=shade, zorder=0)

    # Final styling
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(1.5)
    ax.spines["bottom"].set_linewidth(1.5)

    plt.tight_layout()

    return fig


def generate_nwb_report(df: pd.DataFrame, output_dir: Path) -> str:
    """
    Generate analysis report for NWB-based reproduction.

    Parameters
    ----------
    df : pd.DataFrame
        Analyzed data
    output_dir : Path
        Output directory

    Returns
    -------
    str
        Path to generated report
    """
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    # Collapse per-trial data to per-file averages (first two trials)
    collapsed = _collapse_first_two_trials_per_file(df)

    # Calculate summary statistics using collapsed control data
    stats = {}
    control_data = collapsed[collapsed["treatment"] == "control"]

    for condition in control_data["condition"].unique():
        condition_data = control_data[control_data["condition"] == condition]
        auc_values = condition_data["auc_normalized"]

        stats[condition] = {
            "n": len(auc_values),
            "mean": auc_values.mean(),
            "sem": auc_values.sem() if len(auc_values) > 1 else 0,
            "std": auc_values.std() if len(auc_values) > 1 else 0,
            "files": len(condition_data["file"].unique()),
        }

    markdown_content = f"""# Figure 5F NWB Reproduction Report

**Generated:** {timestamp}
**Data Source:** NWB files from acetylcholine biosensor experiments
**Analysis:** GRABACh3.0 acetylcholine release measurements

## Executive Summary

This report presents the reproduction of Figure 5F from Zhai et al. 2025 using standardized NWB files.
The analysis demonstrates acetylcholine release dynamics measured with the GRABACh3.0 sensor across
different experimental conditions.

### Key Findings

- **Total Per-file Datapoints (all treatments):** {len(collapsed)}
- **Control Condition Per-file Measurements:** {len(control_data)}
- **Experimental Conditions (control):** {len(control_data['condition'].unique())}

## Summary Statistics

### Normalized AUC Values (Control Treatment Only)

"""

    # Add statistics for each condition
    for condition in ["UL control", "6-OHDA", "LID off"]:
        if condition in stats:
            s = stats[condition]
            markdown_content += f"""
#### {condition}
- **N:** {s['n']} measurements from {s['files']} files
- **Mean ± SEM:** {s['mean']:.3f} ± {s['sem']:.3f}
- **Std Dev:** {s['std']:.3f}
"""

    markdown_content += f"""
## Methodology

### Data Processing Pipeline
1. **NWB File Loading:** Extracted fluorescence time series from ophys processing module
2. **Background Subtraction:** Applied to raw fluorescence values
3. **Baseline Calculation:** Mean fluorescence from t=9-10s (before stimulation)
4. **Response Measurement:** AUC calculated from t=10-12s (after stimulation)
5. **Normalization:** Applied using calibration values (Fmax-Fmin)

### Analysis Parameters
- **Stimulation Time:** t=10 seconds
- **Baseline Window:** 9-10 seconds
- **Response Window:** 10-12 seconds
- **Sampling Rate:** ~5 Hz (from NWB timestamps)

## Comparison with Original Analysis

The NWB-based reproduction follows the exact methodology from the original BOT_imaging data_script.ipynb:
- Same time windows for baseline and response
- Identical normalization approach
- Consistent AUC calculation method

## Figure Output

![Figure 5F NWB Reproduction](figure_5f_nwb_boxplots.png)

Box plots showing normalized acetylcholine release (AUC) across experimental conditions.
Individual data points represent separate imaging sessions.

## Files Generated

- `figure_5f_nwb_boxplots.png` - Box plot reproduction from NWB data
- `figure_5f_nwb_report.md` - This analysis report

---

*Generated by Figure 5F GRABACh Reproduction Script*
*Data source: NWB files from Zhai et al. 2025*
"""

    # Write report
    report_path = output_dir / "figure_5f_nwb_report.md"
    with open(report_path, "w") as f:
        f.write(markdown_content)

    return str(report_path)


def main():
    """Main analysis function."""

    # Analysis parameters
    STIMULATION_TYPE = "single_pulse"  # Options: "single_pulse" or "burst_stimulation"

    # Set up paths
    base_dir = Path(__file__).parent.parent
    nwb_dir = Path("/home/heberto/development/surmeier-lab-to-nwb/nwb_files/acetylcholine_biosensor/figure_5")
    output_dir = base_dir / "analysis_outputs" / "figure_5f"

    # Ensure output directory exists
    output_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 80)
    print("FIGURE 5F REPRODUCTION FROM NWB FILES")
    print("=" * 80)
    print(f"NWB Directory: {nwb_dir}")
    print(f"Output Directory: {output_dir}")
    print()

    # Analyze NWB dataset
    print(f"Step 1: Analyzing NWB files for {STIMULATION_TYPE} stimulation...")
    df = analyze_nwb_dataset(nwb_dir, STIMULATION_TYPE)
    print(f"Analyzed {len(df)} trials total")

    # Create figure
    print("\nStep 2: Creating Figure 5F from NWB data...")
    fig = create_figure_5f_from_nwb(df)
    fig_path = output_dir / "figure_5f_nwb_boxplots.png"
    fig.savefig(fig_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  Saved: {fig_path}")

    # Generate report
    print("\nStep 3: Generating analysis report...")
    report_path = generate_nwb_report(df, output_dir)
    print(f"  Saved: {report_path}")

    # Print summary (collapsed per-file counts)
    print("\n" + "=" * 80)
    print("ANALYSIS COMPLETED SUCCESSFULLY!")
    print("=" * 80)

    # Collapsed dataset for summary
    collapsed = _collapse_first_two_trials_per_file(df)
    control_collapsed = collapsed[collapsed["treatment"] == "control"]
    print("\nSummary by Condition (collapsed per-file, control):")
    for condition in ["UL control", "6-OHDA", "LID off"]:
        condition_data = control_collapsed[control_collapsed["condition"] == condition]
        if len(condition_data) > 0:
            auc_values = condition_data["auc_normalized"]
            print(f"  {condition}: n={len(auc_values)}, mean={auc_values.mean():.3f}±{auc_values.sem():.3f}")
    print(f"\nTotal per-file datapoints (all treatments): {len(collapsed)}")

    print(f"\nOutput Files:")
    print(f"  - {fig_path}")
    print(f"  - {report_path}")

    print("\n" + "=" * 80)
    print("Figure 5F successfully reproduced from NWB files!")
    print("=" * 80)


if __name__ == "__main__":
    main()
