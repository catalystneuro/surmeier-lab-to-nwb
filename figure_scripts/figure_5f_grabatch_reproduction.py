#!/usr/bin/env python3
"""
Figure 5F GRABACh3.0 Acetylcholine Reproduction Script

This script reproduces Figure 5F from Zhai et al. 2025 using NWB files containing
GRABACh3.0 acetylcholine sensor BOT (Brightness Over Time) data. It loads two-photon
imaging data, calculates normalized ΔF/(Fmax-Fmin) and AUC values, and generates
publication-quality box plots.

The analysis adapts the original BOT_imaging data_script.ipynb to work with NWB
structured data, providing acetylcholine release measurements across different
pharmacological conditions and experimental groups.

Figure 5F shows normalized AUC values for:
- UL control (unlesioned control)
- PD (6-OHDA lesioned)
- LID off (dyskinetic mice, off-state)

With treatments: control, +quinpirole, +sulpiride

Usage:
    python figure_scripts/figure_5f_grabatch_reproduction.py

Output:
    - analysis_outputs/figure_5f/figure_5f_reproduction_report.md
    - analysis_outputs/figure_5f/figure_5f_grabatch_boxplots.png
    - analysis_outputs/figure_5f/figure_5f_traces.png
"""

import warnings
from datetime import datetime
from pathlib import Path
from typing import Dict, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pynwb import NWBHDF5IO
from tqdm import tqdm

# Suppress expected warnings when calculating statistics
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


def load_grabatch_data(nwb_dir: Path) -> pd.DataFrame:
    """
    Load GRABACh3.0 BOT data from NWB files.

    Parameters
    ----------
    nwb_dir : Path
        Directory containing Figure 5 GRABACh3.0 NWB files

    Returns
    -------
    pd.DataFrame
        Combined dataframe with all trials and fluorescence data
    """
    nwb_files = list(nwb_dir.glob("*.nwb"))

    if not nwb_files:
        raise FileNotFoundError(
            f"No GRABACh3.0 NWB files found in {nwb_dir}. "
            "Please run figure_5_acetylcholine_grabatch.py first to generate NWB files."
        )

    all_data = []
    print(f"Loading data from {len(nwb_files)} NWB files...")

    for nwb_file in tqdm(nwb_files, desc="Loading NWB files"):
        with NWBHDF5IO(str(nwb_file), "r") as io:
            nwb = io.read()

            # Check if file has stimulus table - if not, it's calibration and should be skipped
            if "stimulus_table" not in nwb.stimulus:
                continue

            # Extract metadata from session info
            session_id = nwb.session_id
            subject_id = nwb.subject.subject_id

            # Extract condition from session_id (now includes condition at beginning)
            # Session ID format: condition_YYYY-MM-DD_sliceXROIX_treatment_stimulation_sN
            parts = session_id.split("_")
            if len(parts) >= 6:
                # Handle multi-word conditions
                if parts[0] == "LID" and parts[1] == "off":
                    condition = "LID off"
                elif parts[0] == "UL" and parts[1] == "control":
                    condition = "UL control"
                else:
                    condition = parts[0]

                # Treatment is at position 4 (after date and slice info)
                treatment = parts[4]

                # Handle multi-word treatments
                if treatment == "50nM" and len(parts) > 5 and parts[5] == "dopamine":
                    treatment = "50nM_dopamine"
                elif treatment == "acetylcholine" and len(parts) > 5 and parts[5] == "calibration":
                    treatment = "acetylcholine_calibration"
                elif treatment == "TTX" and len(parts) > 5 and parts[5] == "calibration":
                    treatment = "TTX_calibration"

                # Get stimulation type from later in the parts
                if "single" in parts and "pulse" in parts:
                    stimulation = "single_pulse"
                elif "burst" in parts and "stimulation" in parts:
                    stimulation = "burst_stimulation"
                else:
                    stimulation = "unknown"

                # Map treatment names to match original analysis
                treatment_mapping = {
                    "control": "control",
                    "50nM_dopamine": "50nMDA",
                    "quinpirole": "quinpirole",
                    "sulpiride": "sulpiride",
                    "acetylcholine_calibration": "ACh",
                    "TTX_calibration": "TTX",
                }

                treatment = treatment_mapping.get(treatment, treatment)
            else:
                condition = "unknown"
                treatment = "unknown"
                stimulation = "unknown"

            # Filter for single pulse trials only using stimulus table
            stimulus_table = nwb.stimulus["stimulus_table"]
            if len(stimulus_table) == 0:
                continue

            # Check if this is a single pulse trial
            stimulus_type = stimulus_table["stimulus_type"].data[0]
            if stimulus_type != "single_pulse":
                continue

            # Skip non-single-pulse experimental trials (but include calibration)
            if treatment not in ["ACh", "TTX"] and stimulus_type != "single_pulse":
                continue

            # Access fluorescence data from acetylcholine channel (Ch2)
            if "TwoPhotonSeriesCh2" in nwb.acquisition:
                ophys_ts = nwb.acquisition["TwoPhotonSeriesCh2"]

                # Get fluorescence data and timestamps using the proper method
                fluorescence_data = ophys_ts.data[:]
                timestamps = ophys_ts.get_timestamps()

                # Extract acetylcholine ROI fluorescence (GRABACh3.0 sensor)
                # Note: Originally called "Region 1" in Prairie View, now acetylcholine channel
                if len(fluorescence_data.shape) > 1:
                    # If multi-dimensional, take mean over spatial dimensions (ROI average)
                    acetylcholine_fluorescence = np.mean(fluorescence_data, axis=(1, 2))
                else:
                    acetylcholine_fluorescence = fluorescence_data

                all_data.append(
                    {
                        "session_id": session_id,
                        "subject_id": subject_id,
                        "condition": condition,
                        "treatment": treatment,
                        "stimulation": stimulation,
                        "timestamps": timestamps,
                        "acetylcholine_fluorescence": acetylcholine_fluorescence,
                        "nwb_file": nwb_file.name,
                    }
                )

    # Special handling for calibration files that might not have stimulus tables
    # Load ACh and TTX calibration files separately
    for nwb_file in nwb_files:
        if "_ACh-" in nwb_file.name or "_TTX-" in nwb_file.name:
            with NWBHDF5IO(str(nwb_file), "r") as io:
                nwb = io.read()

                # Skip if already processed (has stimulus table)
                if "stimulus_table" in nwb.stimulus:
                    continue

                # Extract metadata
                session_id = nwb.session_id
                subject_id = nwb.subject.subject_id
                parts = session_id.split("_")

                # Parse condition
                if parts[0] == "LID" and parts[1] == "off":
                    condition = "LID off"
                elif parts[0] == "UL" and parts[1] == "control":
                    condition = "UL control"
                else:
                    condition = parts[0]

                # Determine treatment from filename
                if "_ACh-" in nwb_file.name:
                    treatment = "ACh"
                elif "_TTX-" in nwb_file.name:
                    treatment = "TTX"
                else:
                    continue

                stimulation = "calibration"

                # Access fluorescence data
                if "TwoPhotonSeriesCh2" in nwb.acquisition:
                    ophys_ts = nwb.acquisition["TwoPhotonSeriesCh2"]
                    fluorescence_data = ophys_ts.data[:]
                    timestamps = ophys_ts.get_timestamps()

                    if len(fluorescence_data.shape) > 1:
                        acetylcholine_fluorescence = np.mean(fluorescence_data, axis=(1, 2))
                    else:
                        acetylcholine_fluorescence = fluorescence_data

                    all_data.append(
                        {
                            "session_id": session_id,
                            "subject_id": subject_id,
                            "condition": condition,
                            "treatment": treatment,
                            "stimulation": stimulation,
                            "timestamps": timestamps,
                            "acetylcholine_fluorescence": acetylcholine_fluorescence,
                            "nwb_file": nwb_file.name,
                        }
                    )

    if not all_data:
        raise ValueError(
            f"No GRABACh3.0 fluorescence data could be loaded from NWB files. "
            f"Checked {len(nwb_files)} files in {nwb_dir}. "
            f"Looking for 'TwoPhotonSeriesCh2' in acquisition and filtering for single_pulse trials. "
            f"First few filenames: {[f.name for f in nwb_files[:3]]}"
        )

    print(f"Loaded {len(all_data)} trials from {len(set([d['condition'] for d in all_data]))} conditions")
    return pd.DataFrame(all_data)


def calculate_bot_metrics(
    timestamps: np.ndarray,
    acetylcholine_fluorescence: np.ndarray,
    baseline_start: float = 9.0,
    baseline_end: float = 10.0,
    analysis_start: float = 10.0,
    analysis_end: float = 12.0,
    pmt_background: float = 162.0,
) -> Dict[str, float]:
    """
    Calculate BOT metrics following the original notebook analysis.

    Adapted from the BOT_imaging data_script.ipynb methodology:
    1. Subtract PMT background
    2. Calculate F0 from baseline period (9-10s)
    3. Calculate ΔF/F0
    4. Calculate AUC from analysis window (10-12s)

    Parameters
    ----------
    timestamps : np.ndarray
        Time points in seconds
    acetylcholine_fluorescence : np.ndarray
        Raw GRABACh3.0 acetylcholine sensor fluorescence values
    baseline_start : float, default=9.0
        Start of baseline period (seconds)
    baseline_end : float, default=10.0
        End of baseline period (seconds)
    analysis_start : float, default=10.0
        Start of analysis window (seconds)
    analysis_end : float, default=12.0
        End of analysis window (seconds)
    pmt_background : float, default=162.0
        PMT background to subtract

    Returns
    -------
    Dict[str, float]
        Dictionary with F0, ΔF/F0 trace, AUC, and peak amplitude
    """
    # Subtract PMT background from acetylcholine channel (matching original analysis)
    acetylcholine_bg_subtracted = acetylcholine_fluorescence - pmt_background

    # Find indices for baseline and analysis windows
    baseline_mask = (timestamps >= baseline_start) & (timestamps <= baseline_end)
    analysis_mask = (timestamps >= analysis_start) & (timestamps <= analysis_end)

    # Calculate F0 (baseline acetylcholine fluorescence)
    acetylcholine_f0 = np.mean(acetylcholine_bg_subtracted[baseline_mask])

    # Calculate ΔF/F0 for acetylcholine channel
    acetylcholine_delta_f_over_f0 = (acetylcholine_bg_subtracted - acetylcholine_f0) / acetylcholine_f0

    # Calculate AUC in analysis window for acetylcholine channel
    acetylcholine_analysis_trace = acetylcholine_delta_f_over_f0[analysis_mask]
    analysis_times = timestamps[analysis_mask]

    if len(acetylcholine_analysis_trace) > 0:
        time_interval = np.mean(np.diff(analysis_times))
        acetylcholine_auc = np.sum(acetylcholine_analysis_trace) * time_interval
        acetylcholine_peak_amplitude = np.max(acetylcholine_analysis_trace)
    else:
        acetylcholine_auc = np.nan
        acetylcholine_peak_amplitude = np.nan

    return {
        "acetylcholine_f0": acetylcholine_f0,
        "acetylcholine_delta_f_over_f0": acetylcholine_delta_f_over_f0,
        "acetylcholine_auc": acetylcholine_auc,
        "acetylcholine_peak_amplitude": acetylcholine_peak_amplitude,
        "acetylcholine_bg_subtracted": acetylcholine_bg_subtracted,
        "timestamps": timestamps,
    }


def calculate_normalized_metrics(
    df: pd.DataFrame, fmax_treatment: str = "ACh", fmin_treatment: str = "TTX"
) -> pd.DataFrame:
    """
    Calculate normalized ΔF/(Fmax-Fmin) metrics using calibration trials.

    Follows the original notebook approach for normalization using
    ACh (Fmax calibration) and TTX (Fmin calibration) trials.

    Parameters
    ----------
    df : pd.DataFrame
        Trial data with calculated metrics
    fmax_treatment : str, default="ACh"
        Treatment name for Fmax calibration
    fmin_treatment : str, default="TTX"
        Treatment name for Fmin calibration

    Returns
    -------
    pd.DataFrame
        Data with normalized metrics added
    """
    results = []

    # Group by condition (UL control, PD, LID off) and session
    for (condition, session_id), session_data in df.groupby(["condition", "session_id"]):

        # Find calibration trials for this session
        fmax_trials = session_data[session_data["treatment"] == fmax_treatment]
        fmin_trials = session_data[session_data["treatment"] == fmin_treatment]

        # Calculate Fmax and Fmin if calibration trials exist
        if not fmax_trials.empty and not fmin_trials.empty:
            # Use the mean F0 from calibration trials (following notebook approach)
            fmax_final = fmax_trials["acetylcholine_f0"].mean()
            fmin_final = fmin_trials["acetylcholine_f0"].mean()

            # Calculate normalization factor
            normalization_factor = fmax_final - fmin_final
        else:
            # If no calibration trials, skip normalization for this session
            print(f"Warning: No calibration trials found for {condition} {session_id}")
            normalization_factor = None
            fmax_final = np.nan
            fmin_final = np.nan

        # Process all non-calibration trials for this session
        experimental_trials = session_data[~session_data["treatment"].isin([fmax_treatment, fmin_treatment])]

        for _, trial in experimental_trials.iterrows():
            trial_dict = trial.to_dict()

            # Add calibration info
            trial_dict["fmax_final"] = fmax_final
            trial_dict["fmin_final"] = fmin_final

            # Calculate normalized metrics if calibration is available
            if normalization_factor is not None and normalization_factor != 0:
                # Normalized ΔF/(Fmax-Fmin) following notebook methodology
                acetylcholine_bg_subtracted = trial["acetylcholine_fluorescence"] - 162.0  # PMT background
                acetylcholine_delta_f_normalized = (
                    acetylcholine_bg_subtracted - trial["acetylcholine_f0"]
                ) / normalization_factor

                # Calculate normalized AUC for acetylcholine channel
                timestamps = trial["timestamps"]
                analysis_mask = (timestamps >= 10.0) & (timestamps <= 12.0)

                if np.any(analysis_mask):
                    acetylcholine_analysis_trace = acetylcholine_delta_f_normalized[analysis_mask]
                    analysis_times = timestamps[analysis_mask]
                    time_interval = np.mean(np.diff(analysis_times))

                    acetylcholine_normalized_auc = np.sum(acetylcholine_analysis_trace) * time_interval
                    acetylcholine_normalized_peak = np.max(acetylcholine_analysis_trace)
                else:
                    acetylcholine_normalized_auc = np.nan
                    acetylcholine_normalized_peak = np.nan

                trial_dict["acetylcholine_normalized_auc"] = acetylcholine_normalized_auc
                trial_dict["acetylcholine_normalized_peak"] = acetylcholine_normalized_peak
                trial_dict["acetylcholine_delta_f_normalized"] = acetylcholine_delta_f_normalized
            else:
                trial_dict["acetylcholine_normalized_auc"] = np.nan
                trial_dict["acetylcholine_normalized_peak"] = np.nan
                trial_dict["acetylcholine_delta_f_normalized"] = None

            results.append(trial_dict)

    return pd.DataFrame(results)


def create_grabatch_boxplots(df: pd.DataFrame) -> Tuple[plt.Figure, plt.Axes]:
    """
    Create box plots for GRABACh3.0 normalized AUC data matching Figure 5F.

    Parameters
    ----------
    df : pd.DataFrame
        Trial data with normalized metrics

    Returns
    -------
    Tuple[plt.Figure, plt.Axes]
        Figure and axes objects
    """
    # Filter for main experimental treatments (exclude calibration)
    experimental_data = df[df["treatment"].isin(["control", "quinpirole", "sulpiride"])].copy()

    if experimental_data.empty:
        raise ValueError("No experimental data found for box plots")

    # Remove invalid data
    experimental_data = experimental_data.dropna(subset=["acetylcholine_normalized_auc"])

    fig, axes = plt.subplots(1, 3, figsize=(12, 5), sharey=True)

    # Define condition order and styling
    conditions = ["UL control", "PD", "LID off"]
    condition_titles = ["UL control", "PD", "LID off-state"]
    treatments = ["control", "quinpirole", "sulpiride"]
    treatment_labels = ["control", "+quinpirole", "+sulpiride"]

    # Colors to match paper style
    colors = ["white", "lightgray", "gray"]

    for i, (condition, title) in enumerate(zip(conditions, condition_titles)):
        ax = axes[i]

        condition_data = experimental_data[experimental_data["condition"] == condition]

        if condition_data.empty:
            print(f"Warning: No data for condition {condition}")
            ax.set_title(title, fontsize=12, fontweight="bold")
            ax.set_ylabel("acetylcholine normalized AUC" if i == 0 else "", fontsize=12)
            continue

        # Prepare data for box plot
        plot_data = []
        labels = []

        for treatment, label in zip(treatments, treatment_labels):
            treatment_data = condition_data[condition_data["treatment"] == treatment][
                "acetylcholine_normalized_auc"
            ].dropna()

            if not treatment_data.empty:
                plot_data.append(treatment_data.values)
                labels.append(label)
            else:
                plot_data.append([])
                labels.append(label)

        # Create box plot
        box_plot = ax.boxplot(
            plot_data,
            tick_labels=labels,
            patch_artist=True,
            boxprops=dict(linewidth=1),
            whiskerprops=dict(linewidth=1),
            capprops=dict(linewidth=1),
            medianprops=dict(color="black", linewidth=2),
            flierprops=dict(marker="o", markersize=4, markeredgecolor="black", alpha=0.7),
        )

        # Color the boxes
        for patch, color in zip(box_plot["boxes"], colors):
            patch.set_facecolor(color)
            patch.set_edgecolor("black")

        # Add individual data points
        for j, data in enumerate(plot_data):
            if len(data) > 0:
                x_pos = np.random.normal(j + 1, 0.04, len(data))
                ax.scatter(x_pos, data, color="black", alpha=0.6, s=15, zorder=3)

        # Styling
        ax.set_title(title, fontsize=12, fontweight="bold")
        ax.set_ylabel("acetylcholine normalized AUC" if i == 0 else "", fontsize=12)
        ax.set_ylim(-0.1, 0.8)
        ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8])

        # Rotate x-axis labels if needed
        ax.tick_params(axis="x", rotation=45)

    plt.suptitle("Figure 5F: GRABACh3.0 Acetylcholine Release", fontsize=14, fontweight="bold")
    plt.tight_layout()

    return fig, axes


def create_example_traces(df: pd.DataFrame) -> Tuple[plt.Figure, plt.Axes]:
    """
    Create example ΔF/F0 traces showing the analysis methodology.

    Parameters
    ----------
    df : pd.DataFrame
        Trial data with calculated metrics

    Returns
    -------
    Tuple[plt.Figure, plt.Axes]
        Figure and axes objects
    """
    fig, ax = plt.subplots(1, 1, figsize=(8, 6))

    # Find example traces from UL control condition
    control_data = df[
        (df["condition"] == "UL control") & (df["treatment"].isin(["control", "quinpirole", "sulpiride"]))
    ]

    if control_data.empty:
        ax.text(0.5, 0.5, "No example traces available", ha="center", va="center", transform=ax.transAxes)
        return fig, ax

    # Plot representative traces for each treatment
    colors = {"control": "black", "quinpirole": "green", "sulpiride": "orange"}

    for treatment in ["control", "quinpirole", "sulpiride"]:
        treatment_data = control_data[control_data["treatment"] == treatment]

        if not treatment_data.empty:
            # Take the first available trace
            example = treatment_data.iloc[0]
            timestamps = example["timestamps"]
            acetylcholine_delta_f_over_f0 = example["acetylcholine_delta_f_over_f0"]

            ax.plot(timestamps, acetylcholine_delta_f_over_f0, color=colors[treatment], linewidth=1.5, label=treatment)

    # Add analysis windows
    ax.axvspan(9, 10, alpha=0.2, color="blue", label="baseline (F0)")
    ax.axvspan(10, 12, alpha=0.2, color="red", label="analysis (AUC)")

    # Styling
    ax.set_xlabel("Time (s)", fontsize=12)
    ax.set_ylabel("acetylcholine ΔF/F0", fontsize=12)
    ax.set_xlim(8, 14)
    ax.legend(loc="upper right", frameon=False)
    ax.set_title("Example GRABACh3.0 Traces - Analysis Windows", fontsize=14)

    plt.tight_layout()
    return fig, ax


def generate_markdown_report(df: pd.DataFrame, output_dir: Path) -> str:
    """
    Generate comprehensive markdown report.

    Parameters
    ----------
    df : pd.DataFrame
        Trial data with calculated metrics
    output_dir : Path
        Output directory

    Returns
    -------
    str
        Path to generated report
    """
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    # Calculate summary statistics
    experimental_data = df[
        df["treatment"].isin(["control", "quinpirole", "sulpiride"]) & df["acetylcholine_normalized_auc"].notna()
    ]

    markdown_content = f"""# Figure 5F GRABACh3.0 Acetylcholine Analysis Report

**Generated:** {timestamp}
**Data Source:** NWB files with GRABACh3.0 BOT imaging data
**Analysis:** Normalized ΔF/(Fmax-Fmin) and AUC measurements for acetylcholine release

## Executive Summary

This report presents the reproduction of Figure 5F from Zhai et al. 2025, demonstrating
acetylcholine release dynamics measured with the GRABACh3.0 sensor across different
experimental conditions and pharmacological treatments.

### Key Findings

- **Total Trials Analyzed:** {len(experimental_data)}
- **Experimental Conditions:** {len(experimental_data['condition'].unique())}
- **Pharmacological Treatments:** {len(experimental_data['treatment'].unique())}
- **Sessions:** {len(experimental_data['session_id'].unique())}

## Methodology

### Experimental Protocol
GRABACh3.0 acetylcholine sensor imaging with electrical stimulation:
- **Sensor:** Genetically encoded fluorescent acetylcholine sensor
- **Imaging:** Two-photon microscopy at 920nm excitation
- **Stimulation:** Electrical pulses (single or burst) at 3s post-baseline
- **Analysis:** BOT (Brightness Over Time) quantification

### Analysis Approach
Following the original BOT_imaging data_script.ipynb methodology:

1. **Background Subtraction:** PMT background (162 counts) subtracted
2. **Baseline Calculation:** F0 from 9-10s baseline period
3. **ΔF/F0 Calculation:** (F-F0)/F0 for raw measurements
4. **Calibration:** Fmax (ACh) and Fmin (TTX) for normalization
5. **Normalization:** ΔF/(Fmax-Fmin) for sensor saturation correction
6. **AUC Calculation:** Area under curve from 10-12s analysis window

## Results by Condition

"""

    # Add condition-specific results
    for condition in experimental_data["condition"].unique():
        condition_data = experimental_data[experimental_data["condition"] == condition]

        markdown_content += f"""### {condition}

- **Sample Size:** {len(condition_data)} trials from {len(condition_data['session_id'].unique())} sessions
- **Treatments:** {', '.join(condition_data['treatment'].unique())}

**Treatment Summary:**
"""

        for treatment in condition_data["treatment"].unique():
            treatment_data = condition_data[condition_data["treatment"] == treatment]
            auc_data = treatment_data["acetylcholine_normalized_auc"].dropna()

            if not auc_data.empty:
                markdown_content += f"""
- **{treatment}:** n={len(auc_data)}, AUC={auc_data.mean():.3f}±{auc_data.sem():.3f} (mean±SEM)"""

    markdown_content += f"""

## Calibration Analysis

### Fmax/Fmin Calibration
- **Fmax (ACh):** Sensor saturation with 100μM acetylcholine
- **Fmin (TTX):** Activity blockade with 10μM tetrodotoxin
- **Purpose:** Normalize for sensor expression and background fluorescence

## Figures

### Box Plot Analysis
![GRABACh3.0 Box Plots](figure_5f_grabatch_boxplots.png)

**Normalized AUC Comparison:** Statistical comparison of acetylcholine release
across experimental conditions and pharmacological treatments.

### Example Traces
![Example Traces](figure_5f_traces.png)

**Analysis Methodology:** Representative ΔF/F0 traces showing baseline and
analysis windows used for AUC calculations.

## Data Quality

- **Valid Trials:** {len(experimental_data)} trials with complete datasets
- **Calibration:** Fmax/Fmin normalization applied per session
- **Time Resolution:** {'N/A (no valid data)' if experimental_data.empty else f"{experimental_data.iloc[0]['timestamps'][1] - experimental_data.iloc[0]['timestamps'][0]:.3f}s"} sampling
- **Analysis Window:** 2-second post-stimulus period (10-12s)

## Statistical Notes

- **Normalization:** ΔF/(Fmax-Fmin) accounts for sensor expression variability
- **AUC Calculation:** Numerical integration over 2-second analysis window
- **Background Correction:** PMT dark current subtracted (162 counts)
- **Baseline:** F0 calculated from 1-second pre-stimulus period (9-10s)

## Files Generated

- `figure_5f_grabatch_boxplots.png` - Box plot statistical comparison
- `figure_5f_traces.png` - Example traces with analysis windows
- `figure_5f_reproduction_report.md` - This detailed report

---

*Report generated by Figure 5F GRABACh3.0 Reproduction Script*
*Analysis follows BOT_imaging data_script.ipynb methodology*
"""

    # Write report
    report_path = output_dir / "figure_5f_reproduction_report.md"
    with open(report_path, "w") as f:
        f.write(markdown_content)

    return str(report_path)


def main():
    """Main function to orchestrate the Figure 5F reproduction analysis."""

    # Set up paths
    base_dir = Path(__file__).parent.parent
    nwb_dir = base_dir / "nwb_files" / "acetylcholine_biosensor" / "figure_5"
    output_dir = base_dir / "analysis_outputs" / "figure_5f"

    # Ensure output directory exists
    output_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 80)
    print("FIGURE 5F GRABATCH3.0 ACETYLCHOLINE REPRODUCTION ANALYSIS")
    print("=" * 80)
    print(f"NWB Directory: {nwb_dir}")
    print(f"Output Directory: {output_dir}")
    print()

    # Load GRABACh3.0 data
    print("Step 1: Loading GRABACh3.0 BOT data from NWB files...")
    raw_df = load_grabatch_data(nwb_dir)

    # Calculate basic BOT metrics for all trials
    print("Step 2: Calculating BOT metrics (F0, ΔF/F0, AUC)...")
    bot_results = []

    for _, trial in raw_df.iterrows():
        metrics = calculate_bot_metrics(trial["timestamps"], trial["acetylcholine_fluorescence"])

        trial_dict = trial.to_dict()
        trial_dict.update(metrics)
        bot_results.append(trial_dict)

    bot_df = pd.DataFrame(bot_results)

    # Calculate normalized metrics using calibration trials
    print("Step 3: Calculating normalized metrics using Fmax/Fmin calibration...")
    df = calculate_normalized_metrics(bot_df)

    if df.empty:
        raise ValueError("No valid data after normalization")

    # Create box plots
    print("Step 4: Creating GRABACh3.0 box plots...")
    fig1, _ = create_grabatch_boxplots(df)
    fig1_path = output_dir / "figure_5f_grabatch_boxplots.png"
    fig1.savefig(fig1_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig1)
    print(f"  Saved: {fig1_path}")

    # Create example traces
    print("Step 5: Creating example traces...")
    fig2, _ = create_example_traces(df)
    fig2_path = output_dir / "figure_5f_traces.png"
    fig2.savefig(fig2_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig2)
    print(f"  Saved: {fig2_path}")

    # Generate report
    print("Step 6: Generating analysis report...")
    report_path = generate_markdown_report(df, output_dir)
    print(f"  Saved: {report_path}")

    # Print summary
    print("\n" + "=" * 80)
    print("ANALYSIS COMPLETED SUCCESSFULLY!")
    print("=" * 80)

    experimental_data = df[
        df["treatment"].isin(["control", "quinpirole", "sulpiride"]) & df["acetylcholine_normalized_auc"].notna()
    ]

    print(f"Total Trials Analyzed: {len(experimental_data)}")
    print(f"Experimental Conditions: {len(experimental_data['condition'].unique())}")
    print(f"Sessions: {len(experimental_data['session_id'].unique())}")

    print("\nCondition Summary:")
    for condition in experimental_data["condition"].unique():
        condition_data = experimental_data[experimental_data["condition"] == condition]
        print(f"  {condition}: n={len(condition_data)} trials")

        for treatment in condition_data["treatment"].unique():
            treatment_data = condition_data[condition_data["treatment"] == treatment]
            auc_data = treatment_data["acetylcholine_normalized_auc"].dropna()
            if not auc_data.empty:
                print(f"    {treatment}: n={len(auc_data)}, " f"AUC={auc_data.mean():.3f}±{auc_data.sem():.3f}")

    print(f"\nOutput Files:")
    print(f"  - {fig1_path}")
    print(f"  - {fig2_path}")
    print(f"  - {report_path}")

    print("\n" + "=" * 80)
    print("Figure 5F reproduction demonstrates NWB-powered GRABACh3.0 analysis!")
    print("Normalized acetylcholine release measurements with calibration correction.")
    print("=" * 80)


if __name__ == "__main__":
    main()
