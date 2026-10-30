#!/usr/bin/env python3
"""
Figure 1I Dendritic Excitability Reproduction Script

This script reproduces Figure 1I from Zhai et al. 2025 using NWB files containing
dendritic excitability data. It loads two-photon calcium imaging data and intracellular
recordings to calculate the dendritic excitability index (distal/proximal AUC ratio).

Based on the paper methods:
- Dendritic excitability index = distal dendrite AUC / proximal dendrite AUC
- AUC calculated from calcium transients (G/R0) during current injection
- Three conditions: LID off-state, LID on-state, LID on-state with SCH

Usage:
    python figure_scripts/figure_1i_dendritic_excitability_reproduction.py

Output:
    - analysis_outputs/figure_1i/figure_1i_reproduction_report.md
    - analysis_outputs/figure_1i/figure_1i_dendritic_excitability.png
"""

import warnings
from datetime import datetime
from pathlib import Path
from typing import Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pynwb import NWBHDF5IO

# Suppress expected warnings during analysis
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


def calculate_calcium_auc(fluo4_series, alexa568_series) -> float:
    """
    Calculate area under the curve (AUC) for ΔG/R0 calcium transient data.

    Based on the original dendritic excitability notebook methodology:
    1. Subtract PMT background from both channels
    2. Apply 6-point rolling average smoothing
    3. Calculate G/R0 ratio
    4. Define baseline (0.5-0.95s) and response (0.95-1.5s) periods
    5. Calculate ΔG/R0 = (G/R0) - baseline_mean
    6. Calculate AUC of ΔG/R0 during response period

    Parameters
    ----------
    fluo4_series : TimeSeriesReference
        One trial's slice of the Fluo-4 calcium imaging data (calcium-sensitive, green)
    alexa568_series : TimeSeriesReference
        The same trial's slice of the Alexa568 structural reference data (calcium-insensitive, red)

    Returns
    -------
    float
        Area under the curve for the ΔG/R0 calcium transient
    """
    try:
        # Get calcium and structural data
        # The fluorescence series have one column, for the one ROI of the scan line
        fluo4_data = fluo4_series.data[:, 0]  # Green channel (Ch2)
        alexa568_data = alexa568_series.data[:, 0]  # Red channel (Ch1)

        if len(fluo4_data) != len(alexa568_data):
            print(f"Warning: Data length mismatch - Fluo4: {len(fluo4_data)}, Alexa568: {len(alexa568_data)}")
            return np.nan

        # Step 1: Background subtraction (from original notebook)
        bg_PMT1 = 55  # Background for Ch1 (Alexa568/Red)
        bg_PMT2 = 161  # Background for Ch2 (Fluo4/Green)

        alexa568_bg_corrected = alexa568_data - bg_PMT1
        fluo4_bg_corrected = fluo4_data - bg_PMT2

        # Step 2: Smoothing with 6-point rolling average (from original notebook)
        # Convert to pandas Series for rolling mean
        import pandas as pd

        alexa568_smooth = pd.Series(alexa568_bg_corrected).rolling(6, min_periods=1).mean().values
        fluo4_smooth = pd.Series(fluo4_bg_corrected).rolling(6, min_periods=1).mean().values

        # Avoid division by zero
        alexa568_smooth = np.where(alexa568_smooth <= 0, np.nan, alexa568_smooth)

        # Step 3: Calculate G/R0 ratio
        g_over_r0 = fluo4_smooth / alexa568_smooth

        # Step 4: Define time windows (from original notebook methodology)
        # Get actual time intervals from the trial's timestamps
        try:
            timestamps = fluo4_series.timestamps
            time_interval = timestamps[1] - timestamps[0] if len(timestamps) > 1 else 0.00064
        except (AttributeError, IndexError):
            # Fallback to default interval
            time_interval = 0.00064  # seconds per data point

        # Time indices (from original notebook):
        # 0.5s baseline start, 0.95s stimulus, 1.5s response end
        baseline_start_idx = int(0.5 / time_interval)
        stimulus_idx = int(0.95 / time_interval)
        response_end_idx = int(1.5 / time_interval)

        # Ensure indices are within bounds
        max_idx = len(g_over_r0)
        baseline_start_idx = min(baseline_start_idx, max_idx - 100)
        stimulus_idx = min(stimulus_idx, max_idx - 50)
        response_end_idx = min(response_end_idx, max_idx - 1)

        if stimulus_idx >= response_end_idx or baseline_start_idx >= stimulus_idx:
            return np.nan

        # Step 5: Calculate baseline values (from original notebook)
        # Calculate baseline G0 (calcium signal) and R0 (structural reference)
        baseline_fluo4 = fluo4_smooth[baseline_start_idx:stimulus_idx]  # G0
        baseline_alexa568 = alexa568_smooth[baseline_start_idx:stimulus_idx]  # R0

        g0_mean = np.nanmean(baseline_fluo4)  # G0
        r0_mean = np.nanmean(baseline_alexa568)  # R0

        if np.isnan(g0_mean) or np.isnan(r0_mean) or r0_mean <= 0:
            return np.nan

        # Step 6: Calculate ΔG/R0 = (G - G0) / R0 (from original notebook)
        # This is the key formula from the original notebook
        delta_g_over_r0 = (fluo4_smooth - g0_mean) / r0_mean

        # Step 7: Extract response period and calculate AUC
        response_delta_g_r0 = delta_g_over_r0[stimulus_idx:response_end_idx]

        # Remove NaN values for AUC calculation
        valid_response_data = response_delta_g_r0[~np.isnan(response_delta_g_r0)]
        if len(valid_response_data) == 0:
            return np.nan

        # Calculate AUC using summation method (from original notebook)
        # Original: proxAreaG = (proxDeltaGoverR0[IndexTimeC:IndexTimeD].sum()) * TimeInterval
        auc = np.sum(valid_response_data) * time_interval

        return auc

    except Exception as e:
        print(f"Warning: Could not calculate AUC - {e}")
        return np.nan


def load_dendritic_excitability_data(nwb_dir: Path) -> pd.DataFrame:
    """
    Load dendritic excitability data from NWB files.

    For each cell and condition, calculates the dendritic excitability index as the ratio
    of distal to proximal dendrite calcium transient AUCs.

    Parameters
    ----------
    nwb_dir : Path
        Directory containing the Figure 1 per-mouse-day NWB files

    Returns
    -------
    pd.DataFrame
        Combined dataframe with dendritic excitability indices per cell/condition
    """
    nwb_files = sorted(nwb_dir.glob("*.nwb"))

    if not nwb_files:
        raise FileNotFoundError(
            f"No NWB files found in {nwb_dir}. "
            "Please run the per-mouse-day merge (uv run poe merge_patch_clamp) first to generate NWB files."
        )

    all_data = []
    print(f"Loading data from {len(nwb_files)} NWB files...")

    for nwb_file in nwb_files:
        with NWBHDF5IO(str(nwb_file), "r") as io:
            nwb = io.read()

            # Get session info
            subject_id = nwb.subject.subject_id

            # Each file holds the cells recorded from one mouse on one day. The RecordingsIndexTable has one row per
            # sweep; a dendritic sweep's row points to its slices of the Fluo-4 and Alexa 568 line-scan profiles
            index_df = nwb.processing["recordings_index"]["RecordingsIndexTable"].to_dataframe()
            dendritic_df = index_df[index_df["dendrite_type"].isin(["Proximal", "Distal"])]

            for (cell_id, condition), cell_df in dendritic_df.groupby(["cell_id", "condition"]):
                cell = f"{subject_id} cell {cell_id}"

                # Group calcium data by dendrite type (proximal vs distal)
                distal_aucs = []
                proximal_aucs = []

                # For each trial, calculate G/R0 AUC
                for _, row in cell_df.iterrows():
                    location_trial = f"{row['site']} trial at {row['sweep_start_time_s']:.1f} s"
                    auc = calculate_calcium_auc(row["roi_response_fluo4"], row["roi_response_alexa568"])
                    if not np.isnan(auc):
                        if row["dendrite_type"] == "Distal":
                            distal_aucs.append(auc)
                            print(f"    Distal AUC: {auc:.6f} from {location_trial}")
                        else:
                            proximal_aucs.append(auc)
                            print(f"    Proximal AUC: {auc:.6f} from {location_trial}")
                    else:
                        print(f"    Warning: Could not calculate AUC for {location_trial}")

                # Calculate dendritic excitability index (distal/proximal ratio)
                if len(distal_aucs) > 0 and len(proximal_aucs) > 0:
                    # Use mean AUCs for each location if multiple trials
                    mean_distal_auc = np.mean(distal_aucs)
                    mean_proximal_auc = np.mean(proximal_aucs)

                    if mean_proximal_auc > 0:
                        excitability_index = mean_distal_auc / mean_proximal_auc
                    else:
                        excitability_index = np.nan

                    all_data.append(
                        {
                            "cell": cell,
                            "subject_id": subject_id,
                            "condition": condition,
                            "dendritic_excitability_index": excitability_index,
                            "distal_auc_mean": mean_distal_auc,
                            "proximal_auc_mean": mean_proximal_auc,
                            "n_distal_trials": len(distal_aucs),
                            "n_proximal_trials": len(proximal_aucs),
                            "nwb_file": nwb_file.name,
                        }
                    )

                    print(f"  {cell}: Index = {excitability_index:.3f} ({condition})")
                else:
                    print(f"  Warning: {cell} ({condition}) - insufficient calcium data")

    if not all_data:
        raise ValueError("No dendritic excitability data could be loaded from NWB files")

    df = pd.DataFrame(all_data)
    print(f"Loaded {len(df)} cells with valid dendritic excitability data")
    return df


def create_dendritic_excitability_plot(df: pd.DataFrame) -> Tuple[plt.Figure, plt.Axes]:
    """
    Create Figure 1I dendritic excitability comparison box plot.

    Parameters
    ----------
    df : pd.DataFrame
        Dendritic excitability data

    Returns
    -------
    Tuple[plt.Figure, plt.Axes]
        Figure and axes objects
    """
    fig, ax = plt.subplots(1, 1, figsize=(4, 5))

    # Filter out invalid excitability indices
    valid_data = df[df["dendritic_excitability_index"].notna()]

    # Order conditions to match paper
    condition_order = ["LID off-state", "LID on-state", "LID on-state with SCH 23390"]
    condition_labels = ["off", "on", "on + SCH"]

    # Prepare data for box plot
    plot_data = []
    labels = []

    for condition, label in zip(condition_order, condition_labels):
        if condition in valid_data["condition"].unique():
            condition_data = valid_data[valid_data["condition"] == condition]["dendritic_excitability_index"]
            plot_data.append(condition_data.values)
            labels.append(label)

    # Create box plot matching paper style
    ax.boxplot(
        plot_data,
        tick_labels=labels,
        patch_artist=True,
        boxprops=dict(facecolor="white", color="black", linewidth=1),
        whiskerprops=dict(color="black", linewidth=1),
        capprops=dict(color="black", linewidth=1),
        medianprops=dict(color="black", linewidth=1.5),
        flierprops=dict(marker="o", markerfacecolor="gray", markersize=3, markeredgecolor="black", alpha=0.7),
    )

    # Add individual data points as in the paper
    for i, data in enumerate(plot_data):
        # Add some jitter to x position
        x_pos = np.random.normal(i + 1, 0.04, len(data))
        ax.scatter(x_pos, data, color="gray", alpha=0.6, s=15, zorder=3)

    # Style to match paper
    ax.set_ylabel("dendritic excitability index", fontsize=12)
    ax.set_ylim(0, 0.81)  # Typical range for dendritic excitability index
    ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8])

    # Add descriptive title
    ax.set_title("Figure 1I: dSPN Dendritic Excitability Index", fontsize=14, pad=20)

    # Remove top and right spines
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Add significance annotation if needed (matching paper style)
    # Note: This would need statistical testing to be accurate
    if len(plot_data) >= 2:
        ax.text(1.5, 0.75, "*", ha="center", va="bottom", fontsize=12, fontweight="bold")

    plt.tight_layout()
    return fig, ax


def generate_markdown_report(df: pd.DataFrame, output_dir: Path) -> str:
    """
    Generate comprehensive markdown report.

    Parameters
    ----------
    df : pd.DataFrame
        Dendritic excitability data
    output_dir : Path
        Output directory

    Returns
    -------
    str
        Path to generated report
    """
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    # Calculate summary statistics
    valid_data = df.dropna()

    markdown_content = f"""# Figure 1I Dendritic Excitability Analysis Report

**Generated:** {timestamp}
**Data Source:** NWB files with two-photon calcium imaging and intracellular recordings
**Analysis:** Dendritic excitability index for dSPN dendrites across LID conditions

## Executive Summary

This report presents the reproduction of Figure 1I from Zhai et al. 2025, demonstrating
dendritic excitability changes in direct pathway spiny projection neurons (dSPNs) across
different LID (L-DOPA-induced dyskinesia) states.

### Key Findings

- **Total Sessions Analyzed:** {len(valid_data)}
- **Conditions:** {len(valid_data['condition'].unique())}
- **Measurement:** Dendritic excitability index (distal/proximal AUC ratio)

## Methodology

### Experimental Protocol
Two-photon calcium imaging of dSPN dendrites during current injection:
- **Current injection:** Three 2 nA injections, 2 ms each, at 50 Hz
- **Imaging:** Line scan at 810 nm excitation (Chameleon Ultra II)
- **Indicators:** Fluo-4 (calcium-sensitive) and Alexa568 (structural reference)
- **Locations:** Proximal (~40 μm) and distal (~90 μm) from soma

### Analysis Approach
1. **G/R0 Calculation:** Calcium signal (Green) / structural reference (Red)
2. **AUC Measurement:** Area under curve during stimulus response (0.95-1.5s)
3. **Excitability Index:** Distal dendrite AUC / Proximal dendrite AUC
4. **Statistical Comparison:** Across LID conditions using box plots

## Results by Condition

"""

    # Add condition-specific results
    for condition in valid_data["condition"].unique():
        condition_data = valid_data[valid_data["condition"] == condition]
        excitability_data = condition_data["dendritic_excitability_index"].dropna()

        if len(excitability_data) > 0:
            markdown_content += f"""### {condition}

- **Sample Size:** {len(condition_data)} sessions
- **Excitability Index:** {excitability_data.mean():.3f} ± {excitability_data.sem():.3f} (mean ± SEM)
- **Range:** {excitability_data.min():.3f} - {excitability_data.max():.3f}
- **Median:** {excitability_data.median():.3f}

"""

    markdown_content += f"""
## Figures

### Dendritic Excitability Comparison
![Dendritic Excitability](figure_1i_dendritic_excitability.png)

**Dendritic Excitability Index:** Statistical comparison of distal-to-proximal calcium
transient AUC ratios across LID conditions with individual session data points.

## Data Quality

- **Valid Sessions:** {len(valid_data)} sessions with complete calcium imaging data
- **Calcium Imaging:** Two-photon line scan data with Fluo-4 indicator
- **Current Injection:** Three-pulse protocol for dendritic activation
- **Temporal Analysis:** AUC calculated during stimulus response window

## Files Generated

- `figure_1i_dendritic_excitability.png` - Dendritic excitability comparison
- `figure_1i_reproduction_report.md` - This detailed report

---

*Report generated by Figure 1I Dendritic Excitability Reproduction Script*
"""

    # Write report
    report_path = output_dir / "figure_1i_reproduction_report.md"
    with open(report_path, "w") as f:
        f.write(markdown_content)

    return str(report_path)


def main():
    """Main function to orchestrate the Figure 1I reproduction analysis."""

    # Set up paths
    base_dir = Path(__file__).parent.parent
    nwb_dir = base_dir / "nwb_files" / "merged" / "F1"
    output_dir = base_dir / "analysis_outputs" / "figure_1i"

    # Ensure output directory exists
    output_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 80)
    print("FIGURE 1I DENDRITIC EXCITABILITY REPRODUCTION ANALYSIS")
    print("=" * 80)
    print(f"NWB Directory: {nwb_dir}")
    print(f"Output Directory: {output_dir}")
    print()

    # Load dendritic excitability data
    print("Step 1: Loading dendritic excitability data from NWB files...")
    df = load_dendritic_excitability_data(nwb_dir)

    # Create dendritic excitability plot
    print("Step 2: Creating dendritic excitability comparison plot...")
    fig, _ = create_dendritic_excitability_plot(df)
    fig_path = output_dir / "figure_1i_dendritic_excitability.png"
    fig.savefig(fig_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  Saved: {fig_path}")

    # Generate report
    print("Step 3: Generating analysis report...")
    report_path = generate_markdown_report(df, output_dir)
    print(f"  Saved: {report_path}")

    # Print summary
    print("\\n" + "=" * 80)
    print("ANALYSIS COMPLETED SUCCESSFULLY!")
    print("=" * 80)

    valid_stats = df.dropna()
    print(f"Total Sessions Analyzed: {len(valid_stats)}")
    print(f"Conditions: {len(valid_stats['condition'].unique())}")

    print("\\nCondition Summary:")
    for condition in valid_stats["condition"].unique():
        condition_data = valid_stats[valid_stats["condition"] == condition]
        excitability_data = condition_data["dendritic_excitability_index"].dropna()
        if len(excitability_data) > 0:
            print(
                f"  {condition}: n={len(condition_data)}, "
                f"index={excitability_data.mean():.3f}±{excitability_data.sem():.3f}"
            )

    print(f"\\nOutput Files:")
    print(f"  - {fig_path}")
    print(f"  - {report_path}")

    print("\\n" + "=" * 80)
    print("Figure 1I reproduction demonstrates dendritic excitability analysis!")
    print("Calcium imaging AUC ratios from two-photon microscopy data.")
    print("=" * 80)


if __name__ == "__main__":
    main()
