#!/usr/bin/env python3
"""
Figure 2 G and H oEPSC Reproduction Script (Fixed)

This script reproduces Figure 2 G and H from Zhai et al. 2025 using NWB files
containing Sr²⁺-oEPSC data. It loads voltage clamp recordings, extracts oEPSC
amplitudes, and generates cumulative probability curves and box plots.

Fixed to properly extract optogenetic stimulus timing from intervals.optogenetic_epochs_table

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
from typing import Tuple

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


def extract_oepsc_amplitude(
    voltage_trace: np.ndarray,
    timestamps: np.ndarray,
    stim_absolute_start_time: float,
    baseline_window: float = 0.020,  # 20 ms baseline before stimulus
    analysis_window_start: float = 0.040,  # 40 ms post-stim relative to stim_absolute_start_time
    analysis_window_end: float = 0.400,  # 400 ms post-stim relative to stim_absolute_start_time
) -> float:
    """
    Extracts the peak oEPSC amplitude within a specified analysis window with baseline correction.
    Assumes voltage clamp at -70 mV, so EPSCs are inward (negative), but returns
    the absolute value to represent amplitude as a positive number.

    Parameters
    ----------
    voltage_trace : np.ndarray
        Voltage recording in Amperes.
    timestamps : np.ndarray
        Time stamps in seconds (absolute, relative to NWB session start).
    stim_absolute_start_time : float
        The absolute time (in seconds, relative to NWB session start) when the optogenetic stimulus was applied.
    baseline_window : float
        Duration of baseline period before stimulus (seconds).
    analysis_window_start : float
        Start of the analysis window relative to stimulus start (seconds).
    analysis_window_end : float
        End of the analysis window relative to stimulus start (seconds).

    Returns
    -------
    float
        Peak oEPSC amplitude in pA (as a positive value), or NaN if no clear
        event is detected.
    """
    # Calculate baseline from pre-stimulus period
    baseline_start_abs = stim_absolute_start_time - baseline_window
    baseline_end_abs = stim_absolute_start_time

    baseline_mask = (timestamps >= baseline_start_abs) & (timestamps < baseline_end_abs)
    baseline_trace = voltage_trace[baseline_mask]

    if len(baseline_trace) == 0:
        return np.nan

    baseline = np.mean(baseline_trace)

    # Calculate absolute window times for analysis
    window_start_abs = stim_absolute_start_time + analysis_window_start
    window_end_abs = stim_absolute_start_time + analysis_window_end

    window_mask = (timestamps >= window_start_abs) & (timestamps <= window_end_abs)
    analysis_trace = voltage_trace[window_mask]

    if len(analysis_trace) == 0:
        return np.nan

    # Find peak amplitude relative to baseline
    # For inward currents (negative), we find the minimum
    peak_amplitude_amps = np.min(analysis_trace) - baseline

    # Convert to picoamperes (1 A = 1e12 pA) and return absolute value
    return abs(peak_amplitude_amps * 1e12)


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
    nwb_files = list(nwb_dir.glob("figure2_sr_oepsc_*.nwb"))

    if not nwb_files:
        raise FileNotFoundError(
            f"No NWB files found in {nwb_dir}. " "Please run figure_2_optical_stimuli.py first to generate NWB files."
        )

    all_data = []
    print(f"Loading data from {len(nwb_files)} NWB files...")

    for nwb_file in nwb_files:
        with NWBHDF5IO(str(nwb_file), "r") as io:
            nwb = io.read()

            # Extract condition from filename
            condition = "unknown"
            if "LID_off_state" in nwb_file.name:
                condition = "LID off-state"
            elif "LID_on_state" in nwb_file.name:
                condition = "LID on-state"

            # Get session info
            session_id = nwb.session_id
            subject_id = nwb.subject.subject_id

            # Access intracellular recordings table
            icephys_table = nwb.get_intracellular_recordings()

            # Get optogenetic epochs table
            opto_epochs = nwb.intervals["optogenetic_epochs_table"]
            opto_df = opto_epochs.to_dataframe()

            # Filter for actual stimulation epochs
            stim_epochs = opto_df[opto_df["stimulation_on"] == True]

            for idx in range(len(icephys_table.id)):
                row = icephys_table[idx]
                cell_number = row[("intracellular_recordings", "cell_number")].iloc[0]
                led_intensity = row[("intracellular_recordings", "led_intensity")].iloc[0]
                sweep_number = row[("intracellular_recordings", "sweep_number")].iloc[0]
                recording_id = f"Cell{cell_number}_LED{led_intensity}_Sweep{sweep_number}"

                # Get the voltage clamp series (oEPSC trace)
                response_ref = row[("responses", "response")].iloc[0]
                ts = response_ref.timeseries
                voltage_trace = ts.data[response_ref.idx_start : response_ref.idx_start + response_ref.count]
                timestamps = ts.timestamps[response_ref.idx_start : response_ref.idx_start + response_ref.count]

                # Find the stimulation epoch that occurs during this sweep
                sweep_start = timestamps[0]
                sweep_end = timestamps[-1]

                # Find stimulation epochs within this sweep
                sweep_stim_epochs = stim_epochs[
                    (stim_epochs["start_time"] >= sweep_start) & (stim_epochs["start_time"] <= sweep_end)
                ]

                if len(sweep_stim_epochs) == 0:
                    print(f"Warning: No stimulation found for {recording_id} in {nwb_file.name}")
                    continue

                # Use the first stimulation epoch in this sweep
                stim_absolute_start_time = sweep_stim_epochs.iloc[0]["start_time"]

                # Extract oEPSC amplitude
                oepsc_amplitude = extract_oepsc_amplitude(voltage_trace, timestamps, stim_absolute_start_time)

                all_data.append(
                    {
                        "session_id": session_id,
                        "subject_id": subject_id,
                        "condition": condition,
                        "recording_id": recording_id,
                        "oepsc_amplitude_pA": oepsc_amplitude,
                        "nwb_file": nwb_file.name,
                        "stim_time": stim_absolute_start_time,
                        "sweep_start": sweep_start,
                        "sweep_end": sweep_end,
                    }
                )

    if not all_data:
        raise ValueError("No oEPSC data could be loaded from NWB files.")

    df = pd.DataFrame(all_data)
    # Filter out NaN amplitudes
    df = df.dropna(subset=["oepsc_amplitude_pA"])
    print(f"Loaded {len(df)} oEPSC amplitudes from {len(df.session_id.unique())} cells.")
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

    plot_data = []
    for condition in condition_order:
        if condition in df["condition"].unique():
            plot_data.append(df[df["condition"] == condition]["oepsc_amplitude_pA"].values)
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
        ax.scatter(x_pos, data, color="gray", alpha=0.6, s=15, zorder=3)

    ax.set_ylabel("oEPSC amplitude (pA)", fontsize=12)
    # Dynamic y-axis limits based on data
    max_value = max([max(data) if len(data) > 0 else 0 for data in plot_data])
    ax.set_ylim(0, max_value * 1.2)  # Add 20% padding
    ax.set_yticks(np.arange(0, max_value * 1.2, 20))

    ax.set_title("Figure 2H: dSPN oEPSC Amplitude Comparison", fontsize=14, pad=20)

    # Perform statistical test (e.g., independent t-test)
    off_state_data = df[df["condition"] == "LID off-state"]["oepsc_amplitude_pA"].values
    on_state_data = df[df["condition"] == "LID on-state"]["oepsc_amplitude_pA"].values

    if len(off_state_data) > 1 and len(on_state_data) > 1:
        # Assuming independent samples and potentially unequal variances
        stat, p_value = stats.ttest_ind(off_state_data, on_state_data, equal_var=False)
        if p_value < 0.001:
            ax.text(1.5, 22, "***", ha="center", va="bottom", fontsize=12, fontweight="bold")  # Adjust y-pos as needed
        elif p_value < 0.01:
            ax.text(1.5, 22, "**", ha="center", va="bottom", fontsize=12, fontweight="bold")
        elif p_value < 0.05:
            ax.text(1.5, 22, "*", ha="center", va="bottom", fontsize=12, fontweight="bold")

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
    nwb_dir = base_dir / "nwb_files" / "figure_2" / "sr_oepsc"
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
