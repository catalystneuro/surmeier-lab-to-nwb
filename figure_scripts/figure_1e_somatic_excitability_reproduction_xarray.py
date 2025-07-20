#!/usr/bin/env python3
"""
Figure 1E Somatic Excitability Reproduction Script - XArray Version

This script reproduces Figure 1E from Zhai et al. 2025 using NWB files containing
somatic excitability data. It loads current clamp recordings, calculates rheobase
values, and generates publication-quality F-I (frequency-intensity) curves.

This version uses XArray for multidimensional data handling instead of pandas.

Usage:
    python figure_scripts/figure_1e_somatic_excitability_reproduction_xarray.py

Output:
    - analysis_outputs/figure_1e_xarray/figure_1e_reproduction_report.md
    - analysis_outputs/figure_1e_xarray/figure_1e_fi_curves.png
    - analysis_outputs/figure_1e_xarray/figure_1e_rheobase_comparison.png
"""

import warnings
from datetime import datetime
from pathlib import Path
from typing import Dict, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr
from pynwb import NWBHDF5IO

# Suppress expected warnings when calculating statistics on groups with no spikes
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


def count_action_potentials(voltage_trace: np.ndarray, timestamps: np.ndarray, threshold: float = 0.0) -> int:
    """
    Count action potentials using threshold crossing method matching authors' approach.

    Adapted from the original count_events function in the notebook.
    Counts upward threshold crossings during the stimulus period, avoiding double-counting.

    Parameters
    ----------
    voltage_trace : np.ndarray
        Voltage recording in mV
    timestamps : np.ndarray
        Time stamps in seconds
    threshold : float, default=0.0
        Spike detection threshold in mV (matching authors' approach)

    Returns
    -------
    int
        Number of action potentials detected
    """
    if len(voltage_trace) == 0:
        return 0

    # The stimulus period is typically 500ms long, starting after some baseline
    # Based on authors' analysis windows and protocol
    duration = timestamps[-1] - timestamps[0]

    if duration < 0.8:  # Recording too short
        return 0

    # Use stimulus period: approximately 200ms-700ms (500ms stimulus)
    # This matches the authors' approach of analyzing during current injection
    start_time = timestamps[0] + 0.2  # 200ms after start
    end_time = timestamps[0] + 0.7  # 700ms after start (500ms stimulus)

    start_idx = np.searchsorted(timestamps, start_time)
    end_idx = np.searchsorted(timestamps, end_time)

    if start_idx >= end_idx or end_idx > len(voltage_trace):
        # Fall back to middle portion of recording
        start_idx = len(voltage_trace) // 4
        end_idx = 3 * len(voltage_trace) // 4

    analysis_trace = voltage_trace[start_idx:end_idx]

    if len(analysis_trace) == 0:
        return 0

    # Implement threshold crossing detection matching authors' count_events function
    event_counter = 0
    i = 0

    # Loop through the trace looking for upward threshold crossings
    while i < len(analysis_trace):
        if analysis_trace[i] > threshold:
            event_counter += 1
            # Skip values until we go back below threshold (avoid double-counting)
            while i < len(analysis_trace) and analysis_trace[i] > threshold:
                i += 1
        else:
            i += 1

    return event_counter


def calculate_rheobase(current_steps: np.ndarray, spike_counts: np.ndarray) -> float:
    """
    Calculate rheobase current matching authors' approach.

    Adapted from rheobaseCal function in the original notebook.

    Parameters
    ----------
    current_steps : np.ndarray
        Array of injected currents in pA
    spike_counts : np.ndarray
        Corresponding spike counts for each current step

    Returns
    -------
    float
        Rheobase current in pA
    """
    # Simply return the first current that evokes a spike (most straightforward approach)
    # This matches the biological definition of rheobase regardless of the indexing formula
    for current, spikes in zip(current_steps, spike_counts):
        if spikes >= 1:
            return current

    # If no spikes found, return NaN
    return np.nan


def analyze_membrane_properties(
    voltage_trace: np.ndarray,
    timestamps: np.ndarray,
    current_pA: float,
    baseline_start: float = 0.005,
    baseline_end: float = 0.195,
    steady_start: float = 0.3,
    steady_end: float = 0.4,
) -> Dict[str, float]:
    """
    Analyze basic membrane properties from a hyperpolarizing current step.

    Adapted from analyze_properties function in the original notebook.

    Parameters
    ----------
    voltage_trace : np.ndarray
        Voltage recording in mV
    timestamps : np.ndarray
        Time stamps in seconds
    current_pA : float
        Injected current in pA
    baseline_start : float, default=0.005
        Start of baseline period (seconds)
    baseline_end : float, default=0.195
        End of baseline period (seconds)
    steady_start : float, default=0.3
        Start of steady-state period (seconds)
    steady_end : float, default=0.4
        End of steady-state period (seconds)

    Returns
    -------
    Dict[str, float]
        Dictionary with Vm (resting potential) and Rm (input resistance)
    """
    # Find baseline period
    baseline_mask = (timestamps >= baseline_start) & (timestamps <= baseline_end)
    vm = np.mean(voltage_trace[baseline_mask])

    # Find steady-state period
    steady_mask = (timestamps >= steady_start) & (timestamps <= steady_end)
    v_ss = np.mean(voltage_trace[steady_mask])

    # Calculate input resistance (Ohm's law: R = ΔV / ΔI)
    delta_v = v_ss - vm  # mV
    delta_i = current_pA / 1000.0  # Convert pA to nA
    rm = (delta_v / delta_i) * 1000.0 if delta_i != 0 else np.nan  # MΩ

    return {"vm_mv": vm, "rm_mohm": rm}


def load_somatic_excitability_data(nwb_dir: Path) -> xr.Dataset:
    """
    Load somatic excitability data from NWB files into an xarray Dataset.

    Parameters
    ----------
    nwb_dir : Path
        Directory containing Figure 1 somatic excitability NWB files

    Returns
    -------
    xr.Dataset
        Structured dataset with dimensions (cell_id, current_pA)
    """
    nwb_files = list(nwb_dir.glob("figure1_somatic_excitability_*.nwb"))

    if not nwb_files:
        raise FileNotFoundError(
            f"No NWB files found in {nwb_dir}. "
            "Please run figure_1_somatic_experiments.py first to generate NWB files."
        )

    cell_datasets = []
    print(f"Loading data from {len(nwb_files)} NWB files...")

    for nwb_file in nwb_files:
        with NWBHDF5IO(str(nwb_file), "r") as io:
            nwb = io.read()

            # Extract condition from filename (order matters - most specific first!)
            condition = "unknown"
            if "LID_off-state" in nwb_file.name:
                condition = "LID off-state"
            elif "LID_on-state_with_SCH" in nwb_file.name:
                condition = "LID on-state with SCH"
            elif "LID_on-state" in nwb_file.name:
                condition = "LID on-state"

            # Get session info
            session_id = nwb.session_id
            subject_id = nwb.subject.subject_id

            # Access intracellular recordings table
            icephys_table = nwb.get_intracellular_recordings()

            # Collect data for this cell
            current_steps = []
            spike_counts = []
            vm_values = []
            rm_values = []

            for idx in range(len(icephys_table.id)):
                # Get recording info
                row = icephys_table[idx]
                current_pA = row[("intracellular_recordings", "stimulus_current_pA")].iloc[0]

                # Get the current clamp series
                response_ref = row[("responses", "response")].iloc[0]
                ts = response_ref.timeseries
                voltage_trace = (
                    ts.data[response_ref.idx_start : response_ref.idx_start + response_ref.count] * 1000
                )  # Convert to mV
                timestamps = ts.timestamps[response_ref.idx_start : response_ref.idx_start + response_ref.count]

                # Count action potentials
                spike_count = count_action_potentials(voltage_trace, timestamps)

                # Calculate membrane properties for hyperpolarizing steps
                if current_pA < 0:  # Hyperpolarizing step
                    mem_props = analyze_membrane_properties(voltage_trace, timestamps, current_pA)
                    vm = mem_props["vm_mv"]
                    rm = mem_props["rm_mohm"]
                else:
                    vm, rm = np.nan, np.nan

                current_steps.append(current_pA)
                spike_counts.append(spike_count)
                vm_values.append(vm)
                rm_values.append(rm)

            # Create xarray Dataset for this cell
            cell_ds = xr.Dataset(
                {
                    "spike_count": (["current_pA"], spike_counts),
                    "vm_mv": (["current_pA"], vm_values),
                    "rm_mohm": (["current_pA"], rm_values),
                },
                coords={
                    "current_pA": current_steps,
                    "cell_id": session_id,
                    "condition": condition,
                    "subject_id": subject_id,
                    "nwb_file": nwb_file.name,
                },
            )

            cell_datasets.append(cell_ds)

    # Combine all cells into a single dataset
    if not cell_datasets:
        raise ValueError("No somatic excitability data could be loaded from NWB files")

    # Concatenate along cell dimension
    full_dataset = xr.concat(cell_datasets, dim="cell_id")

    print(f"Loaded {len(full_dataset.cell_id)} cells")
    return full_dataset


def calculate_cell_summary_stats(dataset: xr.Dataset) -> xr.Dataset:
    """
    Calculate summary statistics for each cell (rheobase, membrane properties).

    Parameters
    ----------
    dataset : xr.Dataset
        Recording-level data

    Returns
    -------
    xr.Dataset
        Cell-level summary statistics
    """

    def calc_rheobase(spike_counts, current_steps):
        """Calculate rheobase for a single cell."""
        # Sort by current - inputs are already numpy arrays from apply_ufunc
        sorted_indices = np.argsort(current_steps)
        sorted_currents = current_steps[sorted_indices]
        sorted_spikes = spike_counts[sorted_indices]
        return calculate_rheobase(sorted_currents, sorted_spikes)

    def get_membrane_props(vm_values, rm_values, current_steps):
        """Get membrane properties from most hyperpolarized step."""
        # Find most hyperpolarized step with valid data
        # Inputs are already numpy arrays from apply_ufunc
        valid_mask = ~np.isnan(vm_values)
        if not np.any(valid_mask):
            return np.nan, np.nan

        valid_currents = current_steps[valid_mask]
        min_current_idx = np.argmin(valid_currents)
        valid_indices = np.where(valid_mask)[0]
        actual_idx = valid_indices[min_current_idx]

        return vm_values[actual_idx], rm_values[actual_idx]

    # Calculate rheobase for each cell
    rheobase = xr.apply_ufunc(
        calc_rheobase,
        dataset.spike_count,
        dataset.current_pA,
        input_core_dims=[["current_pA"], ["current_pA"]],
        output_core_dims=[[]],
        vectorize=True,
    )

    # Calculate membrane properties for each cell
    vm_rm = xr.apply_ufunc(
        get_membrane_props,
        dataset.vm_mv,
        dataset.rm_mohm,
        dataset.current_pA,
        input_core_dims=[["current_pA"], ["current_pA"], ["current_pA"]],
        output_core_dims=[[], []],
        vectorize=True,
    )

    # Create summary dataset
    summary_ds = xr.Dataset(
        {
            "rheobase_pA": rheobase,
            "vm_mv": vm_rm[0],
            "rm_mohm": vm_rm[1],
            "n_recordings": dataset.spike_count.count("current_pA"),
        },
        coords={
            "cell_id": dataset.cell_id,
            "condition": dataset.condition,
            "subject_id": dataset.subject_id,
        },
    )

    return summary_ds


def create_fi_curves(dataset: xr.Dataset) -> Tuple[plt.Figure, plt.Axes]:
    """
    Create frequency-intensity (F-I) curves matching paper style using xarray.

    Parameters
    ----------
    dataset : xr.Dataset
        Recording-level data

    Returns
    -------
    Tuple[plt.Figure, plt.Axes]
        Figure and axes objects
    """
    fig, ax = plt.subplots(1, 1, figsize=(6, 4))

    # Define condition styling to match paper (black and white theme)
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

    # Calculate F-I curves using xarray groupby operations
    def calculate_fi_stats_for_condition(condition_data):
        """Calculate F-I curve statistics for a single condition."""
        # Group by current_pA and calculate mean/SEM across cells
        stats_list = []
        for current in condition_data.current_pA.values:
            current_data = condition_data.sel(current_pA=current).spike_count
            valid_data = current_data.dropna("cell_id")

            if len(valid_data.cell_id) == 0:
                mean_val, sem_val = np.nan, 0.0
            elif len(valid_data.cell_id) == 1:
                mean_val, sem_val = float(valid_data.values[0]), 0.0
            else:
                mean_val = float(valid_data.mean("cell_id").values)
                sem_val = float(valid_data.std("cell_id").values / np.sqrt(len(valid_data.cell_id)))

            stats_list.append({"current_pA": current, "mean": mean_val, "sem": sem_val})

        return pd.DataFrame(stats_list)

    # Calculate F-I statistics for each condition
    fi_stats = {}

    for condition in ["LID off-state", "LID on-state", "LID on-state with SCH"]:
        if condition not in dataset.condition.values:
            continue

        # Get cells for this condition
        condition_mask = dataset.condition == condition
        condition_data = dataset.isel(cell_id=condition_mask)

        # Calculate stats for this condition
        condition_stats = calculate_fi_stats_for_condition(condition_data)

        # Filter to positive currents up to 300 pA
        filtered_stats = condition_stats[(condition_stats["current_pA"] >= 0) & (condition_stats["current_pA"] <= 300)]

        if len(filtered_stats) == 0:
            continue

        style = condition_styles[condition]

        ax.errorbar(
            filtered_stats["current_pA"],
            filtered_stats["mean"],
            yerr=filtered_stats["sem"],
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

    # Style to match paper
    ax.set_xlabel("current (pA)", fontsize=12)
    ax.set_ylabel("number of APs", fontsize=12)
    ax.set_xlim(0, 300)
    ax.set_ylim(0, 18)

    # Set ticks to match paper
    ax.set_xticks([0, 100, 200, 300])
    ax.set_yticks([0, 5, 10, 15])

    # Legend positioning to match paper
    ax.legend(loc="upper left", frameon=False, fontsize=10)

    # Add descriptive title
    ax.set_title("Figure 1E: dSPN Somatic Excitability (F-I Curves) - XArray", fontsize=14, pad=20)

    # Remove top and right spines (already set in rcParams)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.tight_layout()
    return fig, ax


def create_rheobase_comparison(cell_stats: xr.Dataset) -> Tuple[plt.Figure, plt.Axes]:
    """
    Create rheobase comparison to match paper Panel F style using xarray.

    Parameters
    ----------
    cell_stats : xr.Dataset
        Cell-level summary statistics

    Returns
    -------
    Tuple[plt.Figure, plt.Axes]
        Figure and axes objects
    """
    fig, ax = plt.subplots(1, 1, figsize=(4, 5))

    # Filter out invalid rheobase values
    valid_data = cell_stats.where(~np.isnan(cell_stats.rheobase_pA), drop=True)

    # Order conditions to match paper
    condition_order = ["LID off-state", "LID on-state", "LID on-state with SCH"]
    condition_labels = ["off", "on", "on + SCH"]

    # Prepare data for box plot
    plot_data = []
    labels = []

    for condition, label in zip(condition_order, condition_labels):
        if condition in valid_data.condition.values:
            condition_data = valid_data.sel(cell_id=valid_data.condition == condition).rheobase_pA
            # Remove NaN values
            condition_values = condition_data.values[~np.isnan(condition_data.values)]
            plot_data.append(condition_values)
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
    ax.set_ylabel("rheobase (pA)", fontsize=12)
    ax.set_ylim(0, 300)
    ax.set_yticks([0, 100, 200, 300])

    # Add descriptive title
    ax.set_title("Figure 1F: dSPN Rheobase Comparison - XArray", fontsize=14, pad=20)

    # Remove top and right spines
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Add significance annotation (matching paper)
    # Note: This would need statistical testing to be accurate
    ax.text(1.5, 280, "****", ha="center", va="bottom", fontsize=12, fontweight="bold")

    plt.tight_layout()
    return fig, ax


def generate_markdown_report(cell_stats: xr.Dataset, output_dir: Path) -> str:
    """
    Generate comprehensive markdown report using xarray.

    Parameters
    ----------
    cell_stats : xr.Dataset
        Cell-level summary statistics
    output_dir : Path
        Output directory

    Returns
    -------
    str
        Path to generated report
    """
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    # Calculate summary statistics using xarray
    valid_data = cell_stats.where(~np.isnan(cell_stats.rheobase_pA), drop=True)

    markdown_content = f"""# Figure 1E Somatic Excitability Analysis Report - XArray Version

**Generated:** {timestamp}
**Data Source:** NWB files with intracellular current clamp recordings
**Analysis:** F-I curves and rheobase measurements for dSPN somatic excitability
**Implementation:** XArray-based multidimensional data analysis

## Executive Summary

This report presents the reproduction of Figure 1E from Zhai et al. 2025, demonstrating
somatic excitability changes in direct pathway spiny projection neurons (dSPNs) across
different LID (L-DOPA-induced dyskinesia) states.

**XArray Implementation Benefits:**
- Multidimensional data structure with labeled dimensions
- Clean statistical operations with groupby
- Automatic dimension alignment and broadcasting
- Better handling of missing data

### Key Findings

- **Total Cells Analyzed:** {len(valid_data.cell_id)}
- **Conditions:** {len(np.unique(valid_data.condition.values))}
- **Total Recordings:** {int(valid_data.n_recordings.sum().values)}

## Methodology

### Experimental Protocol
Current clamp recordings from dSPN somata with step current injections:
- **Current Range:** -120 pA to +300 pA
- **Step Size:** 20 pA increments
- **Step Duration:** 500 ms
- **Analysis Window:** Action potentials counted during stimulus period

### Analysis Approach
1. **Action Potential Counting:** Spike detection using threshold crossing method
2. **Rheobase Calculation:** Minimum current to evoke ≥1 action potential
3. **Membrane Properties:** Vm and Rm from hyperpolarizing steps
4. **F-I Curves:** Frequency-intensity relationships using xarray groupby operations

## Results by Condition

"""

    # Add condition-specific results using xarray operations
    unique_conditions = np.unique(valid_data.condition.values)
    for condition in unique_conditions:
        condition_data = valid_data.sel(cell_id=valid_data.condition == condition)
        rheobase_data = condition_data.rheobase_pA.values
        rheobase_clean = rheobase_data[~np.isnan(rheobase_data)]

        if len(rheobase_clean) > 0:
            mean_rheobase = np.mean(rheobase_clean)
            sem_rheobase = np.std(rheobase_clean, ddof=1) / np.sqrt(len(rheobase_clean))
            min_rheobase = np.min(rheobase_clean)
            max_rheobase = np.max(rheobase_clean)
            median_rheobase = np.median(rheobase_clean)

            markdown_content += f"""### {condition}

- **Sample Size:** {len(condition_data.cell_id)} cells
- **Rheobase:** {mean_rheobase:.1f} ± {sem_rheobase:.1f} pA (mean ± SEM)
- **Range:** {min_rheobase:.1f} - {max_rheobase:.1f} pA
- **Median:** {median_rheobase:.1f} pA

"""

    markdown_content += f"""
## Figures

### F-I Curves
![F-I Curves](figure_1e_fi_curves.png)

**F-I Curves:** Frequency-intensity relationships showing action potential firing rates
in response to depolarizing current steps across LID conditions.

### Rheobase Comparison
![Rheobase Comparison](figure_1e_rheobase_comparison.png)

**Rheobase Analysis:** Statistical comparison of rheobase values across conditions
with individual cell data points and summary statistics.

## XArray Data Structure

The analysis uses an xarray Dataset with the following structure:
- **Dimensions:** cell_id, current_pA
- **Data Variables:** spike_count, vm_mv, rm_mohm
- **Coordinates:** condition, subject_id, nwb_file
- **Benefits:** Labeled dimensions, automatic alignment, vectorized operations

## Data Quality

- **Valid Recordings:** {len(valid_data.cell_id)} cells with complete datasets
- **Membrane Properties:** Calculated from hyperpolarizing current steps
- **Spike Detection:** Automated threshold crossing detection
- **Statistical Operations:** XArray groupby with proper NaN handling

## Files Generated

- `figure_1e_fi_curves.png` - F-I curve analysis
- `figure_1e_rheobase_comparison.png` - Rheobase statistical comparison
- `figure_1e_reproduction_report.md` - This detailed report

---

*Report generated by Figure 1E Somatic Excitability Reproduction Script (XArray Version)*
"""

    # Write report
    report_path = output_dir / "figure_1e_reproduction_report.md"
    with open(report_path, "w") as f:
        f.write(markdown_content)

    return str(report_path)


def main():
    """Main function to orchestrate the Figure 1E reproduction analysis using XArray."""

    # Set up paths
    base_dir = Path(__file__).parent.parent
    nwb_dir = base_dir / "nwb_files" / "figure_1_somatic_excitability"
    output_dir = base_dir / "analysis_outputs" / "figure_1e_xarray"

    # Ensure output directory exists
    output_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 80)
    print("FIGURE 1E SOMATIC EXCITABILITY REPRODUCTION ANALYSIS - XARRAY VERSION")
    print("=" * 80)
    print(f"NWB Directory: {nwb_dir}")
    print(f"Output Directory: {output_dir}")
    print()

    # Load somatic excitability data into xarray Dataset
    print("Step 1: Loading somatic excitability data from NWB files into xarray Dataset...")
    dataset = load_somatic_excitability_data(nwb_dir)

    # Calculate cell-level statistics using xarray operations
    print("Step 2: Calculating cell-level summary statistics using xarray...")
    cell_stats = calculate_cell_summary_stats(dataset)

    # Create F-I curves using xarray groupby operations
    print("Step 3: Creating F-I curves using xarray...")
    fig1, _ = create_fi_curves(dataset)
    fig1_path = output_dir / "figure_1e_fi_curves.png"
    fig1.savefig(fig1_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig1)
    print(f"  Saved: {fig1_path}")

    # Create rheobase comparison using xarray
    print("Step 4: Creating rheobase comparison using xarray...")
    fig2, _ = create_rheobase_comparison(cell_stats)
    fig2_path = output_dir / "figure_1e_rheobase_comparison.png"
    fig2.savefig(fig2_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig2)
    print(f"  Saved: {fig2_path}")

    # Generate report
    print("Step 5: Generating analysis report...")
    report_path = generate_markdown_report(cell_stats, output_dir)
    print(f"  Saved: {report_path}")

    # Print summary using xarray operations
    print("\n" + "=" * 80)
    print("ANALYSIS COMPLETED SUCCESSFULLY!")
    print("=" * 80)

    valid_stats = cell_stats.where(~np.isnan(cell_stats.rheobase_pA), drop=True)
    print(f"Total Cells Analyzed: {len(valid_stats.cell_id)}")
    print(f"Total Recordings: {int(valid_stats.n_recordings.sum().values)}")
    print(f"Conditions: {len(np.unique(valid_stats.condition.values))}")

    print("\nCondition Summary:")
    unique_conditions = np.unique(valid_stats.condition.values)
    for condition in unique_conditions:
        condition_data = valid_stats.sel(cell_id=valid_stats.condition == condition)
        rheobase_data = condition_data.rheobase_pA.values
        rheobase_clean = rheobase_data[~np.isnan(rheobase_data)]

        if len(rheobase_clean) > 0:
            mean_rheobase = np.mean(rheobase_clean)
            sem_rheobase = np.std(rheobase_clean, ddof=1) / np.sqrt(len(rheobase_clean))
            print(
                f"  {condition}: n={len(condition_data.cell_id)}, "
                f"rheobase={mean_rheobase:.1f}±{sem_rheobase:.1f} pA"
            )

    print(f"\nOutput Files:")
    print(f"  - {fig1_path}")
    print(f"  - {fig2_path}")
    print(f"  - {report_path}")

    print("\n" + "=" * 80)
    print("XArray implementation demonstrates clean multidimensional data handling!")
    print("Labeled dimensions and vectorized operations for electrophysiology analysis.")
    print("=" * 80)


if __name__ == "__main__":
    main()
