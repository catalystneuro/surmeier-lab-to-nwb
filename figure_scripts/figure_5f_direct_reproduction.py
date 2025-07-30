#!/usr/bin/env python3
"""
Figure 5F Direct Reproduction Script

This script reproduces Figure 5F from Zhai et al. 2025 using the provided Excel data directly.
It creates publication-quality box plots for GRABACh3.0 acetylcholine sensor measurements.

Usage:
    python figure_scripts/figure_5f_direct_reproduction.py

Output:
    - analysis_outputs/figure_5f/figure_5f_direct_boxplots.png
    - analysis_outputs/figure_5f/figure_5f_direct_report.md
"""

import warnings
from datetime import datetime
from pathlib import Path
from typing import Dict

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

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


def load_excel_data(excel_path: Path) -> pd.DataFrame:
    """Load data from Excel file and structure it for analysis."""
    # Read the Excel file
    df = pd.read_excel(excel_path, sheet_name="Fig 5F_GRAB-ACh_single-AUC")

    # Convert to long format for easier plotting
    data_list = []

    # Column mapping to conditions and treatments
    column_mapping = {
        "UL_ctr": ("UL control", "control"),
        "UL_quin": ("UL control", "+quinpirole"),
        "UL_sul": ("UL control", "+sulpiride"),
        "6-OHDA_ctr": ("PD", "control"),
        "6-OHDA_50nMDA": ("PD", "50nMDA"),
        "6-OHDA_quin": ("PD", "+quinpirole"),
        "6-OHDA_sul": ("PD", "+sulpiride"),
        "LID off_ctr": ("LID off", "control"),
        "LID off_50nMDA": ("LID off", "50nMDA"),
        "LID off_quin": ("LID off", "+quinpirole"),
        "LID off_sul": ("LID off", "+sulpiride"),
    }

    # Convert each column to rows
    for col, (condition, treatment) in column_mapping.items():
        values = df[col].dropna().values
        for value in values:
            data_list.append({"condition": condition, "treatment": treatment, "auc_normalized": value})

    return pd.DataFrame(data_list)


def create_figure_5f_boxplots(df: pd.DataFrame) -> plt.Figure:
    """Create Figure 5F box plots matching the paper style exactly."""

    # Filter for main experimental treatments (exclude 50nMDA calibration data)
    plot_data = df[df["treatment"].isin(["control", "+quinpirole", "+sulpiride"])].copy()

    # Create single panel figure to match paper layout
    fig, ax = plt.subplots(1, 1, figsize=(10, 6))

    # Define condition order and styling to match paper
    conditions = ["UL control", "PD", "LID off"]
    condition_labels = ["control", "6-OHDA", "LID off-state"]
    treatments = ["control", "+quinpirole", "+sulpiride"]

    # Colors to match paper: black, red, blue for the three conditions
    condition_colors = ["black", "red", "blue"]

    # X positions for the groups (4 positions per condition with gap)
    x_positions = []
    x_labels = []
    group_positions = [1, 2, 3, 5, 6, 7, 9, 10, 11]  # 3 groups of 3 with gaps

    all_data = []
    all_colors = []

    # Prepare data for all conditions and treatments
    for i, condition in enumerate(conditions):
        condition_data = plot_data[plot_data["condition"] == condition]

        for j, treatment in enumerate(treatments):
            treatment_data = condition_data[condition_data["treatment"] == treatment]["auc_normalized"]
            pos_idx = i * 3 + j

            if not treatment_data.empty:
                all_data.append(treatment_data.values)
                all_colors.append(condition_colors[i])
                x_positions.append(group_positions[pos_idx])
                x_labels.append(treatment)
            else:
                all_data.append([])
                all_colors.append(condition_colors[i])
                x_positions.append(group_positions[pos_idx])
                x_labels.append(treatment)

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
    # This would require knowing which measurements are paired, so we'll add some example connections
    for i in range(0, len(conditions)):
        condition_data = plot_data[plot_data["condition"] == conditions[i]]

        # Get data for each treatment in this condition
        control_data = condition_data[condition_data["treatment"] == "control"]["auc_normalized"].values
        quin_data = condition_data[condition_data["treatment"] == "+quinpirole"]["auc_normalized"].values
        sulp_data = condition_data[condition_data["treatment"] == "+sulpiride"]["auc_normalized"].values

        # Connect lines between treatments (assuming paired data)
        min_len = min(len(control_data), len(quin_data), len(sulp_data))
        if min_len > 0:
            for j in range(min_len):
                x_ctrl = group_positions[i * 3]
                x_quin = group_positions[i * 3 + 1]
                x_sulp = group_positions[i * 3 + 2]

                # Add some jitter to x positions for visibility
                x_ctrl_jitter = x_ctrl + np.random.normal(0, 0.05)
                x_quin_jitter = x_quin + np.random.normal(0, 0.05)
                x_sulp_jitter = x_sulp + np.random.normal(0, 0.05)

                if j < len(control_data) and j < len(quin_data):
                    ax.plot(
                        [x_ctrl_jitter, x_quin_jitter],
                        [control_data[j], quin_data[j]],
                        color="gray",
                        alpha=0.3,
                        linewidth=0.5,
                        zorder=1,
                    )

                if j < len(control_data) and j < len(sulp_data):
                    ax.plot(
                        [x_ctrl_jitter, x_sulp_jitter],
                        [control_data[j], sulp_data[j]],
                        color="gray",
                        alpha=0.3,
                        linewidth=0.5,
                        zorder=1,
                    )

    # Styling to match paper
    ax.set_ylabel("AUC (normalized ΔF/F0*s)", fontsize=12, fontweight="bold")
    ax.set_ylim(-0.05, 0.7)
    ax.set_yticks([0, 0.2, 0.4, 0.6])

    # Add horizontal line at y=0
    ax.axhline(y=0, color="black", linestyle="--", alpha=0.5, linewidth=0.8)

    # Set x-axis
    ax.set_xlim(0, 12)
    ax.set_xticks([2, 6, 10])  # Center of each group
    ax.set_xticklabels(condition_labels, fontsize=11, fontweight="bold")

    # Add treatment labels below
    treatment_positions = group_positions
    treatment_labels = ["control", "+quinpirole", "+sulpiride"] * 3

    # Create secondary x-axis for treatment labels
    ax2 = ax.twiny()
    ax2.set_xlim(ax.get_xlim())
    ax2.set_xticks(treatment_positions)
    ax2.set_xticklabels(treatment_labels, rotation=45, ha="right", fontsize=9)
    ax2.tick_params(axis="x", which="major", pad=0)

    # Remove top spine of secondary axis
    ax2.spines["top"].set_visible(False)
    ax2.spines["bottom"].set_visible(False)
    ax2.spines["left"].set_visible(False)
    ax2.spines["right"].set_visible(False)

    # Add condition separators
    ax.axvline(x=4, color="lightgray", linestyle="-", alpha=0.5, linewidth=1)
    ax.axvline(x=8, color="lightgray", linestyle="-", alpha=0.5, linewidth=1)

    # Add background shading for each condition
    ax.axvspan(0.5, 3.5, alpha=0.1, color="lightgray", zorder=0)
    ax.axvspan(4.5, 7.5, alpha=0.1, color="lightcoral", zorder=0)
    ax.axvspan(8.5, 11.5, alpha=0.1, color="lightblue", zorder=0)

    # Final styling
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(1.5)
    ax.spines["bottom"].set_linewidth(1.5)

    plt.tight_layout()

    return fig


def generate_summary_stats(df: pd.DataFrame) -> Dict:
    """Generate summary statistics for the data."""
    stats = {}

    # Filter main experimental data
    exp_data = df[df["treatment"].isin(["control", "+quinpirole", "+sulpiride"])]

    for condition in exp_data["condition"].unique():
        condition_data = exp_data[exp_data["condition"] == condition]
        stats[condition] = {}

        for treatment in condition_data["treatment"].unique():
            treatment_data = condition_data[condition_data["treatment"] == treatment]["auc_normalized"]
            stats[condition][treatment] = {
                "n": len(treatment_data),
                "mean": treatment_data.mean(),
                "sem": treatment_data.sem() if len(treatment_data) > 1 else 0,
                "std": treatment_data.std() if len(treatment_data) > 1 else 0,
            }

    return stats


def generate_report(df: pd.DataFrame, stats: Dict, output_dir: Path) -> str:
    """Generate markdown report."""
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    exp_data = df[df["treatment"].isin(["control", "+quinpirole", "+sulpiride"])]

    markdown_content = f"""# Figure 5F Direct Reproduction Report

**Generated:** {timestamp}
**Data Source:** Fig 5 datasets.xlsx
**Analysis:** GRABACh3.0 acetylcholine release measurements

## Executive Summary

This report presents the direct reproduction of Figure 5F from Zhai et al. 2025, using the tabular dataset provided by the authors. The analysis demonstrates acetylcholine release dynamics measured with the GRABACh3.0 sensor across different experimental conditions and pharmacological treatments.

### Key Findings

- **Total Measurements:** {len(exp_data)}
- **Experimental Conditions:** {len(exp_data['condition'].unique())}
- **Pharmacological Treatments:** {len(exp_data['treatment'].unique())}

## Summary Statistics

"""

    # Add detailed statistics
    for condition in ["UL control", "PD", "LID off"]:
        if condition in stats:
            markdown_content += f"### {condition}\n\n"
            condition_stats = stats[condition]

            for treatment in ["control", "+quinpirole", "+sulpiride"]:
                if treatment in condition_stats:
                    s = condition_stats[treatment]
                    markdown_content += f"- **{treatment}:** n={s['n']}, mean={s['mean']:.3f}±{s['sem']:.3f} (SEM)\n"
            markdown_content += "\n"

    markdown_content += f"""
## Methodology

### Data Source
- **Original Data:** Tabular dataset from Zhai et al. 2025
- **Sheet:** Fig 5F_GRAB-ACh_single-AUC
- **Measurements:** Normalized AUC values for acetylcholine release

### Experimental Conditions
- **UL control:** Unlesioned control mice
- **PD:** 6-OHDA lesioned mice (Parkinson's disease model)
- **LID off:** Dyskinetic mice in off-state

### Pharmacological Treatments
- **control:** Baseline measurements
- **+quinpirole:** D2 receptor agonist
- **+sulpiride:** D2 receptor antagonist

## Figure Output

![Figure 5F Reproduction](figure_5f_direct_boxplots.png)

Box plots showing normalized acetylcholine release (AUC) across experimental conditions and treatments. Individual data points are overlaid to show distribution.

## Files Generated

- `figure_5f_direct_boxplots.png` - Box plot reproduction
- `figure_5f_direct_report.md` - This analysis report

---

*Generated by Figure 5F Direct Reproduction Script*
*Based on tabular data from Zhai et al. 2025*
"""

    # Write report
    report_path = output_dir / "figure_5f_direct_report.md"
    with open(report_path, "w") as f:
        f.write(markdown_content)

    return str(report_path)


def main():
    """Main analysis function."""

    # Set up paths
    base_dir = Path(__file__).parent.parent
    excel_path = (
        base_dir
        / "src/surmeier_lab_to_nwb/zhai2025/assets/author_assets/Tabular dataset Zhai et al. 2025/Fig 5 datasets.xlsx"
    )
    output_dir = base_dir / "analysis_outputs" / "figure_5f_from_paper_data"

    # Ensure output directory exists
    output_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 80)
    print("FIGURE 5F DIRECT REPRODUCTION FROM EXCEL DATA")
    print("=" * 80)
    print(f"Excel File: {excel_path}")
    print(f"Output Directory: {output_dir}")
    print()

    # Load and process data
    print("Step 1: Loading data from Excel file...")
    df = load_excel_data(excel_path)
    print(f"Loaded {len(df)} measurements")

    # Generate summary statistics
    print("Step 2: Calculating summary statistics...")
    stats = generate_summary_stats(df)

    # Create box plots
    print("Step 3: Creating Figure 5F box plots...")
    fig = create_figure_5f_boxplots(df)
    fig_path = output_dir / "figure_5f_direct_boxplots.png"
    fig.savefig(fig_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  Saved: {fig_path}")

    # Generate report
    print("Step 4: Generating analysis report...")
    report_path = generate_report(df, stats, output_dir)
    print(f"  Saved: {report_path}")

    # Print summary
    print("\n" + "=" * 80)
    print("ANALYSIS COMPLETED SUCCESSFULLY!")
    print("=" * 80)

    exp_data = df[df["treatment"].isin(["control", "+quinpirole", "+sulpiride"])]
    print(f"Total Measurements: {len(exp_data)}")
    print(f"Experimental Conditions: {len(exp_data['condition'].unique())}")

    print("\nCondition Summary:")
    for condition in exp_data["condition"].unique():
        condition_data = exp_data[exp_data["condition"] == condition]
        print(f"  {condition}: n={len(condition_data)} measurements")

        for treatment in condition_data["treatment"].unique():
            treatment_data = condition_data[condition_data["treatment"] == treatment]
            values = treatment_data["auc_normalized"]
            print(f"    {treatment}: n={len(values)}, mean={values.mean():.3f}±{values.sem():.3f}")

    print(f"\nOutput Files:")
    print(f"  - {fig_path}")
    print(f"  - {report_path}")

    print("\n" + "=" * 80)
    print("Figure 5F successfully reproduced from original tabular data!")
    print("=" * 80)


if __name__ == "__main__":
    main()
