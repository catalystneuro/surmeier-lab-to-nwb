# Surmeier Lab to NWB: Zhai 2025 Conversion

NWB conversion scripts for Surmeier lab data from the **Zhai et al. 2025** paper on levodopa-induced dyskinesia research. This project converts multimodal neuroscience data to the [Neurodata Without Borders](https://nwb-overview.readthedocs.io/) format for standardized data sharing and analysis.

## Project Overview

This repository contains conversion scripts for data from the **Zhai et al. 2025** paper studying levodopa-induced dyskinesia (LID) in a mouse model of Parkinson's disease. The study examines cellular and molecular mechanisms underlying dyskinesia through comprehensive analysis of direct and indirect pathway spiny projection neurons (dSPNs and iSPNs) in the dorsolateral striatum.

### Paper Details
- **Title**: M1 muscarinic receptor-mediated dendritic excitability underlies striatal dysfunction in levodopa-induced dyskinesia
- **Preprint**: [bioRxiv 2025.01.02.631090](https://www.biorxiv.org/content/10.1101/2025.01.02.631090v1.full)
- **Key findings**: M1 muscarinic receptors drive dendritic hyperexcitability in indirect pathway neurons, contributing to dyskinetic behaviors

### Data Types Converted

- **Electrophysiology**: Patch-clamp recordings (MultiClamp 700B)
- **Two-photon imaging**: Calcium imaging, acetylcholine biosensor (GRABACh3.0), spine density analysis
- **Confocal microscopy**: High-resolution spine density measurements
- **Optogenetics**: Light-evoked postsynaptic currents (oEPSCs)
- **Behavioral data**: AIM scores, contralateral rotations, video recordings
- **Pharmacology**: Drug treatments and receptor manipulations

### Experimental Models
- **6-OHDA lesion model**: Parkinson's disease simulation
- **Chronic levodopa treatment**: Dyskinesia induction
- **Cell-type specific analysis**: dSPNs vs iSPNs
- **Genetic manipulations**: CDGI knockout, M1R antagonists

## Installation

### Installation from GitHub

This package is currently only available from GitHub (not yet released on PyPI). This project requires Python 3.12 or higher.

#### Option 1: Using uv (Recommended)
The project uses `uv` for fast and reliable dependency management. Install `uv` first ([installation instructions](https://docs.astral.sh/uv/getting-started/installation/)), then:

```bash
git clone https://github.com/catalystneuro/surmeier-lab-to-nwb
cd surmeier-lab-to-nwb
uv sync
```

This will create a virtual environment and install all dependencies automatically.

#### Option 2: Using conda
If you prefer conda for environment management ([installation instructions](https://docs.conda.io/en/latest/miniconda.html)):

```bash
git clone https://github.com/catalystneuro/surmeier-lab-to-nwb
cd surmeier-lab-to-nwb
conda create -n surmeier-lab-to-nwb python=3.12
conda activate surmeier-lab-to-nwb
pip install --editable .
```

#### Option 3: Using pip with virtual environment
Alternatively, you can use standard Python tools:

```bash
git clone https://github.com/catalystneuro/surmeier-lab-to-nwb
cd surmeier-lab-to-nwb
python -m venv .venv
source .venv/bin/activate  # On Windows: .venv\Scripts\activate
pip install --editable .
```

### Running the full conversion

All conversion tasks are wired up as `poethepoet` tasks defined in `pyproject.toml`. The canonical "convert everything for DANDI" command is:

```bash
uv run poe convert_all_full
```

This runs the seven-step pipeline in order: per-mouse-day patch-clamp merge, then spine density, confocal, AIM behavior, videos, biosensor, and optical stimulation conversions. Total runtime is roughly 30-45 minutes; outputs land in `nwb_files/`.

Other useful entry points:

```bash
# Patch-clamp only (the merged per-mouse-day files: Figures 1, 3, 6, 7, 8)
uv run poe merge_patch_clamp

# Single modality
uv run poe convert_spine_full
uv run poe convert_biosensor_full
uv run poe convert_behavior_full
uv run poe convert_videos_full
uv run poe convert_optical_full
uv run poe convert_confocal_full

# Quick stub-test variants (small subset of data, for development)
uv run poe convert_all
```

### Running a specific conversion script

Each figure-specific script under `src/surmeier_lab_to_nwb/zhai2025/conversion_scripts/` can also be invoked directly. The poe tasks above are thin wrappers around these; they exist as a fallback for fine-grained re-runs.

```bash
# Convert Figure 5 acetylcholine biosensor data
python src/surmeier_lab_to_nwb/zhai2025/conversion_scripts/acetylcholine_biosensor/figure_5_acetylcholine_biosensor.py --stub-test=False
```

### Per-mouse-day file architecture (patch-clamp)

For Figures 1, 3, 6, 7, and 8 (somatic + dendritic excitability), the conversion produces **one NWB file per mouse-day**. Each file contains all cells patched in that mouse on that day, with:

- One `CurrentClampSeries` per (cell, modality) holding all sweeps concatenated end-to-end (Pattern B). Individual sweeps are sliced via `idx_start + count` columns on the `IntracellularRecordingsTable`.
- The full 5-level icephys hierarchical chain (`IntracellularRecordings` -> `SimultaneousRecordings` -> `SequentialRecordings` -> `RepetitionsTable` -> `ExperimentalConditionsTable`).
- Paired line-scan kymographs for dendritic sweeps, consolidated per (cell, channel, dendrite location).
- A denormalized `RecordingsIndexTable` in `processing/recordings_index/` that pre-joins the chain for fast pandas-style queries.

The `session_id` follows a 5-token format: `<cellType>++<state>++<pharm>++<genotype>++<date>` (e.g., `dSPN++OffState++none++WT++20170202`).

### Reading from the dandiset

The notebooks under `notebooks/` stream files directly from DANDI without downloading. Start with `how_to_use_this_dataset.ipynb` for a tour of the new per-mouse-day patch-clamp files, the `RecordingsIndexTable`, line-scan access, AIM `DynamicTable`, optogenetic Sr-oEPSC, and acetylcholine biosensor data. The per-figure reproduction notebooks (`figure_1E_dspn_somatic_excitability.ipynb`, `figure_2GH_oepsc_analysis.ipynb`, `figure_5F_acetylcholine_biosensor.ipynb`) demonstrate downstream analysis on the published figures.

**Note:** All installation methods above install the repository in [editable mode](https://pip.pypa.io/en/stable/cli/pip_install/#editable-installs), allowing you to modify the source code if needed.


## Repository structure
Each conversion is organized in a directory of its own in the `src` directory:

    surmeier-lab-to-nwb/
    ├── LICENSE
    ├── README.md
    ├── pyproject.toml
    ├── notebooks/                       # Streaming-based reproduction notebooks
    │   ├── how_to_use_this_dataset.ipynb
    │   ├── figure_1E_dspn_somatic_excitability.ipynb
    │   ├── figure_2GH_oepsc_analysis.ipynb
    │   └── figure_5F_acetylcholine_biosensor.ipynb
    ├── scripts/                         # Conversion runners and validation
    │   ├── full_per_mouse_day_merge.py  # Production patch-clamp merger
    │   ├── prototype_per_mouse_day_merger.py  # Merger core (library)
    │   ├── validate_all_fi_curves.py    # Paper-panel F-I validation
    │   ├── validate_figure_7c.py
    │   ├── validate_against_dandi_streaming.py  # Cross-validation vs DANDI
    │   └── wipe_dandiset.py             # Bulk-delete DANDI assets (re-upload)
    └── src/
        └── surmeier_lab_to_nwb/
            └── zhai2025/
                ├── conversion_notes_folder/   # Detailed conversion documentation
                ├── conversion_scripts/        # Per-figure conversion scripts
                │   ├── acetylcholine_biosensor/
                │   ├── aim_behavior/
                │   ├── confocal_spine_density/
                │   ├── dendritic_excitability/
                │   ├── optical_stimulation/
                │   ├── somatic_excitability/
                │   ├── spine_density/
                │   └── videos/
                ├── interfaces/                # Custom NeuroConv interfaces (Bruker BOT, line-scan, behavior)
                ├── mouse_day_merger.py        # Per-mouse-day bundle builder (raw folders → bundles)
                ├── subject_registry.py        # Data Connections spreadsheet parser
                └── utils/                     # Utility functions

## Figure-Specific Conversions

The `zhai2025` conversion includes specialized scripts for each figure in the paper. **For the dandiset deposit, the patch-clamp scripts under Dendritic Excitability and Somatic Excitability are superseded by the per-mouse-day merger** (`uv run poe merge_patch_clamp`). The per-figure scripts below are still functional and used as the lower-level converters that feed the merger; they remain the right entry point for ad-hoc per-condition output.

### Dendritic Excitability (Figures 1, 3, 6, 7) — merged into per-mouse-day files
- Two-photon calcium imaging with patch-clamp electrophysiology
- Back-propagating action potential measurements
- Cell-type specific analysis (dSPNs vs iSPNs)
- **Scripts**:
  - [Figure 1](src/surmeier_lab_to_nwb/zhai2025/conversion_scripts/dendritic_excitability/figure_1_dendritic_excitability.py)
  - [Figure 3](src/surmeier_lab_to_nwb/zhai2025/conversion_scripts/dendritic_excitability/figure_3_dendritic_excitability.py)
  - [Figure 6](src/surmeier_lab_to_nwb/zhai2025/conversion_scripts/dendritic_excitability/figure_6_dendritic_excitability.py)
  - [Figure 7](src/surmeier_lab_to_nwb/zhai2025/conversion_scripts/dendritic_excitability/figure_7_dendritic_excitability.py)
  - [Figure 7 OxoM](src/surmeier_lab_to_nwb/zhai2025/conversion_scripts/dendritic_excitability/figure_7_oxoM_dendritic_excitability.py)

### Somatic Excitability (Figures 1, 3, 6, 7, 8) — merged into per-mouse-day files
- Patch-clamp electrophysiology recordings
- Current injection protocols and F-I curves
- **Scripts**:
  - [Figure 1](src/surmeier_lab_to_nwb/zhai2025/conversion_scripts/somatic_excitability/figure_1_somatic_excitability.py)
  - [Figure 3](src/surmeier_lab_to_nwb/zhai2025/conversion_scripts/somatic_excitability/figure_3_somatic_excitability.py)
  - [Figure 6](src/surmeier_lab_to_nwb/zhai2025/conversion_scripts/somatic_excitability/figure_6_somatic_excitability.py)
  - [Figure 7](src/surmeier_lab_to_nwb/zhai2025/conversion_scripts/somatic_excitability/figure_7_somatic_excitability.py)
  - [Figure 8](src/surmeier_lab_to_nwb/zhai2025/conversion_scripts/somatic_excitability/figure_8_somatic_excitability.py)

### Spine Density Analysis (Figures 2, 4, 6, 7, 8)
- Two-photon and confocal microscopy
- 3D spine reconstruction and counting
- Multiple imaging platforms (Bruker, Nikon, Olympus)
- **Scripts**:
  - [Figure 2](src/surmeier_lab_to_nwb/zhai2025/conversion_scripts/spine_density/figure_2_spine_density.py)
  - [Figure 4](src/surmeier_lab_to_nwb/zhai2025/conversion_scripts/spine_density/figure_4_spine_density.py)
  - [Figure 4 Confocal (Nikon)](src/surmeier_lab_to_nwb/zhai2025/conversion_scripts/confocal_spine_density/figure_4_confocal_spine_density_nikon.py)
  - [Figure 4 Confocal (Olympus)](src/surmeier_lab_to_nwb/zhai2025/conversion_scripts/confocal_spine_density/figure_4h_confocal_spine_density_olympus.py)
  - [Figure 6](src/surmeier_lab_to_nwb/zhai2025/conversion_scripts/spine_density/figure_6_spine_density.py)
  - [Figure 7](src/surmeier_lab_to_nwb/zhai2025/conversion_scripts/spine_density/figure_7_spine_density.py)
  - [Figure 8](src/surmeier_lab_to_nwb/zhai2025/conversion_scripts/spine_density/figure_8_spine_density.py)

### Acetylcholine Biosensor (Figure 5)
- GRABACh3.0 fluorescent sensor recordings
- Two-photon line scans and full-field imaging
- Real-time acetylcholine dynamics
- **Script**: [Figure 5](src/surmeier_lab_to_nwb/zhai2025/conversion_scripts/acetylcholine_biosensor/figure_5_acetylcholine_biosensor.py)

### Behavioral Analysis (Figures 7, 8)
- AIM (Abnormal Involuntary Movement) scoring
- Contralateral rotation counting
- Video recording integration
- **Scripts**:
  - [Figure 7 AIM](src/surmeier_lab_to_nwb/zhai2025/conversion_scripts/aim_behavior/figure_7_behavioral_aim_experiments.py)
  - [Figure 7 Videos](src/surmeier_lab_to_nwb/zhai2025/conversion_scripts/videos/figure_7_behavioral_videos.py)
  - [Figure 8 AIM](src/surmeier_lab_to_nwb/zhai2025/conversion_scripts/aim_behavior/figure_8_behavioral_aim_experiments.py)
  - [Supplementary Figure 3 Videos](src/surmeier_lab_to_nwb/zhai2025/conversion_scripts/videos/supplementary_figure_3_behavioral_videos.py)

### Optogenetics (Figures 2, 4)
- Light-evoked postsynaptic currents (oEPSCs)
- Stimulus timing and metadata
- **Scripts**:
  - [Figure 2](src/surmeier_lab_to_nwb/zhai2025/conversion_scripts/optical_stimulation/figure_2_optical_stimuli.py)
  - [Figure 4](src/surmeier_lab_to_nwb/zhai2025/conversion_scripts/optical_stimulation/figure_4_optical_stimuli.py)
