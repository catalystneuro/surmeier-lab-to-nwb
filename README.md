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

### Running a specific conversion
Once you have installed the package, you can run any of the conversion scripts in a notebook or a python file:

Figure-specific conversion scripts are located in `src/surmeier_lab_to_nwb/zhai2025/conversion_scripts/`. Each script handles specific data types and experimental conditions:

```bash
# Example: Convert Figure 1 dendritic excitability data
python src/surmeier_lab_to_nwb/zhai2025/conversion_scripts/dendritic_excitability/figure_1_dendritic_excitability.py

# Example: Convert Figure 5 acetylcholine biosensor data
python src/surmeier_lab_to_nwb/zhai2025/conversion_scripts/acetylcholine_biosensor/figure_5_acetylcholine_biosensor.py
```

**Note:** All methods above install the repository in [editable mode](https://pip.pypa.io/en/stable/cli/pip_install/#editable-installs), allowing you to modify the source code if needed.


## Repository structure
Each conversion is organized in a directory of its own in the `src` directory:

    surmeier-lab-to-nwb/
    ├── LICENSE
    ├── README.md
    ├── pyproject.toml
    └── src/
        └── surmeier_lab_to_nwb/
            └── zhai2025/
                ├── conversion_notes_folder/  # Detailed conversion documentation
                ├── conversion_scripts/  # Figure-specific conversion scripts
                │   ├── acetylcholine_biosensor/
                │   ├── aim_behavior/
                │   ├── confocal_spine_density/
                │   ├── dendritic_excitability/
                │   ├── optical_stimulation/
                │   ├── somatic_excitability/
                │   ├── spine_density/
                │   └── videos/
                ├── interfaces/         # Custom NWB interfaces
                └── utils/             # Utility functions

## Figure-Specific Conversions

The `zhai2025` conversion includes specialized scripts for each figure in the paper:

### Dendritic Excitability (Figures 1, 3, 6, 7)
- Two-photon calcium imaging with patch-clamp electrophysiology
- Back-propagating action potential measurements
- Cell-type specific analysis (dSPNs vs iSPNs)
- **Scripts**:
  - [Figure 1](src/surmeier_lab_to_nwb/zhai2025/conversion_scripts/dendritic_excitability/figure_1_dendritic_excitability.py)
  - [Figure 3](src/surmeier_lab_to_nwb/zhai2025/conversion_scripts/dendritic_excitability/figure_3_dendritic_excitability.py)
  - [Figure 6](src/surmeier_lab_to_nwb/zhai2025/conversion_scripts/dendritic_excitability/figure_6_dendritic_excitability.py)
  - [Figure 7](src/surmeier_lab_to_nwb/zhai2025/conversion_scripts/dendritic_excitability/figure_7_dendritic_excitability.py)
  - [Figure 7 OxoM](src/surmeier_lab_to_nwb/zhai2025/conversion_scripts/dendritic_excitability/figure_7_oxoM_dendritic_excitability.py)

### Somatic Excitability (Figures 1, 3, 6, 7, 8)
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
