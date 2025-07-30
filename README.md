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

### Basic installation

You can install the latest release of the package with pip:

```
pip install surmeier-lab-to-nwb
```

We recommend that you install the package inside a [virtual environment](https://docs.python.org/3/tutorial/venv.html). A simple way of doing this is to use a [conda environment](https://docs.conda.io/projects/conda/en/latest/user-guide/concepts/environments.html) from the `conda` package manager ([installation instructions](https://docs.conda.io/en/latest/miniconda.html)). Detailed instructions on how to use conda environments can be found in their [documentation](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html).

### Running a specific conversion
Once you have installed the package with pip, you can run any of the conversion scripts in a notebook or a python file:

Figure-specific conversion scripts are located in `src/surmeier_lab_to_nwb/zhai2025/conversion_scripts/`. Each script handles specific data types and experimental conditions:

```bash
# Example: Convert Figure 1 dendritic excitability data
python src/surmeier_lab_to_nwb/zhai2025/conversion_scripts/dendritic_excitability/figure_1_dendritic_excitability.py

# Example: Convert Figure 5 acetylcholine biosensor data
python src/surmeier_lab_to_nwb/zhai2025/conversion_scripts/acetylcholine_biosensor/figure_5_acetylcholine_biosensor.py
```

## Installation from GitHub
Another option is to install the package directly from GitHub. This option has the advantage that the source code can be modified if you need to amend some of the code we originally provided to adapt to future experimental differences. To install the conversion from GitHub you will need to use `git` ([installation instructions](https://github.com/git-guides/install-git)). We also recommend the installation of `conda` ([installation instructions](https://docs.conda.io/en/latest/miniconda.html)) as it contains all the required machinery in a single and simple install.

From a terminal (note that conda should install one in your system) you can do the following:

```
git clone https://github.com/catalystneuro/surmeier-lab-to-nwb
cd surmeier-lab-to-nwb
conda env create --file make_env.yml
conda activate surmeier-lab-to-nwb-env
```

This creates a [conda environment](https://docs.conda.io/projects/conda/en/latest/user-guide/concepts/environments.html) which isolates the conversion code from your system libraries.  We recommend that you run all your conversion related tasks and analysis from the created environment in order to minimize issues related to package dependencies.

Alternatively, if you want to avoid conda altogether (for example if you use another virtual environment tool) you can install the repository with the following commands using only pip:

```
git clone https://github.com/catalystneuro/surmeier-lab-to-nwb
cd surmeier-lab-to-nwb
pip install --editable .
```

Note:
both of the methods above install the repository in [editable mode](https://pip.pypa.io/en/stable/cli/pip_install/#editable-installs).


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
