# Protocol-Figure Mapping for Zhai et al. 2025

This document provides a comprehensive mapping between experimental protocols and figure-specific data in Zhai et al. 2025, showing detailed relationships between protocols, conversion scripts, and experimental procedures.

## Protocol Files (Renamed for Clarity)

| New Name | Description | Original Filename |
|----------|-------------|-------------------|
| `dendritic_calcium_2PLSM.md` | Two-photon imaging of dendritic calcium and spines | Ex_vivo_imaging_of_dendritic_calcium_transients_and_spines_with_two-photon_laser_scanning_microscopy_2PLSM.md |
| `genetically_encoded_sensors_2PLSM.md` | Two-photon imaging of biosensors (GRABACh3.0) | Ex_vivo_imaging_of_genetically_encoded_sensors_with_two-photon_laser_scanning_microscopy_2PLSM.md |
| `patch_clamp_optogenetics.md` | Patch clamp with optogenetic stimulation | Ex_vivo_mouse_brain_patch_clamp_recordings_combined_with_optogenetic_stimulation.md |
| `confocal_microscopy.md` | Confocal imaging of spine density | High-resolution_confocal_microscopy_of_sparsely_labeled_neurons.md |
| `stereotaxic_viral_injection_v1.md` | Viral injection surgery (version 1) | Mouse_Stereotaxic_Surgeries_for_Intracranial_Viral_Injection_v1.md |
| `stereotaxic_viral_injection_v2.md` | Viral injection surgery (version 2) | Mouse_Stereotaxic_Surgeries_for_Intracranial_Viral_Injection_v2.md |
| `levodopa_dyskinesia_rating.md` | AIM behavioral scoring | Rating_of_levodopa-induced_dyskinesia_in_6-OHDA_lesioned_mice.md |
| `6OHDA_lesion_model.md` | Parkinson's disease model creation | Unilateral_6-hydroxydopamine_lesion_mouse_model_of_Parkinson_disease.md |

## Detailed Figure-Protocol-Script Mapping

### Figure 1: Dendritic and Somatic Excitability (dSPNs)

#### Figure 1A-E: Dendritic Excitability
- **Conversion Script**: `dendritic_excitability/figure_1_dendritic_excitability.py`
- **Experimental Methods**: Dendritic patch clamp recordings in current clamp mode
- **Related Protocols**:
  - `patch_clamp_optogenetics.md` - **PRIMARY** - Dendritic patch clamp methodology
  - `stereotaxic_viral_injection_v1.md` or `v2.md` - For Drd1-Tdtomato reporter expression
  - `6OHDA_lesion_model.md` - For creating PD model
- **Data Type**: Dendritic current clamp recordings with current injection protocols
- **Cell Type**: Direct pathway SPNs (dSPNs)
- **Treatment Groups**: Saline vs L-DOPA treated mice

#### Figure 1F-J: Somatic Excitability
- **Conversion Script**: `somatic_excitability/figure_1_somatic_excitability.py`
- **Experimental Methods**: Whole-cell patch clamp recordings in current clamp mode
- **Related Protocols**:
  - `patch_clamp_optogenetics.md` - **PRIMARY** - Somatic patch clamp methodology
  - `stereotaxic_viral_injection_v1.md` or `v2.md` - For Drd1-Tdtomato reporter expression
  - `6OHDA_lesion_model.md` - For creating PD model
- **Data Type**: Somatic current clamp recordings with F-I protocol (-120 to +300 pA steps)
- **Cell Type**: Direct pathway SPNs (dSPNs)
- **Treatment Groups**: Saline vs L-DOPA treated mice

### Figure 2: Spine Density and Optical Stimulation (dSPNs)

#### Figure 2A-D: Spine Density Analysis
- **Conversion Script**: `spine_density/figure_2_spine_density.py`
- **Experimental Methods**: Two-photon microscopy for spine density quantification
- **Related Protocols**:
  - `dendritic_calcium_2PLSM.md` - **PRIMARY** - Two-photon imaging methodology for spine imaging
  - `stereotaxic_viral_injection_v1.md` or `v2.md` - For sparse neuronal labeling
  - `6OHDA_lesion_model.md` - For PD model
- **Data Type**: Z-stack images (0.15 μm pixels, 0.3 μm z-steps) for 3D spine reconstruction
- **Cell Type**: Direct pathway SPNs (dSPNs)
- **Treatment Groups**: Saline vs L-DOPA treated mice

#### Figure 2E-H: Optical Stimulation
- **Conversion Script**: `optical_stimulation/figure_2_optical_stimuli.py`
- **Experimental Methods**: Two-photon imaging with optical stimulation protocols
- **Related Protocols**:
  - `dendritic_calcium_2PLSM.md` - **PRIMARY** - Two-photon imaging with optical stimulation
  - `patch_clamp_optogenetics.md` - For optogenetic stimulation protocols
  - `stereotaxic_viral_injection_v1.md` or `v2.md` - For ChR2 expression and reporter labeling
  - `6OHDA_lesion_model.md` - For PD model
- **Data Type**: Calcium imaging data with optogenetic stimulation protocols
- **Cell Type**: Direct pathway SPNs (dSPNs)
- **Stimulation**: ChR2-mediated optical activation

### Figure 3: Dendritic and Somatic Excitability (iSPNs)

#### Figure 3A-E: Dendritic Excitability (iSPNs)
- **Conversion Script**: `dendritic_excitability/figure_3_dendritic_excitability.py`
- **Experimental Methods**: Dendritic patch clamp recordings
- **Related Protocols**:
  - `patch_clamp_optogenetics.md` - **PRIMARY** - Dendritic patch clamp methodology
  - `stereotaxic_viral_injection_v1.md` or `v2.md` - For Drd2-EGFP reporter expression
  - `6OHDA_lesion_model.md` - For PD model
- **Data Type**: Dendritic current clamp recordings
- **Cell Type**: Indirect pathway SPNs (iSPNs)
- **Treatment Groups**: Saline vs L-DOPA treated mice

#### Figure 3F-J: Somatic Excitability (iSPNs)
- **Conversion Script**: `somatic_excitability/figure_3_somatic_excitability.py`
- **Experimental Methods**: Whole-cell patch clamp recordings in current clamp mode
- **Related Protocols**:
  - `patch_clamp_optogenetics.md` - **PRIMARY** - Somatic patch clamp methodology
  - `stereotaxic_viral_injection_v1.md` or `v2.md` - For Drd2-EGFP reporter expression
  - `6OHDA_lesion_model.md` - For PD model
- **Data Type**: Somatic current clamp recordings with F-I protocol
- **Cell Type**: Indirect pathway SPNs (iSPNs)
- **Treatment Groups**: Saline vs L-DOPA treated mice

### Figure 4: Comprehensive Spine Density Analysis

#### Figure 4A-D: Two-Photon Spine Density
- **Conversion Script**: `spine_density/figure_4_spine_density.py`
- **Experimental Methods**: Two-photon microscopy spine density analysis
- **Related Protocols**:
  - `dendritic_calcium_2PLSM.md` - **PRIMARY** - Two-photon spine imaging
  - `stereotaxic_viral_injection_v1.md` or `v2.md` - For sparse labeling
  - `6OHDA_lesion_model.md` - For PD model
- **Data Type**: High-resolution z-stack images for spine quantification
- **Cell Types**: Both dSPNs and iSPNs
- **Treatment Groups**: Saline vs L-DOPA treated mice

#### Figure 4E-F: Optical Stimulation Analysis
- **Conversion Script**: `optical_stimulation/figure_4_optical_stimuli.py`
- **Experimental Methods**: Two-photon imaging with optical stimulation
- **Related Protocols**:
  - `dendritic_calcium_2PLSM.md` - **PRIMARY** - Two-photon imaging methodology
  - `patch_clamp_optogenetics.md` - For optogenetic protocols
  - `stereotaxic_viral_injection_v1.md` or `v2.md` - For ChR2 and reporter expression
  - `6OHDA_lesion_model.md` - For PD model
- **Data Type**: Optical stimulation response data
- **Cell Types**: Both dSPNs and iSPNs

#### Figure 4G: Confocal Spine Density (Nikon)
- **Conversion Script**: `confocal_spine_density/figure_4_confocal_spine_density_nikon.py`
- **Experimental Methods**: High-resolution confocal microscopy (Nikon system)
- **Related Protocols**:
  - `confocal_microscopy.md` - **PRIMARY** - Confocal imaging methodology
  - `stereotaxic_viral_injection_v1.md` or `v2.md` - For sparse neuronal labeling
  - `6OHDA_lesion_model.md` - For PD model
- **Data Type**: High-resolution confocal images for spine analysis
- **Microscope**: Nikon confocal system
- **Cell Types**: Both dSPNs and iSPNs

#### Figure 4H: Confocal Spine Density (Olympus)
- **Conversion Script**: `confocal_spine_density/figure_4h_confocal_spine_density_olympus.py`
- **Experimental Methods**: High-resolution confocal microscopy (Olympus system)
- **Related Protocols**:
  - `confocal_microscopy.md` - **PRIMARY** - Confocal imaging methodology
  - `stereotaxic_viral_injection_v1.md` or `v2.md` - For sparse neuronal labeling
  - `6OHDA_lesion_model.md` - For PD model
- **Data Type**: High-resolution confocal images for spine analysis
- **Microscope**: Olympus confocal system
- **Cell Types**: Both dSPNs and iSPNs

### Figure 5: Acetylcholine Biosensor Imaging

#### Figure 5A-H: GRABACh3.0 Biosensor Analysis
- **Conversion Script**: `acetylcholine_biosensor/figure_5_acetylcholine_biosensor.py`
- **Experimental Methods**: GRABACh3.0 biosensor imaging with two-photon microscopy
- **Related Protocols**:
  - `genetically_encoded_sensors_2PLSM.md` - **PRIMARY** - Biosensor imaging methodology
  - `stereotaxic_viral_injection_v1.md` or `v2.md` - For AAV5-hSyn-ACh3.0 injection
  - `6OHDA_lesion_model.md` - For PD model
- **Data Type**: Brightness-over-time (BOT) measurements of acetylcholine release dynamics
- **Biosensor**: GRABACh3.0 genetically encoded acetylcholine sensor
- **Brain Region**: Dorsal striatum
- **Treatment Groups**: Saline vs L-DOPA treated mice

### Figure 6: Combined Excitability and Spine Analysis

#### Figure 6A-C: Dendritic Excitability
- **Conversion Script**: `dendritic_excitability/figure_6_dendritic_excitability.py`
- **Experimental Methods**: Dendritic patch clamp recordings
- **Related Protocols**:
  - `patch_clamp_optogenetics.md` - **PRIMARY** - Dendritic patch clamp methodology
  - `stereotaxic_viral_injection_v1.md` or `v2.md` - For cell identification
  - `6OHDA_lesion_model.md` - For PD model
- **Data Type**: Dendritic current clamp recordings
- **Cell Types**: Both dSPNs and iSPNs

#### Figure 6D-F: Somatic Excitability
- **Conversion Script**: `somatic_excitability/figure_6_somatic_excitability.py`
- **Experimental Methods**: Somatic patch clamp recordings
- **Related Protocols**:
  - `patch_clamp_optogenetics.md` - **PRIMARY** - Somatic patch clamp methodology
  - `stereotaxic_viral_injection_v1.md` or `v2.md` - For cell identification
  - `6OHDA_lesion_model.md` - For PD model
- **Data Type**: Somatic current clamp recordings
- **Cell Types**: Both dSPNs and iSPNs

#### Figure 6G-I: Spine Density Analysis
- **Conversion Script**: `spine_density/figure_6_spine_density.py`
- **Experimental Methods**: Two-photon spine density measurement
- **Related Protocols**:
  - `dendritic_calcium_2PLSM.md` - **PRIMARY** - Two-photon spine imaging
  - `stereotaxic_viral_injection_v1.md` or `v2.md` - For sparse labeling
  - `6OHDA_lesion_model.md` - For PD model
- **Data Type**: Spine density quantification
- **Cell Types**: Both dSPNs and iSPNs

### Figure 7: Behavioral Assessment and Physiology

#### Figure 7A-C: AIM Behavioral Scoring
- **Conversion Script**: `aim_behavior/figure_7_behavioral_aim_experiments.py`
- **Experimental Methods**: AIM (Abnormal Involuntary Movement) behavioral scoring
- **Related Protocols**:
  - `levodopa_dyskinesia_rating.md` - **PRIMARY** - AIM behavioral scoring methodology
  - `6OHDA_lesion_model.md` - For PD model creation
- **Data Type**: AIM scores (locomotion, axial, limb, orolingual dyskinesia)
- **Assessment**: Behavioral rating scales for dyskinesia severity

#### Figure 7D: Behavioral Video Analysis
- **Conversion Script**: `videos/figure_7_behavioral_videos.py`
- **Experimental Methods**: Video analysis of dyskinetic behaviors
- **Related Protocols**:
  - `levodopa_dyskinesia_rating.md` - **PRIMARY** - Video-based behavioral assessment
  - `6OHDA_lesion_model.md` - For PD model
- **Data Type**: Video recordings and automated movement analysis
- **Analysis**: Computer-assisted behavioral quantification

#### Figure 7E-G: Dendritic Excitability
- **Conversion Script**: `dendritic_excitability/figure_7_dendritic_excitability.py`
- **Experimental Methods**: Dendritic patch clamp recordings
- **Related Protocols**:
  - `patch_clamp_optogenetics.md` - **PRIMARY** - Dendritic patch clamp methodology
  - `stereotaxic_viral_injection_v1.md` or `v2.md` - For reporter expression
  - `6OHDA_lesion_model.md` - For PD model
- **Data Type**: Dendritic excitability measurements correlated with behavior
- **Cell Types**: Both dSPNs and iSPNs

#### Figure 7H-J: OxoM Dendritic Excitability
- **Conversion Script**: `dendritic_excitability/figure_7_oxoM_dendritic_excitability.py`
- **Experimental Methods**: Dendritic patch clamp with oxotremorine-M (oxoM) application
- **Related Protocols**:
  - `patch_clamp_optogenetics.md` - **PRIMARY** - Dendritic patch clamp with pharmacology
  - `stereotaxic_viral_injection_v1.md` or `v2.md` - For cell identification
  - `6OHDA_lesion_model.md` - For PD model
- **Data Type**: Dendritic excitability under muscarinic receptor activation
- **Pharmacology**: OxoM (muscarinic agonist) treatment
- **Cell Types**: Both dSPNs and iSPNs

#### Figure 7K-M: Somatic Excitability
- **Conversion Script**: `somatic_excitability/figure_7_somatic_excitability.py`
- **Experimental Methods**: Somatic patch clamp recordings
- **Related Protocols**:
  - `patch_clamp_optogenetics.md` - **PRIMARY** - Somatic patch clamp methodology
  - `stereotaxic_viral_injection_v1.md` or `v2.md` - For cell identification
  - `6OHDA_lesion_model.md` - For PD model
- **Data Type**: Somatic excitability correlated with behavioral measures
- **Cell Types**: Both dSPNs and iSPNs

#### Figure 7N-P: Spine Density Analysis
- **Conversion Script**: `spine_density/figure_7_spine_density.py`
- **Experimental Methods**: Spine density measurement in behaviorally characterized animals
- **Related Protocols**:
  - `dendritic_calcium_2PLSM.md` - **PRIMARY** - Two-photon spine imaging
  - `stereotaxic_viral_injection_v1.md` or `v2.md` - For sparse labeling
  - `6OHDA_lesion_model.md` - For PD model
- **Data Type**: Spine density measurements correlated with AIM scores
- **Cell Types**: Both dSPNs and iSPNs

### Figure 8: Comprehensive Multi-Modal Analysis

#### Figure 8A-C: AIM Behavioral Assessment
- **Conversion Script**: `aim_behavior/figure_8_behavioral_aim_experiments.py`
- **Experimental Methods**: Comprehensive AIM behavioral scoring
- **Related Protocols**:
  - `levodopa_dyskinesia_rating.md` - **PRIMARY** - AIM behavioral scoring
  - `6OHDA_lesion_model.md` - For PD model
- **Data Type**: Multi-session AIM scores and behavioral progression
- **Assessment**: Longitudinal behavioral tracking

#### Figure 8D-F: Somatic Excitability
- **Conversion Script**: `somatic_excitability/figure_8_somatic_excitability.py`
- **Experimental Methods**: Somatic patch clamp in behaviorally characterized animals
- **Related Protocols**:
  - `patch_clamp_optogenetics.md` - **PRIMARY** - Somatic patch clamp methodology
  - `stereotaxic_viral_injection_v1.md` or `v2.md` - For cell identification
  - `6OHDA_lesion_model.md` - For PD model
- **Data Type**: Somatic excitability in animals with known dyskinesia severity
- **Cell Types**: Both dSPNs and iSPNs

#### Figure 8G-I: Spine Density Correlation
- **Conversion Script**: `spine_density/figure_8_spine_density.py`
- **Experimental Methods**: Spine density analysis correlated with behavior and physiology
- **Related Protocols**:
  - `dendritic_calcium_2PLSM.md` - **PRIMARY** - Two-photon spine imaging
  - `confocal_microscopy.md` - For high-resolution spine analysis
  - `stereotaxic_viral_injection_v1.md` or `v2.md` - For sparse labeling
  - `6OHDA_lesion_model.md` - For PD model
- **Data Type**: Comprehensive spine density analysis with multi-modal correlations
- **Cell Types**: Both dSPNs and iSPNs

## Supplementary Figures

### Supplementary Figure 3: Extended Behavioral Video Analysis
- **Conversion Script**: `videos/supplementary_figure_3_behavioral_videos.py`
- **Experimental Methods**: Extended video analysis of dyskinetic behaviors
- **Related Protocols**:
  - `levodopa_dyskinesia_rating.md` - **PRIMARY** - Video-based behavioral assessment
  - `6OHDA_lesion_model.md` - For PD model
- **Data Type**: Extended video datasets and automated analysis
- **Analysis**: Computer vision-based movement quantification

## Protocol Usage Matrix

| Protocol | Fig 1 | Fig 2 | Fig 3 | Fig 4 | Fig 5 | Fig 6 | Fig 7 | Fig 8 | SF3 |
|----------|-------|-------|-------|-------|-------|-------|-------|-------|-----|
| `6OHDA_lesion_model.md` | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| `stereotaxic_viral_injection_v1/v2.md` | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| `patch_clamp_optogenetics.md` | ✓ | ○ | ✓ | ○ | - | ✓ | ✓ | ✓ | - |
| `dendritic_calcium_2PLSM.md` | - | ✓ | - | ✓ | - | ✓ | ✓ | ✓ | - |
| `genetically_encoded_sensors_2PLSM.md` | - | - | - | - | ✓ | - | - | - | - |
| `confocal_microscopy.md` | - | - | - | ✓ | - | - | - | ✓ | - |
| `levodopa_dyskinesia_rating.md` | - | - | - | - | - | - | ✓ | ✓ | ✓ |

**Legend**: ✓ = Primary use, ○ = Secondary use (optical stimulation), - = Not used

## Conversion Script Organization

### By Experimental Category
```
conversion_scripts/
├── dendritic_excitability/       # Dendritic patch clamp experiments
│   ├── figure_1_dendritic_excitability.py
│   ├── figure_3_dendritic_excitability.py
│   ├── figure_6_dendritic_excitability.py
│   ├── figure_7_dendritic_excitability.py
│   └── figure_7_oxoM_dendritic_excitability.py
├── somatic_excitability/         # Somatic patch clamp experiments
│   ├── figure_1_somatic_excitability.py
│   ├── figure_3_somatic_excitability.py
│   ├── figure_6_somatic_excitability.py
│   ├── figure_7_somatic_excitability.py
│   └── figure_8_somatic_excitability.py
├── spine_density/                # Two-photon spine density analysis
│   ├── figure_2_spine_density.py
│   ├── figure_4_spine_density.py
│   ├── figure_6_spine_density.py
│   ├── figure_7_spine_density.py
│   └── figure_8_spine_density.py
├── optical_stimulation/          # Optical stimulation experiments
│   ├── figure_2_optical_stimuli.py
│   └── figure_4_optical_stimuli.py
├── confocal_spine_density/       # Confocal microscopy spine analysis
│   ├── figure_4_confocal_spine_density_nikon.py
│   └── figure_4h_confocal_spine_density_olympus.py
├── acetylcholine_biosensor/      # Biosensor imaging experiments
│   └── figure_5_acetylcholine_biosensor.py
├── aim_behavior/                 # Behavioral assessment experiments
│   ├── figure_7_behavioral_aim_experiments.py
│   └── figure_8_behavioral_aim_experiments.py
└── videos/                       # Video analysis experiments
    ├── figure_7_behavioral_videos.py
    └── supplementary_figure_3_behavioral_videos.py
```

## Key Experimental Variables

### Cell Type Identification
- **dSPNs**: Drd1-Tdtomato reporter mice
- **iSPNs**: Drd2-EGFP reporter mice

### Treatment Groups
- **Saline**: Control group (vehicle injection)
- **L-DOPA**: L-DOPA + carbidopa treated animals (dyskinetic)

### Brain Regions
- **Dorsal Striatum**: Primary region of interest
- **Lesioned Side**: Ipsilateral to 6-OHDA injection
- **Control Side**: Contralateral hemisphere

### Recording Conditions
- **Ex vivo**: Brain slice preparations (300-350 μm thickness)
- **Temperature**: 32-34°C
- **ACSF**: Standard artificial cerebrospinal fluid

## Missing/Inferred Protocols

The following protocols are referenced in experimental methods but not available as standalone files:

1. **Brain slice preparation** - Standard ex vivo slice protocol
2. **Perfusion and fixation** - For histological processing
3. **Immunofluorescence staining** - Post-hoc cell identification
4. **Drug preparation and application** - L-DOPA, oxoM, etc.
5. **Data analysis pipelines** - Statistical and imaging analysis methods

## Notes

1. **Protocol Versions**: Two versions of viral injection protocols exist - usage may vary by experiment
2. **Adaptation Required**: Protocols may need adaptation for specific experimental configurations
3. **Combination Protocols**: Most experiments require multiple protocols in sequence
4. **Quality Control**: All experiments include appropriate controls and validation steps
5. **Reproducibility**: Protocol standardization ensures experimental reproducibility across figures

---

*This comprehensive mapping integrates information from conversion scripts, experimental methods, and protocol files to provide complete traceability between protocols and figure-specific data in Zhai et al. 2025.*
