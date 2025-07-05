# Figure 6 Conversion Notes

## Overview
Figure 6: "M1R mediated the dendritic, but not the somatic, alterations in iSPNs in LID off-state"

This figure demonstrates that M1 muscarinic receptors (M1Rs) specifically mediate dendritic, but not somatic, alterations in indirect pathway spiny projection neurons (iSPNs) during the off-state of levodopa-induced dyskinesia (LID). The data includes:

- **Somatic excitability**: Current injection protocols testing M1R contribution to somatic changes
- **Dendritic excitability**: Combined patch clamp + two-photon line scan imaging with M1R pharmacological manipulation
- **Spine density**: Two-photon microscopy Z-stacks examining M1R effects on dendritic spine morphology
- **Pharmacological intervention**: Systemic M1R antagonist treatment during LID off-state

## Data Structure

**Raw Data Location**: `/home/heberto/development/surmeier-lab-to-nwb/link_to_raw_data/Figure 6/`

### Overview

```
Figure 6/
├── Dendritic excitability/
│   ├── control/ (9 animals)
│   └── M1R antagonist/ (7 animals)
├── Somatic excitability/
│   ├── control/ (8 animals)
│   └── M1R antagonist/ (10 animals)
└── Spine density/
    ├── control/ (8 animals)
    └── M1R antagonist/ (9 animals)
```

### Detailed Structure

#### Somatic Excitability

**File Organization**: Current injection experiments with 21 current steps per animal

```
Somatic excitability/
├── control/
│   ├── 0217a/ → cell1-001 through cell1-021 (-120 pA to +300 pA, 20 pA steps)
│   ├── 0217b/ → cell2-001 through cell2-021
│   ├── 0225a, 0225b, 0225c/
│   └── 0811a, 0812a, 0812b/
└── M1R antagonist/
    ├── 0218b/ → cell2-001 through cell2-021 (+ THP treatment)
    ├── 0218c/ → cell3-001 through cell3-021 (+ THP treatment)
    ├── 0220a, 0220c, 0312a, 0312b/
    └── 0721a, 0721b, 0721c, 0722b/
```

**File Bundle per Current Level**:
```
cell1-001/
├── cell1-001.xml                                    # Experiment metadata
├── cell1-001_Cycle00001_VoltageOutput_001.xml       # Stimulus protocol
├── cell1-001_Cycle00001_VoltageRecording_001.csv    # Voltage recording data
└── cell1-001_Cycle00001_VoltageRecording_001.xml    # Recording metadata
```

**Animal-Specific Cell Numbering**:
- **Animal "a"**: cell1-XXX, **Animal "b"**: cell2-XXX, **Animal "c"**: cell3-XXX

#### Dendritic Excitability

**File Organization**: Line scan imaging + electrophysiology with location-specific recordings

```
Dendritic excitability/
├── control/
│   ├── 0217a/ → 02172020_Cell1_[dist1|dist2|prox1|prox2]_trio-[001|002|003]
│   ├── 0217b, 0225a, 0225b/
│   └── 0425a, 0425b, 0811a, 0812a, 0812b/
└── M1R antagonist/
    ├── 0218b/ → 02182020_Cell2_[dist1|dist2|prox1|prox2]_trio-[001|002|003] (+ THP)
    ├── 0218c/ → 02182020_Cell3_[dist1|dist2|prox1|prox2]_trio-[001|002|003] (+ THP)
    └── 0220c, 0312a, 0312b, 0430a, 0430b/
```

**Naming Convention**: `[MMDDYYYY]_Cell[N]_[location][number]_trio-[trial]`
- **Locations**: dist (distal ~90μm), prox (proximal ~40μm)
- **Trials**: 3 repetitions per location (001, 002, 003)

**File Bundle per Recording**:
```
02182020_Cell2_dist1_trio-001/
├── *_Cycle00001_Ch1_000001.ome.tif          # Fluo-4 calcium channel
├── *-Cycle00001_Ch1Source.tif               # Field of view with scan line
├── *_Cycle00001_Ch2_000001.ome.tif          # Alexa Fluor 568 structural channel
├── *-Cycle00001_Ch2Source.tif               # Structural channel field of view
├── *_Cycle00001_LineProfileData.csv         # Averaged fluorescence traces
├── *_Cycle00001_VoltageOutput_001.xml       # Stimulus protocol
├── *_Cycle00001_VoltageRecording_001.csv    # Electrophysiology data
├── *_Cycle00001_VoltageRecording_001.xml    # Recording metadata
├── *.env                                    # Environment configuration
└── *.xml                                    # Master experiment file
```

#### Spine Density

**File Organization**: Z-stack imaging for spine morphology analysis

```
Spine density/
├── control/
│   ├── 02172020a/ → Decon_20200217_Cell1_[dist12|prox1|prox2]/
│   ├── 02172020b, 02252020a, 02252020b, 02252020c/
│   └── 08112020a, 08122020a, 08122020b/
└── M1R antagonist/
    ├── 02182020a/ → Decon_20200218_Cell1_[dist1|dist2|prox12]/ (+ THP)
    ├── 02182020c, 02202020c, 03122020a, 03122020b/
    └── 07212020a, 07212020b, 07212020c, 07222020b/
```

**Z-Stack Structure**:
```
Decon_20200217_Cell1_dist12/
├── 20_ZSeries-*_#.ome_Z01.tif    # Z-plane 1
├── 20_ZSeries-*_#.ome_Z02.tif    # Z-plane 2
├── 20_ZSeries-*_#.ome_Z03.tif    # Z-plane 3
└── ...                           # Additional Z-planes
```

**Imaging Parameters**:
- **Pixel size**: 0.15 μm, **Z-step**: 0.3 μm
- **Analysis**: Semi-automated spine counting with NeuronStudio

## Experimental Protocols

### Somatic Excitability Protocol

**Methodology**: Whole-cell patch clamp recording in current clamp mode from iSPNs in LID off-state mice

**Experimental Design**:
- **Test Current**: 200 pA current injections (500 ms duration) used for comparison
- **Animal Model**: Dyskinetic Drd2 BAC mice (expressing markers for iSPN identification)
- **Pharmacological Treatment**: M1R antagonist trihexyphenidyl hydrochloride (THP, 3 mg/kg, i.p.) vs saline control
- **Timing**: THP administered at beginning of off-state (~2 hr after levodopa), sacrifice 16 hours later

**Current Steps** (based on directory structure):
- Range: -120 pA to +300 pA
- Step size: 20 pA
- Duration: 500 ms each
- Total steps: 21 (from cell1-001 to cell1-021)

**Analysis**:
- F-I relationship curves comparing control vs M1R antagonist treatment
- **Key Finding**: M1R antagonism did NOT significantly impact somatic excitability elevation in off-state iSPNs
- Statistical analysis: Mann-Whitney test (n=8-10 cells from 4-5 mice per group)

### Dendritic Excitability Protocol

**Methodology**: Combination of patch clamp electrophysiology and two-photon laser scanning microscopy (2PLSM) from iSPNs in LID off-state mice

**Experimental Design**:
- **Animal Model**: Dyskinetic Drd2 BAC mice (iSPNs identified by lack of Drd1-Tdtomato expression)
- **Pharmacological Treatment**:
  - **In vivo**: THP (3 mg/kg, i.p.) at beginning of off-state, sacrifice 16 hours later
  - **Ex vivo**: Continuous bath application of selective M1R antagonist VU 0255035 (5 μM)
- **Recording Conditions**: Sagittal brain slices (280 μm), maintained at 31-32°C

**Methodology**: Same as Figures 1 and 3 - iSPNs patch clamped, filled with Ca2+-sensitive dye Fluo-4 and Ca2+-insensitive dye Alexa Fluor 568, and stimulated with brief current steps (three 2 nA injections, 2 ms each, at ~50 Hz).

**Stimulus Protocol**:
```
Current (nA)
  2 |   ___      ___      ___
    |  |   |    |   |    |   |
  0 |__|   |____|   |____|   |____
      <-2ms->    <-2ms->    <-2ms->
      <---20ms--><---20ms--->
```

**Physiological Rationale**:
- **Back-propagating action potentials (bAPs)**: Somatically delivered current steps evoke spikes that back-propagate into iSPN dendrites
- **M1R dependency testing**: Assessing whether M1R signaling contributes to LID-induced dendritic adaptations
- **Dual pharmacological approach**: In vivo THP + ex vivo VU 0255035 for comprehensive M1R blockade

**Analysis**:
- Ca2+ signals measured at proximal (~40 μm from soma) and distal (~90 μm from soma) dendritic locations
- **Dendritic excitability index**: Distal/proximal ratio used to quantify dendritic excitability
- **Statistical Analysis**: Mann-Whitney test (n=7-9 cells from 4-5 mice per group)

**Key Finding**: M1R antagonism prevented the off-state elevation in dendritic excitability (p = 0.012), indicating M1R-dependent dendritic adaptations in LID

### Spine Density Protocol

**Methodology**: Two-photon laser scanning microscopy with deconvolution from iSPNs in LID off-state mice

**Experimental Design**:
- **Animal Model**: Same as dendritic excitability experiments (Dyskinetic Drd2 BAC mice)
- **Pharmacological Treatment**: THP (3 mg/kg, i.p.) vs saline control at beginning of off-state
- **Timing**: 16-hour treatment period before sacrifice and slice preparation

**Imaging Parameters**: Same as Figure 2:
- **Pixel size**: 0.15 μm
- **Z-step size**: 0.3 μm
- **Analysis**: Semi-automated spine counting with NeuronStudio

**Statistical Analysis**:
- Mann-Whitney test (n=8-9 cells from 4-5 mice per group)
- **Key Finding**: M1R antagonism prevented the apparent elevation in spine density (p = 0.0003)

**Interpretation**: Suggests that M1R signaling contributes to LID-induced morphological adaptations in iSPN dendrites

### Pharmacological Manipulation

**M1R Antagonist Treatment Protocol**:

**Primary Drug - Trihexyphenidyl hydrochloride (THP)**:
- **Dose**: 3 mg/kg body weight
- **Route**: Intraperitoneal injection
- **Timing**: Administered at beginning of LID off-state (~2 hr after levodopa administration)
- **Duration**: 16-hour treatment period before sacrifice
- **Rationale**: Good brain bioavailability and favorable pharmacokinetics for peripheral administration
- **Limitation**: Not fully selective for M1Rs (broad anticholinergic effects)

**Secondary Drug - VU 0255035**:
- **Concentration**: 5 μM in bath solution
- **Application**: Continuous superfusion during ex vivo slice recordings
- **Selectivity**: Selective M1R antagonist for more specific receptor blockade
- **Purpose**: Complement THP treatment with highly selective M1R antagonism

**Control Treatment**:
- **Vehicle**: Saline injection (matched volume and timing)
- **Comparison**: Direct within-subject comparison of M1R-dependent vs independent effects


## Key Metadata

### Stimulus Parameters
Same as Figure 1 - three 2 nA current injections, 2 ms each, at ~50 Hz frequency.

### Channel Assignments
- **Channel 1 (Ch1)**: Fluo-4 calcium-sensitive dye (green, 490-560 nm)
- **Channel 2 (Ch2)**: Alexa Fluor 568 structural dye (red, 580-620 nm)

### Line Scan Parameters
- **Pixels per line**: 64
- **Dwell time**: 10 μs/pixel
- **Typical acquisition**: 2500 lines (time points)

## Acquisition Parameters

### Electrophysiology
> All the electrophysiological recordings were made using a MultiClamp 700B amplifier (Axon Instrument, USA), and signals were filtered at 2 kHz and digitized at 10 kHz. Voltage protocols and data acquisition were performed by PrairieView 5.3 (Bruker).

### Two-Photon Imaging
> The recorded SPN was visualized using 810 nm excitation laser (Chameleon Ultra II, Coherent, Santa Clara, USA). Dendritic structure was visualized by the red signal of Alexa Fluor 568 detected by a Hamamatsu R3982 side-on photomultiplier tube (PMT, 580-620 nm). Calcium transients, as signals in the green channel, were detected by a Hamamatsu H7422P-40 GaAsP PMT (490-560 nm, Hamamatsu Photonics, Japan).

## Key Findings and Clinical Significance

### Primary Results Summary

**Differential M1R Effects**:
1. **Somatic Excitability**: M1R antagonism did NOT significantly impact the elevation in iSPN somatic excitability during LID off-state
2. **Dendritic Excitability**: M1R antagonism prevented the off-state elevation in dendritic excitability (p = 0.012)
3. **Spine Density**: M1R antagonism prevented the apparent elevation in spine density (p = 0.0003)

**Key Interpretation**: M1R signaling specifically contributes to LID-induced dendritic, but not somatic, adaptations in iSPNs

### Mechanistic Insights

**M1R-Dependent Plasticity**:
- M1R signaling is required for pathological dendritic adaptations in LID off-state
- Dendritic compartment shows selective vulnerability to M1R-mediated changes
- Spine morphology changes are M1R-dependent, suggesting structural plasticity involvement

**Therapeutic Implications**:
- M1R signaling in iSPNs represents potential therapeutic target for LID treatment
- Selective targeting of dendritic M1R effects might preserve beneficial somatic functions
- **Limitation**: Global M1R antagonism would cause unacceptable side effects, requiring targeted approaches

**Research Context**:
- Complements Figure 3 (D2R antagonist effects) by examining M1R contribution
- Led to follow-up studies with CDGI knockout mice and CRISPR-based M1R deletion for more specificity
- Part of comprehensive analysis of receptor systems involved in LID pathophysiology

## Experimental Design Notes

### Cell Type Identification
- **Animal Model**: Drd2 BAC transgenic mice (indirect pathway marker expression)
- **iSPN Identification**: Lack of Drd1-Tdtomato expression (negative selection for indirect pathway)
- **Confirmation**: Electrophysiological properties and morphological characteristics
- **LID Model**: 6-OHDA lesioning + chronic levodopa treatment to induce dyskinesia

### Experimental Timeline and Controls
- **Lesioning**: Unilateral 6-OHDA injection in medial forebrain bundle
- **LID Induction**: Chronic levodopa treatment (dyskinesiogenic doses)
- **Treatment Window**: M1R antagonist during off-state (critical timing)
- **Controls**: Saline vehicle injection with matched timing and volume
- **Recording Conditions**: Ex vivo brain slices (280 μm) at 31-32°C

### Comparison with Other Figures
- **Figure 3**: D2R antagonist (sulpiride) effects on iSPNs - tests dopaminergic modulation
- **Figure 6**: M1R antagonist effects on iSPNs - tests cholinergic modulation
- **Complementary Analysis**: Different neurotransmitter systems in same cell type and disease model
- **Mechanistic Progression**: Led to Figure 7 (CDGI knockout) for more specific targeting

### Statistical Approach
- **Design**: Between-subjects comparison (control vs M1R antagonist treatment)
- **Analysis**: Non-parametric Mann-Whitney test (appropriate for small sample sizes)
- **Sample Sizes**: n=7-10 cells from 4-5 mice per group (adequate power for effect sizes observed)
- **Significance**: p < 0.05 threshold

## Related Files
- [Overview](../conversion_notes_overview.md)
- [Figure 1 Notes](figure_1_conversion_notes.md) - Similar experimental design
- [Figure 3 Notes](figure_3_conversion_notes.md) - D2R antagonist effects on iSPNs
- [Figure 7 Notes](figure_7_conversion_notes.md) - CDGI knockout effects
