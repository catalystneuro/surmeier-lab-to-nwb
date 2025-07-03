# Figure 6 Conversion Notes

## Overview
Figure 6 examines the role of M1 muscarinic receptors in iSPN dendritic excitability changes. The data includes:
- **Somatic excitability**: Current injection protocols with and without M1R antagonist
- **Dendritic excitability**: Combined patch clamp + two-photon line scan imaging
- **Spine density**: Two-photon microscopy Z-stacks of dendritic segments
- **Pharmacological manipulation**: Effects of M1R antagonist on iSPN properties

## Data Structure

```
Figure 6/
├── Dendritic excitability/
│   ├── control/
│   └── M1R antagonist/
├── Somatic excitability/
│   ├── control/
│   └── M1R antagonist/
└── Spine density/
    ├── control/
    └── M1R antagonist/
```

## Experimental Protocols

### Somatic Excitability Protocol

**Methodology**: Whole-cell patch clamp recording in current clamp mode

Similar to Figures 1 and 3, voltage changes recorded in response to current injection steps of 500 ms duration.

**Current Steps**:
- Range: -120 pA to +300 pA
- Step size: 20 pA
- Duration: 500 ms each
- Total steps: 21 (from cell1-001 to cell1-021)

**Analysis**: F-I relationship curves and rheobase measurements

### Dendritic Excitability Protocol

**Methodology**: Combination of patch clamp electrophysiology and two-photon laser scanning microscopy (2PLSM)

> To assess the role of M1 muscarinic receptors in iSPN dendritic excitability, the same methodology as Figures 1 and 3 was applied. iSPNs were patch clamped, filled with Ca2+-sensitive dye Fluo-4 and Ca2+-insensitive dye Alexa Fluor 568, **and injected with brief current steps (three 2 nA injections, 2 ms each, at 50 Hz)**.

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
- **M1R modulation**: Testing whether M1 muscarinic receptor signaling affects dendritic calcium transients
- **Pharmacological approach**: M1R antagonist blocks muscarinic signaling to reveal receptor-dependent effects

**Analysis**:
- Ca2+ signals measured at proximal (~40 μm from soma) and distal (~90 μm from soma) dendritic locations
- **Distal/proximal ratio**: Used as index of dendritic excitability, normalized for experimental variables
- **Control vs. M1R antagonist**: Comparison reveals M1R-dependent component of dendritic excitability

**Key Finding**: M1R antagonist blocks specific components of iSPN dendritic excitability

### Spine Density Protocol

**Methodology**: Two-photon laser scanning microscopy with deconvolution

Same as Figure 2:
- **Pixel size**: 0.15 μm
- **Z-step size**: 0.3 μm
- **Analysis**: Semi-automated spine counting with NeuronStudio

### Pharmacological Manipulation

**M1R Antagonist Treatment**:
- **Purpose**: Test whether M1 muscarinic receptor signaling contributes to iSPN dendritic changes
- **Application**: Bath application during experiments
- **Comparison**: Control vs. M1R antagonist conditions

## Data Organization

### Somatic Excitability Data Structure

**Directory Structure**:
```
control/
├── 0217a/          # February 17, 2020, animal "a"
│   ├── cell1-001/  # -120 pA injection
│   ├── cell1-002/  # -100 pA injection
│   ...
│   └── cell1-021/  # +300 pA injection
├── 0217b/
├── 0218a/
...

M1R antagonist/
├── 0217a/
│   ├── cell1-001/
│   ├── cell1-002/
│   ...
│   └── cell1-021/
...
```

**File Bundle per Current Level**:
```
cell1-001/
├── cell1-001_Cycle00001_VoltageOutput_001.xml       # Stimulus protocol
├── cell1-001_Cycle00001_VoltageRecording_001.csv    # Voltage recording
├── cell1-001_Cycle00001_VoltageRecording_001.xml    # Recording metadata
├── cell1-001.xml                                    # Experiment metadata
└── References/                                      # Calibration files
```

### Dendritic Excitability Data Structure

**Raw Data Location**: `/home/heberto/development/surmeier-lab-to-nwb/link_to_raw_data/Figure 6/Dendritic excitability/`

**Naming Convention**: `[date]_Cell[cell number]_[location][location number]_trio-[trial number]`

**Complete Directory Structure**:
```
Figure 6/Dendritic excitability/
├── control/
│   ├── 0217a/          # February 17, 2020, animal "a"
│   │   ├── 02172020_Cell1_dist1_trio-001
│   │   ├── 02172020_Cell1_dist1_trio-002
│   │   ├── 02172020_Cell1_dist1_trio-003
│   │   ├── 02172020_Cell1_prox1_trio-001
│   │   ├── 02172020_Cell1_prox1_trio-002
│   │   └── 02172020_Cell1_prox1_trio-003
│   ├── 0217b/          # February 17, 2020, animal "b"
│   ├── 0218a/          # February 18, 2020, animal "a"
│   ├── 0218b/          # February 18, 2020, animal "b"
│   ├── 0219a/          # February 19, 2020, animal "a"
│   ├── 0219b/          # February 19, 2020, animal "b"
│   └── [Additional control experiments]
└── M1R antagonist/
    ├── 0217a/          # February 17, 2020, animal "a" + M1R antagonist
    │   ├── 02172020_Cell1_dist1_trio-001
    │   ├── 02172020_Cell1_dist1_trio-002
    │   └── 02172020_Cell1_dist1_trio-003
    ├── 0217b/          # February 17, 2020, animal "b" + M1R antagonist
    ├── 0218a/          # February 18, 2020, animal "a" + M1R antagonist
    ├── 0218b/          # February 18, 2020, animal "b" + M1R antagonist
    └── [Additional M1R antagonist experiments]
```

**Experimental Design**:
- **Paired experiments**: Same animals tested under control and M1R antagonist conditions
- **Within-subject design**: Allows direct comparison of M1R-dependent effects
- **Temporal sequence**: Control recordings often performed before drug application

**File Bundle per Recording**:
```
02172020_Cell1_dist1_trio-001/
├── *_Cycle00001_Ch1_000001.ome.tif          # Fluo-4 calcium channel kymograph
├── *-Cycle00001_Ch1Source.tif               # Field of view with scan line overlay
├── *_Cycle00001_Ch2_000001.ome.tif          # Alexa Fluor 568 structural channel
├── *-Cycle00001_Ch2Source.tif               # Structural channel field of view
├── *_Cycle00001_LineProfileData.csv         # Averaged fluorescence traces
├── *_Cycle00001_VoltageOutput_001.xml       # Stimulus protocol definition
├── *_Cycle00001_VoltageRecording_001.csv    # Raw electrophysiology data
├── *_Cycle00001_VoltageRecording_001.xml    # Recording metadata
├── *.env                                    # Environment configuration
├── *.xml                                    # Master experiment file
└── References/                              # Calibration files
```

### Spine Density Data Structure

**Directory Structure**:
```
control/
├── 02172020a/      # February 17, 2020, animal "a"
│   ├── Decon_20200217_Cell1_dist12/
│   ├── Decon_20200217_Cell1_prox1/
│   └── Decon_20200217_Cell1_prox2/
├── 02182020a/
...

M1R antagonist/
├── 02172020a/
│   ├── Decon_20200217_Cell1_dist1/
│   ├── Decon_20200217_Cell1_dist2/
│   └── Decon_20200217_Cell1_prox12/
...
```

**Z-Stack Files**:
```
Decon_20200217_Cell1_dist12/
├── 20_ZSeries-20200217_Cell1_dist12-001_Cycle00001_Ch1_#.ome_Z01.tif
├── 20_ZSeries-20200217_Cell1_dist12-001_Cycle00001_Ch1_#.ome_Z02.tif
├── 20_ZSeries-20200217_Cell1_dist12-001_Cycle00001_Ch1_#.ome_Z03.tif
...
├── 20_ZSeries-20200217_Cell1_dist12-001_Cycle00001_Ch1_#.ome_Z04.tif
```

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

## Key Findings

### M1R Role in iSPN Excitability
- **Control**: Baseline iSPN excitability and spine density
- **M1R antagonist**: Blocks M1 muscarinic receptor signaling
- **Effect**: Tests whether M1R contributes to dendritic excitability changes

### Mechanistic Insights
- M1R signaling pathway involvement in iSPN dendritic plasticity
- Relationship between muscarinic signaling and spine morphology
- Potential therapeutic target for LID treatment

## Experimental Design Notes

### Cell Type Identification
- iSPNs identified by lack of Drd1-Tdtomato expression (negative selection)
- Confirmed by electrophysiological properties and morphology

### Pharmacological Controls
- **Control**: Normal ACSF
- **M1R antagonist**: Specific M1 muscarinic receptor blocker
- **Timing**: Drug applied during recording session

### Comparison with Other Figures
- **Figure 3**: D2R antagonist (sulpiride) effects on iSPNs
- **Figure 6**: M1R antagonist effects on iSPNs
- **Complementary**: Tests different receptor systems in same cell type

## Related Files
- [Overview](../conversion_notes_overview.md)
- [Figure 1 Notes](figure_1_conversion_notes.md) - Similar experimental design
- [Figure 3 Notes](figure_3_conversion_notes.md) - D2R antagonist effects on iSPNs
- [Figure 7 Notes](figure_7_conversion_notes.md) - CDGI knockout effects
