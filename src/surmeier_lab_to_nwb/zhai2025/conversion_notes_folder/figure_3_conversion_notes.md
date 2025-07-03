# Figure 3 Conversion Notes

## Overview
Figure 3 examines iSPN (indirect pathway spiny projection neurons) excitability changes between LID on- and off-states, and the role of D2 receptor signaling. The data includes:
- **Somatic excitability**: Current injection protocols measuring rheobase and F-I relationships
- **Dendritic excitability**: Combined patch clamp + two-photon line scan imaging measuring calcium transients
- **Pharmacological manipulation**: Effects of sulpiride (D2R antagonist) on iSPN excitability

## Data Structure

```
Figure 3/
├── Dendritic excitability/
│   ├── LID off-state/
│   ├── LID on-state/
│   └── LID on-state with sul/
└── Somatic excitability/
    ├── LID off-state/
    ├── LID on-state/
    └── LID on-state with sul (iSPN)/
```

## Experimental Protocols

### Somatic Excitability Protocol

**Methodology**: Whole-cell patch clamp recording in current clamp mode

Similar to Figure 1, voltage changes recorded in response to current injection steps of 500 ms duration.

**Current Steps**:
- Range: -120 pA to +300 pA
- Step size: 20 pA
- Duration: 500 ms each
- Total steps: 21 (from cell1-001 to cell1-021)

**Analysis**: F-I relationship curves and rheobase measurements

### Dendritic Excitability Protocol

**Methodology**: Combination of patch clamp electrophysiology and two-photon laser scanning microscopy (2PLSM)

> To assess iSPN dendritic excitability in on- and off-states, and the role of D2 receptor signaling, the same methodology as Figure 1 was used. iSPNs in ex vivo brain slices were patch clamped, filled with Ca2+-sensitive dye Fluo-4 and Ca2+-insensitive dye Alexa Fluor 568, **and injected with brief current steps (three 2 nA injections, 2 ms each, at 50 Hz)**.

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
- **Calcium signaling**: bAPs trigger transient opening of voltage-dependent Ca2+ channels, producing measurable calcium transients
- **Dendritic excitability assessment**: Magnitude of Ca2+ signals serves as a surrogate estimate of dendritic depolarization extent

**Analysis**:
- Ca2+ signals measured at proximal (~40 μm from soma) and distal (~90 μm from soma) dendritic locations
- **Distal/proximal ratio**: Used as index of dendritic excitability, normalized for experimental variables
- **iSPN vs dSPN comparison**: Opposite effects observed - iSPNs show decreased excitability in on-state (vs. increased in dSPNs)

**Key Finding**: iSPN dendritic excitability decreased in on-state, reversed by sulpiride (D2R antagonist)

### Pharmacological Manipulation

**Sulpiride Treatment**:
- **Drug**: Sulpiride (selective D2 dopamine receptor antagonist)
- **Purpose**: Test whether observed changes in iSPN excitability depend on ongoing D2R signaling
- **Application**: Bath application during "LID on-state with sul" experiments

## Data Organization

### Somatic Excitability Data Structure

**Directory Structure**:
```
LID off-state/
├── 05232016_1/     # May 23, 2016, cell 1
├── 05232016_2/     # May 23, 2016, cell 2
├── 05242016_1/
├── 05242016_2/
...
```

Each of this is a session folder containing the recordings for a specific cell.

```
├── 05232016_1
│   ├── cell1-001
│   ├── cell1-002
│   ├── cell1-003
│   ├── cell1-004
│   ├── cell1-005
│   ├── cell1-006
│   ├── cell1-007
│   ├── cell1-008
│   ├── cell1-009
│   ├── cell1-010
│   ├── cell1-011
│   ├── cell1-012
│   ├── cell1-013
│   ├── cell1-014
│   ├── cell1-015
│   ├── cell1-016
│   ├── cell1-017
│   ├── cell1-018
│   ├── cell1-019
│   ├── cell1-020
│   └── cell1-021
├── 05232016_2
│   ├── cell2-001
│   ├── cell2-002
│   ├── cell2-003
│   ├── cell2-004
│   ├── cell2-005

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

**Raw Data Location**: `/home/heberto/development/surmeier-lab-to-nwb/link_to_raw_data/Figure 3/Dendritic excitability/`

**Naming Convention**: `[date]_Cell[cell number]_[location][location number]_trio-[trial number]`

**Complete Directory Structure**:
```
Figure 3/Dendritic excitability/
├── LID off-state/
│   ├── 0523a/          # May 23, 2016, animal "a"
│   │   ├── 05232016_Cell1_dist1_trio-001
│   │   ├── 05232016_Cell1_dist1_trio-002
│   │   ├── 05232016_Cell1_dist1_trio-003
│   │   ├── 05232016_Cell1_prox1_trio-001
│   │   ├── 05232016_Cell1_prox1_trio-002
│   │   └── 05232016_Cell1_prox1_trio-003
│   ├── 0523b/          # May 23, 2016, animal "b"
│   │   ├── 05232016_Cell2_dist1_trio-001
│   │   ├── 05232016_Cell2_dist1_trio-002
│   │   └── 05232016_Cell2_dist1_trio-003
│   ├── 0524a/          # May 24, 2016, animal "a"
│   ├── 0524b/          # May 24, 2016, animal "b"
│   ├── 0525a/          # May 25, 2016, animal "a"
│   └── [Additional off-state experiments]
├── LID on-state/
│   ├── 0523a/          # May 23, 2016, animal "a" (on-state)
│   ├── 0523b/          # May 23, 2016, animal "b" (on-state)
│   ├── 0524a/          # May 24, 2016, animal "a" (on-state)
│   ├── 0524b/          # May 24, 2016, animal "b" (on-state)
│   └── [Additional on-state experiments]
└── LID on-state with sul/  # "sul" = sulpiride (D2R antagonist)
    ├── 0523a/          # May 23, 2016, animal "a" + sulpiride
    ├── 0523b/          # May 23, 2016, animal "b" + sulpiride
    ├── 0524a/          # May 24, 2016, animal "a" + sulpiride
    └── [Additional sulpiride experiments]
```

**Experimental Design Notes**:
- **Cell type focus**: iSPNs (indirect pathway SPNs) vs. dSPNs in Figure 1
- **D2R modulation**: Sulpiride blocks D2 dopamine receptors to test receptor-dependent effects
- **Comparative design**: Same animals tested across different pharmacological conditions when possible

**File Bundle per Recording**:
```
05232016_Cell1_dist1_trio-001/
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

**Imaging Strategy**:
- **Dendritic targeting**: Line scans positioned along iSPN dendrites at specific distances from soma
- **Dual-channel acquisition**: Simultaneous recording of calcium (Fluo-4) and structural (Alexa Fluor 568) signals
- **Spatial calibration**: Distances measured from soma to scan locations using structural channel
- **Signal validation**: Structural channel confirms dendrite integrity throughout recording

## Key Findings

### iSPN Excitability Changes
- **Off-state**: Baseline iSPN excitability
- **On-state**: Decreased iSPN excitability (opposite to dSPN changes in Figure 1)
- **Sulpiride effect**: D2R antagonist blocks the on-state decrease in excitability

### Mechanistic Insights
- iSPN excitability changes are D2R-dependent
- Sulpiride application reverses the on-state phenotype
- Demonstrates opposing regulation of direct vs. indirect pathway neurons

## Experimental Design Notes

### Cell Type Identification
- iSPNs identified by lack of Drd1-Tdtomato expression (negative selection)
- Confirmed by electrophysiological properties and morphology

### Pharmacological Controls
- **Control**: Normal ACSF
- **Sulpiride**: D2R antagonist application
- **Timing**: Drug applied during recording session

## Related Files
- [Overview](../conversion_notes_overview.md)
- [Figure 1 Notes](figure_1_conversion_notes.md) - Similar experimental design for dSPNs
- [Figure 6 Notes](figure_6_conversion_notes.md) - M1R antagonist effects
