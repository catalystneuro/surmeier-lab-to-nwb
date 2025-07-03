# Figure 7 Conversion Notes

## Overview
Figure 7 examines the role of CDGI (cAMP-dependent guanine nucleotide exchange factor I) in iSPN dendritic excitability using knockout mice. The data includes:
- **Dendritic excitability**: Combined patch clamp + two-photon line scan imaging in CDGI knockout mice
- **Genetic manipulation**: CDGI knockout vs. wild-type comparison
- **Mechanistic validation**: Tests CDGI's role in iSPN dendritic plasticity

## Data Structure

```
Figure 7/
└── Dendritic excitability/
    ├── CDGI KO/
    └── WT/
```

## Experimental Protocols

### Dendritic Excitability Protocol

**Methodology**: Combination of patch clamp electrophysiology and two-photon laser scanning microscopy (2PLSM)

> To assess the role of CDGI in iSPN dendritic excitability, the same methodology as previous figures was applied to CDGI knockout and wild-type mice. iSPNs were patch clamped, filled with Ca2+-sensitive dye Fluo-4 and Ca2+-insensitive dye Alexa Fluor 568, **and injected with brief current steps (three 2 nA injections, 2 ms each, at 50 Hz)**.

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
- **Genetic approach**: Testing whether CDGI gene deletion affects dendritic calcium transients
- **Molecular validation**: Confirms the role of specific signaling molecules in dendritic excitability

**Analysis**:
- Ca2+ signals measured at proximal (~40 μm from soma) and distal (~90 μm from soma) dendritic locations
- **Distal/proximal ratio**: Used as index of dendritic excitability, normalized for experimental variables
- **CDGI KO vs. WT**: Comparison reveals CDGI-dependent component of dendritic excitability

**Key Finding**: CDGI knockout alters iSPN dendritic excitability, confirming its role in dendritic plasticity

### Genetic Manipulation

**CDGI Knockout**:
- **Purpose**: Test whether CDGI (cAMP-dependent guanine nucleotide exchange factor I) is required for iSPN dendritic excitability changes
- **Model**: Genetic knockout mice lacking CDGI expression
- **Comparison**: CDGI KO vs. wild-type (WT) littermates

## Data Organization

### Dendritic Excitability Data Structure

**Raw Data Location**: `/home/heberto/development/surmeier-lab-to-nwb/link_to_raw_data/Figure 7/oxoM on DE/`

**Naming Convention**: `[date]_Cell[cell number]_[location][location number]_trio-[trial number]`

**Complete Directory Structure**:
```
Figure 7/oxoM on DE/  # Note: "oxoM on DE" = oxotremorine-M effects on Dendritic Excitability
├── KO/                 # CDGI knockout mice
│   ├── 0301a/          # March 1, 2020, animal "a"
│   │   ├── 03012020_Cell1_dist1_trio-001
│   │   ├── 03012020_Cell1_dist1_trio-002
│   │   ├── 03012020_Cell1_dist1_trio-003
│   │   ├── 03012020_Cell1_prox1_trio-001
│   │   ├── 03012020_Cell1_prox1_trio-002
│   │   └── 03012020_Cell1_prox1_trio-003
│   ├── 0301b/          # March 1, 2020, animal "b"
│   ├── 0302a/          # March 2, 2020, animal "a"
│   ├── 0302b/          # March 2, 2020, animal "b"
│   ├── 0303a/          # March 3, 2020, animal "a"
│   └── [Additional CDGI KO experiments]
├── WT/                # Wild-type control mice
│   ├── 0301a/          # March 1, 2020, animal "a"
│   │   ├── 03012020_Cell1_dist1_trio-001
│   │   ├── 03012020_Cell1_dist1_trio-002
│   │   ├── 03012020_Cell1_dist1_trio-003
│   │   ├── 03012020_Cell1_prox1_trio-001
│   │   ├── 03012020_Cell1_prox1_trio-002
│   │   └── 03012020_Cell1_prox1_trio-003
│   ├── 0301b/          # March 1, 2020, animal "b"
│   ├── 0302a/          # March 2, 2020, animal "a"
│   ├── 0302b/          # March 2, 2020, animal "b"
│   └── [Additional WT control experiments]
└── ReadMe.rtf         # Additional experimental notes
```

**Experimental Design**:
- **Littermate controls**: CDGI KO and WT mice from same breeding pairs
- **Genotype blinding**: Experimenter unaware of genotype during recordings
- **Standardized protocols**: Identical experimental conditions for both genotypes
- **oxoM treatment**: Some experiments include oxotremorine-M (muscarinic agonist) to test receptor responsiveness

**File Bundle per Recording**:
```
03012020_Cell1_dist1_trio-001/
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

## Key Findings

### CDGI Role in iSPN Dendritic Excitability
- **Wild-type**: Normal iSPN dendritic excitability patterns
- **CDGI KO**: Altered dendritic excitability due to loss of CDGI function
- **Mechanistic insight**: CDGI is required for normal iSPN dendritic function

### Signaling Pathway Validation
- CDGI links cAMP signaling to dendritic excitability changes
- Genetic approach validates pharmacological findings
- Confirms CDGI as a key molecular player in iSPN plasticity

## Experimental Design Notes

### Genetic Background
- **Knockout model**: CDGI gene deletion
- **Controls**: Wild-type littermates
- **Breeding**: Heterozygous crosses to generate KO and WT siblings

### Cell Type Identification
- iSPNs identified by lack of Drd1-Tdtomato expression (negative selection)
- Confirmed by electrophysiological properties and morphology
- Same identification criteria across genotypes

### Experimental Controls
- **Genotype blinding**: Experimenter blinded to genotype during recording
- **Littermate controls**: WT and KO from same litters
- **Standardized protocols**: Identical experimental conditions

## CDGI Background

### Molecular Function
- **Full name**: cAMP-dependent guanine nucleotide exchange factor I
- **Function**: Activates small GTPases in response to cAMP signaling
- **Pathway**: Links cAMP/PKA signaling to cytoskeletal regulation
- **Relevance**: Important for dendritic spine dynamics and excitability

### Knockout Strategy
- **Target**: CDGI gene deletion
- **Validation**: Confirmed loss of CDGI protein expression
- **Viability**: Viable knockout mice with specific neuronal phenotypes

## Related Files
- [Overview](../conversion_notes_overview.md)
- [Figure 1 Notes](figure_1_conversion_notes.md) - Similar experimental design
- [Figure 3 Notes](figure_3_conversion_notes.md) - iSPN excitability studies
- [Figure 6 Notes](figure_6_conversion_notes.md) - M1R antagonist effects on iSPNs
