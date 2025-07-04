# Figure 7 Conversion Notes

## Overview
Figure 7 examines the role of CalDAG-GEFI (CDGI) in iSPN function using knockout mice across multiple experimental paradigms and includes behavioral assessment. The data includes:
- **Somatic excitability**: Current clamp electrophysiology in CDGI knockout vs wildtype mice
- **Dendritic excitability**: Combined patch clamp + two-photon line scan imaging with oxotremorine-M (oxoM) testing
- **Spine density**: Two-photon image stacks for morphological analysis
- **Behavioral assessment**: AIM scoring and contralateral rotation analysis
- **Genetic validation**: CDGI knockout vs. wild-type comparison across multiple domains

## Complete Data Structure

```
Figure 7/
├── AIM rating (Figure 7J)/
│   └── AIM testing_CDGI KO.xlsx                    # AIM scores behavioral data
├── contralateral rotations (Figure 7I)/
│   └── CDGI KO videos/                             # Video recordings for rotation analysis
├── KO DE on vs off/                                # Figure 7F - CDGI KO dendritic excitability
│   ├── KO off-state/                               # LID off-state condition
│   └── KO on-state/                                # LID on-state condition
├── KO SE on vs off/                                # Figure 7C - CDGI KO somatic excitability
│   ├── KO off-state/                               # LID off-state condition
│   └── KO on-state/                                # LID on-state condition
├── KO spine density/                               # Figure 7H - CDGI KO spine morphology
│   ├── KO off/                                     # LID off-state condition
│   └── KO on/                                      # LID on-state condition
└── oxoM on DE (Figure 7E)/                        # Oxotremorine-M effects on dendritic excitability
    ├── KO/                                         # CDGI KO oxo-M dendritic response
    └── WT/                                         # Wildtype oxo-M dendritic response
```

## Abbreviations and Terminology
- **KO**: Knockout (CDGI gene deletion)
- **WT**: Wildtype (control animals)
- **DE**: Dendritic Excitability
- **SE**: Somatic Excitability
- **oxoM**: Oxotremorine-M (muscarinic M1 receptor agonist for testing receptor function)
- **CDGI**: CalDAG-GEFI (Calcium and DAG-regulated Guanine nucleotide Exchange Factor I)
- **AIM**: Abnormal Involuntary Movement scale (dyskinesia severity measure)

## CDGI Background
From the paper:
> The M1R-mediated modulation of iSPN dendrites is dependent upon CalDAG-GEFI (CDGI), a striatum-enriched, Ca²⁺-activated guanine nucleotide exchange factor

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

**oxoM Dendritic Excitability Directory Structure**:
```
Figure 7/oxoM on DE/  # Oxotremorine-M effects on Dendritic Excitability
├── KO/                 # CDGI knockout mice
│   ├── 0301a/          # March 1, 2020, animal "a"
│   │   ├── 03012020_Cell1_dist1_trio-001    # BEFORE oxo-M (baseline)
│   │   ├── 03012020_Cell1_dist1_trio-002    # BEFORE oxo-M (baseline)
│   │   ├── 03012020_Cell1_dist1_trio-003    # BEFORE oxo-M (baseline)
│   │   ├── 03012020_Cell1_dist1_trio-004    # AFTER oxo-M application
│   │   ├── 03012020_Cell1_dist1_trio-005    # AFTER oxo-M application
│   │   ├── 03012020_Cell1_dist1_trio-006    # AFTER oxo-M application
│   │   ├── 03012020_Cell1_prox1_trio-001    # BEFORE oxo-M (baseline)
│   │   ├── 03012020_Cell1_prox1_trio-002    # BEFORE oxo-M (baseline)
│   │   ├── 03012020_Cell1_prox1_trio-003    # BEFORE oxo-M (baseline)
│   │   ├── 03012020_Cell1_prox1_trio-004    # AFTER oxo-M application
│   │   ├── 03012020_Cell1_prox1_trio-005    # AFTER oxo-M application
│   │   └── 03012020_Cell1_prox1_trio-006    # AFTER oxo-M application
│   ├── 0301b/          # March 1, 2020, animal "b"
│   ├── 0302a/          # March 2, 2020, animal "a"
│   ├── 0302b/          # March 2, 2020, animal "b"
│   ├── 0303a/          # March 3, 2020, animal "a"
│   └── [Additional CDGI KO experiments]
├── WT/                # Wild-type control mice
│   ├── 0301a/          # March 1, 2020, animal "a"
│   │   ├── 03012020_Cell1_dist1_trio-001    # BEFORE oxo-M (baseline)
│   │   ├── 03012020_Cell1_dist1_trio-002    # BEFORE oxo-M (baseline)
│   │   ├── 03012020_Cell1_dist1_trio-003    # BEFORE oxo-M (baseline)
│   │   ├── 03012020_Cell1_dist1_trio-004    # AFTER oxo-M application
│   │   ├── 03012020_Cell1_dist1_trio-005    # AFTER oxo-M application
│   │   ├── 03012020_Cell1_dist1_trio-006    # AFTER oxo-M application
│   │   ├── 03012020_Cell1_prox1_trio-001    # BEFORE oxo-M (baseline)
│   │   ├── 03012020_Cell1_prox1_trio-002    # BEFORE oxo-M (baseline)
│   │   ├── 03012020_Cell1_prox1_trio-003    # BEFORE oxo-M (baseline)
│   │   ├── 03012020_Cell1_prox1_trio-004    # AFTER oxo-M application
│   │   ├── 03012020_Cell1_prox1_trio-005    # AFTER oxo-M application
│   │   └── 03012020_Cell1_prox1_trio-006    # AFTER oxo-M application
│   ├── 0301b/          # March 1, 2020, animal "b"
│   ├── 0302a/          # March 2, 2020, animal "a"
│   ├── 0302b/          # March 2, 2020, animal "b"
│   └── [Additional WT control experiments]
└── ReadMe.rtf         # Additional experimental notes
```

## oxoM Protocol Convention
From `ReadMe.rtf`:
```
-001, -002, -003    # BEFORE oxo-M application (baseline measurements)
-004, -005, -006    # AFTER oxo-M application (drug response)
-007, -008, -009    # AFTER oxo-M application (additional recordings)
```

**Experimental Design**:
This follows a paired before/after design where:
* **Baseline recordings (-001 to -003)**: Measure dendritic excitability in normal ACSF
* **Drug application**: Add oxotremorine-M to test M1 receptor responsiveness
* **Post-drug recordings (-004 to -006/009)**: Measure dendritic excitability after oxo-M application
* **Comparison**: KO vs WT response to M1 receptor stimulation reveals CDGI-dependent component

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

## Experimental Paradigms

### 1. Somatic Excitability (Figure 7C)
**Location**: `Figure 7/KO SE on vs off/`
- **KO off-state/**: CDGI knockout mice during LID off-state
- **KO on-state/**: CDGI knockout mice during LID on-state
- **Protocol**: Current clamp electrophysiology with F-I curve measurements
- **Purpose**: Test CDGI role in somatic excitability changes during dyskinesia

### 2. Dendritic Excitability - LID States (Figure 7F)
**Location**: `Figure 7/KO DE on vs off/`
- **KO off-state/**: CDGI knockout mice during LID off-state
- **KO on-state/**: CDGI knockout mice during LID on-state
- **Protocol**: Line scan + electrophysiology during dyskinesia states
- **Purpose**: Test CDGI role in dendritic excitability during dyskinesia

### 3. Dendritic Excitability - oxoM Response (Figure 7E)
**Location**: `Figure 7/oxoM on DE/`
- **KO/**: CDGI knockout mice with oxoM before/after protocol
- **WT/**: Wildtype mice with oxoM before/after protocol
- **Protocol**: Paired before/after design testing M1 receptor responsiveness
- **Purpose**: Test whether CDGI is required for M1 receptor-mediated dendritic changes

### 4. Spine Density (Figure 7H)
**Location**: `Figure 7/KO spine density/`
- **KO off/**: CDGI knockout spine morphology during LID off-state
- **KO on/**: CDGI knockout spine morphology during LID on-state
- **Protocol**: Two-photon image stacks for morphological analysis
- **Purpose**: Test CDGI role in structural plasticity during dyskinesia

### 5. Behavioral Assessment (Figure 7I, 7J)
**Location**: `Figure 7/AIM rating/` and `Figure 7/contralateral rotations/`
- **AIM testing_CDGI KO.xlsx**: Abnormal Involuntary Movement scores
- **CDGI KO videos/**: Video recordings for rotation quantification
- **Purpose**: Assess behavioral phenotype of CDGI knockout mice

## Key Findings

### CDGI Role Across Multiple Domains
- **Somatic excitability**: CDGI affects iSPN intrinsic properties during dyskinesia
- **Dendritic excitability**: CDGI required for dendritic calcium responses during dyskinesia and M1 stimulation
- **Spine morphology**: CDGI influences structural plasticity during dyskinesia
- **Behavior**: CDGI knockout affects dyskinesia-related behaviors

### Mechanistic Insights
- **M1 receptor pathway**: CDGI is downstream effector of M1 receptor signaling in dendrites
- **State-dependent effects**: CDGI function varies between LID on/off states
- **Genetic validation**: Knockout approach confirms CDGI's causal role in iSPN plasticity

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
