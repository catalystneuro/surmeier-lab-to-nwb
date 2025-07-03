# Figure 8 Conversion Notes

## Overview
Figure 8 examines the effects of M1R CRISPR deletion on iSPN properties and dyskinetic behaviors. The data includes:
- **M1R CRISPR somatic excitability**: Current injection protocols in M1R-deleted vs. control iSPNs
- **M1R CRISPR spine density**: Two-photon microscopy Z-stacks comparing spine morphology
- **AIM behavioral scoring**: Abnormal involuntary movement assessments in M1R CRISPR mice

## Data Structure

```
Figure 8/
├── M1R CRISPR AIMs/
│   └── AIM raw score_M1R CRISPR.xlsx
├── M1R CRISPR SE/
│   ├── interleaved control/
│   └── M1R CRISPR/
└── M1R CRISPR spine density/
    ├── control/
    └── M1R CRISPR/
```

## Experimental Protocols

### M1R CRISPR Gene Editing

**Methodology**: CRISPR-Cas9 gene editing to delete M1 muscarinic receptors specifically from iSPNs

**Experimental Groups**:
- **M1R CRISPR group**: AAV-Cas9 + AAV-gRNA-FusionRed → M1R deletion
- **Control group**: Saline + AAV-gRNA-FusionRed → No M1R deletion

From the paper's methods:
> Control mice were injected with only the gRNA-FR vector

### Somatic Excitability Protocol

**Methodology**: Whole-cell patch clamp recording in current clamp mode

Similar to Figures 1 and 3, voltage changes recorded in response to current injection steps of 500 ms duration.

**Current Steps**:
- Range: -120 pA to +300 pA
- Step size: 20 pA
- Duration: 500 ms each
- Total steps: 21 (from cell1-001 to cell1-021)

**Analysis**: F-I relationship curves and rheobase measurements

### Spine Density Protocol

**Methodology**: Two-photon laser scanning microscopy with deconvolution

Same as Figure 2:
- **Pixel size**: 0.15 μm
- **Z-step size**: 0.3 μm
- **Analysis**: Semi-automated spine counting with NeuronStudio

### AIM Behavioral Assessment

**Methodology**: Abnormal Involuntary Movement scoring

**Categories**:
- **Ax**: Axial AIMs (abnormal trunk/neck movements)
- **Li**: Limb AIMs (abnormal forelimb movements)
- **Ol**: Orolingual AIMs (abnormal mouth/tongue movements)
- **Total score**: Sum of all three categories

**Time Points**: 20min, 40min, 60min, 80min, 100min, 120min during levodopa sessions

## Data Organization

### M1R CRISPR Somatic Excitability Data Structure

**Directory Structure**:
```
interleaved control/
├── 20221012a/          # October 12, 2022, animal "a"
│   ├── cell1-001/      # -120 pA injection
│   ├── cell1-002/      # -100 pA injection
│   ...
│   └── cell1-021/      # +300 pA injection
├── 20221012b/
├── 20221012c/
├── 20221130a/
├── 20221130b/
└── 20221130c/

M1R CRISPR/
├── 20221004b/
│   ├── cell2-001/
│   ├── cell2-002/
│   ...
│   └── cell2-021/
├── 20221004c/
├── 20221005a/
├── 20221005c/
├── 20221005d/
├── 20221129a/
├── 20221129b/
└── 20221129c/
```

**File Bundle per Current Level**:
```
cell2-001/
├── cell2-001_Cycle00001_VoltageOutput_001.xml       # Stimulus protocol
├── cell2-001_Cycle00001_VoltageRecording_001.csv    # Voltage recording
├── cell2-001_Cycle00001_VoltageRecording_001.xml    # Recording metadata
├── cell2-001.env                                    # Environment config
├── cell2-001.xml                                    # Experiment metadata
└── References/                                      # Calibration files
```

### M1R CRISPR Spine Density Data Structure

**Directory Structure**:
```
control/
├── 20221012a/          # October 12, 2022, animal "a"
│   ├── Decon_20221012_cell1_dist1/
│   ├── Decon_20221012_cell1_prox1/
│   └── Decon_20221012_cell1_medium12/
├── 20221012b/
├── 20221012c/
├── 20221130a/
├── 20221130b/
└── 20221130c/

M1R CRISPR/
├── 20221004a/
│   ├── Decon_20221004_cell1_dist1/
│   ├── Decon_20221004_Cell1_medium12/
│   └── Decon_20221004_cell1_prox1/
├── 20221004b/
│   ├── Decon_20221004_cell2_dist1/
│   ├── Decon_20221004_cell2_dist2/
│   ├── Decon_20221004_cell2_prox1/
│   └── Decon_20221004_cell2_prox23/
...
```

**Z-Stack Files**:
```
Decon_20221004_cell1_dist1/
├── 20_ZSeries-20221004_cell1_dist1-001_Cycle00001_Ch1_#.ome_Z01.tif
├── 20_ZSeries-20221004_cell1_dist1-001_Cycle00001_Ch1_#.ome_Z02.tif
├── 20_ZSeries-20221004_cell1_dist1-001_Cycle00001_Ch1_#.ome_Z03.tif
...
├── 20_ZSeries-20221004_cell1_dist1-001_Cycle00001_Ch1_#.ome_Z07.tif
```

### AIM Scoring Data Structure

**File**: `AIM raw score_M1R CRISPR.xlsx`

**Structure**:
- **Columns B-G**: Time points (20min, 40min, 60min, 80min, 100min, 120min)
- **Column H**: Final total or summary score
- **Column I**: "60min rotations/min" - contralateral rotation counts

**Mouse Information**:
- Mouse IDs with weights (21g, 22g, 22.5g, 23g - typical adult mouse weights)
- Categories per mouse: Ax, Li, Ol, Total score

## Key Metadata

### Stimulus Parameters
Same as Figure 1 - current injection steps from -120 pA to +300 pA in 20 pA increments.

### Channel Assignments
- **Channel 1 (Ch1)**: Alexa Fluor 568 structural dye (red, 580-620 nm)

### Z-Stack Parameters
- **Pixel size**: 0.15 μm
- **Z-step size**: 0.3 μm
- **Typical stack**: 7-15 Z-planes per dendritic segment

## Acquisition Parameters

### Electrophysiology
> All the electrophysiological recordings were made using a MultiClamp 700B amplifier (Axon Instrument, USA), and signals were filtered at 2 kHz and digitized at 10 kHz. Voltage protocols and data acquisition were performed by PrairieView 5.3 (Bruker).

### Two-Photon Imaging
> The recorded SPN was visualized using 810 nm excitation laser (Chameleon Ultra II, Coherent, Santa Clara, USA). Dendritic structure was visualized by the red signal of Alexa Fluor 568 detected by a Hamamatsu R3982 side-on photomultiplier tube (PMT, 580-620 nm).

## Key Findings

### M1R CRISPR Effects on iSPN Excitability
- **Control**: Normal iSPN somatic excitability patterns
- **M1R CRISPR**: Altered excitability due to loss of M1R function in iSPNs
- **Mechanistic insight**: M1Rs are required for normal iSPN excitability

### M1R CRISPR Effects on Spine Density
- **Control**: Baseline spine density in off-state conditions
- **M1R CRISPR**: Changes in spine morphology due to M1R deletion
- **Comparison**: Effects specific to M1R loss vs. general genetic manipulation

### Behavioral Validation
- **AIM scores**: Reduced dyskinetic behaviors in M1R CRISPR mice
- **Therapeutic potential**: M1R deletion as potential treatment strategy
- **Cell-type specificity**: iSPN-specific M1R deletion effects

## Experimental Design Notes

### CRISPR-Cas9 Strategy
- **Target**: M1 muscarinic receptor gene deletion in iSPNs
- **Delivery**: AAV-mediated Cas9 and gRNA delivery
- **Controls**: gRNA-FusionRed only (no Cas9)
- **Validation**: Confirmed M1R protein loss

### Cell Type Identification
- iSPNs identified by lack of Drd1-Tdtomato expression (negative selection)
- Confirmed by electrophysiological properties and morphology
- Same identification criteria across experimental groups

### Experimental Controls
- **Interleaved controls**: Control and CRISPR mice processed together
- **Genetic background**: Same mouse strain and breeding
- **Standardized protocols**: Identical experimental conditions

## Related Files
- [Overview](conversion_notes_overview.md)
- [Figure 1 Notes](figure_1_conversion_notes.md) - Similar somatic excitability design
- [Figure 2 Notes](figure_2_conversion_notes.md) - Similar spine density methods
- [Figure 6 Notes](figure_6_conversion_notes.md) - M1R antagonist effects on iSPNs
- [Figure 7 Notes](figure_7_conversion_notes.md) - CDGI knockout behavioral studies
