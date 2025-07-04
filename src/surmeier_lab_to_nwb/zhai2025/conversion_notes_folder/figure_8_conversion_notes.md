# Figure 8 Conversion Notes

## Overview
Figure 8 examines the effects of M1R CRISPR deletion on iSPN properties and dyskinetic behaviors using targeted gene editing. The data includes:
- **M1R CRISPR somatic excitability**: Current injection protocols in M1R-deleted vs. control iSPNs
- **M1R CRISPR spine density**: Two-photon microscopy Z-stacks comparing spine morphology in off-state
- **AIM behavioral scoring**: Abnormal involuntary movement assessments in M1R CRISPR mice

## Complete Data Structure

```
Figure 8/
├── M1R CRISPR AIMs/
│   └── AIM raw score_M1R CRISPR.xlsx        # AIM scores like previous figures
├── M1R CRISPR SE/                           # Figure 8C-D - Somatic excitability
│   ├── interleaved control/                 # Control iSPNs somatic excitability data
│   └── M1R CRISPR/                          # M1R-deleted iSPNs somatic excitability data
└── M1R CRISPR spine density/                # Figure 8E-F - Spine morphology
    ├── control/                             # Control iSPN spine morphology - off-state
    └── M1R CRISPR/                          # M1R-deleted iSPN spine morphology - off-state
```

## CRISPR-Cas9 Gene Editing Approach
**M1R CRISPR** refers to the CRISPR-Cas9 gene editing technique used to delete M1 muscarinic receptors (M1Rs) specifically from iSPNs.

**Experimental Groups**:
- **M1R CRISPR group**: AAV-Cas9 + AAV-gRNA-FusionRed → M1R deletion
- **Control group**: Saline + AAV-gRNA-FusionRed → No M1R deletion

From the paper's methods:
> Control mice were injected with only the gRNA-FR vector

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
These are somatic excitability experiments with M1R CRISPR mice. The structure is like Figure 1 with multiple sessions across recording dates:

```
interleaved control/     # Control group (no M1R deletion)
├── 20221012a/          # October 12, 2022, animal "a"
├── 20221012b/          # October 12, 2022, animal "b"
├── 20221012c/          # October 12, 2022, animal "c"
├── 20221130a/          # November 30, 2022, animal "a"
├── 20221130b/          # November 30, 2022, animal "b"
└── 20221130c/          # November 30, 2022, animal "c"

M1R CRISPR/             # M1R-deleted group
├── 20221004b/          # October 4, 2022, animal "b"
├── 20221004c/          # October 4, 2022, animal "c"
├── 20221005a/          # October 5, 2022, animal "a"
├── 20221005c/          # October 5, 2022, animal "c"
├── 20221005d/          # October 5, 2022, animal "d"
├── 20221129a/          # November 29, 2022, animal "a"
├── 20221129b/          # November 29, 2022, animal "b"
└── 20221129c/          # November 29, 2022, animal "c"

Total: 17 directories (6 control + 8 M1R CRISPR sessions)
```

**Session Organization**:
Each session folder contains current injection steps organized by stimulus intensity:
```
20221004b/              # Single session example
├── cell2-001/          # -120 pA injection (step 1)
├── cell2-002/          # -100 pA injection (step 2)
├── cell2-003/          # -80 pA injection (step 3)
├── cell2-004/          # -60 pA injection (step 4)
├── cell2-005/          # -40 pA injection (step 5)
├── cell2-006/          # -20 pA injection (step 6)
├── cell2-007/          # 20 pA injection (step 7)
├── cell2-008/          # 40 pA injection (step 8)
├── cell2-009/          # 60 pA injection (step 9)
├── cell2-010/          # 80 pA injection (step 10)
├── cell2-011/          # 100 pA injection (step 11)
├── cell2-012/          # 120 pA injection (step 12)
├── cell2-013/          # 140 pA injection (step 13)
├── cell2-014/          # 160 pA injection (step 14)
├── cell2-015/          # 180 pA injection (step 15)
├── cell2-016/          # 200 pA injection (step 16)
├── cell2-017/          # 220 pA injection (step 17)
├── cell2-018/          # 240 pA injection (step 18)
├── cell2-019/          # 260 pA injection (step 19)
├── cell2-020/          # 280 pA injection (step 20)
└── cell2-021/          # 300 pA injection (step 21)
```

**File Bundle per Current Level**:
Each session folder (one per stimulus level) contains only patch clamp data:
```
cell2-001/
├── cell2-001_Cycle00001_VoltageOutput_001.xml       # Stimulus protocol definition
├── cell2-001_Cycle00001_VoltageRecording_001.csv    # Raw voltage recording data
├── cell2-001_Cycle00001_VoltageRecording_001.xml    # Recording metadata
├── cell2-001.env                                    # Environment configuration
├── cell2-001.xml                                    # Master experiment metadata
└── References/                                      # Calibration files
```

This structure contains sessions of patch clamp data, with every session divided to account for the stimuli like in Figure 1.

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
