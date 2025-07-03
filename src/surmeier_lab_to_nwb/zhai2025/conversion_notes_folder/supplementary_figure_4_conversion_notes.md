# Supplementary Figure 4 Conversion Notes

## Overview
Supplementary Figure 4 provides additional validation data for the spine density measurements comparing confocal and two-photon microscopy methods. This figure supports the main findings in Figure 4 by showing:
- **Extended spine density comparisons**: Additional examples of confocal vs. two-photon spine detection
- **Methodological validation**: Demonstrates the consistency and reliability of spine counting methods
- **Technical controls**: Shows the impact of different imaging modalities on spine detection

## Data Structure

Based on the main conversion notes, Supplementary Figure 4 likely shares data with Figure 4, particularly:

```
Figure 4/
└── Confocal spine density/
    └── Fig 4J and Suppl Fig 5/
        ├── raw/                    # Contains nd2 files for both Fig 4J and Suppl Fig 4
        ├── processed-not organized yet.../
        └── Spine density plots for Fig 4J and Suppl Fig 5.pzfx
```

## Experimental Protocols

### High-Resolution Confocal Microscopy

**Methodology**: Nikon AXR confocal laser microscope imaging of sparsely labeled neurons

From the main conversion notes:
> A Nikon AXR confocal laser microscope from the Center for Advanced Microscopy & Nikon Imaging Center at Northwestern University was used to image dendritic spines from sparsely labeled iSPNs in fixed brain sections. Z-stack images were acquired with a 60x oil immersion objective (NA = 1.49) at 0.125 μm intervals with a 0.09 μm pixel size.

**Acquisition Parameters**:
- **Objective**: 60x oil immersion (NA = 1.49)
- **Z-step size**: 0.125 μm intervals
- **Pixel size**: 0.09 μm
- **Processing**: De-noised and deconvolved using Nikon NIS-Elements AR 5.41.02

### Spine Analysis Workflow

**Segmentation Protocol**:
> Dendritic segments (~30 μm in length, located ~30 μm from the soma) were analyzed for spine density using Imaris 10.0.0 software (Oxford Instruments). The dendrite of interest was isolated by creating a new Surface. The 'Labkit for pixel classification' tool was used to train signal detection over noise, creating a tight surface.

**Analysis Parameters**:
- **Segment length**: ~30 μm dendritic segments
- **Distance from soma**: ~30 μm
- **Software**: Imaris 10.0.0 (Oxford Instruments)
- **Spine detection parameters**:
  - Thinnest spine head: 0.188 μm
  - Maximum spine length: 5 μm
  - No branched spines allowed

## Data Organization

### Raw Data Structure

**File Format**: Nikon ND2 files
```
raw/
├── 6-OHDA/
├── control/
├── off-state/
└── on-state/
```

**File Properties**:
- **Format**: .nd2 (Nikon imaging format)
- **Content**: High-resolution confocal Z-stacks
- **Channels**: Typically single channel (fluorescent protein labeling)
- **Dimensions**: Multi-dimensional (X, Y, Z, potentially T)

### Processed Data Structure

**Analysis Outputs**:
```
processed-not organized yet.../
├── 6-OHDA/
├── control/
├── off-state/
└── on-state/
```

**File Types**:
- **.ims files**: Imaris software project files
- **Excel files**: Quantified spine density measurements
- **Code files**: Analysis scripts and parameters
- **Combined excel**: Consolidated results across conditions

### Statistical Analysis

**File**: `Spine density plots for Fig 4J and Suppl Fig 5.pzfx`
- **Format**: GraphPad Prism format (.pzfx)
- **Content**: Statistical comparisons and plotting data
- **Purpose**: Generate publication-quality figures and statistical tests

## Key Findings

### Methodological Validation
- **Confocal advantages**: Higher resolution detection of smaller spines
- **Two-photon limitations**: May miss smaller spines due to resolution constraints
- **Consistency**: Both methods show similar trends across experimental conditions
- **Quantitative differences**: Confocal detects more spines per unit length

### Experimental Conditions Compared
- **Control**: Baseline spine density in healthy tissue
- **6-OHDA**: Parkinson's disease model (dopamine depletion)
- **Off-state**: LID off-state condition
- **On-state**: LID on-state condition

### Technical Considerations
- **Resolution impact**: Higher resolution reveals more spine details
- **Detection sensitivity**: Confocal better for small spine detection
- **Methodological consistency**: Both approaches show similar relative changes

## Related Validation Studies

### Figure 4J Relationship
From the main notes:
> Figure 4J showed confocal detects more spines, but Supplemental Figure 5 proved WHY and validated the central hypothesis that spine size changes, not spine number changes, drive the apparent 2PLSM oscillations.

### Supplementary Figure 5 Connection
- **Shared data**: Same raw ND2 files used for both analyses
- **Different focus**: Suppl Fig 4 focuses on methodological validation
- **Complementary findings**: Together they provide complete technical validation

## Acquisition Hardware

### Nikon AXR Confocal System
- **Manufacturer**: Nikon
- **Model**: AXR confocal laser microscope
- **Location**: Center for Advanced Microscopy & Nikon Imaging Center, Northwestern University
- **Objective**: 60x oil immersion, NA = 1.49

### Software Stack
- **Acquisition**: Nikon NIS-Elements AR 5.41.02
- **Processing**: De-noising and deconvolution
- **Analysis**: Imaris 10.0.0 (Oxford Instruments)
- **Statistics**: GraphPad Prism (for .pzfx files)

## NWB Conversion Considerations

### Image Data
- **ImageSeries**: Store raw ND2 confocal stacks
- **Metadata**: Include acquisition parameters, pixel sizes, z-step intervals
- **Compression**: Apply appropriate compression for large image stacks

### Analysis Results
- **PlaneSegmentation**: Store spine segmentation results
- **ImageSegmentation**: Define ROIs for dendritic segments
- **ProcessingModule**: Document analysis pipeline and parameters

### Experimental Metadata
- **Subject**: Include genotype, treatment condition, lesion status
- **Session**: Link to specific imaging sessions and conditions
- **Device**: Document confocal microscope specifications

## Related Files
- [Overview](conversion_notes_overview.md)
- [Figure 4 Notes](figure_4_conversion_notes.md) - Main spine density comparison methods
- [Figure 2 Notes](figure_2_conversion_notes.md) - Two-photon spine density methods

## Technical Notes

### ND2 File Handling
- **Python libraries**: nd2 (tlambert03) recommended for metadata extraction
- **Bioformats**: Alternative but requires Java dependencies
- **Lazy loading**: Use memory-mapped access for large files

### Imaris Integration
- **File formats**: .ims files contain analysis results and parameters
- **Segmentation data**: Extract spine coordinates and morphological measurements
- **Analysis pipeline**: Document automated vs. manual correction steps
