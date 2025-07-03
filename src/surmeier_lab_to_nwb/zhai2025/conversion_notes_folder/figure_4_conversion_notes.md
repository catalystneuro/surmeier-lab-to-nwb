# Figure 4 Conversion Notes

## Overview
Figure 4 compares spine density measurements between two-photon and confocal microscopy methods, demonstrating that confocal microscopy detects more spines due to higher resolution. The data includes:
- **Confocal spine density**: High-resolution confocal images with Imaris analysis
- **Two-photon spine density**: Standard two-photon Z-stacks with NeuronStudio analysis
- **Sr2+-oEPSC**: Optogenetic stimulation recordings (similar to Figure 2)

## Data Structure

```
Figure 4/
├── Confocal spine density/
│   ├── Fig 4H/
│   ├── Fig 4I/
│   └── Fig 4J and Suppl Fig 5/
├── Spine density/
│   ├── control iSPN/
│   ├── LID off-state iSPN/
│   ├── LID on-state iSPN/
│   └── PD iSPN/
└── Sr-oEPSC/
    ├── LID off-state/
    └── LID on-state/
```

## Experimental Protocols

### High-Resolution Confocal Microscopy

> A Nikon AXR confocal laser microscope from the Center for Advanced Microscopy & Nikon Imaging Center at Northwestern University was used to image dendritic spines from sparsely labeled iSPNs in fixed brain sections. Z-stack images were acquired with a 60x oil immersion objective (NA = 1.49) at 0.125 μm intervals with a 0.09 μm pixel size. The images were de-noised and deconvolved using Nikon NIS-Elements AR 5.41.02 software.

**Acquisition Parameters**:
- **Objective**: 60x oil immersion (NA = 1.49)
- **Z-step size**: 0.125 μm
- **Pixel size**: 0.09 μm
- **Processing**: De-noised and deconvolved with Nikon NIS-Elements

### Confocal Analysis with Imaris

> Dendritic segments (~30 μm in length, located ~30 μm from the soma) were analyzed for spine density using Imaris 10.0.0 software (Oxford Instruments). The dendrite of interest was isolated by creating a new Surface. The 'Labkit for pixel classification' tool was used to train signal detection over noise, creating a tight surface. This surface was further refined to include only the dendritic segment of interest. A Filament for the dendritic segment was generated with the aid of an embedded supervised learning function. Spine detection was performed with the following parameters: thinnest spine head set at 0.188 μm, maximum spine length set at 5 μm, and no branched spines allowed. Spines were further refined using an embedded supervised learning function and manually corrected if necessary.

**Analysis Parameters**:
- **Segment length**: ~30 μm
- **Distance from soma**: ~30 μm
- **Minimum spine head**: 0.188 μm
- **Maximum spine length**: 5 μm
- **Branched spines**: Not allowed
- **Software**: Imaris 10.0.0 with Labkit pixel classification

### Two-Photon Spine Density (Standard Method)

Same as Figure 2:
- **Pixel size**: 0.15 μm
- **Z-step size**: 0.3 μm
- **Analysis**: Semi-automated spine counting with NeuronStudio

## Data Organization

### Confocal Spine Density Data Structure

#### Figure 4H
```
Fig 4H/
├── processed/      # Contains processed TIFF and JPG files
└── raw/           # Contains Olympus OIF files
    └── 5104-3 60x str_Cycle/
        ├── Image_01_01_01_01.oif
        ├── Image_01_01_01_01.oif.files/
        ├── Image_01_01_02_01.oif
        ├── Image_01_01_02_01.oif.files/
        ...
        ├── MATL_01_01.log
        └── TileConfiguration.txt
```

#### Figure 4I
```
Fig 4I/
├── processed/
│   ├── 6-OHDA/
│   ├── control/
│   ├── off-state/
│   └── on-state/
└── raw/
    ├── 6-OHDA/
    │   └── 8040-slide 1-slice 2-cell 1 proxi.nd2
    ├── control/
    │   └── 3824-slide 2-slice 2-cell 2 proxi.nd2
    ├── off-state/
    │   └── 8041-slide 2-slice 3-cell 2 proxi.nd2
    └── on-state/
        └── 8939-slide 1-slice 3-cell 2 den2 proxi.nd2
```

#### Figure 4J and Supplemental Figure 5
```
Fig 4J and Suppl Fig 5/
├── processed-not organized yet.../
├── raw/
│   ├── 6-OHDA/
│   ├── control/
│   ├── off-state/
│   └── on-state/
└── Spine density plots for Fig 4J and Suppl Fig 5.pzfx
```

### File Format Details

#### Olympus OIF Files (Figure 4H)
- **OIF files**: Olympus Image Files containing metadata
- **OIF.files folders**: Associated image slices and data
- **Contents per folder**:
  ```
  Image_01_01_01_01.oif.files/
  ├── s_238358488-1826303.roi    # Region of interest file
  ├── s_C001Z001.pty             # Parameter file for Z-slice 1
  ├── s_C001Z001.tif             # TIFF image for Z-slice 1
  ├── s_C001Z002.pty             # Parameter file for Z-slice 2
  ├── s_C001Z002.tif             # TIFF image for Z-slice 2
  ...
  ```

#### Nikon ND2 Files (Figures 4I, 4J)
- **ND2 files**: Nikon native format containing multi-dimensional image data
- **Metadata**: Embedded acquisition parameters, timestamps, channel information
- **Processing**: Can be read with nd2 Python library or bioformats

### Two-Photon Spine Density (Standard)

Same structure as Figure 2:
```
LID on-state iSPN/
├── 0411a2019/
│   └── Decon_20190411_Cell1_dist1/
├── 04122019a/
│   └── Decon_20190412_Cell1_prox12/
...
```

## Key Metadata

### Olympus OIF Metadata Example
```
[Acquisition Parameters Common]
AFAE Coefficient=1.0
Acquisition Device="FV10i"
AcquisitionMode=1
ImageCaputreDate='2024-09-28 12:53:12'
LaserTransmissivity01=3.5
LaserWavelength01=473
Number of Acquisition Channel=1
PinholeCoefficient_60X=2.0
PinholeDiameter=107000
```

### PTY File Parameters
```
[Acquisition Parameters Common]
Confocal="ON"
Magnification=60.0
ObjectiveLens NAValue=1.35
ObjectiveLens Name="UPLSAP60xO"
PMTDetectingMode="Analog"
PMTVoltage=599
PinholeDiameter=214000
```

## Analysis Methods

### Confocal vs. Two-Photon Comparison

**Resolution Differences**:
- **Confocal**: 0.09 μm pixels, 0.125 μm Z-steps
- **Two-photon**: 0.15 μm pixels, 0.3 μm Z-steps
- **Result**: Confocal detects ~2x more spines due to higher resolution

**Software Differences**:
- **Confocal**: Imaris with supervised learning and manual correction
- **Two-photon**: NeuronStudio semi-automated counting

### GraphPad Prism Analysis
- **File**: `Spine density plots for Fig 4J and Suppl Fig 5.pzfx`
- **Software**: GraphPad Prism statistical analysis
- **Content**: Statistical comparisons and plotting data

## Key Findings

### Methodological Validation
- **Higher resolution reveals more spines**: Confocal microscopy detects approximately twice as many spines as two-photon microscopy
- **Spine size changes**: Supplemental Figure 5 validates that spine size (not number) drives apparent changes in two-photon measurements
- **Technical consideration**: Resolution limits affect spine detection sensitivity

### Biological Implications
- Previous two-photon studies may underestimate absolute spine density
- Relative changes between conditions remain valid
- Importance of consistent methodology for comparative studies

## Acquisition Systems

### Nikon AXR Confocal System
- **Location**: Center for Advanced Microscopy & Nikon Imaging Center, Northwestern University
- **Objective**: 60x oil immersion (NA = 1.49)
- **Software**: Nikon NIS-Elements AR 5.41.02

### Olympus FV10i System
- **Confocal mode**: ON
- **Laser wavelength**: 473 nm
- **PMT detection**: Analog mode

## Related Files
- [Overview](../conversion_notes_overview.md)
- [Figure 2 Notes](figure_2_conversion_notes.md) - Standard two-photon spine density methods
- [Supplementary Figure 5](supplementary_figure_5_conversion_notes.md) - Spine size validation
