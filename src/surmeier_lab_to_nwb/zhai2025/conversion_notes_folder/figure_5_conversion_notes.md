# Figure 5 Conversion Notes

## Overview
Figure 5 examines acetylcholine (ACh) release dynamics using the genetically encoded fluorescent sensor GRABACh3.0. The data includes:
- **Two-photon imaging**: Time series of GRABACh3.0 fluorescence changes
- **Electrical stimulation**: Single pulse and burst stimulation protocols
- **Pharmacological manipulation**: Effects of dopamine, quinpirole (D2R agonist), and sulpiride (D2R antagonist)

## Data Structure

```
Figure 5/
├── LID off (Figure 5E)/
├── PD (6-OHDA lesioned Figure 5D)/
└── UL control (Unlesioned Figure 5C)/
```

## Experimental Protocols

### GRABACh3.0 Imaging Protocol

> ACh release was assessed by imaging GRABACh3.0, a genetically encoded fluorescent sensor of ACh, using 2PLSM. Acute slices with striatal expression of GRABACh3.0 were prepared as described above, transferred to a recording chamber, and continuously perfused with normal ACSF at 32–34°C. A two-photon laser (Chameleon Ultra II, Coherent, Santa Clara, CA) tuned to 920 nm was used to excite GRABACh3.0. Fluorescence was imaged using an Ultima In Vitro Multiphoton Microscope system (Bruker, Billerica, MA) with an Olympus 60x/0.9 NA water-immersion objective lens and a Hamamatsu H7422P-40 GaAsP PMT (490 nm to 560 nm, Hamamatsu Photonics, Hamamatsu, Japan).

**Acquisition Parameters**:
- **Excitation**: 920 nm (Chameleon Ultra II laser)
- **Objective**: Olympus 60x/0.9 NA water-immersion
- **Detection**: Hamamatsu H7422P-40 GaAsP PMT (490-560 nm)
- **Pixel size**: 0.388 μm × 0.388 μm
- **Dwell time**: 8 μs
- **Frame rate**: 21.26 fps

### Stimulation Protocols

> Time series images of the GRABACh3.0 were acquired with 0.388 μm × 0.388 μm pixels, 8-μs dwell time and a frame rate of 21.26 fps. After 3-s baseline acquisition, synchronous ACh release was evoked by delivering a single (1 ms x 0.3 mA) or a train of 20 electrical stimuli (1 ms x 0.3 mA at 20 Hz) by a concentric bipolar electrode (CBAPD75, FHC) placed at 200 μm ventral to the region of interest. Imaging was continued for at least another 5 s. Two trials were performed for each stimulation protocol and data averaged.

**Stimulation Parameters**:
- **Single pulse**: 1 ms × 0.3 mA
- **Burst stimulation**: 20 pulses at 20 Hz (1 ms × 0.3 mA each)
- **Electrode**: Concentric bipolar electrode (CBAPD75, FHC)
- **Placement**: 200 μm ventral to imaging region
- **Baseline**: 3 s before stimulation
- **Post-stimulation**: ≥5 s imaging
- **Trials**: 2 per protocol, averaged

### Pharmacological Treatments

> The slices were imaged for (in this sequence) control, 50 nM DA (only for lesioned mice), 10 μM quinpirole and 10 μM sulpiride treatment conditions, with at least 5 min of perfusion for each treatment.

**Treatment Sequence**:
1. **Control**: Baseline ACSF
2. **50 nM DA**: Dopamine (only for 6-OHDA lesioned mice)
3. **10 μM quinpirole**: D2R agonist
4. **10 μM sulpiride**: D2R antagonist

**Perfusion**: ≥5 min between treatments

### Calibration Protocol

> The minimal (Fmin) and maximal fluorescence intensity (Fmax) were determined by applying 10 μM TTX (to block any basal transmission) and 100 μM acetylcholine chloride (to saturate GRABACh3.0 signal), respectively.

**Calibration Conditions**:
- **Fmin**: 10 μM TTX (blocks sodium channels)
- **Fmax**: 100 μM acetylcholine chloride (saturates sensor)

## Data Organization

### Directory Structure

```
LID off/
├── 05242024slice1/
├── 05242024slice2/
├── 05292024slice1/
├── 05292024slice2/
├── 05312024slice1/
├── 05312024slice2/
├── 06032024slice1/
├── 06042024slice1/
├── 06042024slice2/
├── 06042024slice3/
└── 06042024slice4/

PD/
├── 04122024slice1ROI1/
├── 04122024slice1ROI2/
├── 04122024slice2ROI1/
├── 04122024slice2ROI2/
├── 04162024slice1ROI1/
├── 04162024slice1ROI2/
├── 04162024slice2ROI1/
├── 04162024slice2ROI2/
├── 05092024slice1/
└── 05092024slice2/

UL control/
├── 04022024slice1ROI1/
├── 04022024slice1ROI2/
├── 04022024slice2ROI1/
├── 04022024slice2ROI2/
├── 04022024slice3ROI1/
├── 04022024slice3ROI2/
├── 04052024slice1ROI1/
├── 04052024slice1ROI2/
├── 04052024slice2ROI1/
├── 04102024slice1ROI1/
├── 04102024slice1ROI2/
└── 04102024slice2ROI1/
```

### Experimental Naming Convention

**Format**: `BOT_[date]_slice[#]ROI[#]_[treatment]_[stimulation]-[trial]`

**Example**: `BOT_04162024_slice2ROI1_50nMDA_burst-001`

**Treatment Abbreviations**:
- **ctr**: Control (baseline ACSF)
- **50nMDA**: 50 nM dopamine
- **quin**: 10 μM quinpirole (D2R agonist)
- **sul**: 10 μM sulpiride (D2R antagonist)
- **ACh**: 100 μM acetylcholine (calibration)
- **TTX**: 10 μM tetrodotoxin (calibration)

**Stimulation Types**:
- **single**: Single pulse stimulation
- **burst**: Burst stimulation (20 pulses at 20 Hz)

### File Bundle per Experiment

```
BOT_04162024_slice2ROI1_sul_single-001/
├── BOT_04162024_slice2ROI1_sul_single-001_Cycle00001-botData.csv
├── BOT_04162024_slice2ROI1_sul_single-001_Cycle00001_Ch2_000001.ome.tif
├── BOT_04162024_slice2ROI1_sul_single-001_Cycle00001_Ch2_000002.ome.tif
├── BOT_04162024_slice2ROI1_sul_single-001_Cycle00001_Ch2_000003.ome.tif
...
├── BOT_04162024_slice2ROI1_sul_single-001_Cycle00001_Ch2_000121.ome.tif
├── BOT_04162024_slice2ROI1_sul_single-001_Cycle00001_Ch3_000001.ome.tif
├── BOT_04162024_slice2ROI1_sul_single-001_Cycle00001_Ch3_000002.ome.tif
├── BOT_04162024_slice2ROI1_sul_single-001_Cycle00001_Ch3_000003.ome.tif
...
├── BOT_04162024_slice2ROI1_sul_single-001_Cycle00001_VoltageOutput_001.xml
├── BOT_04162024_slice2ROI1_sul_single-001.env
└── BOT_04162024_slice2ROI1_sul_single-001.xml
```

## Key Metadata

### Channel Assignments
- **Channel 2 (Ch2)**: GRABACh3.0 fluorescence (PMT gain 900)
- **Channel 3 (Ch3)**: Dodt contrast/transmitted light (PMT gain 320)

### BOT Data Analysis

> The whole image was the region of interest (ROI) used for analysis. Fluorescent intensity data were analyzed by custom Python code (accessible upon request). Briefly, the fluorescence intensity values were first background-subtracted (the background resulted from PMT was measured by imaging with same PMT voltage but zero laser power). Baseline fluorescence F0 was the average fluorescence over the 1 s period right before stimulation. ΔF = F − F0 was normalized by (Fmax − Fmin) and then analyzed.

**Analysis Steps**:
1. **Background subtraction**: PMT dark current removal
2. **Baseline calculation**: F0 = average 1 s before stimulation
3. **ΔF calculation**: ΔF = F - F0
4. **Normalization**: ΔF/(Fmax - Fmin)

### Acquisition Timeline
- **Baseline**: 3 s pre-stimulation
- **Stimulation**: At t = 3 s
- **Post-stimulation**: ≥5 s recording
- **Total duration**: ≥8 s per trial

## Acquisition Parameters

### Two-Photon System
- **Microscope**: Ultima In Vitro Multiphoton Microscope (Bruker)
- **Laser**: Chameleon Ultra II (Coherent)
- **Wavelength**: 920 nm
- **Objective**: Olympus 60x/0.9 NA water-immersion

### Detection
- **PMT**: Hamamatsu1 H7422P-40 GaAsP
- **Detection range**: 490-560 nm
- **Gain settings**: Ch2 = 900, Ch3 = 320

### Stimulation Equipment
- **Electrode**: Concentric bipolar electrode (CBAPD75, FHC)
- **Stimulator**: Electrical stimulator for precise timing
- **Parameters**: 1 ms pulses, 0.3 mA amplitude

## Key Findings

### ACh Release Dynamics
- **Control**: Baseline ACh release patterns
- **6-OHDA lesion**: Altered ACh dynamics in PD model
- **LID off-state**: Specific ACh release profile
- **Dopamine effects**: D2R-mediated modulation of ACh release

### Pharmacological Responses
- **Quinpirole**: D2R activation reduces ACh release
- **Sulpiride**: D2R blockade increases ACh release
- **State-dependent**: Different responses in off vs. on states

## Related Files
- [Overview](../conversion_notes_overview.md)
- [Figure 1 Notes](figure_1_conversion_notes.md) - Two-photon imaging methods
- [Figure 2 Notes](figure_2_conversion_notes.md) - Electrical stimulation protocols
