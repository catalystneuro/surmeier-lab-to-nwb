# Figure 2 Conversion Notes

## Overview
Figure 2 examines dSPN spine morphology and synaptic strength changes between LID on- and off-states. The data includes:
- **Spine density**: Two-photon microscopy Z-stacks of dendritic segments
- **Sr2+-oEPSC**: Optogenetic stimulation with voltage clamp recordings measuring synaptic strength

## Data Structure

```
Figure 2/
├── Spine density/
│   ├── control dSPN/
│   ├── LID off-state dSPN/
│   ├── LID on-state dSPN/
│   └── PD dSPN/
└── Sr-oEPSC/
    ├── LID off-state/
    └── LID on-state/
```

## Experimental Protocols

### Spine Density Protocol

**Methodology**: Two-photon laser scanning microscopy with deconvolution

> For assessment of dendritic spine density, images of dendritic segments (proximal: ~40 μm from soma; distal: > 80 μm from soma) were acquired with 0.15 μm pixels with 0.3 μm z-steps. Images were deconvolved in AutoQuant X3.0.4 (MediaCybernetics, Rockville, MD) and semi-automated spine counting was performed using 3D reconstructions in NeuronStudio (CNIC, Mount Sinai School of Medicine, New York). On average, two proximal and two distal dendrites were imaged and analyzed per neuron.

**Acquisition Parameters**:
- **Pixel size**: 0.15 μm
- **Z-step size**: 0.3 μm
- **Analysis**: Semi-automated spine counting with NeuronStudio

### Sr2+-oEPSC Protocol

**Methodology**: Optogenetic stimulation with voltage clamp recording

> For Sr2+-oEPSC experiments, patch pipettes (3-4 MΩ resistance) were loaded with 120 CsMeSO3, 5 NaCl, 0.25 EGTA, 10 HEPES, 4 Mg-ATP, 0.3 Na-GTP, 10 TEA, 5 QX-314 (pH 7.25, osmolarity 280-290 mOsm/L). SPNs were held at -70 mV in the voltage-clamp configuration. After patching, recording solution was changed to Ca2+-free ACSF containing 3 mM SrCl2 and 10 μM gabazine (10 μM, to suppress GABAA-mediated currents). Slices were incubated with this Ca2+-free solution for 25 min before recording. EPSCs were evoked every 30 s by whole-field LED illumination (single 0.3-ms pulses).

**Key Parameters**:
- **Holding potential**: -70 mV
- **Recording solution**: Ca2+-free ACSF + 3 mM SrCl2 + 10 μM gabazine
- **Stimulation**: Blue LED, 0.3 ms pulses every 30 s
- **Analysis window**: 40-400 ms post-stimulation

## Data Organization

### Spine Density Data Structure

**Directory Structure**:
```
LID on-state dSPN/
├── 0411a2019/      # April 11, 2019, session "a"
├── 04122019a/      # April 12, 2019, session "a"
├── 05012019a/      # May 1, 2019, session "a"
├── 0706a/          # July 6, 2017, session "a"
├── 0706b/          # July 6, 2017, session "b"
├── 0707a/          # July 7, 2017, session "a"
...
```

**File Structure per Cell**:
```
0706a/
├── Decon_07062017_Cell1_dist1/     # Distal dendrite 1
├── Decon_07062017_Cell1_dist2/     # Distal dendrite 2
├── Decon_20170706_Cell1_prox12/    # Proximal dendrites 1&2
└── Decon_20170706_Cell1_prox2/     # Proximal dendrite 2
```

**Z-Stack Files**:
```
Decon_20170720_Cell1_dist1/
├── 20_ZSeries-20170720_Cell1_dist1-001_Cycle00001_Ch1_#.ome_Z001.tif
├── 20_ZSeries-20170720_Cell1_dist1-001_Cycle00001_Ch1_#.ome_Z002.tif
├── 20_ZSeries-20170720_Cell1_dist1-001_Cycle00001_Ch1_#.ome_Z003.tif
...
├── 20_ZSeries-20170720_Cell1_dist1-001_Cycle00001_Ch1_#.ome_Z013.tif
```

### Sr2+-oEPSC Data Structure

**Directory Structure**:
```
LID off-state/
├── 07052023a/      # July 5, 2023, animal "a"
│   ├── cell1_LED14-001/    # Cell 1, trial 1
│   ├── cell1_LED14-002/    # Cell 1, trial 2
│   ...
│   └── cell1_LED14-010/    # Cell 1, trial 10
├── 07052023b/      # July 5, 2023, animal "b"
│   ├── cell2_LED12-006/    # Cell 2, trial 6
│   ...
│   └── cell2_LED12-015/    # Cell 2, trial 15
```

**File Bundle per Trial**:
```
cell1_LED14-001/
├── cell1_LED14-001_Cycle00001_VoltageOutput_001.xml     # LED stimulus protocol
├── cell1_LED14-001_Cycle00001_VoltageRecording_001.csv  # Voltage clamp recording
├── cell1_LED14-001_Cycle00001_VoltageRecording_001.xml  # Recording metadata
├── cell1_LED14-001.env                                  # Environment settings
├── cell1_LED14-001.xml                                  # Experiment metadata
└── References/                                          # Calibration files
```

## Optogenetic Stimulation

### Viral Vector Injection
> To selectively activate corticostriatal synapses, an adeno-associated virus (AAV) carrying a channelrhodopsin 2 (ChR2) expression construct (AAV5-hSyn-hChR2(H134R)-EYFP) was injected into the motor cortex ipsilateral to the 6-OHDA lesion.

### LED Stimulation Parameters

From VoltageOutput XML files:
```json
"WaveformComponent_PulseTrain": {
    "Name": "Waveform Component #1",
    "PulseCount": "1",
    "PulseWidth": "0.3",
    "PulseSpacing": "49.7",
    "PulsePotentialStart": "5",
    "PulsePotentialDelta": "0",
    "RestPotential": "0",
    "FirstPulseDelay": "20",
    "Repetitions": "1",
    "DelayBetweenReps": "2000"
}
```

**Parameter Interpretation**:
- **PulseCount**: 1 → Single LED pulse
- **PulseWidth**: 0.3 ms → LED on for 0.3 milliseconds
- **PulsePotentialStart**: 5 V → LED pulse voltage
- **FirstPulseDelay**: 20 ms → Pulse begins 20 ms after trial start
- **DelayBetweenReps**: 2000 ms → 2 second interval between trials (matches "every 30 s" protocol)

### Device Information
- **LED Device**: Blue LED driver for ChR2 activation
- **Excitation wavelength**: ~470 nm (blue light for ChR2)
- **Stimulation site**: Corticostriatal terminals in dorsal striatum

## Analysis Methods

### Spine Density Analysis
1. **Image Acquisition**: Z-stacks with 0.15 μm pixels, 0.3 μm z-steps
2. **Deconvolution**: AutoQuant X3.0.4 software
3. **Segmentation**: NeuronStudio semi-automated 3D reconstruction
4. **Quantification**: Spine counts per μm of dendrite length

### Sr2+-oEPSC Analysis
> Amplitudes of Sr2+-oEPSCs were analyzed automatically using TaroTools Event Detection in Igor Pro 8 (WaveMetrics, Portland, Oregon) followed by visual verification. Events were measured between 40 ms to 400 ms after photostimulation. The threshold for detection of an event was greater than 5 SD above the noise.

**Analysis Parameters**:
- **Detection window**: 40-400 ms post-stimulation
- **Threshold**: >5 SD above noise
- **Software**: TaroTools Event Detection in Igor Pro 8
- **Verification**: Manual visual confirmation

## Detailed Sweep Timeline for Sr2+-oEPSC Experiments

Based on the paper and XML metadata, here's how each sweep was structured:

### Single Sweep Timeline Diagram

**Sweep Timeline (per 0.3-ms LED pulse):**
```
┌─────────────────────────────────────────────────────────────────┐
│ Voltage Recording (continuous throughout)                       │
└─────────────────────────────────────────────────────────────────┘

Time:    0ms    20ms   20.3ms              420ms        30,000ms
         │       │      │                   │              │
         │   ┌───┐│      │                   │              │
Recording│   │LED││      │                   │              │
Start    │   │ON ││      │                   │              │
         │   └───┘│      │                   │              │
         │       │      │                   │              │
         ▼       ▼      ▼                   ▼              ▼
    ┌─────────┬──┬────────────────────────┬─────────────────┐
    │Baseline │ Stim │  EPSC Detection    │   Wait Period   │
    │ (20ms)  │0.3ms │    Window          │  (until 30sec)  │
    │         │      │   (40-400ms        │                 │
    │         │      │   post-stimulus)   │                 │
    └─────────┴──┴────────────────────────┴─────────────────┘
```

### Key Timing Elements

**Before Stimulation (0-20ms):**
- Continuous voltage recording for baseline measurement
- Cell held at -70mV in voltage clamp
- 20ms of pre-stimulus baseline

**During Stimulation (20-20.3ms):**
- 0.3ms blue LED pulse activates ChR2 in cortical terminals
- Voltage recording continues throughout
- Light activates glutamate release

**After Stimulation (20.3-420ms):**
- Critical detection window: 40-400ms post-stimulus
- Strontium causes asynchronous glutamate release
- Individual EPSCs detected as discrete events
- Voltage recording captures all synaptic responses

**Inter-Sweep Interval:**
- Total cycle: 30 seconds between LED pulses
- Allows complete clearance of strontium effects
- Prevents synaptic depression from over-stimulation

### What Gets Analyzed

From the paper: "Events were measured between 40 ms to 400 ms after photostimulation"

So the voltage recording is continuous, but the analysis focuses on the 40-400ms window where strontium-mediated asynchronous EPSCs occur, allowing detection of individual synaptic events rather than summed responses.

## Key Findings

### Spine Morphology Changes
- **Off-state**: Reduced spine density in both proximal and distal dendrites
- **On-state**: Increased mushroom spine density (larger, mature spines)
- **Frequency vs. Amplitude**: Sr2+-oEPSC frequency unchanged, amplitude increased in on-state

### Synaptic Strength Correlation
- Structural changes (spine enlargement) correlate with functional changes (increased oEPSC amplitude)
- Suggests potentiation of existing synapses rather than formation of new synapses

## Related Files
- [Overview](../conversion_notes_overview.md)
- [Figure 1 Notes](figure_1_conversion_notes.md) - Related electrophysiology methods
- [Figure 4 Notes](figure_4_conversion_notes.md) - Additional spine density data
