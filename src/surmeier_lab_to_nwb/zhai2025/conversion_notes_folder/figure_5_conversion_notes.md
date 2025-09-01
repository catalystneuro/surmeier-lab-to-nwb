# Figure 5 Conversion Notes

## Overview
Figure 5 examines acetylcholine (ACh) release dynamics using the genetically encoded fluorescent sensor GRABACh3.0 in a mouse model of Parkinson's disease and levodopa-induced dyskinesia (LID). This study reveals state-dependent modulation of cholinergic signaling and its role in LID pathophysiology.

### Key Findings
- **ACh release elevation**: Dramatically increased in parkinsonian and LID off-states due to loss of dopaminergic inhibition
- **D2R-mediated control**: Dopamine continues to robustly inhibit ACh release in dyskinetic mice
- **State-dependent oscillations**: ACh levels oscillate counter to dopamine levels between on- and off-states
- **Therapeutic implications**: Disrupting M1R signaling in iSPNs attenuates LID while enhancing levodopa benefits

## Experimental Context
This figure demonstrates that cholinergic interneurons (ChIs) play a critical role in LID pathophysiology, particularly during the off-state when striatal dopamine levels are low and ACh release is disinhibited. The study challenges previous assumptions about cholinergic involvement in dyskinesia by revealing state-specific changes.

## Data Structure

```
Figure 5/
├── LID off (Figure 5E) - Dyskinetic mice 24-48h after last levodopa dose
├── PD (6-OHDA lesioned Figure 5D) - Parkinsonian mice without levodopa treatment
└── UL control (Unlesioned Figure 5C) - Healthy control mice
```

### Subject Groups
- **UL control**: Unlesioned control mice (healthy striatum)
- **PD**: 6-OHDA lesioned mice (parkinsonian state, no levodopa)
- **LID off**: Dyskinetic mice in off-state (24-48h post-levodopa)

### ROI Designation
Some experiments include ROI1 and ROI2 designations, representing different fields of view within the same slice to increase sampling or find optimal GRABACh3.0 expression areas. This is relevant for subject matching during analysis.

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
- **Temperature**: 32-34°C
- **Perfusion**: Continuous ACSF flow

### Sensor Properties
**GRABACh3.0** is an optimized genetically encoded acetylcholine sensor that provides:
- High sensitivity and selectivity for ACh detection
- Fast kinetics suitable for monitoring rapid neurotransmitter release
- Minimal interference with endogenous cholinergic signaling
- Robust signal-to-noise ratio for quantitative measurements

### Stimulation Protocols

> Time series images of the GRABACh3.0 were acquired with 0.388 μm × 0.388 μm pixels, 8-μs dwell time and a frame rate of 21.26 fps. After 3-s baseline acquisition, synchronous ACh release was evoked by delivering a single (1 ms x 0.3 mA) or a train of 20 electrical stimuli (1 ms x 0.3 mA at 20 Hz) by a concentric bipolar electrode (CBAPD75, FHC) placed at 200 μm ventral to the region of interest. Imaging was continued for at least another 5 s. Two trials were performed for each stimulation protocol and data averaged.

**Stimulation Parameters**:
- **Single pulse**: 1 ms × 0.3 mA (single electrical stimulus)
- **Burst stimulation**: 20 pulses at 20 Hz (1 ms × 0.3 mA each pulse)
- **Electrode**: Concentric bipolar electrode (CBAPD75, FHC)
- **Placement**: 200 μm ventral to imaging region
- **Baseline**: 3 s pre-stimulation recording
- **Post-stimulation**: ≥5 s continued imaging
- **Trials**: 2 per protocol, data averaged
- **Stimulation timing**: Delivered at t = 3 s after baseline start

**Stimulation Rationale**:
- **Single pulse**: Mimics physiological ACh release from single action potentials
- **Burst stimulation**: Simulates high-frequency firing patterns of cholinergic interneurons
- Both protocols reveal different aspects of ACh release dynamics and pharmacological modulation

**Important Note on Metadata Discrepancies**:
XML metadata in the raw data files shows different stimulation parameters than those reported in the paper:
- **XML metadata**: LED stimulator with 5V pulses, 0.1ms pulse width, ~204Hz burst frequency
- **Paper specifications**: Concentric bipolar electrode, 0.3mA amplitude, 1ms pulse width, 20Hz burst frequency
- **Conversion decision**: Paper specifications are used in NWB conversion as they represent the actual experimental parameters used

### Pharmacological Treatments

> The slices were imaged for (in this sequence) control, 50 nM DA (only for lesioned mice), 10 μM quinpirole and 10 μM sulpiride treatment conditions, with at least 5 min of perfusion for each treatment.

**Treatment Sequence**:
1. **Control**: Baseline ACSF (artificial cerebrospinal fluid)
2. **50 nM DA**: Dopamine (only for 6-OHDA lesioned mice)
3. **10 μM quinpirole**: D2R agonist
4. **10 μM sulpiride**: D2R antagonist

**Perfusion Protocol**: ≥5 min equilibration between treatments

**Pharmacological Rationale**:
- **50 nM DA**: Mimics on-state dopamine levels in vivo (73-75 nM detected in striatum during levodopa treatment)
- **Quinpirole**: D2 receptor agonist that activates D2Rs on cholinergic interneurons, leading to ACh release suppression
- **Sulpiride**: D2 receptor antagonist that blocks dopamine's inhibitory effect on ACh release, revealing baseline dopaminergic tone

**Drug Effects**:
- **+DA**: Mimics natural ligand effect—activating dopamine receptors (primarily D2Rs on ChIs) and inhibiting ACh release
- **+Quinpirole**: Activates D2 receptors, suppressing ACh release even in absence of endogenous dopamine
- **+Sulpiride**: Blocks D2 receptors, preventing dopamine inhibition and resulting in increased evoked ACh signal

### Calibration Protocol

> The minimal (Fmin) and maximal fluorescence intensity (Fmax) were determined by applying 10 μM TTX (to block any basal transmission) and 100 μM acetylcholine chloride (to saturate GRABACh3.0 signal), respectively.

**Calibration Conditions**:
- **Fmin**: 10 μM TTX (tetrodotoxin - blocks sodium channels, eliminates neural activity)
- **Fmax**: 100 μM acetylcholine chloride (saturates GRABACh3.0 sensor response)

**Calibration Purpose**:
- **TTX application**: Determines minimal fluorescence by blocking all neural transmission and endogenous ACh release
- **ACh saturation**: Establishes maximal sensor response for normalization of evoked signals
- **Normalization formula**: ΔF/(Fmax - Fmin) provides standardized measurement across experiments

## Data Organization

**Raw Data Location**: `/home/heberto/development/surmeier-lab-to-nwb/link_to_raw_data/Figure 5_SF2/`

### Directory Structure

```
Figure 5_SF2/ (Raw data directory)
├── LID off/ (Figure 5E - Dyskinetic mice, off-state)
│   ├── 05242024slice1/
│   ├── 05242024slice2/
│   ├── 05292024slice1/
│   ├── 05292024slice2/
│   ├── 05312024slice1/
│   ├── 05312024slice2/
│   ├── 06032024slice1/
│   ├── 06042024slice1/
│   ├── 06042024slice2/
│   ├── 06042024slice3/
│   └── 06042024slice4/
├── PD/ (Figure 5D - 6-OHDA lesioned, parkinsonian)
│   ├── 04122024slice1ROI1/
│   ├── 04122024slice1ROI2/
│   ├── 04122024slice2ROI1/
│   ├── 04122024slice2ROI2/
│   ├── 04162024slice1ROI1/
│   ├── 04162024slice1ROI2/
│   ├── 04162024slice2ROI1/
│   ├── 04162024slice2ROI2/
│   ├── 05092024slice1/
│   └── 05092024slice2/
└── UL control/ (Figure 5C - Unlesioned control)
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
    └── 04102024slice2ROI2/
```

### Naming Convention
- **Date format**: MMDDYYYY (e.g., 04162024 = April 16, 2024)
- **Slice numbering**: Sequential numbering within each experimental day
- **ROI designation**: Region of interest within slice (ROI1, ROI2) for multiple imaging fields
- **LID off naming**: Uses "slice1A" format instead of "slice1ROI1" (e.g., `05242024slice1A`)
- **PD and UL control**: Use standard "sliceXROIY" format for multiple fields per slice
- **Experimental variation**: Some slices have single ROI, others have multiple ROIs per slice

### Key Structural Differences Observed
- **LID off group**: Uses "sliceXA" naming convention instead of "sliceXROIY"
- **Missing TTX in some experiments**: Not all slice directories contain TTX calibration trials
- **Variable ACh trials**: Some LID off experiments have 2 ACh calibration trials instead of 1
- **Consistent PD/UL structure**: PD and UL control groups follow standard ROI1/ROI2 naming

### Experimental Naming Convention

**Format**: `BOT_[date]_slice[#]ROI[#]_[treatment]_[stimulation]-[trial]`

**Example**: `BOT_04162024_slice2ROI1_50nMDA_burst-001`

**Treatment Abbreviations**:
- **ctr**: Control (baseline ACSF)
- **50nMDA**: 50 nM dopamine application (mimics on-state DA levels)
- **quin**: 10 μM quinpirole (D2R agonist)
- **sul**: 10 μM sulpiride (D2R antagonist)
- **ACh**: 100 μM acetylcholine chloride (for sensor saturation/calibration, Fmax)
- **TTX**: 10 μM tetrodotoxin (blocks sodium channels, determines Fmin)

**Stimulation Types**:
- **single**: Single pulse electrical stimulation (1 ms × 0.3 mA)
- **burst**: Burst stimulation (20 pulses at 20 Hz, 1 ms × 0.3 mA each)

**Trial Numbering**:
- **-001, -002**: Two trials performed per stimulation protocol, data averaged for analysis
- Ensures reproducibility and reduces noise in measurements

### File Bundle per Experiment

Each experimental trial contains multiple file types organized in a structured directory:

```
BOT_04162024_slice2ROI1_sul_single-001/
├── BOT_04162024_slice2ROI1_sul_single-001_Cycle00001-botData.csv
├── BOT_04162024_slice2ROI1_sul_single-001_Cycle00001_Ch2_000001.ome.tif
├── BOT_04162024_slice2ROI1_sul_single-001_Cycle00001_Ch2_000002.ome.tif
├── BOT_04162024_slice2ROI1_sul_single-001_Cycle00001_Ch2_000003.ome.tif
... (continues to ~121 frames)
├── BOT_04162024_slice2ROI1_sul_single-001_Cycle00001_Ch2_000121.ome.tif
├── BOT_04162024_slice2ROI1_sul_single-001_Cycle00001_Ch3_000001.ome.tif
├── BOT_04162024_slice2ROI1_sul_single-001_Cycle00001_Ch3_000002.ome.tif
├── BOT_04162024_slice2ROI1_sul_single-001_Cycle00001_Ch3_000003.ome.tif
... (continues to ~121 frames)
├── BOT_04162024_slice2ROI1_sul_single-001_Cycle00001_VoltageOutput_001.xml
├── BOT_04162024_slice2ROI1_sul_single-001.env
└── BOT_04162024_slice2ROI1_sul_single-001.xml
```

**File Types and Contents**:
- **botData.csv**: Time series fluorescence data (brightness over time) for GRABACh3.0 sensor
  - Format: `BOT_[experiment]_Cycle00001-botData.csv`
- **Ch2_*.ome.tif**: GRABACh3.0 fluorescence channel (PMT gain 900) - main data channel
  - Format: `BOT_[experiment]_Cycle00001_Ch2_000XXX.ome.tif`
  - Typical count: ~106 frames per trial (varies slightly)
- **Ch3_*.ome.tif**: Dodt contrast/transmitted light channel (PMT gain 320) - reference/background
  - Format: `BOT_[experiment]_Cycle00001_Ch3_000XXX.ome.tif`
  - Typical count: ~91 frames per trial (fewer than Ch2)
- **VoltageOutput_001.xml**: Stimulation protocol metadata (timing, amplitude of electrical stimuli)
  - Format: `BOT_[experiment]_Cycle00001_VoltageOutput_001.xml`
- **.env**: Experiment environment metadata (microscope settings, PMT gains, other parameters)
  - Format: `BOT_[experiment].env`
- **.xml**: Session metadata including start time and acquisition parameters
  - Format: `BOT_[experiment].xml`

**Frame Count Observations**:
- **Ch2 frames**: Typically 106 frames per trial (at 21.26 fps ≈ 5 seconds total)
- **Ch3 frames**: Typically 91 frames per trial (shorter acquisition for reference channel)
- **Cycle designation**: All experiments use "Cycle00001" indicating single acquisition cycle per trial

**Complete Experimental Set per Slice/ROI**:

**Standard Set (PD and UL control groups)** - 18 experimental conditions:
```
BOT_[date]_slice[#]ROI[#]_ctr_single-001        # Control + single pulse - trial 1
BOT_[date]_slice[#]ROI[#]_ctr_single-002        # Control + single pulse - trial 2
BOT_[date]_slice[#]ROI[#]_ctr_burst-001         # Control + burst - trial 1
BOT_[date]_slice[#]ROI[#]_ctr_burst-002         # Control + burst - trial 2
BOT_[date]_slice[#]ROI[#]_50nMDA_single-001     # 50 nM DA + single - trial 1
BOT_[date]_slice[#]ROI[#]_50nMDA_single-002     # 50 nM DA + single - trial 2
BOT_[date]_slice[#]ROI[#]_50nMDA_burst-001      # 50 nM DA + burst - trial 1
BOT_[date]_slice[#]ROI[#]_50nMDA_burst-002      # 50 nM DA + burst - trial 2
BOT_[date]_slice[#]ROI[#]_quin_single-001       # Quinpirole + single - trial 1
BOT_[date]_slice[#]ROI[#]_quin_single-002       # Quinpirole + single - trial 2
BOT_[date]_slice[#]ROI[#]_quin_burst-001        # Quinpirole + burst - trial 1
BOT_[date]_slice[#]ROI[#]_quin_burst-002        # Quinpirole + burst - trial 2
BOT_[date]_slice[#]ROI[#]_sul_single-001        # Sulpiride + single - trial 1
BOT_[date]_slice[#]ROI[#]_sul_single-002        # Sulpiride + single - trial 2
BOT_[date]_slice[#]ROI[#]_sul_burst-001         # Sulpiride + burst - trial 1
BOT_[date]_slice[#]ROI[#]_sul_burst-002         # Sulpiride + burst - trial 2
BOT_[date]_slice[#]ROI[#]_ACh-001               # ACh saturation (Fmax calibration)
BOT_[date]_slice[#]ROI[#]_TTX-001               # TTX application (Fmin calibration)
```

**LID off Group Variation** - Uses "sliceXA" naming and may have additional ACh trials:
```
BOT_[date]_slice[#]A_ctr_single-001             # Control + single pulse - trial 1
BOT_[date]_slice[#]A_ctr_single-002             # Control + single pulse - trial 2
... (same pattern as above)
BOT_[date]_slice[#]A_ACh-001                    # ACh saturation - trial 1
BOT_[date]_slice[#]A_ACh-002                    # ACh saturation - trial 2 (sometimes present)
```

**Note**: Not all experiments include TTX calibration trials. Some LID off experiments have duplicate ACh calibration trials for improved Fmax determination.

## Key Metadata

### Channel Assignments
- **Channel 2 (Ch2)**: GRABACh3.0 fluorescence (PMT gain 900) - primary data channel
- **Channel 3 (Ch3)**: Dodt contrast/transmitted light (PMT gain 320) - reference channel

### Technical Specifications
**Dodt Contrast**: Gradient contrast transmitted light microscopy technique that provides high-contrast visualization of unstained living tissue without interfering with fluorescence detection. Uses principles similar to oblique contrast to produce images comparable to Differential Interference Contrast (DIC) without phase contrast artifacts.

**BOT (Brightness Over Time)**: Method to visualize changes in fluorescence intensity over time, used in imaging experiments to track dynamic processes like neurotransmitter release. The CSV files contain this time series data for quantitative analysis.

### BOT Data Analysis

> The whole image was the region of interest (ROI) used for analysis. Fluorescent intensity data were analyzed by custom Python code (accessible upon request). Briefly, the fluorescence intensity values were first background-subtracted (the background resulted from PMT was measured by imaging with same PMT voltage but zero laser power). Baseline fluorescence F0 was the average fluorescence over the 1 s period right before stimulation. ΔF = F − F0 was normalized by (Fmax − Fmin) and then analyzed.

**Analysis Pipeline**:
1. **Background subtraction**: PMT dark current removal (measured with same PMT voltage but zero laser power)
2. **Baseline calculation**: F0 = average fluorescence over 1 s period immediately before stimulation
3. **ΔF calculation**: ΔF = F - F0 (change from baseline)
4. **Normalization**: ΔF/(Fmax - Fmin) using calibration values
5. **ROI analysis**: Whole image used as region of interest for fluorescence quantification

**Data Processing Notes**:
- Custom Python analysis code available upon request from authors
- Standardized normalization allows comparison across experiments and conditions
- Time series data in CSV files represents processed fluorescence values
- Analysis accounts for PMT noise and baseline drift

### Acquisition Timeline
- **Baseline**: 3 s pre-stimulation recording
- **Stimulation**: Delivered at t = 3 s
- **Post-stimulation**: ≥5 s continued recording
- **Total duration**: ≥8 s per trial
- **Frame count**: ~121 frames per trial (at 21.26 fps)
- **Temporal resolution**: ~47 ms per frame

### Data Availability
**Stimulus Data**: No direct electrical stimulus waveform data available. Stimulation parameters are specified in methods but actual stimulus traces are not recorded in the dataset. The VoltageOutput XML files contain stimulation protocol metadata but not the actual voltage/current waveforms.

**File Extensions Found**:
```bash
.csv  # Time series fluorescence data
.env  # Environment/settings metadata
.tif  # Image stack frames (OME-TIFF format)
.xml  # Metadata and protocol information
```

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
- **Control (Unlesioned)**: Normal baseline ACh release with intact dopaminergic modulation
- **6-OHDA lesion**: Dramatically elevated ACh release due to loss of dopaminergic inhibition
- **LID off-state**: Sustained elevation in ACh release similar to parkinsonian state
- **State-dependent modulation**: ACh levels oscillate counter to dopamine between on/off states

### Pharmacological Responses
- **Quinpirole (D2R agonist)**: Strongly suppresses ACh release in all conditions
- **Sulpiride (D2R antagonist)**: Increases ACh release in control, no effect in lesioned tissue
- **50 nM Dopamine**: Mimics on-state conditions, strongly inhibits ACh release in lesioned mice
- **D2R signaling preservation**: Dopamine continues to robustly modulate ACh in dyskinetic tissue

### Mechanistic Insights
- **Cholinergic disinhibition**: Loss of dopaminergic input leads to ChI hyperactivity
- **Preserved D2R function**: D2 receptors on ChIs remain functional in parkinsonian state
- **Counter-oscillation**: ACh and dopamine levels oscillate in opposite directions between LID states
- **Off-state pathophysiology**: Elevated ACh during off-state drives aberrant synaptic plasticity in SPNs
- **M1R-mediated effects**: Cholinergic signaling through M1Rs on iSPNs contributes to LID severity
- **Therapeutic target**: Disrupting M1R/CalDAG-GEFI signaling reduces dyskinesia while preserving levodopa benefits

### Clinical Implications
- **State-specific treatment**: Understanding on/off-state differences may guide timing of interventions
- **Cholinergic modulation**: Targeting specific muscarinic receptor subtypes could provide therapeutic benefits
- **Preserved dopaminergic control**: D2R-mediated ACh regulation remains intact, suggesting potential for targeted therapies
- **Gene therapy potential**: Cell-type specific targeting of M1Rs or CalDAG-GEFI in iSPNs represents a novel therapeutic approach

## NWB Conversion Considerations

### Primary Data Types
- **Fluorescence time series**: GRABACh3.0 sensor data (botData.csv files)
- **Image stacks**: Two-photon microscopy frames (Ch2 and Ch3 OME-TIFF files)
- **Stimulation metadata**: Electrical stimulation protocols (VoltageOutput XML files)
- **Session metadata**: Experimental parameters and settings (XML and ENV files)

### Data Structure for NWB
- **OpticalChannel**: Define Ch2 (GRABACh3.0) and Ch3 (Dodt) channels with appropriate metadata
- **TwoPhotonSeries**: Store image stack data with proper temporal and spatial resolution
- **RoiResponseSeries**: Store processed fluorescence time series from botData.csv
- **Stimulation**: Represent electrical stimulation events and protocols
- **Pharmacology**: Document drug applications and concentrations

### Metadata Requirements
- **Subject information**: Lesion status, treatment group, experimental timeline
- **Imaging parameters**: Laser settings, PMT gains, pixel size, frame rate
- **Pharmacological treatments**: Drug concentrations, application timing, perfusion protocols
- **Stimulation protocols**: Electrode placement, stimulus parameters, timing
- **Analysis parameters**: ROI definition, normalization methods, calibration values

### Challenges and Solutions
- **Multiple file formats**: Integrate CSV, TIFF, XML, and ENV data into unified NWB structure
- **State designation**: Properly annotate on-state vs off-state experiments
- **ROI handling**: Account for whole-image ROI analysis approach
- **Calibration data**: Include Fmin/Fmax calibration for proper normalization
- **Trial averaging**: Represent both individual trials and averaged data
- **Data orientation warnings**: NWB inspector may flag TwoPhotonSeries orientation (time vs spatial dimensions) and report fewer frames during validation. This is expected for acetylcholine biosensor experiments (~45 frames, 8.78s) with high spatial resolution (128×128 pixels). The actual NWB files contain all available frames (45 per channel) as confirmed by direct TIFF file analysis. The orientation is scientifically correct for capturing biosensor responses.

### Actual vs Protocol Acquisition Parameters

**Protocol specifications vs actual acquisition:**
- **Protocol frame rate**: 21.26 fps (likely incorrect or refers to different mode)
- **Hardware-limited frame rate**: 5.05 fps (scan line rate 646.2 Hz ÷ 128 pixels)
- **Actual frame rate**: 5.12 fps (45 frames ÷ 8.78s) - matches hardware limit
- **Protocol duration**: ≥8s total
- **Actual duration**: 8.78s total
- **Frames available**: 45 frames per channel (Ch2 and Ch3)
- **File format**: Individual OME-TIFF files per frame per channel

**Technical Analysis**: The actual acquisition rate (5.12 fps) matches the hardware-limited maximum frame rate (5.05 fps) calculated from scan line rate ÷ image height. The protocol specification of 21.26 fps appears to be incorrect or refers to a different imaging mode. The BOT acquisition captured data at the maximum possible rate for 128×128 images on this system, providing complete temporal coverage of acetylcholine release dynamics.

## XML-Based Stimulation Type Detection

### Background
The Figure 5 conversion script initially attempted to determine stimulation type from BOT folder names using complex regex patterns. However, this approach proved unreliable due to inconsistent naming conventions and numerous edge cases. The solution was to extract stimulation information directly from the VoltageOutput XML files.

### Detection Procedure
1. **XML File Location**: Each experimental trial folder contains a VoltageOutput XML file:
   ```
   BOT_[experiment]_Cycle00001_VoltageOutput_001.xml
   ```

2. **Calibration Trial Detection**:
   - **Calibration trials** (TTX, ACh) have **no VoltageOutput XML file**
   - Absence of this file indicates no electrical stimulation was applied
   - Function returns `"calibration"` when XML file is missing

3. **Stimulation Type Parsing**: For trials with electrical stimulation, the XML contains:
   ```xml
   <Experiment>
     <Name>1pulse</Name>  <!-- or "20pulses@20Hz" -->
     <Waveform>
       <Name>LED</Name>
       <Enabled>true</Enabled>
       <WaveformComponent_PulseTrain>
         <PulseCount>1</PulseCount>  <!-- or 20 for burst -->
       </WaveformComponent_PulseTrain>
     </Waveform>
   </Experiment>
   ```

### XML Structure Analysis

#### Single Pulse Stimulation
```xml
<Experiment>
  <Name>1pulse</Name>
  <Waveform>
    <Name>LED</Name>
    <Enabled>true</Enabled>
    <WaveformComponent_PulseTrain>
      <PulseCount>1</PulseCount>
      <PulseWidth>1</PulseWidth>
      <PulseSpacing>499</PulseSpacing>
      <PulsePotentialStart>5</PulsePotentialStart>
      <FirstPulseDelay>10000</FirstPulseDelay>
    </WaveformComponent_PulseTrain>
  </Waveform>
</Experiment>
```

#### Burst Stimulation (20 pulses @ 20Hz)
```xml
<Experiment>
  <Name>20pulses@20Hz</Name>
  <Waveform>
    <Name>LED</Name>
    <Enabled>true</Enabled>
    <WaveformComponent_PulseTrain>
      <PulseCount>20</PulseCount>
      <PulseWidth>1</PulseWidth>
      <PulseSpacing>49</PulseSpacing>
      <PulsePotentialStart>5</PulsePotentialStart>
      <FirstPulseDelay>10000</FirstPulseDelay>
    </WaveformComponent_PulseTrain>
  </Waveform>
</Experiment>
```

### Key Detection Logic
- **LED Channel**: The stimulation waveform uses the "LED" channel name, which corresponds to the electrical stimulator output
- **PulseCount**: Distinguishes between single pulse (`PulseCount=1`) and burst stimulation (`PulseCount=20`)
- **Experiment Name**: Provides human-readable confirmation (`"1pulse"` vs `"20pulses@20Hz"`)

### Implementation Functions
```python
def get_stimulation_type_from_xml(bot_folder: Path) -> str:
    """Determine stimulation type from VoltageOutput XML file.

    Returns:
        - "single" for PulseCount=1
        - "burst" for PulseCount=20
        - "calibration" for missing XML (TTX/ACh trials)
    """
```

### Advantages of XML-Based Detection
1. **Eliminates 280+ edge cases** from folder name parsing
2. **Direct source of truth** - actual stimulation parameters used
3. **Simplified code** - reduced from ~150 lines to ~50 lines
4. **Reliable detection** - no assumptions about naming conventions
5. **Future-proof** - works regardless of folder naming changes

### Legacy Folder-Name Approach Issues
The previous approach used complex regex patterns to infer stimulation type from folder names:
- Required handling numerous edge cases and exceptions
- Made assumptions about naming consistency
- Failed when folder names deviated from expected patterns
- Created maintenance burden with growing exception lists

The XML-based approach provides definitive stimulation information directly from the acquisition system metadata, eliminating the need for folder name inference entirely.

## Related Files
- [Overview](../conversion_notes_overview.md)
- [Figure 1 Notes](figure_1_conversion_notes.md) - Two-photon imaging methods
- [Figure 2 Notes](figure_2_conversion_notes.md) - Electrical stimulation protocols
