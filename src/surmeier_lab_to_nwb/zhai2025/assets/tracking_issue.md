# Surmeier Lab Conversion Progress

# Pre-Conversion

- [x] Repo Setup
- [x] Initial Inspection and Notes: Delineation of *Projects*, *Experiments* and *Data Streams*: [notes](https://github.com/catalystneuro/surmeier-lab-to-nwb/blob/main/src/surmeier_lab_to_nwb/zhai2025/conversion_notes.md)
- [x] Identify and request missing data/metadata/READMEs
- [x] Acquire **all** data needed for conversion

# Scope of Work points

These are the points that we wrote in the scope of work

## Build interfaces to convert the following data streams to NWB format:
- [x] Convert ABF format: this data is not available and it will not be converted
- [x] Convert Bruker optical sensor recordings with appropriate metadata:
    - [x] Line Scans [Line Scan Interface](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/2)
    - [x] Fluorescence traces [Prairie View BOT interface](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/12) [fix Bruker roiextractors](https://github.com/catalystneuro/roiextractors/pull/438) [fix Bruker neuroconv](https://github.com/catalystneuro/neuroconv/pull/1375)
    - [x] Prairie View Optogenetical Stimuli [PrairieViewOptogeneticsInterface](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/14)
- [x] Convert electrophysiological data acquired with the Bruker system [Pairie View Interface](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/5)
- [x] Build an interface for manual segmentation data in Surmeier lab format
    - [x] Two-photon laser scanning microscopy image stacks for spine density [Figure 2](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/7)   [Neuroconv Interface Bug](https://github.com/catalystneuro/neuroconv/pull/1365)
    - [x] Confocal microscopy image stacks for spine density Nikon  [Figure 4](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/13)
- [ ] Integrate behavioral annotations from Surmeier lab custom format
- [ ] Include behavioral video recordings
- [ ] Include electrical stimulation signals and metadata

# Conversion of figures data

The data, however, was packaged per figure so I am gonna show progress per figure as well here.

## Figure 1:
Intracellular electrophysiology (VoltageClamp), Image Stacks and LineScans.

### Somatic Excitability
- [x] Build needed interfaces [Pairie View Interface](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/5)
- [x] Conversion script [PR](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/1)
- [ ] Time alignment
- [ ] Add icephys table hierarchical structure
- [ ] Double check metadata
- [ ] Inspector
- [ ] Upload data

### Dendritic Excitability
- [x] Build needed interfaces [Pairie View Interface](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/5) [Line Scan Interface](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/2)
- [x] Conversion script [PR](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/3)
- [x] Time alignment [PR](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/19)
- [x] Add icephys table hierarchical structure [PR](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/20)
- [ ] Double check metadata
- [ ] Inspector
- [ ] Upload data


## Figure 2:
This has the image stacks (spine density) and intracellular electrophysiology (CurrentClamp) data with optogenetic stimulation (Sr-oEPSC). The image stacks here are from two-photon laser scanning microscopy.

### Spine density image stacks
- [x] Build needed interfaces [Neuroconv PR](https://github.com/catalystneuro/neuroconv/pull/1365)
- [x] Conversion script [PR](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/7)
- [ ] Time alignment
- [ ] Double check metadata
- [ ] Inspector
- [ ] Upload data

### Current clamp with optogenetic stimulation
- [x] Build needed interfaces [Pairie View Interface](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/5) [PrairieViewOptogeneticsInterface](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/14)
- [x] Conversion script [PR](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/14)
- [ ] Time alignment
- [ ] Add icephys table hierarchical structure
- [ ] Double check metadata
- [ ] Inspector
- [ ] Upload data

## Figure 3:
Intracellular electrophysiology (VoltageClamp) and LineScan data.

### Somatic Excitability
- [x] Build needed interfaces [Pairie View Interface](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/5)
- [x] Conversion script [PR](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/15)
- [ ] Time alignment of the data
- [ ] Add icephys table hierarchical structure
- [ ] Double check metadata
- [ ] Inspector
- [ ] Upload data

### Dendritic Excitability
- [x] Build needed interfaces [Pairie View Interface](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/5)
- [x] Conversion script [PR](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/16)
- [x] Time alignment [PR](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/17)
- [x] Add icephys table hierarchical structure [PR](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/18)
- [ ] Double check metadata
- [ ] Inspector
- [ ] Upload data

## Figure 4:

### Current clamp with optogenetic stimulation
- [x] Build needed interfaces [Pairie View Interface](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/5)
- [ ] Conversion script
- [ ] Time alignment
- [ ] Add icephys table hierarchical structure
- [ ] Double check metadata
- [ ] Inspector
- [ ] Upload data

## Spine Density
- [ ] Build needed interfaces
- [ ] Conversion script
- [ ] Time alignment
- [ ] Add icephys table hierarchical structure (N/A - no electrophysiology)
- [ ] Double check metadata
- [ ] Inspector
- [ ] Upload data

## Confocal Spine Density
Image stacks (spine density) like figure 2 but they come in a different format as they come from confocal microscopy.
- [x] Build needed interfaces [Nikon Image Stack Interface](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/13)
- [ ] Conversion script
- [ ] Time alignment
- [ ] Add icephys table hierarchical structure (N/A - no electrophysiology)
- [ ] Double check metadata
- [ ] Inspector
- [ ] Upload data

## Figure 5:
Fluorescence imaging data from two-photon laser scanning microscopy.
- [x] Build needed interfaces [Prairie View BOT interface](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/12)
- [ ] Conversion script
- [ ] Time alignment
- [ ] Add icephys table hierarchical structure (N/A - no electrophysiology)
- [ ] Double check metadata
- [ ] Inspector
- [ ] Upload data

## Figure 6:
Intracellular electrophysiology (VoltageClamp), two-photon laser stacks for spine density and line scans.
Structure: Three experiment types (Dendritic excitability, Somatic excitability, Spine density), each with control and M1R antagonist conditions.

### Dendritic Excitability
- [x] Build needed interfaces [Pairie View Interface](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/5) [Line Scan Interface](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/2)
- [x] Conversion script (handle both control and M1R antagonist conditions) [PR](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/21)
- [x] Time alignment [PR](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/21)
- [x] Add icephys table hierarchical structure [PR](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/21)
- [ ] Double check metadata
- [ ] Inspector
- [ ] Upload data

### Somatic Excitability
- [x] Build needed interfaces [Pairie View Interface](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/5)
- [ ] Conversion script (handle both control and M1R antagonist conditions)
- [ ] Time alignment
- [ ] Add icephys table hierarchical structure
- [ ] Double check metadata
- [ ] Inspector
- [ ] Upload data

### Spine Density
- [x] Build needed interfaces [Neuroconv Interface](https://github.com/catalystneuro/neuroconv/pull/1365)
- [ ] Conversion script (handle both control and M1R antagonist conditions)
- [ ] Double check metadata
- [ ] Inspector
- [ ] Upload data

## Figure 7:
CDGI knockout effects on iSPN function and behavior. Includes somatic excitability, dendritic excitability, spine density, and behavioral assessment comparing CDGI knockout vs wildtype mice.

### Somatic Excitability (CDGI KO on vs off)
- [x] Build needed interfaces [Pairie View Interface](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/5)
- [ ] Conversion script (handle CDGI KO off-state and on-state conditions)
- [ ] Time alignment
- [ ] Add icephys table hierarchical structure
- [ ] Double check metadata
- [ ] Inspector
- [ ] Upload data

### Dendritic Excitability (CDGI KO on vs off)
- [x] Build needed interfaces [Pairie View Interface](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/5) [Line Scan Interface](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/2)
- [x] Conversion script (handle CDGI KO off-state and on-state conditions) [PR](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/22)
- [x] Time alignment [PR](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/22)
- [x] Add icephys table hierarchical structure [PR](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/22)
- [ ] Double check metadata
- [ ] Inspector
- [ ] Upload data

### Dendritic Excitability (oxoM response - CDGI KO vs WT)
- [x] Build needed interfaces [Pairie View Interface](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/5) [Line Scan Interface](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/2)
- [ ] Conversion script (handle paired before/after oxoM protocol for both KO and WT)
- [ ] Time alignment
- [ ] Add icephys table hierarchical structure
- [ ] Double check metadata
- [ ] Inspector
- [ ] Upload data

### Spine Density (CDGI KO on vs off)
- [x] Build needed interfaces [Neuroconv Interface](https://github.com/catalystneuro/neuroconv/pull/1365)
- [ ] Conversion script (handle CDGI KO off-state and on-state conditions)
- [ ] Time alignment
- [ ] Add icephys table hierarchical structure (N/A - no electrophysiology)
- [ ] Double check metadata
- [ ] Inspector
- [ ] Upload data

### Behavioral Assessment (AIM scoring and contralateral rotations)
- [ ] Build needed interfaces (behavioral scoring interface needed)
- [ ] Conversion script (AIM scores and rotation videos)
- [ ] Time alignment
- [ ] Add icephys table hierarchical structure (N/A - no electrophysiology)
- [ ] Double check metadata
- [ ] Inspector
- [ ] Upload data

## Figure 8:
M1R CRISPR deletion effects on iSPN properties and dyskinetic behaviors using targeted gene editing. Includes somatic excitability, spine density, and behavioral assessment comparing M1R CRISPR vs control mice.

### Somatic Excitability (M1R CRISPR vs interleaved control)
- [x] Build needed interfaces [Pairie View Interface](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/5)
- [ ] Conversion script (handle M1R CRISPR and interleaved control conditions)
- [ ] Time alignment
- [ ] Add icephys table hierarchical structure
- [ ] Double check metadata
- [ ] Inspector
- [ ] Upload data

### Spine Density (M1R CRISPR vs control - off-state)
- [x] Build needed interfaces [Neuroconv Interface](https://github.com/catalystneuro/neuroconv/pull/1365)
- [ ] Conversion script (handle M1R CRISPR and control conditions)
- [ ] Time alignment
- [ ] Add icephys table hierarchical structure (N/A - no electrophysiology)
- [ ] Double check metadata
- [ ] Inspector
- [ ] Upload data

### Behavioral Assessment (AIM scoring)
- [ ] Build needed interfaces (behavioral scoring interface needed)
- [ ] Conversion script (AIM scores for M1R CRISPR mice)
- [ ] Time alignment
- [ ] Add icephys table hierarchical structure (N/A - no electrophysiology)
- [ ] Double check metadata
- [ ] Inspector
- [ ] Upload data

## Nice to have

First ensure that the raw data is available and on dandi and then we can do these improvements:
- [ ] Extract the times of the intracellular events so we can write a single time series instead of many.
- [ ] Add an extension for line scan data
- [ ] Add an extension for image stacks that keeps the metadata of the microscopy
- [ ] Separate the nwbfiles per subject

# Data uploading and conversion packaging
- [ ] Setup Dandiset(s):
- [ ] README/Documentation:
- [ ] Example Notebooks:
- [ ] Middle Meeting
- [ ] Final Meeting
