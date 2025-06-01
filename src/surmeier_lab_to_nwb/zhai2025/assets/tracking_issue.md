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
- [x] Convert electrophysiological data acquired with the Bruker system [Pairie View Interface](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/5)
- [ ] Build an interface for manual segmentation data in Surmeier lab format
    - [x] Two-photon laser scanning microscopy image stacks for spine density [Figure 2](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/7)   [Neuroconv Interface](https://github.com/catalystneuro/neuroconv/pull/1365)
    - [ ] Confocal microscopy image stacks for spine density
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
- [ ] Metadata
- [ ] Inspector
- [ ] Upload data

### Dendritic Excitability
- [x] Build needed interfaces [Pairie View Interface](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/5) [Line Scan Interface](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/2)
- [x] Conversion script [PR](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/3)
- [ ] Metadata
- [ ] Inspector
- [ ] Upload data


## Figure 2:
This has the image stacks (spine density) and intracellular electrophysiology (CurrentClamp) data with optogenetic stimulation (Sr-oEPSC). The image stacks here are from two-photon laser scanning microscopy.

### Spine density image stacks
- [x] Build needed interfaces [Neuroconv PR](https://github.com/catalystneuro/neuroconv/pull/1365)
- [x] Conversion script [PR](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/7)
- [ ] Metadata
- [ ] Inspector
- [ ] Upload data

### CurrentClamp with optogenetic stimulation
- [x] Build needed interfaces [Pairie View Interface](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/5)
- [ ] Conversion script
- [ ] Metadata
- [ ] Inspector
- [ ] Upload data

## Figure 3:
Intracellular electrophysiology (VoltageClamp) and LineScan data.

## Somatic Excitability
- [x] Build needed interfaces [Pairie View Interface](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/5)
- [ ] Conversion script
- [ ] Metadata
- [ ] Inspector
- [ ] Upload data

### Dendritic Excitability
- [x] Build needed interfaces [Pairie View Interface](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/5)
- [ ] Conversion script
- [ ] Metadata
- [ ] Inspector
- [ ] Upload data

## Figure 4:
Image stacks (spine density) like figure 2 but they come in a different format as they come from confocal microscopy.
- [ ] Build needed interfaces
- [ ] Conversion script
- [ ] Metadata
- [ ] Inspector
- [ ] Upload data

## Figure 5:
Fluorescence imaging data from two-photon laser scanning microscopy.
- [x] Build needed interfaces [Prairie View BOT interface](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/12)
- [ ] Conversion script
- [ ] Metadata
- [ ] Inspector
- [ ] Upload data

## Figure 6:
Intracellular electrophysiology (VoltageClamp), two-photon laser stacks for spine density and line scans.
- [ ] Build needed interfaces
- [ ] Conversion script
- [ ] Metadata
- [ ] Inspector
- [ ] Upload data

## Figure 7:
Intracellular electrophysiology (VoltageClamp), two-photon laser stacks for spine density, and line scans. Behavior scoring data (AIM) and behavioral videos.
- [ ] Build needed interfaces
- [ ] Conversion script
- [ ] Metadata
- [ ] Inspector
- [ ] Upload data

## Figure 8:
Behavioral scoring, intracellular electrophysiology (VoltageClamp), two-photon laser stacks for spine density, and line scans.
- [ ] Build needed interfaces
- [ ] Conversion script
- [ ] Metadata
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
