# Surmeier Lab Conversion Progress

# Pre-Conversion

- [x] Repo Setup
- [x] Initial Inspection and Notes: Delineation of *Projects*, *Experiments* and *Data Streams*: [notes](https://github.com/catalystneuro/surmeier-lab-to-nwb/blob/main/src/surmeier_lab_to_nwb/zhai2025/conversion_notes.md)
- [x] Identify and request missing data/metadata/READMEs
- [x] Acquire **all** data needed for conversion

# Scope of Work points

[Scope of Work Link (Private)](https://docs.google.com/document/d/1ITmVEEOQ1TbBC8hvZvkdQ-VcdW9_ynUWFwfRxBBrAGQ/edit?tab=t.0#heading=h.6ol4x2q1uqwd)

These are the points that we wrote in the scope of work

## Build interfaces to convert the following data streams to NWB format:
- [x] Convert ABF format: this data is not available and it will not be converted
- [x] Convert Bruker optical sensor recordings with appropriate metadata:
    - [x] Line Scans [Line Scan Interface](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/2)
    - [x] Fluorescence traces [Prairie View BOT interface](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/12) | [fix Bruker roiextractors](https://github.com/catalystneuro/roiextractors/pull/438) | [fix Bruker neuroconv](https://github.com/catalystneuro/neuroconv/pull/1375)
    - [x] Prairie View Optogenetical Stimuli [PrairieViewOptogeneticsInterface](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/14)
- [x] Convert electrophysiological data acquired with the Bruker system [Pairie View Interface](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/5)
- [x] Build an interface for manual segmentation data in Surmeier lab format
    - [x] Two-photon laser scanning microscopy image stacks for spine density [Neuroconv Interface Improvement 1](https://github.com/catalystneuro/neuroconv/pull/1365) | [Neuroconv Interface Improvement 2](https://github.com/catalystneuro/neuroconv/pull/1439) | [Neuroconv Interface Improvement 3](https://github.com/catalystneuro/neuroconv/pull/1441)
    - [x] Confocal microscopy image stacks for spine density Nikon  [Figure 4](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/13)
- [x] Integrate behavioral annotations from Surmeier lab custom format [Behavioral AIM Scoring Interface](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/32)
- [x] Include behavioral video recordings [Conversion Scripts for Video](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/35)
- [x] Include electrical stimulation signals and metadata [All the Prairie View Interfaces Extract Stimulation Metadata]



# Conversion of figures data

The data, however, was packaged per figure so I am gonna show progress per figure as well here.

## Figure 1:
Intracellular electrophysiology (VoltageClamp), Image Stacks and LineScans.

### Somatic Excitability
- [x] Build needed interfaces [Pairie View Interface](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/5)
- [x] Conversion script [PR](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/1)
- [x] Time alignment [PR](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/24)
- [x] Add icephys table hierarchical structure [PR](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/24)
- [x] Inspector [Commit](https://github.com/catalystneuro/surmeier-lab-to-nwb/commit/79930c77c568c646de43c4f4d66268e022c35602)
- [x] Upload data

### Dendritic Excitability
- [x] Build needed interfaces [Pairie View Interface](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/5) | [Line Scan Interface](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/2)
- [x] Conversion script [PR](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/3)
- [x] Time alignment [PR](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/19)
- [x] Add icephys table hierarchical structure [PR](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/20)
- [x] Inspector [Commit](https://github.com/catalystneuro/surmeier-lab-to-nwb/commit/79930c77c568c646de43c4f4d66268e022c35602)
- [x] Upload data




## Figure 2:
This has the image stacks (spine density) and intracellular electrophysiology (CurrentClamp) data with optogenetic stimulation (Sr-oEPSC). The image stacks here are from two-photon laser scanning microscopy.

### Spine density image stacks
- [x] Build needed interfaces [Neuroconv PR](https://github.com/catalystneuro/neuroconv/pull/1365)
- [x] Conversion script [PR](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/7)
- [x] [Commit](https://github.com/catalystneuro/surmeier-lab-to-nwb/commit/d6242917671c51a21b693241edb5774eecafe364)
- [x] Upload data

### Voltage clamp with optogenetic stimulation
- [x] Build needed interfaces [Pairie View Interface](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/5) | [PrairieViewOptogeneticsInterface](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/14)
- [x] Conversion script [PR](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/14)
- [x] Time alignment [PR](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/28)
- [x] Add icephys table hierarchical structure [PR](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/28)
- [x] Use ndx-optogenetics extension [PR](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/29)
- [x] Inspector [Commit](https://github.com/catalystneuro/surmeier-lab-to-nwb/commit/d6242917671c51a21b693241edb5774eecafe364)
- [x] Upload data

## Figure 3:
Intracellular electrophysiology (VoltageClamp) and LineScan data.

### Somatic Excitability
- [x] Build needed interfaces [Pairie View Interface](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/5)
- [x] Conversion script [PR](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/15)
- [x] Time alignment of the data [PR](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/23)
- [x] Add icephys table hierarchical structure [PR](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/23)
- [x] Inspector [Commit](https://github.com/catalystneuro/surmeier-lab-to-nwb/commit/d6242917671c51a21b693241edb5774eecafe364)
- [x] Upload data

### Dendritic Excitability
- [x] Build needed interfaces [Pairie View Interface](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/5)
- [x] Conversion script [PR](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/16)
- [x] Time alignment [PR](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/17)
- [x] Add icephys table hierarchical structure [PR](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/18)
- [x] Inspector [Commit](https://github.com/catalystneuro/surmeier-lab-to-nwb/commit/d6242917671c51a21b693241edb5774eecafe364)
- [x] Upload data

## Figure 4:

### Voltage clamp with optogenetic stimulation
- [x] Build needed interfaces [Pairie View Interface](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/5)
- [x] Conversion script [PR](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/30)
- [x] Time alignment [PR](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/30)
- [x] Add icephys table hierarchical structure [PR](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/30)
- [x] Use ndx-optogenetics extension [PR](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/30)
- [x] Inspector [Commit](https://github.com/catalystneuro/surmeier-lab-to-nwb/commit/b77f8a7dc53c62091597f9975d04d4ab26a96eef)
- [x] Upload data


### Spine Density
- [x] Build needed interfaces [Neuroconv PR](https://github.com/catalystneuro/neuroconv/pull/1365)
- [x] Conversion script [PR](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/36)
- [x] Inspector [Commit](https://github.com/catalystneuro/surmeier-lab-to-nwb/commit/b77f8a7dc53c62091597f9975d04d4ab26a96eef)
- [x] Upload data

### Confocal Spine Density
Image stacks (spine density) like figure 2 but they come in a different format as they come from confocal microscopy.
- [x] Build needed interfaces [Nikon Image Stack Interface](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/13)
- [x] Conversion script [PR](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/37)
- [x] Inspector [Commit](https://github.com/catalystneuro/surmeier-lab-to-nwb/commit/b77f8a7dc53c62091597f9975d04d4ab26a96eef)
- [x] Upload data

## Figure 5:
Fluorescence imaging data from two-photon laser scanning microscopy.
- [x] Build needed interfaces [Prairie View BOT interface](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/12)
- [x] Conversion script [PR](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/38)
- [x] Inspector [Commit](https://github.com/catalystneuro/surmeier-lab-to-nwb/commit/76858abb498f606b42804a49168619dd2c6ac333)
- [x] Upload data

## Figure 6:
Intracellular electrophysiology (VoltageClamp), two-photon laser stacks for spine density and line scans.
Structure: Three experiment types (Dendritic excitability, Somatic excitability, Spine density), each with control and M1R antagonist conditions.

### Dendritic Excitability
- [x] Build needed interfaces [Pairie View Interface](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/5) | [Line Scan Interface](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/2)
- [x] Conversion script [PR](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/21)
- [x] Time alignment [PR](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/21)
- [x] Add icephys table hierarchical structure [PR](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/21)
- [x] Inspector [Commit](https://github.com/catalystneuro/surmeier-lab-to-nwb/commit/3edb76c5195ac9c071d91065eb455831aeb55549)
- [x] Upload data


### Somatic Excitability
- [x] Build needed interfaces [Pairie View Interface](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/5)
- [x] Conversion script [PR](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/26)
- [x] Time alignment [PR](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/26)
- [x] Add icephys table hierarchical structure [PR](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/26)
- [x] Inspector [Commit](https://github.com/catalystneuro/surmeier-lab-to-nwb/commit/3edb76c5195ac9c071d91065eb455831aeb55549)
- [x] Upload data

### Spine Density
- [x] Build needed interfaces [Neuroconv Interface](https://github.com/catalystneuro/neuroconv/pull/1365)
- [x] Conversion script [PR](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/39)
- [x] Inspector [Commit](https://github.com/catalystneuro/surmeier-lab-to-nwb/commit/3edb76c5195ac9c071d91065eb455831aeb55549)
- [x] Upload data

## Figure 7:
CDGI knockout effects on iSPN function and behavior. Includes somatic excitability, dendritic excitability, spine density, and behavioral assessment comparing CDGI knockout vs wildtype mice.

### Somatic Excitability (CDGI KO on vs off)
- [x] Build needed interfaces [Pairie View Interface](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/5)
- [x] Conversion script (handle CDGI KO off-state and on-state conditions) [PR](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/25)
- [x] Time alignment [PR](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/25)
- [x] Add icephys table hierarchical structure [PR](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/25)
- [x] Inspector [Commit](https://github.com/catalystneuro/surmeier-lab-to-nwb/commit/76858abb498f606b42804a49168619dd2c6ac333)
- [x] Upload data

### Dendritic Excitability (CDGI KO on vs off)
- [x] Build needed interfaces [Pairie View Interface](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/5) | [Line Scan Interface](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/2)
- [x] Conversion script (handle CDGI KO off-state and on-state conditions) [PR](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/22)
- [x] Time alignment [PR](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/22)
- [x] Add icephys table hierarchical structure [PR](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/22)
- [x] Inspector [Commit](https://github.com/catalystneuro/surmeier-lab-to-nwb/commit/76858abb498f606b42804a49168619dd2c6ac333)
- [x] Upload data

### Dendritic Excitability (oxoM response - CDGI KO vs WT)
- [x] Build needed interfaces [Pairie View Interface](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/5) | [Line Scan Interface](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/2)
- [x] Conversion script (handle paired before/after oxoM protocol for both KO and WT) [PR](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/27)
- [x] Time alignment [PR](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/27)
- [x] Add icephys table hierarchical structure [PR](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/27)
- [x] Inspector [Commit](https://github.com/catalystneuro/surmeier-lab-to-nwb/commit/76858abb498f606b42804a49168619dd2c6ac333)
- [x] Upload data

### Spine Density (CDGI KO on vs off)
- [x] Build needed interfaces [Neuroconv Interface](https://github.com/catalystneuro/neuroconv/pull/1365)
- [x] Conversion script [PR](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/39)
- [x] Inspector [Commit](https://github.com/catalystneuro/surmeier-lab-to-nwb/commit/76858abb498f606b42804a49168619dd2c6ac333)
- [x] Upload data

### Behavioral Assessment (AIM scoring and contralateral rotations)
- [x] Build needed interfaces (behavioral scoring interface needed) [PR](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/32)
- [x] Conversion script (AIM scores and rotation videos) [PR](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/32)
- [x] Inspector [Commit](https://github.com/catalystneuro/surmeier-lab-to-nwb/commit/76858abb498f606b42804a49168619dd2c6ac333)
- [x] Upload data

### Behavioral Videos (Contralateral rotation documentation)
- [x] Build video interfaces (ExternalVideoInterface)
- [x] Conversion script for Figure 7 behavioral videos [PR](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/35)
- [x] Inspector [Commit](https://github.com/catalystneuro/surmeier-lab-to-nwb/commit/76858abb498f606b42804a49168619dd2c6ac333)
- [x] Upload data
**Data**: 434 video files documenting contralateral rotation behavior in CDGI KO vs WT mice across multiple experimental sessions (2017-2020). Videos recorded at 20-120 min intervals post L-DOPA injection to assess therapeutic benefit and dyskinesia severity.

## Figure 8:
M1R CRISPR deletion effects on iSPN properties and dyskinetic behaviors using targeted gene editing. Includes somatic excitability, spine density, and behavioral assessment comparing M1R CRISPR vs control mice.

### Somatic Excitability (M1R CRISPR vs interleaved control)
- [x] Build needed interfaces [Pairie View Interface](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/5)
- [x] Conversion script (handle M1R CRISPR and interleaved control conditions) [PR](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/31)
- [x] Time alignment [PR](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/31)
- [x] Add icephys table hierarchical structure [PR](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/31)
- [x] Inspector [Commit](f4b989da3b66620c0164a834650acc26d82fa7d0)
- [x] Upload data

### Spine Density (M1R CRISPR vs control - off-state)
- [x] Build needed interfaces [Neuroconv PR](https://github.com/catalystneuro/neuroconv/pull/1365)
- [x] Conversion script [PR](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/39)
- [x] Inspector [Commit](f4b989da3b66620c0164a834650acc26d82fa7d0)
- [x] Upload data

### Behavioral Assessment (AIM scoring)
- [x] Build needed interfaces [PR](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/32)
- [x] Conversion script (AIM scores for M1R CRISPR mice) [PR](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/33)
- [x] Inspector [Commit](f4b989da3b66620c0164a834650acc26d82fa7d0)
- [x] Upload data

## Supplementary Figure 3:
M1R CRISPR behavioral video documentation showing dyskinesia assessment in gene-edited mice.

### Behavioral Videos (M1R CRISPR dyskinesia assessment)
- [x] Build video interfaces (ExternalVideoInterface)
- [x] Conversion script for Supplementary Figure 3 behavioral videos [PR](https://github.com/catalystneuro/surmeier-lab-to-nwb/pull/35)
- [x] Inspector [Commit](https://github.com/catalystneuro/surmeier-lab-to-nwb/commit/76858abb498f606b42804a49168619dd2c6ac333)
- [x] Upload data
**Data**: 464 video files documenting dyskinetic behavior in M1R CRISPR vs Control mice across multiple experimental sessions (2023-2024). Videos recorded at 20-120 min intervals post L-DOPA injection to assess effects of M1R gene editing on dyskinesia severity.

## Nice to have

First ensure that the raw data is available and on dandi and then we can do these improvements:
- [x] Extract the times of the intracellular events so we can write a single time series instead of many.
- [ ] Add an extension for line scan data
- [ ] Add an extension for image stacks that keeps the metadata of the microscopy
- [x] Separate the nwbfiles per subject
- [x] use ndx-optogenetics

# Data uploading and conversion packaging
- [x] Create Dandiset(s): [Link](https://dandiarchive.org/dandiset/001538)
- [x] Add protocols.io to dandiset
- [x] Add anatomy metadata to Dandiset
- [x] Add disease metadata to Dandiset
- [x] README/Documentation here in the repo.
- [x] Example Notebooks:
- [x] Final Meeting

Example of neurosift access with the modalities:


- Surmeier Presentation
https://neurosift.app/dandiset/001538
Dataset (ecephys):
https://neurosift.app/nwb?url=https://api.dandiarchive.org/api/assets/b341d3b5-8413-49d4-b2c9-f93bd0985778/download/&dandisetId=001538&dandisetVersion=draft
Dataset (video)
https://neurosift.app/nwb?url=https://api.dandiarchive.org/api/assets/c574b342-ee3b-4deb-8bcd-0104c64b78a2/download/&dandisetId=001538&dandisetVersion=draft
Datset (optogenetics)
https://neurosift.app/nwb?url=https://api.dandiarchive.org/api/assets/276d2b36-18cc-4f68-be2f-9fa1ea288abc/download/&dandisetId=001538&dandisetVersion=draft
Dataset (ophys)
https://neurosift.app/nwb?url=https://api.dandiarchive.org/api/assets/812b3bae-285d-481e-9001-97f0247407a1/download/&dandisetId=001538&dandisetVersion=draft
Dataset (figure 5)
https://neurosift.app/nwb?url=https://api.dandiarchive.org/api/assets/c8ecbe04-6a51-4dae-a3fa-43ba3fa6c5f8/download/&dandisetId=001538&dandisetVersion=draft
