# Zhai 2025 Conversion Notes - Overview

## Paper
Here is the paper:
https://www.biorxiv.org/content/10.1101/2025.01.02.631090v1.full

## Scope of Work (Private Access)
https://docs.google.com/document/d/1ITmVEEOQ1TbBC8hvZvkdQ-VcdW9_ynUWFwfRxBBrAGQ/edit?usp=sharing

## Glossary

### Biological Terms
* **DLS** dorsolateral striatum
* **6-OHDA** 6-hydroxydopamine, a neurotoxin used to induce Parkinson's disease-like symptoms in animal models by selectively destroying dopaminergic neurons.
* **levodopa** the treatment for Parkinson's disease
* **levodopa-induced dyskinesia (LID)** a set of erratic movements induced by levodopa after long-term treatment
* **rheobase** the minimal electric current required to excite a tissue (as nerve or muscle) given an indefinitely long time during which the current is applied,
* **SPNs** spiny projection neurons
* **dSPNs** direct pathway SPNs
* **iSPNs** indirect pathway SPNs
* **AIM** lower abnormal involuntary movement (AIM)
* **Contralateral rotations** are a behavioral measure used in rodent models of Parkinson's disease. When one side of the brain is lesioned (for example, with 6-OHDA), the imbalance in motor circuitry causes the animal to rotate predominantly toward the side opposite the lesion. By counting these rotations after a treatment like levodopa, researchers can assess changes in motor function and dyskinetic behavior.
* **M1Rs**: Muscarinic acetylcholine receptors type 1. M1 receptor activation occurs when acetylcholine (the natural neurotransmitter) or other muscarinic agonists (such as muscarine, which is where these receptors get their name) bind to the receptor. M1Rs are G-protein coupled receptors that primarily signal through the Gq/11 pathway, leading to increased neuronal excitability in many cell types including striatal indirect pathway spiny projection neurons.
* **SCH** SCH23390, a pharmacological compound that acts as a selective antagonist at the dopamine D₁ receptor (D1R). In other words, SCH23390 blocks D1-type dopaminergic signaling in neurons.
* **SUL/sulpiride** sulpiride, a selective D2 dopamine receptor antagonist. By blocking D2-type receptors, sulpiride is used to test whether observed changes (in this case, decreased iSPN excitability in the on-state) depend on ongoing D2R signaling.
* **quinpirole or DA** D2R agonist(s)
* **CDGI** stands for CalDAG-GEFI, a calcium-activated guanine nucleotide exchange factor that is highly expressed in the striatum. It plays a critical role in linking M1 muscarinic receptor activation to intracellular signaling pathways that regulate dendritic excitability and synaptic plasticity in indirect pathway spiny projection neurons (iSPNs). In the paper, disrupting CDGI was found to blunt the dendritic adaptations in iSPNs and reduce dyskinetic behaviors, suggesting its important role in the pathophysiology of levodopa-induced dyskinesia.
* **D1**: Direct pathway spiny projection neurons (dSPNs) that express dopamine D1 receptors
* **D2**: Indirect pathway spiny projection neurons (iSPNs) that express dopamine D2 receptors
* **ACSF** stands for Artificial Cerebrospinal Fluid

### Technical Terms
* **2PLSM** Two Photon Laser Scanning Microscopy
* **Ex vivo** refers to a procedure or experiment performed on tissue taken from a living organism, but studied outside that organism's normal biological context (for example, in a chamber or petri dish)
* **oEPSC** optogenetically evoked postsynaptic current
* **GRABACh3.0** is a genetically encoded fluorescent sensor designed to detect acetylcholine (ACh) release. It is used to monitor ACh dynamics in real-time.
* **AAV5-hSyn-GRAB** refers to an adeno‐associated virus (serotype 5) that uses the human synapsin promoter to drive neuronal expression of the genetically encoded acetylcholine sensor GRABACh3.0. The authors injected this viral vector into the dorsolateral striatum so that they could optically monitor acetylcholine release in brain slices. This tool was critical for assessing how cholinergic signaling changes under conditions of dopamine depletion and during levodopa-induced dyskinesia.

### File Format Notes
The STK metadata in this context likely refers to metadata extracted from a MetaMorph STK file, a format used in microscopy for storing image sequences. STK files are essentially TIFF files with additional MetaMorph-specific metadata embedded in the tags. This metadata typically contains acquisition parameters, calibration details, and channel information.

These are tags appearing in the dendritic experiments of the Figures 1 and 3.

### Key Device
MultiClamp 700B - Electrophysiology amplifier used for patch clamp recordings

## Scope of Work (Conversion Part)

### Conversion of Data Streams

Build interfaces to convert the following data streams to NWB format:

- Convert ABF format: I have not found .abf data
- Convert Bruker optical sensor recordings with appropriate metadata: these are the line scans and the full-field timelapse image of figure 5.
- Convert electrophysiological data acquired with the Bruker system: the prairie view data.
- Build an interface for manual segmentation data in Surmeier lab format.
- Integrate behavioral annotations from Surmeier lab custom format. AIM scoring and contralateral rotations.
- Include behavioral video recordings. Videos are present.
- Include electrical stimulation signals and metadata. Those are the .xml data from PrairieView and the *VolageOutput*.xml files that contain a textual description of the stimulus.

Each conversion will utilize compression for efficient storage and implement chunking strategies to optimize for cloud storage. Detailed documentation will be provided for the installation and usage of the conversion software, including scripts to handle data from various protocols:
1. Ex-vivo brain slices with optogenetics: combine ABF recordings with optogenetic stimulus information.
2. Electrophysiology combined with optical imaging. Bruker optical sensor data with electrical recordings.
3. Two-photon laser scanning microscopy with electrical stimulation: process imaging data with stimulation events
4. Behavioral pharmacology: convert videos, tracking data, and drug administration metadata

The conversion of each protocol will handle the time alignment of data from each data stream.

## Figure-to-Data Stream Mapping

* **Figure 1**: Intracellular electrophysiology (VoltageClamp), Image Stacks and LineScans.
* **Figure 2**: Image stacks (spine density) and intracellular electrophysiology (CurrentClamp) data with optogenetic stimulation (Sr-oEPSC). The image stacks here are from two-photon laser scanning microscopy.
* **Figure 3**: Intracellular electrophysiology (VoltageClamp), LineScan data
* **Figure 4**: Image stacks (spine density) like figure 2 but they come in a different format as they come from confocal microscopy.
* **Figure 5**: Fluorescence imaging data from two-photon laser scanning microscopy.
* **Figure 6**: Intracellular electrophysiology (VoltageClamp), two-photon laser stacks for spine density and line scans.
* **Figure 7**: Intracellular electrophysiology (VoltageClamp), two-photon laser stacks for spine density, and line scans. Behavior scoring data (AIM) and behavioral videos.
* **Figure 8**: Behavioral scoring, intracellular electrophysiology (VoltageClamp), two-photon laser stacks for spine density, and line scans.

## Data Organization - Raw Data for Figs

### Data Files Overview

1. **Word File**: [Zhai et al. SA_manuscript_finalv2] - Most recent version of the submitted manuscript to Science Advances (adv8224) on 8 JAN 2025.
2. **PDF File**: [adv8224_SupplementalMaterial_v1] - Associated Supplemental Materials and Figures.
3. **Excel File**: [Key resources table_Zhai_v1] - Key resources table for relevant software (incl. Notebooks), protocols, antibodies, viruses, animals, chemicals, and hardware.
4. **Excel File**: [Data Connections] - Master sheets for raw data connections to panel figures, resources table, and some experiment metadata.
5. **Folder**: [Tabular dataset Zhai et al. 2025] - Includes readme file inside folder; spreadsheets for each Figure with separate pages for each Figure panel's data.
6. **Folder**: [Raw data for Figs] - Folders arranged by Figure panels.

## Key Devices and Software

### Acquisition Hardware
- **Bruker Ultima In Vitro Multiphoton Microscope System**: Two-photon imaging of ACh sensor
- **Chameleon Ultra II Laser (Coherent)**: Two-photon excitation laser
- **Nikon AXR Confocal Laser Microscope**: High-resolution confocal imaging
- **MultiClamp 700B amplifier (Axon Instrument)**: Electrophysiology recordings

### Software
- **PrairieView 5.3 (Bruker)**: Voltage protocols and data acquisition
- **AutoQuant X3.0.4**: For deconvolving two-photon images
- **Nikon NIS-Elements AR 5.41.02**: For de-noising and deconvolving high-resolution confocal images
- **NeuronStudio**: For semi-automated 3D reconstruction and spine counting from two-photon images
- **Imaris 10.0.0**: For detailed segmentation and 3D reconstruction of dendritic spines in confocal datasets

## General Acquisition Parameters

### Electrophysiology
> All the electrophysiological recordings were made using a MultiClamp 700B amplifier (Axon Instrument, USA), and signals were filtered at 2 kHz and digitized at 10 kHz. Voltage protocols and data acquisition were performed by PrairieView 5.3 (Bruker). The amplifier command voltage and all light source shutter and modulator signals were sent via the PCI-NI6713 analog-to-digital converter card (National Instruments, Austin, TX).

### Two-Photon Imaging
> The recorded SPN was visualized using 810 nm excitation laser (Chameleon Ultra II, Coherent, Santa Clara, USA). Dendritic structure was visualized by the red signal of Alexa Fluor 568 detected by a Hamamatsu R3982 side-on photomultiplier tube (PMT, 580-620 nm). Calcium transients, as signals in the green channel, were detected by a Hamamatsu H7422P-40 GaAsP PMT (490-560 nm, Hamamatsu Photonics, Japan).

## Analysis Scripts
* [Python script for analyzing Prairie View 5-generated .csv files of somatic excitability data](https://doi.org/10.5281/zenodo.14145776)
* [Python script for analyzing Prairie View 5-generated .csv files of dendritic excitability data](https://doi.org/10.5281/zenodo.14163502)
* [Python script for analyzing Prairie View 5-generated .csv files of BOT imaging data](https://doi.org/10.5281/zenodo.14163731)

## Related Figure-Specific Notes
- [Figure 1 Conversion Notes](figure_1_conversion_notes.md)
- [Figure 2 Conversion Notes](figure_2_conversion_notes.md)
- [Figure 3 Conversion Notes](figure_3_conversion_notes.md)
- [Figure 4 Conversion Notes](figure_4_conversion_notes.md)
- [Figure 5 Conversion Notes](figure_5_conversion_notes.md)
- [Figure 6 Conversion Notes](figure_6_conversion_notes.md)
- [Figure 7 Conversion Notes](figure_7_conversion_notes.md)
- [Figure 8 Conversion Notes](figure_8_conversion_notes.md)
- [Supplementary Figures Conversion Notes](supplementary_figures_conversion_notes.md)
