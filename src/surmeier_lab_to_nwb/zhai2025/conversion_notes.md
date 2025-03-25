# Zhai 2025 Conversion Notes

## Paper
Here is the paper:
https://www.biorxiv.org/content/10.1101/2025.01.02.631090v1.full

## Scope of Work (Private Access)
https://docs.google.com/document/d/1ITmVEEOQ1TbBC8hvZvkdQ-VcdW9_ynUWFwfRxBBrAGQ/edit?usp=sharing


### Some glossary for better reading the paper
* **DLS** dorsolateral striatum
* **levodopa** the treatment for Parkinson's disease
* **levodopa-induced dyskinesia (LID)** a set of erratic movements induced by levodopa after long-term treatment
* **rheobase** the minimal electric current required to excite a tissue (as nerve or muscle) given an indefinitely long time during which the current is applied,
* **SPNs** spiny projection neurons
* **dSPNs** direct pathway SPNs
* **iSPNs** indirect pathway SPNs
* ** AIM ** lower abnormal involuntary movement (AIM)
* **M1Rs** Muscarinic acetylcholine receptors. M1 receptor activation occurs when acetylcholine (the natural neurotransmitter) or other muscarinic agonists (such as muscarine, which is where these receptors get their name).
* **SCH** SCH23390, a pharmacological compound that acts as a selective antagonist at the dopamine D₁ receptor (D1R). In other words, SCH23390 blocks D1-type dopaminergic signaling in neurons.
* **2PLSM** Two Photon Laser Scanning Microscopy
* **SUL/sulpiride** sulpiride, a selective D2 dopamine receptor antagonist. By blocking D2-type receptors, sulpiride is used to test whether observed changes (in this case, decreased iSPN excitability in the on-state) depend on ongoing D2R signaling.
* **Ex vivo** refers to a procedure or experiment performed on tissue taken from a living organism, but studied outside that organism’s normal biological context (for example, in a chamber or petri dish)
* **oEPSC** optogenetically evoked postsynaptic current
* **quinpirole or DA** D2R agonist(s) (
* **GRABACh3.0** is a genetically encoded fluorescent sensor designed to detect acetylcholine (ACh) release. It is used to monitor ACh dynamics in real-time.
* **AAV5-hSyn-GRAB** refers to an adeno‐associated virus (serotype 5) that uses the human synapsin promoter to drive neuronal expression of the genetically encoded acetylcholine sensor GRABACh3.0. The authors injected this viral vector into the dorsolateral striatum so that they could optically monitor acetylcholine release in brain slices. This tool was critical for assessing how cholinergic signaling changes under conditions of dopamine depletion and during levodopa-induced dyskinesia.

* **Contralateral rotations** are a behavioral measure used in rodent models of Parkinson’s disease. When one side of the brain is lesioned (for example, with 6-OHDA), the imbalance in motor circuitry causes the animal to rotate predominantly toward the side opposite the lesion. By counting these rotations after a treatment like levodopa, researchers can assess changes in motor function and dyskinetic behavior.

* **CDGI** stands for CalDAG-GEFI, a calcium-activated guanine nucleotide exchange factor that is highly expressed in the striatum. It plays a critical role in linking M1 muscarinic receptor activation to intracellular signaling pathways that regulate dendritic excitability and synaptic plasticity in indirect pathway spiny projection neurons (iSPNs). In the paper, disrupting CDGI was found to blunt the dendritic adaptations in iSPNs and reduce dyskinetic behaviors, suggesting its important role in the pathophysiology of levodopa-induced dyskinesia.

Corroborate this:
The STK metadata in this context likely refers to metadata extracted from a MetaMorph STK file, a format used in microscopy for storing image sequences. STK files are essentially TIFF files with additional MetaMorph-specific metadata embedded in the tags. This metadata typically contains acquisition parameters, calibration details, and channel information.

These are tags appearing in the dendritic experiments of the Figures 1 and 3.

What is this Device:
 MultiClamp 700B




# Scope of Work (the conversion part)

## Conversion of Data Streams

Build interfaces to convert the following data streams to NWB format:

- Convert ABF format: I have not found .abf data
- Convert Bruker optical sensor recordings with appropriate metadata.
- Convert electrophysiological data acquired with the Bruker system
- Build an interface for manual segmentation data in Surmeier lab format.
- Integrate behavioral annotations from Surmeier lab custom format. Is this AIM scoring?
- Include behavioral video recordings. Are those
- Include electrical stimulation signals and metadata. I think this is on figure 5

Each conversion will utilize compression for efficient storage and implement chunking strategies to optimize for cloud storage. Detailed documentation will be provided for the installation and usage of the conversion software, including scripts to handle data from various protocols:
1. Ex-vivo brain slices with optogenetics: combine ABF recordings with optogenetic stimulus information.
2. Electrophysiology combined with optical imaging. Bruker optical sensor data with electrical recordings.
3. Two-photon laser scanning microscopy with electrical stimulation: process imaging data with stimulation events
4. Behavioral pharmacology: convert videos, tracking data, and drug administration metadata

The conversion of each protocol will handle the time alignment of data from each data stream.

## Which figures are related to which streams:
* Figure 1 contains recordings in xml data and the currents are on the folder organization. It also has some confocal imaging.
* Figure 2 has optogenetics and image stacks but no segmentation data as far as I can tell.
* Figure 3 also electrophysiology
* Figure 4 like figure 2 but also has "confocal spine density" data.
* Figure 5: should have fluoresence traces. This looks like it has two photon time
* Figure

# Data Organization - Raw data for Figs

## Data Files

1. **Word File**: [Zhai et al. SA_manuscript_finalv2] - Most recent version of the submitted manuscript to Science Advances (adv8224) on 8 JAN 2025.
2. **PDF File**: [adv8224_SupplementalMaterial_v1] - Associated Supplemental Materials and Figures.
3. **Excel File**: [Key resources table_Zhai_v1] - Key resources table for relevant software (incl. Notebooks), protocols, antibodies, viruses, animals, chemicals, and hardware.
4. **Excel File**: [Data Connections] - Master sheets for raw data connections to panel figures, resources table, and some experiment metadata.
5. **Folder**: [Tabular dataset Zhai et al. 2025] - Includes readme file inside folder; spreadsheets for each Figure with separate pages for each Figure panel’s data.
6. **Folder**: [Raw data for Figs] - Folders arranged by Figure panels.

## Figure 1

```
[256K]  .
├── [256K]  Dendritic excitability
│   ├── [256K]  LID off-state
│   ├── [256K]  LID on-state
│   └── [256K]  LID on-state with SCH
├── [256K]  Immunostaining images (Fig. 1B) confocal images
└── [256K]  Somatic excitability
    ├── [256K]  LID off-state
    ├── [256K]  LID on-state
    └── [256K]  LID on-state with SCH
```

### Somatic Excitability
It has a readme with the following info:

Read me:
Somatic excitability data were generated by recording voltage changes in response to current injection steps that are 500 ms in duration.
> Each folder with format “XXXXXXXX_X” represents data from a cell.
> Within the cell, “cell*-00*” represents a voltage recording trace.
> “cell*-001”: injection of -120pA current…
> “cell*-006”: injection of -20pA current…
> “cell*-007”: injection of 20pA current… and so forth.

15 directories on LID on-state

```
Command: 	 du -sh ./* | pbcopyfull
27M	./04112019_1
27M	./04122019_1
27M	./05012019_1
27M	./05032019_1
27M	./05072019_1
27M	./05072019_2
27M	./07062017_1
27M	./07062017_2
27M	./07072017_1
21M	./07072017_2
27M	./07192017_1
27M	./07192017_2
27M	./07202017_1
27M	./07202017_2
27M	./07212017_1
```

Let's look at one of them:

04112019_1
```
1,3M    ./cell1-001
1,3M    ./cell1-002
1,3M    ./cell1-003
1,3M    ./cell1-004
1,3M    ./cell1-005
1,3M    ./cell1-006
1,3M    ./cell1-007
1,3M    ./cell1-008
1,3M    ./cell1-009
1,3M    ./cell1-010
1,3M    ./cell1-011
1,3M    ./cell1-012
1,3M    ./cell1-013
1,3M    ./cell1-014
1,3M    ./cell1-015
1,3M    ./cell1-016
1,3M    ./cell1-017
1,3M    ./cell1-018
1,3M    ./cell1-019
1,3M    ./cell1-020
1,3M    ./cell1-021
```

And they look like this:

```bash
[256K]  .
├── [256K]  cell1-001
│   ├── [3.7K]  cell1-001_Cycle00001_VoltageOutput_001.xml
│   ├── [240K]  cell1-001_Cycle00001_VoltageRecording_001.csv
│   ├── [ 13K]  cell1-001_Cycle00001_VoltageRecording_001.xml
│   └── [ 554]  cell1-001.xml
├── [256K]  cell1-002
│   ├── [3.7K]  cell1-002_Cycle00001_VoltageOutput_001.xml
│   ├── [242K]  cell1-002_Cycle00001_VoltageRecording_001.csv
│   ├── [ 13K]  cell1-002_Cycle00001_VoltageRecording_001.xml
│   └── [ 554]  cell1-002.xml
├── [256K]  cell1-003
```

How do they look:

```xml
───────┬────────────────────────────────────────────────────────────────────────────────
       │ File: cell1-001.xml
───────┼────────────────────────────────────────────────────────────────────────────────
   1   │ <PVScan version="5.3.64.300" date="4/11/2019 3:17:18 PM" notes="">
   2   │   <Sequence type="VoltageRecording" cycle="1" time="15:17:18.1161356">
   3   │     <VoltageRecording name=" LIVE" triggerMode="None" cycle="1" index="1" confi
       │ gurationFile="cell1-001_Cycle00001_VoltageRecording_001.xml" dataFile="cell1-00
       │ 1_Cycle00001_VoltageRecording_001.csv" relativeTime="0" absoluteTime="0" />
   4   │     <VoltageOutput name="Step 1_-120 pA" triggerMode="None" filename="cell1-001
       │ _Cycle00001_VoltageOutput_001.xml" relativeTime="0" absoluteTime="0" />
   5   │   </Sequence>
   6   │ </PVScan>
```

- cell1-001_Cycle00001_VoltageOutput_001.xml
These is just a TimeSeries with two columns (time in ms) and Primary which from the xml should be in (mVThis is probably the )

- cell1-001_Cycle00001_VoltageRecording_001.xml
Most likely  XML describing how the raw data are acquired and scaled, which signals are recorded vs. disabled, and how the software’s user interface plots them.

- cell1-001_Cycle00001_VoltageOutput_001.xml

LID-off-state seems have the same structure.


### Dendritic Excitability
Missing a readme discussing the format but same top level structure as Somatic Excitability.

```bash
Command: 	 du -shc ./*
585M	./LID off-state
525M	./LID on-state
353M	./LID on-state with SCH
1,5G	total
```

LID on state

```bash
├── 0706a
│   ├── 0706a1
│   └── 0706a2
├── 0706b
│   └── 0706b1
├── 0707a
│   ├── 0707a1
│   └── 0707a2
├── 0707b
│   ├── 07072017_Cell2_dist2_trio-001
│   ├── 07072017_Cell2_dist2_trio-002
│   ├── 07072017_Cell2_dist2_trio-003
│   ├── 07072017_Cell2_prox2_trio-001
│   ├── 07072017_Cell2_prox2_trio-002
│   └── 07072017_Cell2_prox2_trio-003
├── 0719a
│   ├── 0719a1
│   └── 0719a2
├── 0719b
│   ├── 0719b1
│   └── 0719b2
├── 0720a
│   └── 0720a1
├── 0720b
│   ├── 0720b1
│   └── 0720b2
├── 0721a
│   ├── 0721a1
│   └── 0721a2
└── 0721b
    ├── 0721b1
    └── 0721b2
```
And the tree:


```bash
 heberto  …  LID on-state  0706a  0706a1  1  tree -shL 2
[256K]  .
├── [256K]  07062017_Cell1_dist1_trio-001
│   ├── [156K]  07062017_Cell1_dist1_trio-001_Cycle00001_Ch1_000001.ome.tif
│   ├── [ 33K]  07062017_Cell1_dist1_trio-001-Cycle00001_Ch1Source.tif
│   ├── [155K]  07062017_Cell1_dist1_trio-001_Cycle00001_Ch2_000001.ome.tif
│   ├── [ 33K]  07062017_Cell1_dist1_trio-001-Cycle00001_Ch2Source.tif
│   ├── [173K]  07062017_Cell1_dist1_trio-001_Cycle00001_LineProfileData.csv
│   ├── [5.4K]  07062017_Cell1_dist1_trio-001_Cycle00001_VoltageOutput_001.xml
│   ├── [997K]  07062017_Cell1_dist1_trio-001_Cycle00001_VoltageRecording_001.csv
│   ├── [ 13K]  07062017_Cell1_dist1_trio-001_Cycle00001_VoltageRecording_001.xml
│   ├── [286K]  07062017_Cell1_dist1_trio-001.env
│   ├── [6.4K]  07062017_Cell1_dist1_trio-001.xml
│   └── [256K]  References
├── [256K]  07062017_Cell1_dist1_trio-002
│   ├── [156K]  07062017_Cell1_dist1_trio-002_Cycle00001_Ch1_000001.ome.tif
│   ├── [ 33K]  07062017_Cell1_dist1_trio-002-Cycle00001_Ch1Source.tif
│   ├── [155K]  07062017_Cell1_dist1_trio-002_Cycle00001_Ch2_000001.ome.tif
│   ├── [ 33K]  07062017_Cell1_dist1_trio-002-Cycle00001_Ch2Source.tif
│   ├── [172K]  07062017_Cell1_dist1_trio-002_Cycle00001_LineProfileData.csv
│   ├── [5.4K]  07062017_Cell1_dist1_trio-002_Cycle00001_VoltageOutput_001.xml
│   ├── [996K]  07062017_Cell1_dist1_trio-002_Cycle00001_VoltageRecording_001.csv
│   ├── [ 13K]  07062017_Cell1_dist1_trio-002_Cycle00001_VoltageRecording_001.xml
│   ├── [286K]  07062017_Cell1_dist1_trio-002.env
│   ├── [6.4K]  07062017_Cell1_dist1_trio-002.xml
│   └── [256K]  References

```

There are tiffs here.

The references are also tiff fies

* Question: what are the references?

checking the metadata for one channel
```bash
  <Pixels DimensionOrder="XYZCT" ID="Pixels:1" PhysicalSizeX="0.194576" PhysicalSizeY="0.194576" PhysicalSizeZ="1" SizeC="2" SizeT="1" SizeX="31" SizeY="2500" SizeZ="1" Type="uint16">
```

Indicates that the tiffs at the bottom level of the tree are not time series. Images with two channels.

Question: what is ht

#### Figure 1 Someiatc Excitability Readme
Somatic excitability data were generated by recording voltage changes in response to current injection steps that are 500 ms in duration.

Each folder with format “XXXXXXXX_X” represents data from a cell.
Within the cell, “cell*-00*” represents a voltage recording trace.
“cell*-001”: injection of -120pA current…
“cell*-006”: injection of -20pA current…
“cell*-007”: injection of 20pA current… and so forth.


## Figure 2

The data of figure 2
```
.
├── Spine density
│   ├── control dSPN
│   ├── LID off-state dSPN
│   ├── LID on-state dSPN
│   └── PD dSPN  # PD most likely means Parkinson's disease
└── Sr-oEPSC
    ├── LID off-state
    └── LID on-state

### Spine Density
This is data that counts the density of the spines and is presented as image stacks.

```

LID on-state dSPN:

```bash
0411a2019  05012019a  0706a  0707a  0719a  0720a
04122019a  05032019a  0706b  0707b  0719b  0721a
```
Not clear why the folder have such a different name convention.

Here the bottom are confocal stacks:

```bash
.
├── 0411a2019
│   └── Decon_20190411_Cell1_dist1
├── 04122019a
│   └── Decon_20190412_Cell1_prox12
├── 05012019a
│   └── Decon_20190501_Cell1_prox12
├── 05032019a
│   ├── Decon-20190503_Cell1_dist12-001
│   └── Decon20190503_Cell1_Prox1-001
├── 0706a
│   ├── Decon_07062017_Cell1_dist1
│   ├── Decon_07062017_Cell1_dist2
│   ├── Decon_20170706_Cell1_prox12
│   └── Decon_20170706_Cell1_prox2
├── 0706b
│   ├── Decon_20170706_Cell2_dist1
│   ├── Decon_20170706_Cell2_dist2
│   ├── Decon_20170706_Cell2_prox1
│   └── Decon_20170706_Cell2_prox2
├── 0707a
│   ├── Decon_20170707_Cell1_dist1
│   ├── Decon_20170707_Cell1_dist2
│   ├── Decon_20170707_Cell1_prox1
│   └── Decon_20170707_Cell1_prox23
├── 0707b
│   ├── Decon_20170707_Cell2_prox12
│   └── Decon_20170707_Cell2_realprox12
```

Inside there seems to be stacks:

```bash
├── Decon_20170720_Cell1_dist1
│   ├── 20_ZSeries-20170720_Cell1_dist1-001_Cycle00001_Ch1_#.ome_Z001.tif
│   ├── 20_ZSeries-20170720_Cell1_dist1-001_Cycle00001_Ch1_#.ome_Z002.tif
│   ├── 20_ZSeries-20170720_Cell1_dist1-001_Cycle00001_Ch1_#.ome_Z003.tif
│   ├── 20_ZSeries-20170720_Cell1_dist1-001_Cycle00001_Ch1_#.ome_Z004.tif
│   ├── 20_ZSeries-20170720_Cell1_dist1-001_Cycle00001_Ch1_#.ome_Z005.tif
│   ├── 20_ZSeries-20170720_Cell1_dist1-001_Cycle00001_Ch1_#.ome_Z006.tif
│   ├── 20_ZSeries-20170720_Cell1_dist1-001_Cycle00001_Ch1_#.ome_Z007.tif
│   ├── 20_ZSeries-20170720_Cell1_dist1-001_Cycle00001_Ch1_#.ome_Z008.tif
│   ├── 20_ZSeries-20170720_Cell1_dist1-001_Cycle00001_Ch1_#.ome_Z009.tif
│   ├── 20_ZSeries-20170720_Cell1_dist1-001_Cycle00001_Ch1_#.ome_Z010.tif
│   ├── 20_ZSeries-20170720_Cell1_dist1-001_Cycle00001_Ch1_#.ome_Z011.tif
│   ├── 20_ZSeries-20170720_Cell1_dist1-001_Cycle00001_Ch1_#.ome_Z012.tif
│   ├── 20_ZSeries-20170720_Cell1_dist1-001_Cycle00001_Ch1_#.ome_Z013.tif

```
Some info on a specific stack
```bash
python -m tifffile 20_ZSeries-20170720_Cell1_dist1-001_Cycle00001_Ch1_#.ome_Z001.tif

Reading TIFF header: 0.316915 s
Reading image data: 0.000305 s
Generating report:   0.004349 s

TiffFile '20_ZSeries-20170720_Cell1_d…cle00001_Ch1_#.ome_Z001.tif'  129.46 KiB  mediacy

TiffPageSeries 0  256x256  uint16  YX  uniform  1 Pages  @8

TiffPage 0 @132276  256x256  uint16  minisblack memmappable  mediacy
```

### Sr-oEPSC

Here is the Optical Stimulation data

```bash
.
├── 07052023a
│   ├── cell1_LED14-001
│   ├── cell1_LED14-002
│   ├── cell1_LED14-003
│   ├── cell1_LED14-004
│   ├── cell1_LED14-005
│   ├── cell1_LED14-006
│   ├── cell1_LED14-007
│   ├── cell1_LED14-008
│   ├── cell1_LED14-009
│   └── cell1_LED14-010
├── 07052023b
│   ├── cell2_LED12-006
│   ├── cell2_LED12-007
│   ├── cell2_LED12-008
│   ├── cell2_LED12-009
│   ├── cell2_LED12-010
│   ├── cell2_LED12-011
│   ├── cell2_LED12-012
│   ├── cell2_LED12-013
│   ├── cell2_LED12-014
│   └── cell2_LED12-015
```

And inside there

```bash
├── cell1_LED14-001
│   ├── cell1_LED14-001_Cycle00001_VoltageOutput_001.xml
│   ├── cell1_LED14-001_Cycle00001_VoltageRecording_001.csv
│   ├── cell1_LED14-001_Cycle00001_VoltageRecording_001.xml
│   ├── cell1_LED14-001.env
│   ├── cell1_LED14-001.xml
│   └── References
├── cell1_LED14-002
│   ├── cell1_LED14-002_Cycle00001_VoltageOutput_001.xml
│   ├── cell1_LED14-002_Cycle00001_VoltageRecording_001.csv
│   ├── cell1_LED14-002_Cycle00001_VoltageRecording_001.xml
│   ├── cell1_LED14-002.env
│   ├── cell1_LED14-002.xml
```


### File Descriptions for Optical Stimulation Data

- **cell1_LED14-001_Cycle00001_VoltageOutput_001.xml**: XML file specifying the voltage output protocol for the experiment. It defines timing, amplitude, and waveforms for LED stimulation, current/voltage steps, or other output signals.

- **cell1_LED14-001_Cycle00001_VoltageRecording_001.csv**: CSV file containing raw time-series data for the recorded signals, such as membrane currents/voltages or auxiliary signals captured during the experiment.

- **cell1_LED14-001_Cycle00001_VoltageRecording_001.xml**: XML file describing metadata for the CSV recording, including sampling rates, scaling factors (e.g., pA, mV), channel mappings, and amplifier settings.

- **cell1_LED14-001.env**: Environment/configuration file containing information about the microscope or acquisition system settings (e.g., laser powers, PMT gains, etc.).

- **cell1_LED14-001.xml**: A higher-level XML file referencing the specific acquisition cycles and linking them to their respective output/recording metadata.

- **References (directory)**: A folder possibly containing calibration data or other supplementary files that further document the experimental setup or parameters.

## Figure 3

Very similar, also for dendritic experiments we have this at the bottom:

```bash
 05232016_Cell1_dist1_trio-001-Cycle00001_Ch1Source.tif              References
 05232016_Cell1_dist1_trio-001-Cycle00001_Ch2Source.tif
 05232016_Cell1_dist1_trio-001.env
 05232016_Cell1_dist1_trio-001.xml
 05232016_Cell1_dist1_trio-001_Cycle00001_Ch1_000001.ome.tif
 05232016_Cell1_dist1_trio-001_Cycle00001_Ch2_000001.ome.tif
 05232016_Cell1_dist1_trio-001_Cycle00001_LineProfileData.csv
 05232016_Cell1_dist1_trio-001_Cycle00001_VoltageOutput_001.xml
 05232016_Cell1_dist1_trio-001_Cycle00001_VoltageRecording_001.csv
 05232016_Cell1_dist1_trio-001_Cycle00001_VoltageRecording_001.xml
```

## Figure 4

```bash
.
├── Confocal spine density
│   ├── Fig 4H
│   ├── Fig 4I
│   └── Fig 4J and Suppl Fig 5
├── Spine density (These are like figure 2)
│   ├── control iSPN
│   ├── LID off-state iSPN
│   ├── LID on-state iSPN
│   └── PD iSPN
└── Sr-oEPSC (these are like figure 2)
    ├── LID off-state
    └── LID on-state
```

Extensions found in this figure are the most varied:

```
 find . -type f -name "*.*" | rev | cut -d. -f1 | rev | sort | uniq
bmp
csv
env
ims
jpg
log
lut
nd2
oif
png
pty
pzfx
roi
swc
tif
txt
xml
```


- **.oif files**
    Olympus Image Files – a proprietary format used by Olympus microscopy systems to store multi-dimensional (e.g., multi-channel, z-stack) raw image data.

- **.oif.files folder**
    For each `.oif` file, there is an accompanying folder (named `.oif.files`) that contains the individual components of the image. This folder holds the individual image slices, metadata, and auxiliary files related to that acquisition.

- **.roi files**
    Region of Interest files that store information about specific areas or contours drawn on the image for analysis or segmentation.

- **.pty files**
    These files likely contain parameter or pointer data associated with each image slice. They might store settings such as acquisition parameters, focus information, or metadata for that particular slice.

- **.tif files**
    TIFF (Tagged Image File Format) images – these are the actual image slices exported from the microscope system. They are widely used for storing high-quality microscopy data.

- **.lut files**
    Lookup Table files that define color maps or intensity scaling. They’re used by imaging software to render grayscale data in false color for better visualization.

- **.bmp files**
    Bitmap image files (in this case, a thumbnail image). These provide a quick, low-resolution preview of the image.

- **.log files (e.g., MATL_01_01.log)**
    Log files document processing steps, parameters, or any errors during image acquisition or subsequent processing, often generated by analysis software such as MATLAB.

- **TileConfiguration.txt**
    This text file likely contains information about how individual image tiles are arranged and stitched together to form a mosaic. It details the layout used during the image acquisition if a large area was imaged in tiles.

The raw data for Figure 4 comes from an Olympus system and includes the native multi-dimensional `.oif` files along with a folder of associated image slices (`.tif`), metadata (`.pty`), ROI definitions (`.roi`), and other auxiliary files (lookup tables, thumbnails, logs, and tile configuration). These files together allow for subsequent processing steps—such as deconvolution, segmentation, and 3D reconstruction—to accurately quantify dendritic spine density.

### Spine Density

```bash
tree -L 3 -I "*.xml|*.env|*.csv"
```

The only files are xml env and csv

### Confocal spine density

```bash
├── Fig 4H
│   ├── processed  # Contains tiffs and jpgs
│   └── raw
├── Fig 4I
│   ├── processed
│   └── raw
└── Fig 4J and Suppl Fig 5
    ├── processed-not organized yet, see if needed-denoise-decon,ims,excel files,code,combined excel
    ├── raw
    └── Spine density plots for Fig 4J and Suppl Fig 5.pzfx
```

What is pzfx?
Apparently [Graph Pad](https://www.graphpad.com/)

This is a proprietary format for Prism, a software used for scientific graphing and statistical analysis. The .pzfx file format is a compressed file that contains data, graphs, and analysis results created in GraphPad Prism.

```bash
├── Fig 4H
│   ├── processed
│   │   ├── 5104-3 60x str_Cycle fused
│   │   │   ├── Fused_5104-3 60x str_Cycle stitched.oif - C=0_Fused_5104-3 60x str_Cycle stitched.oif - C=0_Image_01_01_02_010000.tif
│   │   │   ├── Fused_5104-3 60x str_Cycle stitched.oif - C=0_Fused_5104-3 60x str_Cycle stitched.oif - C=0_Image_01_01_02_010001.tif
...
...
...
│   │   ├── MAX_Fused_5104-3 60x str_Cycle stitched.jpg
│   │   └── MAX_Fused_5104-3 60x str_Cycle stitched with scale bar.jpg
│   └── raw
│       └── 5104-3 60x str_Cycle
│           ├── Image_01_01_01_01.oif
│           ├── Image_01_01_01_01.oif.files
│           ├── Image_01_01_02_01.oif
│           ├── Image_01_01_02_01.oif.files
│           ├── Image_01_01_03_01.oif
│           ├── Image_01_01_03_01.oif.files
│           ├── Image_01_01_04_01.oif
│           ├── Image_01_01_04_01.oif.files
│           ├── Image_01_01_05_01.oif
│           ├── Image_01_01_05_01.oif.files
│           ├── Image_01_01_06_01.oif
│           ├── Image_01_01_06_01.oif.files
│           ├── Image_01_01_07_01.oif
│           ├── Image_01_01_07_01.oif.files
│           ├── Image_01_01_08_01.oif
│           ├── Image_01_01_08_01.oif.files
│           ├── Image_01_01_09_01.oif
│           ├── Image_01_01_09_01.oif.files
│           ├── MATL_01_01.log
│           └── TileConfiguration.txt
├── Fig 4I
│   ├── processed
│   │   ├── 6-OHDA
│   │   │   ├── 8040-slide 1-slice 2-cell 1 proxi - Denoised decon_2024-07-27T13-10-58.187.png
│   │   │   ├── 8040-slide 1-slice 2-cell 1 proxi - Denoised decon.ims
│   │   │   └── 8040-slide 1-slice 2-cell 1 proxi - Denoised decon.nd2
│   │   ├── control
│   │   │   ├── 3824-slide 2-slice 2-cell 2 proxi - Denoised decon_2024-11-15T13-58-34.120.png
│   │   │   ├── 3824-slide 2-slice 2-cell 2 proxi - Denoised decon.ims
│   │   │   └── 3824-slide 2-slice 2-cell 2 proxi - Denoised decon.nd2
│   │   ├── off-state
│   │   │   ├── 8041-slide 2-slice 3-cell 2 proxi - Denoised decon_2024-07-27T16-23-35.514.png
│   │   │   ├── 8041-slide 2-slice 3-cell 2 proxi - Denoised decon.ims
│   │   │   └── 8041-slide 2-slice 3-cell 2 proxi - Denoised decon.nd2
│   │   └── on-state
│   │       ├── 8939-slide 1-slice 3-cell 2 den2 proxi - Denoised decon_2024-07-28T13-53-52.753.png
│   │       ├── 8939-slide 1-slice 3-cell 2 den2 proxi - Denoised decon.ims
│   │       └── 8939-slide 1-slice 3-cell 2 den2 proxi - Denoised decon.nd2
│   └── raw
│       ├── 6-OHDA
│       │   └── 8040-slide 1-slice 2-cell 1 proxi.nd2
│       ├── control
│       │   └── 3824-slide 2-slice 2-cell 2 proxi.nd2
│       ├── off-state
│       │   └── 8041-slide 2-slice 3-cell 2 proxi.nd2
│       └── on-state
│           └── 8939-slide 1-slice 3-cell 2 den2 proxi.nd2
└── Fig 4J and Suppl Fig 5
    ├── processed-not organized yet, see if needed-denoise-decon,ims,excel files,code,combined excel
    │   ├── 6-OHDA
    │   ├── control
    │   ├── off-state
    │   └── on-state
    ├── raw
    │   ├── 6-OHDA
    │   │   ├── 3823-slide 2-slice 4-cell 2 proxi.nd2
    │   │   ├── 3825-slide 2-slice 3-cell 1 proxi.nd2
    │   │   ├── 3825-slide 2-slice 3-cell 3 proxi.nd2
    │   │   ├── 8040-slide 1-slice 2-cell 1 proxi.nd2
    │   │   ├── 8040-slide 1-slice 2-cell 2 proxi.nd2
    │   │   ├── 8040-slide 1-slice 2-cell 3 proxi.nd2
    │   │   ├── 8040-slide 1-slice 3-cell 1 proxi.nd2
    │   │   ├── 8040-slide 1-slice 3-cell 2 den2 proxi.nd2
    │   │   ├── 8040-slide 1-slice 3-cell 2 proxi.nd2
    │   │   ├── 8941-slide 1-slice 2-cell 1 den 2 proxi.nd2
    │   │   ├── 8941-slide 1-slice 2-cell 1 proxi.nd2
    │   │   ├── 8941-slide 1-slice 3-cell 1 proxi.nd2
    │   │   ├── 8941-slide 1-slice 3-cell 2 proxi.nd2
    │   │   ├── 8941-slide 1-slice 3-cell 3 proxi.nd2
    │   │   ├── 8941-slide 1-slice 3-cell 4 den2 proxi.nd2
    │   │   ├── 8941-slide 1-slice 3-cell 4 proxi.nd2
    │   │   ├── 8941-slide 2-slice 2-cell 1 proxi.nd2
    │   │   ├── 8941-slide 2-slice 2-cell 2 proxi.nd2
    │   │   ├── 8941-slide 2-slice 2-cell 3 proxi.nd2
    │   │   ├── 8941-slide 2-slice 2-cell 4 proxi.nd2
    │   │   ├── 8941-slide 2-slice 2-cell 5 den2 proxi.nd2
    │   │   └── 8941-slide 2-slice 2-cell 5 proxi.nd2
    │   ├── control
    │   │   ├── 3824-slide 2-slice 2-cell 1 proxi.nd2
    │   │   ├── 3824-slide 2-slice 2-cell 2 proxi.nd2
    │   │   ├── 5104-slice 3-cell 1 proxi3.nd2
    │   │   ├── 5104-slide 2-slice 2-cell 1 den 2 proxi.nd2
    │   │   ├── 5104-slide 2-slice 2-cell 1 proxi.nd2
    │   │   ├── 5104-slide 2-slice 2-cell 2 proxi.nd2
    │   │   ├── 8038-slide 1-slice 2-cell 1 den2 proxi.nd2
    │   │   ├── 8038-slide 1-slice 2-cell 1 proxi_0001.nd2
    │   │   ├── 8038-slide 1-slice 2-cell 2 den2 proxi.nd2
    │   │   ├── 8038-slide 1-slice 2-cell 2 proxi.nd2
    │   │   ├── 8038-slide 1-slice 2-cell 3 proxi.nd2
    │   │   ├── 8940-slide 1-slice 1-cell 1 den2 proxi.nd2
    │   │   ├── 8940-slide 1-slice 1-cell 1 proxi.nd2
    │   │   ├── 8940-slide 1-slice 2-cell 1 den2 proxi.nd2
    │   │   ├── 8940-slide 1-slice 2-cell 1 proxi.nd2
    │   │   ├── 8940-slide 1-slice 2-cell 2 den2 proxi.nd2
    │   │   ├── 8940-slide 1-slice 2-cell 2 proxi.nd2
    │   │   └── 8940-slide 1-slice 2-cell 3 proxi.nd2
    │   ├── off-state
    │   │   ├── 3826-slide 2-slice 3-cell 1 proxi.nd2
    │   │   ├── 3826-slide 2-slice 3-cell 2 proxi.nd2
    │   │   ├── 5105-slice 3-cell 1 proxi3.nd2
    │   │   ├── 5105-slice 3-cell 2 den2 proxi.nd2
    │   │   ├── 5105-slice 3-cell 2 proxi.nd2
    │   │   ├── 5105-slide 2-slice 2-cell 1 den 2 proxi.nd2
    │   │   ├── 5105-slide 2-slice 2-cell 1 proxi.nd2
    │   │   ├── 5105-slide 2-slice 3-cell 1 proxi.nd2
    │   │   ├── 5105-slide 2-slice 3-cell 2 den 2 proxi.nd2
    │   │   ├── 5105-slide 2-slice 3-cell 2 proxi.nd2
    │   │   ├── 5105-slide 2-slice 3-cell 3 proxi.nd2
    │   │   ├── 8041-slide 1-slice 4-cell 1 proxi.nd2
    │   │   ├── 8041-slide 2-slice 2-cell 1 proxi.nd2
    │   │   ├── 8041-slide 2-slice 3-cell 1den2 proxi.nd2
    │   │   ├── 8041-slide 2-slice 3-cell 1 proxi.nd2
    │   │   ├── 8041-slide 2-slice 3-cell 2 proxi.nd2
    │   │   ├── 8041-slide 3-slice 3-cell 1 proxi.nd2
    │   │   ├── 8041-slide 3-slice 3-cell 3 proxi.nd2
    │   │   ├── 8041-slide 3-slice 3-cell 4 den2 proxi.nd2
    │   │   ├── 8041-slide 4-slice 3-cell 2 proxi.nd2
    │   │   ├── 8943-slide 1-slice 1-cell 1 den2 proxi.nd2
    │   │   ├── 8943-slide 1-slice 1-cell 1 proxi.nd2
    │   │   ├── 8943-slide 1-slice 1-cell 2 proxi.nd2
    │   │   ├── 8943-slide 1-slice 2-cell 1 den2 proxi.nd2
    │   │   └── 8943-slide 1-slice 2-cell 1 proxi.nd2
    │   └── on-state
    │       ├── 5107-slide 2-slice 2-cell 1 den 2 proxi.nd2
    │       ├── 5107-slide 2-slice 2-cell 1 proxi.nd2
    │       ├── 5107-slide 2-slice 2-cell 2 proxi.nd2
    │       ├── 5107-slide 2-slice 3-cell 3 proxi.nd2
    │       ├── 5108-slice 2-cell 1 den2 proxi.nd2
    │       ├── 5108-slice 2-cell 1 den3 proxi.nd2
    │       ├── 5108-slice 2-cell 1 proxi.nd2
    │       ├── 5108-slide 2-slice 1-cell 1 den2 proxi.nd2
    │       ├── 5108-slide 2-slice 1-cell 1 proxi.nd2
    │       ├── 5108-slide 2-slice 1-cell 2 proxi.nd2
    │       ├── 5108-slide 2-slice 2-cell 1 den 2 proxi.nd2
    │       ├── 8042-slide 1-slice 2-cell 1 proxi.nd2
    │       ├── 8042-slide 1-slice 2-cell 2 den2 proxi.nd2
    │       ├── 8042-slide 1-slice 2-cell 2 proxi.nd2
    │       ├── 8042-slide 1-slice 2-cell 3 proxi.nd2
    │       ├── 8042-slide 1-slice 2-cell 4 proxi.nd2
    │       ├── 8939-slide 1-slice 2-cell 1 den2 proxi.nd2
    │       ├── 8939-slide 1-slice 2-cell 1 proxi.nd2
    │       ├── 8939-slide 1-slice 2-cell 2 proxi.nd2
    │       ├── 8939-slide 1-slice 3-cell 1 proxi.nd2
    │       ├── 8939-slide 1-slice 3-cell 2 den2 proxi.nd2
    │       └── 8939-slide 1-slice 3-cell 2 proxi.nd2
    └── Spine density plots for Fig 4J and Suppl Fig 5.pzfx
```

What are:
* nd2 files:  raw image data files produced by Nikon imaging system?
* ims files: Imari software?
* oif files: Olympus Microscope?
* pyt files:


## Figure 5

```bash
.
├── LID off (Figure 5 E)
├── PD (6-OHDA) (Figure 5 D)
└── UL control (Figure 5 C)
```
Note that each of them correspond to a paper figure

• Quinpirole is a D2 receptor agonist. When applied, it activates D2 receptors on cholinergic interneurons, leading to a suppression of ACh release.

• +DA indicates the bath application of dopamine. Adding dopamine mimics the natural ligand's effect—activating dopamine receptors (primarily D2 receptors on ChIs) and thus inhibiting ACh release.

• +Sulpiride refers to the application of a D2 receptor antagonist. By blocking D2 receptors, sulpiride prevents dopamine from exerting its inhibitory effect on ACh release, resulting in an increased evoked ACh signal.


```bash
.
├── LID off
│   ├── 05242024slice1
│   ├── 05242024slice2
│   ├── 05292024slice1
│   ├── 05292024slice2
│   ├── 05312024slice1
│   ├── 05312024slice2
│   ├── 06032024slice1
│   ├── 06042024slice1
│   ├── 06042024slice2
│   ├── 06042024slice3
│   └── 06042024slice4
├── PD
│   ├── 04122024slice1ROI1
│   ├── 04122024slice1ROI2
│   ├── 04122024slice2ROI1
│   ├── 04122024slice2ROI2
│   ├── 04162024slice1ROI1
│   ├── 04162024slice1ROI2
│   ├── 04162024slice2ROI1
│   ├── 04162024slice2ROI2
│   ├── 05092024slice1
│   └── 05092024slice2
└── UL control
    ├── 04022024slice1ROI1
    ├── 04022024slice1ROI2
    ├── 04022024slice2ROI1
    ├── 04022024slice2ROI2
    ├── 04022024slice3ROI1
    ├── 04022024slice3ROI2
    ├── 04052024slice1ROI1
    ├── 04052024slice1ROI2
    ├── 04052024slice2ROI1
    ├── 04102024slice1ROI1
    ├── 04102024slice1ROI2
    ├── 04102024slice2ROI1
    └── 04102024slice2ROI2
```

These folders contain what looks like Two Photon time series data:

```bash
├── BOT_04122024_slice1ROI1_sul_single-002
│   ├── BOT_04122024_slice1ROI1_sul_single-002_Cycle00001-botData.csv
│   ├── BOT_04122024_slice1ROI1_sul_single-002_Cycle00001_Ch2_000001.ome.tif
│   ├── BOT_04122024_slice1ROI1_sul_single-002_Cycle00001_Ch2_000002.ome.tif
│   ├── BOT_04122024_slice1ROI1_sul_single-002_Cycle00001_Ch2_000003.ome.tif
│   ├── BOT_04122024_slice1ROI1_sul_single-002_Cycle00001_Ch2_000004.ome.tif
│   ├── BOT_04122024_slice1ROI1_sul_single-002_Cycle00001_Ch2_000005.ome.tif
│   ├── BOT_04122024_slice1ROI1_sul_single-002_Cycle00001_Ch2_000006.ome.tif
│   ├── BOT_04122024_slice1ROI1_sul_single-002_Cycle00001_Ch2_000007.ome.tif
```

## Figure 6

```bash
 heberto  (e) work  …  LID_paper_Zhai_2025  Raw data for Figs  Figure 6  1  tree -L 2
.
├── Dendritic excitability
│   ├── control
│   └── M1R antagonist
├── Somatic excitability
│   ├── control
│   └── M1R antagonist
└── Spine density
    ├── control
    └── M1R antagonist

```

## Figure 7

```bash

.
├── AIM rating
│   └── AIM testing_CDGI KO.xlsx  # This is an excel file
├── contralateral rotations
│   └── CDGI KO videos  # These are the videos
├── KO DE on vs off
│   ├── KO off-state
│   └── KO on-state
├── KO SE on vs off
│   ├── KO off-state
│   └── KO on-state
├── KO spine density
│   ├── KO off
│   └── KO on
└── oxoM on DE
    ├── KO
    ├── ReadMe.rtf
    └── WT
```
How does the excel looks like:

![aim score](./assets/aim_score.png)

## Figure 8

```bash
 heberto  (e) work  …  LID_paper_Zhai_2025  Raw data for Figs  Figure 8  1  tree -L 2
.
├── M1R CRISPR AIMs
│   └── AIM raw score_M1R CRISPR.xlsx
├── M1R CRISPR SE
│   ├── interleaved control
│   └── M1R CRISPR
└── M1R CRISPR spine density
    ├── control
    └── M1R CRISPR
```

### M1R CRISPR SE

```bash
tree -L 3 -I "*.xml|*.env|*.csv"
```
Only contains this three types of files that look like this:

```bash
.
├── cell1-001
│   ├── cell1-001_Cycle00001_VoltageOutput_001.xml
│   ├── cell1-001_Cycle00001_VoltageRecording_001.csv
│   ├── cell1-001_Cycle00001_VoltageRecording_001.xml
│   ├── cell1-001.env
│   ├── cell1-001.xml
│   └── References

```


### M1R CRISPR spine density
I expected to see spine density calculations on the M1R CRISPR spie density but there is only xml files

```bash
tree -L 4 -I "*.tif*"
.
├── control
│   ├── 20221012a
│   │   ├── Decon_20221012_cell1_dist1
│   │   ├── Decon_20221012_cell1_dist2
│   │   ├── Decon_20221012_cell1_prox1
│   │   └── Decon_20221012_cell1_prox2
│   ├── 20221012b
│   │   ├── Decon_20221012_cell2_dist1
│   │   ├── Decon_20221012_cell2_dist2
│   │   ├── Decon_20221012_cell2_prox1
│   │   └── Decon_20221012_cell2_prox23
│   ├── 20221012c
│   │   ├── Decon_20221012_cell3_dist1
│   │   ├── Decon_20221012_cell3_prox1
│   │   └── Decon_20221012_cell3_prox2
│   ├── 20221130a
│   │   ├── Decon_20221130_cell1_dist1
│   │   ├── Decon_20221130_cell1-dist2
│   │   ├── Decon_20221130_cell1_prox1
│   │   └── Decon_20221130_cell1_prox2
│   ├── 20221130b
│   │   ├── Decon_20221130_cell2_dist12
│   │   ├── Decon_20221130_cell2_prox1
│   │   └── Decon_20221130_cell2_prox2
│   └── 20221130c
│       ├── Decon_20221130_cell3_dist1
│       ├── Decon_20221130_cell3_dist2
│       └── Decon_20221130_cell3_prox12
└── M1R CRISPR
    ├── 20221004a
    │   ├── Decon_20221004_cell1_dist1
    │   ├── Decon_20221004_Cell1_medium12
    │   └── Decon_20221004_cell1_prox1
    ├── 20221004b
    │   ├── Decon_20221004_cell2_dist1
    │   ├── Decon_20221004_cell2_dist2
    │   ├── Decon_20221004_cell2_prox1
    │   └── Decon_20221004_cell2_prox23
    ├── 20221004c
    │   ├── Decon_20221004_cell3_dist1
    │   ├── Decon_20221004_cell3_dist2
    │   ├── Decon_20221004_cell3_prox1
    │   └── Decon-20221004_cell3_prox2
    ├── 20221004d
    │   ├── Decon_20221004_cell4_dist1
    │   └── Decon_20221004_cell4_prox1
    ├── 20221005a
    │   ├── Decon_20221004_cell1_dist1
    │   ├── Decon_20221004_cell1_dist2
    │   ├── Decon-20221005_cell1_prox1
    │   └── Decon_20221005_cell1_prox2
    ├── 20221005b
    │   ├── Decon-20221005_cell2_dist1
    │   └── Decon_20221005_cell2_prox1
    ├── 20221005c
    │   ├── Decon_20221005_cell3_dist1
    │   └── Decon-20221005_cell3_prox12
    ├── 20221129a
    │   ├── Decon_20221129_cell1_dist1
    │   ├── Decon_20221129_cell1_dist2
    │   ├── Decon_20221129_cell1_prox12
    │   └── Decon_20221129_cell1_prox3
    ├── 20221129b
    │   ├── Decon_20221129_cell2_dist1
    │   ├── Decon_20221129_cell2_dist2
    │   ├── Decon_20221129_cell2_prox1
    │   └── Decon_20221129_cell2_prox2
    └── 20221129c
        ├── Decon-20221129_cell3_dist12
        ├── Decon_20221129_cell3_prox1
        └── Decon_20221129_cell3_prox2
```

## Supplementary Figure 3

This seems to contain a lot of videos:


```bash
.
└── M1R CRISPR videos
    ├── 01172024-
    │   ├── 1st session
    │   ├── 2nd session
    │   ├── 3rd session
    │   ├── 4th session
    │   └── 5th session
    ├── 02232024-
    │   ├── 1st session
    │   ├── 2nd session
    │   ├── 3rd session
    │   ├── 4th session
    │   └── 5th session
    ├── 05032024-
    │   ├── 1st session
    │   ├── 2nd session
    │   ├── 3rd session
    │   ├── 4th session
    │   └── 5th session
    ├── 09272023-
    │   ├── 1st session
    │   ├── 2nd session
    │   ├── 3rd session
    │   ├── 4th session
    │   └── 5th session
    └── 12082023-
        ├── 1st session
        ├── 2nd session
        ├── 3rd session
        ├── 4th session
        └── 5th session
```

# Data Connections
This is a table shared by the authors:

| TECHNIQUE                   | FIG     | CONT      | PD        | LID OFF   | LID OFF  | LID ON    | ON + agon | OXO     | OXO     | MSN    | MOUSE           | Age     | Virus Inj  | Ldopa   | Slicing | Patching  | Imaging  | Stimulation                     | Analysis               | Analysis        | Analysis       | Analysis       | Analysis         |
|-----------------------------|---------|-----------|-----------|-----------|----------|-----------|-----------|---------|---------|--------|----------------|---------|------------|---------|---------|-----------|----------|--------------------------------|------------------------|----------------|---------------|---------------|-----------------|
| Somatic Excitability        | 1F      | -         | -         | 15(8)8    | -        | 15(10)10  | 10(5)5    | -       | -       | D1     | tdTomato       | 7-12 wks | -          | -       | I Clamp | 2P GEI    | Neg120 to 300pA, 0.5s | pA separate traces       | 14145776       | RMP [5, 195]ms | [200, 700]ms   | #APs in 0.5s     |
| Somatic Excitability        | 3D      | -         | -         | 12(5)5    | -        | 11(8)8    | 7(4)4     | -       | -       | D2     | eGFP           | 7-12 wks | -          | -       | I Clamp | 2P GEI    | Neg120 to 300pA, 0.5s | -                      | 14145776       | Not shown       | -               | #APs in 0.5s     |
| VEHICLE CONTROL FOR FIG6    | 6C      | 8(4)4     | -         | 12(5)5    | 10(5)5   | -         | -         | -       | -       | D2     | eGFP           | 7-12 wks | M1R ant    | -       | I Clamp | 2P GEI    | Neg120 to 300pA, 0.5s | -                      | 14145776       | -               | -               | #APs in 0.5s     |
| VEHICLE CONTROL FOR FIG6    | 7C      | -         | -         | 13(5)5    | -        | 8(4)4     | -         | -       | -       | D2     | CDGIko X eGFP  | ~8 wks   | -          | -       | I Clamp | 2P GEI    | Neg120 to 300pA, 0.5s | -                      | 14145776       | -               | -               | #APs in 0.5s     |
| VEHICLE CONTROL FOR FIG8    | 8D      | 6(3)2     | -         | -         | -        | 8(4)3     | -         | -       | -       | D2     | Adora2-Cre     | ~8 wks   | 3 viruses  | -       | I Clamp | 2P GEI    | Neg120 to 300pA, 0.5s | -                      | 14145776       | -               | -               | #APs in 0.5s     |
| VEHICLE CONTROL FOR FIG8    | SF4     | -         | -         | -         | -        | -         | 7(3)3     | 6(3)3   | -       | D2     | Adora2-Cre     | ~8 wks   | 3 viruses  | -       | I Clamp | 2P GEI    | Neg120 to 300pA, 0.5s | -                      | 14145776       | -               | -               | #APs in 0.5s     |
| Dendritic Excitability      | 1I      | -         | -         | 10(7)7    | -        | 10(5)5    | 9(9)5     | -       | -       | D1     | tdTomato       | 7-12 wks | -          | -       | I Clamp | 2P Fluo-4 | 2nA, 2ms X3 @50Hz      | 4r3I294nqv1y/v1        | 14163502       | F0 [0.5, 0.95]s | F [1.0, 1.5]s  | (G-G0)/G0R0 AUC  |
| Dendritic Excitability      | 3F      | -         | -         | 10(5)5    | -        | 10(5)5    | 8(6)6     | -       | -       | D2     | eGFP           | 7-12 wks | -          | -       | I Clamp | 2P Fluo-4 | 2nA, 2ms X3 @50Hz      | 4r3I294nqv1y/v1        | 14163502       | Dist/Prox        | Same dendrite   | (G-G0)/G0R0 AUC  |
| ADD 25 apr 24, REM PT 25 FEB 20 | 6D  | 9(5)5     | -         | -         | -        | 7(4)4     | -         | -       | -       | D2     | eGFP           | 7-12 wks | M1R ant    | -       | I Clamp | 2P Fluo-4 | 2nA, 2ms X3 @50Hz      | 4r3I294nqv1y/v1        | 14163502       | G0/R0 ratio      | -               | (G-G0)/G0R0 AUC  |
| Spine Density               | 2BC     | 10(6)6    | 11(5)5    | 13(8)8    | -        | 8(5)5     | 10d, 11p  | -       | -       | D1     | tdTomato       | 7-12 wks | -          | -       | V Clamp | 2P Alexa  | Fluor feature morph   | Autoquant               | -              | -               | -               | #spines/um       |
| Presynaptic Excitability    | 2H      | SF1A      | -         | 10(4)3    | -        | 8(4)4     | -         | -       | -       | D1     | tdTomato       | 7-12 wks | ChR2-eYFP  | -       | V Clamp | 2P GEI    | ?%, 0.3ms (LED), 10 runs | 5qpvokx914o/v1        | TaroTools      | V0 [0, 20]ms   | 40 to 400ms     | EPSCs >5 stdev   |
| Ach Dynamics (BOT)          | 5F      | SF2       | 13/7(3)3  | 10/6(3)3  | 11/11(5)5 | -         | -         | -       | -       | D1 & D2 | C57bl/6J       | ~8 wks   | ACh3.0     | -       | 2P GEI  | 0.3mA, 1ms (electrode) | 21fps                    | 14163731       | F0 [9,10]s      | F[10,12]s       | (G-G0)/G0 AUC    |
| AIM Score                   | 7I [7J] | SF3A      | (8) [10]  | -         | -        | (9) [11]  | -         | -       | -       | D2     | CDGIko X eGFP  | ~8 wks   | -          | -       | -       | -         | -                             | 5jy18dm8rg2w/v1        | Spreadsheet    | -               | -               | AIM scores       |

## Working comments:

How to look for files:
```bash
find . -type f \( -iname "*.mp4" -o -iname "*.mov" -o -iname "*.avi" -o -iname "*.mkv" \)
```
```bash
ls -R | rg ".MP4"
```

###  Acquisition

> All the electrophysiological recordings were made using a MultiClamp 700B amplifier (Axon Instrument, USA), and signals were filtered at 2 kHz and digitized at 10 kHz. Voltage protocols and data acquisition were performed by PraireView 5.3 (Bruker). The amplifier command voltage and all light source shutter and modulator signals were sent via the PCI-NI6713 analog-to-digital converter card (National Instruments, Austin, TX).


[Ultima In Vitro Multiphoton
Microscope system](https://www.bruker.com/en/products-and-solutions/fluorescence-microscopy/multiphoton-microscopes/ultima-in-vitro.html)


#### Spine Density Calculation

> For assessment of dendritic spine density, images of dendritic segments (proximal: ~40 μm from soma; distal: > 80 μm from soma) were acquired with 0.15 μm pixels with 0.3 μm z-steps.
Images were deconvolved in AutoQuant X3.0.4 (MediaCybernetics, Rockville, MD) and semi-automated spine counting was performed using 3D reconstructions in NeuronStudio (CNIC, Mount Sinai School of Medicine, New York). On average, two proximal and two distal dendrites were imaged and analyzed per neuron.

But then there is also this section

#### High resolution confocal microscopy of sparsely labeled neurons

A Nikon AXR confocal laser microscope from the Center for Advanced Microscopy & Nikon Imaging Center at Northwestern University was used to image dendritic spines from sparsely labeled iSPNs in fixed brain sections. Z-stack images were acquired with a 60x oil immersion objective (NA = 1.49) at 0.125 μm intervals with a 0.09 μm pixel size. The images were de-noised and deconvolved using Nikon NIS-Elements AR 5.41.02 software.

Dendritic segments (~30 μm in length, located ~30 μm from the soma) were analyzed for spine density using Imaris 10.0.0 software (Oxford Instruments). The dendrite of interest was isolated by creating a new Surface. The ‘Labkit for pixel classification’ tool was used to train signal detection over noise, creating a tight surface. This surface was further refined to include only the dendritic segment of interest. A Filament for the dendritic segment was generated with the aid of an embedded supervised learning function. Spine detection was performed with the following parameters: thinnest spine head set at 0.188 μm, maximum spine length set at 5 μm, and no branched spines allowed. Spines were further refined using an embedded supervised learning function and manually corrected if necessary.


### Amplitudes of Sr2-oEPSC
> Amplitudes of Sr2-oEPSCs were
analyzed automatically using TaroTools Event Detection in Igor Pro 8 (WaveMetrics, Portland,
Oregon) followed by visual verification. Events were measured between 40 ms to 400 ms after


https://www.wavemetrics.com/forum/general/igorpro-neuromatic-event-detection

### Microscope

> Fluorescence was imaged using an Ultima In Vitro Multiphoton Microscope system (Bruker, Billerica, MA) with an Olympus 60x/0.9 NA water-immersion objective lens and a Hamamatsu H7422P-40 GaAsP PMT (490 nm to 560 nm, Hamamatsu Photonics,
