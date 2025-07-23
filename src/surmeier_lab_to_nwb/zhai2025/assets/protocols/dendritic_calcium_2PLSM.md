# Ex vivo imaging of dendritic calcium transients and spines with two-photon laser scanning microscopy (2PLSM)

**Date:** December 20, 2024

**Forked from:** Ex vivo imaging of genetically encoded sensors with two-photon laser scanning microscopy (2PLSM)

**DOI:** [dx.doi.org/10.17504/protocols.io.4r3l294nqv1y/v1](https://dx.doi.org/10.17504/protocols.io.4r3l294nqv1y/v1)

## Authors

**Enrico Zampese¹, Shenyu Zhai¹**

¹Northwestern University, Feinberg School of Medicine

**Corresponding Author:** Shenyu Zhai, Northwestern University

## Citation

**Protocol Citation:** Enrico Zampese, Shenyu Zhai 2024. Ex vivo imaging of dendritic calcium transients and spines with two-photon laser scanning microscopy (2PLSM). protocols.io https://dx.doi.org/10.17504/protocols.io.4r3l294nqv1y/v1

## License

This is an open access protocol distributed under the terms of the Creative Commons Attribution License, which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.

## Protocol Information

- **Status:** Working - We use this protocol and it's working
- **Created:** November 21, 2024
- **Last Modified:** December 20, 2024
- **Protocol ID:** 112567

## Keywords

2P, ex vivo slices, two-photon, imaging, dendritic spines, spiny projection neurons, SPNs, imaging of dendritic calcium transient, imaging dendritic ca2, dendritic calcium transient, dendritic spines of striatal spiny projection neuron, dendritic spine, striatal spiny projection neuron, ex vivo brain slice, scanning microscopy, imaging, spine, photon laser, ca2

## Funding

**Aligning Science Across Parkinson's (ASAP)**
Grant ID: ASAP-020551

## Abstract

This protocol is adapted from: [http://dx.doi.org/10.17504/protocols.io.n92ldry68g5b/v1](http://dx.doi.org/10.17504/protocols.io.n92ldry68g5b/v1) and modified for imaging dendritic Ca²⁺ transients and dendritic spines of striatal spiny projection neurons (SPNs) in ex vivo brain slices with two-photon laser scanning microscopy (2PLSM).

## Guidelines

Animal use procedures were approved by the Northwestern Institutional Animal Care and Use committee.

## Materials

### Equipment and Software
- 2PLSM optical workstation and computer with data acquisition software
- Artificial cerebrospinal fluid (aCSF)
- Internal solution
- Holding chamber
- Carbogen (95% O₂/5% CO₂) tank connected to regulator and bubblers
- Peristaltic pump with tubing and connectors, including inlet and outlet to microscope's imaging chamber
- Waste solution collector
- Patch pipettes (3-4 MΩ resistance) pulled from thick-wall borosilicate glass capillaries

### Reagents
- 2 mM stock solution of Fluo-4 (Thermo Fisher Scientific, F14200), dissolved in internal solution
- 1 mM stock solution of Alexa Fluor 568 hydrazide (Thermo Fisher Scientific, A10437), dissolved in internal solution
- Syringe filters (e.g. Millipore SLGV013SL)
- Sterile 1 ml syringes without needles (e.g. BD 309659)
- MicroFil Pipette Filler
- Slice hold-down (e.g. Warner Instrument, 64-1418)
- Ex vivo brain slices in which iSPNs or dSPNs are fluorescently labeled

### 2PLSM Optical Workstation and Computer Software

The laser scanning optical workstation embodies an Ultima dual-excitation-channel scan head (Bruker Nano Fluorescence Microscopy Unit). The foundation of the system is the Olympus BX-51WIF upright microscope with a 60X/0.9NA water-dipping objective lens. The automation of the XY stage motion, lens focus, and manipulator XYZ movement was provided by FM-380 shifting stage, axial focus module for Olympus scopes, and manipulators (Luigs & Neumann).

Cell visualization and patching were made possible by a ½" CCD video camera (Hitachi) imaged through a Dodt contrast tube, a 2x magnification changer (Bruker), and MicroManager software. Electrophysiological signals were sent and collected with a 700B patch clamp amplifier and MultiClamp Commander software with computer input and output signals were controlled by Prairie View 5.3-5.5 using a National Instruments PCI6713 output card and PCI6052e input card.

A two-photon laser (Chameleon Ultra II, Coherent, Santa Clara, CA) is used to excite the fluorescent dyes (Alexa and Fluo4). The excitation wavelength 810 nm is selected based on the dyes being imaged. Laser power attenuation is achieved with two Pockels' cell electro-optic modulators (models M350-80-02-BK and M350-50-02-BK, Con Optics) controlled by Prairie View 5.3–5.5 software. The two modulators are aligned in series to provide enhanced modulation range for fine control of the excitation dose, to limit the sample maximum power, and to serve as a rapid shutter during line-scan acquisitions.

Fluorescence is imaged using an Ultima In Vitro Multiphoton Microscope system (Bruker, Billerica, MA) with a Hamamatsu H7422P-40 GaAsP photomultiplier tube (PMT, 490 nm to 560 nm, Hamamatsu Photonics, Hamamatsu, Japan) and a Hamamatsu R3982 side-on PMT(580-620 nm). Dodt-tube-based transmission detector with Hamamatsu R3982 PMT allowed cell visualization during laser scanning. Scanning signals were sent and received by the NI PCI-6110 analog-to-digital converter card in the system computer (Bruker Nano Fluorescence). Data acquisition was done by Praire View 5.3-5.5 software (Bruker).

## Solutions

### aCSF

Different types of aCSF are adopted by different groups and optimized for different preparations. The procedure described here is modified for experiments focused on the dorsal striatum. The aCSF adopted for these experiments has the following composition (in mM):

- 124 NaCl
- 3 KCl
- 1 NaH₂PO₄
- 2.0 CaCl₂
- 1.0 MgCl₂
- 26 NaHCO₃
- 13.89 glucose

All aCSF solutions are constantly bubbled with carbogen (95% O₂/5% CO₂).

### Internal Solution

For current-clamp recording, we use a potassium gluconate-based internal solution that contains (in mM):

- 115 K-gluconate
- 20 KCl
- 1.5 MgCl₂
- 5 HEPES
- 0.2 EGTA
- 2 Mg-ATP
- 0.5 Na-GTP
- 10 Na-phosphocreatine

pH 7.25, osmolarity 280-290 mOsm/L. This recipe can be modified depending on animal age, cell type, and experimental design.

For Ca²⁺ imaging and visualization of dendritic structures, the internal solution is supplemented with:
- 100 μM Fluo-4 (Thermo Fisher Scientific, F14200)
- 50 μM Alexa Fluor 568 hydrazide (Thermo Fisher Scientific, A10437)

## Before Start

Stocks of internal solution and fluorescent dyes are prepared prior to the experimental day, aliquoted, and stored at -20°C.

On the day of experiment, ex vivo brain slices are prepared from mice in which either direct-pathway SPNs (dSPNs) or indirect-pathway SPNs (iSPNs) are fluorescently labeled (e.g. D1-TdTomato or D2-eGFP mouse lines), following protocols like this one: [https://dx.doi.org/10.17504/protocols.io.dm6gp328jvzp/v2](https://dx.doi.org/10.17504/protocols.io.dm6gp328jvzp/v2).

## Protocol

### Setting up for electrophysiology and two-photon imaging

1. After being prepared, acute brain slices are kept at room temperature in the holding chamber containing aCSF continuously bubbled with carbogen gas (95% O₂/5% CO₂) until use.

2. Thaw an aliquot of internal solution, an aliquot of Fluo-4, and an aliquot of Alexa 568.

3. Add Fluo-4 and Alexa 568 into the internal solution to make a final concentration of 100 μM Fluo-4 and 50 μM Alexa 568.

4. Fill a 1 ml syringe with the internal solution containing Fluo-4 and Alexa 568. Place a syringe filter on the tip of the syringe and attach a MicroFil Pipette Filler.

5. Turn on the 2PLSM workstation, including the laser, modulators and heated stage. Also turn on MultiClamp 700B amplifier, Axon Digidata 1550B digitizer and micromanipulators.

6. Open the software Prairie View and MicroManager.

7. Start running aCSF through the peristaltic pump, into the recording/imaging chamber. The aCSF reservoir is continuously bubbled with carbogen. Adjust the flow rate to at least 2 ml/min. Make sure the chamber outlet removes the solution from the chamber at the same rate, collecting it into a waste solution collector. A vacuum-based outlet is also recommended because this can help prevent overflow, if available.

### Cell type identification and cell selection for patch clamp

8. Transfer one slice from the holding chamber into the recording/imaging chamber. Gently place a slice hold-down on it.

9. Aided by the camera feed through MicroManager, with the 10X objective, adjust the stage position so that the brain region of interest is in the center of the view.

10. Switch to the 60X immersion objective and bring the focus to ~40-60 microns deep in the slice. Find a region with healthy neurons for patching. Switch to 2PLSM imaging.

11. Within Prairie View, adjust the power and image acquisition settings used for cell type identification (it is recommended to start from lower settings and increase laser power and/or gain if needed). To visualize most fluorescent proteins, we use an excitation wavelength of 920 nm.

12. Start imaging in "Live" mode in both the fluorescence channel (red or green, depending on the fluorescent label of SPNs) and the Dodt channel. Have a third image window open to show the merged image of fluorescence and Dodt signal.

13. Pause 2PLSM imaging and go back to the camera view.

14. Choose a healthy neuron of the desired cell type (i.e. having the correct fluorescence labeling) for patching.

### Whole-cell patch clamp

15. Follow the steps described in this protocol to form a whole-cell patch with the internal solution containing Fluo-4 and Alexa 568: [dx.doi.org/10.17504/protocols.io.rm7vzx1w2gx1/v2](https://dx.doi.org/10.17504/protocols.io.rm7vzx1w2gx1/v2).

### Imaging of dendritic calcium transients

16. After whole-cell recording configuration is established, the patched cell is held in current-clamp mode and allowed to equilibrate with dyes for at least 15 min before imaging.

17. Change the laser wavelength to 810 nm, the preferred wavelength for imaging Fluo-4 and Alexa 568.

18. In "Live" imaging mode, visualize the dendritic structure by the red signal of Alexa Fluor 568 (detected by a Hamamatsu R3982 side-on photomultiplier tube/PMT, 580-620 nm). Fluo-4 signals in the green channel were detected by a Hamamatsu H7422P-40 GaAsP PMT (490-560 nm, Hamamatsu Photonics, Japan).

19. Using the Ruler function of Prairie View, find proximal (~40 microns to soma) and distal (~90 microns to soma) locations that are on the same dendrites.

20. Zoom in on the dendrite with the following imaging parameters: zoom = 8, resolution = 128 pixels per line, dwell time = 10 μs/pixel.

21. Stop "Live" imaging.

22. Open the "Linescan" window in Prairie View. Draw a small line along the dendrite (either "Freehand" or "Straight line", 17 pixels) to specify the location of line-scan imaging. Also, specify the length of imaging or number of scans.

23. Open the "Voltage Output" window and choose the voltage output protocol to be used (for example, to evoke back-propagating action potentials, we use three 2 nA current injections, 2 ms each, at 50 Hz).

24. Within the "Linescan" window, "Synchronize" the line scan with both voltage output and voltage recording.

25. Specify the folder and filename for saving the line-scan data.

26. Start line-scan. Repeat three trials for each dendritic location.

27. Both Fluo-4 and Alexa 568 signals are acquired. Only the Fluo-4 signals are responsive to Ca²⁺-elevating stimulations whereas Alexa 568 signals should be stable and provide an internal control (e.g. for physical vibration or laser stability). PMT voltage settings, once optimized, are fixed for all the experiments.

28. Laser power is adjusted so that the Fluo-4 fluorescence at baseline is bright but far from reaching saturation of the signal. In our conditions, signal saturation is experienced above 3000 fluorescence units; the baseline fluorescence of line-scan for the Fluo-4 signal is normally adjusted to average at 600-1000 units (background caused by PMTs are around 60-160 units). This should allow us to easily measure fluorescence increases without saturation.

### Imaging of dendritic spines

29. Close the "Linescan" window and open the "Z-stack" window.

30. Switch to voltage-clamp mode and hold the patched SPN at -80 mV.

31. Turn off green PMT but leave red PMT open.

32. At zoom=1, select a proximal or distal dendritic segment for imaging.

33. Set the imaging parameters for dendritic spine imaging: zoom = 5.2, resolution = 256 x 256, frame average = 4, dwell time = 10 μs/pixel.

34. In "Live" imaging mode, set the positions of the upper and lower limits of the Z-stack in "Z-stack" window.

35. Specify Z-stack step size (0.3 micron or less).

36. Adjust laser power so that even the images from deeper layers of the stack are bright enough. A laser power gradient can be used for optimal imaging.

37. Define the folder and filename for saving the z-stack image.

38. Start the Z-stack imaging and wait until it finishes.

## After the experiment

39. Make sure the PMTs are turned off and properly shuttered before turning on the room light.

40. Discard slices and waste solutions according to institutional protocols.

41. Carefully wash tubing, recording/imaging chamber, holding chamber, Microfil pipette filler, slice holder, and 60X objective.

42. Turn off all equipment.

43. Export data and proceed with data analysis.
