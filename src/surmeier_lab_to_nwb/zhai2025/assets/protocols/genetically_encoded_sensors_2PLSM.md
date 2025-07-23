# Ex vivo imaging of genetically encoded sensors with two-photon laser scanning microscopy (2PLSM)

*Forked from Ex vivo Ca²⁺ 2PLSM measurements with genetically encoded probes*

**DOI:** [dx.doi.org/10.17504/protocols.io.n92ldry68g5b/v1](https://dx.doi.org/10.17504/protocols.io.n92ldry68g5b/v1)

## Authors
- **Enrico Zampese**¹ - Northwestern University, Feinberg School of Medicine
- **Shenyu Zhai**¹ - Northwestern University

## Protocol Information

**DOI:** [dx.doi.org/10.17504/protocols.io.n92ldry68g5b/v1](https://dx.doi.org/10.17504/protocols.io.n92ldry68g5b/v1)

**Protocol Citation:** Enrico Zampese, Shenyu Zhai 2024. Ex vivo imaging of genetically encoded sensors with two-photon laser scanning microscopy (2PLSM). protocols.io https://dx.doi.org/10.17504/protocols.io.n92ldry68g5b/v1

**Manuscript citation:** [https://www.science.org/doi/10.1126/sciadv.abp8701](https://www.science.org/doi/10.1126/sciadv.abp8701)

**License:** This is an open access protocol distributed under the terms of the Creative Commons Attribution License, which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited

**Protocol status:** Working - We use this protocol and it's working

**Created:** November 14, 2024
**Last Modified:** December 20, 2024
**Protocol Integer ID:** 112135

**Keywords:** 2P, ex vivo slices, two-photon, genetically encoded sensors, imaging, ex vivo brain slice, encoded probe, photon laser, scanning microscopy, 2plsm measurement, encoded sensor, imaging, 2plsm, ca2, ex vivo

## Funding Acknowledgements
**Aligning Science Across Parkinson's (ASAP)**
Grant ID: ASAP-020551

## Abstract
This protocol is adapted from "Ex vivo Ca²⁺ 2PLSM measurements with genetically encoded probes" ([dx.doi.org/10.17504/protocols.io.5jyl8pby7g2w/v1](https://dx.doi.org/10.17504/protocols.io.5jyl8pby7g2w/v1)) and modified for imaging genetically encoded sensors (e.g. GRABACh3.0) in ex vivo brain slices with two-photon laser scanning microscopy (2PLSM).

## Guidelines
Animal use procedures were approved by the Northwestern Institutional Animal Care and Use committee.

## Materials

### Equipment and Software
- 2PLSM optical workstation and computer with imaging softwares (see below)
- Artificial cerebrospinal fluid (aCSF) (see below)
- Holding chamber
- LED light (e.g. CoolLED pE-100-470)
- Carbogen (95% O₂, 5% CO₂) tank connected to regulator and bubblers
- Slice hold-down (e.g. Warner Instrument, 64-1418)
- Ex vivo brain slices expressing genetically encoded sensors (e.g.: GRABACh3.0) (see protocols on stereotaxic surgeries for viral injections: [https://dx.doi.org/10.17504/protocols.io.81wgby191vpk/v1](https://dx.doi.org/10.17504/protocols.io.81wgby191vpk/v1) and ex vivo slice preparation: [https://dx.doi.org/10.17504/protocols.io.dm6gp328jvzp/v2](https://dx.doi.org/10.17504/protocols.io.dm6gp328jvzp/v2))
- Peristaltic pump with tubing and connectors, including inlet and outlet to microscope's imaging chamber
- Microscope heating system with probe
- Freshly prepared neuromodulator solution (e.g. 100 mM acetylcholine chloride stock solution for GRABACh3.0) for measurement of Fₘₐₓ
- Tetrodotoxin (TTX) solution (10 mM) for measurement of Fₘᵢₙ
- Waste solution collector
- Concentric bipolar electrode (CBAPD75, FHC)

### 2PLSM Optical Workstation and Computer Software

The laser scanning optical workstation embodies an Ultima dual-excitation-channel scan head (Bruker Nano Fluorescence Microscopy Unit). The foundation of the system is the Olympus BX-51WIF upright microscope with a 60X/0.9NA water-dipping objective lens. The automation of the XY stage motion, lens focus, and manipulator XYZ movement was provided by FM-380 shifting stage, axial focus module for Olympus scopes, and manipulators (Luigs & Neumann). Cell visualization were made possible by a ½" CCD video camera (Hitachi) imaged through a Dodt contrast tube, a 2x magnification changer (Bruker), and microManager software.

A two-photon laser (Chameleon Ultra II, Coherent, Santa Clara, CA) was used to excite the sensor (e.g. GRABACh3.0). The excitation wavelength was selected based on the probe being imaged. Laser power attenuation was achieved with two Pockels' cell electro-optic modulators (models M350-80-02-BK and M350-50-02-BK, Con Optics) controlled by Prairie View 5.3–5.5 software. The two modulators were aligned in series to provide enhanced modulation range for fine control of the excitation dose, to limit the sample maximum power, and to serve as a rapid shutter during line scan or time series acquisitions.

Fluorescence was imaged using an Ultima In Vitro Multiphoton Microscope system (Bruker, Billerica, MA) with a Hamamatsu H7422P-40 GaAsP PMT (490 nm to 560 nm, Hamamatsu Photonics, Hamamatsu, Japan). Dodt-tube-based transmission detector with Hamamatsu R3982 PMT allowed cell visualization during laser scanning. Scanning signals were sent and received by the NI PCI-6110 analog-to-digital converter card in the system computer (Bruker Nano Fluorescence). Imaging data were acquired using Praire View 5.3-5.5 software (Bruker).

### aCSF Solution

Different types of aCSF are adopted by different groups and optimized for different preparations. The procedure described here is modified for experiments focused on the dorsal striatum. The aCSF adopted for these experiments has the following composition (in mM): 124 NaCl, 3 KCl, 1 NaH₂PO₄, 2.0 CaCl₂, 1.0 MgCl₂, 26 NaHCO₃ and 13.89 glucose. All aCSF solutions are constantly bubbled with carbogen (95% O₂/5% CO₂).

## Before Start
An AAV that expresses a genetically encoded sensor (e.g. GRABACh3.0) is injected into the brain region of interest through stereotaxic surgery, following the protocol: [https://dx.doi.org/10.17504/protocols.io.dm6gp328jvzp/v2](https://dx.doi.org/10.17504/protocols.io.dm6gp328jvzp/v2) or a similar protocol. A few weeks after viral injection, ex vivo slices are prepared following protocols like: [https://dx.doi.org/10.17504/protocols.io.dm6gp328jvzp/v2](https://dx.doi.org/10.17504/protocols.io.dm6gp328jvzp/v2).

## Standard Experimental Procedure

1. **Prepare brain slices**: Acute brain slices expressing genetically encoded sensors are prepared and kept at room temperature in a chamber containing aCSF continuously bubbled with carbogen gas (95% O₂/5% CO₂) until the start of the imaging experiment.

2. **Setup equipment**: Turn on the 2PLSM workstation, including the laser, modulators, heated stage and computer softwares.

3. **Start perfusion**: Start running aCSF through the peristaltic pump, into the microscope chamber. The reservoir of aCSF is also continuously bubbled with carbogen. Make sure the chamber outlet removes the solution from the chamber at the same rate, collecting it into a waste solution collector. A vacuum-based outlet is also recommended because this can help prevent overflow, if available.

4. **Temperature equilibration**: Temperature probe should be inserted in the solution in the chamber. As the system is turned on, the temperature in the microscope chamber should reach 32-33°C.

5. **Transfer slice**: Transfer one slice from the holding chamber into the imaging chamber. Gently place a slice hold-down on it and move it to the center of the viewing field.

6. **Initial verification**: With the eyepieces and using the LED as the light source, verify with a low magnification objective the expression of the sensor in the brain region of interest and adjust the stage position so that the sensor-expressing region is in the center of the view.

7. **Position electrode**: Place the tip of the concentric bipolar electrode onto the slice. The location and intensity of the stimulation depend on the experimental design but should remain the same for all the experiments.

8. **Switch to high magnification**: Switch to the 60X immersion objective and focus on the slice surface. Switch to 2PLSM imaging.

9. **Adjust imaging settings**: Within the imaging software, adjust power and image acquisition settings (it is recommended to start from lower settings and increase laser power and/or gain if needed) and start imaging in "Live" mode. The 2P excitation wavelength for GFP-based sensors is 920nm. Other kinds of sensors might require different excitation wavelengths.

10. **Optimize parameters**: Identify a cell/region to image from, optimize imaging settings including zoom, field of view, resolution, dwell time, frame rate. For imaging GRABACh3.0, our preferred settings are: 0.388 µm × 0.388 µm pixels, 8-µs dwell time, and a frame rate of 21.26 fps.

11. **Signal optimization**: Laser power and PMTs gain are adjusted so that the fluorescence at baseline is bright but far from reaching saturation of the signal. In our conditions, signal saturation is experienced above 3000 fluorescence units; the baseline fluorescence for the region of interest (ROI) is normally adjusted to average at ~600 units (background caused by PMTs are around 60-160 units). This should allow us to easily measure fluorescence increases without saturation.

12. **Equilibration period**: It is recommended to wait at least 10-15 mins after placing the slice in the chamber and lowering the 60x objective before starting any experiment. This should give sufficient time to the slice to stabilize and equilibrate properly with the working temperature of the chamber. Not waiting a sufficient time might result in changes in focus/movement and instable fluorescent baseline.

13. **Imaging modes**: 2PLSM imaging with genetically encoded sensors can be performed either as time-lapse experiments or as continuous acquisitions, depending on the timescale of the phenomenon under observation.
    - **13.1**: For slow pharmacological effects, time-lapse acquisitions are preferred. In this case, a series of frames are acquired at regular intervals. Each series of frames is then averaged to represent one time point.
    - **13.2**: For faster effects (e.g. acute stimulation), continuous acquisitions are preferred. In Prairie View, this is performed as a Brightness Over Time (BOT) acquisition.

14. **Calibration (optional)**: If desired, for a semi-quantitative approach, it is recommended to perform an in situ calibration of the dynamic range of each probe for each ROI examined. This is obtained by estimating minimum and maximum fluorescence levels (Fₘᵢₙ and Fₘₐₓ) for each ROI.

## Calibration Protocol

15. **Minimum fluorescence (Fₘᵢₙ)**: To obtain the minimum fluorescence (Fₘᵢₙ), aCSF containing 10 µM TTX is washed in for at least 10 min. BOT acquisition is performed until a stable level is reached.

16. **Maximum fluorescence (Fₘₐₓ)**: To obtain the maximum fluorescence (Fₘₐₓ), aCSF containing a saturating concentration of the neurotransmitter/neuromodulator (e.g. 100 µM acetylcholine chloride for GRABACh3.0) is washed in for at least 10 min. BOT acquisition is performed until a stable maximum level is obtained.

17. **Cleanup**: Slices are discarded and the imaging chamber is cleaned after Fₘᵢₙ and Fₘₐₓ measurements.

## After the Experiment

18. **Disposal**: Discard slices and waste solutions according to institutional protocols.

19. **Equipment cleaning**: Carefully wash tubing, microscope chamber, stimulating electrode, slice holder, and 60X objective.

20. **Equipment shutdown**: Turn off all equipments.

21. **Data analysis**: Export data and proceed with image analysis.
