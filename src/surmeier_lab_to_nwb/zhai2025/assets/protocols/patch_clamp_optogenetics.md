# Ex vivo mouse brain patch clamp recordings combined with optogenetic stimulation

*Forked from Ex vivo mouse brain patch clamp recordings combined with uncaging and optogenetic stimulation*

**DOI:** [dx.doi.org/10.17504/protocols.io.5qpvokxp9l4o/v1](https://dx.doi.org/10.17504/protocols.io.5qpvokxp9l4o/v1)

## Authors
- **Enrico Zampese**¹ - Northwestern University
- **DeNard V Simmons**¹ - Northwestern University

## Protocol Information

**DOI:** [dx.doi.org/10.17504/protocols.io.5qpvokxp9l4o/v1](https://dx.doi.org/10.17504/protocols.io.5qpvokxp9l4o/v1)

**Protocol Citation:** Enrico Zampese, DeNard V Simmons 2024. Ex vivo mouse brain patch clamp recordings combined with optogenetic stimulation. protocols.io https://dx.doi.org/10.17504/protocols.io.5qpvokxp9l4o/v1

**License:** This is an open access protocol distributed under the terms of the Creative Commons Attribution License, which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited

**Protocol status:** Working - We use this protocol and it's working

**Created:** March 20, 2024
**Last Modified:** March 22, 2024
**Protocol Integer ID:** 97016

**Keywords:** ASAPCRN, electrophysiology, optogenetics, uncaging, optogenetic stimulation, optogenetic stimulation in this protocol, mouse brain patch clamp recording, brain

## Funding Acknowledgements
- **NIH** - Grant ID: NS34696
- **CHDI Foundation** - Grant ID: n/a
- **JPB Foundation** - Grant ID: n/a
- **ASAP** - Grant ID: 020551

## Disclaimer
The protocols.io team notes that research involving animals and humans must be conducted according to internationally-accepted standards and should always have prior approval from an Institutional Ethics Committee or Board.

## Abstract
In this protocol we detail the steps to perform ex-vivo brain slice electrophysiology and optogenetic stimulation.

## Guidelines
Please adhere to institutional guidelines and regulations.

## Materials

### Solutions

#### Internal Solution
These are prepared prior to the experiment day and aliquoted in 1.5 ml tubes, and stored at -20°C until the day of experiment.

**Internal solution composition (in mM):**
- 120 potassium-D-gluconate
- 13 KCl
- 10 HEPES
- 0.05 EGTA
- 4 ATP-Mg²⁺
- 0.5 GTP-Na
- 10 phosphocreatine-di (tris)

pH was adjusted to 7.25 with KOH and osmolarity to 275-280 mOsm.

#### Recording Solution
**(in mM):**
- 125 NaCl
- 3 KCl
- 1 MgCl₂
- 2 CaCl₂
- 25 NaHCO₃
- 1.25 NaH₂PO₄
- 10 glucose

Saturated with 95% O₂-5% CO₂; pH 7.4; 300 mOsm/l

#### Low Ca²⁺ Optogenetics Recording Solution
**(in mM):**
- 124 NaCl
- 4.5 KCl
- 1.8 MgCl₂
- 1.2 CaCl₂
- 25 NaHCO₃
- 1.25 NaH₂PO₄
- 10 glucose

Saturated with 95% O₂-5% CO₂; pH 7.4; 300 mOsm/l

### Hardware and Miscellaneous
- Targeted focal spot blue LED (473 nm pE-100, CoolLED)
- Blood-gas mixture (95% O₂, 5% CO₂) tank connected to bubblers
- Slice holder
- Brain slices expressing genetically encoded opsins in holding chamber with aCSF
- Peristaltic pump or gravity flow perfusion with tubing and connectors, including inlet and outlet to microscope's imaging chamber
- Stage heating system with probe
- Waste solution collector
- 10% Ethanol in water (wash solution)
- Image analysis software

## Safety Warnings
Please follow institutional safety guidelines and chemical safety datasheets.

## Before Start
Please note that optogenetic experiments require the expression of a genetically encoded opsin in neurons of interest. Mice genetically engineered to express opsin in specific neuronal populations are available. Alternatively, opsin expression can be obtained via injections of appropriate viral vectors.

Please refer to existing protocols for brain slices preparation and/or viral injections.

## Procedure

### Prepare Patch Pipettes

1. **Setup puller**: Turn on the Sutter P-1000 puller and enter the desired pull protocol.

2. **Pull pipettes**: Insert a thick-walled borosilicate glass capillary and press pull.

3. **Check resistance**: Pipette resistance must be of 3 to 5 megaohms.

### Setting up Patch Rig and Environment

4. **Turn on equipment**: Turn on the MultiClamp 700B Amplifier, Axon Digidata 1550B digitizer, micromanipulator, computer tower and the associated software. **Note:** amplifier and digitizer must be turned on prior to opening software.

5. **Start gas flow**: Turn on O₂/CO₂ tank and bubble aCSF solution.

6. **Prepare internal solution**: Take an aliquot of internal solution from the -20°C fridge.

7. **Fill syringe**: Fill syringe with internal solution, place a filter on the end of the syringe, and place a MicroFil Pipette Filler on the end of the filter.

8. **Start perfusion**: Start perfusing microscope chamber with aCSF (either with a gravity system or a peristaltic pump; a gravity system is normally preferred because of the lower electrical noise).

9. **Adjust flow rate**: Adjust the flow rate to the desired value (recommended at least 2 mL/min).

10. **Temperature control**: Depending on the experimental protocol, turn on stage heater and set to temperature.

### Examine Slices and Patching Cells

11. **Secure slice**: Secure down slice with a harp (slice anchor).

12. **Locate region**: Locate and focus the desired brain region under the 4x objective.

13. **Switch objective**: Change the microscope lens to the 60x objective.

14. **Focus on cells**: Focus on healthy neurons in slices for patching.

15. **Fill pipette**: For whole cell configuration, fill a glass micropipette one-third full of internal solution. For cell attached configuration fill it with internal solution or aCSF. Ensure there is no residual internal solution on exterior of glass micropipette, as this may introduce salts into the micromanipulator and add additional noise to recordings. Remove any air bubbles by gently flicking the glass micropipette.

16. **Mount pipette**: Gently place the glass micropipette onto the electrode holder.

17. **Position electrode**: Position the electrode using a micromanipulator.

18. **Approach slice**: Under the 60x objective, bring the tip of the glass pipette above the slice.

19. **Apply pressure**: Apply a positive pressure and maintain it. Zero the voltage offset.

20. **Approach cell**: Approach the cell diagonally. The positive pressure should create a small dimple on the cell.

21. **Form seal**: Once a dimple is formed, release the positive pressure, and apply a small amount of negative pressure. The resistance should begin to increase rapidly.

22. **Voltage clamp**: As the resistance increases, clamp the cell at your resting potential of interest (typically -60 mV).

23. **Break in**: After a giga-ohm seal is formed, apply a few quick pulses of negative pressure to break into the cell to record in whole cell configuration.

24. **Start recording**: Start recording.

### Recordings and Optogenetic Stimulation

25. **Equilibration (15 min)**: If in whole cell configuration, allow time for the internal solution to fill the cell.

26. **Acquisition**: Electrophysiological recordings and stimulation are acquired and synchronized via the pClamp software.

27. **Optogenetic stimulation**: Single-photon optogenetic stimulation is obtained with a 473 nm LED, whose pulse number, intensity, and duration are controlled via pClamp.

## After the Experiments

28. **Remove pipette**: Gently remove and discard the recording pipette.

29. **Remove slice**: Remove slice from the recording chamber. If desired, the slice can be gently transferred into a plate and bathed with a fixative solution for further processing.
    - **29.1**: **Note:** use appropriate PPE when working with fixatives, including operating under a chemical hood.

30. **Disposal**: Discard waste according to institutional protocols.

31. **Clean equipment**: Rinse perfusion tubing microscope chamber with water and/or a 10% ethanol in water.

32. **Save data**: Save and export files.
