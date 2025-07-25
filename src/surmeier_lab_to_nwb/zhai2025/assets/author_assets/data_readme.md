# GENERAL INFORMATION
This readme file was generated on [2025-01-23 and 2025-02-04] by [DAVID WOKOSIN]

Title of Dataset Folder: LID paper Zhai et al. 2025:
Description of Dataset: Components supporting raw data and context for research article Zhai et al. 2025
Principal Investigator: D. James Surmeier, j-surmeier@northwestern.edu, ORCID: 0000-0002-6376-5225
Date of Data Collection: 2015-11 to 2024-06

# FILE OVERVIEW
1. Word File, [Zhai et al. SA_manuscript_finalv2] most recent version of the submitted manuscript submitted to  Science Advances (adv8224) on 8 JAN 2025
2. pdf File, [adv8224_SupplementalMaterial_v1] associated Supplemental Materials and Figures
3. Excel File, [Key resources table_Zhai_v1] key resources table for relevant software (incl. Notebooks), protocols, antibodies, viruses, animals, chemicals, and hardware
4. Excel File, [Data Connections] master sheets for raw data connections to panel figures, resources table, and some experiment metadata
5. FOLDER, [Tabular dataset Zhai et al. 2025], includes readme file inside folder; spreadsheets for each Figure with separate pages for each Figure panel’s data
6. FOLDER, [Raw data for Figs], folders arranged by Figure panels

# SUMMARY INFORMATION
Submission spans 8 Figures and 5 Supplemental Figures.
Primary experiment threads are somatic excitability (AP number & rheobase, electrophysiology), dendritic excitability (calcium imaging distal/proximal ratio along individual dendrites, electrophysiology and 2P imaging), dendritic spine density (morphology via 2P and confocal imaging), somatic EPSC minis (optical stimulation and electrophysiology), LID AIM behavior scoring.

All non-behavior experiments were performed on a single Bruker Ultima 2P system (S/N 4503) running Prairie View 5.x software which uses three different A-to-D cards to synchronize (shared clock and DMM memory access) the recorded images with protocol defined experimental input and output voltages from the optical workstation system.  The system experimental configurations for dendritic excitability, somatic EPSC minis, and the ACh BOT imaging from Fig 5 all create .csv files, during acquisition, as the output (raw) data from the images saved within the folder.  All instrument system setting metadata is included in the Prairie View raw data folders.  The confocal spine density images were acquired on a Core Facility Nikon A1R confocal system.  The AIM behavior was all scored manually, in real time, by Dr. Zhai with the videos providing potential refinement for later review.

Brain slice regions pursued are Spiny Projection Neurons: direct, D1, Figs 1 & 2, td-tomato; indirect, D2, Figs 3 & 4, 6 � 8, eGFP.

Main animal manipulations are 6-OHDA (PD) lesioning, LevaDopa treatments with different delays (LID OFF, ON).
