# Supplementary Figure 3 Conversion Notes

## Overview
Supplementary Figure 3 contains behavioral video recordings of M1R CRISPR mice during levodopa-induced dyskinesia testing. The data includes:
- **Behavioral videos**: Video recordings of contralateral rotations and abnormal involuntary movements
- **Longitudinal assessment**: Multiple sessions across treatment progression
- **M1R CRISPR validation**: Behavioral phenotype of M1R-deleted mice

## Data Structure

```
Supplementary Figure 3/
└── M1R CRISPR videos/
    ├── 01172024-/
    │   ├── 1st session/
    │   ├── 2nd session/
    │   ├── 3rd session/
    │   ├── 4th session/
    │   └── 5th session/
    ├── 02232024-/
    ├── 05032024-/
    ├── 09272023-/
    └── 12082023-/
```

## Experimental Protocols

### Behavioral Video Recording Protocol

**Methodology**: Video recording of mouse behavior during levodopa treatment sessions

**Recording Schedule**:
- **Sessions**: 5 sessions total, every other day
- **Session progression**:
  - Day 1: 6 mg/kg L-DOPA + 12 mg/kg benserazide
  - Day 3: 6 mg/kg L-DOPA + 12 mg/kg benserazide
  - Day 5: 12 mg/kg L-DOPA + 12 mg/kg benserazide
  - Day 7: 12 mg/kg L-DOPA + 12 mg/kg benserazide
  - Day 9: 12 mg/kg L-DOPA + 12 mg/kg benserazide

**Recording Time Points**: 20min, 40min, 60min, 80min post-injection

### M1R CRISPR Experimental Design

**Groups**:
- **M1R CRISPR mice**: AAV-Cas9 + AAV-gRNA-FusionRed → M1R deletion in iSPNs
- **Control mice**: Saline + AAV-gRNA-FusionRed → No M1R deletion

**Behavioral Measures**:
- **Contralateral rotations**: Rotations toward lesioned side
- **AIM scoring**: Abnormal involuntary movements (axial, limb, orolingual)

## Data Organization

### Video File Structure

**Naming Convention**: `[MouseID]_[SessionNumber]_[TimePoint].mov`

**Directory Structure**:
```
01172024-/                    # Date: Jan 17, 2024 - Experimental cohort 1
├── 1st session/              # Day 1: 6 mg/kg L-DOPA + 12 mg/kg benserazide
│   ├── [MouseID]_1_20min.mov # Mouse ID, Session 1, 20min post-injection
│   ├── [MouseID]_1_40min.mov # Mouse ID, Session 1, 40min post-injection
│   ├── [MouseID]_1_60min.mov # Mouse ID, Session 1, 60min post-injection
│   └── [MouseID]_1_80min.mov # Mouse ID, Session 1, 80min post-injection
├── 2nd session/              # Day 3: 6 mg/kg L-DOPA + 12 mg/kg benserazide
│   ├── [MouseID]_2_20min.mov # Mouse ID, Session 2, 20min post-injection
│   ├── [MouseID]_2_40min.mov # Mouse ID, Session 2, 40min post-injection
│   ├── [MouseID]_2_60min.mov # Mouse ID, Session 2, 60min post-injection
│   └── [MouseID]_2_80min.mov # Mouse ID, Session 2, 80min post-injection
├── 3rd session/              # Day 5: 12 mg/kg L-DOPA + 12 mg/kg benserazide
├── 4th session/              # Day 7: 12 mg/kg L-DOPA + 12 mg/kg benserazide
└── 5th session/              # Day 9: 12 mg/kg L-DOPA + 12 mg/kg benserazide

02232024-/                    # Date: Feb 23, 2024 - Experimental cohort 2
├── 1st session/              # Sessions every other day
├── 2nd session/
├── 3rd session/
├── 4th session/
└── 5th session/

05032024-/                    # Date: May 3, 2024 - Experimental cohort 3
├── 1st session/
├── 2nd session/
├── 3rd session/
├── 4th session/
└── 5th session/

09272023-/                    # Date: Sep 27, 2023 - Experimental cohort 4
├── 1st session/
├── 2nd session/
├── 3rd session/
├── 4th session/
└── 5th session/

12082023-/                    # Date: Dec 8, 2023 - Experimental cohort 5
├── 1st session/
├── 2nd session/
├── 3rd session/
├── 4th session/
└── 5th session/
```

### Video File Properties

**Format**: .mov (QuickTime Movie format)
**Content**: Behavioral recordings of individual mice
**Duration**: Typically several minutes per recording
**Purpose**: Document contralateral rotations and dyskinetic movements

## Key Metadata

### Recording Parameters
- **Time points**: 20, 40, 60, 80 minutes post-injection
- **Session frequency**: Every other day (5 sessions total)
- **Drug escalation**: Increasing L-DOPA doses across sessions

### Behavioral Measures
- **Contralateral rotations**: Quantified rotations toward lesioned side
- **AIM categories**: Axial, limb, and orolingual abnormal movements
- **Scoring intervals**: Discrete time windows during each session

### Mouse Information
- **Strain**: M1R CRISPR mice and controls
- **Lesion**: 6-OHDA unilateral lesion model
- **Treatment**: L-DOPA + benserazide (carbidopa analog)

## Experimental Design Notes

### Longitudinal Design
- **Baseline**: Pre-treatment behavioral assessment
- **Progression**: 5 sessions over 9 days
- **Escalation**: Increasing L-DOPA doses to induce dyskinesia
- **Documentation**: Video evidence of behavioral changes

### M1R CRISPR Validation
- **Hypothesis**: M1R deletion reduces dyskinetic behaviors
- **Comparison**: CRISPR vs. control mice behavioral responses
- **Outcome**: Video documentation of reduced AIMs in M1R CRISPR mice

### Video Analysis
- **Manual scoring**: Trained observers score videos for AIMs
- **Rotation counting**: Quantification of contralateral rotations
- **Blinded analysis**: Scorers blinded to experimental group

## Behavioral Assessment Protocol

### L-DOPA Treatment Schedule
1. **Session 1 (Day 1)**: 6 mg/kg L-DOPA + 12 mg/kg benserazide
2. **Session 2 (Day 3)**: 6 mg/kg L-DOPA + 12 mg/kg benserazide
3. **Session 3 (Day 5)**: 12 mg/kg L-DOPA + 12 mg/kg benserazide
4. **Session 4 (Day 7)**: 12 mg/kg L-DOPA + 12 mg/kg benserazide
5. **Session 5 (Day 9)**: 12 mg/kg L-DOPA + 12 mg/kg benserazide

### Video Recording Timeline
- **Pre-injection**: Baseline recording
- **20 min post**: Early response phase
- **40 min post**: Peak response phase
- **60 min post**: Sustained response phase
- **80 min post**: Late response phase

## Related Files
- [Overview](conversion_notes_overview.md)
- [Figure 7 Notes](figure_7_conversion_notes.md) - CDGI knockout behavioral videos
- [Figure 8 Notes](figure_8_conversion_notes.md) - M1R CRISPR AIM scoring data

## Technical Notes

### Video File Handling
- **Storage**: Large file sizes require efficient storage solutions
- **Compression**: Consider video compression for NWB storage
- **Metadata**: Extract recording timestamps and session information
- **Synchronization**: Align video timing with behavioral scoring data

### NWB Conversion Considerations
- **BehavioralTimeSeries**: Store quantified behavioral measures
- **ImageSeries**: Store video data with appropriate compression
- **TimeIntervals**: Define behavioral epochs and scoring windows
- **Subject metadata**: Include mouse ID, genotype, and treatment history
