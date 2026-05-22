"""
Utility functions for extracting conditions from NWB session IDs.

This module provides hardcoded dictionaries for experimental condition mappings
and shared utility functions used in the Zhai 2025 paper conversion scripts.
"""

import argparse
import re
from typing import Optional

# Mapping from folder path names to paper condition names
FOLDER_TO_PAPER_CONDITION = {
    # Figure 1 conditions
    "LID off-state": "LID off-state",
    "LID on-state": "LID on-state",
    "LID on-state with SCH": "LID on-state with SCH",
    # Figure 3 conditions
    "LID on-state with sul": "LID on-state with sul (iSPN)",
    "LID on-state with sul (iSPN)": "LID on-state with sul (iSPN)",
    # Figure 6 conditions
    "control": "control",
    "M1R antagonist": "M1R antagonist",
    # Figure 7 conditions
    "KO off-state": "KO off-state",
    "KO on-state": "KO on-state",
    # Figure 8 conditions
    "M1R CRISPR": "M1R CRISPR",
    "interleaved control": "interleaved control",
}

# Main condition mapping table
format_condition = {
    # Figure 1 conditions
    "LID off-state": {
        "CamelCase": "LIDOffState",
        "underscore": "lid_off_state",
        "human_readable": "LID off-state",
        "description_and_purpose": "Levodopa-induced dyskinesia animals in off-state (24h after last L-DOPA injection)",
    },
    "LID on-state": {
        "CamelCase": "LIDOnState",
        "underscore": "lid_on_state",
        "human_readable": "LID on-state",
        "description_and_purpose": "Levodopa-induced dyskinesia animals in on-state (1h after L-DOPA injection)",
    },
    "LID on-state with SCH": {
        "CamelCase": "LIDOnStateWithSchD1Antagonist",
        "underscore": "lid_on_state_with_sch",
        "human_readable": "LID on-state with SCH D1 antagonist",
        "description_and_purpose": "LID on-state with SCH-23390 D1 receptor antagonist to block direct pathway",
    },
    # Figure 3 conditions
    "LID on-state with sul (iSPN)": {
        "CamelCase": "LIDOnStateWithSulpiride",
        "underscore": "lid_on_state_with_sul",
        "human_readable": "LID on-state with sulpiride",
        "description_and_purpose": "LID on-state with sulpiride D2 receptor antagonist to block indirect pathway",
    },
    # Figure 5 conditions
    "6-OHDA": {
        "CamelCase": "SixOHDA",
        "underscore": "6_ohda",
        "human_readable": "6-OHDA",
        "description_and_purpose": "6-hydroxydopamine lesioned animals (Parkinson's disease model)",
    },
    "off-state": {
        "CamelCase": "OffState",
        "underscore": "off_state",
        "human_readable": "off-state",
        "description_and_purpose": "Animals in off-state (without L-DOPA treatment)",
    },
    # Figure 6 conditions
    "control": {
        "CamelCase": "Control",
        "underscore": "control",
        "human_readable": "control",
        "description_and_purpose": "Control condition for M1R antagonist experiments",
    },
    "M1R antagonist": {
        "CamelCase": "M1RAntagonist",
        "underscore": "m1r_antagonist",
        "human_readable": "M1R antagonist",
        "description_and_purpose": "Muscarinic M1 receptor antagonist treatment",
    },
    # Figure 7 conditions
    "KO off-state": {
        "CamelCase": "KOOffState",
        "underscore": "ko_off_state",
        "human_readable": "KO off-state",
        "description_and_purpose": "ChI-M1R knockout animals in off-state",
    },
    "KO on-state": {
        "CamelCase": "KOOnState",
        "underscore": "ko_on_state",
        "human_readable": "KO on-state",
        "description_and_purpose": "ChI-M1R knockout animals in on-state",
    },
    # Figure 7 oxoM treatment conditions
    "WT oxoM treatment": {
        "CamelCase": "WTOxoMTreatment",
        "underscore": "wt_oxom_treatment",
        "human_readable": "WT oxoM treatment",
        "description_and_purpose": "Wildtype animals with oxotremorine-M muscarinic agonist treatment",
    },
    "CDGI KO oxoM treatment": {
        "CamelCase": "CDGIKOOxoMTreatment",
        "underscore": "cdgi_ko_oxom_treatment",
        "human_readable": "CDGI KO oxoM treatment",
        "description_and_purpose": "CalDAG-GEFI knockout animals with oxotremorine-M muscarinic agonist treatment",
    },
    # Figure 8 conditions
    "M1R CRISPR": {
        "CamelCase": "M1RCRISPR",
        "underscore": "m1r_crispr",
        "human_readable": "M1R CRISPR",
        "description_and_purpose": "M1R knockout using CRISPR-Cas9 system",
    },
    "interleaved control": {
        "CamelCase": "InterleavedControl",
        "underscore": "interleaved_control",
        "human_readable": "interleaved control",
        "description_and_purpose": "Control condition interleaved with M1R CRISPR experiments",
    },
    # Spine density conditions
    "control dSPN": {
        "CamelCase": "ControlDSPN",
        "underscore": "control_dspn",
        "human_readable": "control dSPN",
        "description_and_purpose": "Control direct pathway spiny projection neurons",
    },
    "LID off-state dSPN": {
        "CamelCase": "LIDOffStateDSPN",
        "underscore": "lid_off_state_dspn",
        "human_readable": "LID off-state dSPN",
        "description_and_purpose": "LID off-state direct pathway spiny projection neurons",
    },
    "LID on-state dSPN": {
        "CamelCase": "LIDOnStateDSPN",
        "underscore": "lid_on_state_dspn",
        "human_readable": "LID on-state dSPN",
        "description_and_purpose": "LID on-state direct pathway spiny projection neurons",
    },
    "PD dSPN": {
        "CamelCase": "PDDSN",
        "underscore": "pd_dspn",
        "human_readable": "PD dSPN",
        "description_and_purpose": "Parkinson's disease direct pathway spiny projection neurons",
    },
    "control iSPN": {
        "CamelCase": "ControlISPN",
        "underscore": "control_ispn",
        "human_readable": "control iSPN",
        "description_and_purpose": "Control indirect pathway spiny projection neurons",
    },
    "LID off-state iSPN": {
        "CamelCase": "LIDOffStateISPN",
        "underscore": "lid_off_state_ispn",
        "human_readable": "LID off-state iSPN",
        "description_and_purpose": "LID off-state indirect pathway spiny projection neurons",
    },
    "LID on-state iSPN": {
        "CamelCase": "LIDOnStateISPN",
        "underscore": "lid_on_state_ispn",
        "human_readable": "LID on-state iSPN",
        "description_and_purpose": "LID on-state indirect pathway spiny projection neurons",
    },
    "PD iSPN": {
        "CamelCase": "PDISPN",
        "underscore": "pd_ispn",
        "human_readable": "PD iSPN",
        "description_and_purpose": "Parkinson's disease indirect pathway spiny projection neurons",
    },
    "KO off": {
        "CamelCase": "KOOff",
        "underscore": "ko_off",
        "human_readable": "KO off",
        "description_and_purpose": "ChI-M1R knockout animals in off-state for spine density analysis",
    },
    "KO on": {
        "CamelCase": "KOOn",
        "underscore": "ko_on",
        "human_readable": "KO on",
        "description_and_purpose": "ChI-M1R knockout animals in on-state for spine density analysis",
    },
}

# Create reverse mapping for quick lookup
CAMELCASE_TO_PAPER_CONDITION = {formats["CamelCase"]: condition for condition, formats in format_condition.items()}

# Canonical session ID mappings for 6-token format
FIG_MAPPING = {
    "Figure1": "F1",  # dSPN excitability
    "Figure2": "F2",  # dSPN spines & oEPSCs
    "Figure3": "F3",  # iSPN excitability
    "Figure4": "F4",  # iSPN spines & oEPSCs
    "Figure5": "F5",  # GRAB-ACh photometry + behaviour
    "Figure6": "F6",  # M1R antagonist in iSPNs
    "Figure7": "F7",  # CDGI knock-out
    "Figure8": "F8",  # iSPN-specific M1R KO
    "SupplementaryFigure3": "FSup3",  # High-res confocal spines
    "SupplementaryFigure": "FSup",  # Alternative naming
}

COMPARTMENT_MAPPING = {
    "SomaticExcitability": "som",  # Whole-cell recording at soma
    "DendriticExcitability": "dend",  # Dendritic patch or imaging
    "AcetylcholineGRAB": "stri",  # In-vivo striatal signal
    "BehavioralAIM": "behav",  # Whole-animal behaviour
    "BehavioralVideos": "behav",  # Video recording
    "SpineDensity": "dend",  # Dendritic imaging
    "ConfocalSpineDensity": "dend",  # Confocal dendritic imaging
    "OpticalStimuli": "dend",  # Dendritic stimulation
}

MEASUREMENT_MAPPING = {
    "SomaticExcitability": None,  # Default for intrinsic excitability
    "DendriticExcitability": None,  # Default for intrinsic excitability
    "SpineDensity": "spine",  # 2-photon spine density
    "ConfocalSpineDensityNikon": "confSpine",  # Confocal spine density
    "ConfocalSpineDensityOlympus": "confSpine",  # Confocal spine density
    "OpticalStimuli": "oEPSC",  # Sr²⁺ oEPSC amplitude
    "AcetylcholineGRAB": "AChFP",  # GRAB-ACh photometry
    "BehavioralAIM": "AIMs",  # AIM score
    "BehavioralVideos": "video",  # Raw video
}

SPN_MAPPING = {
    "dSPN": "dspn",  # Direct pathway
    "iSPN": "ispn",  # Indirect pathway
    "pan": "whole",  # Non cell-specific
    None: "NA",  # NA
}

STATE_MAPPING = {
    "LID on-state": "ON",
    "LID off-state": "OFF",
    "LID on-state with SCH": "ON",  # Still on-state, just with D1 antagonist
    "LID on-state with sul (iSPN)": "ON",  # Still on-state, just with D2 antagonist
    "KO on-state": "ON",
    "KO off-state": "OFF",
    "control": "CTRL",
    "6-OHDA": "PD",
    "PD": "PD",
    "M1R CRISPR": "OFF",
    "interleaved control": "OFF",
    "control dSPN": "CTRL",
    "control iSPN": "CTRL",
    "LID on-state dSPN": "ON",
    "LID off-state dSPN": "OFF",
    "LID on-state iSPN": "ON",
    "LID off-state iSPN": "OFF",
    "PD dSPN": "PD",
    "PD iSPN": "PD",
    "KO on": "ON",
    "KO off": "OFF",
    "WT oxoM treatment": "ON",  # oxoM induces on-like state
    "CDGI KO oxoM treatment": "ON",
    "off-state": "OFF",
}

PHARM_MAPPING = {
    "none": "none",
    "LID on-state with SCH": "D1RA",
    "LID on-state with sul (iSPN)": "D2RA",
    "M1R antagonist": "M1RA",
    "control": "none",
    "WT oxoM treatment": "none",  # oxoM is part of the treatment, not pharmacology
    "CDGI KO oxoM treatment": "none",
    "interleaved control": "none",
    "M1R CRISPR": "none",
}

GENO_MAPPING = {
    "WT": "WT",
    "KO": "CDGIKO",
    "CDGI-KO": "CDGIKO",
    "CDGI KO": "CDGIKO",
    "M1R CRISPR": "iSPN-M1RKO",
    "iSPN-M1R-KO": "iSPN-M1RKO",
    "control": "WT",
    "interleaved control": "WT",
}

# Token descriptions for revised session ID schema
TOKEN_DESCRIPTIONS = {
    # Figure descriptions
    "fig": {
        "F1": "Figure 1 - dSPN excitability",
        "F2": "Figure 2 - dSPN spines & oEPSCs",
        "F3": "Figure 3 - iSPN excitability",
        "F4": "Figure 4 - iSPN spines & oEPSCs",
        "F5": "Figure 5 - GRAB-ACh photometry + behaviour",
        "F6": "Figure 6 - M1R antagonist in iSPNs",
        "F7": "Figure 7 - CDGI knock-out",
        "F8": "Figure 8 - iSPN-specific M1R KO",
        "Fconf": "Supplementary confocal spines",
    },
    # Measurement + Compartment descriptions (CamelCase merged tokens)
    "meas_comp": {
        "SomExc": "somatic excitability",
        "DendExc": "dendritic excitability",
        "DendSpine": "dendritic spine density",
        "DendConfSpine": "dendritic confocal spine density",
        "DendOEPSC": "dendritic optically-evoked EPSC",
        "StriAChFP": "striatal acetylcholine fluorescence photometry",
        "BehavAIMs": "behavioral abnormal involuntary movements",
        "BehavRot": "behavioral rotation",
        "BehavVideo": "behavioral video recording",
    },
    # Cell type descriptions
    "cell_type": {
        "dSPN": "direct pathway spiny projection neuron",
        "iSPN": "indirect pathway spiny projection neuron",
        "WholeStriatum": "whole striatum bulk signal",
        "NonCell": "non-cellular (behavioral)",
    },
    # State descriptions (CamelCase)
    "state": {
        "OnState": "dyskinetic on-state",
        "OffState": "parkinsonian off-state",
        "LesionedControl": "lesioned control",
        "UnlesionedControl": "unlesioned control",
    },
    # Pharmacology descriptions (CamelCase with drug tags)
    "pharm": {
        "none": "no pharmacological intervention",
        "D1RaSch": "D1 receptor antagonist (SCH-23390)",
        "D2RaSul": "D2 receptor antagonist (sulpiride)",
        "M1RaTri": "M1 muscarinic receptor antagonist (trihexyphenidyl)",
    },
    # Genotype descriptions (CamelCase, no separators)
    "geno": {
        "WT": "wild-type",
        "CDGIKO": "CalDAG-GEFI knockout",
        "iSPNM1RKO": "iSPN-specific M1 receptor knockout",
    },
}

# Legacy genotype prefix mapping for electrode descriptions (backward compatibility)
GENOTYPE_DESCRIPTION_MAPPING = {
    "WT": "Wild Type",  # No prefix for wild-type background
    "CDGIKO": "CDGI knockout ",  # Global CalDAG-GEFI knock-out – removes M1R effector
    "iSPN-M1RKO": "",  # iSPN-restricted M1R CRISPR knock-out – already described in condition
    "M1RCRISPR": "",  # Alternative naming for M1R CRISPR
}

# Pharmacology additions mapping - keys match session ID pharm tokens
PHARMACOLOGY_ADDITIONS = {
    "D1RaSch": "SCH-23390: D1 receptor antagonist (bath application) to block direct pathway during LID on-state.",
    "D2RaSul": "Sulpiride: D2 receptor antagonist (bath application) to investigate the role of D2 receptor signaling in LID-induced excitability changes in iSPNs.",
    "M1RaTri": "M1 muscarinic receptor antagonist: Trihexyphenidyl hydrochloride (THP, 3 mg/kg i.p.) administered to assess M1R contribution to LID-induced changes.",
}


def extract_condition_from_session_id(session_id: str) -> Optional[str]:
    """
    Extract experimental condition from NWB session ID.

    Supports the mouse-day session-id schema (current) plus two legacy schemas:
    - Mouse-day (current): "dSPN++OnState++D1RaSch++WT++20190411"
    - Per-session revised: "F1++SomExc++dSPN++OnState++D1RaSch++WT++20170202153933"
    - Per-session legacy:  "Figure1++SomaticExcitability++LIDOffState++20170202153933"

    Parameters
    ----------
    session_id : str
        Session ID from NWB file

    Returns
    -------
    Optional[str]
        Paper condition name (e.g., "LID off-state"), or None if not found
    """
    condition_mapping = {
        ("OffState", "none"): "LID off-state",
        ("OnState", "none"): "LID on-state",
        ("OnState", "D1RaSch"): "LID on-state with SCH",
        ("OnState", "D2RaSul"): "LID on-state with sul (iSPN)",
        ("LesionedControl", "none"): "control",
        ("LesionedControl", "M1RaTri"): "M1R antagonist",
    }

    # Current schema: cellType++state++pharm++geno++YYYYMMDD (5 tokens)
    mouse_day_pattern = r"^(\w+)\+\+(\w+)\+\+(\w+)\+\+(\w+)\+\+(\d{8})$"
    match = re.match(mouse_day_pattern, session_id)
    if match:
        state_token = match.group(2)
        pharm_token = match.group(3)
        if (state_token, pharm_token) in condition_mapping:
            return condition_mapping[(state_token, pharm_token)]

    # Legacy per-session revised schema: fig++measComp++cellType++state++pharm++geno++timestamp
    revised_pattern = r"F\d\+\+\w+\+\+\w+\+\+(\w+)\+\+(\w+)\+\+\w+\+\+\d{14}"
    match = re.search(revised_pattern, session_id)
    if match:
        condition_key = (match.group(1), match.group(2))
        if condition_key in condition_mapping:
            return condition_mapping[condition_key]

    # Legacy per-session schema: Figure[N]++ExperimentType++Condition++Timestamp
    pattern1 = r"Figure\d\+\+\w+\+\+(\w+)\+\+\d{14}"
    match = re.search(pattern1, session_id)
    if match:
        return CAMELCASE_TO_PAPER_CONDITION.get(match.group(1))

    # Most flexible fallback: any known camel-case condition substring
    for camel_case_condition in CAMELCASE_TO_PAPER_CONDITION.keys():
        if camel_case_condition in session_id:
            return CAMELCASE_TO_PAPER_CONDITION[camel_case_condition]

    return None


def generate_canonical_session_id(
    cell_type: str,
    state: str,
    pharm: str,
    geno: str,
    date: str,
) -> str:
    """
    Generate a canonical session ID at the mouse-day level.

    The dataset uses per-mouse-day NWB files (one file collects all cells patched
    from one mouse on one experimental day across all modalities), so the
    session_id reflects mouse-day attributes only. Per-cell, per-modality, and
    per-figure information lives inside the file via the icephys hierarchical
    tables and the electrode column on IntracellularRecordingsTable.

    Schema: cellType++state++pharm++geno++date

    Parameters
    ----------
    cell_type : str
        Cell type: dSPN, iSPN, pan
    state : str
        Levodopa/disease state: OnState, OffState, LesionedControl, UnlesionedControl
    pharm : str
        Pharmacology (CamelCase): none, D1RaSch, D2RaSul, M1RaTri
        For paired-recording cells (baseline + post-drug on the same cell), use
        the most distinguishing token (e.g., D1RaSch for ON+SCH days). The
        ExperimentalConditions table inside the file carries both rows.
    geno : str
        Genotype (CamelCase): WT, CDGIKO, iSPNM1RKO
    date : str
        Recording day in YYYYMMDD format (matches the subject_id date).

    Returns
    -------
    str
        Canonical session ID with format: cellType++state++pharm++geno++date

    Examples
    --------
    >>> generate_canonical_session_id("dSPN", "OnState", "D1RaSch", "WT", "20190411")
    'dSPN++OnState++D1RaSch++WT++20190411'

    >>> generate_canonical_session_id("iSPN", "OffState", "none", "WT", "20170202")
    'iSPN++OffState++none++WT++20170202'
    """
    return f"{cell_type}++{state}++{pharm}++{geno}++{date}"


def str_to_bool(v):
    """
    Convert string or boolean input to boolean value for argparse.

    This function is used consistently across all conversion scripts for handling
    boolean command line arguments like --stub-test.

    Parameters
    ----------
    v : str or bool
        Input value to convert to boolean

    Returns
    -------
    bool
        Converted boolean value

    Raises
    ------
    argparse.ArgumentTypeError
        If input cannot be converted to boolean
    """
    if isinstance(v, bool):
        return v
    if v.lower() in ("yes", "true", "t", "y", "1"):
        return True
    elif v.lower() in ("no", "false", "f", "n", "0"):
        return False
    else:
        raise argparse.ArgumentTypeError("Boolean value expected.")
