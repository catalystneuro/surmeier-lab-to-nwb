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


def extract_condition_from_session_id(session_id: str) -> Optional[str]:
    """
    Extract experimental condition from NWB session ID.

    Parameters
    ----------
    session_id : str
        Session ID from NWB file (e.g., "Figure1++SomaticExcitability++LIDOffState++20170202153933")

    Returnsm
    -------
    Optional[str]
        Paper condition name (e.g., "LID off-state"), or None if not found
    """
    # Pattern 1: Figure[N]++ExperimentType++Condition++Timestamp (somatic/dendritic excitability)
    pattern1 = r"Figure\d\+\+\w+\+\+(\w+)\+\+\d{14}"
    match = re.search(pattern1, session_id)

    if match:
        camel_case_condition = match.group(1)
        return CAMELCASE_TO_PAPER_CONDITION.get(camel_case_condition)

    # Pattern 2: More flexible pattern to catch variations
    # Look for any known camel case condition in the session ID
    for camel_case_condition in CAMELCASE_TO_PAPER_CONDITION.keys():
        if camel_case_condition in session_id:
            return CAMELCASE_TO_PAPER_CONDITION[camel_case_condition]

    return None


def generate_canonical_session_id(
    fig: str,
    compartment: str,
    measurement: str,
    spn_type: str,
    state: str,
    pharmacology: str,
    genotype: str,
    timestamp: str,
) -> str:
    """
    Generate a canonical session ID following the 6-token format plus timestamp.

    This function creates a standardized session ID that encodes all experimental
    parameters in a machine-parsable format. Each token represents a specific
    experimental dimension, allowing for systematic organization and retrieval
    of data across the entire dataset.

    Parameters
    ----------
    fig : str
        Figure code from allowed values:
        - "F1": Fig 1 – dSPN excitability
        - "F2": Fig 2 – dSPN spines & oEPSCs
        - "F3": Fig 3 – iSPN excitability
        - "F4": Fig 4 – iSPN spines & oEPSCs
        - "F5": Fig 5 – GRAB-ACh photometry + behaviour
        - "F6": Fig 6 – M1R antagonist in iSPNs
        - "F7": Fig 7 – CDGI knock-out
        - "F8": Fig 8 – iSPN-specific M1R KO
        - "Fconf": Suppl. high-res confocal spines

    compartment : str
        Biological compartment from allowed values:
        - "som": Whole-cell recording at soma – intrinsic excitability
        - "dend": Dendritic patch or imaging in slice
        - "stri": In-vivo bulk striatal signal (fiber-photometry)
        - "behav": Whole-animal behaviour / video recording

    measurement : str
        Measurement type from allowed values:
        - "None": Default for intrinsic excitability
        - "spine": Spine density from two-photon stacks
        - "confSpine": Spine density from fixed-tissue confocal stacks
        - "oEPSC": Sr²⁺-evoked asynchronous EPSCs – synaptic strength
        - "AChFP": ΔF/F trace from GRAB-ACh fiber-photometry
        - "AIMs": Abnormal Involuntary Movement score – dyskinesia severity
        - "Rot": Ipsilateral rotation count – motor asymmetry
        - "video": Raw behaviour video (unscored)

    spn_type : str
        SPN pathway type from allowed values:
        - "dspn": Direct-pathway SPN (D1-positive)
        - "ispn": Indirect-pathway SPN (D2-positive)
        - "pan": Non cell-specific / bulk signal

    state : str
        Levodopa/disease state from allowed values:
        - "ON": 30 min post-levodopa (hyper-dopaminergic on-state)
        - "OFF": 24–48 h post-levodopa (hypo-dopaminergic off-state)
        - "CTRL": Unlesioned control – never received levodopa
        - "PD": Parkinsonian baseline (6-OHDA lesion, pre-levodopa)

    pharmacology : str
        Acute pharmacology from allowed values:
        - "none": No acute drug relevant to the assay
        - "D1RA": Bath SCH-23390 – isolates D1-signalling
        - "D2RA": Bath sulpiride – isolates D2-signalling
        - "M1RA": Systemic M1-receptor antagonist during off-state

    genotype : str
        Genetic background from allowed values:
        - "WT": Wild-type background
        - "CDGIKO": Global CalDAG-GEFI knock-out – removes M1R effector
        - "iSPN-M1RKO": iSPN-restricted M1R CRISPR knock-out – tests cell-autonomous M1 role

    timestamp : str
        Session timestamp (typically YYYYMMDDHHMMSS format)

    Returns
    -------
    str
        Canonical session ID with format: fig++compartment++measurement++spn_type++state++pharmacology++genotype++timestamp

    Examples
    --------
    >>> generate_canonical_session_id("F1", "som", "None", "dspn", "OFF", "none", "WT", "20240115143022")
    'F1++som++None++dspn++OFF++none++WT++20240115143022'

    >>> generate_canonical_session_id("F7", "behav", "AIMs", "pan", "OFF", "none", "CDGIKO", "20240115143022")
    'F7++behav++AIMs++pan++OFF++none++CDGIKO++20240115143022'

    Notes
    -----
    - Use "None" for any token that is not applicable (e.g., measurement="None" for intrinsic excitability)
    - All tokens are required and must be from the allowed values listed above
    - The session ID uses "++" as delimiter with no underscores or dashes
    - This format enables systematic queries across experimental dimensions
    """
    # Validate inputs match expected values (optional but helpful for debugging)
    valid_figs = {"F1", "F2", "F3", "F4", "F5", "F6", "F7", "F8", "Fconf"}
    valid_compartments = {"som", "dend", "stri", "behav"}
    valid_measurements = {"None", "spine", "confSpine", "oEPSC", "AChFP", "AIMs", "Rot", "video"}
    valid_spn_types = {"dspn", "ispn", "pan"}
    valid_states = {"ON", "OFF", "CTRL", "PD"}
    valid_pharmacology = {"none", "D1RA", "D2RA", "M1RA"}
    valid_genotypes = {"WT", "CDGIKO", "iSPN-M1RKO"}

    # Note: We don't enforce validation to allow flexibility, but these are the expected values

    # Build the canonical session ID
    session_id = f"{fig}++{compartment}++{measurement}++{spn_type}++{state}++{pharmacology}++{genotype}++{timestamp}"

    return session_id


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
