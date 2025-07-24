"""
Utility functions for extracting conditions from NWB session IDs.

This module provides hardcoded dictionaries for experimental condition mappings
used in the Zhai 2025 paper conversion scripts.
"""

import re
from typing import Optional

# Mapping from folder path names to paper condition names
FOLDER_TO_PAPER_CONDITION = {
    # Figure 1 conditions
    "LID off-state": "LID off-state",
    "LID on-state": "LID on-state",
    "LID on-state with SCH": "LID on-state with SCH",
    # Figure 3 conditions
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
}

# Create reverse mapping for quick lookup
CAMELCASE_TO_PAPER_CONDITION = {formats["CamelCase"]: condition for condition, formats in format_condition.items()}


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
