"""
Centralized condition mappings for all Zhai 2025 conversion scripts.

This module contains standardized condition mappings used across different
experiment types to ensure consistency in session IDs, table names, and descriptions.
"""

# Centralized condition mapping for all experiments (somatic, dendritic, etc.)
# Maps folder names (as they appear in raw data) to standardized formats
CONDITION_MAPPING = {
    # Figure 1 conditions
    "LID off-state": {
        "camel_case": "LIDOffState",
        "underscore": "lid_off_state",
        "human_readable": "LID off-state",
        "description": "Levodopa-induced dyskinesia animals in off-state (24h after last L-DOPA injection)",
    },
    "LID on-state": {
        "camel_case": "LIDOnState",
        "underscore": "lid_on_state",
        "human_readable": "LID on-state",
        "description": "Levodopa-induced dyskinesia animals in on-state (1h after L-DOPA injection)",
    },
    "LID on-state with SCH": {
        "camel_case": "LIDOnStateWithSchD1Antagonist",
        "underscore": "lid_on_state_with_sch",
        "human_readable": "LID on-state with SCH D1 antagonist",
        "description": "LID on-state with SCH-23390 D1 receptor antagonist to block direct pathway",
    },
    # Figure 3 conditions
    "LID on-state with sul (iSPN)": {
        "camel_case": "LIDOnStateWithSulpiride",
        "underscore": "lid_on_state_with_sul",
        "human_readable": "LID on-state with sulpiride",
        "description": "LID on-state with sulpiride D2 receptor antagonist to investigate indirect pathway",
    },
    # Figure 6 conditions
    "control": {
        "camel_case": "Control",
        "underscore": "control",
        "human_readable": "control",
        "description": "Control condition for M1R antagonist experiments",
    },
    "M1R antagonist": {
        "camel_case": "M1RAntagonist",
        "underscore": "m1r_antagonist",
        "human_readable": "M1R antagonist",
        "description": "M1 muscarinic receptor antagonist (trihexyphenidyl) to test M1R contribution to LID",
    },
    # Figure 7 conditions
    "KO off-state": {
        "camel_case": "KOOffState",
        "underscore": "ko_off_state",
        "human_readable": "CDGI knockout off-state",
        "description": "CDGI knockout mice in off-state to test role of cholinergic interneurons",
    },
    "KO on-state": {
        "camel_case": "KOOnState",
        "underscore": "ko_on_state",
        "human_readable": "CDGI knockout on-state",
        "description": "CDGI knockout mice in on-state to test role of cholinergic interneurons",
    },
    # Figure 8 conditions
    "M1R CRISPR": {
        "camel_case": "M1RCRISPR",
        "underscore": "m1r_crispr",
        "human_readable": "M1R CRISPR knockout",
        "description": "Cell-specific M1R knockout using CRISPR-Cas9 in striatal neurons",
    },
    "interleaved control": {
        "camel_case": "InterleavedControl",
        "underscore": "interleaved_control",
        "human_readable": "interleaved control",
        "description": "Control condition interleaved with M1R CRISPR experiments",
    },
}


def get_condition_mapping(condition: str, format_type: str) -> str:
    """
    Get standardized condition name for experiments (somatic, dendritic, etc.).

    Parameters
    ----------
    condition : str
        The condition name as it appears in the raw data folder
    format_type : str
        The desired format: "camel_case" for session IDs, "underscore" for tables,
        "human_readable" for descriptions, "description" for experimental context

    Returns
    -------
    str
        The standardized condition name in the requested format

    Raises
    ------
    ValueError
        If condition or format_type is not found in the mapping
    """
    if condition not in CONDITION_MAPPING:
        raise ValueError(f"Unknown condition '{condition}'. Available conditions: {list(CONDITION_MAPPING.keys())}")

    if format_type not in ["camel_case", "underscore", "human_readable", "description"]:
        raise ValueError(
            f"Unknown format_type '{format_type}'. Must be 'camel_case', 'underscore', 'human_readable', or 'description'"
        )

    return CONDITION_MAPPING[condition][format_type]
