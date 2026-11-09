"""
Condition names and shared utility functions used in the Zhai 2025 paper conversion scripts.
"""

import argparse

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

# Pharmacology additions mapping - keys match session ID pharm tokens
PHARMACOLOGY_ADDITIONS = {
    "D1RaSch": "SCH-23390: D1 receptor antagonist (bath application) to block direct pathway during LID on-state.",
    "D2RaSul": "Sulpiride: D2 receptor antagonist (bath application) to investigate the role of D2 receptor signaling in LID-induced excitability changes in iSPNs.",
    "M1RaTri": "M1 muscarinic receptor antagonist: Trihexyphenidyl hydrochloride (THP, 3 mg/kg i.p.) administered to assess M1R contribution to LID-induced changes.",
}


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


def surgery_and_pharmacology(state: str, general_metadata: dict) -> tuple[str, str]:
    """Surgery and pharmacology text for a group of mice, given its session ID state token.

    general_metadata.yaml describes the lesioned, levodopa-treated mice. Unlesioned controls had neither; lesioned
    controls had the lesion but no levodopa, so their pharmacology stops before the levodopa priming.
    """
    surgery = general_metadata["NWBFile"]["surgery"]
    pharmacology = general_metadata["NWBFile"]["pharmacology"]
    if state == "UnlesionedControl":
        return "No 6-OHDA lesion (unlesioned control).", "No 6-OHDA and no levodopa (unlesioned control)."
    if state == "LesionedControl":
        return surgery, pharmacology[: pharmacology.index("L-DOPA priming")].rstrip() + " No levodopa."
    return surgery, pharmacology


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
