"""
Shared utility functions used in the Zhai 2025 paper conversion scripts.
"""

import argparse


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

    metadata.yaml describes the lesioned, levodopa-treated mice. Unlesioned controls had neither; lesioned
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
