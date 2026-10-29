"""Mouse-day merger: walks raw-data folders and produces per-mouse-day NWB files.

For now the merger only knows Figure 1 (dSPN somatic + dendritic excitability across
LID off, LID on, and LID on+SCH conditions). The folder conventions for other patch-clamp
figures differ in surface detail (some use MMDDYYYY_N, others MMDD<letter>, Figure 8 uses
YYYYMMDD<letter>); add per-figure handlers in `_figure_specs` to extend coverage.

See `obsidian_docs/plan_per_cell_merge_with_full_icephys_chain.md` for the design.
"""

from __future__ import annotations

import re
from collections import defaultdict
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Iterable

from surmeier_lab_to_nwb.zhai2025.subject_registry import (
    SubjectEntry,
    index_by_key,
    load_registry,
)

REPO_ROOT = Path(__file__).resolve().parents[3]
RAW_DATA_ROOT = REPO_ROOT / "link_to_raw_data"

# Map condition subfolder names (as they appear in the lab's tree) to the
# Data-Connections category token. Different figures use slightly different
# strings; this normalization keeps the matcher figure-agnostic.
_CONDITION_TO_CATEGORY = {
    "LID off-state": "LID OFF",
    "LID on-state": "LID ON",
    "LID on-state with SCH": "LID ON",  # baseline category; SCH is per-recording
    "LID on-state with sul (iSPN)": "LID ON",  # iSPN sulpiride paired-recording
    "LID on-state with sul": "LID ON",
    "KO off-state": "LID OFF",
    "KO on-state": "LID ON",
    "control": "CONT",
    "M1R antagonist": "M1R ant",
    "M1R CRISPR": "M1R CRSPR",
    "interleaved control": "CONT",
    "UL control": "CONT",  # GRAB-ACh, unlesioned
    "PD": "6-OHDA",  # GRAB-ACh, lesioned without levodopa
    "LID off": "LID OFF",  # GRAB-ACh
}

# Per-figure metadata: the Data Connections "animal" reporter token and the
# directory under link_to_raw_data where this figure's data lives.
_FIGURE_SPECS = {
    "F1": {
        "animal": "tD-Tomato",
        "cell_type": "dSPN",
        "raw_root": "Figure 1",
        "somatic_subdir": "Somatic excitability",
        "dendritic_subdir": "Dendritic excitability",
        "date_format": "MMDDYYYY_N",
        "trust_folder_category": True,  # Some 2019 cells not in spreadsheet
    },
    "F3": {
        "animal": "eGFP",
        "cell_type": "iSPN",
        "raw_root": "Figure 3",
        "somatic_subdir": "Somatic excitability",
        "dendritic_subdir": "Dendritic excitability",
        "date_format": "MMDDYYYY_N",
        "trust_folder_category": True,  # Some 2019/2023+ cells not in spreadsheet
    },
    "F6": {
        "animal": "eGFP",
        "cell_type": "iSPN",
        "raw_root": "Figure 6",
        "somatic_subdir": "Somatic excitability",
        "dendritic_subdir": "Dendritic excitability",
        "date_format": "MMDD_letter",
        "trust_folder_category": True,  # 2023 M1Rant cells not in spreadsheet
    },
    "F7": {
        "animal": "CDGIko X eGFP",
        "cell_type": "iSPN",
        "raw_root": "Figure 7",
        "somatic_subdir": "KO SE on vs off",
        "dendritic_subdir": "KO DE on vs off",
        "date_format": "MMDD_letter",
        # The Data Connections spreadsheet's category labels for the CDGI-KO entries
        # disagree with the paper's published state labels. Cross-checked against the
        # paper's Fig 7C panel n's:
        #   - Paper says "off-state, n = 13 cells from 5 mice" → matches the
        #     `KO off-state/` folder (5 mice × 13 cells: May 25, 26, 27, Jun 16, 19)
        #   - Paper says "on-state, n = 8 cells from 4 mice" → matches the
        #     `KO on-state/` folder (4 mice × 8 cells: Aug 16, Sep 1, 2, 13)
        # The folder names are authoritative for state. Spreadsheet rows for these
        # dates carry inverted category labels (likely a lab-internal convention
        # for the KO cohort that we don't need to follow). The matcher resolves
        # subject_id from (date, animal) alone for these figures rather than
        # requiring the spreadsheet's category to match the folder.
        "trust_folder_category": True,
    },
    "F8": {
        "animal": "Adora2-Cre",
        "cell_type": "iSPN",
        "raw_root": "Figure 8",
        "somatic_subdir": "M1R CRISPR SE",
        "dendritic_subdir": None,  # Fig 8 has no dendritic excitability
        "date_format": "YYYYMMDD_letter",
    },
    "F2": {
        "animal": "tD-Tomato",
        "cell_type": "dSPN",
        "raw_root": "Figure 2_SF1A",
        "somatic_subdir": None,
        "dendritic_subdir": None,
        "oepsc_subdir": "Sr-oEPSC",
    },
    "F4": {
        "animal": "eGFP",
        "cell_type": "iSPN",
        "raw_root": "Figure 4_SF1B_SF5",
        "somatic_subdir": None,
        "dendritic_subdir": None,
        "oepsc_subdir": "Sr-oEPSC",
    },
    "F5": {
        "animal": "C57bl/6J",
        "raw_root": "Figure 5_SF2",
        "grab_ach": True,  # condition folders sit directly under raw_root
    },
}


@dataclass
class SweepRecord:
    """One sweep folder on disk, parsed into structured fields."""

    modality: str  # "soma_excitability", "dend_excitability", etc.
    cell_index: int | None  # 1, 2, 3, 4 within the mouse-day; None for a GRAB-ACh recording, which has a region
    condition_subfolder: str  # raw lab string, e.g., "LID off-state"
    recording_path: Path  # the leaf folder containing the .xml + raw files
    pv_scan_date: datetime  # absolute wall-clock from PVScan @date
    dendrite_location: str | None = None  # "prox", "dist", or None for somatic
    location_number: int | None = None  # e.g., 1 for "prox1", 2 for "dist2"
    site_label: str | None = None  # site as written in the folder name, e.g., "prox1", "realprox1", "dist1real"
    suffix: str | None = None  # e.g., "SCH" for cell1_SCH-NNN folders, "quin_single" for a GRAB-ACh recording
    led_setting: float | None = None  # LED setting written in a Sr2+-oEPSC folder name, e.g., 20 for cell1_LED20-001
    region: tuple[int, int] | None = None  # (slice, ROI) of a GRAB-ACh recording


@dataclass
class MouseDayBundle:
    """All sweeps from one mouse-day, grouped by cell."""

    subject_entry: SubjectEntry
    # cell_index -> list of SweepRecord; (slice, ROI) -> list of SweepRecord for GRAB-ACh regions
    cells: dict[int | tuple[int, int], list[SweepRecord]] = field(default_factory=lambda: defaultdict(list))

    @property
    def date_compact(self) -> str:
        return self.subject_entry.date_compact

    def cell_indices(self) -> list[int]:
        return sorted(self.cells.keys())

    def sweeps_for(self, cell_index: int, modality: str | None = None) -> list[SweepRecord]:
        sweeps = self.cells.get(cell_index, [])
        if modality is not None:
            sweeps = [s for s in sweeps if s.modality == modality]
        return sorted(sweeps, key=lambda s: s.pv_scan_date)

    @property
    def session_start_time(self) -> datetime:
        all_sweeps = [s for sweeps in self.cells.values() for s in sweeps]
        if not all_sweeps:
            raise ValueError(f"MouseDayBundle for {self.subject_entry.subject_id} has no sweeps")
        return min(s.pv_scan_date for s in all_sweeps)

    def has_paired_conditions(self) -> bool:
        """True if any cell has sweeps from multiple condition subfolders (e.g., ON + ON-with-SCH)."""
        for sweeps in self.cells.values():
            condition_set = {s.condition_subfolder for s in sweeps}
            if len(condition_set) > 1:
                return True
        return False


# --------- folder name parsers ---------------------------------------------------

# Somatic Figure 1/3 convention: MMDDYYYY_N (date then underscore + cell index)
_RE_SOMATIC_MMDDYYYY = re.compile(r"^(\d{2})(\d{2})(\d{4})_(\d+)$")

# Letter convention: MMDD<letter>[<digit>] — used for dendritic everywhere, plus
# Fig 6/7 somatic. We only match the OUTER MMDD<letter> (the dendritic-location
# digit, if any, is inside the inner folder).
_RE_LETTER_OUTER = re.compile(r"^(\d{2})(\d{2})([a-d])$")
# Variant with full year embedded: MMDDYYYY<letter> (Fig 3 ON LID 2019 cohort).
_RE_LETTER_OUTER_FULLYEAR = re.compile(r"^(\d{2})(\d{2})(\d{4})([a-d])$")

# Dendritic LOCATION folder inside a cell-day folder: MMDD<letter><digit>
_RE_DENDRITIC_LOCATION = re.compile(r"^(\d{2})(\d{2})[a-d](\d+)$")

# Dendritic LEAF recording: MMDDYYYY_Cell<N>_<location><N>_trio-<NNN>
_RE_DENDRITIC_LEAF = re.compile(r"^(\d{2})(\d{2})(\d{4})_Cell(\d+)_(prox|dist)(\d+)_trio-(\d+)$")
# Dendritic LEAF recording variant: YYYYMMDD_Cell<N>_<location><N>_trio-<NNN>
_RE_DENDRITIC_LEAF_YYYYMMDD = re.compile(r"^(\d{4})(\d{2})(\d{2})_Cell(\d+)_(prox|dist)(\d+)_trio-(\d+)$")

# Fig 8 outer: YYYYMMDD<letter>
_RE_YYYYMMDD_LETTER_OUTER = re.compile(r"^(\d{4})(\d{2})(\d{2})([a-d])$")


# Flexible dendritic leaf pattern: handles all observed lab variations.
#   - Date in MMDDYYYY or YYYYMMDD order
#   - "Cell" or "cell" (case-insensitive)
#   - "trio" or "Trio" (case-insensitive)
#   - Location: prox, dist, with "real" before or after the site (realprox1, dist1real) for a
#     second patch of dendrite the lab recorded at that site after moving the microscope
#   - Optional trailing pharm suffix: `_trio_SCH-NNN`, `_trio_sul-NNN`
#   - "Cell<N>_<loc><N>_trio[_pharm]" OR "cell<N>_trio_<loc><N>" (lab swapped order in 2023)
#   - Optional trailing digit on the location (some leaves have `trio_dist` without digit)
#   - `_trio` missing when nothing else follows the site (2016-06-08 cell 3 prox1, 2017-06-19 cell 2 dist1)
#   - `_bsl` instead of `_trio`: a background scan, the line moved off the dendrite and no stimulus
#     (2017-05-25 and 2017-05-26 cell 2 only)
_RE_DENDRITIC_LEAF_FLEX = re.compile(
    # Date: 8 or 9 digits (9 allows for lab typos like `007182016`).
    r"^(?P<date>\d{8,9})_[Cc]ell(?P<cell>\d+)_"
    r"(?:"
    r"(?P<site1>(?P<loc1>real(?:prox|dist)|prox|dist)(?P<num1>\d+)(?:real)?)(?:(?P<background>_bsl)|_[Tt]rio(?:_[A-Za-z]+)?)?"  # Cell1_dist1_trio, _trio_SCH or _bsl
    r"|[Tt]rio_(?P<site2>(?P<loc2>real(?:prox|dist)|prox|dist)(?P<num2>\d*))"  # Cell1_trio_dist1
    r")"
    r"-(?P<seq>\d+)$"
)


def _match_dendritic_leaf(name: str):
    """Match a dendritic leaf folder name across all observed lab naming variations.

    Returns a tuple (cell_in_name, loc_type, loc_num, site_label, modality) or None if no match;
    modality is "dend_background" for a `_bsl` background scan and "dend_excitability" otherwise.
    The loc_type is normalized: "realprox" -> "prox", "realdist" -> "dist". The lab's
    "real" marks a second patch of dendrite recorded at the site after moving the
    microscope; it is kept in site_label, the site as written in the folder name.
    """
    # The lab misspelled prox2 as proxt2 in the folders of 2019-05-08 cell 3 (Figure 3, on-state with sulpiride).
    m = _RE_DENDRITIC_LEAF_FLEX.match(name.replace("_proxt2_", "_prox2_"))
    if m:
        cell = int(m.group("cell"))
        raw_loc = m.group("loc1") or m.group("loc2")
        # Normalize "realprox"/"realdist" -> "prox"/"dist"
        loc_type = raw_loc.replace("real", "")
        num_str = m.group("num1") or m.group("num2") or "1"
        loc_num = int(num_str) if num_str else 1
        site_label = m.group("site1") or m.group("site2")
        modality = "dend_background" if m.group("background") else "dend_excitability"
        return cell, loc_type, loc_num, site_label, modality
    return None


# Somatic LEAF recording: cell<N>-<NNN> or cell<N>_<suffix>-<NNN> (e.g., SCH, sul, oxoM)
_RE_SOMATIC_LEAF = re.compile(r"^cell(\d+)(_[A-Za-z]+)?-(\d+)$")

# Sr2+-oEPSC LEAF recording: cell<N>_LED<setting>-<NNN>. Observed variations: no cell
# (08162023a/LED20-NNN), no underscore (cell1LED20), "int" before the setting (cell1_LEDint14),
# and "b", "_1mspulse" or "_real" after it (cell3_LED100b, cell1_LED100_1mspulse, cell2_LED4_real).
_RE_OEPSC_LEAF = re.compile(r"^(?:cell(\d+)_?)?LED(?:int)?(\d+(?:\.\d+)?)(?:b|_1mspulse|_real)?-(\d+)$")

# GRAB-ACh region folder: MMDDYYYYslice<N>ROI<M>, or MMDDYYYYslice<N> on the May and June days, which filmed one
# region per slice.
_RE_GRAB_ACH_REGION = re.compile(r"^(\d{8})slice(\d+)(?:ROI(\d+))?$")
# GRAB-ACh recording: BOT_MMDDYYYY_<label>_<drug and stimulus>-NNN, e.g. BOT_04022024_slice1ROI1_quin_single-001.
_RE_GRAB_ACH_LEAF = re.compile(r"^BOT_(\d{8})_([^_]+)_(.+)-(\d+)$")
# The label names the slice and ROI (slice1ROI2), ROI "A" (slice1A, the only ROI of the slice) or the slice alone.
_RE_GRAB_ACH_LABEL = re.compile(r"^slice(\d+)(?:ROI(\d+)|(A))?$")
# Labels that disagree with their folder, placed with the lab's Data Connections sheet and the stage position in
# each recording's XML.
_GRAB_ACH_LABEL_FIXES = {
    # 2024-04-12, a TTX recording in the slice 1 ROI 2 folder, at the ROI 2 position.
    "slice1ROI22": "slice1ROI2",
    # 2024-05-31, two 50 nM dopamine recordings at the position of slice1A, 3 minutes after its last one; slice1A has
    # no others of that kind, and the sheet counts 2 regions that day.
    "slice1realA": "slice1A",
}
# GRAB-ACh recordings left out on purpose.
_GRAB_ACH_EXCLUDED = {
    # A false start: filed under ROI 1 but recorded at the ROI 2 position, stopped after 2.9 s, before its stimulus
    # at 10 s, and restarted 8 s later at the same spot as BOT_04022024_slice1ROI2_sul_sinlge-001. Only baseline.
    "BOT_04022024_slice1ROI_sul_sinlge-001",
}

_LETTER_TO_INDEX = {"a": 1, "b": 2, "c": 3, "d": 4}


def _parse_pv_scan_date(xml_path: Path) -> datetime | None:
    """Extract the PVScan @date attribute from a master XML at 1-second precision.

    The attribute format is "M/D/YYYY H:MM:SS AM/PM" (e.g., "2/2/2017 3:57:24 PM").
    Returns None if the file can't be parsed.
    """
    try:
        # Read only enough to find the PVScan tag; full parse is wasteful.
        text = xml_path.read_text(errors="replace", encoding="utf-8")[:2048]
    except OSError:
        return None
    m = re.search(r'<PVScan[^>]*\bdate="([^"]+)"', text)
    if not m:
        return None
    raw = m.group(1)
    for fmt in ("%m/%d/%Y %I:%M:%S %p", "%m/%d/%Y %H:%M:%S"):
        try:
            return datetime.strptime(raw, fmt)
        except ValueError:
            continue
    return None


def _walk_somatic_mmddyyyy_n(
    figure_root: Path,
    somatic_subdir: str,
    condition_subfolders: Iterable[str],
) -> list[SweepRecord]:
    """Walk somatic excitability folders using the MMDDYYYY_N convention.

    Used by Figures 1 and 3. Convention:
        <figure_root>/<somatic_subdir>/<condition>/MMDDYYYY_N/cellN-NNN/
    """
    out: list[SweepRecord] = []
    somatic_root = figure_root / somatic_subdir
    if not somatic_root.exists():
        return out
    for condition in condition_subfolders:
        cond_dir = somatic_root / condition
        if not cond_dir.exists():
            continue
        for cell_day_dir in sorted(cond_dir.iterdir()):
            if not cell_day_dir.is_dir():
                continue
            m = _RE_SOMATIC_MMDDYYYY.match(cell_day_dir.name)
            if not m:
                raise ValueError(f"Unrecognized somatic cell folder: {cell_day_dir}")
            mm, dd, yyyy, cell_index = m.groups()
            cell_index = int(cell_index)
            for leaf in sorted(cell_day_dir.iterdir()):
                if not leaf.is_dir():
                    continue
                leaf_match = _RE_SOMATIC_LEAF.match(leaf.name)
                if not leaf_match:
                    raise ValueError(f"Unrecognized somatic recording folder: {leaf}")
                _leaf_cell, suffix_part, _seq = leaf_match.groups()
                suffix = suffix_part.lstrip("_") if suffix_part else None
                xml_path = leaf / f"{leaf.name}.xml"
                pv_date = _parse_pv_scan_date(xml_path) if xml_path.exists() else None
                if pv_date is None:
                    raise ValueError(f"No recording time in {xml_path}")
                out.append(
                    SweepRecord(
                        modality="soma_excitability",
                        cell_index=cell_index,
                        condition_subfolder=condition,
                        recording_path=leaf,
                        pv_scan_date=pv_date,
                        suffix=suffix,
                    )
                )
    return out


def _walk_somatic_mmdd_letter(
    figure_root: Path,
    somatic_subdir: str,
    condition_subfolders: Iterable[str],
) -> list[SweepRecord]:
    """Walk somatic excitability folders using the MMDD<letter>/cellN-NNN convention.

    Used by Figures 6 and 7. Convention:
        <figure_root>/<somatic_subdir>/<condition>/MMDD<letter>/cellN-NNN/

    The LEAF `cellN` is authoritative for the cell index. The outer letter is a
    folder-ordering convention that usually but not always matches (e.g., 0616b
    contains `cell3-NNN` in Fig 7 KO data).
    """
    out: list[SweepRecord] = []
    somatic_root = figure_root / somatic_subdir
    if not somatic_root.exists():
        return out
    for condition in condition_subfolders:
        cond_dir = somatic_root / condition
        if not cond_dir.exists():
            continue
        for cell_day_dir in sorted(cond_dir.iterdir()):
            if not cell_day_dir.is_dir():
                continue
            if not _RE_LETTER_OUTER.match(cell_day_dir.name):
                raise ValueError(f"Unrecognized somatic cell folder: {cell_day_dir}")
            for leaf in sorted(cell_day_dir.iterdir()):
                if not leaf.is_dir():
                    continue
                leaf_match = _RE_SOMATIC_LEAF.match(leaf.name)
                if not leaf_match:
                    raise ValueError(f"Unrecognized somatic recording folder: {leaf}")
                leaf_cell, suffix_part, _seq = leaf_match.groups()
                cell_index = int(leaf_cell)
                suffix = suffix_part.lstrip("_") if suffix_part else None
                xml_path = leaf / f"{leaf.name}.xml"
                pv_date = _parse_pv_scan_date(xml_path) if xml_path.exists() else None
                if pv_date is None:
                    raise ValueError(f"No recording time in {xml_path}")
                out.append(
                    SweepRecord(
                        modality="soma_excitability",
                        cell_index=cell_index,
                        condition_subfolder=condition,
                        recording_path=leaf,
                        pv_scan_date=pv_date,
                        suffix=suffix,
                    )
                )
    return out


def _walk_somatic_yyyymmdd_letter(
    figure_root: Path,
    somatic_subdir: str,
    condition_subfolders: Iterable[str],
) -> list[SweepRecord]:
    """Walk somatic excitability folders using the YYYYMMDD<letter>/cellN-NNN convention.

    Used by Figure 8 (the only figure with 8-digit date prefixes).
    """
    out: list[SweepRecord] = []
    somatic_root = figure_root / somatic_subdir
    if not somatic_root.exists():
        return out
    for condition in condition_subfolders:
        cond_dir = somatic_root / condition
        if not cond_dir.exists():
            continue
        for cell_day_dir in sorted(cond_dir.iterdir()):
            if not cell_day_dir.is_dir():
                continue
            outer = _RE_YYYYMMDD_LETTER_OUTER.match(cell_day_dir.name)
            if not outer:
                raise ValueError(f"Unrecognized somatic cell folder: {cell_day_dir}")
            _yyyy, _mm, _dd, letter = outer.groups()
            cell_index = _LETTER_TO_INDEX.get(letter)
            if cell_index is None:
                continue
            for leaf in sorted(cell_day_dir.iterdir()):
                if not leaf.is_dir():
                    continue
                leaf_match = _RE_SOMATIC_LEAF.match(leaf.name)
                if not leaf_match:
                    raise ValueError(f"Unrecognized somatic recording folder: {leaf}")
                leaf_cell, suffix_part, _seq = leaf_match.groups()
                if int(leaf_cell) != cell_index:
                    raise ValueError(f"{leaf} names cell {leaf_cell} but is in the folder of cell {cell_index}")
                suffix = suffix_part.lstrip("_") if suffix_part else None
                xml_path = leaf / f"{leaf.name}.xml"
                pv_date = _parse_pv_scan_date(xml_path) if xml_path.exists() else None
                if pv_date is None:
                    raise ValueError(f"No recording time in {xml_path}")
                out.append(
                    SweepRecord(
                        modality="soma_excitability",
                        cell_index=cell_index,
                        condition_subfolder=condition,
                        recording_path=leaf,
                        pv_scan_date=pv_date,
                        suffix=suffix,
                    )
                )
    return out


def _walk_oepsc_mmddyyyy_letter(
    figure_root: Path,
    oepsc_subdir: str,
    condition_subfolders: Iterable[str],
) -> list[SweepRecord]:
    """Walk Sr2+-oEPSC folders using the MMDDYYYY<letter>/cellN_LED<setting>-NNN convention.

    Used by Figures 2 and 4. The letter gives the cell (a = cell 1); a leaf that names a cell must agree.
    """
    out: list[SweepRecord] = []
    oepsc_root = figure_root / oepsc_subdir
    for condition in condition_subfolders:
        cond_dir = oepsc_root / condition
        if not cond_dir.exists():
            continue
        for cell_day_dir in sorted(cond_dir.iterdir()):
            if not cell_day_dir.is_dir():
                continue
            outer = _RE_LETTER_OUTER_FULLYEAR.match(cell_day_dir.name)
            if not outer:
                raise ValueError(f"Unrecognized Sr2+-oEPSC cell folder: {cell_day_dir}")
            cell_index = _LETTER_TO_INDEX[outer.group(4)]
            for leaf in sorted(cell_day_dir.iterdir()):
                if not leaf.is_dir():
                    continue
                leaf_match = _RE_OEPSC_LEAF.match(leaf.name)
                if not leaf_match:
                    raise ValueError(f"Unrecognized Sr2+-oEPSC recording folder: {leaf}")
                leaf_cell, led_setting, _seq = leaf_match.groups()
                if leaf_cell is not None and int(leaf_cell) != cell_index:
                    raise ValueError(f"{leaf} names cell {leaf_cell} but is in the folder of cell {cell_index}")
                xml_path = leaf / f"{leaf.name}.xml"
                pv_date = _parse_pv_scan_date(xml_path) if xml_path.exists() else None
                if pv_date is None:
                    raise ValueError(f"No recording time in {xml_path}")
                out.append(
                    SweepRecord(
                        modality="oepsc",
                        cell_index=cell_index,
                        condition_subfolder=condition,
                        recording_path=leaf,
                        pv_scan_date=pv_date,
                        led_setting=float(led_setting),
                    )
                )
    return out


def _walk_grab_ach(figure_root: Path) -> list[SweepRecord]:
    """Walk GRAB-ACh folders: <condition>/MMDDYYYYslice<N>[ROI<M>]/BOT_MMDDYYYY_<label>_<drug and stimulus>-NNN.

    Used by Figure 5. The folder gives the region (slice, ROI), ROI 1 when it names only the slice; the label in a
    recording's name must agree with it.
    """
    out: list[SweepRecord] = []
    for cond_dir in sorted(p for p in figure_root.iterdir() if p.is_dir()):
        for region_dir in sorted(p for p in cond_dir.iterdir() if p.is_dir()):
            outer = _RE_GRAB_ACH_REGION.match(region_dir.name)
            if not outer:
                raise ValueError(f"Unrecognized GRAB-ACh region folder: {region_dir}")
            folder_date, slice_number, roi_number = outer.groups()
            region = (int(slice_number), int(roi_number or 1))
            for leaf in sorted(p for p in region_dir.iterdir() if p.is_dir()):
                if leaf.name in _GRAB_ACH_EXCLUDED:
                    continue
                leaf_match = _RE_GRAB_ACH_LEAF.match(leaf.name)
                if not leaf_match:
                    raise ValueError(f"Unrecognized GRAB-ACh recording folder: {leaf}")
                leaf_date, label, drug_and_stimulus, _seq = leaf_match.groups()
                label_match = _RE_GRAB_ACH_LABEL.match(_GRAB_ACH_LABEL_FIXES.get(label, label))
                if not label_match:
                    raise ValueError(f"Unrecognized GRAB-ACh recording folder: {leaf}")
                leaf_slice, leaf_roi, leaf_letter = label_match.groups()
                # ROI "A" is ROI 1; a label without a ROI takes the folder's.
                leaf_region = (int(leaf_slice), int(leaf_roi) if leaf_roi else 1 if leaf_letter else region[1])
                if leaf_date != folder_date:
                    raise ValueError(f"{leaf} is dated {leaf_date} but is in a folder dated {folder_date}")
                if leaf_region != region:
                    raise ValueError(f"{leaf} names {label} but is in the folder of slice {region[0]} ROI {region[1]}")
                xml_path = leaf / f"{leaf.name}.xml"
                pv_date = _parse_pv_scan_date(xml_path) if xml_path.exists() else None
                if pv_date is None:
                    raise ValueError(f"No recording time in {xml_path}")
                out.append(
                    SweepRecord(
                        modality="grab_ach",
                        cell_index=None,
                        condition_subfolder=cond_dir.name,
                        recording_path=leaf,
                        pv_scan_date=pv_date,
                        suffix=drug_and_stimulus,
                        region=region,
                    )
                )
    return out


def _walk_dendritic_mmdd_letter(
    figure_root: Path,
    dendritic_subdir: str,
    condition_subfolders: Iterable[str],
) -> list[SweepRecord]:
    """Walk dendritic excitability folders using the MMDD<letter> convention.

    Used by Figures 1, 3, 6, 7. Two layouts supported:
      Layout A (Fig 1):  <root>/<cond>/MMDD<letter>/MMDD<letter><digit>/MMDDYYYY_Cell<N>_<loc><N>_trio-<NNN>/
      Layout B (Fig 3+): <root>/<cond>/MMDD<letter>/MMDDYYYY_Cell<N>_<loc><N>_trio-<NNN>/
    The walker auto-detects which layout each (condition, cell_day) folder uses.
    """
    out: list[SweepRecord] = []
    dend_root = figure_root / dendritic_subdir
    if not dend_root.exists():
        return out
    for condition in condition_subfolders:
        cond_dir = dend_root / condition
        if not cond_dir.exists():
            continue
        for cell_day_dir in sorted(cond_dir.iterdir()):
            if not cell_day_dir.is_dir():
                continue
            # Layout C: outer name is MMDD<letter><digit> (e.g., 0411a2 = April 11,
            # Cell a, location 2). The leaves are dendritic recordings directly inside.
            # The lab uses this shortcut for some paired-recording cells (Fig 1 SCH, etc.).
            outer_loc = _RE_DENDRITIC_LOCATION.match(cell_day_dir.name)
            if outer_loc:
                # The cell-day is implicit; walk leaves directly via the same _emit_leaf
                # logic used below. Skip the standard _RE_LETTER_OUTER check.
                def _emit_leaf(leaf_dir: Path):
                    parsed = _match_dendritic_leaf(leaf_dir.name)
                    if parsed is None:
                        raise ValueError(f"Unrecognized dendritic recording folder: {leaf_dir}")
                    cell_in_name, loc_type, loc_num, site_label, modality = parsed
                    xml_path = leaf_dir / f"{leaf_dir.name}.xml"
                    pv_date = _parse_pv_scan_date(xml_path) if xml_path.exists() else None
                    if pv_date is None:
                        raise ValueError(f"No recording time in {xml_path}")
                    out.append(
                        SweepRecord(
                            modality=modality,
                            cell_index=cell_in_name,
                            condition_subfolder=condition,
                            recording_path=leaf_dir,
                            pv_scan_date=pv_date,
                            dendrite_location=loc_type,
                            location_number=loc_num,
                            site_label=site_label,
                        )
                    )

                for grandchild in sorted(cell_day_dir.iterdir()):
                    if grandchild.is_dir():
                        _emit_leaf(grandchild)
                continue  # don't fall through to the LETTER_OUTER branch
            # Match either MMDD<letter> or MMDDYYYY<letter> for the cell-day folder.
            outer = _RE_LETTER_OUTER.match(cell_day_dir.name) or _RE_LETTER_OUTER_FULLYEAR.match(cell_day_dir.name)
            if not outer:
                raise ValueError(f"Unrecognized dendritic cell folder: {cell_day_dir}")
            # The letter is always the LAST captured group regardless of which regex matched.
            letter = outer.groups()[-1]
            cell_index = _LETTER_TO_INDEX.get(letter)
            if cell_index is None:
                continue

            # Two layouts: with intermediate MMDD<letter><digit> folders (Fig 1)
            # or without (Fig 3+). We walk both possibilities by inspecting each child:
            # if it matches the leaf pattern directly, treat it as a leaf;
            # if it matches the location pattern, recurse into it.
            def _emit_leaf(leaf_dir: Path):
                parsed = _match_dendritic_leaf(leaf_dir.name)
                if parsed is None:
                    raise ValueError(f"Unrecognized dendritic recording folder: {leaf_dir}")
                cell_in_name, loc_type, loc_num, site_label, modality = parsed
                xml_path = leaf_dir / f"{leaf_dir.name}.xml"
                pv_date = _parse_pv_scan_date(xml_path) if xml_path.exists() else None
                if pv_date is None:
                    raise ValueError(f"No recording time in {xml_path}")
                out.append(
                    SweepRecord(
                        modality=modality,
                        cell_index=cell_in_name,
                        condition_subfolder=condition,
                        recording_path=leaf_dir,
                        pv_scan_date=pv_date,
                        dendrite_location=loc_type,
                        location_number=loc_num,
                        site_label=site_label,
                    )
                )

            for child in sorted(cell_day_dir.iterdir()):
                if not child.is_dir():
                    continue
                if _match_dendritic_leaf(child.name) is not None:
                    # Layout B (Fig 3+): leaf directly under cell-day folder.
                    _emit_leaf(child)
                elif _RE_DENDRITIC_LOCATION.match(child.name):
                    # Layout A (Fig 1): intermediate location folder; recurse.
                    for grandchild in sorted(child.iterdir()):
                        if grandchild.is_dir():
                            _emit_leaf(grandchild)
                else:
                    raise ValueError(f"Unrecognized dendritic recording folder: {child}")
    return out


def build_figure_bundles(figure_code: str = "F1") -> list[MouseDayBundle]:
    """Build per-mouse-day bundles for one figure.

    Walks raw-data folders, parses PVScan @date for each sweep, joins against the
    Data Connections subject registry to assign a Subject per mouse-day.
    """
    spec = _FIGURE_SPECS[figure_code]
    fig_root = RAW_DATA_ROOT / spec["raw_root"]
    if not fig_root.exists():
        raise FileNotFoundError(f"Figure raw-data root not found: {fig_root}")

    # All condition subfolders found under either the somatic or dendritic dir.
    condition_subfolders: list[str] = []
    for modality_dir in (spec.get("somatic_subdir"), spec.get("dendritic_subdir"), spec.get("oepsc_subdir")):
        if modality_dir is None:
            continue
        modality_root = fig_root / modality_dir
        if modality_root.exists():
            for sub in modality_root.iterdir():
                if sub.is_dir() and sub.name not in condition_subfolders:
                    condition_subfolders.append(sub.name)

    sweeps: list[SweepRecord] = []
    somatic_subdir = spec.get("somatic_subdir")
    dendritic_subdir = spec.get("dendritic_subdir")
    date_format = spec.get("date_format", "MMDDYYYY_N")
    if somatic_subdir is not None:
        if date_format == "MMDDYYYY_N":
            sweeps.extend(_walk_somatic_mmddyyyy_n(fig_root, somatic_subdir, condition_subfolders))
        elif date_format == "MMDD_letter":
            sweeps.extend(_walk_somatic_mmdd_letter(fig_root, somatic_subdir, condition_subfolders))
        elif date_format == "YYYYMMDD_letter":
            sweeps.extend(_walk_somatic_yyyymmdd_letter(fig_root, somatic_subdir, condition_subfolders))
        else:
            raise NotImplementedError(f"date_format={date_format!r} not yet supported")
    if dendritic_subdir is not None:
        sweeps.extend(_walk_dendritic_mmdd_letter(fig_root, dendritic_subdir, condition_subfolders))
    oepsc_subdir = spec.get("oepsc_subdir")
    if oepsc_subdir is not None:
        sweeps.extend(_walk_oepsc_mmddyyyy_letter(fig_root, oepsc_subdir, condition_subfolders))
    if spec.get("grab_ach"):
        sweeps.extend(_walk_grab_ach(fig_root))

    # Load the subject registry and index by (date, animal, category).
    entries = load_registry()
    indexed = index_by_key(entries)
    # For figures where the spreadsheet's category labels are unreliable (see
    # `trust_folder_category` in `_FIGURE_SPECS`), also build a (date, animal)
    # index so we can resolve subject_id without requiring category agreement.
    indexed_by_date_animal: dict[tuple[str, str], SubjectEntry] = {}
    for entry in entries:
        if entry.section == "VIDEOS":
            continue
        key = (entry.date_compact, entry.animal)
        if key not in indexed_by_date_animal:
            indexed_by_date_animal[key] = entry

    trust_folder = spec.get("trust_folder_category", False)

    # Group sweeps by (date_compact, cell_index), then look up the Subject.
    bundles_by_date: dict[str, MouseDayBundle] = {}
    unmatched: list[SweepRecord] = []
    for sw in sweeps:
        date_compact = sw.pv_scan_date.strftime("%Y%m%d")
        folder_category = _CONDITION_TO_CATEGORY.get(sw.condition_subfolder)
        if folder_category is None:
            unmatched.append(sw)
            continue
        if trust_folder:
            # Use (date, animal) lookup; ignore spreadsheet's category column.
            # If the spreadsheet's PVScan date doesn't match any row (e.g., the
            # spreadsheet has year-discrepant entries for Fig 7's June 2016 vs 2017
            # CDGI-KO cells), fall back to constructing the Subject from the folder
            # data directly. The folder name + PVScan date + animal genotype are
            # sufficient to identify the mouse-day; the spreadsheet is a useful
            # source of cells_per_experiment metadata when present but not required.
            from surmeier_lab_to_nwb.zhai2025.subject_registry import _build_subject_id

            base_entry = indexed_by_date_animal.get((date_compact, spec["animal"]))
            subject_id = _build_subject_id(date_compact, folder_category, spec["animal"])
            if base_entry is not None:
                subject = SubjectEntry(
                    subject_id=subject_id,
                    year=base_entry.year,
                    month=base_entry.month,
                    day=base_entry.day,
                    animal=base_entry.animal,
                    category=folder_category,
                    section=base_entry.section,
                    cells_per_experiment=base_entry.cells_per_experiment,
                )
            else:
                # Construct a Subject from folder data alone.
                subject = SubjectEntry(
                    subject_id=subject_id,
                    year=sw.pv_scan_date.year,
                    month=sw.pv_scan_date.month,
                    day=sw.pv_scan_date.day,
                    animal=spec["animal"],
                    category=folder_category,
                    section=f"FOLDER_ONLY_{spec['raw_root']}",
                    cells_per_experiment={},
                )
        else:
            subject_key = (date_compact, spec["animal"], folder_category)
            subject = indexed.get(subject_key)
            if subject is None:
                unmatched.append(sw)
                continue
        bundle = bundles_by_date.get(date_compact)
        if bundle is None:
            bundle = MouseDayBundle(subject_entry=subject)
            bundles_by_date[date_compact] = bundle
        bundle.cells[sw.cell_index if sw.region is None else sw.region].append(sw)

    bundles = sorted(bundles_by_date.values(), key=lambda b: b.date_compact)
    if unmatched:
        listing = "\n".join(f"  {sw.recording_path}" for sw in unmatched)
        raise ValueError(f"{len(unmatched)} recording(s) did not match a Data Connections subject:\n{listing}")
    return bundles


if __name__ == "__main__":
    bundles = build_figure_bundles("F1")
    print(f"\nFigure 1 bundles: {len(bundles)} mouse-days")
    for b in bundles[:5]:
        cells = b.cell_indices()
        modalities_per_cell = {c: sorted({s.modality for s in b.sweeps_for(c)}) for c in cells}
        paired = " [PAIRED]" if b.has_paired_conditions() else ""
        print(f"  {b.subject_entry.subject_id}{paired}")
        for c in cells:
            sweeps = b.sweeps_for(c)
            modality_counts = defaultdict(int)
            condition_counts = defaultdict(int)
            for sw in sweeps:
                modality_counts[sw.modality] += 1
                condition_counts[sw.condition_subfolder] += 1
            print(f"    Cell{c}: {dict(modality_counts)}  conditions={dict(condition_counts)}")
    print(f"\nTotal cells across all mouse-days: {sum(len(b.cells) for b in bundles)}")
    print(f"Mouse-days with paired conditions: {sum(1 for b in bundles if b.has_paired_conditions())}")
