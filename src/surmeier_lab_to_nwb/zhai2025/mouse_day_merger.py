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
}


@dataclass
class SweepRecord:
    """One sweep folder on disk, parsed into structured fields."""

    modality: str  # "soma_excitability", "dend_excitability", etc.
    cell_index: int  # 1, 2, 3, 4 within the mouse-day
    condition_subfolder: str  # raw lab string, e.g., "LID off-state"
    recording_path: Path  # the leaf folder containing the .xml + raw files
    pv_scan_date: datetime  # absolute wall-clock from PVScan @date
    dendrite_location: str | None = None  # "prox", "dist", or None for somatic
    location_number: int | None = None  # e.g., 1 for "prox1", 2 for "dist2"
    suffix: str | None = None  # e.g., "SCH" for cell1_SCH-NNN folders


@dataclass
class MouseDayBundle:
    """All sweeps from one mouse-day, grouped by cell."""

    subject_entry: SubjectEntry
    # cell_index -> list of SweepRecord
    cells: dict[int, list[SweepRecord]] = field(default_factory=lambda: defaultdict(list))

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
#   - Location: prox, dist, or realprox/realdist (Fig 3 has "realprox" — re-imaged proximal location)
#   - Optional trailing pharm suffix: `_trio_SCH-NNN`, `_trio_sul-NNN`
#   - "Cell<N>_<loc><N>_trio[_pharm]" OR "cell<N>_trio_<loc><N>" (lab swapped order in 2023)
#   - Optional trailing digit on the location (some leaves have `trio_dist` without digit)
_RE_DENDRITIC_LEAF_FLEX = re.compile(
    # Date: 8 or 9 digits (9 allows for lab typos like `007182016`).
    r"^(?P<date>\d{8,9})_[Cc]ell(?P<cell>\d+)_"
    r"(?:"
    r"(?P<loc1>real(?:prox|dist)|prox|dist)(?P<num1>\d+)_[Tt]rio(?:_[A-Za-z]+)?"  # Cell1_dist1_trio or Cell1_dist1_trio_SCH
    r"|[Tt]rio_(?P<loc2>real(?:prox|dist)|prox|dist)(?P<num2>\d*)"  # Cell1_trio_dist1
    r")"
    r"-(?P<seq>\d+)$"
)


def _match_dendritic_leaf(name: str):
    """Match a dendritic leaf folder name across all observed lab naming variations.

    Returns a tuple (cell_in_name, loc_type, loc_num) or None if no match.
    The loc_type is normalized: "realprox" -> "prox", "realdist" -> "dist" (the lab's
    "real" prefix indicates re-imaged after objective re-positioning; same compartment).
    """
    m = _RE_DENDRITIC_LEAF_FLEX.match(name)
    if m:
        cell = int(m.group("cell"))
        raw_loc = m.group("loc1") or m.group("loc2")
        # Normalize "realprox"/"realdist" -> "prox"/"dist"
        loc_type = raw_loc.replace("real", "")
        num_str = m.group("num1") or m.group("num2") or "1"
        loc_num = int(num_str) if num_str else 1
        return cell, loc_type, loc_num
    return None


# Somatic LEAF recording: cell<N>-<NNN> or cell<N>_<suffix>-<NNN> (e.g., SCH, sul, oxoM)
_RE_SOMATIC_LEAF = re.compile(r"^cell(\d+)(_[A-Za-z]+)?-(\d+)$")

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
                continue
            mm, dd, yyyy, cell_index = m.groups()
            cell_index = int(cell_index)
            for leaf in sorted(cell_day_dir.iterdir()):
                if not leaf.is_dir():
                    continue
                leaf_match = _RE_SOMATIC_LEAF.match(leaf.name)
                if not leaf_match:
                    continue
                _leaf_cell, suffix_part, _seq = leaf_match.groups()
                suffix = suffix_part.lstrip("_") if suffix_part else None
                xml_path = leaf / f"{leaf.name}.xml"
                pv_date = _parse_pv_scan_date(xml_path) if xml_path.exists() else None
                if pv_date is None:
                    continue
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
                continue
            for leaf in sorted(cell_day_dir.iterdir()):
                if not leaf.is_dir():
                    continue
                leaf_match = _RE_SOMATIC_LEAF.match(leaf.name)
                if not leaf_match:
                    continue
                leaf_cell, suffix_part, _seq = leaf_match.groups()
                cell_index = int(leaf_cell)
                suffix = suffix_part.lstrip("_") if suffix_part else None
                xml_path = leaf / f"{leaf.name}.xml"
                pv_date = _parse_pv_scan_date(xml_path) if xml_path.exists() else None
                if pv_date is None:
                    continue
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
                continue
            _yyyy, _mm, _dd, letter = outer.groups()
            cell_index = _LETTER_TO_INDEX.get(letter)
            if cell_index is None:
                continue
            for leaf in sorted(cell_day_dir.iterdir()):
                if not leaf.is_dir():
                    continue
                leaf_match = _RE_SOMATIC_LEAF.match(leaf.name)
                if not leaf_match:
                    continue
                leaf_cell, suffix_part, _seq = leaf_match.groups()
                if int(leaf_cell) != cell_index:
                    continue
                suffix = suffix_part.lstrip("_") if suffix_part else None
                xml_path = leaf / f"{leaf.name}.xml"
                pv_date = _parse_pv_scan_date(xml_path) if xml_path.exists() else None
                if pv_date is None:
                    continue
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
                        return
                    cell_in_name, loc_type, loc_num = parsed
                    xml_path = leaf_dir / f"{leaf_dir.name}.xml"
                    pv_date = _parse_pv_scan_date(xml_path) if xml_path.exists() else None
                    if pv_date is None:
                        return
                    out.append(
                        SweepRecord(
                            modality="dend_excitability",
                            cell_index=cell_in_name,
                            condition_subfolder=condition,
                            recording_path=leaf_dir,
                            pv_scan_date=pv_date,
                            dendrite_location=loc_type,
                            location_number=loc_num,
                        )
                    )

                for grandchild in sorted(cell_day_dir.iterdir()):
                    if grandchild.is_dir():
                        _emit_leaf(grandchild)
                continue  # don't fall through to the LETTER_OUTER branch
            # Match either MMDD<letter> or MMDDYYYY<letter> for the cell-day folder.
            outer = _RE_LETTER_OUTER.match(cell_day_dir.name) or _RE_LETTER_OUTER_FULLYEAR.match(cell_day_dir.name)
            if not outer:
                continue
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
                    return
                cell_in_name, loc_type, loc_num = parsed
                xml_path = leaf_dir / f"{leaf_dir.name}.xml"
                pv_date = _parse_pv_scan_date(xml_path) if xml_path.exists() else None
                if pv_date is None:
                    return
                out.append(
                    SweepRecord(
                        modality="dend_excitability",
                        cell_index=cell_in_name,
                        condition_subfolder=condition,
                        recording_path=leaf_dir,
                        pv_scan_date=pv_date,
                        dendrite_location=loc_type,
                        location_number=loc_num,
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
    for modality_dir in (spec.get("somatic_subdir"), spec.get("dendritic_subdir")):
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
        bundle.cells[sw.cell_index].append(sw)

    bundles = sorted(bundles_by_date.values(), key=lambda b: b.date_compact)
    if unmatched:
        # Print rather than raise: matcher should surface orphans without aborting.
        sample = unmatched[:5]
        print(f"WARNING: {len(unmatched)} sweep(s) did not match a Data Connections subject:")
        for sw in sample:
            print(
                f"  {sw.modality:20s}  cell={sw.cell_index}  cond={sw.condition_subfolder!r}  path={sw.recording_path.name}"
            )
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
