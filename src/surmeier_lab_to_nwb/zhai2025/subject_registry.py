"""Subject registry derived from the lab's Data Connections spreadsheet.

The per-mouse-day merge requires a single source-of-truth mapping from
(date, animal_reporter, category) to a stable subject_id. The spreadsheet
`assets/author_assets/Data Connections.xlsx` provides this: each row in
the D1/D2 sheets corresponds to one mouse, and each row in VIDEOS carries
an explicit lab-assigned MOUSE id alongside the same fields.

This module:

- Parses all four sheets (OVERVIEW, D1, D2, VIDEOS)
- Returns a `(date, animal, category) -> SubjectEntry` registry
- Encodes subject_id as `Subject{YYYYMMDD}{category_slug}{animal_slug}`
- Exposes a `cells_for(date, animal, category, modality)` lookup that returns
  the cell indices the lab patched for a given experiment

See `obsidian_docs/plan_per_cell_merge_with_full_icephys_chain.md` Decision 8
for the design rationale.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[3]
DATA_CONNECTIONS_XLSX = (
    REPO_ROOT / "src" / "surmeier_lab_to_nwb" / "zhai2025" / "assets" / "author_assets" / "Data Connections.xlsx"
)

_MONTH_MAP = {
    "JAN": "01",
    "FEB": "02",
    "MAR": "03",
    "APR": "04",
    "MAY": "05",
    "JUN": "06",
    "JUL": "07",
    "AUG": "08",
    "SEP": "09",
    "OCT": "10",
    "NOV": "11",
    "DEC": "12",
}

# Map column "soma_ex_F1F" -> ("modality", "figure_panel") so the per-experiment
# cell lookup stays explicit. Keys are sheet+section-scoped.
_MODALITY_FIGURE_KEYS = {
    # D1 patch-clamp section (Figures 1, 2)
    "D1_patch.soma_ex_F1F": ("soma_excitability", "F1F"),
    "D1_patch.dend_ex_F1I": ("dend_excitability", "F1I"),
    "D1_patch.spine_dens_F2BC": ("spine_density", "F2BC"),
    "D1_patch.soma_ex_addSCH_F1F": ("soma_excitability_sch", "F1F"),
    "D1_patch.dend_ex_addSCH_F1I": ("dend_excitability_sch", "F1I"),
    "D1_patch.minis_F2H": ("oepsc", "F2H"),
    # D1 biosensor section (Figure 5)
    "D1_biosensor.cells_F5F": ("ach_biosensor", "F5F"),
    # D2 patch-clamp section (Figures 3, 4)
    "D2_patch.soma_ex_F3D": ("soma_excitability", "F3D"),
    "D2_patch.dend_ex_F3F": ("dend_excitability", "F3F"),
    "D2_patch.spine_dens_F4BC": ("spine_density", "F4BC"),
    "D2_patch.soma_ex_addSUL_F3D": ("soma_excitability_sul", "F3D"),
    "D2_patch.dend_ex_addSUL_F3F": ("dend_excitability_sul", "F3F"),
    "D2_patch.minis_F4G": ("oepsc", "F4G"),
    "D2_patch.cf_spine_F4J": ("cf_spine_density", "F4J"),
    # D2 M1R CRISPR / vehicle control section (Figure 6)
    "D2_M1RCRISPR.soma_ex_F6C": ("soma_excitability", "F6C"),
    "D2_M1RCRISPR.dend_ex_F6D": ("dend_excitability", "F6D"),
    "D2_M1RCRISPR.spine_dens_F6F": ("spine_density", "F6F"),
}


@dataclass(frozen=True)
class SubjectEntry:
    """One row of the registry: a single mouse identified by (date, animal, category)."""

    subject_id: str
    year: int
    month: int  # 1-12
    day: int
    animal: str
    category: str
    section: str
    mouse_id_explicit: int | None = None  # only set for VIDEOS sheet rows
    # Maps modality_key -> list of cell indices the lab patched for this experiment.
    # Modality keys are namespaced like "soma_excitability" / "dend_excitability".
    cells_per_experiment: dict[str, tuple[int, ...]] = field(default_factory=dict)

    @property
    def date_iso(self) -> str:
        """YYYY-MM-DD."""
        return f"{self.year:04d}-{self.month:02d}-{self.day:02d}"

    @property
    def date_compact(self) -> str:
        """YYYYMMDD."""
        return f"{self.year:04d}{self.month:02d}{self.day:02d}"


def _slug(value: str) -> str:
    """Normalize a free-text token (category or animal) into a CamelCase-ish slug.

    The lab uses tokens like "LID OFF", "6-OHDA", "M1R CRSPR", "CDGIko X eGFP".
    Strip punctuation/spaces, keep alphanumerics.
    """
    if value is None:
        return ""
    s = str(value).strip()
    s = re.sub(r"[\s/\\-]+", "", s)  # remove whitespace, slashes, dashes
    s = re.sub(r"[^\w]+", "", s)
    return s


def _parse_cells(value) -> tuple[int, ...]:
    """Parse a Data Connections cell-index entry into a tuple of ints.

    Examples of input formats observed in the spreadsheet:
        "1,2,3"         -> (1, 2, 3)
        "2,3"           -> (2, 3)
        "1"             -> (1,)
        "10/5(3)2"      -> (10, 5, 3, 2)   (Fig 7 with slashes/parens; treat all numbers as cell hints)
        "1AB, 2A"       -> (1, 2)          (videos: session labels — extract trailing-digit groups)
        NaN/None/""     -> ()
    """
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return ()
    s = str(value).strip()
    if not s or s.lower() == "nan":
        return ()
    nums = re.findall(r"\d+", s)
    if not nums:
        return ()
    return tuple(int(n) for n in nums)


def _year_int_or_none(value) -> int | None:
    """Coerce a value to a four-digit year, or None if not parseable."""
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return None
    try:
        n = int(str(value).replace(",", "").replace(".0", ""))
    except (ValueError, TypeError):
        return None
    return n if 2000 < n < 2030 else None


def _build_subject_id(date_compact: str, category: str, animal: str) -> str:
    """Build the canonical Subject ID from (date, category, animal)."""
    return f"Subject{date_compact}{_slug(category)}{_slug(animal)}"


def _extract_section_rows(
    df: pd.DataFrame,
    header_row: int,
    columns: dict[str, int],
    extra_cols: dict[str, int],
    section_label: str,
) -> list[SubjectEntry]:
    """Parse one horizontal section of a sheet (one set of YEAR/MONTH/DAY/.../experiment columns).

    `columns` maps "year", "month", "day", "animal", "category" -> their column indices.
    `extra_cols` maps a free-form name -> column index for cell-list columns.
    """
    rows: list[SubjectEntry] = []
    for i in range(header_row + 1, len(df)):
        row = df.iloc[i]
        year = _year_int_or_none(row.iloc[columns["year"]])
        if year is None:
            continue
        month_raw = row.iloc[columns["month"]]
        day_raw = row.iloc[columns["day"]]
        animal = row.iloc[columns["animal"]]
        category = row.iloc[columns["category"]]
        if any(pd.isna(v) for v in (month_raw, day_raw, animal, category)):
            continue
        month_str = str(month_raw).strip().upper()
        month_num = _MONTH_MAP.get(month_str)
        if month_num is None:
            continue
        try:
            day = int(day_raw)
        except (ValueError, TypeError):
            continue
        animal_str = str(animal).strip()
        category_str = str(category).strip()
        date_compact = f"{year:04d}{month_num}{day:02d}"
        subject_id = _build_subject_id(date_compact, category_str, animal_str)
        cells_per_experiment: dict[str, tuple[int, ...]] = {}
        for col_name, col_idx in extra_cols.items():
            full_key = f"{section_label}.{col_name}"
            modality_figure = _MODALITY_FIGURE_KEYS.get(full_key)
            if modality_figure is None:
                continue
            modality_key, _figure_panel = modality_figure
            cells = _parse_cells(row.iloc[col_idx])
            if cells:
                cells_per_experiment[modality_key] = cells
        rows.append(
            SubjectEntry(
                subject_id=subject_id,
                year=year,
                month=int(month_num),
                day=day,
                animal=animal_str,
                category=category_str,
                section=section_label,
                cells_per_experiment=cells_per_experiment,
            )
        )
    return rows


def _load_videos_rows(df: pd.DataFrame) -> list[SubjectEntry]:
    """Parse the VIDEOS sheet, which uses explicit MOUSE IDs."""
    rows: list[SubjectEntry] = []
    for i in range(5, len(df)):
        row = df.iloc[i]
        mouse_id = row.iloc[0]
        if pd.isna(mouse_id):
            continue
        try:
            mid = int(mouse_id)
        except (ValueError, TypeError):
            continue
        year = _year_int_or_none(row.iloc[1])
        if year is None:
            continue
        month_str = str(row.iloc[2]).strip().upper()
        month_num = _MONTH_MAP.get(month_str)
        if month_num is None:
            continue
        try:
            day = int(row.iloc[3])
        except (ValueError, TypeError):
            continue
        animal_str = str(row.iloc[4]).strip()
        category_str = str(row.iloc[5]).strip()
        date_compact = f"{year:04d}{month_num}{day:02d}"
        # For videos, append the lab's explicit MOUSE ID to disambiguate.
        # Multiple mice can share (date, animal, category) in VIDEOS because the
        # lab ran several mice through the AIM protocol per day.
        subject_id = f"Subject{date_compact}{_slug(category_str)}{_slug(animal_str)}M{mid}"
        rows.append(
            SubjectEntry(
                subject_id=subject_id,
                year=year,
                month=int(month_num),
                day=day,
                animal=animal_str,
                category=category_str,
                section="VIDEOS",
                mouse_id_explicit=mid,
                cells_per_experiment={},
            )
        )
    return rows


def load_registry(xlsx_path: Path = DATA_CONNECTIONS_XLSX) -> list[SubjectEntry]:
    """Load the Data Connections spreadsheet and return all subject entries.

    Combines all four sheets. Returns a flat list rather than a dict because some
    sheets (videos) can have multiple subjects for the same (date, animal, category)
    distinguished only by lab-assigned MOUSE id.
    """
    if not xlsx_path.exists():
        raise FileNotFoundError(f"Data Connections spreadsheet not found at {xlsx_path}")

    entries: list[SubjectEntry] = []

    # D1 sheet: header at row 3, two horizontal sections
    df_d1 = pd.read_excel(xlsx_path, sheet_name="D1(FIGs1,2,5)", header=None)
    entries.extend(
        _extract_section_rows(
            df_d1,
            header_row=3,
            columns={"year": 1, "month": 2, "day": 3, "animal": 4, "category": 5},
            extra_cols={
                "soma_ex_F1F": 6,
                "dend_ex_F1I": 7,
                "spine_dens_F2BC": 8,
                "soma_ex_addSCH_F1F": 9,
                "dend_ex_addSCH_F1I": 10,
                "minis_F2H": 11,
            },
            section_label="D1_patch",
        )
    )
    entries.extend(
        _extract_section_rows(
            df_d1,
            header_row=3,
            columns={"year": 13, "month": 14, "day": 15, "animal": 16, "category": 17},
            extra_cols={"cells_F5F": 18},
            section_label="D1_biosensor",
        )
    )

    # D2 sheet: header at row 4, two horizontal sections
    df_d2 = pd.read_excel(xlsx_path, sheet_name="D2(FIGs3,4,6,7,8)", header=None)
    entries.extend(
        _extract_section_rows(
            df_d2,
            header_row=4,
            columns={"year": 1, "month": 2, "day": 3, "animal": 4, "category": 5},
            extra_cols={
                "soma_ex_F3D": 6,
                "dend_ex_F3F": 7,
                "spine_dens_F4BC": 8,
                "soma_ex_addSUL_F3D": 9,
                "dend_ex_addSUL_F3F": 10,
                "minis_F4G": 11,
                "cf_spine_F4J": 12,
            },
            section_label="D2_patch",
        )
    )
    entries.extend(
        _extract_section_rows(
            df_d2,
            header_row=4,
            columns={"year": 14, "month": 15, "day": 16, "animal": 17, "category": 18},
            extra_cols={"soma_ex_F6C": 19, "dend_ex_F6D": 20, "spine_dens_F6F": 21},
            section_label="D2_M1RCRISPR",
        )
    )

    # Videos sheet
    df_v = pd.read_excel(xlsx_path, sheet_name="VIDEOS", header=None)
    entries.extend(_load_videos_rows(df_v))

    return entries


def index_by_key(entries: Iterable[SubjectEntry]) -> dict[tuple[str, str, str], SubjectEntry]:
    """Return a `(date_compact, animal, category) -> SubjectEntry` dict.

    Raises ValueError if any two entries share the same key (excluding VIDEOS
    where the explicit MOUSE id disambiguates).
    """
    out: dict[tuple[str, str, str], SubjectEntry] = {}
    duplicates: list[tuple[tuple[str, str, str], list[SubjectEntry]]] = []
    grouped: dict[tuple[str, str, str], list[SubjectEntry]] = {}
    for e in entries:
        if e.section == "VIDEOS":
            continue
        key = (e.date_compact, e.animal, e.category)
        grouped.setdefault(key, []).append(e)
    for key, group in grouped.items():
        if len(group) > 1:
            duplicates.append((key, group))
        else:
            out[key] = group[0]
    if duplicates:
        msg_lines = [f"Found {len(duplicates)} duplicate (date, animal, category) keys:"]
        for key, group in duplicates[:5]:
            msg_lines.append(f"  {key}: {[e.section for e in group]}")
        raise ValueError("\n".join(msg_lines))
    return out


def summarize(entries: Iterable[SubjectEntry]) -> dict:
    """Return a short summary of the registry: counts per section, category, animal."""
    from collections import Counter

    section_counts = Counter(e.section for e in entries)
    category_counts = Counter(e.category for e in entries)
    animal_counts = Counter(e.animal for e in entries)
    return {
        "total": sum(section_counts.values()),
        "by_section": dict(section_counts),
        "by_category": dict(category_counts),
        "by_animal": dict(animal_counts),
    }


if __name__ == "__main__":
    # Smoke test: load the registry, summarize, print a few entries.
    entries = load_registry()
    summary = summarize(entries)
    print(f"Total subject entries: {summary['total']}")
    print(f"By section: {summary['by_section']}")
    print(f"By category: {summary['by_category']}")
    print(f"By animal: {summary['by_animal']}")

    # Verify the uniqueness invariant for icephys+biosensor sections
    indexed = index_by_key(entries)
    print(f"\nUnique (date, animal, category) keys (excluding VIDEOS): {len(indexed)}")

    # Print one example with cells-per-experiment populated
    examples_with_cells = [e for e in entries if e.cells_per_experiment]
    print(f"\nSample entry with cells_per_experiment:")
    if examples_with_cells:
        e = examples_with_cells[0]
        print(f"  subject_id   : {e.subject_id}")
        print(f"  date         : {e.date_iso}")
        print(f"  animal       : {e.animal}")
        print(f"  category     : {e.category}")
        print(f"  section      : {e.section}")
        print(f"  cells        : {e.cells_per_experiment}")

    # Probe one specific mouse: 2017-02-02 LID OFF eGFP (the Fig 1 example mouse)
    feb2_key = ("20170202", "eGFP", "LID OFF")
    if feb2_key in indexed:
        e = indexed[feb2_key]
        print(f"\nFigure 1 example mouse 2017-02-02 LID OFF eGFP:")
        print(f"  subject_id   : {e.subject_id}")
        print(f"  cells_per_exp: {e.cells_per_experiment}")
    else:
        # Fallback to first LID OFF entry
        lidoff = [e for e in entries if e.category == "LID OFF"]
        if lidoff:
            e = lidoff[0]
            print(f"\nFirst LID OFF mouse: {e.subject_id} ({e.date_iso}), cells: {e.cells_per_experiment}")
