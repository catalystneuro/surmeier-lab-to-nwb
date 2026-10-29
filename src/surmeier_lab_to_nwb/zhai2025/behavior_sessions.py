"""Levodopa sessions of the behaviour mice (Figures 7 and 8, Supplementary Figure 3): one NWB file per mouse and session.

Each mouse was lesioned and then given levodopa every other day, five times. In each session the lab rated its
abnormal involuntary movements (AIMs) every 20 min, in one workbook per figure, and filmed it at the same times to count
contralateral rotations. A session's file holds its AIM scores and its videos.
"""

from __future__ import annotations

import math
import os
import re
import uuid
from dataclasses import dataclass, field
from datetime import date, datetime, time
from pathlib import Path
from zoneinfo import ZoneInfo

import cv2
import pandas as pd
from neuroconv.utils import load_dict_from_file
from pynwb import NWBHDF5IO, NWBFile
from pynwb.epoch import TimeIntervals
from pynwb.file import Subject
from pynwb.image import ImageSeries

from surmeier_lab_to_nwb.zhai2025.subjects import SUBJECTS, SubjectEntry

REPO_ROOT = Path(__file__).resolve().parents[3]
RAW_DATA_ROOT = REPO_ROOT / "link_to_raw_data"
GENERAL_METADATA = Path(__file__).parent / "general_metadata.yaml"
TIMEZONE = ZoneInfo("America/Chicago")

AIM_WORKBOOKS = {
    "Figure 7": RAW_DATA_ROOT / "Figure 7" / "AIM rating" / "AIM testing_CDGI KO.xlsx",
    "Figure 8": RAW_DATA_ROOT / "Figure 8" / "M1R CRISPR AIMs" / "AIM raw score_M1R CRISPR.xlsx",
}
VIDEO_FOLDERS = {
    "Figure 7": RAW_DATA_ROOT / "Figure 7" / "contralateral rotations" / "CDGI KO videos",
    "Figure 8": RAW_DATA_ROOT / "Sup. Figure 3" / "M1R CRISPR videos",
}

_RE_MOUSE = re.compile(r"ET#\s*0*(\d+)")
_RE_MINUTES = re.compile(r"^\s*(\d+)\s*min\s*$")
_RE_WEIGHT = re.compile(r"^\s*([\d.]+)\s*g\s*$")
_RE_COMPACT_DATE = re.compile(r"(?<!\d)(\d{2})(\d{2})(\d{4})(?!\d)")
_RE_SLASH_DATE = re.compile(r"(\d{1,2})/(\d{1,2})(?:/(\d{4}))?")
_RE_ORDINAL = re.compile(r"(\d)[a-z]{2}\b")  # "1st", "2nd", ...
_LABELS = {"Ax": "axial", "Li": "limb", "Ol": "orolingual"}
# Cells the lab used as column headers or left empty, which say nothing about a session
_HEADER_CELLS = {"cage", "cage#", "date:", "(wt)", "(ko)", "g"}

_RE_COHORT_FOLDER = re.compile(r"^(\d{2})(\d{2})(\d{4})-$")
_RE_SESSION_FOLDER = re.compile(r"^(\d)[a-z]{2} session$")  # "1st session", ...
_VIDEO_NAMES = (
    re.compile(r"^0*(?P<mouse>\d+)_(?P<session>\d)_(?P<minutes>\d+)(?:min|imn|mi)$"),  # "imn" and "mi" are typos
    re.compile(r"^(?P<date>\d{8})_ET0*(?P<mouse>\d+)_(?P<minutes>\d+)min$"),
    re.compile(r"^MVI_\d{4}_ET0*(?P<mouse>\d+)_(?P<minutes>\d+)min$"),
)
# Video names that give the wrong mouse. There is no mouse 449 or 548: the 2017-12-15 cohort is 465, 466, 467 and 549,
# the camera filmed 449 in 549's turn of the AIM sheet's order, and 548's one video is the 40 min one 549 lacks.
_MOUSE_NUMBER_TYPOS = {("12152017-", 449): 549, ("12152017-", 548): 549}
# A video named for session 5 in the session 4 folder, which has no other 40 min video of the mouse
_SESSION_NUMBER_TYPOS = {"12152017-/4th session/0465_5_40min.mp4"}

_LEVODOPA_DOSE_MG_PER_KG = {1: 6, 2: 6, 3: 12, 4: 12, 5: 12}
_GROUP = {"CONT": "control", "ko": "CDGI knockout", "M1R CRSPR": "M1R CRISPR"}
_SESSION_ID_GENOTYPE = {"CONT": "WT", "ko": "CDGIKO", "M1R CRSPR": "M1RCRISPR"}
_CRISPR_INJECTION = (
    " About a month before levodopa, {viruses} injected into the dorsolateral striatum on the side of the lesion "
    "(0.5 µl) at six sites, at AP 0.9, ML -2.3, DV -3.4 and -2.8; AP 0.6, ML -1.5, DV -3.4 and -2.8; and AP 0.24, "
    "ML -1.9, DV -3.3 and -2.7 mm relative to Bregma."
)
_CRISPR_VIRUSES = {
    "M1R CRSPR": "AAV9-hSyn-DIO-saCAS-minWPRE and AAV9-U6-gRNA(M1R)_U6-gRNA2-hSyn-FusionRed (1:1), which delete M1R "
    "from the Cre-expressing iSPNs, were",
    "CONT": "the gRNA virus AAV9-U6-gRNA(M1R)_U6-gRNA2-hSyn-FusionRed mixed with saline (control, no Cas9) was",
}


@dataclass
class AimScores:
    """The AIM scores of one mouse in one session, as one block of the lab's workbook."""

    figure: str
    sheet: str
    mouse_number: int
    session_number: int
    date: date
    scores: dict[str, dict[int, float]]  # "axial", "limb", "orolingual" -> minutes after levodopa -> score
    totals: dict[int, float]  # minutes -> the lab's total
    weight_g: float | None = None
    notes: list[str] = field(default_factory=list)  # what the lab wrote next to the cohort, the session or the mouse
    remarks: list[str] = field(default_factory=list)  # where the sheet contradicts itself


@dataclass
class Video:
    """One video of a mouse in a session."""

    path: Path
    mouse_number: int
    session_number: int
    minutes: int  # after the levodopa injection
    named_date: str | None = None  # MMDDYYYY, when the file name has one
    note: str = ""


@dataclass
class BehaviorSession:
    subject: SubjectEntry
    aim: AimScores
    videos: list[Video]


def _number(value) -> float:
    """A score cell: a number, or NaN when empty."""
    if value is None or (isinstance(value, str) and not value.strip()):
        return math.nan
    if isinstance(value, (int, float)):
        return float(value)
    raise ValueError(f"Not a score: {value!r}")


def _session_date(value, sheet_year: int) -> date | None:
    """The date of a session header cell, or None if the cell is not one."""
    if isinstance(value, datetime):
        return value.date()
    if not isinstance(value, str):
        return None
    text = value.strip()
    if not (text.lower().startswith("date") or re.match(r"^\d{1,2}/\d{1,2}/\d{4}", text)):
        return None
    match = _RE_COMPACT_DATE.search(text)
    if match:
        return date(int(match.group(3)), int(match.group(1)), int(match.group(2)))
    match = _RE_SLASH_DATE.search(text)
    if match:  # the 2017-12-15 sheet writes no year
        return date(int(match.group(3) or sheet_year), int(match.group(1)), int(match.group(2)))
    return None


def read_aim_workbook(figure: str) -> list[AimScores]:
    """Read every mouse and session of a figure's AIM workbook; raises on a block it cannot read.

    A sheet is one cohort. A session starts at a row with its date; a row with "20min", "40min", ... gives the columns
    of the times. A mouse is a block of four rows: its number ("ET#2586") with the axial scores, then the limb and the
    orolingual scores, then the lab's totals.
    """
    sessions = []
    workbook = pd.ExcelFile(AIM_WORKBOOKS[figure])
    for sheet in workbook.sheet_names:
        frame = pd.read_excel(workbook, sheet_name=sheet, header=None)
        rows = frame.astype(object).where(frame.notna(), None).values.tolist()
        sheet_year = int(sheet[-4:])
        time_columns, session_number, session_date = {}, 0, None
        sheet_notes, session_notes = [], []
        index = 0
        while index < len(rows):
            row = rows[index]
            times = {col: int(_RE_MINUTES.match(v).group(1)) for col, v in enumerate(row) if _is_minutes(v)}
            if times:
                time_columns = times
            found = _session_date(row[0], sheet_year)
            mouse = _RE_MOUSE.search(row[0]) if isinstance(row[0], str) else None
            if found is not None:
                session_number += 1
                ordinal = _RE_ORDINAL.search(str(row[0]))
                if ordinal and int(ordinal.group(1)) != session_number:
                    raise ValueError(f"{figure} sheet {sheet} row {index}: session {session_number} says {row[0]!r}")
                session_date, session_notes = found, []
            elif mouse is not None:
                if session_date is None:
                    raise ValueError(f"{figure} sheet {sheet} row {index}: mouse before any session date")
                sessions.append(
                    _read_mouse_block(rows, index, time_columns, figure, sheet, session_number, session_date)
                )
                sessions[-1].notes[:0] = sheet_notes + session_notes
                index += 4
                continue
            elif isinstance(row[0], str) and row[0].strip().lower() not in _HEADER_CELLS:
                # A cohort's notes ("cage 1050713", "both WT") before its first session, a session's ("Dose: 6+12")
                # after its date; the 2024-01-17 cohort's "Start date" repeats its first session's date.
                if not row[0].strip().lower().startswith("start date"):
                    (sheet_notes if session_date is None else session_notes).append(row[0].strip())
            index += 1
    return sessions


def _is_minutes(value) -> bool:
    return isinstance(value, str) and _RE_MINUTES.match(value) is not None


def _read_mouse_block(rows, index, time_columns, figure, sheet, session_number, session_date) -> AimScores:
    where = f"{figure} sheet {sheet} row {index}"
    mouse_number = int(_RE_MOUSE.search(rows[index][0]).group(1))
    label_column = next(col for col, v in enumerate(rows[index]) if isinstance(v, str) and v.strip().startswith("Ax"))
    scores, notes, remarks, extra_labels, weight = {}, [], [], [], None
    for offset, (label, name) in enumerate(_LABELS.items()):
        cell = str(rows[index + offset][label_column]).strip()
        if not cell.startswith(label):
            raise ValueError(f"{where}: expected {label} in row {index + offset}, found {cell!r}")
        if cell != label:
            extra_labels.append(cell)
        scores[name] = {minutes: _number(rows[index + offset][col]) for col, minutes in time_columns.items()}
    if extra_labels:
        remarks.append(
            f"The lab wrote {', '.join(repr(cell) for cell in extra_labels)} in the label cells, scores at no time, "
            "which its totals leave out."
        )
    for offset in range(1, 4):
        text = rows[index + offset][0]
        if not isinstance(text, str) or text.strip().lower() in _HEADER_CELLS:
            continue
        weight_match = _RE_WEIGHT.match(text)
        if weight_match:
            weight = float(weight_match.group(1))
        else:
            notes.append(text.strip())

    totals = {}
    for col, minutes in time_columns.items():
        values = [scores[name][minutes] for name in _LABELS.values() if not math.isnan(scores[name][minutes])]
        lab_total = _number(rows[index + 3][col])
        if not values:  # not scored at this time; the lab's sheet sums the empty cells to 0
            if not (math.isnan(lab_total) or lab_total == 0):
                raise ValueError(f"{where}: total {lab_total} at {minutes} min without scores")
            continue
        totals[minutes] = lab_total
        if lab_total != sum(values):
            remarks.append(
                f"At {minutes} min the lab's total ({lab_total:g}) is not the sum of its scores ({sum(values):g})."
            )
    return AimScores(figure, sheet, mouse_number, session_number, session_date, scores, totals, weight, notes, remarks)


def find_videos(figure: str) -> list[Video]:
    """Every video of a figure, by mouse, session and time; raises on a name or folder it cannot read."""
    videos = []
    root = VIDEO_FOLDERS[figure]
    for cohort in sorted(p for p in root.iterdir() if p.is_dir()):
        if not _RE_COHORT_FOLDER.match(cohort.name):
            raise ValueError(f"Unrecognized video folder: {cohort}")
        for session_folder in sorted(p for p in cohort.iterdir() if p.is_dir()):
            folder_match = _RE_SESSION_FOLDER.match(session_folder.name)
            if folder_match is None:
                raise ValueError(f"Unrecognized session folder: {session_folder}")
            session_number = int(folder_match.group(1))
            for path in sorted(p for p in session_folder.iterdir() if p.name != ".DS_Store"):
                match = next((pattern.match(path.stem) for pattern in _VIDEO_NAMES if pattern.match(path.stem)), None)
                if match is None:
                    raise ValueError(f"Unrecognized video name: {path}")
                fields = match.groupdict()
                relative = path.relative_to(root).as_posix()
                named_mouse = int(fields["mouse"])
                mouse_number = _MOUSE_NUMBER_TYPOS.get((cohort.name, named_mouse), named_mouse)
                note = ""
                if mouse_number != named_mouse:
                    note = f"The file name says mouse {named_mouse}, a typo for {mouse_number}."
                named_session = fields.get("session")
                if named_session is not None and int(named_session) != session_number:
                    if relative not in _SESSION_NUMBER_TYPOS:
                        raise ValueError(
                            f"{path} is named for session {named_session} in a session {session_number} folder"
                        )
                    note = f"The file name says session {named_session}; it is in the lab's session {session_number} folder."
                videos.append(
                    Video(path, mouse_number, session_number, int(fields["minutes"]), fields.get("date"), note)
                )
    return videos


def build_behavior_sessions(figure: str) -> list[BehaviorSession]:
    """One session per mouse and levodopa dose, with its AIM scores and its videos; raises if a video has no session."""
    by_number = {subject.mouse_number: subject for subject in SUBJECTS if subject.mouse_number is not None}
    sessions = {}
    for aim in read_aim_workbook(figure):
        key = (aim.mouse_number, aim.session_number)
        if key in sessions:
            raise ValueError(f"{figure}: mouse {aim.mouse_number} has session {aim.session_number} twice")
        if aim.mouse_number not in by_number:
            raise ValueError(f"{figure}: mouse {aim.mouse_number} is not in the subjects table")
        sessions[key] = BehaviorSession(by_number[aim.mouse_number], aim, [])
    for video in find_videos(figure):
        session = sessions.get((video.mouse_number, video.session_number))
        if session is None:
            raise ValueError(f"{video.path} has no AIM session of mouse {video.mouse_number}")
        if video.named_date is not None and video.named_date != session.aim.date.strftime("%m%d%Y"):
            raise ValueError(f"{video.path} is dated {video.named_date}, its session {session.aim.date}")
        if any(other.minutes == video.minutes for other in session.videos):
            raise ValueError(f"{video.path}: two videos at {video.minutes} min")
        session.videos.append(video)
    return [sessions[key] for key in sorted(sessions)]


def build_session_id(session: BehaviorSession) -> str:
    genotype = _SESSION_ID_GENOTYPE[session.subject.category]
    return f"pan++OnState++none++{genotype}++{session.aim.date.strftime('%Y%m%d')}"


def write_behavior_session(session: BehaviorSession, output_path: Path) -> Path:
    """Write one mouse's session: its AIM scores as intervals and its videos as files the NWB file links to, by a path
    relative to it."""
    general = load_dict_from_file(GENERAL_METADATA)
    nwb_meta, subject_meta = general["NWBFile"], general["Subject"]
    subject, aim = session.subject, session.aim
    group = _GROUP[subject.category]
    dose = _LEVODOPA_DOSE_MG_PER_KG[aim.session_number]

    description = (
        f"Levodopa session {aim.session_number} of 5 of mouse {subject.mouse_number} ({group}, {aim.figure}): AIM "
        "(abnormal involuntary movement) scores and the videos filmed to count contralateral rotations. After an "
        f"i.p. injection of levodopa ({dose} mg/kg) with benserazide (12 mg/kg), axial, limb and orolingual "
        "movements were each observed for one minute every 20 min and rated from 0 to 4. Times in this file count "
        "from the levodopa injection, whose clock time the lab did not record; the session start is the date."
    )
    if aim.notes:
        description += " The lab's notes in its AIM sheet: " + ", ".join(f"'{note}'" for note in aim.notes) + "."
    description = " ".join([description, *aim.remarks])
    surgery = nwb_meta["surgery"]
    if aim.figure == "Figure 8":
        surgery += _CRISPR_INJECTION.format(viruses=_CRISPR_VIRUSES[subject.category])
    nwbfile = NWBFile(
        session_description=description,
        identifier=str(uuid.uuid4()),
        session_start_time=datetime.combine(aim.date, time(0), tzinfo=TIMEZONE),
        session_id=build_session_id(session),
        experiment_description=(
            "Levodopa-induced dyskinesia in 6-OHDA-lesioned mice: CDGI knockout and wild-type mice (Figure 7), and "
            "mice with M1R deleted from iSPNs by CRISPR-Cas9 and their controls (Figure 8, Supplementary Figure 3)."
        ),
        experimenter=nwb_meta["experimenter"],
        institution=nwb_meta["institution"],
        lab=nwb_meta["lab"],
        keywords=nwb_meta["keywords"],
        related_publications=nwb_meta["related_publications"],
        surgery=surgery,
        pharmacology=nwb_meta["pharmacology"],
    )
    nwbfile.subject = Subject(
        subject_id=subject.subject_id,
        species=subject_meta["species"],
        sex=subject_meta["sex"],
        age=subject_meta["age"],
        strain=subject.strain,
        genotype=subject.genotype,
        weight=None if aim.weight_g is None else f"{aim.weight_g:g} g",
        description=(
            f"Mouse {subject.mouse_number} in the lab's numbering, {group} group. Unilateral 6-OHDA lesion, then five "
            "levodopa sessions every other day."
        ),
    )

    aim_table = TimeIntervals(
        name="aim_scores",
        description=(
            "AIM scores, one row per one-minute observation, starting at the time after the levodopa injection the "
            "lab's sheet gives. Each score is from 0 to 4; empty where the lab did not score."
        ),
    )
    for name in _LABELS.values():
        aim_table.add_column(name=name, description=f"The {name} AIM score, from 0 to 4.")
    aim_table.add_column(name="total", description="The lab's total of the three scores, as its sheet writes it.")
    for minutes, total in sorted(aim.totals.items()):
        aim_table.add_interval(
            start_time=minutes * 60.0,
            stop_time=minutes * 60.0 + 60.0,
            total=total,
            **{name: aim.scores[name][minutes] for name in _LABELS.values()},
        )
    nwbfile.add_time_intervals(aim_table)

    for video in sorted(session.videos, key=lambda video: video.minutes):
        capture = cv2.VideoCapture(str(video.path))
        rate, frames = capture.get(cv2.CAP_PROP_FPS), int(capture.get(cv2.CAP_PROP_FRAME_COUNT))
        capture.release()
        if not rate or not frames:
            raise ValueError(f"Cannot read the frame rate and frame count of {video.path}")
        nwbfile.add_acquisition(
            ImageSeries(
                name=f"Video{video.minutes}Min",
                description=(
                    f"Video of the mouse from {video.minutes} min after the levodopa injection, filmed to count its "
                    f"contralateral rotations. {video.note}"
                ).strip(),
                external_file=[os.path.relpath(video.path, output_path.parent)],
                format="external",
                starting_frame=[0],
                num_samples=frames,
                rate=rate,
                starting_time=video.minutes * 60.0,
            )
        )

    with NWBHDF5IO(output_path, "w") as io:
        io.write(nwbfile)
    return output_path
