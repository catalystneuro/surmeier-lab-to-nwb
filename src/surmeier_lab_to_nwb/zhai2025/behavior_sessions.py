"""Levodopa sessions of the behaviour mice (Figures 7 and 8, Supplementary Figure 3): one NWB file per mouse and session.

Each mouse was lesioned and then given levodopa every other day, five times. In each session the lab rated its
abnormal involuntary movements (AIMs) every 20 min, in one workbook per figure, and filmed it at the same times to count
contralateral rotations. A session's file holds its AIM scores and its videos.
"""

from __future__ import annotations

import math
import re
import subprocess
import uuid
from dataclasses import dataclass, field
from datetime import date, datetime, time
from pathlib import Path
from zoneinfo import ZoneInfo

import pandas as pd
from hdmf.common import MeaningsTable
from hed.errors import ErrorHandler
from ndx_ethogram import Ethogram
from ndx_hed import HedLabMetaData, HedTags, HedValueVector
from ndx_hed.utils.hed_nwb_validator import HedNWBValidator
from neuroconv.datainterfaces import ExternalVideoInterface
from neuroconv.tools import configure_and_write_nwbfile
from neuroconv.tools.nwb_helpers import make_nwbfile_from_metadata

from surmeier_lab_to_nwb.zhai2025.metadata import session_metadata, shared_metadata
from surmeier_lab_to_nwb.zhai2025.subjects import SUBJECTS, SubjectEntry

REPO_ROOT = Path(__file__).resolve().parents[3]
RAW_DATA_ROOT = REPO_ROOT / "link_to_raw_data"
TIMEZONE = ZoneInfo("America/Chicago")

AIM_WORKBOOKS = {
    "Figure 7": RAW_DATA_ROOT / "Figure 7" / "AIM rating" / "AIM testing_CDGI KO.xlsx",
    "Figure 8": RAW_DATA_ROOT / "Figure 8" / "M1R CRISPR AIMs" / "AIM raw score_M1R CRISPR.xlsx",
}
VIDEO_FOLDERS = {
    "Figure 7": RAW_DATA_ROOT / "Figure 7" / "contralateral rotations" / "CDGI KO videos",
    "Figure 8": RAW_DATA_ROOT / "Sup. Figure 3" / "M1R CRISPR videos",
}
_BEHAVIOR = session_metadata("behavior")
_PUBLISHED_FIGURE = _BEHAVIOR["published_figure"]

_RE_MOUSE = re.compile(r"ET#\s*0*(\d+)")
_RE_MINUTES = re.compile(r"^\s*(\d+)\s*min\s*$")
_RE_WEIGHT = re.compile(r"^\s*([\d.]+)\s*g\s*$")
_RE_COMPACT_DATE = re.compile(r"(?<!\d)(\d{2})(\d{2})(\d{4})(?!\d)")
_RE_SLASH_DATE = re.compile(r"(\d{1,2})/(\d{1,2})(?:/(\d{4}))?")
_RE_ORDINAL = re.compile(r"(\d)[a-z]{2}\b")  # "1st", "2nd", ...
_LABELS = {"Ax": "axial", "Li": "limb", "Ol": "orolingual"}
# The AIM rating scheme, one entry per scored category and the lab's total, and the HED tags of each category
_AIM_CATEGORIES = {name: text.format(scale=_BEHAVIOR["scale"]) for name, text in _BEHAVIOR["categories"].items()}
_HED_SCHEMA_VERSION = _BEHAVIOR["HED"]["schema_version"]
_AIM_HED = _BEHAVIOR["HED"]["categories"]
_AIM_SCORE_HED = _BEHAVIOR["HED"]["score"]
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
# The lab's 80 min video of mouse 1590's session 1 has the same video as its 20 min video of session 2, and shows the
# cage label of session 2: the session 1 label also gives the weight of that day, "BW24.5". So session 1 has no 80 min
# video.
_COPIED_VIDEOS = {"02232018-/1st session/1590_1_80min.mov"}
_MISSING_VIDEO_NOTES = {
    (mouse, session): note
    for mouse, notes in _BEHAVIOR["missing_video_notes"].items()
    for session, note in notes.items()
}

_LEVODOPA_DOSE_MG_PER_KG = _BEHAVIOR["levodopa_dose_mg_per_kg_by_session"]
_GROUP = _BEHAVIOR["group_by_category"]
_SESSION_ID_GENOTYPE = {"CONT": "WT", "ko": "CDGIKO", "M1R CRSPR": "M1RCRISPR"}
_CRISPR_INJECTION = " " + _BEHAVIOR["NWBFile"]["crispr_injection"]
_CRISPR_VIRUSES = _BEHAVIOR["NWBFile"]["crispr_viruses"]


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
            _BEHAVIOR["remarks"]["extra_labels"].format(labels=", ".join(repr(cell) for cell in extra_labels))
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
                _BEHAVIOR["remarks"]["total_not_sum"].format(minutes=minutes, total=lab_total, sum=sum(values))
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
                if relative in _COPIED_VIDEOS:
                    continue
                named_mouse = int(fields["mouse"])
                mouse_number = _MOUSE_NUMBER_TYPOS.get((cohort.name, named_mouse), named_mouse)
                note = ""
                if mouse_number != named_mouse:
                    note = _BEHAVIOR["Video"]["wrong_mouse_note"].format(
                        named_mouse=named_mouse, mouse_number=mouse_number
                    )
                named_session = fields.get("session")
                if named_session is not None and int(named_session) != session_number:
                    if relative not in _SESSION_NUMBER_TYPOS:
                        raise ValueError(
                            f"{path} is named for session {named_session} in a session {session_number} folder"
                        )
                    note = _BEHAVIOR["Video"]["wrong_session_note"].format(
                        named_session=named_session, session_number=session_number
                    )
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


def _frame_count(video_path: Path) -> int:
    """The frame count the video's container header states; raises when it states none."""
    result = subprocess.run(
        ["ffprobe", "-v", "error", "-select_streams", "v:0", "-show_entries", "stream=nb_frames", "-of", "csv=p=0"]
        + [str(video_path)],
        check=True,
        capture_output=True,
        text=True,
    )
    return int(result.stdout)


def make_web_playable_copy(video_path: Path, copy_path: Path) -> Path:
    """Copy a video as an `.mp4` a browser can stream: the same video stream, without re-encoding, with its index at the
    start of the file and without the audio, which may carry the lab room's conversations. A copy that exists is kept;
    raises unless the copy has the frame count of the original."""
    if not copy_path.exists():
        copy_path.parent.mkdir(parents=True, exist_ok=True)
        # Written under another name first, so an interrupted copy is never taken for a finished one
        partial_path = copy_path.with_name(f"partial_{copy_path.name}")
        subprocess.run(
            ["ffmpeg", "-nostdin", "-v", "error", "-y", "-i", str(video_path), "-map", "0:v:0", "-c:v", "copy", "-an"]
            + ["-movflags", "+faststart", str(partial_path)],
            check=True,
        )
        partial_path.rename(copy_path)
    if _frame_count(copy_path) != _frame_count(video_path):
        raise ValueError(f"{copy_path} does not have the frame count of {video_path}")
    return copy_path


def write_behavior_session(session: BehaviorSession, output_path: Path, video_folder: Path) -> Path:
    """Write one mouse's session: its AIM scores as intervals and its videos as files the NWB file links to, by a path
    relative to it. The videos are web-playable copies without audio, made in `video_folder` when missing."""
    metadata = shared_metadata()
    subject, aim = session.subject, session.aim
    group = _GROUP[subject.category]
    dose = _LEVODOPA_DOSE_MG_PER_KG[aim.session_number]
    text = _BEHAVIOR["NWBFile"]

    description = text["session_description"].format(
        session_number=aim.session_number,
        mouse_number=subject.mouse_number,
        group=group,
        figure=_PUBLISHED_FIGURE[aim.figure],
        dose=dose,
    )
    if aim.notes:
        description += " " + text["lab_notes"].format(notes=", ".join(f"'{note}'" for note in aim.notes))
    description = " ".join([description, *aim.remarks])
    if (subject.mouse_number, aim.session_number) in _MISSING_VIDEO_NOTES:
        description += " " + _MISSING_VIDEO_NOTES[(subject.mouse_number, aim.session_number)]
    surgery = metadata["NWBFile"]["surgery"]
    if aim.figure == "Figure 8":
        surgery += _CRISPR_INJECTION.format(viruses=_CRISPR_VIRUSES[subject.category])
    metadata["NWBFile"].update(
        session_description=description,
        identifier=str(uuid.uuid4()),
        session_start_time=datetime.combine(aim.date, time(0), tzinfo=TIMEZONE),
        session_id=build_session_id(session),
        experiment_description=text["experiment_description"],
        surgery=surgery,
    )
    # The behaviour was recorded in the animal, so the slice preparation does not apply
    del metadata["NWBFile"]["slices"]
    metadata["Subject"].update(
        subject_id=subject.subject_id,
        strain=subject.strain,
        genotype=subject.genotype,
        weight=None if aim.weight_g is None else f"{aim.weight_g:g} g",
        description=_BEHAVIOR["Subject"]["description"].format(mouse_number=subject.mouse_number, group=group),
    )
    # The subjects table gives no strain or genotype for some lines, which the file then leaves out
    metadata["Subject"] = {key: value for key, value in metadata["Subject"].items() if value is not None}
    nwbfile = make_nwbfile_from_metadata(metadata)

    # The rating scheme as an ethogram catalogue and each rating as an event, the layout of NeuroConv's manual
    # scoring interfaces: a rating is taken over a one-minute observation, and a 0 is a rating too.
    ethogram_text, events_text = _BEHAVIOR["Ethogram"], _BEHAVIOR["Events"]
    ethogram = Ethogram(name=ethogram_text["name"], description=ethogram_text["description"], exclusive=False)
    for behavior, definition in _AIM_CATEGORIES.items():
        ethogram.add_row(behavior=behavior, definition=definition, behavior_type="point")
    nwbfile.create_processing_module(name="behavior", description=ethogram_text["processing_module_description"]).add(
        ethogram
    )
    aim_scores = nwbfile.create_events_table(
        name=events_text["name"],
        description=events_text["description"],
        source_description=events_text["source_description"],
    )
    aim_scores.add_column(name="behavior", description=events_text["columns"]["behavior"])
    aim_scores.add_column(
        name="score", description=events_text["columns"]["score"], col_cls=HedValueVector, hed=_AIM_SCORE_HED
    )
    for minutes, total in sorted(aim.totals.items()):
        ratings = {name: aim.scores[name][minutes] for name in _LABELS.values()}
        for behavior, score in [*ratings.items(), ("total", total)]:
            if not math.isnan(score):
                aim_scores.add_row(timestamp=minutes * 60.0, duration=60.0, behavior=behavior, score=score)
    behavior_meanings = MeaningsTable(target=aim_scores["behavior"], description=events_text["meanings_description"])
    behavior_meanings.add_column(name="HED", description=events_text["hed_column_description"], col_cls=HedTags)
    for behavior, definition in _AIM_CATEGORIES.items():
        behavior_meanings.add_row(value=behavior, meaning=definition, HED=_AIM_HED[behavior])
    aim_scores.add_meanings_table(behavior_meanings)
    hed_metadata = HedLabMetaData(hed_schema_version=_HED_SCHEMA_VERSION)
    nwbfile.add_lab_meta_data(hed_metadata)
    issues = HedNWBValidator(hed_metadata).validate_file(nwbfile, ErrorHandler(check_for_warnings=True))
    if issues:
        raise ValueError(f"The AIM HED tags do not validate against HED {_HED_SCHEMA_VERSION}: {issues}")

    # One video interface per video, placed at its minute after the injection. The lab's files name no camera, so the
    # videos link to no device.
    for video in sorted(session.videos, key=lambda video: video.minutes):
        relative_path = video.path.relative_to(VIDEO_FOLDERS[aim.figure])
        copy_path = make_web_playable_copy(video.path, video_folder / relative_path.with_suffix(".mp4"))
        interface = ExternalVideoInterface(file_paths=[copy_path], metadata_key=f"video_{video.minutes}_min")
        if not all(interface.get_header_frame_rates()) or not all(interface.get_header_frame_counts()):
            raise ValueError(f"Cannot read the frame rate and frame count of {copy_path}")
        (video_key,) = interface.temporal_alignment.keys()
        interface.temporal_alignment[video_key].move_start_to(video.minutes * 60.0)
        metadata = interface.get_metadata()
        video_metadata = metadata["Behavior"]["ExternalVideos"][interface.metadata_key]
        del video_metadata["device_metadata_key"]
        video_metadata.update(
            name=_BEHAVIOR["Video"]["name"].format(minutes=video.minutes),
            description=_BEHAVIOR["Video"]["description"].format(minutes=video.minutes, note=video.note).strip(),
        )
        interface.add_to_nwbfile(nwbfile, metadata=metadata)

    # Writes each video's path relative to the NWB file
    configure_and_write_nwbfile(nwbfile, nwbfile_path=output_path, metadata=shared_metadata())
    return output_path
