"""Tests for the behaviour sessions. They read the lab's raw data through link_to_raw_data."""

from collections import Counter
from pathlib import Path

import pytest
from pynwb import NWBHDF5IO

from surmeier_lab_to_nwb.zhai2025.behavior_sessions import (
    build_behavior_sessions,
    write_behavior_session,
)

REPO_ROOT = Path(__file__).resolve().parents[1]

pytestmark = pytest.mark.skipif(not (REPO_ROOT / "link_to_raw_data").exists(), reason="needs the lab's raw data")


def test_every_mouse_and_session_of_the_aim_sheets_is_read():
    # The published files kept 9 mice: files were named by day, so mice scored the same day overwrote each other, the
    # 2017-12-15 sheet was skipped and the Figure 8 parser found three mice.
    mice = {}
    for figure, groups in [("Figure 7", {"CONT": 10, "ko": 11}), ("Figure 8", {"CONT": 11, "M1R CRSPR": 10})]:
        sessions = build_behavior_sessions(figure)
        mice[figure] = {session.subject.mouse_number: session.subject.category for session in sessions}
        assert Counter(mice[figure].values()) == groups
        assert Counter(session.subject.mouse_number for session in sessions) == {mouse: 5 for mouse in mice[figure]}
    assert {465, 466, 467, 549, 4123, 4124} <= set(mice["Figure 7"])


def test_every_video_belongs_to_its_mouse_and_session():
    # The published videos were filed under 7 mice; mouse 549's are named 449 and 548 in the share.
    sessions = build_behavior_sessions("Figure 7")
    assert sum(len(session.videos) for session in sessions) == 434
    first = next(s for s in sessions if s.subject.mouse_number == 549 and s.aim.session_number == 1)
    assert sorted(video.minutes for video in first.videos) == [20, 40, 60, 80, 100]
    assert sum(len(session.videos) for session in build_behavior_sessions("Figure 8")) == 464


def test_session_file_has_the_mouse_its_scores_and_videos(tmp_path):
    # Every AIM file said wild type, knockouts included, and the videos all said knockout.
    session = next(
        s for s in build_behavior_sessions("Figure 7") if s.subject.mouse_number == 1944 and s.aim.session_number == 1
    )
    nwbfile_path = write_behavior_session(session, tmp_path / "session.nwb")

    with NWBHDF5IO(nwbfile_path, "r") as io:
        nwbfile = io.read()
        assert nwbfile.subject.subject_id == "Subject20180223koCDGIkoXeGFPM1944"
        assert nwbfile.subject.genotype.startswith("CDGI knockout")
        assert nwbfile.subject.weight == "26.6 g"
        assert nwbfile.session_id == "pan++OnState++none++CDGIKO++20180223"
        scores = nwbfile.intervals["aim_scores"].to_dataframe()
        assert list(scores["start_time"]) == [1200.0, 2400.0, 3600.0, 4800.0]
        assert (scores["total"] == scores[["axial", "limb", "orolingual"]].sum(axis=1)).all()
        assert sorted(nwbfile.acquisition) == ["Video20Min", "Video40Min", "Video60Min", "Video80Min"]
        assert nwbfile.acquisition["Video40Min"].starting_time == 2400.0
