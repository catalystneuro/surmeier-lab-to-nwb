"""Tests for the subjects table."""

from surmeier_lab_to_nwb.zhai2025.subjects import SUBJECTS


def test_every_mouse_has_its_own_subject_id():
    subject_ids = [subject.subject_id for subject in SUBJECTS]
    assert len(subject_ids) == len(set(subject_ids))
