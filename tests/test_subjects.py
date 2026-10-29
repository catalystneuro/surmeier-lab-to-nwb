"""Tests for the subjects table."""

from surmeier_lab_to_nwb.zhai2025.subjects import SUBJECTS


def test_every_mouse_has_its_own_subject_id():
    subject_ids = [subject.subject_id for subject in SUBJECTS]
    assert len(subject_ids) == len(set(subject_ids))


def test_every_numbered_mouse_has_its_own_number():
    numbers = [subject.mouse_number for subject in SUBJECTS if subject.mouse_number is not None]
    assert len(numbers) == len(set(numbers)) == 58
