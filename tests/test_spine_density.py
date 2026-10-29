"""Tests for the two-photon spine-density conversions. They read the lab's raw data through link_to_raw_data."""

from pathlib import Path

import pytest

from surmeier_lab_to_nwb.zhai2025.conversion_scripts.spine_density.figure_2_spine_density import (
    convert_data_to_nwb,
)

REPO_ROOT = Path(__file__).resolve().parents[1]
FIGURE_2 = REPO_ROOT / "link_to_raw_data/Figure 2_SF1A/Spine density"

pytestmark = pytest.mark.skipif(not (REPO_ROOT / "link_to_raw_data").exists(), reason="needs the lab's raw data")


def test_control_and_6ohda_stacks_have_their_own_group():
    # Both groups used to be labelled OffState and described as lesioned and levodopa-treated.
    control = convert_data_to_nwb(FIGURE_2 / "control dSPN/0110b", "control dSPN")
    assert control.session_id == "dSPN++UnlesionedControl++none++WT++20170110"
    assert control.surgery.startswith("No 6-OHDA lesion")
    assert "levodopa" not in control.subject.description.replace("no levodopa", "")

    lesioned = convert_data_to_nwb(FIGURE_2 / "PD dSPN/0524a", "PD dSPN")
    assert lesioned.session_id == "dSPN++LesionedControl++none++WT++20170524"
    assert "6-hydroxydopamine" in lesioned.surgery
    assert "L-DOPA" not in lesioned.pharmacology
    assert lesioned.subject.description.count("levodopa") == 1 and "without levodopa" in lesioned.subject.description
