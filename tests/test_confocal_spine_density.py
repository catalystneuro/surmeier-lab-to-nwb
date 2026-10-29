"""Tests for the confocal spine-density conversions. They read the lab's raw data through link_to_raw_data."""

from datetime import datetime
from pathlib import Path
from zoneinfo import ZoneInfo

import pytest

from surmeier_lab_to_nwb.zhai2025.conversion_scripts.confocal_spine_density.figure_4_confocal_spine_density_nikon import (
    convert_session_to_nwbfile,
)
from surmeier_lab_to_nwb.zhai2025.conversion_scripts.confocal_spine_density.figure_4h_confocal_spine_density_olympus import (
    convert_session_to_nwbfile as convert_olympus_session_to_nwbfile,
)

REPO_ROOT = Path(__file__).resolve().parents[1]
FIGURE_4H = REPO_ROOT / "link_to_raw_data/Figure 4_SF1B_SF5/Confocal spine density/Fig 4H/raw"
FIGURE_4J = REPO_ROOT / "link_to_raw_data/Figure 4_SF1B_SF5/Confocal spine density/Fig 4J and Suppl Fig 5/raw"

pytestmark = pytest.mark.skipif(not (REPO_ROOT / "link_to_raw_data").exists(), reason="needs the lab's raw data")


def test_dendrites_of_one_cell_get_their_own_session_ids():
    # Cell 1 of slide 1, slice 2 of mouse 8038 has two imaged dendrites.
    first = convert_session_to_nwbfile(FIGURE_4J / "control/8038-slide 1-slice 2-cell 1 proxi_0001.nd2", "control")
    second = convert_session_to_nwbfile(FIGURE_4J / "control/8038-slide 1-slice 2-cell 1 den2 proxi.nd2", "control")
    assert first.session_id.endswith("_Slide1_Slice2_Cell1_Dendrite1")
    assert second.session_id.endswith("_Slide1_Slice2_Cell1_Dendrite2")


def test_session_id_leaves_out_a_slide_the_file_name_does_not_give():
    nwbfile = convert_session_to_nwbfile(FIGURE_4J / "on-state/5108-slice 2-cell 1 den3 proxi.nd2", "on-state")
    assert nwbfile.session_id == "iSPN++OnState++none++WT++20230901_Slice2_Cell1_Dendrite3"


def test_olympus_session_starts_when_the_mosaic_was_captured():
    # The session date used to be the day of conversion.
    nwbfile = convert_olympus_session_to_nwbfile(FIGURE_4H / "5104-3 60x str_Cycle")
    assert nwbfile.session_start_time == datetime(2024, 9, 28, 12, 53, 12, tzinfo=ZoneInfo("America/Chicago"))
    assert nwbfile.session_id == "iSPN++UnlesionedControl++none++WT++20240928"


def test_control_stacks_are_described_as_unlesioned():
    # Every confocal mouse used to be described as lesioned and levodopa-treated, controls included.
    nwbfile = convert_session_to_nwbfile(FIGURE_4J / "control/8038-slide 1-slice 2-cell 1 proxi_0001.nd2", "control")
    assert "++UnlesionedControl++" in nwbfile.session_id
    assert nwbfile.surgery.startswith("No 6-OHDA lesion")
    assert "levodopa" not in nwbfile.subject.description.replace("no levodopa", "")
