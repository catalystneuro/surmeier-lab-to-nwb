"""Tests for the confocal spine-density conversions. They read the lab's raw data through link_to_raw_data."""

from datetime import datetime
from pathlib import Path
from zoneinfo import ZoneInfo

import pytest

from surmeier_lab_to_nwb.zhai2025.conversion_scripts.confocal_spine_density.figure_4_confocal_spine_density_nikon import (
    convert_mouse_to_nwbfile,
    group_stacks_by_mouse,
)
from surmeier_lab_to_nwb.zhai2025.conversion_scripts.confocal_spine_density.figure_4h_confocal_spine_density_olympus import (
    convert_session_to_nwbfile as convert_olympus_session_to_nwbfile,
)

REPO_ROOT = Path(__file__).resolve().parents[1]
FIGURE_4H = REPO_ROOT / "link_to_raw_data/Figure 4_SF1B_SF5/Confocal spine density/Fig 4H/raw"
FIGURE_4J = REPO_ROOT / "link_to_raw_data/Figure 4_SF1B_SF5/Confocal spine density/Fig 4J and Suppl Fig 5/raw"

pytestmark = pytest.mark.skipif(not (REPO_ROOT / "link_to_raw_data").exists(), reason="needs the lab's raw data")


def test_stacks_grouped_by_mouse():
    # Each stack was its own file; the paper counts 18, 22, 25 and 22 dendrites per group.
    stacks = {condition: group_stacks_by_mouse(FIGURE_4J / condition) for condition in ("control", "6-OHDA")}
    assert {mouse: len(files) for mouse, files in stacks["control"].items()} == {3824: 2, 5104: 4, 8038: 5, 8940: 7}
    assert sum(len(files) for files in stacks["6-OHDA"].values()) == 22


def test_mouse_file_holds_its_stacks_under_the_sheet_mouse():
    # The files said Drd2-EGFP and named the mouse by number only, apart from its Data Connections row.
    nwbfile = convert_mouse_to_nwbfile(group_stacks_by_mouse(FIGURE_4J / "control")[8038], "control")
    assert nwbfile.subject.subject_id == "Subject20240116CONTAdora2CreM8038"
    assert nwbfile.subject.genotype == "Tg(Adora2a-cre)KG139Gsat (RRID:MMRRC_031168-UCD)"
    assert nwbfile.subject.strain is None
    assert nwbfile.session_id == "iSPN++UnlesionedControl++none++WT++20240116"
    # Figure 4J of the preprint is Figure 4K of the published paper
    assert "the data of Figure 4, J and K, and fig. S5" in nwbfile.session_description
    assert nwbfile.surgery.startswith("No 6-OHDA lesion") and "AAV9-pCAG-flex-eGFP" in nwbfile.surgery
    assert "levodopa" not in nwbfile.subject.description.replace("no levodopa", "")
    # Cell 1 of slide 1, slice 2 has two imaged dendrites.
    assert {"ImagesSlide1Slice2Cell1Dendrite1", "ImagesSlide1Slice2Cell1Dendrite2"} <= set(nwbfile.acquisition)
    assert len(nwbfile.acquisition) == 5


def test_stack_name_without_a_slide_leaves_it_out():
    nwbfile = convert_mouse_to_nwbfile(group_stacks_by_mouse(FIGURE_4J / "on-state")[5108], "on-state")
    assert "ImagesSlice2Cell1Dendrite3" in nwbfile.acquisition
    assert nwbfile.session_start_time.date().isoformat() == "2023-09-01"


def test_stack_whose_clocks_disagree_has_no_acquisition_time():
    # Its frame clock says 2024-06-04, months after the data share's date for the file.
    nwbfile = convert_mouse_to_nwbfile(group_stacks_by_mouse(FIGURE_4J / "6-OHDA")[8040], "6-OHDA")
    assert "uncertain" in nwbfile.acquisition["ImagesSlide1Slice3Cell2Dendrite1"].description
    assert "2024-06-04" not in nwbfile.session_description.split("on ")[-1]


def test_olympus_session_starts_when_the_mosaic_was_captured():
    # The session date used to be the day of conversion.
    nwbfile = convert_olympus_session_to_nwbfile(FIGURE_4H / "5104-3 60x str_Cycle")
    assert nwbfile.session_start_time == datetime(2024, 9, 28, 12, 53, 12, tzinfo=ZoneInfo("America/Chicago"))
    assert nwbfile.session_id == "iSPN++UnlesionedControl++none++WT++20240928"
    assert "the image of Figure 4I" in nwbfile.session_description
    assert nwbfile.subject.subject_id == "Subject20230820CONTAdora2CreM5104"
    assert nwbfile.subject.genotype.startswith("Tg(Adora2a-cre)")
