"""Tests for the two-photon spine-density conversions. They read the lab's raw data through link_to_raw_data."""

from pathlib import Path

import pytest

from surmeier_lab_to_nwb.zhai2025.conversion_scripts.spine_density import (
    figure_8_spine_density,
)
from surmeier_lab_to_nwb.zhai2025.conversion_scripts.spine_density.figure_2_spine_density import (
    convert_data_to_nwb,
)
from surmeier_lab_to_nwb.zhai2025.conversion_scripts.spine_density.spine_density_utils import (
    cell_date,
    group_cells_by_mouse_day,
)

REPO_ROOT = Path(__file__).resolve().parents[1]
FIGURE_2 = REPO_ROOT / "link_to_raw_data/Figure 2_SF1A/Spine density"
FIGURE_4 = REPO_ROOT / "link_to_raw_data/Figure 4_SF1B_SF5/Spine density"
FIGURE_8 = REPO_ROOT / "link_to_raw_data/Figure 8/M1R CRISPR spine density"

pytestmark = pytest.mark.skipif(not (REPO_ROOT / "link_to_raw_data").exists(), reason="needs the lab's raw data")


def test_control_and_6ohda_stacks_have_their_own_group():
    # Both groups used to be labelled OffState and described as lesioned and levodopa-treated.
    control = convert_data_to_nwb([FIGURE_2 / "control dSPN/0110b"], "control dSPN")
    assert control.session_id == "dSPN++UnlesionedControl++none++WT++20170110"
    assert control.surgery.startswith("No 6-OHDA lesion")
    assert "levodopa" not in control.subject.description.replace("no levodopa", "")

    lesioned = convert_data_to_nwb([FIGURE_2 / "PD dSPN/0524a", FIGURE_2 / "PD dSPN/0524b"], "PD dSPN")
    assert lesioned.session_id == "dSPN++LesionedControl++none++WT++20170524"
    assert "6-hydroxydopamine" in lesioned.surgery
    assert "L-DOPA" not in lesioned.pharmacology
    assert lesioned.subject.description.count("levodopa") == 1 and "without levodopa" in lesioned.subject.description


def test_every_cell_of_a_mouse_day_is_in_its_file():
    # The cells of one day used to overwrite each other: Figure 2 has 46 cell folders, and only 28 files were left.
    conditions = [condition for condition in FIGURE_2.iterdir() if condition.is_dir()]
    mouse_days = [cells for condition in conditions for cells in group_cells_by_mouse_day(condition).values()]
    assert len(mouse_days) == 28 and sum(len(cells) for cells in mouse_days) == 46

    # 2016-12-16 has cells b, c and d; folder b's stacks say Cell1, but the folder and the lab's sheet say cell 2.
    nwbfile = convert_data_to_nwb(sorted((FIGURE_2 / "PD dSPN").glob("1216?")), "PD dSPN")
    assert nwbfile.subject.subject_id == "Subject201612166OHDAtDTomato"
    stacks = sorted(name for name in nwbfile.acquisition if name.startswith("Images"))
    assert {name[len("ImagesCell")] for name in stacks} == {"2", "3", "4"}
    # Cell 3 has dist2 and dist2_real, two spots of one site.
    assert {"ImagesCell3Distal2Spot1", "ImagesCell3Distal2Spot2"} <= set(stacks)


def test_cell_day_comes_from_the_folder_when_file_names_disagree():
    # The image names of this cell say 2014-04-26; the folder says 2019-04-26.
    assert cell_date(FIGURE_4 / "LID on-state iSPN/04262019a") == "20190426"
    # The image names of these cells say 06-02; the folder and the stack names say 05-02.
    assert cell_date(FIGURE_4 / "PD iSPN/0502b") == "20160502"


def test_crispr_mice_have_the_cre_line_as_genotype():
    # The CRISPR deletion was made by virus; the files called it the genotype and the background C57BL/6, which the
    # paper does not give for these mice.
    nwbfile = figure_8_spine_density.convert_data_to_nwb([FIGURE_8 / "M1R CRISPR/20221005a"], "M1R CRISPR")
    assert nwbfile.subject.subject_id == "Subject20221005M1RCRSPRAdora2Cre"
    assert nwbfile.subject.genotype == "Tg(Adora2a-cre)KG139Gsat (RRID:MMRRC_031168-UCD)"
    assert nwbfile.subject.strain is None


def test_descriptions_have_no_literal_quote_marks():
    # The YAML wrapped the default descriptions in quote marks inside folded blocks, and the files kept them.
    nwbfile = convert_data_to_nwb([FIGURE_2 / "LID off-state dSPN/0202a"], "LID off-state dSPN")
    assert nwbfile.subject.description.startswith("Adult Drd1-Tdtomato transgenic mouse")
    assert nwbfile.session_description.startswith("Dendritic spine density assessment")
    assert '"' not in nwbfile.subject.description + nwbfile.session_description


def test_descriptions_end_without_a_newline():
    # The folded YAML blocks kept their final line break, so every description ended with a newline.
    lid = convert_data_to_nwb([FIGURE_2 / "LID off-state dSPN/0202a"], "LID off-state dSPN")
    control = convert_data_to_nwb([FIGURE_2 / "control dSPN/0110b"], "control dSPN")
    for text in (lid.session_description, lid.subject.description, control.subject.description):
        assert text == text.strip()
