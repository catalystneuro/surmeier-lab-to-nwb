"""Every mouse of the dataset, one row each: the day and the line and group the lab lists it under, its strain and its
genotype.

The rows come from the lab's Data Connections spreadsheet (`assets/author_assets/Data Connections.xlsx`), one per mouse
of its D1, D2 and VIDEOS sheets, and are edited here. `note` says where a row departs from the spreadsheet and why.
Two kinds of rows were left out: sheet rows with no recordings in the data share (2016-07-17 Drd2-eGFP LID ON, and six
Adora2-Cre rows of September 2021 that list no cells).

Strain is the background colony and genotype the modified loci relative to it, with the line stock each allele came
from. The paper (Methods, Animals) says the Drd1-tdTomato and Drd2-eGFP lines were backcrossed to C57BL/6 and the CDGI
knockout is on C57BL/6J; it gives no background for Adora2-Cre. The stock RRIDs are those of the lab's key resources
table; the registered stocks are on other backgrounds (FVB for Drd1-tdTomato, mixed for Drd2-eGFP and Adora2-Cre).
"""

from __future__ import annotations

import re
from dataclasses import dataclass

C57BL6 = "C57BL/6"
C57BL6J = "C57BL/6J"

DRD1_TDTOMATO = "Tg(Drd1-tdTomato)5Calak, hemizygous (RRID:MMRRC_030512-UNC)"
DRD2_EGFP = "Tg(Drd2-EGFP)S118Gsat, hemizygous (RRID:MMRRC_000230-UNC)"
CDGI_KO = "CDGI knockout"
CDGI_KO_DRD2_EGFP = f"{CDGI_KO}; {DRD2_EGFP}"
# The paper does not say how many copies of the Cre transgene these mice carry.
ADORA2A_CRE = "Tg(Adora2a-cre)KG139Gsat (RRID:MMRRC_031168-UCD)"
WILD_TYPE = "wild type"


@dataclass(frozen=True)
class SubjectEntry:
    """One mouse."""

    date_compact: str  # YYYYMMDD: recording day, a behaviour mouse's first session, or a confocal mouse's imaging day
    animal: str | None  # the line, as Data Connections writes it ("tD-Tomato", "eGFP", "CDGIko X eGFP", ...)
    category: str  # the group, as Data Connections writes it ("CONT", "LID OFF", "M1R CRSPR", ...)
    strain: str | None
    genotype: str
    mouse_number: int | None = None  # the lab's number, for the behaviour and confocal mice
    note: str = ""

    @property
    def subject_id(self) -> str:
        """Subject, the date, the group and the line, and for a numbered mouse M and its number."""
        subject_id = f"Subject{self.date_compact}{_slug(self.category)}{_slug(self.animal)}"
        return subject_id if self.mouse_number is None else f"{subject_id}M{self.mouse_number}"


def _slug(value: str | None) -> str:
    """Drop spaces, slashes, dashes and other punctuation: "CDGIko X eGFP" -> "CDGIkoXeGFP"."""
    return re.sub(r"\W+", "", value or "")


def find_subject(date_compact: str, animal: str, category: str) -> SubjectEntry:
    """The mouse the lab recorded on a day in a line and group; raises unless the table has exactly one."""
    matches = [
        subject
        for subject in SUBJECTS
        if (subject.date_compact, subject.animal, subject.category) == (date_compact, animal, category)
    ]
    if len(matches) != 1:
        raise ValueError(f"The subjects table has {len(matches)} mice for {date_compact}, {animal}, {category}")
    return matches[0]


SUBJECTS = (
    # Drd1-tdTomato mice: dSPN patch clamp, spines and Sr2+-oEPSCs (Figures 1 and 2)
    SubjectEntry("20160714", "tD-Tomato", "CONT", C57BL6, DRD1_TDTOMATO),
    SubjectEntry("20160727", "tD-Tomato", "CONT", C57BL6, DRD1_TDTOMATO),
    SubjectEntry("20160810", "tD-Tomato", "CONT", C57BL6, DRD1_TDTOMATO),
    SubjectEntry("20160812", "tD-Tomato", "LID OFF", C57BL6, DRD1_TDTOMATO),
    SubjectEntry("20161101", "tD-Tomato", "6-OHDA", C57BL6, DRD1_TDTOMATO),
    SubjectEntry("20161111", "tD-Tomato", "LID OFF", C57BL6, DRD1_TDTOMATO),
    SubjectEntry("20161114", "tD-Tomato", "LID OFF", C57BL6, DRD1_TDTOMATO),
    SubjectEntry("20161116", "tD-Tomato", "LID OFF", C57BL6, DRD1_TDTOMATO),
    SubjectEntry("20161117", "tD-Tomato", "6-OHDA", C57BL6, DRD1_TDTOMATO),
    SubjectEntry("20161121", "tD-Tomato", "6-OHDA", C57BL6, DRD1_TDTOMATO),
    SubjectEntry("20161129", "tD-Tomato", "CONT", C57BL6, DRD1_TDTOMATO),
    SubjectEntry("20161215", "tD-Tomato", "CONT", C57BL6, DRD1_TDTOMATO),
    SubjectEntry("20161216", "tD-Tomato", "6-OHDA", C57BL6, DRD1_TDTOMATO),
    SubjectEntry("20170110", "tD-Tomato", "CONT", C57BL6, DRD1_TDTOMATO),
    SubjectEntry("20170202", "tD-Tomato", "LID OFF", C57BL6, DRD1_TDTOMATO),
    SubjectEntry("20170203", "tD-Tomato", "LID OFF", C57BL6, DRD1_TDTOMATO),
    SubjectEntry("20170204", "tD-Tomato", "LID OFF", C57BL6, DRD1_TDTOMATO),
    SubjectEntry("20170524", "tD-Tomato", "6-OHDA", C57BL6, DRD1_TDTOMATO),
    SubjectEntry("20170525", "tD-Tomato", "LID OFF", C57BL6, DRD1_TDTOMATO),
    SubjectEntry("20170706", "tD-Tomato", "LID ON", C57BL6, DRD1_TDTOMATO),
    SubjectEntry("20170707", "tD-Tomato", "LID ON", C57BL6, DRD1_TDTOMATO),
    SubjectEntry("20170719", "tD-Tomato", "LID ON", C57BL6, DRD1_TDTOMATO),
    SubjectEntry("20170720", "tD-Tomato", "LID ON", C57BL6, DRD1_TDTOMATO),
    SubjectEntry("20170721", "tD-Tomato", "LID ON", C57BL6, DRD1_TDTOMATO),
    SubjectEntry("20190411", "tD-Tomato", "LID ON", C57BL6, DRD1_TDTOMATO),
    SubjectEntry("20190412", "tD-Tomato", "LID ON", C57BL6, DRD1_TDTOMATO),
    SubjectEntry("20190501", "tD-Tomato", "LID ON", C57BL6, DRD1_TDTOMATO),
    SubjectEntry("20190503", "tD-Tomato", "LID ON", C57BL6, DRD1_TDTOMATO),
    SubjectEntry("20190507", "tD-Tomato", "LID ON", C57BL6, DRD1_TDTOMATO),
    SubjectEntry("20230705", "tD-Tomato", "LID ON", C57BL6, DRD1_TDTOMATO),
    SubjectEntry("20230706", "tD-Tomato", "LID OFF", C57BL6, DRD1_TDTOMATO),
    SubjectEntry("20230710", "tD-Tomato", "LID ON", C57BL6, DRD1_TDTOMATO),
    SubjectEntry("20230713", "tD-Tomato", "LID OFF", C57BL6, DRD1_TDTOMATO),
    SubjectEntry("20230802", "tD-Tomato", "LID ON", C57BL6, DRD1_TDTOMATO),
    SubjectEntry("20230816", "tD-Tomato", "LID ON", C57BL6, DRD1_TDTOMATO),
    SubjectEntry("20230817", "tD-Tomato", "LID OFF", C57BL6, DRD1_TDTOMATO),
    # C57BL/6J mice: GRAB-ACh imaging (Figure 5)
    SubjectEntry("20240402", "C57bl/6J", "CONT", C57BL6J, WILD_TYPE),
    SubjectEntry("20240405", "C57bl/6J", "CONT", C57BL6J, WILD_TYPE),
    SubjectEntry("20240410", "C57bl/6J", "CONT", C57BL6J, WILD_TYPE),
    SubjectEntry("20240412", "C57bl/6J", "6-OHDA", C57BL6J, WILD_TYPE),
    SubjectEntry("20240416", "C57bl/6J", "6-OHDA", C57BL6J, WILD_TYPE),
    SubjectEntry("20240509", "C57bl/6J", "6-OHDA", C57BL6J, WILD_TYPE),
    SubjectEntry("20240524", "C57bl/6J", "LID OFF", C57BL6J, WILD_TYPE),
    SubjectEntry("20240529", "C57bl/6J", "LID OFF", C57BL6J, WILD_TYPE),
    SubjectEntry("20240531", "C57bl/6J", "LID OFF", C57BL6J, WILD_TYPE),
    SubjectEntry("20240603", "C57bl/6J", "LID OFF", C57BL6J, WILD_TYPE),
    SubjectEntry("20240604", "C57bl/6J", "LID OFF", C57BL6J, WILD_TYPE),
    # Drd2-eGFP mice: iSPN patch clamp, spines and Sr2+-oEPSCs (Figures 3, 4 and 6), and the Figure 7E wild types
    SubjectEntry("20151119", "eGFP", "CONT", C57BL6, DRD2_EGFP),
    SubjectEntry("20151123", "eGFP", "CONT", C57BL6, DRD2_EGFP),
    SubjectEntry("20151124", "eGFP", "6-OHDA", C57BL6, DRD2_EGFP),
    SubjectEntry("20151202", "eGFP", "6-OHDA", C57BL6, DRD2_EGFP),
    SubjectEntry("20151216", "eGFP", "CONT", C57BL6, DRD2_EGFP),
    SubjectEntry("20151217", "eGFP", "CONT", C57BL6, DRD2_EGFP),
    SubjectEntry("20151222", "eGFP", "6-OHDA", C57BL6, DRD2_EGFP),
    SubjectEntry("20151222", "eGFP", "CONT", C57BL6, DRD2_EGFP),
    SubjectEntry("20160108", "eGFP", "CONT", C57BL6, DRD2_EGFP),
    SubjectEntry("20160215", "eGFP", "6-OHDA", C57BL6, DRD2_EGFP),
    SubjectEntry("20160216", "eGFP", "CONT", C57BL6, DRD2_EGFP),
    SubjectEntry("20160218", "eGFP", "6-OHDA", C57BL6, DRD2_EGFP),
    SubjectEntry("20160502", "eGFP", "6-OHDA", C57BL6, DRD2_EGFP),
    SubjectEntry("20160519", "eGFP", "CONT", C57BL6, DRD2_EGFP),
    SubjectEntry("20160523", "eGFP", "LID OFF", C57BL6, DRD2_EGFP),
    SubjectEntry("20160524", "eGFP", "LID OFF", C57BL6, DRD2_EGFP),
    SubjectEntry("20160608", "eGFP", "LID ON", C57BL6, DRD2_EGFP),
    SubjectEntry("20160609", "eGFP", "LID ON", C57BL6, DRD2_EGFP),
    SubjectEntry("20160718", "eGFP", "LID ON", C57BL6, DRD2_EGFP),
    SubjectEntry("20160719", "eGFP", "LID ON", C57BL6, DRD2_EGFP),
    SubjectEntry("20161116", "eGFP", "LID OFF", C57BL6, DRD2_EGFP),
    SubjectEntry("20161117", "eGFP", "CONT", C57BL6, DRD2_EGFP),
    SubjectEntry("20170525", "eGFP", "LID OFF", C57BL6, DRD2_EGFP),
    SubjectEntry("20170526", "eGFP", "LID OFF", C57BL6, DRD2_EGFP),
    SubjectEntry(
        "20180221",
        "eGFP",
        "NAIVE",
        C57BL6,
        DRD2_EGFP,
        note="Data Connections lists this mouse as LID ON: WT; the paper calls the Figure 7E mice naive.",
    ),
    SubjectEntry(
        "20180305",
        "eGFP",
        "NAIVE",
        C57BL6,
        DRD2_EGFP,
        note="Data Connections lists this mouse as LID ON: WT; the paper calls the Figure 7E mice naive.",
    ),
    SubjectEntry("20190408", "eGFP", "LID ON", C57BL6, DRD2_EGFP),
    SubjectEntry("20190409", "eGFP", "LID ON", C57BL6, DRD2_EGFP),
    SubjectEntry("20190426", "eGFP", "LID ON", C57BL6, DRD2_EGFP),
    SubjectEntry("20190429", "eGFP", "LID ON", C57BL6, DRD2_EGFP),
    SubjectEntry("20190508", "eGFP", "LID ON", C57BL6, DRD2_EGFP),
    SubjectEntry("20190509", "eGFP", "LID ON", C57BL6, DRD2_EGFP),
    SubjectEntry("20200217", "eGFP", "CONT", C57BL6, DRD2_EGFP),
    SubjectEntry("20200218", "eGFP", "M1R ant", C57BL6, DRD2_EGFP),
    SubjectEntry("20200220", "eGFP", "M1R ant", C57BL6, DRD2_EGFP),
    SubjectEntry("20200225", "eGFP", "CONT", C57BL6, DRD2_EGFP),
    SubjectEntry("20200312", "eGFP", "M1R ant", C57BL6, DRD2_EGFP),
    SubjectEntry("20200721", "eGFP", "M1R ant", C57BL6, DRD2_EGFP),
    SubjectEntry("20200722", "eGFP", "M1R ant", C57BL6, DRD2_EGFP),
    SubjectEntry("20200811", "eGFP", "CONT", C57BL6, DRD2_EGFP),
    SubjectEntry("20200812", "eGFP", "CONT", C57BL6, DRD2_EGFP),
    SubjectEntry("20230222", "eGFP", "LID OFF", C57BL6, DRD2_EGFP),
    SubjectEntry("20230228", "eGFP", "LID ON", C57BL6, DRD2_EGFP),
    SubjectEntry("20230301", "eGFP", "LID ON", C57BL6, DRD2_EGFP),
    SubjectEntry("20230330", "eGFP", "LID OFF", C57BL6, DRD2_EGFP),
    SubjectEntry("20230331", "eGFP", "LID ON", C57BL6, DRD2_EGFP),
    SubjectEntry("20230405", "eGFP", "LID ON", C57BL6, DRD2_EGFP),
    SubjectEntry("20230406", "eGFP", "LID OFF", C57BL6, DRD2_EGFP),
    SubjectEntry("20240425", "eGFP", "CONT", C57BL6, DRD2_EGFP),
    SubjectEntry(
        "20240430",
        "eGFP",
        "M1R ant",
        C57BL6,
        DRD2_EGFP,
        note="Data Connections dates this mouse 2024-04-03; the microscope dates its recordings 2024-04-30, and the dendritic cells the sheet lists match.",
    ),
    # CDGI knockout mice crossed with Drd2-eGFP (Figure 7)
    SubjectEntry(
        "20160525",
        "CDGIko X eGFP",
        "LID OFF",
        C57BL6,
        CDGI_KO_DRD2_EGFP,
        note="Data Connections lists this mouse as LID ON; its recordings are in the off-state folder, as the paper's Figure 7C counts require.",
    ),
    SubjectEntry(
        "20160526",
        "CDGIko X eGFP",
        "LID OFF",
        C57BL6,
        CDGI_KO_DRD2_EGFP,
        note="Data Connections lists this mouse as LID ON; its recordings are in the off-state folder, as the paper's Figure 7C counts require.",
    ),
    SubjectEntry(
        "20160527",
        "CDGIko X eGFP",
        "LID OFF",
        C57BL6,
        CDGI_KO_DRD2_EGFP,
        note="Data Connections lists this mouse as LID ON; its recordings are in the off-state folder, as the paper's Figure 7C counts require.",
    ),
    SubjectEntry(
        "20160816",
        "CDGIko X eGFP",
        "LID ON",
        C57BL6,
        CDGI_KO_DRD2_EGFP,
        note="Data Connections lists this mouse as LID OFF; its recordings are in the on-state folder, as the paper's Figure 7C counts require.",
    ),
    SubjectEntry(
        "20160901",
        "CDGIko X eGFP",
        "LID ON",
        C57BL6,
        CDGI_KO_DRD2_EGFP,
        note="Data Connections lists this mouse as LID OFF; its recordings are in the on-state folder, as the paper's Figure 7C counts require.",
    ),
    SubjectEntry(
        "20160902",
        "CDGIko X eGFP",
        "LID ON",
        C57BL6,
        CDGI_KO_DRD2_EGFP,
        note="Data Connections lists this mouse as LID OFF; its recordings are in the on-state folder, as the paper's Figure 7C counts require.",
    ),
    SubjectEntry(
        "20160912",
        "CDGIko X eGFP",
        "LID OFF",
        C57BL6,
        CDGI_KO_DRD2_EGFP,
        note="Data Connections lists this mouse as LID ON; its spine images are in the off-state folder.",
    ),
    SubjectEntry(
        "20160913",
        "CDGIko X eGFP",
        "LID ON",
        C57BL6,
        CDGI_KO_DRD2_EGFP,
        note="Data Connections lists this mouse as LID OFF; its recordings are in the on-state folder, as the paper's Figure 7C counts require.",
    ),
    SubjectEntry(
        "20170616",
        "CDGIko X eGFP",
        "LID OFF",
        C57BL6,
        CDGI_KO_DRD2_EGFP,
        note="Data Connections lists this mouse on 2016-06-16 as LID ON; the microscope dates its recordings 2017-06-16, its recordings are in the off-state folder, and the dendritic and spine cells the sheet lists match.",
    ),
    SubjectEntry(
        "20170619",
        "CDGIko X eGFP",
        "LID OFF",
        C57BL6,
        CDGI_KO_DRD2_EGFP,
        note="Data Connections lists this mouse on 2016-06-19 as LID ON; the microscope dates its recordings 2017-06-19, its recordings are in the off-state folder, and the dendritic and spine cells the sheet lists match.",
    ),
    SubjectEntry(
        "20180228",
        "CDGIko X eGFP",
        "NAIVE",
        C57BL6,
        CDGI_KO_DRD2_EGFP,
        note="Data Connections lists this mouse as LID ON; the paper calls the Figure 7E mice naive.",
    ),
    SubjectEntry(
        "20180301",
        "CDGIko X eGFP",
        "NAIVE",
        C57BL6,
        CDGI_KO_DRD2_EGFP,
        note="Data Connections lists this mouse as LID ON; the paper calls the Figure 7E mice naive.",
    ),
    SubjectEntry(
        "20180314",
        "CDGIko X eGFP",
        "NAIVE",
        C57BL6,
        CDGI_KO_DRD2_EGFP,
        note="Data Connections lists this mouse as LID ON; the paper calls the Figure 7E mice naive.",
    ),
    # Adora2-Cre mice: M1R CRISPR patch clamp and spines (Figure 8) and high-resolution confocal stacks (Figure 4)
    SubjectEntry("20221004", "Adora2-Cre", "M1R CRSPR", None, ADORA2A_CRE),
    SubjectEntry("20221005", "Adora2-Cre", "M1R CRSPR", None, ADORA2A_CRE),
    SubjectEntry("20221012", "Adora2-Cre", "CONT", None, ADORA2A_CRE),
    SubjectEntry("20221129", "Adora2-Cre", "M1R CRSPR", None, ADORA2A_CRE),
    SubjectEntry("20221130", "Adora2-Cre", "CONT", None, ADORA2A_CRE),
    # The confocal mice of Figure 4J. The sheet lists them by group, number of cells and a day they were imaged;
    # the number is the one the lab gave their image files, which match the rows one to one on those three.
    SubjectEntry("20230820", "Adora2-Cre", "CONT", None, ADORA2A_CRE, mouse_number=5104),
    SubjectEntry("20230820", "Adora2-Cre", "LID OFF", None, ADORA2A_CRE, mouse_number=5105),
    SubjectEntry("20230901", "Adora2-Cre", "LID ON", None, ADORA2A_CRE, mouse_number=5108),
    SubjectEntry("20230911", "Adora2-Cre", "LID ON", None, ADORA2A_CRE, mouse_number=5107),
    SubjectEntry("20231022", "Adora2-Cre", "6-OHDA", None, ADORA2A_CRE, mouse_number=3825),
    SubjectEntry("20231022", "Adora2-Cre", "CONT", None, ADORA2A_CRE, mouse_number=3824),
    SubjectEntry("20231026", "Adora2-Cre", "6-OHDA", None, ADORA2A_CRE, mouse_number=3823),
    SubjectEntry("20231029", "Adora2-Cre", "LID OFF", None, ADORA2A_CRE, mouse_number=3826),
    SubjectEntry("20231211", "Adora2-Cre", "6-OHDA", None, ADORA2A_CRE, mouse_number=8040),
    SubjectEntry("20231217", "Adora2-Cre", "LID OFF", None, ADORA2A_CRE, mouse_number=8041),
    SubjectEntry("20240112", "Adora2-Cre", "LID ON", None, ADORA2A_CRE, mouse_number=8042),
    SubjectEntry("20240116", "Adora2-Cre", "CONT", None, ADORA2A_CRE, mouse_number=8038),
    SubjectEntry("20240218", "Adora2-Cre", "6-OHDA", None, ADORA2A_CRE, mouse_number=8941),
    SubjectEntry("20240222", "Adora2-Cre", "LID ON", None, ADORA2A_CRE, mouse_number=8939),
    SubjectEntry("20240227", "Adora2-Cre", "LID OFF", None, ADORA2A_CRE, mouse_number=8943),
    SubjectEntry("20240303", "Adora2-Cre", "CONT", None, ADORA2A_CRE, mouse_number=8940),
    # Behaviour mice, under the lab's numbers: AIM scores and rotation videos (Figures 7 and 8, Supplementary Figure 3)
    SubjectEntry("20171215", "eGFP", "CONT", C57BL6, DRD2_EGFP, mouse_number=465),
    SubjectEntry("20171215", "eGFP", "CONT", C57BL6, DRD2_EGFP, mouse_number=466),
    SubjectEntry("20171215", "eGFP", "CONT", C57BL6, DRD2_EGFP, mouse_number=467),
    SubjectEntry("20180223", "eGFP", "CONT", C57BL6, DRD2_EGFP, mouse_number=1590),
    SubjectEntry("20180320", "eGFP", "CONT", C57BL6, DRD2_EGFP, mouse_number=2586),
    SubjectEntry("20180320", "eGFP", "CONT", C57BL6, DRD2_EGFP, mouse_number=2587),
    SubjectEntry("20180412", "eGFP", "CONT", C57BL6, DRD2_EGFP, mouse_number=2636),
    SubjectEntry("20180412", "eGFP", "CONT", C57BL6, DRD2_EGFP, mouse_number=2638),
    SubjectEntry("20180718", "eGFP", "CONT", C57BL6, DRD2_EGFP, mouse_number=5871),
    SubjectEntry("20180905", "eGFP", "CONT", C57BL6, DRD2_EGFP, mouse_number=5940),
    SubjectEntry("20171215", "CDGIko X eGFP", "ko", C57BL6, CDGI_KO_DRD2_EGFP, mouse_number=549),
    SubjectEntry("20180223", "CDGIko X eGFP", "ko", C57BL6, CDGI_KO_DRD2_EGFP, mouse_number=1944),
    SubjectEntry("20180223", "CDGIko X eGFP", "ko", C57BL6, CDGI_KO_DRD2_EGFP, mouse_number=1945),
    SubjectEntry("20180412", "CDGIko X eGFP", "ko", C57BL6, CDGI_KO_DRD2_EGFP, mouse_number=2828),
    SubjectEntry("20180905", "CDGIko X eGFP", "ko", C57BL6, CDGI_KO_DRD2_EGFP, mouse_number=6571),
    SubjectEntry("20180905", "CDGIko X eGFP", "ko", C57BL6, CDGI_KO_DRD2_EGFP, mouse_number=6572),
    SubjectEntry("20201103", "CDGIko X eGFP", "ko", C57BL6, CDGI_KO_DRD2_EGFP, mouse_number=2516),
    SubjectEntry("20201103", "CDGIko X eGFP", "ko", C57BL6, CDGI_KO_DRD2_EGFP, mouse_number=2517),
    SubjectEntry("20201103", "CDGIko X eGFP", "ko", C57BL6, CDGI_KO_DRD2_EGFP, mouse_number=2519),
    SubjectEntry("20230927", "Adora2-Cre", "CONT", None, ADORA2A_CRE, mouse_number=5844),
    SubjectEntry("20230927", "Adora2-Cre", "M1R CRSPR", None, ADORA2A_CRE, mouse_number=5845),
    SubjectEntry("20230927", "Adora2-Cre", "M1R CRSPR", None, ADORA2A_CRE, mouse_number=5846),
    SubjectEntry("20231208", "Adora2-Cre", "CONT", None, ADORA2A_CRE, mouse_number=7629),
    SubjectEntry("20231208", "Adora2-Cre", "CONT", None, ADORA2A_CRE, mouse_number=7830),
    SubjectEntry("20231208", "Adora2-Cre", "CONT", None, ADORA2A_CRE, mouse_number=7831),
    SubjectEntry("20231208", "Adora2-Cre", "M1R CRSPR", None, ADORA2A_CRE, mouse_number=7630),
    SubjectEntry("20231208", "Adora2-Cre", "M1R CRSPR", None, ADORA2A_CRE, mouse_number=7829),
    SubjectEntry("20240117", "Adora2-Cre", "M1R CRSPR", None, ADORA2A_CRE, mouse_number=8592),
    SubjectEntry("20240117", "Adora2-Cre", "M1R CRSPR", None, ADORA2A_CRE, mouse_number=8594),
    SubjectEntry("20240223", "Adora2-Cre", "CONT", None, ADORA2A_CRE, mouse_number=4790),
    SubjectEntry("20240223", "Adora2-Cre", "CONT", None, ADORA2A_CRE, mouse_number=4793),
    SubjectEntry("20240223", "Adora2-Cre", "CONT", None, ADORA2A_CRE, mouse_number=9366),
    SubjectEntry("20240223", "Adora2-Cre", "CONT", None, ADORA2A_CRE, mouse_number=9367),
    SubjectEntry("20240223", "Adora2-Cre", "CONT", None, ADORA2A_CRE, mouse_number=9368),
    SubjectEntry("20240503", "Adora2-Cre", "CONT", None, ADORA2A_CRE, mouse_number=2422),
    SubjectEntry("20240503", "Adora2-Cre", "CONT", None, ADORA2A_CRE, mouse_number=8240),
    SubjectEntry("20240503", "Adora2-Cre", "M1R CRSPR", None, ADORA2A_CRE, mouse_number=2421),
    SubjectEntry("20240503", "Adora2-Cre", "M1R CRSPR", None, ADORA2A_CRE, mouse_number=2423),
    SubjectEntry("20240503", "Adora2-Cre", "M1R CRSPR", None, ADORA2A_CRE, mouse_number=2424),
    SubjectEntry("20240503", "Adora2-Cre", "M1R CRSPR", None, ADORA2A_CRE, mouse_number=2425),
    SubjectEntry(
        "20171120",
        None,
        "ko",
        C57BL6,
        CDGI_KO,
        mouse_number=4123,
        note="Not in Data Connections, and the AIM workbook does not label it; its five session totals match a knockout column of the lab's Figure 7J table. Its line is not recorded.",
    ),
    SubjectEntry(
        "20171120",
        None,
        "ko",
        C57BL6,
        CDGI_KO,
        mouse_number=4124,
        note="Not in Data Connections, and the AIM workbook does not label it; its five session totals match a knockout column of the lab's Figure 7J table. Its line is not recorded.",
    ),
)
