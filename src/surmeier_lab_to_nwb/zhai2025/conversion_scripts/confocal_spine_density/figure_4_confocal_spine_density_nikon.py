"""
Figure 4 Confocal Spine Density Conversion Script (Nikon) - Zhai et al. 2025
=======================================================================

This script converts the high-resolution confocal z-stacks of Figure 4J and Supplementary Figure 5 into NWB, one file
per mouse with every dendrite segment the lab imaged in it. The mice are Adora2-Cre mice whose iSPNs were sparsely
labeled with a Cre-dependent eGFP virus, in four groups: control, 6-OHDA, LID off-state and LID on-state.

For detailed experimental context, protocols, data structure, and analysis methods,
see: /src/surmeier_lab_to_nwb/zhai2025/conversion_notes_folder/figure_4_conversion_notes.md
"""

import re
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List
from zoneinfo import ZoneInfo

from neuroconv.tools import configure_and_write_nwbfile
from neuroconv.tools.nwb_helpers import make_nwbfile_from_metadata
from neuroconv.utils import dict_deep_update, load_dict_from_file
from pynwb import NWBFile

from surmeier_lab_to_nwb.zhai2025.conversion_scripts.conversion_utils import (
    generate_canonical_session_id,
    surgery_and_pharmacology,
)
from surmeier_lab_to_nwb.zhai2025.interfaces.image_stack_interfaces import (
    NikonImageStackInterface,
)
from surmeier_lab_to_nwb.zhai2025.subjects import SUBJECTS

# "3824-slide 2-slice 2-cell 2 proxi", "5108-slice 2-cell 1 den3 proxi", "5104-slice 3-cell 1 proxi3",
# "8038-slide 1-slice 2-cell 1 proxi_0001", "5105-slide 2-slice 3-cell 1den2 proxi": mouse, slide (some names leave it
# out), slice, cell, and the dendrite when the cell has more than one imaged
_RE_STACK_NAME = re.compile(
    r"^(?P<mouse>\d+)-(?:slide (?P<slide>\d+)-)?slice (?P<slice>\d+)-cell (?P<cell>\d+) ?(?:den ?(?P<dendrite>\d+) )?"
    r"proxi(?:3|_0001)?$"
)
_STATE = {"control": "UnlesionedControl", "6-OHDA": "LesionedControl", "off-state": "OffState", "on-state": "OnState"}
_GROUP_TEXT = {
    "UnlesionedControl": "unlesioned control: no 6-OHDA lesion and no levodopa",
    "LesionedControl": "6-OHDA lesioned, without levodopa",
    "OffState": "6-OHDA lesioned and given levodopa until dyskinetic, perfused in the LID off-state",
    "OnState": "6-OHDA lesioned and given levodopa until dyskinetic, perfused in the LID on-state",
}
TIMEZONE = ZoneInfo("America/Chicago")
# Stacks whose clocks disagree: the text NIS-Elements writes gives 2024-02-23 and the frame clock 2024-06-04, but the
# data share dates the files 2024-01-03, before both, so when they were acquired is uncertain
_UNCERTAIN_ACQUISITION_TIME = {
    "8040-slide 1-slice 3-cell 2 proxi.nd2",
    "8040-slide 1-slice 3-cell 2 den2 proxi.nd2",
    "8041-slide 2-slice 3-cell 2 proxi.nd2",
}


def parse_nd2_filename(nd2_file: Path) -> Dict[str, Any]:
    """
    Read the mouse, slide, slice, cell and dendrite of a stack from its file name; raises on a name it cannot read.

    Returns the mouse number and the name of the stack's Images container, e.g. "ImagesSlide2Slice2Cell2Dendrite1".
    """
    match = _RE_STACK_NAME.match(nd2_file.stem)
    if match is None:
        raise ValueError(f"Unrecognized confocal stack name: {nd2_file}")
    slide, slice_number, cell = match.group("slide"), match.group("slice"), match.group("cell")
    dendrite = match.group("dendrite") or "1"
    location = f"slide {slide}, slice {slice_number}" if slide else f"slice {slice_number}"
    return {
        "mouse_number": int(match.group("mouse")),
        "container_name": f"Images{f'Slide{slide}' if slide else ''}Slice{slice_number}Cell{cell}Dendrite{dendrite}",
        "description": f"Confocal z-stack of proximal dendrite {dendrite} of cell {cell} ({location})",
    }


def get_acquisition_time(nd2_file: Path, interface: NikonImageStackInterface) -> datetime | None:
    """When a stack was acquired, from its frame clock checked against the date in its text; None where uncertain."""
    if nd2_file.name in _UNCERTAIN_ACQUISITION_TIME:
        return None
    acquisition_time = interface.get_session_start_time().astimezone(TIMEZONE)
    text_date = datetime.strptime(interface.file_handle.text_info["date"].split()[0], "%m/%d/%Y").date()
    if acquisition_time.date() != text_date:
        raise ValueError(f"{nd2_file}: the frame clock says {acquisition_time}, the text {text_date}")
    return acquisition_time


def group_stacks_by_mouse(condition_folder: Path) -> Dict[int, List[Path]]:
    """Group the stacks of one condition folder by the mouse their file names give."""
    stacks_per_mouse: Dict[int, List[Path]] = {}
    for nd2_file in sorted(condition_folder.glob("*.nd2")):
        stacks_per_mouse.setdefault(parse_nd2_filename(nd2_file)["mouse_number"], []).append(nd2_file)
    return stacks_per_mouse


def convert_mouse_to_nwbfile(nd2_files: List[Path], condition: str) -> NWBFile:
    """
    Convert every confocal stack of one mouse into one NWB file.

    Parameters
    ----------
    nd2_files : list of Path
        The mouse's ND2 stacks
    condition : str
        Its condition folder: "control", "6-OHDA", "off-state" or "on-state"

    Returns
    -------
    NWBFile
        The populated NWB file object
    """
    state = _STATE[condition]
    stacks = [(nd2_file, parse_nd2_filename(nd2_file), NikonImageStackInterface(nd2_file)) for nd2_file in nd2_files]
    (mouse_number,) = {info["mouse_number"] for _, info, _ in stacks}
    subject = next(subject for subject in SUBJECTS if subject.mouse_number == mouse_number)
    container_names = [info["container_name"] for _, info, _ in stacks]
    if len(set(container_names)) != len(container_names):
        raise ValueError(f"Two stacks of mouse {mouse_number} have the same slide, slice, cell and dendrite")

    acquisition_times = [get_acquisition_time(nd2_file, interface) for nd2_file, _, interface in stacks]
    known_times = [time for time in acquisition_times if time is not None]
    session_start_time = min(known_times)
    imaging_days = sorted({time.strftime("%Y-%m-%d") for time in known_times})
    session_id = generate_canonical_session_id(
        cell_type="iSPN",  # Indirect pathway SPN
        state=state,
        pharm="none",  # No pharmacology
        geno="WT",  # Not a knockout
        date=session_start_time.strftime("%Y%m%d"),
    )

    # Load general and session-specific metadata from YAML files
    general_metadata_path = Path(__file__).parent.parent.parent / "general_metadata.yaml"
    general_metadata = load_dict_from_file(general_metadata_path)

    session_metadata_path = Path(__file__).parent.parent.parent / "session_specific_metadata.yaml"
    session_metadata_template = load_dict_from_file(session_metadata_path)
    script_template = session_metadata_template["figure_4_confocal_spine_density_nikon"]
    surgery, pharmacology = surgery_and_pharmacology(state, general_metadata)
    if state == "OnState":
        pharmacology += " " + script_template["NWBFile"]["on_state_pharmacology"]

    session_specific_metadata = {
        "NWBFile": {
            "session_description": script_template["NWBFile"]["session_description"].format(
                num_stacks=len(stacks),
                mouse_number=mouse_number,
                condition=_GROUP_TEXT[state],
                imaging_days=", ".join(imaging_days),
            ),
            "experiment_description": script_template["NWBFile"]["experiment_description"],
            "session_start_time": session_start_time,
            "session_id": session_id,
            "surgery": surgery + " " + script_template["NWBFile"]["surgery_addition"],
            "pharmacology": pharmacology,
            "slices": script_template["NWBFile"]["slices"],
            "keywords": script_template["NWBFile"]["keywords"],
        },
        "Subject": {
            "subject_id": subject.subject_id,
            "description": script_template["Subject"]["description"].format(
                mouse_number=mouse_number, group=_GROUP_TEXT[state]
            ),
            "strain": subject.strain,
            "genotype": subject.genotype,
        },
    }

    # Merge general metadata with session-specific metadata
    merged_metadata = dict_deep_update(general_metadata, session_specific_metadata)
    if subject.strain is None:  # the paper gives no background for this line
        del merged_metadata["Subject"]["strain"]

    # Create NWB file using neuroconv helper function
    nwbfile = make_nwbfile_from_metadata(merged_metadata)

    for (nd2_file, info, interface), acquisition_time in zip(stacks, acquisition_times):
        if acquisition_time is None:
            acquired = (
                "whose clocks disagree (2024-02-23 in the text NIS-Elements writes, 2024-06-04 from its frames) while "
                "the data share dates it 2024-01-03, so when it was acquired is uncertain; acquired"
            )
        else:
            acquired = f"acquired {acquisition_time:%Y-%m-%d %H:%M} (Chicago time)"
        description = (
            f"{info['description']}, the lab's file '{nd2_file.name}', {acquired} on a Nikon AXR confocal with a 60x "
            "oil immersion objective (NA 1.49), 0.09 µm pixels and 0.125 µm z steps."
        )
        interface_metadata = {
            "Images": {info["container_name"]: {"name": info["container_name"], "description": description}}
        }
        interface.add_to_nwbfile(nwbfile, metadata=interface_metadata)

    return nwbfile


if __name__ == "__main__":
    import logging

    from tqdm import tqdm

    logging.getLogger("nd2").setLevel(logging.ERROR)

    # The four Figure 4I stacks are byte-identical copies of Figure 4J stacks, so only Figure 4J is converted
    data_path = Path("./link_to_raw_data/Figure 4_SF1B_SF5/Confocal spine density/Fig 4J and Suppl Fig 5/raw")

    # Create output directory
    root_dir = Path(__file__).parent.parent.parent.parent.parent.parent  # Go up to repo root
    nwb_files_dir = root_dir / "nwb_files" / "confocal_spine_density" / "figure_4"
    nwb_files_dir.mkdir(parents=True, exist_ok=True)

    for condition in _STATE:
        condition_path = data_path / condition
        assert condition_path.exists(), f"Condition path {condition_path} does not exist"
        stacks_per_mouse = group_stacks_by_mouse(condition_path)
        for nd2_files in tqdm(
            stacks_per_mouse.values(), desc=f"Converting Figure 4J confocal {condition}", unit=" mouse"
        ):
            nwbfile = convert_mouse_to_nwbfile(nd2_files, condition)
            nwbfile_path = nwb_files_dir / f"{nwbfile.subject.subject_id}.nwb"
            configure_and_write_nwbfile(nwbfile, nwbfile_path=nwbfile_path)
