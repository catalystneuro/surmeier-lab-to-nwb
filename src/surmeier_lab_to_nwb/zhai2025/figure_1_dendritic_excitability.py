import logging
import re
import warnings
from pathlib import Path

from neuroconv.tools import configure_and_write_nwbfile

# Suppress specific warnings from pynwb or tifffile
warnings.filterwarnings("ignore", message="invalid value encountered in divide")
warnings.filterwarnings("ignore", message=".*no datetime before year 1.*")
logging.getLogger("tifffile").setLevel(logging.ERROR)


from surmeier_lab_to_nwb.zhai2025 import (
    PrairieViewIntracellularRecordingInterface,
    PrairieViewLineScanInterface,
)

raw_data_folder = Path(
    "/media/heberto/One Touch/Surmeier-CN-data-share/consolidated_data/LID_paper_Zhai_2025/Raw data for Figs"
)
figure_1_folder = raw_data_folder / "Figure 1"
dendritic_excitability_folder = figure_1_folder / "Dendritic excitability"

nwbfile = None

condition_list = ["LID off-state", "LID on-state", "LID on-state with SCH"]
format_condition_name = {
    "LID off-state": "LIDOffState",
    "LID on-state": "LIDOnState",
    "LID on-state with SCH": "LIDOnStateWithSCH",
}
format_proximal_or_distal = {"prox": "Proximal", "dist": "Distal"}

folder_trial_pattern = re.compile(
    r".*_Cell\d+_(?P<location_prefix>real)?(?P<direction>dist|prox)(?P<number>\d+)_trio(?:_SCH|-)+(?P<trial_id>\d+)"
)


for condition in condition_list:
    condition_folder = dendritic_excitability_folder / condition
    assert condition_folder.exists(), f"Condition folder does not exist: {condition_folder}"

    print("=" * 50)
    print(f"Processing condition folder: {condition_folder.name}")
    condition = format_condition_name.get(condition, condition)

    date_folder_path_list = [p for p in condition_folder.iterdir() if p.is_dir() and p.is_dir()]

    for date_folder_path in date_folder_path_list:
        print("=" * 50)
        print(f"Processing date folder: {date_folder_path.name}")
        print("_" * 50)
        cell_folder_path_list = [p for p in date_folder_path.iterdir() if p.is_dir()]

        # Sometimes the recording sessions are not nested
        experiments_are_non_nested = any(p for p in cell_folder_path_list if "Cell" in p.name)
        if experiments_are_non_nested:
            cell_folder_path_list = [date_folder_path]

        for cell_folder_path in cell_folder_path_list:
            print(f"Processing cell folder: {cell_folder_path.name}")
            cell_id = cell_folder_path.name

            non_trial_folders = [".DS_Store", ".DS_Store"]
            trial_folder_path_list = [
                p for p in cell_folder_path.iterdir() if p.is_dir() and p.name not in non_trial_folders
            ]

            for trial_folder_path in trial_folder_path_list:
                match = folder_trial_pattern.match(trial_folder_path.name)
                location_prefix = match.group("location_prefix") or ""
                position = match.group("direction")
                number = match.group("number")
                trial_id = match.group("trial_id")

                xml_metadata_file_path = trial_folder_path / f"{trial_folder_path.name}.xml"
                xml_recording_file_file_path = (
                    trial_folder_path / f"{trial_folder_path.name}_Cycle00001_VoltageRecording_001.xml"
                )

                assert xml_metadata_file_path.exists(), f"Metadata file not found: {xml_metadata_file_path}"
                assert (
                    xml_recording_file_file_path.exists()
                ), f"Recording file not found: {xml_recording_file_file_path}"

                position = format_proximal_or_distal.get(position, position)
                position = f"{location_prefix}{position}{number}"  # e.g., Proximal1, Distal2
                interface_recording = PrairieViewIntracellularRecordingInterface(
                    xml_recording_file_file_path,
                    cell_id=cell_id,
                    trial_id=trial_id,
                    position=position,
                    condition=condition,
                )
                interface_line_scan = PrairieViewLineScanInterface(
                    xml_metadata_file_path,
                    cell_id=cell_id,
                    trial_id=trial_id,
                    position=position,
                    condition=condition,
                )

                if nwbfile is None:
                    nwbfile = interface_recording.create_nwbfile()
                    interface_line_scan.add_to_nwbfile(nwbfile=nwbfile)
                else:
                    interface_recording.add_to_nwbfile(nwbfile=nwbfile)
                    interface_line_scan.add_to_nwbfile(nwbfile=nwbfile)

    print("\n")

nwbfile_path = "nwb_file_test.nwb"
configure_and_write_nwbfile(nwbfile=nwbfile, nwbfile_path=nwbfile_path)
