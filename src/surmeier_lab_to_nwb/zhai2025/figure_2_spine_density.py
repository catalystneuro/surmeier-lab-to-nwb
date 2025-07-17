import re
from datetime import datetime
from pathlib import Path
from typing import Any, Dict
from zoneinfo import ZoneInfo

from neuroconv.datainterfaces import ImageInterface
from neuroconv.tools import configure_and_write_nwbfile
from pynwb import NWBFile


def extract_date_from_tiff_filename(tiff_file: Path) -> datetime:
    """
    Extract date from TIFF filename.

    Example filenames:
    - "20_ZSeries-20170707_Cell2_prox12-001_Cycle00001_Ch1_#.ome_Z026.tif"
    - "20_ZSeries-20190411_Cell1_dist1-001_Cycle00001_Ch1_#.ome_Z01.tif"

    Args:
        tiff_file: Path to TIFF file

    Returns:
        datetime object with extracted date
    """
    filename = tiff_file.name

    # Look for date pattern YYYYMMDD in filename
    date_match = re.search(r"(\d{8})", filename)
    if date_match:
        date_str = date_match.group(1)
        year = int(date_str[:4])
        month = int(date_str[4:6])
        day = int(date_str[6:8])
    else:
        raise ValueError(f"Could not extract date from filename: {filename}")

    # Illinois is in Central Time Zone
    central_tz = ZoneInfo("America/Chicago")

    return datetime(year, month, day, 0, 0, 0, tzinfo=central_tz)


def parse_session_info(session_folder: Path) -> Dict[str, Any]:
    """
        Parse session information, preferring TIFF filename dates over folder names.
    x
    """

    # Look for TIFF files in subfolders
    for subfolder in session_folder.iterdir():
        if subfolder.is_dir():
            tiff_files = [f for f in subfolder.iterdir() if f.suffix.lower() == ".tif"]
            tiff_date = extract_date_from_tiff_filename(tiff_files[0])
            if tiff_date:
                session_start_time = tiff_date
            else:
                raise ValueError(f"Could not extract date from TIFF file in folder: {subfolder.name}")

    # Extract animal ID from folder name
    # TODO: this is a guess in the moment, need to confirm
    folder_name = session_folder.name
    if "2019" in folder_name:
        animal_id = folder_name[8:] if len(folder_name) > 8 else "unknown"
    else:
        animal_id = folder_name[4:] if len(folder_name) > 4 else "unknown"

    return {
        "session_start_time": session_start_time,
        "animal_id": animal_id,
        "date_str": f"{session_start_time.year}-{session_start_time.month:02d}-{session_start_time.day:02d}",
        "session_id": f"{session_start_time.year}{session_start_time.month:02d}{session_start_time.day:02d}_{animal_id}",
    }


def parse_container_info(subfolder_name: str) -> Dict[str, str]:
    """
    Parse container information from subfolder names.

    Examples:
    - "Decon_20190411_Cell1_dist1" -> Cell 1, distal dendrite 1
    - "Decon_20170706_Cell1_prox12" -> Cell 1, proximal dendrites 1&2

    """
    # Extract cell number and location
    cell_match = re.search(r"Cell(\d+)", subfolder_name)
    cell_number = cell_match.group(1) if cell_match else "unknown"

    # Extract location (proximal/distal) and dendrite number
    if "dist" in subfolder_name:
        location_match = re.search(r"dist(\d+)", subfolder_name)
        location = "Distal"
        dendrite_num = location_match.group(1) if location_match else "1"
    elif "prox" in subfolder_name:
        location_match = re.search(r"prox(\d+)", subfolder_name)
        location = "Proximal"
        dendrite_num = location_match.group(1) if location_match else "1"
    elif "Prox" in subfolder_name:
        location_match = re.search(r"Prox_(\d+)", subfolder_name)
        location = "Proximal"
        dendrite_num = location_match.group(1) if location_match else "1"
    else:
        raise ValueError(f"Could not determine location from subfolder name: {subfolder_name}")

    container_name = f"ImagesCell{cell_number}{location}Dendrite{dendrite_num}"

    # TODO: figure out what real means here
    # Most likely the experiment was done twice and this is "correct"
    if "real" in subfolder_name:
        container_name += "Real"

    description = (
        f"Image stack of {location} dendrite {dendrite_num} from dSPN cell {cell_number} for spine density analysis"
    )

    return {
        "container_name": container_name,
        "description": description,
    }


def convert_data_to_nwb(session_folder_path: Path, condition: str) -> NWBFile:
    """
    session_folder_path: Path to the top level folders in the conditions for spiny density figure 2
    """
    # Parse session information
    session_info = parse_session_info(session_folder_path)
    print(f"Session date: {session_info['date_str']}, Animal: {session_info['animal_id']}")

    # Create NWB file with appropriate metadata
    nwbfile = NWBFile(
        session_description=(
            f"Dendritic spine density assessment in direct pathway spiny projection neurons (dSPNs) "
            f"for condition {condition}. Images of dendritic segments "
            "(proximal: ~40 μm from soma; distal: > 80 μm from soma) were acquired with 0.15 μm pixels "
            "with 0.3 μm z-steps using two-photon microscopy. Images were deconvolved in AutoQuant X3.0.4"
        ),
        identifier=f"zhai2025_fig2_spine_density_{session_info['session_id']}",
        session_start_time=session_info["session_start_time"],
        experimenter=[
            "Zhai, Shenyu",
            "Cui, Qiaoling",
            "Wokosin, David",
            "Sun, Linqing",
            "Tkatch, Tatiana",
            "Crittenden, Jill R.",
            "Graybiel, Ann M.",
            "Surmeier, D. James",
        ],
        lab="Surmeier Lab",
        institution="Northwestern University",
        experiment_description=(
            f"Assessment of dendritic spine density changes in dSPNs during condition {condition}. "
            "This experiment is part of Figure 2 from Zhai et al. 2025, investigating how LID affects "
            "dendritic spine density in direct pathway neurons."
        ),
        session_id=session_info["session_id"],
    )

    # TODO: add microscope device

    subfolders = [f for f in session_folder_path.iterdir() if f.is_dir()]

    print(f"Found {len(subfolders)} image stacks in session")

    # Process each image stack
    for subfolder in subfolders:
        print(f"  Processing: {subfolder.name}")

        # Parse container information
        container_info = parse_container_info(subfolder.name)

        # Get all TIFF files in the folder
        tiff_files = sorted([f for f in subfolder.iterdir() if f.suffix.lower() == ".tif"])
        assert len(tiff_files) > 0, f"No TIFF files found in {subfolder.name}"

        print(f"Found {len(tiff_files)} TIFF files in {subfolder.name}")

        interface = ImageInterface(file_paths=tiff_files)
        interface.add_to_nwbfile(nwbfile=nwbfile, container_name=container_info["container_name"])
        print(f"Successfully added image stack: {container_info['container_name']}")

    return nwbfile


if __name__ == "__main__":
    import logging

    from tqdm import tqdm

    # Control verbose output
    verbose = False  # Set to True for detailed output

    logging.getLogger("tifffile").setLevel(logging.ERROR)

    # Define the base path to the data
    base_path = Path(
        "/media/heberto/One Touch/Surmeier-CN-data-share/consolidated_data/LID_paper_Zhai_2025/Raw data for Figs/Figure 2_SF1A/Spine density/"
    )
    conditions = ["PD dSPN", "LID off-state dSPN", "control dSPN", "LID on-state dSPN"]
    for condition in conditions:
        condition_path = base_path / condition
        assert condition_path.exists(), f"Base path does not exist: {base_path}"

        print(f"Processing spine density data for: {condition=}")

        # Get all session folders
        session_folders = [f for f in condition_path.iterdir() if f.is_dir()]
        session_folders.sort()

        print(f"Found {len(session_folders)} session folders")

        # Use tqdm for progress bar when verbose is disabled
        session_iterator = (
            tqdm(session_folders, desc=f"Processing {condition}", disable=verbose) if not verbose else session_folders
        )

        for session_folder in session_iterator:
            if verbose:
                print(f"\nProcessing session: {session_folder.name}")

            # Convert data to NWB format
            nwbfile = convert_data_to_nwb(
                session_folder=session_folder,
                condition=condition,
            )

            nwbfile_path = Path(f"spine_density_{condition}_{session_folder.name}.nwb")
            configure_and_write_nwbfile(nwbfile, nwbfile_path=nwbfile_path)
            print(f"Successfully saved: {nwbfile_path}")
