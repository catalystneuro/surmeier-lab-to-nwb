"""
Figure 7 Behavioral Video Conversion Script

This script processes behavioral video recordings from Figure 7 of the Zhai et al. 2025 paper,
specifically the contralateral rotation behavior for dyskinesia assessment in CDGI knockout mice.

The script handles:
- Behavioral video files (.mov/.mp4) from multiple experimental sessions
- Video metadata extraction and organization by animal and session
- Integration with NWB format using appropriate video interfaces
- Session timing and L-DOPA treatment context
"""

import uuid
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List
from zoneinfo import ZoneInfo

from neuroconv.datainterfaces import ExternalVideoInterface
from neuroconv.tools import configure_and_write_nwbfile
from neuroconv.tools.nwb_helpers import make_nwbfile_from_metadata
from neuroconv.utils import dict_deep_update, load_dict_from_file
from pynwb import NWBFile
from pynwb.epoch import TimeIntervals
from tqdm import tqdm

from surmeier_lab_to_nwb.zhai2025.conversion_scripts.conversion_utils import (
    format_condition,
    generate_canonical_session_id,
)


def find_video_sessions(video_base_path: Path) -> Dict[str, Path]:
    """
    Find all video session directories for Figure 7.

    Parameters
    ----------
    video_base_path : Path
        Base path to CDGI KO video directories

    Returns
    -------
    Dict[str, Path]
        Dictionary mapping session dates to directory paths
    """
    sessions = {}

    for session_dir in video_base_path.iterdir():
        if session_dir.is_dir() and session_dir.name.endswith("-"):
            # Parse date from folder name (format: MMDDYYYY-)
            date_str = session_dir.name[:-1]  # Remove trailing '-'
            if len(date_str) == 8:
                month = date_str[:2]
                day = date_str[2:4]
                year = date_str[4:8]
                session_date = f"{year}-{month}-{day}"
                sessions[session_date] = session_dir

    return sessions


def extract_video_metadata_figure7(video_path: Path) -> Dict[str, Any]:
    """
    Extract metadata from Figure 7 video filename.

    Parameters
    ----------
    video_path : Path
        Path to video file

    Returns
    -------
    Dict[str, Any]
        Dictionary with video metadata including animal_id, session, time_point
    """
    filename = video_path.stem  # Remove extension

    # Common patterns for Figure 7 CDGI videos:
    # ET2586_1_20min, 1944_1_20min, 1590_1_20min, etc.

    metadata = {
        "animal_id": "unknown",
        "session_number": 1,
        "time_point_min": 0,
        "file_format": video_path.suffix.lower(),
        "genotype": "CDGI KO",  # All Figure 7 videos are CDGI knockout
        "treatment": "L-DOPA",
    }

    # Try to parse animal ID and time point
    parts = filename.split("_")

    if len(parts) >= 3:
        # Pattern: AnimalID_Session_TimePoint
        animal_id = parts[0]
        session_num = parts[1]
        time_str = parts[2]

        metadata["animal_id"] = animal_id
        metadata["session_number"] = int(session_num) if session_num.isdigit() else 1

        # Extract time point (remove 'min' suffix)
        if time_str.endswith("min"):
            time_point = time_str[:-3]
            if time_point.isdigit():
                metadata["time_point_min"] = int(time_point)

    elif len(parts) >= 2:
        # Pattern: AnimalID_TimePoint
        animal_id = parts[0]
        time_str = parts[1]

        metadata["animal_id"] = animal_id

        if time_str.endswith("min"):
            time_point = time_str[:-3]
            if time_point.isdigit():
                metadata["time_point_min"] = int(time_point)

    return metadata


def group_videos_by_animal(session_path: Path) -> Dict[str, List[Path]]:
    """
    Group video files by animal ID within a session.

    Parameters
    ----------
    session_path : Path
        Path to session directory

    Returns
    -------
    Dict[str, List[Path]]
        Dictionary mapping animal IDs to lists of video files
    """
    animals = {}

    # Search through all subsession directories
    for subsession_dir in session_path.iterdir():
        if subsession_dir.is_dir():
            for video_file in subsession_dir.iterdir():
                if video_file.suffix.lower() in [".mov", ".mp4"]:
                    metadata = extract_video_metadata_figure7(video_file)
                    animal_id = metadata["animal_id"]

                    if animal_id != "unknown":
                        if animal_id not in animals:
                            animals[animal_id] = []
                        animals[animal_id].append(video_file)

    return animals


def convert_session_to_nwbfile(
    session_date: str,
    animal_id: str,
    video_files: List[Path],
    verbose: bool = False,
) -> NWBFile:
    """
    Convert behavioral video recordings from a single experimental session for a single subject to NWB format.

    This function processes video files from a single behavioral session (one animal on one date) containing
    multiple time points of video recordings during L-DOPA treatment. We organize by session (animal + date)
    because:

    1. **Experimental design**: Each session represents a complete L-DOPA behavioral assessment for one animal
       on one experimental date, with videos recorded at multiple time points (20, 40, 60, 80, 100, 120 min)

    2. **Temporal coherence**: All videos within a session share the same experimental context (same animal,
       same treatment protocol, same experimental conditions) and form a temporally coherent sequence

    3. **NWB structure**: This allows proper temporal alignment of all video recordings within a single NWB file,
       enabling synchronized analysis of behavioral changes across the entire treatment session

    4. **Data integrity**: Each session forms a complete behavioral assessment unit that should be analyzed
       together, making it the natural granularity for NWB file organization

    Parameters
    ----------
    session_date : str
        Session date in YYYY-MM-DD format
    animal_id : str
        Animal identifier (e.g., "ET2586", "1944")
    video_files : List[Path]
        List of video file paths for this animal session
    verbose : bool, default=False
        Whether to print detailed processing information

    Returns
    -------
    NWBFile
        NWB file object containing all video data from the session
    """
    if not video_files:
        raise ValueError(f"No video files found for {animal_id} on {session_date}")

    # Parse session start time
    session_start_time = datetime.strptime(session_date, "%Y-%m-%d")
    session_start_time = session_start_time.replace(tzinfo=ZoneInfo("US/Central"))

    # Map genotype to standardized condition for centralized format_condition dictionary
    genotype = "CDGI KO"  # All Figure 7 videos are CDGI knockout
    genotype_to_condition = {"CDGI KO": "knockout", "KO": "knockout", "WT": "control"}

    # Get standardized condition from genotype
    standardized_condition = genotype_to_condition.get(genotype, "knockout")

    # Use centralized format_condition dictionary
    condition_human_readable = format_condition[standardized_condition]["human_readable"]

    # Create canonical session ID with explicit parameters
    timestamp = session_start_time.strftime("%Y%m%d")

    # Map condition to state - CDGI KO videos capture behavior during ON state
    # All Figure 7 videos are from L-DOPA treated animals (ON state)
    state = "ON"  # All videos are during L-DOPA treatment

    session_id = generate_canonical_session_id(
        fig="F7",
        compartment="behav",  # Whole-animal behaviour
        measurement="video",  # Raw video
        spn_type="pan",  # Non cell-specific
        state=state,
        pharmacology="none",  # No pharmacology
        genotype="CDGIKO",  # CDGI knock-out
        timestamp=timestamp,
    )

    # Load general and session-specific metadata from YAML files
    general_metadata_path = Path(__file__).parent.parent.parent / "general_metadata.yaml"
    general_metadata = load_dict_from_file(general_metadata_path)

    session_metadata_path = Path(__file__).parent.parent.parent / "session_specific_metadata.yaml"
    session_metadata_template = load_dict_from_file(session_metadata_path)
    script_template = session_metadata_template["figure_7_behavioral_videos"]

    # Create session-specific metadata from template with runtime substitutions
    conversion_specific_metadata = {
        "NWBFile": {
            "session_description": script_template["NWBFile"]["session_description"].format(animal_id=animal_id),
            "identifier": str(uuid.uuid4()),
            "session_start_time": session_start_time,
            "session_id": session_id,
            "keywords": script_template["NWBFile"]["keywords"],
        },
        "Subject": {
            "subject_id": f"CDGI_mouse_{animal_id}",
            "description": script_template["Subject"]["description"].format(animal_id=animal_id),
            "genotype": script_template["Subject"]["genotype"],
        },
    }

    # Merge general metadata with session-specific metadata
    metadata = dict_deep_update(general_metadata, conversion_specific_metadata)

    # Create NWB file using neuroconv helper function
    nwbfile = make_nwbfile_from_metadata(metadata)

    # Sort videos by time point for proper temporal organization
    video_metadata_list = [(vid, extract_video_metadata_figure7(vid)) for vid in video_files]
    video_metadata_list.sort(key=lambda x: x[1]["time_point_min"])

    # Add time intervals for behavioral epochs
    behavioral_epochs = TimeIntervals(
        name="behavioral_epochs",
        description="Time intervals for behavioral video recordings during L-DOPA treatment",
    )

    # Add custom columns for behavioral video metadata
    behavioral_epochs.add_column(
        name="animal_id",
        description="Animal identifier",
    )
    behavioral_epochs.add_column(
        name="session_number",
        description="Session number within experimental sequence",
    )
    behavioral_epochs.add_column(
        name="time_point_min",
        description="Time point in minutes post L-DOPA injection",
    )
    behavioral_epochs.add_column(
        name="treatment",
        description="Treatment condition",
    )
    behavioral_epochs.add_column(
        name="video_file",
        description="Path to the video file",
    )

    # Track video names to ensure uniqueness
    video_name_counter = {}

    for video_file, vid_metadata in video_metadata_list:
        # Calculate video start time relative to session start
        time_offset_seconds = float(vid_metadata["time_point_min"] * 60)  # Convert minutes to seconds
        video_start_time = time_offset_seconds

        # Estimate video duration (typical behavioral videos are 1-3 minutes)
        # We'll use a default of 120 seconds, this can be refined with actual video analysis
        video_duration = 120.0  # seconds
        video_stop_time = video_start_time + video_duration

        # Add behavioral epoch
        behavioral_epochs.add_interval(
            start_time=video_start_time,
            stop_time=video_stop_time,
            animal_id=animal_id,
            session_number=vid_metadata["session_number"],
            time_point_min=vid_metadata["time_point_min"],
            treatment="L-DOPA",
            video_file=str(video_file),
        )

        # Create unique video name with counter for duplicates
        base_video_name = f"video_{video_file.stem}"
        if base_video_name in video_name_counter:
            video_name_counter[base_video_name] += 1
            video_name = f"{base_video_name}_{video_name_counter[base_video_name]}"
        else:
            video_name_counter[base_video_name] = 0
            video_name = base_video_name

        video_interface = ExternalVideoInterface(
            file_paths=[video_file],
            video_name=video_name,
            verbose=False,
        )

        # Add video to NWB file
        video_interface.add_to_nwbfile(
            nwbfile=nwbfile,
            metadata={
                "Behavior": {
                    "ExternalVideos": {
                        video_name: {
                            "description": f"Behavioral video at {vid_metadata['time_point_min']} minutes post L-DOPA treatment for contralateral rotation assessment",
                        }
                    }
                }
            },
        )

    # Add behavioral epochs to NWB file
    nwbfile.add_time_intervals(behavioral_epochs)

    return nwbfile


if __name__ == "__main__":
    import argparse

    # Set up argument parser
    parser = argparse.ArgumentParser(description="Convert Figure 7 behavioral videos to NWB format")
    parser.add_argument("--verbose", "-v", action="store_true", help="Enable verbose output")
    args = parser.parse_args()

    # Set up paths
    video_base_path = Path(
        "/media/heberto/One Touch/Surmeier-CN-data-share/consolidated_data/LID_paper_Zhai_2025/Raw data for Figs/Figure 7/contralateral rotations/CDGI KO videos"
    )
    output_base_path = Path("./nwb_files/videos/figure_7")

    # Create output directory
    output_base_path.mkdir(parents=True, exist_ok=True)

    # Setup logging
    verbose = False
    stub_test = True  # Set to True to process only first 2 files per condition for testing
    # Find all video sessions
    sessions = find_video_sessions(video_base_path)

    # Apply stub_test filtering if enabled
    if stub_test:
        sessions_items = list(sessions.items())[:2]
        sessions = dict(sessions_items)

    # Process each session
    for session_date, session_path in tqdm(
        sessions.items(), desc="Converting Figure7 Videos", disable=not verbose, unit=" session"
    ):

        # Group videos by animal
        animals = group_videos_by_animal(session_path)

        # Create NWB files for each animal
        for animal_id, video_files in animals.items():

            # Create output filename
            genotype = "CDGI KO"  # All Figure 7 videos are CDGI knockout
            genotype_safe = genotype.replace(" ", "_").replace("-", "_")
            output_filename = f"zhai2025_figure7_videos_{animal_id}_{session_date.replace('-', '')}_{genotype_safe}.nwb"
            output_path = output_base_path / output_filename

            # Create NWB file
            nwbfile = convert_session_to_nwbfile(
                session_date=session_date,
                animal_id=animal_id,
                video_files=sorted(video_files),
                verbose=verbose,
            )

            # Write the NWB file
            configure_and_write_nwbfile(nwbfile=nwbfile, nwbfile_path=output_path)
