"""
Figure 7 Behavioral Video Conversion Script

This script handles the conversion of behavioral video files from Figure 7 of the Zhai et al. 2025 paper.
Videos show contralateral rotation behavior for dyskinesia assessment in CDGI knockout mice.

The script processes:
- Behavioral video files (.mov/.mp4) organized by session and time point
- Video metadata and timing information
- Integration with NWB format using ExternalVideoInterface
"""

import re
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List
from zoneinfo import ZoneInfo

from neuroconv.datainterfaces import ExternalVideoInterface
from pynwb import NWBHDF5IO, NWBFile
from pynwb.file import Subject


def find_video_sessions(video_base_path: Path) -> Dict[str, Path]:
    """
    Find all video session directories.

    Parameters
    ----------
    video_base_path : Path
        Base path to video directories

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


def extract_video_metadata(video_path: Path) -> Dict[str, Any]:
    """
    Extract metadata from video filename.

    Parameters
    ----------
    video_path : Path
        Path to video file

    Returns
    -------
    Dict[str, Any]
        Dictionary with video metadata
    """
    filename = video_path.name

    # Extract animal ID and time point from filename
    # Patterns: AnimalID_Session_TimePoint.ext or similar

    # Try different patterns
    patterns = [
        r"(\d+)_(\d+)_(\d+)min",  # 1590_1_20min.mov
        r"ET(\d+)_(\d+)min",  # ET5871_20min.mov
        r"(\d+)_(\d+)_(\d+)min",  # 2516_1_20min.MP4
    ]

    metadata = {
        "animal_id": "unknown",
        "session_number": 1,
        "time_point_min": 0,
        "file_format": video_path.suffix.lower(),
    }

    for pattern in patterns:
        match = re.search(pattern, filename)
        if match:
            groups = match.groups()
            if len(groups) >= 2:
                metadata["animal_id"] = groups[0]
                if len(groups) == 3:
                    metadata["session_number"] = int(groups[1])
                    metadata["time_point_min"] = int(groups[2])
                else:
                    metadata["time_point_min"] = int(groups[1])
            break

    return metadata


def create_video_nwb_file(session_date: str, animal_id: str, video_files: List[Path], output_path: Path) -> None:
    """
    Create an NWB file with video data for a single animal session.

    Parameters
    ----------
    session_date : str
        Session date in YYYY-MM-DD format
    animal_id : str
        Animal identifier
    video_files : List[Path]
        List of video file paths for this animal
    output_path : Path
        Path for output NWB file
    """
    if not video_files:
        print(f"No video files found for {animal_id} on {session_date}")
        return

    # Create session start time
    session_start_time = datetime.strptime(session_date, "%Y-%m-%d")
    session_start_time = session_start_time.replace(tzinfo=ZoneInfo("US/Central"))

    # Create NWBFile
    nwbfile = NWBFile(
        session_description=f"Figure 7 Behavioral Videos - Contralateral rotation analysis for CDGI study",
        identifier=f"figure7_videos_{animal_id}_{session_date}",
        session_start_time=session_start_time,
        lab="Surmeier Lab",
        institution="Northwestern University",
        related_publications="https://doi.org/10.1016/j.cell.2024.11.015",
        experiment_description="Behavioral video recordings of contralateral rotation behavior for dyskinesia assessment in CDGI knockout mice following L-DOPA treatment",
    )

    # Add subject information (basic info since we don't have genotype from video files)
    subject = Subject(
        subject_id=f"ET#{animal_id}",
        species="Mus musculus",
        strain="C57BL/6J",
        description=f"CDGI study animal from behavioral video recordings",
    )
    nwbfile.subject = subject

    # Process each video file
    for video_file in video_files:
        metadata = extract_video_metadata(video_file)

        # Create video interface
        try:
            video_interface = ExternalVideoInterface(
                file_paths=[str(video_file)],
                verbose=False,
                video_name=f"behavioral_video_{metadata['time_point_min']}min",
            )

            # Get video metadata
            video_metadata = video_interface.get_metadata()
            video_metadata["NWBFile"].update(session_start_time=session_start_time)

            # Add video to NWB file
            video_interface.add_to_nwbfile(nwbfile=nwbfile, metadata=video_metadata, stub_test=False)

        except Exception as e:
            print(f"Warning: Could not process video {video_file}: {e}")
            continue

    # Write NWB file
    with NWBHDF5IO(str(output_path), "w") as io:
        io.write(nwbfile)

    print(f"Created video NWB file: {output_path}")


def main():
    """
    Main function to process all Figure 7 behavioral video data.
    """
    # Set up paths
    video_base_path = Path(
        "/home/heberto/development/surmeier-lab-to-nwb/link_to_raw_data/Figure 7/contralateral rotations/CDGI KO videos"
    )
    output_base_path = Path("/home/heberto/development/surmeier-lab-to-nwb/nwb_files/figure_7_behavioral_videos")

    # Create output directory
    output_base_path.mkdir(parents=True, exist_ok=True)

    # Find all video sessions
    sessions = find_video_sessions(video_base_path)

    print(f"Found {len(sessions)} video sessions")

    # Process each session
    for session_date, session_path in sessions.items():
        print(f"\nProcessing session: {session_date}")

        # Group videos by animal
        animals = {}

        # Search all subsession directories
        for subsession_dir in session_path.iterdir():
            if subsession_dir.is_dir():
                for video_file in subsession_dir.iterdir():
                    if video_file.suffix.lower() in [".mov", ".mp4"]:
                        metadata = extract_video_metadata(video_file)
                        animal_id = metadata["animal_id"]

                        if animal_id not in animals:
                            animals[animal_id] = []
                        animals[animal_id].append(video_file)

        # Create NWB files for each animal
        for animal_id, video_files in animals.items():
            if animal_id != "unknown":
                print(f"  Animal: {animal_id} ({len(video_files)} videos)")

                # Create NWB file
                output_filename = f"figure7_videos_{animal_id}_{session_date.replace('-', '')}.nwb"
                output_path = output_base_path / output_filename

                create_video_nwb_file(
                    session_date=session_date,
                    animal_id=animal_id,
                    video_files=sorted(video_files),
                    output_path=output_path,
                )

    print(f"\nVideo conversion complete. Files saved to: {output_base_path}")


if __name__ == "__main__":
    main()
