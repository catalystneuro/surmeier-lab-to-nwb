"""
Supplementary Figure 3 Behavioral Video Conversion Script

This script processes behavioral video recordings from Supplementary Figure 3 of the Zhai et al. 2025 paper,
specifically the M1R CRISPR behavioral assessment for dyskinesia in gene-edited mice.

The script handles:
- Behavioral video files (.mov/.mp4) from M1R CRISPR experimental sessions
- Video metadata extraction and organization by animal and session
- Integration with NWB format using appropriate video interfaces
- M1R gene editing context and L-DOPA treatment timing
"""

from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List
from zoneinfo import ZoneInfo

from neuroconv.datainterfaces import ExternalVideoInterface
from neuroconv.tools import configure_and_write_nwbfile
from neuroconv.utils import dict_deep_update, load_dict_from_file
from pynwb import NWBFile
from pynwb.epoch import TimeIntervals
from pynwb.file import Subject
from tqdm import tqdm


def find_video_sessions_supfig3(video_base_path: Path) -> Dict[str, Path]:
    """
    Find all video session directories for Supplementary Figure 3.

    Parameters
    ----------
    video_base_path : Path
        Base path to M1R CRISPR video directories

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


def extract_video_metadata_supfig3(video_path: Path) -> Dict[str, Any]:
    """
    Extract metadata from Supplementary Figure 3 video filename.

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

    # Common patterns for M1R CRISPR videos:
    # 8592_1_20min, 4790_1_40min, 9366_2_60min, etc.

    metadata = {
        "animal_id": "unknown",
        "session_number": 1,
        "time_point_min": 0,
        "file_format": video_path.suffix.lower(),
        "genotype": "M1R CRISPR",  # Default for Sup Fig 3, may need refinement
        "treatment": "L-DOPA",
        "genetic_modification": "M1R deletion via CRISPR-Cas9",
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

    # Infer genotype based on experimental design
    # M1R CRISPR vs Control groups in Supplementary Figure 3
    metadata["genotype"] = "M1R CRISPR"  # Most videos are M1R CRISPR, some may be controls

    return metadata


def group_videos_by_animal_supfig3(session_path: Path) -> Dict[str, List[Path]]:
    """
    Group video files by animal ID within a Supplementary Figure 3 session.

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
                    metadata = extract_video_metadata_supfig3(video_file)
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
    genotype: str,
    verbose: bool = False,
) -> NWBFile:
    """
    Convert M1R CRISPR behavioral video recordings from a single experimental session for a single subject to NWB format.

    This function processes video files from a single behavioral session (one animal on one date) containing
    multiple time points of video recordings during L-DOPA treatment in M1R CRISPR experiments. We organize by
    session (animal + date) because:

    1. **Experimental design**: Each session represents a complete L-DOPA behavioral assessment for one animal
       on one experimental date, with videos recorded at multiple time points (20, 40, 60, 80, 100, 120 min)
       to assess the effects of M1R gene editing on dyskinetic behaviors

    2. **Temporal coherence**: All videos within a session share the same experimental context (same animal,
       same M1R CRISPR treatment status, same behavioral protocol) and form a temporally coherent sequence
       following L-DOPA administration

    3. **NWB structure**: This allows proper temporal alignment of all video recordings within a single NWB file,
       enabling synchronized analysis of behavioral changes across the entire M1R CRISPR treatment session

    4. **Gene editing context**: Each session represents the behavioral outcome of a specific genetic manipulation
       (M1R deletion via CRISPR-Cas9), making the session the natural unit for analyzing gene editing effects

    Parameters
    ----------
    session_date : str
        Session date in YYYY-MM-DD format
    animal_id : str
        Animal identifier (e.g., "8592", "4790")
    video_files : List[Path]
        List of video file paths for this animal session
    verbose : bool, default=False
        Whether to print detailed processing information

    Returns
    -------
    NWBFile
        NWB file object containing all video data from the M1R CRISPR session
    """
    if not video_files:
        raise ValueError(f"No video files found for {animal_id} on {session_date}")

    # Parse session start time
    session_start_time = datetime.strptime(session_date, "%Y-%m-%d")
    session_start_time = session_start_time.replace(tzinfo=ZoneInfo("US/Central"))

    # Load general metadata from YAML file
    metadata_file_path = Path(__file__).parent.parent.parent / "metadata.yaml"
    general_metadata = load_dict_from_file(metadata_file_path)

    # Determine genotype based on video metadata
    first_video_metadata = extract_video_metadata_supfig3(video_files[0])
    genotype = first_video_metadata["genotype"]

    # Create session-specific metadata
    conversion_specific_metadata = {
        "NWBFile": {
            "session_description": f"Supplementary Figure 3 Behavioral Videos - M1R CRISPR behavioral assessment for dyskinesia study",
            "identifier": f"zhai2025_supfig3_videos_{animal_id}_{session_date.replace('-', '')}",
            "session_start_time": session_start_time,
            "session_id": f"supfig3_videos_{animal_id}_{session_date.replace('-', '')}",
            "experiment_description": (
                "Behavioral video recordings for dyskinesia assessment in M1R CRISPR gene-edited mice "
                "following L-DOPA treatment. Videos capture behavioral responses at multiple time points "
                "(20, 40, 60, 80, 100, 120 minutes) post-injection to assess the effects of M1R deletion "
                "on dyskinetic behaviors. Data corresponds to Supplementary Figure 3 from Zhai et al. 2025."
            ),
            "keywords": [
                "M1R CRISPR",
                "gene editing",
                "CRISPR-Cas9",
                "muscarinic receptor",
                "iSPNs",
                "behavioral assessment",
                "AIM",
                "Abnormal Involuntary Movement",
            ],
        },
        "Subject": {
            "subject_id": animal_id,
            "description": (
                f"M1R CRISPR study mouse for behavioral video assessment - Supplementary Figure 3, animal ID: {animal_id}. "
                f"M1R deletion achieved via CRISPR-Cas9 gene editing specifically targeting indirect pathway spiny projection neurons (iSPNs). "
                f"Genotype: {genotype}"
            ),
            "genotype": genotype,
        },
    }

    # Deep merge with general metadata
    metadata = dict_deep_update(general_metadata, conversion_specific_metadata)

    # Create NWB file with merged metadata
    nwbfile = NWBFile(
        session_description=metadata["NWBFile"]["session_description"],
        identifier=metadata["NWBFile"]["identifier"],
        session_start_time=metadata["NWBFile"]["session_start_time"],
        experimenter=metadata["NWBFile"]["experimenter"],
        lab=metadata["NWBFile"]["lab"],
        institution=metadata["NWBFile"]["institution"],
        experiment_description=metadata["NWBFile"]["experiment_description"],
        session_id=metadata["NWBFile"]["session_id"],
        surgery=metadata["NWBFile"]["surgery"],
        pharmacology=metadata["NWBFile"]["pharmacology"],
        keywords=metadata["NWBFile"]["keywords"],
    )

    # Add subject information
    nwbfile.subject = Subject(
        subject_id=metadata["Subject"]["subject_id"],
        species=metadata["Subject"]["species"],
        strain=metadata["Subject"]["strain"],
        genotype=metadata["Subject"]["genotype"],
        sex=metadata["Subject"]["sex"],
        age=metadata["Subject"]["age"],
        description=metadata["Subject"]["description"],
    )

    # Sort videos by time point for proper temporal organization
    video_metadata_list = [(vid, extract_video_metadata_supfig3(vid)) for vid in video_files]
    video_metadata_list.sort(key=lambda x: x[1]["time_point_min"])

    # Add time intervals for behavioral epochs
    behavioral_epochs = TimeIntervals(
        name="behavioral_epochs",
        description="Time intervals for behavioral video recordings during L-DOPA treatment in M1R CRISPR study",
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
        name="genotype",
        description="Animal genotype (M1R CRISPR vs Control)",
    )
    behavioral_epochs.add_column(
        name="genetic_modification",
        description="Type of genetic modification applied",
    )
    behavioral_epochs.add_column(
        name="video_file",
        description="Path to the video file",
    )

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
            genotype=vid_metadata["genotype"],
            genetic_modification=vid_metadata["genetic_modification"],
            video_file=str(video_file),
        )

        # Create video interface for each video
        video_interface = ExternalVideoInterface(
            file_paths=[video_file],
            verbose=False,
        )

        # Add video to NWB file with unique name
        video_name = f"video_{animal_id}_{vid_metadata['session_number']}_{vid_metadata['time_point_min']}min"
        video_interface.add_to_nwbfile(
            nwbfile=nwbfile,
            metadata={
                "name": video_name,
                "description": f"M1R CRISPR behavioral video at {vid_metadata['time_point_min']} minutes post L-DOPA treatment for dyskinesia assessment",
            },
        )

    # Add behavioral epochs to NWB file
    nwbfile.add_time_intervals(behavioral_epochs)

    if verbose:
        print(f"Created NWB file for {animal_id} session {session_date} with {len(video_metadata_list)} videos")

    return nwbfile


if __name__ == "__main__":
    import argparse

    # Set up argument parser
    parser = argparse.ArgumentParser(
        description="Convert Supplementary Figure 3 M1R CRISPR behavioral videos to NWB format"
    )
    parser.add_argument("--verbose", "-v", action="store_true", help="Enable verbose output")
    args = parser.parse_args()

    # Set up paths
    video_base_path = Path(
        "/media/heberto/One Touch/Surmeier-CN-data-share/consolidated_data/LID_paper_Zhai_2025/Raw data for Figs/Sup. Figure 3/M1R CRISPR videos"
    )
    output_base_path = Path("./nwb_files/videos/supplementary_figure_3")

    # Create output directory
    output_base_path.mkdir(parents=True, exist_ok=True)

    # Load metadata if available
    metadata_path = Path(__file__).parent / "metadata" / "supplementary_figure_3_behavioral_videos.yaml"
    metadata = None
    if metadata_path.exists():
        metadata = load_dict_from_file(metadata_path)

    # Control test execution
    stub_test = True  # Set to True to process only first 2 files per condition for testing
    verbose = False
    # Find all video sessions
    sessions = find_video_sessions_supfig3(video_base_path)

    # Apply stub_test filtering if enabled
    if stub_test:
        sessions_items = list(sessions.items())[:2]
        sessions = dict(sessions_items)
        if verbose:
            print(f"stub_test enabled: processing only first {len(sessions)} sessions")

    if verbose:
        print(f"Found {len(sessions)} Supplementary Figure 3 video sessions")

    # Process each session
    for session_date, session_path in tqdm(
        sessions.items(), desc="Converting Supplementary Figure3 Videos", unit=" session"
    ):

        # Group videos by animal
        animals = group_videos_by_animal_supfig3(session_path)

        # Create NWB files for each animal
        for animal_id, video_files in animals.items():
            if verbose:
                print(f"  Animal: {animal_id} ({len(video_files)} videos)")

            # Determine genotype based on video metadata
            first_video_metadata = extract_video_metadata_supfig3(video_files[0])
            genotype = first_video_metadata["genotype"]

            # Create output filename
            genotype_safe = genotype.replace(" ", "_").replace("-", "_")
            output_filename = f"zhai2025_supfig3_videos_{animal_id}_{session_date.replace('-', '')}_{genotype_safe}.nwb"
            output_path = output_base_path / output_filename

            # Create NWB file
            nwbfile = convert_session_to_nwbfile(
                session_date=session_date,
                animal_id=animal_id,
                video_files=sorted(video_files),
                genotype=genotype,
                verbose=verbose,
            )

            # Write the NWB file
            configure_and_write_nwbfile(nwbfile=nwbfile, nwbfile_path=output_path)
            if verbose:
                print(f"    Saved: {output_path}")

    if verbose:
        print(f"\nSupplementary Figure 3 video conversion complete. Files saved to: {output_base_path}")
