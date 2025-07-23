"""
Figure 7 Behavioral Data Conversion Script

This script processes the behavioral data from Figure 7 of the Zhai et al. 2025 paper,
including AIM (Abnormal Involuntary Movement) scoring for CDGI knockout mice.

The script handles:
- AIM scoring data from Excel files with multiple sessions using behavioral_utils
- Genotype information (CDGI KO vs WT)
- Optimized DynamicTable structure for Figure 7J reproduction
"""

import logging
import warnings
from datetime import datetime
from pathlib import Path
from typing import Any, Dict
from zoneinfo import ZoneInfo

from neuroconv.tools import configure_and_write_nwbfile
from neuroconv.utils import dict_deep_update, load_dict_from_file
from pynwb import NWBFile
from pynwb.epoch import TimeIntervals
from pynwb.file import Subject
from tqdm import tqdm

from surmeier_lab_to_nwb.zhai2025.interfaces.behavior_interfaces import (
    AIMBehavioralDynamicTableInterface,
    AIMBehavioralTimeSeriesInterface,
    build_source_data_from_aim_excel_table,
)


def parse_session_info_from_animal_data(session_date: str, animal_id: int, genotype: str) -> Dict[str, Any]:
    """
    Parse essential session information from AIM scoring data.

    Parameters
    ----------
    session_date : str
        Session date in YYYY-MM-DD format
    animal_id : int
        Animal identifier
    genotype : str
        Animal genotype (KO, WT, or unknown)

    Returns
    -------
    Dict[str, Any]
        Dictionary containing session information
    """
    # Parse date
    date_obj = datetime.strptime(session_date, "%Y-%m-%d")

    return {
        "animal_id": str(animal_id),
        "genotype": genotype,
        "session_date": session_date,
        "session_start_time": date_obj.replace(tzinfo=ZoneInfo("US/Central")),
        "original_animal_id": f"ET#{animal_id}",
    }


def convert_session_to_nwbfile(
    session_date: str,
    session_number: int,
    animal_id: int,
    genotype: str,
    processed_data_csv_path: Path,
    verbose: bool = False,
) -> NWBFile:
    """
    Convert AIM scoring data from a single session to NWB format using DataInterface.

    Parameters
    ----------
    session_date : str
        Session date in YYYY-MM-DD format
    session_number : int
        Session number within the day (1, 2, 3, 4, 5)
    animal_id : int
        Animal identifier
    genotype : str
        Animal genotype
    processed_data_csv_path : Path
        Path to the processed CSV file containing AIM data
    verbose : bool, default=False
        Enable verbose output

    Returns
    -------
    NWBFile
        NWB file with the converted data
    """
    # Load general metadata
    metadata_path = Path(__file__).parent / "metadata.yaml"
    metadata = load_dict_from_file(metadata_path)

    # Parse session information
    session_info = parse_session_info_from_animal_data(session_date, animal_id, genotype)

    # Create BIDS-style base session ID with detailed timestamp when available
    session_start_time = session_info["session_start_time"]
    if hasattr(session_start_time, "hour"):
        timestamp = session_start_time.strftime("%Y%m%d_%H%M%S")
    else:
        timestamp = session_start_time.strftime("%Y%m%d")

    base_session_id = f"figure7_BehavioralAIM_{genotype.replace(' ', '_').replace('-', '_')}_{timestamp}"
    script_specific_id = f"Sub{session_info['animal_id']}_Session{session_number}"
    session_id = f"{base_session_id}_{script_specific_id}"

    # Update metadata for behavioral experiment
    session_specific_metadata = {
        "NWBFile": {
            "session_description": "Figure 7 Behavioral Assessment - AIM scoring for dyskinesia analysis in CDGI study",
            "identifier": f"figure7_behavioral_{session_info['animal_id']}_{session_date.replace('-', '')}_s{session_number}",
            "session_id": session_id,
            "session_start_time": session_info["session_start_time"],
            "keywords": ["AIM", "Abnormal Involuntary Movement", "CDGI", "dyskinesia", "L-DOPA"],
            "experiment_description": (
                "Behavioral assessment of CDGI knockout mice using AIM scoring to evaluate "
                "dyskinesia severity following L-DOPA treatment. Data corresponds to Figure 7J "
                "from Zhai et al. 2025. Uses optimized DynamicTable for figure reproduction."
            ),
        },
        "Subject": {
            "subject_id": session_info["original_animal_id"],
            "genotype": session_info["genotype"],
            "description": f"CDGI study animal - {session_info['genotype']}",
            "species": "Mus musculus",
            "strain": "C57BL/6J",
            "sex": "M",
            "age": "P7W/P12W",  # 7-12 weeks old
            "genotype_description": (
                "Hemizygous for BAC transgene (Drd1a-tdTomato or Drd2-eGFP reporter) "
                "back-crossed to C57BL/6 background. "
                + (
                    "Crossed with CDGI-null line maintained on C57BL/6J background."
                    if "KO" in session_info["genotype"]
                    else "Wild-type CDGI."
                )
            ),
        },
    }

    metadata = dict_deep_update(metadata, session_specific_metadata)

    # Create NWBFile
    nwbfile = NWBFile(
        session_description=metadata["NWBFile"]["session_description"],
        identifier=metadata["NWBFile"]["identifier"],
        session_start_time=metadata["NWBFile"]["session_start_time"],
        experimenter=metadata["NWBFile"]["experimenter"],
        lab=metadata["NWBFile"]["lab"],
        institution=metadata["NWBFile"]["institution"],
        experiment_description=metadata["NWBFile"]["experiment_description"],
        session_id=metadata["NWBFile"]["session_id"],
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

    if verbose:
        print(f"  Processing behavioral data for {session_info['original_animal_id']} ({session_info['genotype']})")

    # Create and use AIM Behavioral DynamicTable Interface
    aim_interface = AIMBehavioralDynamicTableInterface(
        processed_data_csv_path=processed_data_csv_path,
        session_date=session_date,
        session_number=session_number,
        animal_id=animal_id,
        verbose=verbose,
    )

    # Add behavioral data to NWB file using the DynamicTable interface
    aim_interface.add_to_nwbfile(nwbfile)

    # Also add behavioral data as BehavioralTimeSeries for time-aligned analysis
    timeseries_interface = AIMBehavioralTimeSeriesInterface(
        processed_data_csv_path=processed_data_csv_path,
        session_date=session_date,
        session_number=session_number,
        animal_id=animal_id,
        verbose=verbose,
    )
    timeseries_interface.add_to_nwbfile(nwbfile)

    # Add experimental epochs for behavioral assessment time intervals
    epochs = TimeIntervals(
        name="behavioral_epochs",
        description="Time intervals for behavioral assessment following L-DOPA administration",
    )

    for time_point in aim_interface.time_points:
        start_time = time_point * 60.0 - 30.0  # 30 seconds before scoring
        stop_time = time_point * 60.0 + 30.0  # 30 seconds after scoring

        epochs.add_interval(
            start_time=max(0, start_time), stop_time=stop_time, tags=["AIM_scoring", f"T{time_point}min"]
        )

    nwbfile.add_time_intervals(epochs)

    if verbose:
        print(f"  Successfully processed behavioral data for {session_info['original_animal_id']}")

    return nwbfile


if __name__ == "__main__":
    import argparse

    # Set up argument parser
    parser = argparse.ArgumentParser(
        description="Convert AIM behavioral data to NWB format with optimized DynamicTable structure"
    )
    parser.add_argument(
        "--verbose", "-v", action="store_true", help="Enable verbose output (show detailed processing for each file)"
    )
    args = parser.parse_args()

    # Suppress warnings
    logging.getLogger("tifffile").setLevel(logging.ERROR)
    warnings.filterwarnings("ignore", message="invalid value encountered in divide")
    warnings.filterwarnings("ignore", message=".*no datetime before year 1.*")
    warnings.filterwarnings("ignore", message="Unknown extension is not supported and will be removed")

    # Set verbose flag
    verbose = args.verbose

    # Set the base path to your data
    base_dir = Path(__file__).parent
    aim_excel_data_file_path = Path(
        "/media/heberto/One Touch/Surmeier-CN-data-share/consolidated_data/LID_paper_Zhai_2025/Raw data for Figs/Figure 7/AIM rating/AIM testing_CDGI KO.xlsx"
    )
    genotype_path = base_dir / "assets" / "data_connections_D2_figures_3_4_6_7_8.csv"
    processed_data_csv_path = base_dir / "assets" / "processed_aim_behavioral_data" / "figure_7_aim_processed_data.csv"

    # Create nwb_files directory at root level
    root_dir = Path(__file__).parent.parent.parent.parent  # Go up to repo root
    nwb_files_dir = root_dir / "nwb_files" / "figure_7" / "behavioral"
    nwb_files_dir.mkdir(parents=True, exist_ok=True)

    if verbose:
        print("Processing Figure 7 Behavioral Assessment data with optimized DynamicTable")

    if verbose:
        print("Verbose mode: ON - Detailed processing information will be shown")
        print("Parsing AIM Excel data...")

    # Parse AIM data using unified pipeline
    pivot_df = build_source_data_from_aim_excel_table(aim_excel_data_file_path, genotype_path, verbose=verbose)

    # Create processed data directory and save processed data for interface use
    processed_data_csv_path.parent.mkdir(parents=True, exist_ok=True)
    pivot_df.to_csv(processed_data_csv_path, index=False)
    if verbose:
        print(f"Saved processed data to {processed_data_csv_path}")

    if verbose:
        print(f"Built source data with {len(pivot_df)} rows ready for NWB")

    # Group by session and animal for NWB file creation
    grouped = pivot_df.groupby(["session_date", "session_number", "animal_id", "genotype"])
    sessions_list = list(grouped)

    if verbose:
        print(f"Found {len(grouped)} unique sessions to process")
        print("\nDetailed processing log:")
        print("=" * 60)

    # Process each session - use tqdm for non-verbose mode
    iterator = (
        sessions_list
        if verbose
        else tqdm(
            sessions_list,
            desc="Converting sessions from figure_7_behavioral_aim_experiments to NWB",
        )
    )

    for (session_date, session_number, animal_id, genotype), session_data in iterator:
        # Sanitize condition string for NWB identifier and session_id
        condition_safe = condition.replace(" ", "_").replace("-", "_").replace("(", "").replace(")", "")

        if verbose:
            print(f"\nProcessing: Animal {animal_id} ({genotype}) on {session_date} session {session_number}")
            print(f"  Session data shape: {session_data.shape}")
            print(f"  Time points: {sorted(session_data['time_minutes'].unique())}")
            print(f"  Score columns: {[col for col in session_data.columns if 'score' in col.lower()]}")

        # Convert session data to NWB format
        nwbfile = convert_session_to_nwbfile(
            session_date=session_date,
            session_number=session_number,
            animal_id=animal_id,
            genotype=genotype,
            processed_data_csv_path=processed_data_csv_path,
            condition=condition_safe,
            verbose=verbose,
        )

        # Create output filename
        nwbfile_path = (
            nwb_files_dir
            / f"figure7_behavioral_{animal_id}_{session_date.replace('-', '')}_s{session_number}_{condition_safe}.nwb"
        )

        # Write NWB file
        configure_and_write_nwbfile(nwbfile=nwbfile, nwbfile_path=nwbfile_path)

        if verbose:
            print(f"  Successfully saved: {nwbfile_path.name}")
            print(f"  File size: {nwbfile_path.stat().st_size / 1024:.1f} KB")
