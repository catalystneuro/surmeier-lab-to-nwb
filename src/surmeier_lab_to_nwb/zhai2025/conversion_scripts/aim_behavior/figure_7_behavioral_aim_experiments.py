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
from neuroconv.tools.nwb_helpers import make_nwbfile_from_metadata
from neuroconv.utils import dict_deep_update, load_dict_from_file
from pynwb import NWBFile
from pynwb.epoch import TimeIntervals
from tqdm import tqdm

from surmeier_lab_to_nwb.zhai2025.conversion_scripts.conversion_utils import (
    format_condition,
    generate_canonical_session_id,
    str_to_bool,
)
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
    # Load general and session-specific metadata from YAML files
    general_metadata_path = Path(__file__).parent.parent.parent / "general_metadata.yaml"
    general_metadata = load_dict_from_file(general_metadata_path)

    session_metadata_path = Path(__file__).parent.parent.parent / "session_specific_metadata.yaml"
    session_metadata_template = load_dict_from_file(session_metadata_path)
    script_template = session_metadata_template["figure_7_behavioral_aim_experiments"]

    # Parse session information
    session_info = parse_session_info_from_animal_data(session_date, animal_id, genotype)

    # Map genotype to standardized condition for centralized format_condition dictionary
    genotype_to_condition = {"KO": "knockout", "WT": "control", "CDGI KO": "knockout", "Wild-type": "control"}

    # Get standardized condition from genotype
    standardized_condition = genotype_to_condition.get(genotype, "control")

    # Use centralized format_condition dictionary
    condition_camel_case = format_condition[standardized_condition]["CamelCase"]
    condition_human_readable = format_condition[standardized_condition]["human_readable"]

    # Create BIDS-style base session ID with detailed timestamp when available
    session_start_time = session_info["session_start_time"]
    if hasattr(session_start_time, "hour"):
        timestamp = session_start_time.strftime("%Y%m%d_%H%M%S")
    else:
        timestamp = session_start_time.strftime("%Y%m%d")

    # Use centralized format_condition dictionary
    condition_human_readable = format_condition[standardized_condition]["human_readable"]

    # Map genotype to state - AIM scoring is done during ON state (L-DOPA treated)
    # All AIM behavioral experiments are performed after L-DOPA administration
    if genotype == "KO" or genotype == "CDGI KO":
        genotype_canonical = "CDGIKO"
    else:
        genotype_canonical = "WT"

    session_id = generate_canonical_session_id(
        fig="F7",
        meas_comp="AIMs",  # AIM scoring
        cell_type="pan",  # Non cell-specific
        state="ON",  # AIM scoring is during L-DOPA treatment (ON state)
        pharm="none",  # No pharmacology
        geno=genotype_canonical,
        timestamp=timestamp,
    )

    # Create session-specific metadata from template with runtime substitutions
    session_specific_metadata = {
        "NWBFile": {
            "session_description": script_template["NWBFile"]["session_description"].format(
                genotype=genotype, animal_id=session_info["animal_id"], session_number=session_number
            ),
            "session_id": session_id,
            "session_start_time": session_info["session_start_time"],
            "keywords": script_template["NWBFile"]["keywords"],
        },
        "Subject": {
            "subject_id": f"CDGI_mouse_{session_info['original_animal_id']}",
            "genotype": (
                script_template["Subject"]["genotype"] if "KO" in session_info["genotype"] else "Wild-type CDGI"
            ),
            "description": script_template["Subject"]["description"].format(
                genotype=session_info["genotype"], animal_id=session_info["original_animal_id"]
            ),
        },
    }

    metadata = dict_deep_update(general_metadata, session_specific_metadata)

    # Create NWB file using neuroconv helper function
    nwbfile = make_nwbfile_from_metadata(metadata)

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
    parser.add_argument(
        "--stub-test",
        type=str_to_bool,
        default=True,
        help="Process only first 2 files per condition for testing (default: True). Use --stub-test=False for full processing.",
    )
    args = parser.parse_args()

    # Suppress warnings
    logging.getLogger("tifffile").setLevel(logging.ERROR)
    warnings.filterwarnings("ignore", message="invalid value encountered in divide")
    warnings.filterwarnings("ignore", message=".*no datetime before year 1.*")
    warnings.filterwarnings("ignore", message="Unknown extension is not supported and will be removed")

    # Set verbose flag
    verbose = args.verbose
    verbose = False
    stub_test = args.stub_test

    # Set the base path to your data
    base_dir = Path(__file__).parent.parent.parent  # Go up to zhai2025 level
    aim_excel_data_file_path = Path(
        "/media/heberto/One Touch/Surmeier-CN-data-share/consolidated_data/LID_paper_Zhai_2025/Raw data for Figs/Figure 7/AIM rating/AIM testing_CDGI KO.xlsx"
    )
    genotype_path = base_dir / "assets" / "data_connections_D2_figures_3_4_6_7_8.csv"
    processed_data_csv_path = base_dir / "assets" / "processed_aim_behavioral_data" / "figure_7_aim_processed_data.csv"

    # Create nwb_files directory at root level
    root_dir = Path(__file__).parent.parent.parent.parent.parent.parent  # Go up to repo root
    nwb_files_dir = root_dir / "nwb_files" / "aim_behavior" / "figure_7"
    nwb_files_dir.mkdir(parents=True, exist_ok=True)

    # Parse AIM data using unified pipeline
    pivot_df = build_source_data_from_aim_excel_table(aim_excel_data_file_path, genotype_path, verbose=verbose)

    # Create processed data directory and save processed data for interface use
    processed_data_csv_path.parent.mkdir(parents=True, exist_ok=True)
    pivot_df.to_csv(processed_data_csv_path, index=False)

    # Group by session and animal for NWB file creation
    grouped = pivot_df.groupby(["session_date", "session_number", "animal_id", "genotype"])
    sessions_list = list(grouped)

    # Apply stub_test filtering if enabled
    if stub_test:
        sessions_list = sessions_list[:2]

    # Process each session - use tqdm for non-verbose mode
    iterator = tqdm(
        sessions_list,
        desc="Converting Figure7 AIMBehavior",
        unit=" session",
    )

    for (session_date, session_number, animal_id, genotype), session_data in iterator:
        # Use genotype as condition for behavioral experiments
        condition = genotype
        condition_safe = condition.replace(" ", "_").replace("-", "_").replace("(", "").replace(")", "")

        # Convert session data to NWB format
        nwbfile = convert_session_to_nwbfile(
            session_date=session_date,
            session_number=session_number,
            animal_id=animal_id,
            genotype=genotype,
            processed_data_csv_path=processed_data_csv_path,
            verbose=verbose,
        )

        # Create output filename using session_id from nwbfile
        nwbfile_path = nwb_files_dir / f"{nwbfile.session_id}.nwb"

        # Write NWB file
        configure_and_write_nwbfile(nwbfile=nwbfile, nwbfile_path=nwbfile_path)
