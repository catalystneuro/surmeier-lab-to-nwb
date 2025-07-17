"""
Figure 7 Behavioral Data Conversion Script

This script processes the behavioral data from Figure 7 of the Zhai et al. 2025 paper,
including AIM (Abnormal Involuntary Movement) scoring for CDGI knockout mice.

The script handles:
- AIM scoring data from Excel files with multiple sessions
- Genotype information (CDGI KO vs WT)
- Time-series behavioral data with L-DOPA response curves
"""

import logging
import warnings
from datetime import datetime
from pathlib import Path
from typing import Any, Dict
from zoneinfo import ZoneInfo

import numpy as np
import pandas as pd
from neuroconv.tools import configure_and_write_nwbfile
from neuroconv.utils import dict_deep_update, load_dict_from_file
from pynwb import NWBFile, TimeSeries
from pynwb.behavior import BehavioralTimeSeries
from pynwb.epoch import TimeIntervals
from pynwb.file import Subject

# Set verbose flag
VERBOSE = False


def parse_session_info_from_excel_data(session_date: str, animal_id: str, genotype: str) -> Dict[str, Any]:
    """
    Parse essential session information from AIM scoring data.

    Parameters
    ----------
    session_date : str
        Session date in YYYY-MM-DD format
    animal_id : str
        Animal identifier (e.g., 'ET#1944')
    genotype : str
        Animal genotype (CDGI_KO, WT, or unknown)

    Returns
    -------
    Dict[str, Any]
        Dictionary containing session information
    """
    # Clean animal ID
    clean_animal_id = animal_id.replace("ET#", "").replace(" ", "").strip()

    # Parse date
    date_obj = datetime.strptime(session_date, "%Y-%m-%d")

    return {
        "animal_id": clean_animal_id,
        "genotype": genotype,
        "session_date": session_date,
        "session_start_time": date_obj.replace(tzinfo=ZoneInfo("US/Central")),
        "original_animal_id": animal_id,
    }


def parse_aim_excel_data(excel_path: Path) -> Dict[str, Any]:
    """
    Parse AIM scoring data from Excel file.

    Parameters
    ----------
    excel_path : Path
        Path to the AIM testing Excel file

    Returns
    -------
    Dict[str, Any]
        Dictionary containing parsed AIM data organized by session and animal
    """
    excel_file = pd.ExcelFile(excel_path)
    aim_data = {}

    # Time points for AIM scoring (in minutes post-L-DOPA)
    time_points = [20, 40, 60, 80, 100, 120]

    for sheet_name in excel_file.sheet_names:
        # Parse date from sheet name (format: MMDDYYYY)
        if len(sheet_name) == 8:
            month = sheet_name[:2]
            day = sheet_name[2:4]
            year = sheet_name[4:8]
            session_date = f"{year}-{month}-{day}"
        else:
            continue

        df = pd.read_excel(excel_path, sheet_name=sheet_name)

        # Extract animal data from the sheet
        animals = {}
        i = 0

        while i < len(df):
            row = df.iloc[i]

            # Look for animal ID (ET# pattern)
            if pd.notna(row.iloc[0]) and "ET#" in str(row.iloc[0]):
                animal_id = str(row.iloc[0])

                # Extract genotype information
                genotype = "unknown"
                if "(KO)" in animal_id:
                    genotype = "CDGI_KO"
                elif "(WT)" in animal_id:
                    genotype = "WT"

                # Clean animal ID
                animal_id = animal_id.replace("(KO)", "").replace("(WT)", "").strip()

                # Extract scoring data for this animal
                if i + 3 < len(df):
                    # Check if this is AIM scoring data (look for 'Ax', 'Li', 'Ol' in next rows)
                    if "Ax" in str(df.iloc[i + 1].iloc[1]) or "Ax" in str(df.iloc[i].iloc[1]):
                        # Determine if categories are in same row or next rows
                        if "Ax" in str(df.iloc[i].iloc[1]):
                            ax_row = i
                            li_row = i + 1
                            ol_row = i + 2
                            total_row = i + 3
                        else:
                            ax_row = i + 1
                            li_row = i + 2
                            ol_row = i + 3
                            total_row = i + 4

                        # Helper function to safely convert values to float
                        def safe_float(val):
                            if pd.notna(val) and str(val).strip() != "":
                                return float(val)
                            return np.nan

                        # Extract scores for each category
                        ax_scores = []
                        li_scores = []
                        ol_scores = []
                        total_scores = []

                        # Get scores from columns 2-7 (time points)
                        for col_idx in range(2, min(8, len(df.columns))):
                            if ax_row < len(df):
                                ax_val = df.iloc[ax_row, col_idx]
                                ax_scores.append(safe_float(ax_val))

                            if li_row < len(df):
                                li_val = df.iloc[li_row, col_idx]
                                li_scores.append(safe_float(li_val))

                            if ol_row < len(df):
                                ol_val = df.iloc[ol_row, col_idx]
                                ol_scores.append(safe_float(ol_val))

                            if total_row < len(df):
                                total_val = df.iloc[total_row, col_idx]
                                total_scores.append(safe_float(total_val))

                        # Look for rotation data in column 9 (index 9)
                        rotation_count = np.nan
                        if len(df.columns) > 9 and ax_row < len(df):
                            rot_val = df.iloc[ax_row, 9]
                            if pd.notna(rot_val):
                                rotation_count = float(rot_val)

                        animals[animal_id] = {
                            "genotype": genotype,
                            "aim_scores": {
                                "axial": ax_scores,
                                "limb": li_scores,
                                "orolingual": ol_scores,
                                "total": total_scores,
                            },
                            "rotation_count": rotation_count,
                            "time_points": time_points[: len(ax_scores)],
                        }

                        i = total_row + 1
                    else:
                        i += 1
                else:
                    i += 1
            else:
                i += 1

        aim_data[session_date] = animals

    return aim_data


def convert_session_to_nwbfile(
    session_date: str, animal_id: str, aim_data: Dict[str, Any], verbose: bool = False
) -> NWBFile:
    """
    Convert AIM scoring data from a single session to NWB format.

    Parameters
    ----------
    session_date : str
        Session date in YYYY-MM-DD format
    animal_id : str
        Animal identifier
    aim_data : Dict[str, Any]
        AIM scoring data for this animal
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
    session_info = parse_session_info_from_excel_data(session_date, animal_id, aim_data["genotype"])

    # Update metadata for behavioral experiment
    behavioral_metadata = dict_deep_update(
        metadata,
        {
            "NWBFile": {
                "session_description": "Figure 7 Behavioral Assessment - AIM scoring for dyskinesia analysis in CDGI study",
                "identifier": f"figure7_behavioral_{session_info['animal_id']}_{session_date.replace('-', '')}",
                "session_start_time": session_info["session_start_time"],
                "experiment_description": "Behavioral assessment of CDGI knockout mice using AIM scoring to evaluate dyskinesia severity following L-DOPA treatment. Data corresponds to Figure 7J from Zhai et al. 2025.",
            },
            "Subject": {
                "subject_id": session_info["original_animal_id"],
                "genotype": session_info["genotype"],
                "description": f"CDGI study animal - {session_info['genotype']}",
            },
        },
    )

    # Create NWBFile
    nwbfile = NWBFile(**behavioral_metadata["NWBFile"])

    # Add subject information
    nwbfile.subject = Subject(**behavioral_metadata["Subject"])

    if verbose:
        print(f"  Processing behavioral data for {session_info['original_animal_id']} ({session_info['genotype']})")

    # Add behavioral time series data
    if "aim_scores" in aim_data:
        aim_scores = aim_data["aim_scores"]
        time_points = aim_data["time_points"]

        # Convert time points to seconds (minutes * 60)
        timestamps = np.array(time_points) * 60.0

        # Create behavioral time series for each AIM category
        behavioral_module = nwbfile.create_processing_module(
            name="behavior", description="Behavioral data including AIM scoring and rotation analysis"
        )

        # Create AIM scoring time series
        aim_categories = ["axial", "limb", "orolingual", "total"]
        aim_ts_data = []

        for category in aim_categories:
            if category in aim_scores:
                scores = np.array(aim_scores[category])
                # Replace NaN values with -1 to indicate missing data
                scores = np.nan_to_num(scores, nan=-1.0)
                aim_ts_data.append(scores)

        if aim_ts_data:
            aim_ts_data = np.column_stack(aim_ts_data)

            aim_time_series = TimeSeries(
                name="AIM_scores",
                data=aim_ts_data,
                timestamps=timestamps,
                unit="score",
                description="Abnormal Involuntary Movement (AIM) scores across different behavioral categories. "
                "Scores range from 0-4 where 0=no abnormal movements, 1=occasional (<50%), "
                "2=frequent (>50%), 3=continuous with interruptions, 4=continuous without interruptions. "
                "Missing values are indicated by -1.",
                comments="AIM scoring performed at regular intervals following L-DOPA administration. "
                "Columns: axial, limb, orolingual, total",
            )

            behavioral_ts = BehavioralTimeSeries(name="AIM_scoring")
            behavioral_ts.add_timeseries(aim_time_series)
            behavioral_module.add(behavioral_ts)

        # Add rotation data if available
        if "rotation_count" in aim_data and not np.isnan(aim_data["rotation_count"]):
            rotation_ts = TimeSeries(
                name="contralateral_rotations",
                data=[aim_data["rotation_count"]],
                timestamps=[timestamps[-1]],  # Use final time point
                unit="rotations/min",
                description="Contralateral rotation count per minute, measured during dyskinesia assessment",
                comments="Quantified from video analysis of rotational behavior",
            )

            rotation_behavioral = BehavioralTimeSeries(name="rotation_analysis")
            rotation_behavioral.add_timeseries(rotation_ts)
            behavioral_module.add(rotation_behavioral)

    # Add experimental epochs
    if "time_points" in aim_data:
        epochs = TimeIntervals(
            name="behavioral_epochs",
            description="Time intervals for behavioral assessment following L-DOPA administration",
        )

        for time_point in aim_data["time_points"]:
            start_time = time_point * 60.0 - 30.0  # 30 seconds before scoring
            stop_time = time_point * 60.0 + 30.0  # 30 seconds after scoring

            epochs.add_interval(
                start_time=max(0, start_time), stop_time=stop_time, tags=[f"AIM_scoring", f"T{time_point}min"]
            )

        nwbfile.add_time_intervals(epochs)

    if verbose:
        print(f"  Successfully processed behavioral data for {session_info['original_animal_id']}")

    return nwbfile


if __name__ == "__main__":
    # Suppress warnings
    logging.getLogger("tifffile").setLevel(logging.ERROR)
    warnings.filterwarnings("ignore", message="invalid value encountered in divide")
    warnings.filterwarnings("ignore", message=".*no datetime before year 1.*")

    # Set the base path to your data
    excel_path = Path("./link_to_raw_data/Figure 7/AIM rating/AIM testing_CDGI KO.xlsx")

    # Create nwb_files directory at root level
    root_dir = Path(__file__).parent.parent.parent.parent  # Go up to repo root
    nwb_files_dir = root_dir / "nwb_files" / "figure_7_behavioral"
    nwb_files_dir.mkdir(parents=True, exist_ok=True)

    print("Processing Figure 7 Behavioral Assessment data")

    # Parse AIM data
    print("Parsing AIM Excel data...")
    aim_data = parse_aim_excel_data(excel_path)

    print(f"Found {len(aim_data)} experimental sessions")

    # Process each session and animal
    for session_date, animals in aim_data.items():
        print(f"\nProcessing session: {session_date}")
        print(f"Found {len(animals)} animals")

        for animal_id, animal_data in animals.items():
            print(f"\nProcessing animal: {animal_id}")

            # Convert session data to NWB format
            nwbfile = convert_session_to_nwbfile(
                session_date=session_date,
                animal_id=animal_id,
                aim_data=animal_data,
                verbose=VERBOSE,
            )

            # Create output filename
            clean_animal_id = animal_id.replace("ET#", "").replace(" ", "_").strip()
            nwbfile_path = nwb_files_dir / f"figure7_behavioral_{clean_animal_id}_{session_date.replace('-', '')}.nwb"

            # Write NWB file
            configure_and_write_nwbfile(nwbfile, nwbfile_path=nwbfile_path)
            print(f"Successfully saved: {nwbfile_path.name}")

    print(f"\nBehavioral data conversion complete. Files saved to: {nwb_files_dir}")
