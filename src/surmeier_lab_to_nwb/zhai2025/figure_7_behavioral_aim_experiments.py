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

import numpy as np
import pandas as pd
from behavioral_utils import parse_aim_excel_to_wide_format
from neuroconv.tools import configure_and_write_nwbfile
from neuroconv.utils import dict_deep_update, load_dict_from_file
from pynwb import NWBFile
from pynwb.core import DynamicTable
from pynwb.epoch import TimeIntervals
from pynwb.file import Subject
from tqdm import tqdm


def transform_wide_to_long_format(wide_df: pd.DataFrame) -> pd.DataFrame:
    """
    Transform wide-format AIM data to long format optimized for NWB DynamicTable.

    This transformation converts data from wide format (one row per animal-session-score_type)
    to long format (one row per animal-session-time_point) to optimize for the DynamicTable
    structure that enables trivial Figure 7J reproduction.

    Parameters
    ----------
    wide_df : pd.DataFrame
        Wide-format DataFrame from behavioral_utils parser with columns:
        ['sheet_name', 'session_date', 'animal_id', 'genotype', 'score_type',
         '20min', '40min', '60min', '80min', '100min', '120min']

    Returns
    -------
    pd.DataFrame
        Long-format DataFrame with one row per time point per animal with columns:
        ['animal_id', 'session_date', 'genotype', 'time_minutes', 'score_type', 'score_value']

    Examples
    --------
    Input (wide format):
    ```
    | animal_id | genotype | score_type | 20min | 40min | 60min |
    |-----------|----------|------------|-------|-------|-------|
    | 1944      | KO       | axial      | 3.5   | 4.0   | 4.0   |
    | 1944      | KO       | limb       | 1.5   | 2.0   | 2.0   |
    | 1944      | KO       | orolingual | 1.5   | 1.5   | 2.0   |
    ```

    Output (long format):
    ```
    | animal_id | genotype | time_minutes | score_type | score_value |
    |-----------|----------|--------------|------------|-------------|
    | 1944      | KO       | 20           | axial      | 3.5         |
    | 1944      | KO       | 40           | axial      | 4.0         |
    | 1944      | KO       | 60           | axial      | 4.0         |
    | 1944      | KO       | 20           | limb       | 1.5         |
    | 1944      | KO       | 40           | limb       | 2.0         |
    | 1944      | KO       | 60           | limb       | 2.0         |
    | 1944      | KO       | 20           | orolingual | 1.5         |
    | 1944      | KO       | 40           | orolingual | 1.5         |
    | 1944      | KO       | 60           | orolingual | 2.0         |
    ```

    Notes
    -----
    This transformation enables the DynamicTable to pivot the data back into the
    optimal structure where all score components (axial, limb, orolingual) are in
    the same row with the same time point, making Figure 7J reproduction trivial:

    Final DynamicTable structure:
    ```
    | animal_id | time_minutes | axial_score | limb_score | orolingual_score | total_score |
    |-----------|--------------|-------------|------------|------------------|-------------|
    | 1944      | 20           | 3.5         | 1.5        | 1.5              | 6.5         |
    | 1944      | 40           | 4.0         | 2.0        | 1.5              | 7.5         |
    | 1944      | 60           | 4.0         | 2.0        | 2.0              | 8.0         |
    ```

    This structure allows Figure 7J reproduction with simple operations:
    - `df.groupby(['time_minutes', 'genotype']).mean()` for plotting
    - `df[['axial_score', 'limb_score', 'orolingual_score']].corr()` for correlations
    - No complex timestamp alignment or data merging required
    """
    long_data = []

    # Time point columns in the wide format
    time_cols = ["20min", "40min", "60min", "80min", "100min", "120min"]

    for _, row in wide_df.iterrows():
        animal_id = row["animal_id"]
        session_date = row["session_date"]
        genotype = row["genotype"]
        score_type = row["score_type"]

        for time_col in time_cols:
            if pd.notna(row[time_col]):  # Only include non-NaN values
                time_minutes = int(time_col.replace("min", ""))

                long_data.append(
                    {
                        "animal_id": animal_id,
                        "session_date": session_date,
                        "genotype": genotype,
                        "time_minutes": time_minutes,
                        "score_type": score_type,
                        "score_value": row[time_col],
                    }
                )

    return pd.DataFrame(long_data)


def create_aim_score_table(long_df: pd.DataFrame) -> DynamicTable:
    """
    Create optimized DynamicTable for AIM scores for Figure 7J reproduction.

    Parameters
    ----------
    long_df : pd.DataFrame
        Long-format DataFrame with AIM scores

    Returns
    -------
    DynamicTable
        NWB DynamicTable optimized for figure reproduction
    """
    # Pivot to get axial, limb, orolingual in separate columns
    pivot_df = long_df.pivot_table(
        index=["animal_id", "session_date", "genotype", "time_minutes"],
        columns="score_type",
        values="score_value",
        aggfunc="first",
    ).reset_index()

    # Fill missing score types with NaN
    for score_type in ["axial", "limb", "orolingual"]:
        if score_type not in pivot_df.columns:
            pivot_df[score_type] = np.nan

    # Calculate total score
    pivot_df["total_score"] = pivot_df[["axial", "limb", "orolingual"]].sum(axis=1, min_count=1)

    # Convert time to timestamps (midpoint of 20-min bins in seconds)
    pivot_df["timestamps"] = pivot_df["time_minutes"] * 60.0

    # Create DynamicTable with explicit data - build row by row
    aim_table = DynamicTable(
        name="AIMScoreTable",
        description="Abnormal Involuntary Movement scores optimized for Figure 7J reproduction. "
        "Each row represents one time point for one animal. Scores range 0-4 where "
        "0=no abnormal movements, 4=continuous without interruptions.",
    )

    # Add columns to the table first (empty)
    aim_table.add_column(name="animal_id", description="Animal identifier")
    aim_table.add_column(name="genotype", description="Animal genotype (CDGI Knockout/Wild Type)")
    aim_table.add_column(name="session_date", description="Session date (YYYY-MM-DD)")
    aim_table.add_column(name="time_minutes", description="Time post L-DOPA (minutes)")
    aim_table.add_column(name="timestamps", description="Time post L-DOPA (seconds)")
    aim_table.add_column(name="axial_score", description="Axial dyskinesia score (0-4)")
    aim_table.add_column(name="limb_score", description="Limb dyskinesia score (0-4)")
    aim_table.add_column(name="orolingual_score", description="Orolingual dyskinesia score (0-4)")
    aim_table.add_column(name="total_score", description="Total AIM score (sum of components)")

    # Add rows
    for _, row in pivot_df.iterrows():
        aim_table.add_row(
            animal_id=int(row["animal_id"]),
            genotype=row["genotype"],
            session_date=row["session_date"],
            time_minutes=int(row["time_minutes"]),
            timestamps=float(row["timestamps"]),
            axial_score=float(row["axial"]) if pd.notna(row["axial"]) else np.nan,
            limb_score=float(row["limb"]) if pd.notna(row["limb"]) else np.nan,
            orolingual_score=float(row["orolingual"]) if pd.notna(row["orolingual"]) else np.nan,
            total_score=float(row["total_score"]) if pd.notna(row["total_score"]) else np.nan,
        )

    return aim_table


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


def parse_aim_data_with_behavioral_utils(excel_path: Path, genotype_csv_path: Path = None) -> pd.DataFrame:
    """
    Parse AIM scoring data using the robust behavioral_utils parser.

    Parameters
    ----------
    excel_path : Path
        Path to the AIM testing Excel file
    genotype_csv_path : Path, optional
        Path to the genotype CSV file

    Returns
    -------
    pd.DataFrame
        Wide-format DataFrame with AIM scores
    """
    return parse_aim_excel_to_wide_format(excel_path=excel_path, genotype_csv_path=genotype_csv_path)


def convert_session_to_nwbfile(
    session_date: str, animal_id: int, genotype: str, session_data: pd.DataFrame, verbose: bool = False
) -> NWBFile:
    """
    Convert AIM scoring data from a single session to NWB format with optimized DynamicTable.

    Parameters
    ----------
    session_date : str
        Session date in YYYY-MM-DD format
    animal_id : int
        Animal identifier
    genotype : str
        Animal genotype
    session_data : pd.DataFrame
        Long-format AIM data for this session/animal
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

    # Update metadata for behavioral experiment
    behavioral_metadata = dict_deep_update(
        metadata,
        {
            "NWBFile": {
                "session_description": "Figure 7 Behavioral Assessment - AIM scoring for dyskinesia analysis in CDGI study",
                "identifier": f"figure7_behavioral_{session_info['animal_id']}_{session_date.replace('-', '')}",
                "session_start_time": session_info["session_start_time"],
                "experiment_description": "Behavioral assessment of CDGI knockout mice using AIM scoring to evaluate dyskinesia severity following L-DOPA treatment. Data corresponds to Figure 7J from Zhai et al. 2025. Uses optimized DynamicTable for figure reproduction.",
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

    # Create behavioral processing module
    behavioral_module = nwbfile.create_processing_module(
        name="behavior", description="Behavioral data including AIM scoring optimized for Figure 7J reproduction"
    )

    # Create optimized AIM score table
    if not session_data.empty:
        aim_table = create_aim_score_table(session_data)
        behavioral_module.add(aim_table)

        # Add experimental epochs
        unique_times = session_data["time_minutes"].unique()
        if len(unique_times) > 0:
            epochs = TimeIntervals(
                name="behavioral_epochs",
                description="Time intervals for behavioral assessment following L-DOPA administration",
            )

            for time_point in unique_times:
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
    excel_path = base_dir / "assets" / "AIM testing_CDGI KO.xlsx"
    genotype_path = base_dir / "assets" / "data_connections_D2_figures_3_4_6_7_8.csv"

    # Create nwb_files directory at root level
    root_dir = Path(__file__).parent.parent.parent.parent  # Go up to repo root
    nwb_files_dir = root_dir / "nwb_files" / "figure_7_behavioral"
    nwb_files_dir.mkdir(parents=True, exist_ok=True)

    print("Processing Figure 7 Behavioral Assessment data with optimized DynamicTable")
    if verbose:
        print("Verbose mode: ON - Detailed processing information will be shown")
    else:
        print("Verbose mode: OFF - Using progress bars for file processing")

    # Parse AIM data using behavioral_utils
    print("Parsing AIM Excel data with behavioral_utils...")
    wide_df = parse_aim_data_with_behavioral_utils(excel_path, genotype_path)

    print(f"Parsed {len(wide_df)} rows of AIM data")

    # Transform to long format
    print("Transforming to long format for NWB optimization...")
    long_df = transform_wide_to_long_format(wide_df)

    print(f"Transformed to {len(long_df)} long-format rows")

    # Group by session and animal for NWB file creation
    grouped = long_df.groupby(["session_date", "animal_id", "genotype"])

    print(f"Found {len(grouped)} unique sessions to process")

    # Process each session and animal
    sessions_list = list(grouped)

    if verbose:
        # Verbose mode: show detailed processing for each file
        print("\nDetailed processing log:")
        print("=" * 60)
        for (session_date, animal_id, genotype), session_data in sessions_list:
            print(f"\nProcessing: Animal {animal_id} ({genotype}) on {session_date}")
            print(f"  Session data shape: {session_data.shape}")
            print(f"  Time points: {sorted(session_data['time_minutes'].unique())}")
            print(f"  Score types: {sorted(session_data['score_type'].unique())}")

            # Convert session data to NWB format
            nwbfile = convert_session_to_nwbfile(
                session_date=session_date,
                animal_id=animal_id,
                genotype=genotype,
                session_data=session_data,
                verbose=verbose,
            )

            # Create output filename
            nwbfile_path = nwb_files_dir / f"figure7_behavioral_{animal_id}_{session_date.replace('-', '')}.nwb"

            # Write NWB file
            configure_and_write_nwbfile(nwbfile, nwbfile_path=nwbfile_path)
            print(f"  ✅ Successfully saved: {nwbfile_path.name}")
            print(f"  File size: {nwbfile_path.stat().st_size / 1024:.1f} KB")
    else:
        # Non-verbose mode: use tqdm progress bar
        print("\nProcessing sessions with progress bar:")
        for (session_date, animal_id, genotype), session_data in tqdm(
            sessions_list, desc="Converting to NWB", unit="session", ncols=80
        ):
            # Convert session data to NWB format (quietly)
            nwbfile = convert_session_to_nwbfile(
                session_date=session_date,
                animal_id=animal_id,
                genotype=genotype,
                session_data=session_data,
                verbose=False,  # Always quiet in non-verbose mode
            )

            # Create output filename
            nwbfile_path = nwb_files_dir / f"figure7_behavioral_{animal_id}_{session_date.replace('-', '')}.nwb"

            # Write NWB file (quietly)
            configure_and_write_nwbfile(nwbfile, nwbfile_path=nwbfile_path)

    print(f"\n✅ Behavioral data conversion complete!")
    print(f"📁 Files saved to: {nwb_files_dir}")
    print(f"📊 Processed {len(sessions_list)} sessions")

    # Show file summary
    nwb_files = list(nwb_files_dir.glob("figure7_behavioral_*.nwb"))
    total_size_mb = sum(f.stat().st_size for f in nwb_files) / (1024 * 1024)
    print(f"💾 Total size: {total_size_mb:.1f} MB")

    print(f"\n🎯 To reproduce Figure 7J:")
    print("1. Load NWB files and extract AIMScoreTable")
    print("2. Convert to DataFrame: table.to_dataframe()")
    print("3. Group by time_minutes and genotype for plotting")
    print(f"\n💡 Run with --verbose for detailed processing information")
    print(f"💡 Example: python {Path(__file__).name} --verbose")
