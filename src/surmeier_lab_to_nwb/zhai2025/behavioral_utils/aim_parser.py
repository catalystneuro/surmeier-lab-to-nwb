"""
AIM (Abnormal Involuntary Movement) data parser for Zhai et al. 2025.

This module provides functions to parse AIM scoring data from Excel files
and convert it to tidy CSV format with proper session date handling.
"""

import re
from datetime import datetime
from pathlib import Path
from typing import List, Optional, Tuple

import numpy as np
import pandas as pd


def parse_date_from_cell(cell_val) -> Optional[datetime]:
    """Parse various date formats from Excel cells."""
    if pd.isna(cell_val):
        return None

    if isinstance(cell_val, datetime):
        return cell_val

    cell_str = str(cell_val).strip()

    # Pattern 1: Date:MMDDYYYY format
    match = re.search(r"Date:\s*(\d{8})", cell_str)
    if match:
        date_str = match.group(1)
        month = date_str[:2]
        day = date_str[2:4]
        year = date_str[4:8]
        try:
            return datetime.strptime(f"{year}-{month}-{day}", "%Y-%m-%d")
        except:
            pass

    # Pattern 2: M/D/YYYY format
    match = re.search(r"(\d{1,2})/(\d{1,2})/(\d{4})", cell_str)
    if match:
        try:
            return datetime.strptime(match.group(0), "%m/%d/%Y")
        except:
            pass

    return None


def find_session_boundaries(df: pd.DataFrame) -> List[Tuple[int, datetime]]:
    """Find all session start rows and their dates in the dataframe."""
    sessions = []

    for i in range(len(df)):
        date = parse_date_from_cell(df.iloc[i, 0])
        if date:
            sessions.append((i, date))

    return sessions


def parse_aim_excel_to_wide_format(
    excel_path: Path, genotype_csv_path: Optional[Path] = None, output_path: Optional[Path] = None
) -> pd.DataFrame:
    """
    Parse AIM scoring data from Excel file with correct session date handling.

    Parameters
    ----------
    excel_path : Path
        Path to the AIM testing Excel file
    genotype_csv_path : Path, optional
        Path to the genotype CSV file for automatic genotype assignment
    output_path : Path, optional
        Path to save the tidy CSV file. If None, CSV is not saved.

    Returns
    -------
    pd.DataFrame
        Tidy DataFrame with columns:
        - sheet_name: Original sheet name as integer
        - session_date: Session date in YYYY-MM-DD format
        - animal_id: Animal identifier as integer
        - genotype: Animal genotype (CDGI Knockout, Wild Type, or unknown)
        - score_type: Type of AIM score (axial, limb, orolingual)
        - 20min, 40min, 60min, 80min, 100min, 120min: Score values at each time point

    Notes
    -----
    This function correctly handles:
    - Multiple sessions per sheet with different dates
    - Proper association of animals with their scoring data
    - Genotype assignment from external CSV file
    - Various date formats in Excel sheets
    """
    excel_file = pd.ExcelFile(excel_path)
    all_data = []

    # Load genotype data if provided
    genotype_map = {}
    if genotype_csv_path and genotype_csv_path.exists():
        genotype_df = pd.read_csv(genotype_csv_path, skiprows=4)  # Skip first 4 rows
        for _, row in genotype_df.iterrows():
            if pd.notna(row["MOUSE"]) and pd.notna(row["Category"]):
                animal_id = str(row["MOUSE"]).replace(",", "")  # Remove commas
                category = row["Category"]
                if category == "ko":
                    genotype_map[animal_id] = "CDGI Knockout"
                elif category in ["CONT", "cont"]:
                    genotype_map[animal_id] = "Wild Type"

    for sheet_name in excel_file.sheet_names:
        # Skip non-date sheets
        if not sheet_name.isdigit() or len(sheet_name) != 8:
            continue

        df = pd.read_excel(excel_path, sheet_name=sheet_name, header=None)

        # Find all session boundaries
        sessions = find_session_boundaries(df)

        if not sessions:
            # No clear date markers, use sheet name as fallback
            month = sheet_name[:2]
            day = sheet_name[2:4]
            year = sheet_name[4:8]
            default_date = datetime.strptime(f"{year}-{month}-{day}", "%Y-%m-%d")
            sessions = [(0, default_date)]

        # Process each session
        for session_idx, (start_row, session_date) in enumerate(sessions):
            # Determine end row for this session
            if session_idx < len(sessions) - 1:
                end_row = sessions[session_idx + 1][0]
            else:
                end_row = len(df)

            # For the first session, start from the beginning of the sheet
            # For subsequent sessions, start after the date marker
            if session_idx == 0:
                search_start = 0
            else:
                search_start = start_row + 1

            # Find all animals in this session
            i = search_start
            while i < end_row:
                row = df.iloc[i]

                # Look for animal ID (ET# pattern)
                if pd.notna(row.iloc[0]) and "ET#" in str(row.iloc[0]):
                    animal_id_str = str(row.iloc[0])

                    # Extract genotype from cell if present
                    genotype = "unknown"
                    if "(KO)" in animal_id_str:
                        genotype = "CDGI Knockout"
                    elif "(WT)" in animal_id_str:
                        genotype = "Wild Type"

                    # Clean animal ID
                    animal_id = animal_id_str.replace("(KO)", "").replace("(WT)", "").replace("ET#", "").strip()

                    # Override with genotype map if available
                    if animal_id in genotype_map:
                        genotype = genotype_map[animal_id]

                    # Find scoring rows - look for the pattern where animal ID and Ax are in same row
                    ax_row = None
                    li_row = None
                    ol_row = None

                    # Check if Ax is in the same row as the animal ID
                    if len(df.iloc[i]) > 1 and "Ax" in str(df.iloc[i].iloc[1]):
                        ax_row = i
                        # Li should be in next row, Ol in row after that
                        if i + 1 < end_row and len(df.iloc[i + 1]) > 1 and "Li" in str(df.iloc[i + 1].iloc[1]):
                            li_row = i + 1
                        if i + 2 < end_row and len(df.iloc[i + 2]) > 1 and "Ol" in str(df.iloc[i + 2].iloc[1]):
                            ol_row = i + 2

                    # Extract scores for each category
                    time_points = ["20min", "40min", "60min", "80min", "100min", "120min"]

                    # Process each score type
                    for score_row, score_type in [(ax_row, "axial"), (li_row, "limb"), (ol_row, "orolingual")]:
                        if score_row is not None:
                            scores = {}

                            for col_idx, time_pt in enumerate(time_points, start=2):
                                if col_idx < len(df.columns):
                                    val = df.iloc[score_row, col_idx]
                                    if pd.notna(val) and str(val).strip() != "":
                                        try:
                                            scores[time_pt] = float(val)
                                        except (ValueError, TypeError):
                                            scores[time_pt] = np.nan
                                    else:
                                        scores[time_pt] = np.nan
                                else:
                                    scores[time_pt] = np.nan

                            # Add row to data
                            row_data = {
                                "sheet_name": int(sheet_name),
                                "session_date": session_date.strftime("%Y-%m-%d"),
                                "animal_id": int(animal_id),
                                "genotype": genotype,
                                "score_type": score_type,
                            }
                            row_data.update(scores)
                            all_data.append(row_data)

                    i = max(ax_row or i, li_row or i, ol_row or i) + 1
                else:
                    i += 1

    # Create DataFrame
    wide_df = pd.DataFrame(all_data)

    # Sort by sheet, date, animal, and score type
    if not wide_df.empty:
        wide_df = wide_df.sort_values(["sheet_name", "session_date", "animal_id", "score_type"])
        wide_df = wide_df.reset_index(drop=True)

    # Save to CSV if output path provided
    if output_path:
        wide_df.to_csv(output_path, index=False)
        print(f"AIM data saved to: {output_path}")

    return wide_df
