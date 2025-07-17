"""
Behavioral Data Interface for Zhai et al. 2025.

This module provides interfaces for parsing and converting behavioral data from the
Zhai et al. 2025 study, including AIM (Abnormal Involuntary Movement) scoring data.
"""

import re
from datetime import datetime
from pathlib import Path
from typing import Dict, Optional

import pandas as pd


def load_genotype_data(data_connections_path: Path) -> Dict[str, Dict[str, str]]:
    """
    Load genotype data from the Data Connections CSV file.

    Parameters
    ----------
    data_connections_path : Path
        Path to the data_connections_D2_figures_3_4_6_7_8.csv file

    Returns
    -------
    Dict[str, Dict[str, str]]
        Dictionary mapping animal_id to genotype information
    """
    df = pd.read_csv(data_connections_path, header=4)
    genotype_map = {}

    for _, row in df.iterrows():
        if pd.notna(row["MOUSE"]) and str(row["MOUSE"]).strip():
            mouse_id = str(row["MOUSE"]).replace(",", "").strip()
            animal_type = row["ANIMAL"] if pd.notna(row["ANIMAL"]) else "unknown"
            category = row["Category"] if pd.notna(row["Category"]) else "unknown"

            # Skip non-numeric mouse IDs
            if mouse_id.isdigit():
                genotype_map[mouse_id] = {"animal_type": animal_type, "category": category}

    return genotype_map


def parse_session_date(cell_value: str) -> Optional[str]:
    """
    Parse session date from various formats found in the first column.

    Parameters
    ----------
    cell_value : str
        Cell value that might contain a session date

    Returns
    -------
    Optional[str]
        Parsed date in YYYY-MM-DD format, or None if not a valid session date
    """
    if not isinstance(cell_value, str):
        cell_value = str(cell_value)

    cell_value = cell_value.strip()

    # Pattern 1: MM/DD/YYYY (Nth session) format
    pattern1 = r"(\d{1,2})/(\d{1,2})/(\d{4})\s*\([^)]*session\)"
    match1 = re.search(pattern1, cell_value)
    if match1:
        month, day, year = match1.groups()
        return f"{year}-{month.zfill(2)}-{day.zfill(2)}"

    # Pattern 2: Date:MMDDYYYY format
    pattern2 = r"Date:(\d{2})(\d{2})(\d{4})"
    match2 = re.search(pattern2, cell_value)
    if match2:
        month, day, year = match2.groups()
        return f"{year}-{month}-{day}"

    # Pattern 3: Date:MM/DD(Nth session) format - infer year from sheet context
    pattern3 = r"Date:(\d{1,2})/(\d{1,2})\([^)]*session\)"
    match3 = re.search(pattern3, cell_value)
    if match3:
        month, day = match3.groups()
        # This pattern needs year context, will be handled in main function
        return f"INCOMPLETE-{month.zfill(2)}-{day.zfill(2)}"

    # Pattern 4: YYYY-MM-DD HH:MM:SS datetime format
    try:
        if len(cell_value) > 10 and ":" in cell_value:
            dt = datetime.fromisoformat(cell_value.replace(" 00:00:00", ""))
            return dt.strftime("%Y-%m-%d")
    except:
        pass

    return None


def is_session_marker(cell_value) -> bool:
    """Check if a cell value looks like a session date marker."""
    if pd.isna(cell_value):
        return False

    cell_str = str(cell_value).strip().lower()

    # Must contain session-related keywords and date-like patterns
    has_session_keyword = any(keyword in cell_str for keyword in ["session", "date:"])
    has_date_pattern = any(pattern in cell_str for pattern in ["/", "2017", "2018", "2019", "2020"])

    return has_session_keyword or has_date_pattern


def extract_genotype_from_text(text: str) -> Optional[str]:
    """Extract genotype from text string."""
    if not isinstance(text, str):
        text = str(text)

    text = text.strip().upper()

    if "(KO)" in text:
        return "KO"
    elif "(WT)" in text:
        return "WT"
    elif text == "(KO)":
        return "KO"
    elif text == "(WT)":
        return "WT"

    return None


def parse_aim_excel_to_tidy_csv(excel_path: Path, output_path: Optional[Path] = None) -> pd.DataFrame:
    """
    Parse AIM scoring data from Excel file and convert to tidy CSV format.

    This function extracts multiple sessions per sheet, with accurate session dates parsed
    from the first column, and improved genotype extraction from various locations.

    Parameters
    ----------
    excel_path : Path
        Path to the AIM testing Excel file
    output_path : Path, optional
        Path to save the tidy CSV file. If None, CSV is not saved.

    Returns
    -------
    pd.DataFrame
        Wide-format DataFrame with columns:
        - sheet_name: Original sheet name (e.g., "02232018")
        - session_date: Actual session date parsed from first column (YYYY-MM-DD format)
        - animal_id: Subject identifier (ET number without prefix)
        - genotype: Animal genotype (KO, WT, or unknown)
        - score_type: Type of AIM score (axial, limb, orolingual)
        - 20min, 40min, 60min, 80min, 100min, 120min: Score values at each time point

    Notes
    -----
    Each sheet contains multiple behavioral testing sessions. Session dates are parsed
    from various formats in the first column. Genotype information is extracted from
    animal IDs or subsequent rows and tracked across sessions for each animal.

    """
    excel_file = pd.ExcelFile(excel_path)
    tidy_data = []

    # Time points for AIM scoring (in minutes post-L-DOPA)
    time_points = [20, 40, 60, 80, 100, 120]

    # Genotype registry to track genotype per animal across sessions
    genotype_registry = {}

    for sheet_name in excel_file.sheet_names:
        print(f"Processing sheet: {sheet_name}")

        # Validate sheet name format (MMDDYYYY)
        if len(sheet_name) != 8 or not sheet_name.isdigit():
            raise ValueError(f'Invalid sheet name format: "{sheet_name}". Expected MMDDYYYY format.')

        df = pd.read_excel(excel_path, sheet_name=sheet_name, header=None)

        # Track current session date and infer year from sheet name for incomplete dates
        sheet_year = sheet_name[4:8]
        current_session_date = None

        # Parse through the sheet row by row
        i = 0
        while i < len(df):
            cell_value = df.iloc[i, 0]

            # Check if this row contains a session date marker
            if pd.notna(cell_value) and is_session_marker(cell_value):
                parsed_date = parse_session_date(str(cell_value))
                if parsed_date:
                    # Handle incomplete dates that need year from sheet context
                    if parsed_date.startswith("INCOMPLETE-"):
                        month_day = parsed_date.replace("INCOMPLETE-", "")
                        current_session_date = f"{sheet_year}-{month_day}"
                    else:
                        current_session_date = parsed_date
                    print(f"  Found session: {current_session_date}")
                i += 1
                continue

            # Look for animal ID (ET# pattern)
            if pd.notna(cell_value) and "ET#" in str(cell_value):
                animal_id = str(cell_value)

                # Clean animal ID
                clean_animal_id = animal_id.replace("(KO)", "").replace("(WT)", "").replace("ET#", "").strip()

                # Extract genotype from animal ID row
                genotype = extract_genotype_from_text(animal_id)

                # If no genotype found in animal ID, check the next row
                if not genotype and i + 1 < len(df):
                    next_row_value = df.iloc[i + 1, 0]
                    if pd.notna(next_row_value):
                        genotype = extract_genotype_from_text(str(next_row_value))

                # Update genotype registry
                if genotype:
                    genotype_registry[clean_animal_id] = genotype
                elif clean_animal_id in genotype_registry:
                    # Use previously found genotype for this animal
                    genotype = genotype_registry[clean_animal_id]
                else:
                    genotype = "unknown"

                # Skip if no current session date (shouldn't happen, but safety check)
                if not current_session_date:
                    print(f"  Warning: No session date found for animal {clean_animal_id}, skipping")
                    i += 1
                    continue

                # Find the scoring rows (Ax, Li, Ol)
                ax_row = None
                li_row = None
                ol_row = None

                # Look in current row and next few rows for scoring categories
                # But stop searching when we hit another animal ID
                for j in range(i, min(i + 5, len(df))):
                    if j < len(df) and len(df.iloc[j]) > 1:
                        # Check if we've hit another animal ID (stop searching)
                        if j > i and pd.notna(df.iloc[j, 0]) and "ET#" in str(df.iloc[j, 0]):
                            break

                        cell_val = str(df.iloc[j].iloc[1]) if pd.notna(df.iloc[j].iloc[1]) else ""
                        if "Ax" in cell_val and ax_row is None:
                            ax_row = j
                        elif "Li" in cell_val and li_row is None:
                            li_row = j
                        elif "Ol" in cell_val and ol_row is None:
                            ol_row = j

                # Extract scores for each time point and create wide format
                if ax_row is not None:
                    # Create base record for each score type
                    for score_type, row_idx in [("axial", ax_row), ("limb", li_row), ("orolingual", ol_row)]:
                        if row_idx is not None and row_idx < len(df):
                            # Initialize record with basic info
                            record = {
                                "sheet_name": sheet_name,
                                "session_date": current_session_date,
                                "animal_id": clean_animal_id,
                                "genotype": genotype,
                                "score_type": score_type,
                            }

                            # Add scores for each time point as separate columns
                            for col_idx, time_point in enumerate(time_points):
                                col_data_idx = col_idx + 2  # Data starts at column 2

                                if col_data_idx < len(df.columns):
                                    score_val = df.iloc[row_idx, col_data_idx]
                                    if pd.notna(score_val) and str(score_val).strip() != "":
                                        try:
                                            record[f"{time_point}min"] = float(score_val)
                                        except ValueError:
                                            record[f"{time_point}min"] = None
                                    else:
                                        record[f"{time_point}min"] = None

                            # Only add record if it has at least one non-null score
                            if any(record.get(f"{tp}min") is not None for tp in time_points):
                                tidy_data.append(record)

                i += 4  # Skip to next animal
            else:
                i += 1

    # Convert to DataFrame
    wide_df = pd.DataFrame(tidy_data)

    # Sort by sheet_name, session_date, animal_id, score_type
    if not wide_df.empty:
        wide_df = wide_df.sort_values(["sheet_name", "session_date", "animal_id", "score_type"])
        wide_df = wide_df.reset_index(drop=True)

    # Save to CSV if output path provided
    if output_path:
        wide_df.to_csv(output_path, index=False)
        print(f"Wide CSV saved to: {output_path}")

    return wide_df
