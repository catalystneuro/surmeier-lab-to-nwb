#!/usr/bin/env python3
"""
Script to extract metadata from DANDI dataset 001538 and export to CSV.
This script iterates over all assets in the dataset and extracts experimental metadata
encoded in the session IDs, then exports everything to a CSV file.
"""

import csv
import os

from dandi.dandiapi import DandiAPIClient
from dotenv import load_dotenv

# Load environment variables from .env file
load_dotenv()


# Session ID Parsing Functions
def get_session_id(asset_path: str) -> str:
    """
    Extract session ID from DANDI asset path.

    DANDI encodes paths as:
    sub-<subject_id>/sub-<subject_id>_ses-<session_id>_[desc-<description>]_<modalities>.nwb

    Example path:
    'sub-SubjectRecordedAt20160523154318/sub-SubjectRecordedAt20160523154318_ses-F3++SomExc++iSPN++OffState++none++WT++20160523154318_icephys.nwb'

    The session_id contains experimental metadata separated by '++':
    F3++SomExc++iSPN++OffState++none++WT++20160523154318
    """
    if not asset_path:
        return ""
    try:
        bottom_level_path = asset_path.split("/")[1]  # Remove top level subject
        session_id_with_ses_prefix = bottom_level_path.split("_")[1]
        session_id = session_id_with_ses_prefix.split("-")[1]
        return session_id
    except (IndexError, AttributeError):
        return ""


def get_figure_number(session_id: str) -> str:
    """Extract which figure this data corresponds to (F1, F2, F3, etc.)"""
    if not session_id:
        return ""
    try:
        return session_id.split("++")[0]
    except IndexError:
        return ""


def get_measurement(session_id: str) -> str:
    """
    Extract measurement type:
    - SomExc: Somatic excitability (patch clamp at cell body)
    - DendExc: Dendritic excitability (patch clamp + 2-photon imaging)
    - DendSpine: Dendritic spine density measurements
    - BehavAIMs: Abnormal involuntary movement behavioral scoring
    - StriAChFP: Striatal acetylcholine fluorescent protein imaging
    """
    if not session_id:
        return ""
    try:
        return session_id.split("++")[1]
    except IndexError:
        return ""


def get_cell_type(session_id: str) -> str:
    """
    Extract cell type:
    - dSPN: Direct pathway spiny projection neurons
    - iSPN: Indirect pathway spiny projection neurons
    - pan: Pan-neuronal (both types)
    """
    if not session_id:
        return ""
    try:
        return session_id.split("++")[2]
    except IndexError:
        return ""


def get_state(session_id: str) -> str:
    """
    Extract experimental state:
    - CTRL: Control condition
    - PD: Parkinson's disease model
    - OffState: No levodopa treatment
    - OnState: With levodopa treatment
    """
    if not session_id:
        return ""
    try:
        return session_id.split("++")[3]
    except IndexError:
        return ""


def get_pharmacology(session_id: str) -> str:
    """
    Extract pharmacological treatment:
    - none: No drugs applied
    - CNO: Clozapine N-oxide (DREADD activator)
    - D1RaSch: D1 receptor agonist SCH23390
    - M1RaOxoM: M1 receptor agonist Oxotremorine-M
    """
    if not session_id:
        return ""
    try:
        return session_id.split("++")[4]
    except IndexError:
        return ""


def get_genotype(session_id: str) -> str:
    """
    Extract animal genotype:
    - WT: Wild type
    - M1RCRISPR: M1 receptor knockout
    - iSPNM1RKO: iSPN-specific M1 receptor knockout
    """
    if not session_id:
        return ""
    try:
        return session_id.split("++")[5]
    except IndexError:
        return ""


def get_timestamp(session_id: str) -> str:
    """Extract recording timestamp (YYYYMMDDHHMMSS format)"""
    if not session_id:
        return ""
    try:
        return session_id.split("++")[6]
    except IndexError:
        return ""


def main():
    """Main function to extract metadata and create CSV"""

    # Load token from environment variable
    token = os.getenv("DANDI_API_TOKEN")
    if not token:
        raise ValueError("DANDI_API_TOKEN environment variable not set. Please set it with your DANDI API token.")

    print("Connecting to DANDI...")
    dandiset_id = "001538"
    client = DandiAPIClient(token=token)
    client.authenticate(token=token)

    dandiset = client.get_dandiset(dandiset_id, "draft")
    assets = dandiset.get_assets()
    assets_list = list(assets)

    print(f"Found {len(assets_list)} assets in dataset {dandiset_id}")

    # Prepare CSV output
    csv_filename = "dandi_001538_metadata.csv"

    # CSV headers
    headers = ["figure_number", "measurement", "cell_type", "state", "pharmacology", "genotype", "asset_path"]

    print(f"Extracting metadata and writing to {csv_filename}...")

    with open(csv_filename, "w", newline="", encoding="utf-8") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(headers)

        for i, asset in enumerate(assets_list):
            print(f"Processing asset {i+1}/{len(assets_list)}: {asset.path}")

            try:
                # Get S3 URL (commented out for performance)
                # s3_url = asset.get_content_url(follow_redirects=1, strip_query=False)

                # Extract metadata
                session_id = get_session_id(asset.path)
                figure_number = get_figure_number(session_id)
                measurement = get_measurement(session_id)
                cell_type = get_cell_type(session_id)
                state = get_state(session_id)
                pharmacology = get_pharmacology(session_id)
                genotype = get_genotype(session_id)
                # timestamp = get_timestamp(session_id)

                # Write row to CSV
                writer.writerow([figure_number, measurement, cell_type, state, pharmacology, genotype, asset.path])

            except Exception as e:
                print(f"Error processing asset {asset.path}: {e}")
                # Write row with error info
                writer.writerow(["", "", "", "", "", "", asset.path])

    print(f"Metadata extraction complete! Results saved to {csv_filename}")

    # Print summary statistics
    print("\nSummary of extracted metadata:")
    with open(csv_filename, "r", encoding="utf-8") as csvfile:
        reader = csv.DictReader(csvfile)
        rows = list(reader)

        print(f"Total assets processed: {len(rows)}")

        # Count unique values for each metadata field
        for field in ["figure_number", "measurement", "cell_type", "state", "pharmacology", "genotype"]:
            unique_values = set(row[field] for row in rows if row[field])
            print(f"Unique {field}: {sorted(unique_values)}")


if __name__ == "__main__":
    main()
