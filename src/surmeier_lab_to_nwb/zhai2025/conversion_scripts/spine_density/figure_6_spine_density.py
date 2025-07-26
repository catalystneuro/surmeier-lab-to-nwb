"""
Figure 6 Spine Density Conversion Script - Zhai et al. 2025
========================================================

This script converts two-photon spine density imaging data from Figure 6 of
Zhai et al. 2025 into NWB (Neurodata Without Borders) format. The data examines
dendritic spine density changes in a mouse model of Parkinson's disease
and levodopa-induced dyskinesia.

For detailed experimental context, protocols, data structure, and analysis methods,
see: /src/surmeier_lab_to_nwb/zhai2025/conversion_notes_folder/figure_6_conversion_notes.md
"""

from pathlib import Path

from neuroconv.tools import configure_and_write_nwbfile
from pynwb import NWBFile

from surmeier_lab_to_nwb.zhai2025.conversion_scripts.conversion_utils import str_to_bool
from surmeier_lab_to_nwb.zhai2025.conversion_scripts.spine_density.spine_density_utils import (
    convert_spine_density_session_to_nwbfile,
)


def convert_data_to_nwb(session_folder_path: Path, condition: str, verbose: bool = False) -> NWBFile:
    """
    Convert Figure 6 spine density data to NWB format.

    This is a wrapper function that calls the shared conversion function with
    Figure 6-specific configuration and session ID parameters.

    Parameters
    ----------
    session_folder_path : Path
        Path to the session folder containing image stack subfolders
    condition : str
        Experimental condition (e.g., "control", "M1R antagonist")
    verbose : bool, default=False
        Enable verbose output showing detailed processing information

    Returns
    -------
    NWBFile
        NWB file with the converted data

    Notes
    -----
    This function handles Figure 6 iSPN spine density data with M1R antagonist
    conditions in the OFF state.
    """
    # Configuration for Figure 6 spine density experiments
    figure_6_config = {
        "metadata_key": "figure_6_spine_density",
    }

    # Local mapping for Figure 6 conditions to revised schema tokens
    figure_6_mappings = {
        "control": {"state": "OffState", "pharm": "none"},
        "M1R antagonist": {"state": "OffState", "pharm": "M1RaThp"},  # Muscarinic type 1 receptor antagonist
    }

    if condition not in figure_6_mappings:
        raise ValueError(f"Unknown condition: {condition}")

    state = figure_6_mappings[condition]["state"]
    pharmacology = figure_6_mappings[condition]["pharm"]

    # Build session ID parameters using revised schema
    session_id_parameters = {
        "fig": "F6",
        "meas_comp": "SpineDens",  # Spine density measurement
        "cell_type": "iSPN",  # Indirect pathway SPN
        "state": state,
        "pharm": pharmacology,
        "geno": "WT",  # Wild-type for all Figure 6
    }

    return convert_spine_density_session_to_nwbfile(
        session_folder_path=session_folder_path,
        condition=condition,
        figure_config=figure_6_config,
        session_id_parameters=session_id_parameters,
        verbose=verbose,
    )


if __name__ == "__main__":
    import argparse
    import logging

    from tqdm import tqdm

    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Convert Figure 6 spine density data to NWB format")
    parser.add_argument(
        "--stub-test",
        type=str_to_bool,
        default=True,
        help="Process only first 2 files per condition for testing (default: True). Use --stub-test=False for full processing.",
    )

    args = parser.parse_args()
    stub_test = args.stub_test

    verbose = False

    logging.getLogger("tifffile").setLevel(logging.ERROR)

    # Define the base path to the data - Figure 6 spine density (two-photon method, iSPNs)
    base_path = Path("/home/heberto/development/surmeier-lab-to-nwb/link_to_raw_data/Figure 6/Spine density/")

    # Create nwb_files directory at root level
    root_dir = Path(__file__).parent.parent.parent.parent.parent.parent  # Go up to repo root
    nwb_files_dir = root_dir / "nwb_files" / "spine_density" / "figure_6"
    nwb_files_dir.mkdir(parents=True, exist_ok=True)

    # Figure 6 conditions use iSPNs (indirect pathway)
    conditions = ["control", "M1R antagonist"]

    for condition in conditions:
        condition_path = base_path / condition
        assert condition_path.exists(), f"Base path does not exist: {condition_path}"

        # Get all session folders
        session_folders = [f for f in condition_path.iterdir() if f.is_dir()]
        session_folders.sort()

        # Apply stub_test filtering if enabled
        if stub_test:
            session_folders = session_folders[:2]

        # Use tqdm for progress bar when verbose is disabled
        session_iterator = tqdm(session_folders, desc=f"Converting Figure6 SpineDensity {condition}", unit=" session")

        for session_folder_path in session_iterator:

            # Convert data to NWB format
            nwbfile = convert_data_to_nwb(
                session_folder_path=session_folder_path,
                condition=condition,
                verbose=verbose,
            )

            # Create output filename using session_id from nwbfile
            nwbfile_path = nwb_files_dir / f"{nwbfile.session_id}.nwb"

            # Write NWB file
            configure_and_write_nwbfile(nwbfile, nwbfile_path=nwbfile_path)
