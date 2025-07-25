# -*- coding: utf-8 -*-
"""
Figure 2 Optical Stimuli Conversion Script - Zhai et al. 2025
============================================================

This script converts Sr²⁺-oEPSC (strontium-substituted optogenetically evoked
postsynaptic current) experimental data from Figure 2 of Zhai et al. 2025 into
NWB (Neurodata Without Borders) format. The data examines corticostriatal
synaptic strength changes in LID on- vs off-states using optogenetic stimulation
with voltage clamp recordings.

For detailed experimental context, protocols, data structure, and analysis methods,
see: /src/surmeier_lab_to_nwb/zhai2025/conversion_notes_folder/figure_2_conversion_notes.md
"""
from pathlib import Path

from neuroconv.tools import configure_and_write_nwbfile
from pynwb import NWBFile

from surmeier_lab_to_nwb.zhai2025.conversion_scripts.conversion_utils import str_to_bool
from surmeier_lab_to_nwb.zhai2025.conversion_scripts.optical_stimulation.optical_stimulation_utils import (
    convert_optical_stimulation_session_to_nwbfile,
)


def convert_session_to_nwbfile(session_folder_path: Path, condition: str, verbose: bool = False) -> NWBFile:
    """
    Convert a single session of Figure 2 Sr²⁺-oEPSC data to NWB format.

    This is a wrapper function that calls the shared conversion function with
    Figure 2-specific configuration and session ID parameters.

    Parameters
    ----------
    session_folder_path : Path
        Path to the session folder containing sweep recordings
    condition : str
        Experimental condition (e.g., "LID off-state", "LID on-state")
    verbose : bool, default=False
        Enable verbose output showing detailed processing information

    Returns
    -------
    NWBFile
        NWB file with the converted data

    Notes
    -----
    This function handles Figure 2 dSPN optical stimulation data (Sr²⁺-oEPSC protocol)
    in LID on- vs off-states.
    """
    # Configuration for Figure 2 optical stimulation experiments
    figure_2_config = {
        "metadata_key": "figure_2_optical_stimuli",
        "cell_type": "dSPN",
        "subject_id_prefix": "dSPN_mouse",
    }

    # Local mapping for Figure 2 conditions to revised schema tokens
    figure_2_mappings = {
        "LID off-state": {"state": "OffState", "pharm": "none"},
        "LID on-state": {"state": "OnState", "pharm": "none"},
    }

    if condition not in figure_2_mappings:
        raise ValueError(f"Unknown condition: {condition}")

    state = figure_2_mappings[condition]["state"]
    pharmacology = figure_2_mappings[condition]["pharm"]

    # Build session ID parameters using revised schema
    session_id_parameters = {
        "fig": "F2",
        "meas_comp": "oEPSC",  # Sr²⁺ oEPSC measurement
        "cell_type": "dSPN",  # Direct pathway SPN
        "state": state,
        "pharm": pharmacology,
        "geno": "WT",  # Wild-type for all Figure 2
    }

    return convert_optical_stimulation_session_to_nwbfile(
        session_folder_path=session_folder_path,
        condition=condition,
        figure_config=figure_2_config,
        session_id_parameters=session_id_parameters,
        verbose=verbose,
    )


if __name__ == "__main__":
    import argparse
    import logging
    import warnings

    from tqdm import tqdm

    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Convert Figure 2 optical stimuli data to NWB format")
    parser.add_argument(
        "--stub-test",
        type=str_to_bool,
        default=True,
        help="Process only first 2 files per condition for testing (default: True). Use --stub-test=False for full processing.",
    )

    args = parser.parse_args()
    stub_test = args.stub_test

    # Control verbose output from here
    verbose = False  # Set to True for detailed output

    # Suppress tifffile warnings
    logging.getLogger("tifffile").setLevel(logging.ERROR)

    # Suppress specific warnings
    warnings.filterwarnings("ignore", message="invalid value encountered in divide")
    warnings.filterwarnings("ignore", message=".*no datetime before year 1.*")

    # Define the base path to the data
    base_path = Path("./link_to_raw_data/Figure 2_SF1A/Sr-oEPSC")
    if not base_path.exists():
        raise FileNotFoundError(f"Base path does not exist: {base_path}")

    # Create nwb_files directory at root level
    root_dir = Path(__file__).parent.parent.parent.parent.parent.parent  # Go up to repo root
    nwb_files_dir = root_dir / "nwb_files" / "optical_stimulation" / "figure_2"
    nwb_files_dir.mkdir(parents=True, exist_ok=True)

    conditions = ["LID on-state", "LID off-state"]

    for condition in conditions:
        condition_path = base_path / condition
        if not condition_path.exists():
            raise FileNotFoundError(f"Expected condition path does not exist: {condition_path}")

        # Get all session folders
        session_folders = [f for f in condition_path.iterdir() if f.is_dir()]
        session_folders.sort()

        # Apply stub_test filtering if enabled
        if stub_test:
            session_folders = session_folders[:2]

        # Use tqdm for progress bar when verbose is disabled
        session_iterator = tqdm(session_folders, desc=f"Converting Figure2 OpticalStimuli {condition}", unit=" session")

        for session_folder in session_iterator:

            # Convert data to NWB format
            nwbfile = convert_session_to_nwbfile(
                session_folder_path=session_folder,
                condition=condition,
                verbose=verbose,
            )

            # Create output filename using session_id from nwbfile
            nwbfile_path = nwb_files_dir / f"{nwbfile.session_id}.nwb"

            configure_and_write_nwbfile(nwbfile, nwbfile_path=nwbfile_path)
