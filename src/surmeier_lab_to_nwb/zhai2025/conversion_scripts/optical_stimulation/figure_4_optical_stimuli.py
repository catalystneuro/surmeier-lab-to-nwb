#!/usr/bin/env python3
"""
Figure 4 Optical Stimuli Conversion Script

This script converts Figure 4 Sr²⁺-oEPSC data from the Zhai et al. 2025 paper.
This data tests iSPN (indirect spiny projection neurons) responses to optogenetic
stimulation of corticostriatal terminals in the LID (L-DOPA-induced dyskinesia) model.

The experiment uses the same Sr²⁺-oEPSC protocol as Figure 2 but on iSPNs instead of dSPNs.
Optogenetic stimulation is delivered via whole-field LED (470nm) targeting ChR2-expressing
corticostriatal terminals, with 0.3ms pulses at 20ms into each 30-second inter-sweep interval.

Data structure:
- Figure 4_SF1B_SF5/Sr-oEPSC/LID on-state/: Sessions with LID-inducing L-DOPA treatment
- Figure 4_SF1B_SF5/Sr-oEPSC/LID off-state/: Sessions without L-DOPA treatment
- Each session contains 10 sweep recordings from iSPN cells

Key differences from Figure 2:
- Cell type: iSPN (indirect spiny projection neurons) instead of dSPN
- Different experimental sessions and dates
- Same optogenetic protocol and LED stimulation parameters
- Same virus (AAV5-hSyn-hChR2(H134R)-EYFP) and injection coordinates

"""

import re
from pathlib import Path

from neuroconv.tools.nwb_helpers import configure_and_write_nwbfile
from pynwb import NWBFile

from surmeier_lab_to_nwb.zhai2025.conversion_scripts.conversion_utils import str_to_bool
from surmeier_lab_to_nwb.zhai2025.conversion_scripts.optical_stimulation.optical_stimulation_utils import (
    convert_optical_stimulation_session_to_nwbfile,
)


def convert_session_to_nwbfile(
    session_folder_path: Path,
    condition: str,
    verbose: bool = False,
) -> NWBFile:
    """
    Convert a single session of Figure 4 Sr²⁺-oEPSC data to NWB format.

    This is a wrapper function that calls the shared conversion function with
    Figure 4-specific configuration and session ID parameters.

    Parameters
    ----------
    session_folder_path : Path
        Path to the session folder containing sweep recordings
    condition : str
        Experimental condition ("LID on-state" or "LID off-state")
    verbose : bool, default=False
        Enable verbose output showing detailed processing information

    Returns
    -------
    NWBFile
        NWB file with the converted data

    Notes
    -----
    This function handles Figure 4 iSPN optical stimulation data (Sr²⁺-oEPSC protocol)
    in LID on- vs off-states.
    """
    # Configuration for Figure 4 optical stimulation experiments
    figure_4_config = {
        "metadata_key": "figure_4_optical_stimuli",
        "cell_type": "iSPN",
        "subject_id_prefix": "iSPN_optical_mouse",
    }

    # Local mapping for Figure 4 conditions to revised schema tokens
    figure_4_mappings = {
        "LID off-state": {"state": "OffState", "pharm": "none"},
        "LID on-state": {"state": "OnState", "pharm": "none"},
    }

    if condition not in figure_4_mappings:
        raise ValueError(f"Unknown condition: {condition}")

    state = figure_4_mappings[condition]["state"]
    pharmacology = figure_4_mappings[condition]["pharm"]

    # Build session ID parameters using revised schema
    session_id_parameters = {
        "fig": "F4",
        "meas_comp": "oEPSC",  # Sr²⁺ oEPSC measurement
        "cell_type": "iSPN",  # Indirect pathway SPN
        "state": state,
        "pharm": pharmacology,
        "geno": "WT",  # Wild-type for all Figure 4
    }

    return convert_optical_stimulation_session_to_nwbfile(
        session_folder_path=session_folder_path,
        condition=condition,
        figure_config=figure_4_config,
        session_id_parameters=session_id_parameters,
        verbose=verbose,
    )


def main():
    """Main conversion function for Figure 4 Sr²⁺-oEPSC data."""
    import argparse

    from tqdm import tqdm

    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Convert Figure 4 optical stimuli data to NWB format")
    parser.add_argument(
        "--stub-test",
        type=str_to_bool,
        default=True,
        help="Process only first 2 files per condition for testing (default: True). Use --stub-test=False for full processing.",
    )

    args = parser.parse_args()
    stub_test = args.stub_test

    # Control verbose output
    verbose = False  # Set to True for detailed output

    # Define raw data and output paths
    raw_data_root = Path("/home/heberto/development/surmeier-lab-to-nwb/link_to_raw_data/Figure 4_SF1B_SF5/Sr-oEPSC")
    output_root = Path("/home/heberto/development/surmeier-lab-to-nwb/nwb_files/optical_stimulation/figure_4")

    # Create output directory
    output_root.mkdir(parents=True, exist_ok=True)

    # Define conditions to process
    conditions = ["LID on-state", "LID off-state"]

    # Process each condition
    for condition in conditions:

        condition_path = raw_data_root / condition
        if not condition_path.exists():
            raise FileNotFoundError(
                f"Condition path does not exist: {condition_path}. Please check the raw data directory."
            )

        # Find all session folders
        session_folders = [
            folder for folder in condition_path.iterdir() if folder.is_dir() and re.match(r"\d{8}[a-z]", folder.name)
        ]

        session_folders.sort()  # Sort by session name

        # Apply stub_test filtering if enabled
        if stub_test:
            session_folders = session_folders[:2]

        # Use tqdm for progress bar when verbose is disabled
        session_iterator = tqdm(session_folders, desc=f"Converting Figure4 OpticalStimuli {condition}", unit=" session")

        # Process each session
        for session_folder in session_iterator:

            # Convert session to NWB
            nwbfile = convert_session_to_nwbfile(
                session_folder_path=session_folder,
                condition=condition,
                verbose=verbose,
            )

            # Create output filename using session_id from nwbfile
            nwbfile_path = output_root / f"{nwbfile.session_id}.nwb"

            # Write NWB file
            configure_and_write_nwbfile(nwbfile, nwbfile_path=nwbfile_path)


if __name__ == "__main__":
    main()
