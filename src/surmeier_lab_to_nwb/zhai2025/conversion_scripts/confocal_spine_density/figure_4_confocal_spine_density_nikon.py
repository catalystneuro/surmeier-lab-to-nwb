import re
import uuid
from pathlib import Path
from typing import Any, Dict

from neuroconv.tools import configure_and_write_nwbfile
from neuroconv.utils import dict_deep_update, load_dict_from_file
from pynwb import NWBFile
from pynwb.file import Subject

from surmeier_lab_to_nwb.zhai2025.interfaces.image_stack_interfaces import (
    NikonImageStackInterface,
)


def parse_nd2_filename(nd2_file: Path) -> Dict[str, Any]:
    """
    Parse ND2 filename to extract metadata.

    Example filenames:
    - "3824-slide 2-slice 2-cell 2 proxi.nd2"
    - "8939-slide 1-slice 3-cell 2 den2 proxi.nd2"
    - "5107-slide 2-slice 2-cell 1 den 2 proxi.nd2"

    Parameters
    ----------
    nd2_file : Path
        Path to ND2 file

    Returns
    -------
    Dict[str, Any]
        Dictionary containing extracted metadata
    """
    filename = nd2_file.stem

    # Extract animal/session ID (first part before dash)
    parts = filename.split("-")
    animal_id = parts[0] if parts else "unknown"

    # Extract slide, slice, cell information
    slide_match = re.search(r"slide\s+(\d+)", filename, re.IGNORECASE)
    slice_match = re.search(r"slice\s+(\d+)", filename, re.IGNORECASE)
    cell_match = re.search(r"cell\s+(\d+)", filename, re.IGNORECASE)

    slide_num = slide_match.group(1) if slide_match else "1"
    slice_num = slice_match.group(1) if slice_match else "1"
    cell_num = cell_match.group(1) if cell_match else "1"

    # Extract dendrite information
    dendrite_info = ""
    if "den2" in filename.lower() or "den 2" in filename.lower():
        dendrite_info = "den2"
    elif "den3" in filename.lower() or "den 3" in filename.lower():
        dendrite_info = "den3"
    elif "den" in filename.lower():
        dendrite_info = "den1"

    # Extract location (proximal/distal)
    location = "proximal" if "proxi" in filename.lower() else "unknown"

    # Create container name - simplified format
    container_name = f"ImagesCell{cell_num}"
    if dendrite_info:
        container_name += f"_{dendrite_info}"

    return {
        "animal_id": animal_id,
        "slide_number": slide_num,
        "slice_number": slice_num,
        "cell_number": cell_num,
        "dendrite_info": dendrite_info,
        "location": location,
        "container_name": container_name,
        "description": (
            f"High-resolution confocal image stack of {location} dendrite from iSPN cell {cell_num} "
            f"({dendrite_info if dendrite_info else 'primary dendrite'}). "
            f"Acquired from slide {slide_num}, slice {slice_num} using Nikon AXR confocal microscope "
            f"with 60x oil immersion objective (NA=1.49), 0.09 μm pixels, 0.125 μm z-steps. "
            f"Images processed with Nikon NIS-Elements and analyzed using Imaris 10.0.0 for "
            f"spine density quantification with supervised learning algorithms."
        ),
    }


def convert_session_to_nwbfile(nd2_file: Path, condition: str, verbose: bool = False) -> NWBFile:
    """
    Convert single confocal ND2 file to NWB format.

    Parameters
    ----------
    nd2_file : Path
        Path to ND2 file
    condition : str
        Experimental condition (e.g., 'control', 'LID on-state')
    verbose : bool, default=False
        Enable verbose output

    Returns
    -------
    NWBFile
        The populated NWB file object
    """
    # Parse filename metadata
    file_info = parse_nd2_filename(nd2_file)
    if verbose:
        print(f"  File: {nd2_file.name}")
        print(f"  Session ID: {file_info['session_id']}")
        print(f"  Container: {file_info['container_name']}")

    # Create NikonImageStackInterface
    interface = NikonImageStackInterface(nd2_file)

    # Get session start time from ND2 metadata
    session_start_time = interface.get_session_start_time()

    # Create BIDS-style base session ID with detailed timestamp when available
    if hasattr(session_start_time, "hour"):
        timestamp = session_start_time.strftime("%Y%m%d_%H%M%S")
    else:
        timestamp = session_start_time.strftime("%Y%m%d")

    base_session_id = f"figure4_ConfocalSpineDensity_{condition.replace(' ', '_').replace('-', '_')}_{timestamp}_Sub{file_info['animal_id']}"

    # Load metadata from YAML file
    metadata_file_path = Path(__file__).parent.parent.parent / "metadata.yaml"
    paper_metadata = load_dict_from_file(metadata_file_path)

    # Create session-specific metadata
    session_specific_metadata = {
        "NWBFile": {
            "session_description": (
                f"High-resolution confocal spine density analysis in indirect pathway spiny projection "
                f"neurons (iSPNs) for condition {condition}. Confocal laser scanning microscopy was used "
                f"to acquire Z-stack images at approximately 30 μm from the soma using a Nikon AXR system. "
                f"Acquisition parameters: 0.09 μm pixels, 0.125 μm z-steps, 60x oil immersion objective "
                f"(NA=1.49). Images were de-noised and deconvolved using Nikon NIS-Elements AR 5.41.02, "
                f"then analyzed using Imaris 10.0.0 with Labkit pixel classification and supervised "
                f"learning for automated spine detection. Manual validation was performed with thinnest "
                f"spine head set at 0.188 μm and maximum spine length of 5 μm. This high-resolution "
                f"confocal data demonstrates detection of approximately 2x more spines compared to "
                f"standard two-photon microscopy, validating methodological differences in spine counting."
            ),
            "identifier": str(uuid.uuid4()),
            "session_start_time": session_start_time,
            "session_id": f"{base_session_id}_Cell{file_info['cell_number']}_Slide{file_info['slide_number']}_Slice{file_info['slice_number']}",
            "keywords": [
                "confocal microscopy",
                "spine density",
                "dendritic spines",
                "iSPNs",
                "high-resolution imaging",
                "Imaris analysis",
                "methodology validation",
            ],
        }
    }

    # Merge paper metadata with session-specific metadata
    merged_metadata = dict_deep_update(paper_metadata, session_specific_metadata)

    # Create NWB file
    nwbfile = NWBFile(
        session_description=merged_metadata["NWBFile"]["session_description"],
        identifier=merged_metadata["NWBFile"]["identifier"],
        session_start_time=merged_metadata["NWBFile"]["session_start_time"],
        experimenter=merged_metadata["NWBFile"]["experimenter"],
        lab=merged_metadata["NWBFile"]["lab"],
        institution=merged_metadata["NWBFile"]["institution"],
        experiment_description=merged_metadata["NWBFile"]["experiment_description"],
        session_id=merged_metadata["NWBFile"]["session_id"],
        keywords=merged_metadata["NWBFile"]["keywords"],
    )

    # Create subject metadata for Figure 4 confocal spine density experiments
    subject = Subject(
        subject_id=f"confocal_mouse_{file_info['animal_id']}",
        species="Mus musculus",
        strain="Drd1-Tdtomato transgenic",
        description=(
            f"Adult Drd1-Tdtomato transgenic mouse with unilateral 6-OHDA lesion (>95% dopamine depletion) "
            f"modeling Parkinson's disease. dSPNs identified by Drd1-Tdtomato expression. "
            f"Confocal microscopy spine density analysis from {file_info['location']} dendrite "
            f"of cell {file_info['cell_number']} on slide {file_info['slide_number']}, slice {file_info['slice_number']}."
        ),
        genotype="Drd1-Tdtomato+",
        sex="M",
        age="P8W/P12W",  # Adult mice, 8-12 weeks in ISO 8601 format
    )
    nwbfile.subject = subject

    # Create metadata for the interface with custom container information
    # Following the same pattern as NeuroConv's ImageInterface
    interface_metadata = {
        "Images": {
            file_info["container_name"]: {"name": file_info["container_name"], "description": file_info["description"]}
        }
    }

    # Add to NWB file using interface with custom metadata
    interface.add_to_nwbfile(nwbfile, metadata=interface_metadata)

    if verbose:
        print(f"  Added Z-stack with {interface.number_of_z_planes} images")
        print(f"  Image dimensions: {interface.width_pixels}x{interface.height_pixels}")

    return nwbfile


if __name__ == "__main__":
    import logging

    from tqdm import tqdm

    # Control verbose output
    verbose = False  # Set to True for detailed output
    stub_test = True  # Set to True to process only first 2 files per condition for testing

    logging.getLogger("nd2").setLevel(logging.ERROR)

    # Process both Figure 4I and Figure 4J data
    datasets = [
        {
            "name": "Fig4I",
            "path": Path(
                "/home/heberto/development/surmeier-lab-to-nwb/link_to_raw_data/Figure 4_SF1B_SF5/Confocal spine density/Fig 4I/raw"
            ),
            "conditions": ["control", "6-OHDA", "off-state", "on-state"],
        },
        {
            "name": "Fig4J",
            "path": Path(
                "/home/heberto/development/surmeier-lab-to-nwb/link_to_raw_data/Figure 4_SF1B_SF5/Confocal spine density/Fig 4J and Suppl Fig 5/raw"
            ),
            "conditions": ["control", "6-OHDA", "off-state", "on-state"],
        },
    ]

    # Create output directory
    root_dir = Path(__file__).parent.parent.parent.parent.parent.parent  # Go up to repo root
    nwb_files_dir = root_dir / "nwb_files" / "confocal_spine_density" / "figure_4"
    nwb_files_dir.mkdir(parents=True, exist_ok=True)

    for dataset in datasets:
        if verbose:
            print(f"\nProcessing {dataset['name']} confocal spine density data...")

        for condition in dataset["conditions"]:
            condition_path = dataset["path"] / condition
            assert condition_path.exists(), f"Condition path {condition_path} does not exist"

            if verbose:
                print(f"Processing {dataset['name']} {condition}...")

            # Get all ND2 files in condition folder
            nd2_files = list(condition_path.glob("*.nd2"))
            nd2_files.sort()

            # Apply stub_test filtering if enabled
            if stub_test:
                nd2_files = nd2_files[:2]
                if verbose:
                    print(f"stub_test enabled: processing only first {len(nd2_files)} ND2 files")

            if verbose:
                print(f"Found {len(nd2_files)} ND2 files in {condition}")

            # Process each ND2 file with progress bar
            file_iterator = tqdm(
                nd2_files,
                desc=f"Converting Figure4 ConfocalSpineDensity {dataset['name']} {condition}",
                disable=verbose,
            )

            for nd2_file in file_iterator:
                if verbose:
                    print(f"  File: {nd2_file.name}")

                # Convert ND2 file to NWB format
                nwbfile = convert_session_to_nwbfile(nd2_file, condition, verbose=verbose)

                # Create safe filename
                nd2_name = nd2_file.stem.replace(" ", "_").replace("-", "_")
                condition_safe = condition.replace("-", "_")
                nwbfile_path = (
                    nwb_files_dir / f"figure4_confocal_{dataset['name'].lower()}_{condition_safe}_{nd2_name}.nwb"
                )

                # Write NWB file
                configure_and_write_nwbfile(nwbfile, nwbfile_path=nwbfile_path)
                if verbose:
                    print(f"Successfully saved: {nwbfile_path.name}")

    if verbose:
        print(f"\nConfocal spine density conversion completed!")
        print(f"Output directory: {nwb_files_dir}")
