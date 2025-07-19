import re
from pathlib import Path
from typing import Any, Dict

from neuroconv.tools import configure_and_write_nwbfile
from neuroconv.utils import dict_deep_update, load_dict_from_file
from pynwb import NWBFile

from surmeier_lab_to_nwb.zhai2025.image_stack_interfaces import NikonImageStackInterface


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

    # Create session identifier
    session_id = f"{animal_id}_slide{slide_num}_slice{slice_num}"

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
        "session_id": session_id,
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

    # Load metadata from YAML file
    metadata_file_path = Path(__file__).parent / "metadata.yaml"
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
            "identifier": f"zhai2025_fig4_confocal_spine_density_{file_info['session_id']}_{condition.replace(' ', '_')}",
            "session_start_time": session_start_time,
            "session_id": f"{condition}_{file_info['session_id']}",
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
    root_dir = Path(__file__).parent.parent.parent.parent  # Go up to repo root
    nwb_files_dir = root_dir / "nwb_files" / "figure_4_confocal_spine_density"
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

            if verbose:
                print(f"Found {len(nd2_files)} ND2 files in {condition}")

            # Process each ND2 file with progress bar
            file_iterator = (
                tqdm(nd2_files, desc=f"{dataset['name']} {condition}", disable=verbose) if not verbose else nd2_files
            )

            for nd2_file in file_iterator:
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
                if verbose or len(nd2_files) == 1:  # Always show success for small datasets or verbose mode
                    print(f"Successfully saved: {nwbfile_path.name}")

    if verbose:
        print(f"\nConfocal spine density conversion completed!")
        print(f"Output directory: {nwb_files_dir}")
