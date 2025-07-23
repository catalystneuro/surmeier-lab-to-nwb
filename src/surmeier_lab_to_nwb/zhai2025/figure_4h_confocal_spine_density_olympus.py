import re
import uuid
from datetime import datetime
from pathlib import Path
from typing import Any, Dict

from neuroconv.datainterfaces import ImageInterface
from neuroconv.tools import configure_and_write_nwbfile
from neuroconv.utils import dict_deep_update, load_dict_from_file
from pynwb import NWBFile
from pynwb.device import Device


def parse_pty_file(pty_file: Path, verbose: bool = False) -> Dict[str, Any]:
    """
    Parse Olympus PTY parameter file to extract acquisition metadata.

    Parameters
    ----------
    pty_file : Path
        Path to PTY parameter file
    verbose : bool, default=False
        Enable verbose output

    Returns
    -------
    Dict[str, Any]
        Dictionary containing extracted parameters
    """
    if verbose:
        print(f"    Parsing PTY file: {pty_file.name}")

    metadata = {}

    # PTY files are Unicode text files with null characters
    with open(pty_file, "r", encoding="utf-16le", errors="ignore") as f:
        content = f.read()

    # Remove null characters and clean up
    content = content.replace("\x00", "")

    # Parse key acquisition parameters
    patterns = {
        "confocal_mode": r'Confocal="([^"]+)"',
        "magnification": r"Magnification=([0-9.]+)",
        "numerical_aperture": r"ObjectiveLens NAValue=([0-9.]+)",
        "objective_name": r'ObjectiveLens Name="([^"]+)"',
        "observation_mode": r'Observation Mode="([^"]+)"',
        "pmt_voltage": r"PMTVoltage=([0-9]+)",
        "pinhole_diameter": r"PinholeDiameter=([0-9]+)",
        "scan_speed": r"ScanSpeed=([0-9.]+)",
        "pixel_size_x": r"WidthConvertValue=([0-9.]+)",
        "pixel_size_y": r"HeightConvertValue=([0-9.]+)",
        "image_width": r"ImageWidth=([0-9]+)",
        "image_height": r"ImageHeight=([0-9]+)",
        "data_type": r'DataType="([^"]+)"',
        "bit_depth": r"ValidBitCounts=([0-9]+)",
        "wavelength": r"AbsPositionValue=([0-9.]+).*nm",
    }

    for key, pattern in patterns.items():
        match = re.search(pattern, content)
        if match:
            value = match.group(1)
            # Convert numeric values
            if key in [
                "magnification",
                "numerical_aperture",
                "scan_speed",
                "pixel_size_x",
                "pixel_size_y",
                "wavelength",
            ]:
                metadata[key] = float(value)
            elif key in ["pmt_voltage", "pinhole_diameter", "image_width", "image_height", "bit_depth"]:
                metadata[key] = int(value)
            else:
                metadata[key] = value

    if verbose:
        print(f"      Extracted {len(metadata)} parameters")

    return metadata


def extract_z_info_from_filename(tiff_file: Path) -> Dict[str, Any]:
    """
    Extract Z-slice information from TIFF filename.

    Example: s_C001Z025.tif -> Z-slice 25

    Parameters
    ----------
    tiff_file : Path
        Path to TIFF file

    Returns
    -------
    Dict[str, Any]
        Dictionary with Z-slice information
    """
    filename = tiff_file.name

    # Extract Z-slice number from filename pattern s_C001Z###.tif
    z_match = re.search(r"Z(\d+)\.tif", filename)
    if not z_match:
        raise ValueError(f"Could not extract Z-slice number from filename: {filename}")

    z_slice_num = int(z_match.group(1))

    # Olympus confocal typically uses smaller z-steps
    z_step_um = 0.125  # From conversion notes: 0.125 μm z-steps
    z_depth_um = (z_slice_num - 1) * z_step_um  # Z001 is depth 0

    return {
        "z_slice_number": z_slice_num,
        "z_depth_um": z_depth_um,
        "image_name": f"ImageZ{z_slice_num:03d}",
        "description": f"Z-slice {z_slice_num} at depth {z_depth_um:.3f} μm. "
        f"Acquired with Olympus confocal microscopy at 0.125 μm z-steps "
        f"for high-resolution spine density analysis.",
    }


def parse_oif_stack_info(oif_folder: Path) -> Dict[str, Any]:
    """
    Parse OIF stack folder to extract stack information.

    Parameters
    ----------
    oif_folder : Path
        Path to OIF.files folder containing TIFF slices

    Returns
    -------
    Dict[str, Any]
        Dictionary with stack information
    """
    folder_name = oif_folder.name

    # Extract image numbers from folder name like "Image_01_01_02_01.oif.files"
    image_match = re.search(r"Image_(\d+)_(\d+)_(\d+)_(\d+)", folder_name)
    if image_match:
        stack_id = f"Stack_{image_match.group(3)}"  # Use the third number as stack identifier
    else:
        stack_id = f"Stack_{folder_name}"

    return {
        "stack_id": stack_id,
        "container_name": f"Images{stack_id}",
        "description": (
            f"High-resolution confocal Z-stack {stack_id} from Olympus FV10i system. "
            f"Acquired at 60x magnification (NA=1.35) with 0.125 μm z-steps and "
            f"~0.207 μm pixels for spine density methodology validation. "
            f"Part of Figure 4H demonstrating confocal vs two-photon resolution differences."
        ),
    }


def create_olympus_confocal_device(nwbfile: NWBFile, pty_metadata: Dict[str, Any]) -> Device:
    """
    Create Olympus confocal microscope device.

    Parameters
    ----------
    nwbfile : NWBFile
        NWB file to add device to
    pty_metadata : Dict[str, Any]
        Metadata from PTY file

    Returns
    -------
    Device
        Device object for Olympus confocal microscope
    """
    objective_name = pty_metadata.get("objective_name", "UPLSAP60xO")
    magnification = pty_metadata.get("magnification", 60.0)
    numerical_aperture = pty_metadata.get("numerical_aperture", 1.35)
    confocal_mode = pty_metadata.get("confocal_mode", "ON")

    device_description = (
        f"Olympus FV10i confocal laser scanning microscope with {objective_name} objective "
        f"(NA={numerical_aperture}, {magnification}x magnification). "
        f"Confocal mode: {confocal_mode}. Used for high-resolution spine density imaging "
        f"with 0.125 μm z-steps and ~0.207 μm pixels. System optimized for detecting "
        f"approximately 2x more spines compared to standard two-photon microscopy."
    )

    device = Device(name="OlympusFV10iConfocal", description=device_description, manufacturer="Olympus")

    nwbfile.add_device(device)
    return device


def convert_session_to_nwbfile(session_folder: Path, verbose: bool = False) -> NWBFile:
    """
    Convert Figure 4H Olympus confocal session to NWB format.

    Parameters
    ----------
    session_folder : Path
        Path to session folder (e.g., "5104-3 60x str_Cycle")
    verbose : bool, default=False
        Enable verbose output

    Returns
    -------
    NWBFile
        The populated NWB file object
    """
    if verbose:
        print(f"Processing session: {session_folder.name}")

    # Get all OIF files in session
    oif_files = sorted(session_folder.glob("*.oif"))

    if not oif_files:
        raise ValueError(f"No OIF files found in {session_folder}")

    if verbose:
        print(f"Found {len(oif_files)} OIF files")

    # Parse session info from folder name
    session_name = session_folder.name
    # Extract animal ID (first part before space or dash)
    animal_id = session_name.split()[0] if " " in session_name else session_name.split("-")[0]

    # Use current date as session start time (OIF files don't contain reliable timestamps)
    session_start_time = datetime.now().astimezone()

    # Create BIDS-style base session ID with detailed timestamp when available
    if hasattr(session_start_time, "hour"):
        timestamp = session_start_time.strftime("%Y%m%d_%H%M%S")
    else:
        timestamp = session_start_time.strftime("%Y%m%d")

    base_session_id = f"figure4h_ConfocalSpineDensityOlympus_methodological_validation_{timestamp}"
    script_specific_id = f"Sub{animal_id}"
    session_id = f"{base_session_id}_{script_specific_id}"

    # Load metadata from YAML file
    metadata_file_path = Path(__file__).parent / "metadata.yaml"
    paper_metadata = load_dict_from_file(metadata_file_path)

    # Create session-specific metadata
    session_specific_metadata = {
        "NWBFile": {
            "session_description": (
                f"High-resolution confocal spine density imaging using Olympus FV10i system "
                f"for methodological validation study. Multiple dendritic segments were imaged "
                f"at 60x magnification (NA=1.35) with 0.125 μm z-steps and ~0.207 μm pixels. "
                f"Data demonstrates that confocal microscopy detects approximately 2x more spines "
                f"compared to standard two-photon microscopy due to superior resolution. "
                f"Acquisition includes {len(oif_files)} image stacks from animal {animal_id} "
                f"processed for Figure 4H spine density comparison analysis."
            ),
            "identifier": str(uuid.uuid4()),
            "session_start_time": session_start_time,
            "session_id": session_id,
            "keywords": [
                "confocal microscopy",
                "spine density",
                "Olympus FV10i",
                "methodological validation",
                "high-resolution imaging",
                "dendritic spines",
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

    # Process each OIF file
    device_created = False

    for oif_file in oif_files:
        if verbose:
            print(f"  Processing OIF: {oif_file.name}")

        # Get corresponding .oif.files folder
        oif_files_folder = session_folder / f"{oif_file.stem}.oif.files"

        if not oif_files_folder.exists():
            raise FileNotFoundError(f"Required .oif.files folder not found: {oif_files_folder}")

        # Get all TIFF files in the folder
        tiff_files = sorted(oif_files_folder.glob("s_C001Z*.tif"))

        if not tiff_files:
            raise FileNotFoundError(f"No TIFF files found in {oif_files_folder}")

        if verbose:
            print(f"    Found {len(tiff_files)} TIFF files")

        # Parse stack information
        stack_info = parse_oif_stack_info(oif_files_folder)

        # Get metadata from first PTY file
        first_pty = oif_files_folder / f"{tiff_files[0].stem}.pty"
        if not first_pty.exists():
            raise FileNotFoundError(f"Required PTY metadata file not found: {first_pty}")

        pty_metadata = parse_pty_file(first_pty, verbose=verbose)

        # Create device (only once)
        if not device_created:
            device = create_olympus_confocal_device(nwbfile, pty_metadata)
            device_created = True
            if verbose:
                print(f"    Created device: {device.name}")

        # Create ImageInterface for this stack
        interface = ImageInterface(file_paths=tiff_files, images_container_metadata_key=stack_info["container_name"])

        # Get base metadata from interface
        interface_metadata = interface.get_metadata()

        # Create custom metadata for individual images
        images_metadata_dict = {}
        for tiff_file in tiff_files:
            z_info = extract_z_info_from_filename(tiff_file)
            file_path_str = str(tiff_file)

            # Calculate resolution from PTY metadata
            pixel_size_um = pty_metadata.get("pixel_size_x", 0.207)  # Default from PTY files
            resolution = 10000.0 / pixel_size_um  # Convert μm/pixel to pixels/cm

            images_metadata_dict[file_path_str] = {
                "name": z_info["image_name"],
                "description": z_info["description"],
                "resolution": resolution,
            }

        # Update interface metadata with custom container info and individual images
        interface_metadata["Images"][stack_info["container_name"]].update(
            {
                "name": stack_info["container_name"],
                "description": stack_info["description"],
                "images": images_metadata_dict,
            }
        )

        # Add to NWB file
        interface.add_to_nwbfile(nwbfile, metadata=interface_metadata)

        if verbose:
            print(f"    Added stack: {stack_info['container_name']} with {len(tiff_files)} images")

    if verbose:
        print(f"Conversion completed for session: {session_name}")

    return nwbfile


if __name__ == "__main__":
    import logging

    from tqdm import tqdm

    # Control verbose output
    verbose = False  # Set to True for detailed output
    stub_test = True  # Set to True to process only first 2 files per condition for testing

    logging.getLogger("PIL").setLevel(logging.ERROR)

    # Define the path to Figure 4H data
    fig4h_path = Path(
        "/home/heberto/development/surmeier-lab-to-nwb/link_to_raw_data/Figure 4_SF1B_SF5/Confocal spine density/Fig 4H/raw"
    )

    # Create output directory
    root_dir = Path(__file__).parent.parent.parent.parent  # Go up to repo root
    nwb_files_dir = root_dir / "nwb_files" / "figure_4" / "confocal_spine_density_olympus"
    nwb_files_dir.mkdir(parents=True, exist_ok=True)

    # Get all session folders
    session_folders = [f for f in fig4h_path.iterdir() if f.is_dir()]
    session_folders.sort()

    # Apply stub_test filtering if enabled
    if stub_test:
        session_folders = session_folders[:2]
        if verbose:
            print(f"stub_test enabled: processing only first {len(session_folders)} session folders")

    if verbose:
        print(f"Found {len(session_folders)} session folders")

    # Process each session folder
    session_iterator = (
        tqdm(
            session_folders,
            desc="Converting sessions from figure_4h_confocal_spine_density_olympus to NWB",
            disable=verbose,
        )
        if not verbose
        else session_folders
    )

    for session_folder in session_iterator:
        if verbose:
            print(f"\nProcessing session: {session_folder.name}")
        elif not verbose:
            session_iterator.set_description(f"Processing {session_folder.name}")

        # Convert session to NWB format
        nwbfile = convert_session_to_nwbfile(session_folder, verbose=verbose)

        # Create output filename
        session_safe = session_folder.name.replace(" ", "_").replace("-", "_")
        nwbfile_path = nwb_files_dir / f"figure4h_olympus_confocal_{session_safe}.nwb"

        # Write NWB file
        configure_and_write_nwbfile(nwbfile, nwbfile_path=nwbfile_path)

        if verbose:
            print(f"Successfully saved: {nwbfile_path.name}")

    if verbose:
        print(f"\nFigure 4H Olympus confocal conversion completed!")
        print(f"Output directory: {nwb_files_dir}")
