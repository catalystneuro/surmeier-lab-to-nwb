from datetime import datetime
from pathlib import Path
from typing import Literal
from uuid import uuid4
from zoneinfo import ZoneInfo

import numpy as np
import pandas as pd
from lxml import etree

# Import pynwb
from pynwb import NWBHDF5IO, NWBFile
from pynwb.icephys import CurrentClampSeries, CurrentClampStimulusSeries


def _get_ephys_vals(element):
    ch_type = element.find(".//PatchclampChannel").text

    if ch_type == "0":
        unit = element.find(".//UnitName").text
        divisor = float(element.find(".//Divisor").text)

        return "primary", {"unit": unit, "divisor": divisor}

    elif ch_type == "1":
        unit = element.find(".//UnitName").text
        divisor = float(element.find(".//Divisor").text)

        return "secondary", {"unit": unit, "divisor": divisor}


def parse_xml(filename):
    """Parse XML metadata from Prairie View electrophysiology files"""
    tree = etree.parse(filename)
    # find all elements associated with enabled channels
    enabled_ch = tree.xpath('.//Enabled[text()="true"]')

    file_attr = {}
    ch_names = []
    for ch in enabled_ch:
        recorded_signal = ch.getparent()
        if recorded_signal.find(".//Type").text == "Physical":
            clamp_device = recorded_signal.find(".//PatchclampDevice").text

            if clamp_device is not None:
                name, ephys_vals = _get_ephys_vals(recorded_signal)
                file_attr[name] = ephys_vals

            else:
                name = recorded_signal.find(".//Name").text

            ch_names.append(name.capitalize())

    file_attr["channels"] = ch_names
    # gets sampling rate
    file_attr["sampling"] = float((tree.find(".//Rate")).text)
    # gets recording time, converts to sec
    file_attr["duration"] = int((tree.find(".//AcquisitionTime")).text)

    # finds the voltage recording csv file name
    datafile = (tree.find(".//DataFile")).text
    # finds the linescan profile file name (if doesn't exist, will be None)
    ls_file = (tree.find(".//AssociatedLinescanProfileFile")).text

    # If ls_file is none this could mean that there is no linescan associated
    # with that voltage recording file or that the file passed to parse_vr is
    # actually a LineScan data file and therefore should be passed to ls_file.
    # In that scenario there is no voltage recording file, so vo_file is None
    if ls_file is None:
        if "LineScan" in datafile:
            ls_file = datafile
            vo_file = None
        elif "LineScan" not in datafile:
            vo_file = datafile
    else:
        vo_file = datafile

    file_attr["voltage recording file"] = vo_file
    file_attr["linescan file"] = ls_file

    return file_attr


def read_voltage_data_as_data_frame(file_path) -> pd.DataFrame:
    """Read the voltage data from CSV file"""
    # Adjust column specifications based on your actual data format
    df = pd.read_csv(file_path, header=0)
    return df


def convert_to_nwb(cell_folder_path: str | Path, output_file: str | Path):
    """
    Convert all recordings to NWB format

    Parameters:
    -----------
    cell_folder_path : Path
        Path to the directory containing the cell recordings
    output_file : Path
        Name of the output NWB file

    """
    # Input current values for each recording
    file_termination_to_current = {
        "001": -120,  # pA
        "002": -100,  # pA
        "003": -80,  # pA
        "004": -60,  # pA
        "005": -40,  # pA
        "006": -20,  # pA
        "007": 20,  # pA
        "008": 40,  # pA
        "009": 60,  # pA
        "010": 80,  # pA
        "011": 100,  # pA
        "012": 120,  # pA
        "013": 140,  # pA
        "014": 160,  # pA
        "015": 180,  # pA
        "016": 200,  # pA
        "017": 220,  # pA
        "018": 240,  # pA
        "019": 260,  # pA
        "020": 280,  # pA
        "021": 300,  # pA
    }

    # Create a new NWB file
    session_id = "CellRecording"
    session_start_time = datetime.now(ZoneInfo("America/Chicago"))
    experiment_description = "State-dependent modulation of spiny projection neurons in LID model (Figure 1)"

    nwbfile = NWBFile(
        session_description="Experiment on spiny projection neurons",
        identifier=str(uuid4()),
        session_start_time=session_start_time,
        experimenter="Zahi et al.",
        lab="Surmeier Lab",
        institution="Northwestern University",
        experiment_description=experiment_description,
        session_id=session_id,
    )

    # Add devices
    amplifier_device = nwbfile.create_device(
        name="MultiClamp 700B",
        description="MultiClamp 700B amplifier (Axon Instruments/Molecular Devices) with PCI-NI6713 analog-to-digital converter, signals filtered at 2 kHz and digitized at 10 kHz",
        manufacturer="Axon Instruments/Molecular Devices",
    )

    digitizer_device = nwbfile.create_device(
        name="PCI-NI6713",
        description="Analog-to-digital converter card, digitizing at 10 kHz",
        manufacturer="National Instruments",
    )

    # Create electrode
    cell_id = cell_folder_path.stem
    name = f"IntracellularElectrode{cell_id}"
    electrode_description = "Recording from dSPN in the dorsolateral striatum"
    electrode = nwbfile.create_icephys_electrode(
        name=name,
        description=electrode_description,
        device=amplifier_device,
        cell_id=cell_id,
    )

    # Process each recording
    recording_indices = []
    for recording_num, input_current in file_termination_to_current.items():
        # Form the folder path for this recording
        recording_folder = cell_folder_path / f"cell1-{recording_num}"

        if not recording_folder.exists():
            raise FileNotFoundError(f"Recording folder {recording_folder} does not exist.")

        # XML file path
        xml_file = recording_folder / f"cell1-{recording_num}_Cycle00001_VoltageRecording_001.xml"

        if not xml_file.exists():
            raise FileNotFoundError(f"XML file {xml_file} does not exist.")

        # Parse XML to get metadata
        metadata = parse_xml(xml_file)

        # Get CSV data file path
        voltage_csv_stem = metadata["voltage recording file"]
        if voltage_csv_stem is None:
            voltage_csv_stem = xml_file.stem
        csv_file_path = recording_folder / f"{voltage_csv_stem}.csv"
        if not csv_file_path.exists():
            raise FileNotFoundError(f"CSV file {csv_file_path} does not exist.")

        # Read voltage recording data
        voltage_df = read_voltage_data_as_data_frame(csv_file_path)

        # Get sampling rate
        sampling_rate = metadata["sampling"]

        # Create stimulus series (current injection)
        constant_current_data = voltage_df.shape[0] * [input_current]  # Create a constant current array
        stimulus = CurrentClampStimulusSeries(
            name=f"current_stimulus_{recording_num}",
            data=constant_current_data,  # Use the constant current array
            rate=sampling_rate,
            electrode=electrode,
            starting_time=np.nan,
        )
        voltage_data = voltage_df[" Primary"].values
        # Create response series (recorded voltage)
        response = CurrentClampSeries(
            name=f"voltage_response_{recording_num}",
            data=voltage_data,  # Use the updated voltage data
            rate=sampling_rate,
            electrode=electrode,
            starting_time=np.nan,
        )

        # Add to NWB file
        nwbfile.add_stimulus(stimulus)
        nwbfile.add_acquisition(response)

        # Add an intracellular recording entry
        recording_index = nwbfile.add_intracellular_recording(
            electrode=electrode, stimulus=stimulus, response=response, id=int(recording_num)
        )

        recording_indices.append(recording_index)
        print(f"Added recording {recording_num} with current {input_current} pA")

    # Write NWB file
    with NWBHDF5IO(output_file, "w") as io:
        io.write(nwbfile)

    print(f"NWB file written to {output_file}")
    return output_file


if __name__ == "__main__":
    # Set the base path to your data
    figure_1_folder = Path(
        "/media/heberto/One Touch/Surmeier-CN-data-share/consolidated_data/LID_paper_Zhai_2025/Raw data for Figs/Figure 1"
    )
    somatic_excitability_folder = figure_1_folder / "Somatic excitability"
    condition: Literal["LID on-state", "LID off-state", "LID on-state with SCH"] = "LID on-state"
    condition_folder = somatic_excitability_folder / f"{condition}"
    cell_id = "04112019_1"
    cell_folder_path = condition_folder / f"{cell_id}"

    # Make sure the path exists
    if not cell_folder_path.exists():
        raise FileNotFoundError(f"Base path {cell_folder_path} does not exist.")

    # Convert data to NWB
    output_file = Path("LID_on_state_cell1_dSPN.nwb").resolve()
    convert_to_nwb(cell_folder_path, output_file)
