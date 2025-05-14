from datetime import datetime
from pathlib import Path
from typing import Literal
from uuid import uuid4
from zoneinfo import ZoneInfo

import numpy as np
import pandas as pd
from lxml import etree
from neuroconv.tools import configure_and_write_nwbfile
from pynwb import NWBFile
from pynwb.device import Device
from pynwb.file import Subject
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


def add_cell_F_I_protcol(
    nwbfile: NWBFile,
    cell_folder_path: Path,
    amplifier_device: Device,
    condition: str,
) -> list:
    """Add cell F-I protocol data to the NWB file"""

    # Input current values for each recording
    protocol_step_to_current = {
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
        "022": 320,  # pA  - some cells like  /LID off-state/02022017_1 go further than 300 pA
        "023": 340,  # pA
        "024": 360,  # pA
        "025": 380,  # pA
        "026": 400,  # pA
    }

    cell_id = cell_folder_path.name
    camel_case_condition = "".join(word.capitalize() for word in condition.replace("-", " ").split(" "))

    # Create electrode for this cell
    name = f"IntracellularElectrodeCell{camel_case_condition}{cell_id}"
    electrode_description = "Recording from dSPN in the dorsolateral striatum"
    electrode = nwbfile.create_icephys_electrode(
        name=name,
        description=electrode_description,
        device=amplifier_device,
        cell_id=cell_id,
    )

    # Process each recording (current step) for this cell
    recording_indices = []  # Store indices for IntracellularRecordingsTable
    recording_folder_path_list = [f for f in cell_folder_path.iterdir() if f.is_dir()]
    for recording_folder_path in recording_folder_path_list:
        cell_name, protocol_step = recording_folder_path.stem.split("-")
        if protocol_step not in protocol_step_to_current:
            raise ValueError(
                f"Protocol step {protocol_step} not found in protocol_step_to_current mapping. for {recording_folder_path}"
            )
        current_in_pA = protocol_step_to_current[protocol_step]
        # Form the folder path for this recording

        # XML file path
        xml_file = recording_folder_path / f"{recording_folder_path.stem}_Cycle00001_VoltageRecording_001.xml"

        if not xml_file.exists():
            raise FileNotFoundError(f"XML file {xml_file} does not exist.")

        # Parse XML to get metadata
        metadata = parse_xml(xml_file)

        # Get CSV data file path
        voltage_csv_stem = metadata["voltage recording file"]
        if voltage_csv_stem is None:
            voltage_csv_stem = xml_file.stem
        csv_file_path = recording_folder_path / f"{voltage_csv_stem}.csv"

        if not csv_file_path.exists():
            raise FileNotFoundError(f"CSV file {csv_file_path} does not exist.")

        # Read voltage recording data
        voltage_df = pd.read_csv(csv_file_path, header=0)

        # Get sampling rate
        sampling_rate = metadata["sampling"]

        # Create stimulus series (current injection)
        stimulus_name = f"CurrentClampStimulusSeries{camel_case_condition}{cell_id}{current_in_pA}pA"
        constant_current_data = np.ones(voltage_df.shape[0]) * current_in_pA
        stimulus = CurrentClampStimulusSeries(
            name=stimulus_name,
            data=constant_current_data,
            rate=sampling_rate,
            electrode=electrode,
            starting_time=np.nan,
        )

        # Get voltage data
        voltage_data = voltage_df[" Primary"].values
        response_name = f"CurrentClampSeries{camel_case_condition}{cell_id}{current_in_pA}pA"

        # Create response series (recorded voltage)
        response = CurrentClampSeries(
            name=response_name,
            data=voltage_data,
            rate=sampling_rate,
            electrode=electrode,
            starting_time=np.nan,
        )

        # Add an intracellular recording entry with a unique ID
        recording_index = nwbfile.add_intracellular_recording(
            electrode=electrode,
            response=response,
            stimulus_current_pA=current_in_pA,
        )

        recording_indices.append(recording_index)
        print(f"  Added recording step {protocol_step} with current {current_in_pA} pA")

    return recording_indices


def convert_data_to_nwb(folder_path: Path, nwbfile_path: Path):
    """
    Convert all cells from a condition to NWB format, with proper organizational structure

    Parameters:
    -----------
    folder_path : Path
        Path to the directory containing all cell recordings
    nwbfile_path : Path
        Path to the output NWB file
    """

    # Create a new NWB file

    # Process just one condition as requested
    session_id = f"SomaticExcitability"
    session_start_time = datetime.now(ZoneInfo("America/Chicago"))
    experiment_description = f"State-dependent modulation of spiny projection neurons in LID model (Figure 1)"

    nwbfile = NWBFile(
        session_description=f"Somatic excitability",
        identifier=str(uuid4()),
        session_start_time=session_start_time,
        experimenter="Zhai et al.",
        lab="Surmeier Lab",
        institution="Northwestern University",
        experiment_description=experiment_description,
        session_id=session_id,
    )

    # Create a subject with the information we have and notes about missing data
    subject = Subject(
        subject_id="Not a single subject",
        species="Mus musculus",
        strain="C57Bl/6",  # Common background strain for the transgenic mice mentioned
        description="Experimental mice with unilateral 6-OHDA lesions in the medial forebrain bundle. Used Drd1-Tdtomato BAC transgenic mice for dSPN experiments. Total number of unique subjects is not specified in the data.",
        genotype="Drd1-Tdtomato bacterial artificial chromosome (BAC) transgenic",
        sex="M",
        age="P49-P84",  # ISO format for 7-12 weeks (postnatal days)
    )

    nwbfile.subject = subject

    # Add devices
    amplifier_device = nwbfile.create_device(
        name="MultiClamp 700B",
        description="MultiClamp 700B amplifier (Axon Instruments/Molecular Devices) with signals filtered at 2 kHz and digitized at 10 kHz",
        manufacturer="Axon Instruments/Molecular Devices",
    )

    # Add a stimulus column to the NWB file
    intracellular_recording_table = nwbfile.get_intracellular_recordings()
    intracellular_recording_table.add_column(
        name="stimulus_current_pA",
        description="The current stimulus applied during the recording in picoamps",
    )

    intracellular_repetitions_table = nwbfile.get_icephys_repetitions()
    intracellular_repetitions_table.add_column(
        name="condition",
        description="Experimental condition for the intracellular recording",
    )

    intracellular_experimental_condition_table = nwbfile.get_icephys_experimental_conditions()
    intracellular_experimental_condition_table.add_column(
        name="condition",
        description="Experimental condition for the intracellular recording",
    )

    condition: Literal["LID on-state", "LID off-state", "LID on-state with SCH"] = "LID on-state"
    conditions = ["LID on-state", "LID off-state", "LID on-state with SCH"]
    for condition in conditions:
        print("" + "-" * 50)
        print(f"\nProcessing condition: {condition}")
        print("" + "-" * 50)
        condition_folder = folder_path / condition

        # Make sure the path exists
        if not condition_folder.exists():
            raise FileNotFoundError(f"Condition folder {condition_folder} does not exist.")

        # Find all cell folders in the condition directory
        cell_folders = [f for f in condition_folder.iterdir() if f.is_dir()]

        if not cell_folders:
            raise ValueError(f"No cell folders found in {condition_folder}")

        print(f"Found {len(cell_folders)} cells in {condition} condition")

        # Process each cell (each cell is a repetition of the protocol)
        sequential_recording_indices = []
        for cell_folder_path in sorted(cell_folders):
            cell_id = cell_folder_path.name
            print(f"\nProcessing cell {cell_id}...")

            recording_indices = add_cell_F_I_protcol(
                nwbfile=nwbfile,
                cell_folder_path=cell_folder_path,
                amplifier_device=amplifier_device,
                condition=condition,
            )

            # Add a simultaneous recording entry (in this case, just one recording per simultaneous recording)
            # There a no simultaneous recordings in this dataset, but we need to add this entry for the NWB structure
            simultaneous_recording_indices = []
            for recording_index in recording_indices:
                simultaneous_index = nwbfile.add_icephys_simultaneous_recording(recordings=[recording_index])
                simultaneous_recording_indices.append(simultaneous_index)

            # Add a sequential recording entry for this cell (grouping all current steps)
            sequential_index = nwbfile.add_icephys_sequential_recording(
                simultaneous_recordings=simultaneous_recording_indices,
                stimulus_type="F-I protocol",
            )

            sequential_recording_indices.append(sequential_index)
            print(f"  Added sequential recording for cell {cell_id} with condition {condition}")

        # Add this cell as a repetition
        repetition_index = nwbfile.add_icephys_repetition(
            sequential_recordings=sequential_recording_indices, condition=condition
        )
        nwbfile.add_icephys_experimental_condition(repetitions=[repetition_index], condition=condition)
        print(f"Added intracellular data for condition {condition} to the NWB file")

    # Write NWB file
    configure_and_write_nwbfile(nwbfile=nwbfile, nwbfile_path=nwbfile_path)

    print(f"NWB file written to {nwbfile_path.resolve()}")

    return nwbfile_path


if __name__ == "__main__":
    # Set the base path to your data
    raw_data_folder = Path(
        "/media/heberto/One Touch/Surmeier-CN-data-share/consolidated_data/LID_paper_Zhai_2025/Raw data for Figs"
    )
    figure_1_folder = raw_data_folder / "Figure 1"
    somatic_excitability_folder = figure_1_folder / "Somatic excitability"

    # Define output file
    nwbfile_path = Path("figure1_somatic_excitability_LID_on_state.nwb")

    # Convert data to NWB
    convert_data_to_nwb(somatic_excitability_folder, nwbfile_path)
