"""GRAB-ACh3.0 sessions (Figure 5): one NWB file per mouse and recording day.

Two-photon movies of the acetylcholine sensor in striatal slices, evoked by electrical stimulation, one movie per
region and zoom with the trials joined in recording order, and a trials table with each trial's region, drug, stimulus
and frames.
"""

from __future__ import annotations

import uuid
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
import pandas as pd
from neuroconv.tools.nwb_helpers import (
    configure_and_write_nwbfile,
    make_nwbfile_from_metadata,
)
from pynwb.base import TimeSeriesReference, TimeSeriesReferenceVectorData
from pynwb.device import Device
from pynwb.image import Image

from surmeier_lab_to_nwb.zhai2025.devices import get_or_create_bruker_ultima_model
from surmeier_lab_to_nwb.zhai2025.metadata import (
    device_metadata,
    session_metadata,
    shared_metadata,
)
from surmeier_lab_to_nwb.zhai2025.mouse_day_merger import MouseDayBundle, SweepRecord
from surmeier_lab_to_nwb.zhai2025.patch_clamp_sessions import (
    _CONDITION_LABEL,
    TIMEZONE,
    _build_session_id,
)

_GRAB_ACH = session_metadata("grab_ach")
# Stimulus words in the recording names, misspellings included; the stimulus is read from the VoltageOutput file.
_GRAB_ACH_STIMULUS_WORDS = {"single", "sinlge", "singl", "burst", "20Hz"}
# Prairie View channel number -> channel name in the file, and what the channel records.
_GRAB_ACH_CHANNELS = {number: channel["name"] for number, channel in _GRAB_ACH["channels"].items()}
_GRAB_ACH_CHANNEL_DESCRIPTION = {channel["name"]: channel["description"] for channel in _GRAB_ACH["channels"].values()}


@dataclass
class SensorTrial:
    """One GRAB-ACh recording, read from its Prairie View files. Times are seconds from the session start."""

    sweep: SweepRecord
    start_time: float
    frame_times: np.ndarray
    frame_period: float
    frames: dict[str, np.ndarray]  # channel -> (frames, x, y)
    brightness: dict[str, np.ndarray]  # channel -> Prairie View's brightness of the whole frame, one value per frame
    references: dict[str, np.ndarray]  # channel -> Prairie View's 16-bit reference image, (x, y)
    zoom: float
    microns_per_pixel: float
    drug: str
    stimulation: str
    stimulus_start_time: float  # NaN without stimulation


def _read_grab_ach_stimulus(voltage_output_path: Path) -> tuple[int, float, float, float]:
    """Pulse count, width (ms), period (ms) and first-pulse delay (ms) of the stimulus in a VoltageOutput file.

    The stimulator is triggered from the output Prairie View calls LED, with 5 V pulses; burst recordings also carry a
    0 V pulse train there, which sends nothing.
    """
    root = ET.parse(voltage_output_path).getroot()
    waveforms = [w for w in root.findall("Waveform") if w.findtext("Enabled") == "true"]
    if [w.findtext("Name") for w in waveforms] != ["LED"]:
        raise ValueError(f"Expected one enabled output, LED, in {voltage_output_path}")
    trains = [
        c
        for c in waveforms[0]
        if c.tag.startswith("WaveformComponent") and float(c.findtext("PulsePotentialStart", "0")) != 0
    ]
    if len(trains) != 1 or trains[0].tag != "WaveformComponent_PulseTrain" or trains[0].findtext("Repetitions") != "1":
        raise ValueError(f"Expected one pulse train on the LED output of {voltage_output_path}")
    train = trains[0]
    width = float(train.findtext("PulseWidth"))
    return (
        int(train.findtext("PulseCount")),
        width,
        width + float(train.findtext("PulseSpacing")),
        float(train.findtext("FirstPulseDelay")),
    )


def _read_grab_ach_trial(sweep: SweepRecord, session_start_time: datetime) -> SensorTrial:
    """Read the frames, brightness, reference images and stimulus of one GRAB-ACh recording."""
    import tifffile

    leaf = sweep.recording_path
    xml_path = leaf / f"{leaf.name}.xml"
    root = ET.parse(xml_path).getroot()
    state = {e.get("key"): e for e in root.find("PVStateShard").findall("PVStateValue")}
    width = int(state["pixelsPerLine"].get("value"))
    height = int(state["linesPerFrame"].get("value"))
    microns = {v.get("index"): float(v.get("value")) for v in state["micronsPerPixel"].findall("IndexedValue")}
    if microns["XAxis"] != microns["YAxis"]:
        raise ValueError(f"Pixels are not square in {xml_path}")
    sequences = root.findall("Sequence")
    if len(sequences) != 1 or sequences[0].get("type") != "BrightnessOverTime":
        raise ValueError(f"Expected one brightness-over-time sequence in {xml_path}")
    sequence = sequences[0]

    # The sequence start, to the microsecond; the PVScan date holds the same moment truncated to the second.
    hours, minutes, seconds = sequence.get("time").split(":")
    started = datetime.combine(sweep.pv_scan_date.date(), datetime.min.time()) + timedelta(
        hours=int(hours), minutes=int(minutes), seconds=float(seconds)
    )
    if not timedelta(0) <= started - sweep.pv_scan_date < timedelta(seconds=1):
        raise ValueError(f"The sequence time of {xml_path} does not match its date")
    start_time = (started - session_start_time).total_seconds()

    frame_elements = sequence.findall("Frame")
    files = {"2": [], "3": []}
    for frame in frame_elements:
        names = {f.get("channel"): (f.get("channelName"), f.get("filename")) for f in frame.findall("File")}
        if names.keys() != files.keys() or names["2"][0] != "Ch2" or names["3"][0] != "Dodt":
            raise ValueError(f"Expected channels Ch2 and Dodt in every frame of {xml_path}")
        for channel, (_name, filename) in names.items():
            files[channel].append(filename)
    relative_times = np.array([float(frame.get("absoluteTime")) for frame in frame_elements])
    frames = {}
    for channel, filenames in files.items():
        # is_ome=False reads only this file's frame; the first file's OME header would pull in the whole series.
        stack = np.stack([tifffile.imread(leaf / filename, is_ome=False) for filename in filenames])
        if stack.shape[1:] != (height, width):
            raise ValueError(f"Frames of {leaf} are {stack.shape[1:]}, the XML says {(height, width)}")
        frames[_GRAB_ACH_CHANNELS[channel]] = stack.transpose(0, 2, 1)

    # Prairie View's brightness-over-time regions: each the whole frame, one per channel, in botData.csv.
    bots = sequence.find("PVBOTs")
    regions = bots.findall("Region")
    if sorted(r.get("channel") for r in regions) != ["2", "3"] or any(
        (r.get("x"), r.get("y"), r.get("width"), r.get("height")) != ("0", "0", str(width), str(height))
        for r in regions
    ):
        raise ValueError(f"Expected one whole-frame brightness region per channel in {xml_path}")
    table = pd.read_csv(leaf / bots.get("botData"))
    if len(table) != len(relative_times) or not np.allclose(table["Timestamp"], relative_times, atol=1e-6):
        raise ValueError(f"The brightness samples of {leaf} do not match its frames")
    brightness = {
        _GRAB_ACH_CHANNELS[r.get("channel")]: table[f"Region {number}"].to_numpy(dtype=np.float32)
        for number, r in enumerate(regions, start=1)
    }

    references = {
        name: tifffile.imread(leaf / "References" / f"{leaf.name}-Cycle00001-Ch{channel}-16bit-Reference.tif").T
        for channel, name in _GRAB_ACH_CHANNELS.items()
    }

    words = sweep.suffix.rsplit("_", 1)
    drug_words = words[0] if len(words) == 2 and words[1] in _GRAB_ACH_STIMULUS_WORDS else sweep.suffix
    if drug_words not in _GRAB_ACH["drug_by_name"]:
        raise ValueError(f"Unknown drug in the name of {leaf}")

    voltage_output = sequence.find("VoltageOutput")
    if voltage_output is None:
        stimulation, stimulus_start_time = "none", float("nan")
    else:
        count, pulse_width, period, delay = _read_grab_ach_stimulus(leaf / voltage_output.get("filename"))
        if count == 1:
            stimulation = f"single {pulse_width:g} ms pulse"
        else:
            stimulation = f"{count} pulses of {pulse_width:g} ms at {1000 / period:g} Hz"
        stimulus_start_time = start_time + float(voltage_output.get("absoluteTime")) + delay * 1e-3

    return SensorTrial(
        sweep=sweep,
        start_time=start_time,
        frame_times=start_time + relative_times,
        frame_period=float(state["framePeriod"].get("value")),
        frames=frames,
        brightness=brightness,
        references=references,
        zoom=float(state["opticalZoom"].get("value")),
        microns_per_pixel=microns["XAxis"],
        drug=_GRAB_ACH["drug_by_name"][drug_words],
        stimulation=stimulation,
        stimulus_start_time=stimulus_start_time,
    )


def merge_grab_ach_mouse_day(bundle: MouseDayBundle, output_path: Path) -> Path:
    """Build one GRAB-ACh NWB file for a mouse-day (Figure 5) and write it to disk."""
    from pynwb.image import Images
    from pynwb.ophys import (
        Fluorescence,
        ImageSegmentation,
        OpticalChannel,
        RoiResponseSeries,
        TwoPhotonSeries,
    )

    session = session_metadata("grab_ach")
    ophys_text = session["Ophys"]
    session_start_time_naive = bundle.session_start_time
    subject_entry = bundle.subject_entry
    category = subject_entry.category
    conditions = {_CONDITION_LABEL[sw.condition_subfolder] for sweeps in bundle.cells.values() for sw in sweeps}
    if len(conditions) != 1:
        raise ValueError(f"{subject_entry.subject_id} has GRAB-ACh recordings under {sorted(conditions)}")
    (condition,) = conditions
    regions = bundle.cell_indices()
    slices = sorted({slice_number for slice_number, _roi in regions})

    metadata = shared_metadata()
    shared = metadata["NWBFile"]
    text = dict(
        surgery=shared["surgery"],
        pharmacology=shared["pharmacology"],
        # The shared text ends with the levodopa treatment, which the 6-OHDA mice did not receive.
        pharmacology_before_levodopa=shared["pharmacology"][: shared["pharmacology"].index("L-DOPA priming")],
        injection=session["injection"].format(side=session["injection_side"][category]),
        bath=session["bath"],
    )
    metadata["NWBFile"].update(
        session_description=session["NWBFile"]["session_description"].format(
            n_regions=len(regions), n_slices=len(slices), condition=condition
        ),
        experiment_description=session["NWBFile"]["experiment_description"],
        identifier=str(uuid.uuid4()),
        session_start_time=session_start_time_naive.replace(tzinfo=TIMEZONE),
        session_id=_build_session_id(bundle),
        surgery=session["NWBFile"]["surgery"][category].format(**text),
        pharmacology=session["NWBFile"]["pharmacology"][category].format(**text),
    )
    metadata["Subject"].update(
        subject_id=subject_entry.subject_id,
        strain=subject_entry.strain,
        genotype=subject_entry.genotype,
        description=session["Subject"]["description"].format(group=session["group_by_category"][category]),
    )
    # The subjects table gives no strain or genotype for some lines, which the file then leaves out
    metadata["Subject"] = {key: value for key, value in metadata["Subject"].items() if value is not None}
    nwbfile = make_nwbfile_from_metadata(metadata)
    microscope = Device(**device_metadata("BrukerUltima"), model=get_or_create_bruker_ultima_model(nwbfile))
    nwbfile.add_device(microscope)
    nwbfile.add_device(Device(**session["StimulatingElectrode"]))
    ophys = nwbfile.create_processing_module(name="ophys", description=ophys_text["processing_module_description"])
    segmentation = ImageSegmentation(name="ImageSegmentation")
    ophys.add(segmentation)
    fluorescence = Fluorescence(name="Fluorescence")
    ophys.add(fluorescence)
    references = Images(**ophys_text["ReferenceImages"])
    nwbfile.add_acquisition(references)

    plane_text, segmentation_text = ophys_text["ImagingPlane"], ophys_text["PlaneSegmentation"]
    movie_text, trace_text = ophys_text["TwoPhotonSeries"], ophys_text["RoiResponseSeries"]
    trial_rows = []  # keyword arguments of add_trial, per trial
    trial_references = []  # (index in trial_rows, number in recording order at the region, trials there, images)
    for slice_number, roi_number in regions:
        region_name = f"Slice{slice_number}Roi{roi_number}"
        trials = sorted(
            (
                _read_grab_ach_trial(sw, session_start_time_naive)
                for sw in bundle.sweeps_for((slice_number, roi_number))
            ),
            key=lambda trial: trial.start_time,
        )
        microns_per_pixel = {trial.microns_per_pixel for trial in trials}
        if len(microns_per_pixel) != 1:
            raise ValueError(f"The pixel size changes between the trials of slice {slice_number} ROI {roi_number}")
        imaging_plane = nwbfile.create_imaging_plane(
            name=plane_text["name"].format(region=region_name),
            description=plane_text["description"].format(slice=slice_number, roi=roi_number),
            optical_channel=[
                OpticalChannel(name=name, description=_GRAB_ACH_CHANNEL_DESCRIPTION[name], emission_lambda=float("nan"))
                for name in _GRAB_ACH_CHANNELS.values()
            ],
            device=microscope,
            excitation_lambda=plane_text["excitation_lambda"],
            indicator=plane_text["indicator"],
            location=plane_text["location"],
            grid_spacing=[microns_per_pixel.pop()] * 2,
            grid_spacing_unit="micrometers",
        )

        zooms = sorted({trial.zoom for trial in trials})
        for zoom in zooms:
            # Trials keep their number in recording order at the region across zooms.
            numbered = [(number, trial) for number, trial in enumerate(trials, start=1) if trial.zoom == zoom]
            group = [trial for _number, trial in numbered]
            if len({trial.frames["GRABACh3"].shape[1:] for trial in group}) != 1:
                raise ValueError(f"The frame size changes at zoom {zoom:g} in slice {slice_number} ROI {roi_number}")
            name = region_name + (f"Zoom{zoom:g}" if len(zooms) > 1 else "")
            frame_shape = group[0].frames["GRABACh3"].shape[1:]
            layout = ophys_text["layout"].format(
                n_trials=len(group), width=frame_shape[0], height=frame_shape[1], rate=1 / group[0].frame_period
            )
            if len(zooms) > 1:
                layout = ophys_text["layout_with_zooms"].format(layout=layout, n_zooms=len(zooms), zoom=zoom)
            timestamps = np.concatenate([trial.frame_times for trial in group])
            movies, traces = {}, {}
            plane_segmentation = segmentation.create_plane_segmentation(
                name=segmentation_text["name"].format(name=name),
                description=segmentation_text["description"].format(slice=slice_number, roi=roi_number),
                imaging_plane=imaging_plane,
            )
            plane_segmentation.add_roi(image_mask=np.ones(frame_shape, dtype=bool))
            for channel in _GRAB_ACH_CHANNELS.values():
                about = dict(slice=slice_number, roi=roi_number, what=_GRAB_ACH_CHANNEL_DESCRIPTION[channel])
                movies[channel] = TwoPhotonSeries(
                    name=movie_text["name"].format(channel=channel, name=name),
                    description=movie_text["description"].format(**about, layout=layout, column=channel.lower()),
                    data=np.concatenate([trial.frames[channel] for trial in group]),
                    imaging_plane=imaging_plane,
                    unit="a.u.",
                    timestamps=timestamps if channel == "GRABACh3" else movies["GRABACh3"],
                )
                nwbfile.add_acquisition(movies[channel])
                traces[channel] = RoiResponseSeries(
                    name=trace_text["name"].format(channel=channel, name=name),
                    description=trace_text["description"].format(**about, channel=channel, name=name),
                    data=np.concatenate([trial.brightness[channel] for trial in group]),
                    rois=plane_segmentation.create_roi_table_region(
                        description=segmentation_text["roi_region_description"], region=[0]
                    ),
                    unit="a.u.",
                    timestamps=movies["GRABACh3"],
                )
                fluorescence.add_roi_response_series(traces[channel])

            first_frame = 0
            for number, trial in numbered:
                count = len(trial.frame_times)
                trial_rows.append(
                    dict(
                        start_time=float(trial.frame_times[0]),
                        stop_time=float(trial.frame_times[-1]) + trial.frame_period,
                        region=region_name,
                        drug=trial.drug,
                        stimulation=trial.stimulation,
                        stimulus_start_time=trial.stimulus_start_time,
                        **{
                            f"{kind}_{channel.lower()}": TimeSeriesReference(first_frame, count, series[channel])
                            for kind, series in (("two_photon", movies), ("roi_response", traces))
                            for channel in _GRAB_ACH_CHANNELS.values()
                        },
                    )
                )
                trial_references.append((len(trial_rows) - 1, number, len(trials), trial.references))
                first_frame += count

    for name, description in session["trials_columns"].items():
        nwbfile.add_trial_column(name=name, description=description)
    for kind in ("two_photon", "roi_response"):
        for channel in _GRAB_ACH_CHANNELS.values():
            series_type = "TwoPhotonSeries" if kind == "two_photon" else "RoiResponseSeries"
            nwbfile.add_trial_column(
                name=f"{kind}_{channel.lower()}",
                description=session["trials_series_column"].format(series_type=series_type, channel=channel),
                col_cls=TimeSeriesReferenceVectorData,
            )

    # The trials table lists the trials of all regions in time order.
    order = sorted(range(len(trial_rows)), key=lambda index: trial_rows[index]["start_time"])
    for index in order:
        nwbfile.add_trial(**trial_rows[index])
    row_of = {index: row for row, index in enumerate(order)}
    image_text = ophys_text["ReferenceImage"]
    for index, number, trials_in_region, images in trial_references:
        row = trial_rows[index]
        for channel, image in images.items():
            references.add_image(
                Image(
                    name=image_text["name"].format(channel=channel, region=row["region"], number=number),
                    description=image_text["description"].format(
                        what=_GRAB_ACH_CHANNEL_DESCRIPTION[channel],
                        region=row["region"],
                        number=number,
                        n_trials=trials_in_region,
                        row=row_of[index],
                        start_time=row["start_time"],
                    ),
                    data=image,
                )
            )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    configure_and_write_nwbfile(nwbfile=nwbfile, nwbfile_path=output_path, backend="hdf5", metadata=shared_metadata())
    return output_path
