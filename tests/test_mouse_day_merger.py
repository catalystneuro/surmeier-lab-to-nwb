"""Tests for the per-mouse-day patch-clamp merger. They read the lab's raw data through link_to_raw_data."""

import sys
from pathlib import Path

import numpy as np
import pytest
import tifffile
from pynwb import NWBHDF5IO

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "scripts"))

from prototype_per_mouse_day_merger import (  # noqa: E402
    merge_grab_ach_mouse_day,
    merge_oepsc_mouse_day,
    merge_one_mouse_day,
)

from surmeier_lab_to_nwb.zhai2025.mouse_day_merger import (  # noqa: E402
    _walk_dendritic_mmdd_letter,
    _walk_grab_ach,
    _walk_oepsc_mmddyyyy_letter,
    build_figure_bundles,
)

pytestmark = pytest.mark.skipif(not (REPO_ROOT / "link_to_raw_data").exists(), reason="needs the lab's raw data")


def test_line_scans_kept_when_the_scan_line_changes_length(tmp_path):
    # Cell 1 Proximal 2 of this mouse-day has a 21-pixel scan line in its first trial and 19 pixels in the next two.
    bundle = next(
        b for b in build_figure_bundles("F1") if b.subject_entry.subject_id == "Subject20161111LIDOFFtDTomato"
    )
    nwbfile_path = merge_one_mouse_day(bundle, tmp_path / "mouse_day.nwb")

    with NWBHDF5IO(nwbfile_path, "r") as io:
        nwbfile = io.read()
        dendrite_type = nwbfile.intracellular_recordings["dendrite_type"].data[:]
        index = nwbfile.processing["recordings_index"]["RecordingsIndexTable"]
        for column in ("line_scan_alexa568", "line_scan_fluo4"):
            references = index[column].data[:]
            missing = [row for row, kind in enumerate(dendrite_type) if kind != "Soma" and references[row][0] < 0]
            assert missing == [], f"dendritic rows without a {column} reference: {missing}"
        assert nwbfile.acquisition["LineScanAlexa568Cell1Proximal2Line1"].data.shape[1] == 21
        assert nwbfile.acquisition["LineScanAlexa568Cell1Proximal2Line2"].data.shape[1] == 19


def test_each_trial_roi_is_the_scan_line_of_its_kymograph(tmp_path):
    # Cells 1 and 3 of this mouse-day each have a location whose scan line moves between trials.
    bundle = next(
        b for b in build_figure_bundles("F1") if b.subject_entry.subject_id == "Subject20161111LIDOFFtDTomato"
    )
    nwbfile_path = merge_one_mouse_day(bundle, tmp_path / "mouse_day.nwb")

    with NWBHDF5IO(nwbfile_path, "r") as io:
        nwbfile = io.read()
        dendrite_type = nwbfile.intracellular_recordings["dendrite_type"].data[:]
        index = nwbfile.processing["recordings_index"]["RecordingsIndexTable"]
        for row, kind in enumerate(dendrite_type):
            if kind == "Soma":
                continue
            for channel in ("alexa568", "fluo4"):
                kymograph = index[f"line_scan_{channel}"][row].timeseries
                rois = index[f"roi_response_{channel}"][row].timeseries.rois
                n_pixels = len(rois.table["pixel_mask"][rois.data[0]])
                assert n_pixels == kymograph.data.shape[1], f"row {row} {channel}: ROI has {n_pixels} pixels"
        # Cell 3 Proximal 2 keeps its 21-pixel length but the line moves, so it gets two ROIs.
        fluorescence = nwbfile.processing["ophys"]["Fluorescence"]
        pixel_masks = []
        for line in ("Line1", "Line2"):
            rois = fluorescence[f"RoiResponseSeriesAlexa568Cell3Proximal2{line}"].rois
            pixel_masks.append([tuple(pixel) for pixel in rois.table["pixel_mask"][rois.data[0]]])
        assert pixel_masks[0] != pixel_masks[1]


def test_real_trials_kept_apart_from_the_first_patch(tmp_path):
    # At Distal 1 of cell 2 the lab recorded three trials, moved the microscope to another patch of the same
    # dendrite and recorded three more in folders marked "dist1real".
    bundle = next(b for b in build_figure_bundles("F3") if b.subject_entry.subject_id == "Subject20160524LIDOFFeGFP")
    nwbfile_path = merge_one_mouse_day(bundle, tmp_path / "mouse_day.nwb")

    with NWBHDF5IO(nwbfile_path, "r") as io:
        nwbfile = io.read()
        assert nwbfile.acquisition["LineScanAlexa568Cell2Distal1Spot1"].data.shape == (7500, 29)
        assert nwbfile.acquisition["LineScanAlexa568Cell2Distal1Spot2"].data.shape == (7500, 24)
        cell_id = nwbfile.intracellular_recordings["cell_id"].data[:]
        site = nwbfile.intracellular_recordings["site"].data[:]
        cell_2_sites = [name for cell, name in zip(cell_id, site) if cell == 2]
        assert cell_2_sites.count("Distal1Spot1") == 3
        assert cell_2_sites.count("Distal1Spot2") == 3


def test_sulpiride_folders_are_one_condition(tmp_path):
    # The lab named the sulpiride folders "LID on-state with sul" (dendritic) and "... with sul (iSPN)" (somatic).
    bundle = next(b for b in build_figure_bundles("F3") if b.subject_entry.subject_id == "Subject20190509LIDONeGFP")
    nwbfile_path = merge_one_mouse_day(bundle, tmp_path / "mouse_day.nwb")

    with NWBHDF5IO(nwbfile_path, "r") as io:
        nwbfile = io.read()
        assert list(nwbfile.icephys_experimental_conditions["condition"].data[:]) == ["LID on-state with sulpiride"]
        assert set(nwbfile.intracellular_recordings["condition"].data[:]) == {"LID on-state with sulpiride"}


def test_one_electrode_per_cell_at_the_cell_body(tmp_path):
    # Each cell had one whole-cell pipette at the cell body; the dendritic sites were only imaged.
    bundle = next(
        b for b in build_figure_bundles("F1") if b.subject_entry.subject_id == "Subject20161111LIDOFFtDTomato"
    )
    nwbfile_path = merge_one_mouse_day(bundle, tmp_path / "mouse_day.nwb")

    with NWBHDF5IO(nwbfile_path, "r") as io:
        nwbfile = io.read()
        cell_ids = nwbfile.intracellular_recordings["cell_id"].data[:]
        assert set(nwbfile.icephys_electrodes) == {f"ElectrodeCell{cell}" for cell in set(cell_ids)}
        electrodes = nwbfile.intracellular_recordings.category_tables["electrodes"]["electrode"].data[:]
        assert all(electrode.name == f"ElectrodeCell{cell}" for electrode, cell in zip(electrodes, cell_ids))


def test_misspelled_proxt2_folders_are_converted(tmp_path):
    # The lab named the Proximal 2 folders of cell 3 "05072019_Cell3_proxt2_trio_sul-NNN".
    bundle = next(b for b in build_figure_bundles("F3") if b.subject_entry.subject_id == "Subject20190508LIDONeGFP")
    nwbfile_path = merge_one_mouse_day(bundle, tmp_path / "mouse_day.nwb")

    with NWBHDF5IO(nwbfile_path, "r") as io:
        nwbfile = io.read()
        cell_id = nwbfile.intracellular_recordings["cell_id"].data[:]
        site = nwbfile.intracellular_recordings["site"].data[:]
        assert [name for cell, name in zip(cell_id, site) if cell == 3].count("Proximal2") == 3


def test_trio_folders_named_without_trio_are_converted(tmp_path):
    # The lab named the Proximal 1 folders of cell 3 "20160608_Cell3_prox1-NNN", without "_trio".
    bundle = next(b for b in build_figure_bundles("F3") if b.subject_entry.subject_id == "Subject20160608LIDONeGFP")
    nwbfile_path = merge_one_mouse_day(bundle, tmp_path / "mouse_day.nwb")

    with NWBHDF5IO(nwbfile_path, "r") as io:
        nwbfile = io.read()
        cell_id = nwbfile.intracellular_recordings["cell_id"].data[:]
        site = nwbfile.intracellular_recordings["site"].data[:]
        assert [name for cell, name in zip(cell_id, site) if cell == 3].count("Proximal1") == 3


def test_background_scans_are_kept_out_of_the_trials(tmp_path):
    # Cell 2 has one "_bsl" scan per site: the line moved off the dendrite, no stimulus.
    bundle = next(b for b in build_figure_bundles("F3") if b.subject_entry.subject_id == "Subject20170525LIDOFFeGFP")
    nwbfile_path = merge_one_mouse_day(bundle, tmp_path / "mouse_day.nwb")

    with NWBHDF5IO(nwbfile_path, "r") as io:
        nwbfile = io.read()
        backgrounds = sorted(name for name in nwbfile.acquisition if name.startswith("LineScanBackground"))
        sites = ["Distal1", "Distal2", "Proximal1", "Proximal2"]
        assert backgrounds == sorted(
            f"LineScanBackground{dye}Cell2{site}" for dye in ("Alexa568", "Fluo4") for site in sites
        )
        cell_id = nwbfile.intracellular_recordings["cell_id"].data[:]
        site = nwbfile.intracellular_recordings["site"].data[:]
        assert [name for cell, name in zip(cell_id, site) if cell == 2].count("Proximal1") == 3


def test_unrecognized_recording_folder_raises(tmp_path):
    # A folder whose name the converter cannot parse used to be skipped without a word.
    (tmp_path / "dendritic" / "LID off-state" / "0524b" / "05242016_Cell2_prox1_notes").mkdir(parents=True)
    with pytest.raises(ValueError, match="Unrecognized dendritic recording folder"):
        _walk_dendritic_mmdd_letter(tmp_path, "dendritic", ["LID off-state"])


def test_recording_without_its_time_raises(tmp_path):
    # A recording folder without the Prairie View file that holds its start time used to be skipped.
    (tmp_path / "dendritic" / "LID off-state" / "0524b" / "05242016_Cell2_prox1_trio-001").mkdir(parents=True)
    with pytest.raises(ValueError, match="No recording time"):
        _walk_dendritic_mmdd_letter(tmp_path, "dendritic", ["LID off-state"])


def test_sr_oepsc_cells_grouped_by_mouse_day():
    # On 2023-08-16 the lab named cell 1's folders "LED20-NNN", without the cell; the folder letter "a" gives it.
    bundle = next(b for b in build_figure_bundles("F2") if b.subject_entry.subject_id == "Subject20230816LIDONtDTomato")
    assert bundle.cell_indices() == [1, 2]
    assert [s.recording_path.parent.name for s in bundle.sweeps_for(1)] == ["08162023a"] * 10
    assert {s.led_setting for s in bundle.sweeps_for(1)} == {20.0}
    assert len(bundle.sweeps_for(2)) == 10


def test_sr_oepsc_recording_naming_another_cell_raises(tmp_path):
    (tmp_path / "Sr-oEPSC" / "LID off-state" / "07062023b" / "cell1_LED20-001").mkdir(parents=True)
    with pytest.raises(ValueError, match="names cell 1 but is in the folder of cell 2"):
        _walk_oepsc_mmddyyyy_letter(tmp_path, "Sr-oEPSC", ["LID off-state"])


def test_sr_oepsc_file_holds_every_cell_of_the_mouse_day(tmp_path):
    # The published files kept one cell per day: for 2023-07-06 only cell 3, whose first sample is -14.83 pA.
    bundle = next(
        b for b in build_figure_bundles("F2") if b.subject_entry.subject_id == "Subject20230706LIDOFFtDTomato"
    )
    nwbfile_path = merge_oepsc_mouse_day(bundle, tmp_path / "mouse_day.nwb")

    with NWBHDF5IO(nwbfile_path, "r") as io:
        nwbfile = io.read()
        irt = nwbfile.intracellular_recordings
        assert sorted(set(irt["cell_id"].data[:])) == [1, 2, 3]
        assert nwbfile.acquisition["VoltageClampSeriesCell3"].data[0] == pytest.approx(-1.48315430e-11)
        # Every sweep links its amplifier command and has one 0.3 ms flash 20 ms after it starts.
        stimuli = irt.category_tables["stimuli"]["stimulus"].data[:]
        assert all(start_index >= 0 for start_index, _count, _series in stimuli)
        pulses = nwbfile.intervals["optogenetic_pulses_table"]
        np.testing.assert_allclose(pulses["start_time"][:] - irt["sweep_start_time"].data[:], 0.020)
        np.testing.assert_allclose(pulses["stop_time"][:] - pulses["start_time"][:], 0.0003)


def test_grab_ach_regions_grouped_by_mouse_day():
    # The published files kept one region per day; the paper counts 13, 10 and 11 regions from 3, 3 and 5 mice.
    bundles = build_figure_bundles("F5")
    assert len(bundles) == 11
    assert sum(len(bundle.cells) for bundle in bundles) == 34
    bundle = next(b for b in bundles if b.subject_entry.subject_id == "Subject20240402CONTC57bl6J")
    assert bundle.cell_indices() == [(1, 1), (1, 2), (2, 1), (2, 2), (3, 1), (3, 2)]
    # The false start filed under slice 1 ROI 1 is left out.
    names = [s.recording_path.name for s in bundle.sweeps_for((1, 1))]
    assert len(names) == 20 and "BOT_04022024_slice1ROI_sul_sinlge-001" not in names
    # On 2024-05-31 the lab named two recordings of slice 1 "slice1realA".
    bundle = next(b for b in bundles if b.subject_entry.subject_id == "Subject20240531LIDOFFC57bl6J")
    assert sum("slice1realA" in s.recording_path.name for s in bundle.sweeps_for((1, 1))) == 2


def test_grab_ach_recording_naming_another_region_raises(tmp_path):
    (tmp_path / "UL control" / "04022024slice1ROI1" / "BOT_04022024_slice1ROI2_ctr_single-001").mkdir(parents=True)
    with pytest.raises(ValueError, match="names slice1ROI2 but is in the folder of slice 1 ROI 1"):
        _walk_grab_ach(tmp_path)


def test_grab_ach_file_holds_every_region_of_the_mouse_day(tmp_path):
    # The published file for 2024-05-31 held only slice 2, which was also filmed at zoom 8 in two trials.
    bundle = next(b for b in build_figure_bundles("F5") if b.subject_entry.subject_id == "Subject20240531LIDOFFC57bl6J")
    nwbfile_path = merge_grab_ach_mouse_day(bundle, tmp_path / "mouse_day.nwb")

    with NWBHDF5IO(nwbfile_path, "r") as io:
        nwbfile = io.read()
        trials = nwbfile.trials.to_dataframe()
        assert trials["region"].value_counts().to_dict() == {"Slice2Roi1": 19, "Slice1Roi1": 18}
        assert nwbfile.acquisition["TwoPhotonSeriesGRABACh3Slice2Roi1Zoom8"].data.shape == (307, 64, 64)
        first = bundle.sweeps_for((1, 1))[0].recording_path
        frame = tifffile.imread(next(first.glob("*_Ch2_000001.ome.tif")), is_ome=False)
        np.testing.assert_array_equal(nwbfile.acquisition["TwoPhotonSeriesGRABACh3Slice1Roi1"].data[0], frame.T)
        # Each trial's frames start at its start time, and its stimulus falls inside it.
        for _, trial in trials.iterrows():
            frames = trial["two_photon_grabach3"]
            assert frames.timeseries.timestamps[frames.idx_start] == trial["start_time"]
        stimulated = trials.dropna(subset=["stimulus_start_time"])
        assert (stimulated["stimulus_start_time"] > stimulated["start_time"]).all()
        assert (stimulated["stimulus_start_time"] < stimulated["stop_time"]).all()
