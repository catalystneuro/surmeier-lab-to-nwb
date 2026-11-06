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
    _build_session_id,
    _general_metadata,
    _parse_voltage_output,
    merge_grab_ach_mouse_day,
    merge_oepsc_mouse_day,
    merge_one_mouse_day,
)

from surmeier_lab_to_nwb.zhai2025.mouse_day_merger import (  # noqa: E402
    _walk_dendritic_mmdd_letter,
    _walk_grab_ach,
    _walk_oepsc_mmddyyyy_letter,
    _walk_oxo_m,
    build_figure_bundles,
)

pytestmark = pytest.mark.skipif(not (REPO_ROOT / "link_to_raw_data").exists(), reason="needs the lab's raw data")


def test_line_scans_kept_when_the_scan_line_changes_length(tmp_path):
    # Cell 2 Proximal 2 of this mouse-day has a 21-pixel scan line in its first trial and 19 pixels in the next two.
    bundle = next(
        b for b in build_figure_bundles("F1") if b.subject_entry.subject_id == "Subject20161111LIDOFFtDTomato"
    )
    nwbfile_path = merge_one_mouse_day(bundle, tmp_path / "mouse_day.nwb")

    with NWBHDF5IO(nwbfile_path, "r") as io:
        nwbfile = io.read()
        dendrite_type = nwbfile.intracellular_recordings["dendrite_type"].data[:]
        index = nwbfile.processing["icephys"]["RecordingsIndexTable"]
        for column in ("line_scan_alexa568", "line_scan_fluo4"):
            references = index[column].data[:]
            missing = [row for row, kind in enumerate(dendrite_type) if kind != "Soma" and references[row][0] < 0]
            assert missing == [], f"dendritic rows without a {column} reference: {missing}"
        assert nwbfile.acquisition["LineScanAlexa568Cell2Proximal2Line1"].data.shape[1] == 21
        assert nwbfile.acquisition["LineScanAlexa568Cell2Proximal2Line2"].data.shape[1] == 19


def test_each_scan_line_gets_its_own_roi(tmp_path):
    # Cell 3 Proximal 2 of this mouse-day keeps its 21-pixel length but the line moves, so it gets two ROIs.
    bundle = next(
        b for b in build_figure_bundles("F1") if b.subject_entry.subject_id == "Subject20161111LIDOFFtDTomato"
    )
    nwbfile_path = merge_one_mouse_day(bundle, tmp_path / "mouse_day.nwb")

    with NWBHDF5IO(nwbfile_path, "r") as io:
        nwbfile = io.read()
        fluorescence = nwbfile.processing["ophys"]["Fluorescence"]
        image_masks = []
        for line in ("Line1", "Line2"):
            rois = fluorescence[f"RoiResponseSeriesAlexa568Cell3Proximal2{line}"].rois
            image_masks.append(np.asarray(rois.table["image_mask"][rois.data[0]]))
        assert not np.array_equal(*image_masks)


def test_fluorescence_roi_is_the_window_prairie_view_averages(tmp_path):
    # Prairie View's line profile is the mean of the line scan over the pixels of its profile window, 0 to 16 of
    # this 19-pixel line, and not over the whole line, which the files used to give as the ROI.
    bundle = next(
        b for b in build_figure_bundles("F1") if b.subject_entry.subject_id == "Subject20161111LIDOFFtDTomato"
    )
    nwbfile_path = merge_one_mouse_day(bundle, tmp_path / "mouse_day.nwb")

    with NWBHDF5IO(nwbfile_path, "r") as io:
        nwbfile = io.read()
        roi_response = nwbfile.processing["ophys"]["Fluorescence"]["RoiResponseSeriesFluo4Cell2Proximal2Line2"]
        rois = roi_response.rois
        assert np.asarray(rois.table["image_mask"][rois.data[0]]).sum() == 17
        line_scan = nwbfile.acquisition["LineScanFluo4Cell2Proximal2Line2"].data[:]
        np.testing.assert_allclose(roi_response.data[:, 0], line_scan[:, :17].mean(axis=1), atol=1)


def test_line_scans_start_with_their_sweeps(tmp_path):
    # The trials of a site used to be stored back to back from the first one's start, so later trials sat too early.
    bundle = next(
        b for b in build_figure_bundles("F1") if b.subject_entry.subject_id == "Subject20161111LIDOFFtDTomato"
    )
    nwbfile_path = merge_one_mouse_day(bundle, tmp_path / "mouse_day.nwb")

    with NWBHDF5IO(nwbfile_path, "r") as io:
        nwbfile = io.read()
        dendrite_type = nwbfile.intracellular_recordings["dendrite_type"].data[:]
        sweep_start_time = nwbfile.intracellular_recordings["sweep_start_time"].data[:]
        index = nwbfile.processing["icephys"]["RecordingsIndexTable"]
        for column in ("line_scan_alexa568", "line_scan_fluo4", "roi_response_alexa568", "roi_response_fluo4"):
            for row, kind in enumerate(dendrite_type):
                if kind != "Soma":
                    assert index[column][row].timestamps[0] == sweep_start_time[row], f"{column} row {row}"


@pytest.mark.parametrize(
    "figure, subject_id, merge",
    [
        ("F3", "Subject20160524LIDOFFeGFP", merge_one_mouse_day),
        ("F2", "Subject20230705LIDONtDTomato", merge_oepsc_mouse_day),
    ],
)
def test_sweeps_start_at_their_recorded_time(tmp_path, figure, subject_id, merge):
    # The sweeps of a cell used to be stored back to back from the first one's start, so later sweeps sat too early.
    # The patch-clamp day has two sweeps recorded for 1 s, whose stimulus outlasts the recording.
    bundle = next(b for b in build_figure_bundles(figure) if b.subject_entry.subject_id == subject_id)
    nwbfile_path = merge(bundle, tmp_path / "mouse_day.nwb")

    with NWBHDF5IO(nwbfile_path, "r") as io:
        nwbfile = io.read()
        irt = nwbfile.intracellular_recordings
        sweep_start_time = irt["sweep_start_time"].data[:]
        for category, column in (("responses", "response"), ("stimuli", "stimulus")):
            for row, (start_index, _, series) in enumerate(irt.category_tables[category][column].data[:]):
                assert series.timestamps[start_index] == sweep_start_time[row], f"{column} row {row}"


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


def test_stimulus_is_the_commanded_current(tmp_path):
    # The dendritic stimulus used to be all zero: only the command's first part, 0 pA for 100 ms, was read. The three
    # 2 ms, 2 nA pulses follow 900 ms after it; a somatic step starts 200 ms in and lasts 500 ms. Two dendritic trials
    # of cell 3 were recorded for only 1 s, before the pulses; their stimulus is the whole 1.152 s command.
    bundle = next(b for b in build_figure_bundles("F3") if b.subject_entry.subject_id == "Subject20160524LIDOFFeGFP")
    nwbfile_path = merge_one_mouse_day(bundle, tmp_path / "mouse_day.nwb")

    with NWBHDF5IO(nwbfile_path, "r") as io:
        nwbfile = io.read()
        irt = nwbfile.intracellular_recordings
        dendrite_type = irt["dendrite_type"].data[:]
        current_pA = irt["stimulus_current_pA"].data[:]
        stimuli = irt.category_tables["stimuli"]["stimulus"].data[:]
        responses = irt.category_tables["responses"]["response"].data[:]
        assert {kind == "Soma" for kind in dendrite_type} == {True, False}
        short_trials = 0
        for row, ((start_index, count, series), (_, response_count, _)) in enumerate(zip(stimuli, responses)):
            samples_per_ms = round(1 / (series.timestamps[start_index + 1] - series.timestamps[start_index])) / 1000
            expected_pA = np.zeros(count)
            if dendrite_type[row] == "Soma":
                expected_pA[round(200 * samples_per_ms) : round(700 * samples_per_ms)] = current_pA[row]
            else:
                assert current_pA[row] == 2000
                for start_ms in (1000, 1020, 1040):
                    expected_pA[round(start_ms * samples_per_ms) : round((start_ms + 2) * samples_per_ms)] = 2000
            if response_count != count:
                short_trials += 1
                assert (response_count, count) == (round(1000 * samples_per_ms), round(1152 * samples_per_ms))
            np.testing.assert_allclose(series.data[start_index : start_index + count], expected_pA * 1e-12, rtol=1e-6)
        assert short_trials == 2


def test_unrecognized_recording_folder_raises(tmp_path):
    # A folder whose name the converter cannot parse used to be skipped without a word.
    (tmp_path / "dendritic" / "LID off-state" / "0524b" / "05242016_Cell2_prox1_notes").mkdir(parents=True)
    with pytest.raises(ValueError, match="Unrecognized dendritic recording folder"):
        _walk_dendritic_mmdd_letter(tmp_path, "dendritic", ["LID off-state"])


def test_dendritic_cell_number_comes_from_the_folder_letter():
    # Folder 1111b holds cell 2 by its letter and the lab's sheet; its recordings' names say Cell1, which made the
    # somatic and dendritic recordings of that cell two cells.
    bundle = next(b for b in build_figure_bundles("F1") if b.date_compact == "20161111")
    assert bundle.cell_indices() == [2, 3]
    assert {s.modality for s in bundle.sweeps_for(2)} == {"soma_excitability", "dend_excitability"}


def test_dendritic_recording_naming_another_cell_raises(tmp_path):
    (tmp_path / "dendritic" / "LID off-state" / "0524b" / "05242016_Cell1_prox1_trio-001").mkdir(parents=True)
    with pytest.raises(ValueError, match="names cell 1 in the folder of cell 2"):
        _walk_dendritic_mmdd_letter(tmp_path, "dendritic", ["LID off-state"])


def test_recordings_of_one_patch_numbered_apart_are_one_cell(tmp_path):
    # The lab numbered these somatic and dendritic recordings differently; they follow each other with nothing in
    # between, sooner than a new cell can be patched and filled with dye.
    bundles = {b.date_compact: b for b in build_figure_bundles("F3")}
    assert bundles["20160524"].cell_indices() == [1, 2, 4]
    assert bundles["20190426"].cell_indices() == [2, 3]
    for bundle in (bundles["20160524"], bundles["20190426"]):
        for cell in bundle.cell_indices():
            assert {s.modality for s in bundle.sweeps_for(cell)} == {"soma_excitability", "dend_excitability"}

    nwbfile_path = merge_one_mouse_day(bundles["20160524"], tmp_path / "mouse_day.nwb")
    with NWBHDF5IO(nwbfile_path, "r") as io:
        nwbfile = io.read()
        dendritic = [
            series for series in nwbfile.acquisition.values() if "dendritic recordings cell 3" in series.description
        ]
        assert len(dendritic) == 1 and dendritic[0].electrode.name == "ElectrodeCell4"


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


def test_oxo_m_dendrites_grouped_by_mouse_day():
    # The published files had no oxo-M recordings; the paper counts 10 dendrites from 5 wild-type cells and 12 from
    # 6 CDGI knockout cells.
    bundles = [b for b in build_figure_bundles("F7") if b.subject_entry.category == "NAIVE"]
    cells = {b.subject_entry.subject_id: b.cell_indices() for b in bundles}
    assert cells == {
        "Subject20180221NAIVEeGFP": [1, 2],
        "Subject20180228NAIVECDGIkoXeGFP": [2, 3],
        "Subject20180301NAIVECDGIkoXeGFP": [2, 3],
        # The third cell's folders say c3, its recordings Cell4.
        "Subject20180305NAIVEeGFP": [1, 2, 4],
        "Subject20180314NAIVECDGIkoXeGFP": [2, 3],
    }
    dendrites = {"WT": set(), "KO": set()}
    for bundle in bundles:
        for cell in bundle.cell_indices():
            sweeps = bundle.sweeps_for(cell)
            dendrites[sweeps[0].condition_subfolder] |= {(bundle.date_compact, cell, s.location_number) for s in sweeps}
            # Every trial after oxo-M was recorded after every trial before it.
            before = [s.pv_scan_date for s in sweeps if s.suffix is None]
            after = [s.pv_scan_date for s in sweeps if s.suffix == "oxoM"]
            assert len(before) == 12 and len(after) >= 12 and max(before) < min(after)
    assert {genotype: len(keys) for genotype, keys in dendrites.items()} == {"WT": 10, "KO": 12}


def test_oxo_m_recording_naming_another_cell_raises(tmp_path):
    (tmp_path / "oxoM on DE" / "WT" / "0221c1d1" / "02212018_Cell2_prox1_trio-001").mkdir(parents=True)
    with pytest.raises(ValueError, match="does not match the day, cell or dendrite of its folder"):
        _walk_oxo_m(tmp_path, "oxoM on DE", ["WT", "KO"])


def test_oxo_m_file_splits_trials_before_and_after_oxo_m(tmp_path):
    bundle = next(b for b in build_figure_bundles("F7") if b.subject_entry.subject_id == "Subject20180221NAIVEeGFP")
    nwbfile_path = merge_one_mouse_day(bundle, tmp_path / "mouse_day.nwb")

    with NWBHDF5IO(nwbfile_path, "r") as io:
        nwbfile = io.read()
        assert nwbfile.session_id == "iSPN++Naive++M1RaOxoM++WT++20180221"
        assert list(nwbfile.icephys_experimental_conditions["condition"].data[:]) == ["naive", "naive with oxo-M"]
        irt = nwbfile.intracellular_recordings
        conditions = irt["condition"].data[:]
        start_times = irt["sweep_start_time"].data[:]
        for cell in (1, 2):
            rows = irt["cell_id"].data[:] == cell
            before = start_times[rows & (conditions == "naive")]
            after = start_times[rows & (conditions == "naive with oxo-M")]
            assert len(before) == len(after) == 12 and before.max() < after.min()
        # The lesion and levodopa text of the other files does not apply; the sheet's label is mentioned.
        assert nwbfile.surgery.startswith("None: naive mouse")
        assert "'LID ON'" in nwbfile.subject.description


def test_knockout_files_carry_the_knockout_genotype(tmp_path):
    # The patch-clamp files had no genotype, which NWB reads as wild type, so the CDGI knockouts looked wild type.
    bundle = next(
        b for b in build_figure_bundles("F7") if b.subject_entry.subject_id == "Subject20180228NAIVECDGIkoXeGFP"
    )
    nwbfile_path = merge_one_mouse_day(bundle, tmp_path / "mouse_day.nwb")

    with NWBHDF5IO(nwbfile_path, "r") as io:
        subject = io.read().subject
        assert subject.genotype == "CDGI knockout; Tg(Drd2-EGFP)S118Gsat, hemizygous (RRID:MMRRC_000230-UNC)"
        assert subject.strain == "C57BL/6"


def test_control_groups_get_the_state_of_their_mice():
    # The lab's sheet files all three as CONT, which the session IDs used to turn into Control.
    states = {}
    for figure in ("F5", "F6", "F8"):
        for bundle in build_figure_bundles(figure):
            conditions = {s.condition_subfolder for sweeps in bundle.cells.values() for s in sweeps}
            for condition in conditions & {"UL control", "control", "interleaved control"}:
                states.setdefault(condition, set()).add(_build_session_id(bundle).split("++")[1])
    # GRAB-ACh unlesioned mice; Figure 6 vehicle and Figure 8 gRNA-only controls, both off-state mice.
    assert states == {"UL control": {"UnlesionedControl"}, "control": {"OffState"}, "interleaved control": {"OffState"}}


def test_session_description_names_the_published_figure_and_paper(tmp_path):
    # The description was an internal note ("Per-mouse-day icephys merge ... category=M1R CRSPR"), the lab's Figure 8
    # recordings are Figure 9 of the published paper, and the files cited only the preprint.
    bundle = next(
        b for b in build_figure_bundles("F8") if b.subject_entry.subject_id == "Subject20221004M1RCRSPRAdora2Cre"
    )
    nwbfile_path = merge_one_mouse_day(bundle, tmp_path / "mouse_day.nwb")

    with NWBHDF5IO(nwbfile_path, "r") as io:
        nwbfile = io.read()
        assert nwbfile.session_description.startswith(
            "Patch-clamp recordings (Figure 9): whole-cell current clamp of 2 indirect-pathway spiny projection "
            "neurons (iSPNs) from one mouse. Condition: LID off-state, M1R deletion. Somatic excitability of 2 cells"
        )
        assert list(nwbfile.related_publications) == ["https://doi.org/10.1126/sciadv.adv8224"]


def test_subject_description_describes_the_mouse(tmp_path):
    # The description was an internal note ("Reporter line Adora2-Cre, category M1R CRSPR, ... See general_metadata.yaml
    # surgery section"), the same for every group.
    bundle = next(
        b for b in build_figure_bundles("F8") if b.subject_entry.subject_id == "Subject20221004M1RCRSPRAdora2Cre"
    )
    with NWBHDF5IO(merge_one_mouse_day(bundle, tmp_path / "patch_clamp.nwb"), "r") as io:
        description = io.read().subject.description
    assert description.startswith(
        "Adora2a-Cre transgenic mouse, whose indirect-pathway spiny projection neurons (iSPNs) express Cre recombinase. "
        "It received a unilateral 6-hydroxydopamine (6-OHDA) lesion of the medial forebrain bundle and, after the "
        "cylinder test, an injection of adeno-associated viruses (AAVs) carrying Cre-dependent Cas9"
    )

    bundle = next(
        b for b in build_figure_bundles("F2") if b.subject_entry.subject_id == "Subject20230706LIDOFFtDTomato"
    )
    with NWBHDF5IO(merge_oepsc_mouse_day(bundle, tmp_path / "sr_oepsc.nwb"), "r") as io:
        description = io.read().subject.description
    assert "brain slices were made 24 to 48 hours after the last dose, in the off-state" in description
    assert "A channelrhodopsin-2 (ChR2) virus was injected into the motor cortex" in description


@pytest.mark.parametrize(
    "figure, subject_id, field, expected",
    [
        ("F6", "Subject20200722M1RanteGFP", "pharmacology", "VU 0255035 (5 µM)"),
        ("F6", "Subject20200811CONTeGFP", "pharmacology", "saline was injected i.p."),
        ("F8", "Subject20221004M1RCRSPRAdora2Cre", "surgery", "AAV9-hSyn-DIO-saCAS-minWPRE"),
    ],
)
def test_figure_6_and_8_files_name_their_treatments(tmp_path, figure, subject_id, field, expected):
    # The files had the general lesion and levodopa text only: no trihexyphenidyl, saline or VU 0255035 for Figure 6,
    # no CRISPR virus for Figure 8.
    bundle = next(b for b in build_figure_bundles(figure) if b.subject_entry.subject_id == subject_id)
    with NWBHDF5IO(merge_one_mouse_day(bundle, tmp_path / "mouse_day.nwb"), "r") as io:
        assert expected in getattr(io.read(), field)


def test_stimulus_current_has_no_float_noise():
    # 0.55 V at 400 pA/V is 220.00000000000003 in floating point, which 88 somatic sweeps had in stimulus_current_pA.
    xml_path = (
        REPO_ROOT / "link_to_raw_data/Figure 1/Somatic excitability/LID off-state/02022017_1/cell1-017"
        "/cell1-017_Cycle00001_VoltageOutput_001.xml"
    )
    assert repr(_parse_voltage_output(xml_path).amplitude_pA) == "220.0"


def test_lesion_is_in_the_left_medial_forebrain_bundle():
    # The surgery and pharmacology text said the right one. The lesion protocol's ML -1.2 mm is left in the paper's
    # convention (it injects "the left DLS" at ML -2.35), and the protocol keeps mice that rarely use the right paw.
    general = _general_metadata()["NWBFile"]
    for field in ("surgery", "pharmacology"):
        text = " ".join(general[field].split())
        assert "into the left medial forebrain bundle" in text
        assert "right medial forebrain bundle" not in text


def test_grna_only_controls_have_the_wild_type_token():
    # The Figure 8 gRNA-only controls got M1RCRISPR, the token of the mice with the deletion, while their behaviour
    # files say WT; without the Cas9 virus their M1R is intact.
    tokens = {}
    for bundle in build_figure_bundles("F8"):
        tokens.setdefault(bundle.subject_entry.category, set()).add(_build_session_id(bundle).split("++")[3])
    assert tokens == {"CONT": {"WT"}, "M1R CRSPR": {"M1RCRISPR"}}
