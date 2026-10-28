"""Tests for the per-mouse-day patch-clamp merger. They read the lab's raw data through link_to_raw_data."""

import sys
from pathlib import Path

import pytest
from pynwb import NWBHDF5IO

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "scripts"))

from prototype_per_mouse_day_merger import merge_one_mouse_day  # noqa: E402

from surmeier_lab_to_nwb.zhai2025.mouse_day_merger import (  # noqa: E402
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
        assert nwbfile.acquisition["LineScanCell1_Proximal2_Alexa568_Line1"].data.shape[1] == 21
        assert nwbfile.acquisition["LineScanCell1_Proximal2_Alexa568_Line2"].data.shape[1] == 19


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
        assert nwbfile.acquisition["LineScanCell2_Distal1_Alexa568"].data.shape == (7500, 29)
        assert nwbfile.acquisition["LineScanCell2_Distal1Real_Alexa568"].data.shape == (7500, 24)
        cell_id = nwbfile.intracellular_recordings["cell_id"].data[:]
        site_label = nwbfile.intracellular_recordings["site_label"].data[:]
        cell_2_sites = [site for cell, site in zip(cell_id, site_label) if cell == 2]
        assert cell_2_sites.count("dist1") == 3
        assert cell_2_sites.count("dist1real") == 3
