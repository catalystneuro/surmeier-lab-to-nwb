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
        for column in ("line_scan_ch1", "line_scan_ch2"):
            references = index[column].data[:]
            missing = [row for row, kind in enumerate(dendrite_type) if kind != "Soma" and references[row][0] < 0]
            assert missing == [], f"dendritic rows without a {column} reference: {missing}"
        assert nwbfile.acquisition["LineScanCell1_Proximal2_Ch1_Line1"].data.shape[1] == 21
        assert nwbfile.acquisition["LineScanCell1_Proximal2_Ch1_Line2"].data.shape[1] == 19
