"""Tests for the figure reproduction scripts."""

import sys
from pathlib import Path

import pandas as pd
import pytest
from scipy import stats

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "figure_scripts"))

from figure_2gh_oepsc_reproduction import generate_markdown_report  # noqa: E402
from figure_5f_acetylcholine_sensor_reproduction import (
    analyze_nwb_dataset,  # noqa: E402
)


def test_figure_2_compares_the_states_with_the_mann_whitney_test(tmp_path):
    # The comparison used Welch's t-test; the paper's Figure 2H uses the Mann-Whitney test.
    off_state, on_state = [8.0, 9.5, 10.1, 11.0], [12.0, 14.5, 15.2, 30.0]
    df = pd.DataFrame(
        {"condition": ["LID off-state"] * 4 + ["LID on-state"] * 4, "oepsc_amplitude_pA": off_state + on_state}
    )
    report = Path(generate_markdown_report(df, tmp_path)).read_text()
    p_value = stats.mannwhitneyu(off_state, on_state, alternative="two-sided").pvalue
    assert "Mann-Whitney U test" in report
    assert f"**p-value:** {p_value:.4f}" in report


def test_figure_5f_raises_on_a_file_it_cannot_read(tmp_path):
    # Every error was caught and printed, and the file left out of the figure.
    (tmp_path / "broken.nwb").write_bytes(b"not an HDF5 file")
    with pytest.raises(OSError):
        analyze_nwb_dataset(tmp_path)
