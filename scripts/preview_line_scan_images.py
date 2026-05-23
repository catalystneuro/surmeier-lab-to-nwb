"""Preview all Alexa 568 source images in a paired-recording dandiset file.

Streams one merged file from DANDI and saves a grid of all its Alexa 568 source images
(one per cell × location × trial) so we can visually pick the one with the clearest
dendrite for the notebook's line-scan visualization.

Output: nwb_files/preview_source_images.png
"""

from __future__ import annotations

import os
import re
from pathlib import Path

import h5py
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import remfile
from dandi.dandiapi import DandiAPIClient
from dotenv import load_dotenv
from pynwb import NWBHDF5IO

REPO_ROOT = Path(__file__).resolve().parents[1]
load_dotenv(REPO_ROOT / ".env")


def main():
    client = DandiAPIClient(token=os.environ["DANDI_API_TOKEN"])
    client.authenticate(token=os.environ["DANDI_API_TOKEN"])
    ds = client.get_dandiset("001538", "draft")

    # Pick the same paired-recording file the notebook uses
    target = next(
        a for a in ds.get_assets() if "dSPN++OnState++D1RaSch" in a.path and a.path.endswith("_icephys+ophys.nwb")
    )
    print(f"Streaming: {target.path}")

    s3_url = target.get_content_url(follow_redirects=1, strip_query=False)
    rf = remfile.File(s3_url)
    h5 = h5py.File(rf, mode="r")
    io = NWBHDF5IO(file=h5, load_namespaces=True)
    nwb = io.read()

    container = nwb.acquisition["ImageLineScanSource"]
    # Filter to Ch1 (Alexa 568) — the structural channel
    alexa_names = sorted(n for n in container.images if "Alexa568" in n)
    print(f"Found {len(alexa_names)} Alexa 568 source images")

    # Plot a grid: rows = (cell, location), cols = trials
    pattern = re.compile(r"ImageAlexa568Cell(\d+)(Proximal|Distal)(\d+)Trial(\d+)")
    groups: dict[tuple[int, str, int], list[tuple[int, str]]] = {}
    for name in alexa_names:
        m = pattern.match(name)
        if not m:
            continue
        cell, loc, loc_num, trial = int(m.group(1)), m.group(2), int(m.group(3)), int(m.group(4))
        groups.setdefault((cell, loc, loc_num), []).append((trial, name))

    n_rows = len(groups)
    n_cols = max(len(v) for v in groups.values())
    cmap = mcolors.LinearSegmentedColormap.from_list("alexa568", ["black", "#FF4500"], N=256)

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(4 * n_cols, 4 * n_rows))
    if n_rows == 1:
        axes = np.array([axes])
    if n_cols == 1:
        axes = axes.reshape(-1, 1)

    for row_idx, (key, entries) in enumerate(sorted(groups.items())):
        cell, loc, loc_num = key
        entries.sort()
        for col_idx, (trial, name) in enumerate(entries):
            ax = axes[row_idx, col_idx]
            img = np.asarray(container.images[name].data, dtype=float)
            ax.imshow(img, cmap=cmap, aspect="equal", origin="upper")
            ax.set_title(f"Cell {cell} {loc}{loc_num} T{trial}", fontsize=10)
            ax.axis("off")
        # Hide unused trial columns in this row
        for col_idx in range(len(entries), n_cols):
            axes[row_idx, col_idx].axis("off")

    fig.suptitle(f"Alexa 568 source images: {target.path.split('/')[-1]}", fontsize=12, y=0.998)
    out_path = REPO_ROOT / "nwb_files" / "preview_source_images.png"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout(rect=(0, 0, 1, 0.98))
    fig.savefig(out_path, dpi=110, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out_path}")


if __name__ == "__main__":
    main()
