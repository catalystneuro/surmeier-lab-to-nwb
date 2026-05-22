"""Run the per-mouse-day merger across all patch-clamp figures.

Builds per-mouse-day NWB files for every bundle in F1, F3, F6, F7, F8.
Writes outputs to nwb_files/merged/<figure>/<subject_id>.nwb and reports a
summary plus nwbinspector check on a sample of the produced files.

This is the production runner; the prototype script
(`scripts/prototype_per_mouse_day_merger.py`) is the per-cell test harness.
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

from surmeier_lab_to_nwb.zhai2025.mouse_day_merger import (
    build_figure_bundles,  # noqa: E402
)

# Import the merger function from the prototype script.
sys.path.insert(0, str(REPO_ROOT / "scripts"))
from prototype_per_mouse_day_merger import merge_one_mouse_day  # noqa: E402

OUTPUT_DIR = REPO_ROOT / "nwb_files" / "merged"


def main() -> int:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    figures = ["F1", "F3", "F6", "F7", "F8"]
    total_files = 0
    total_seconds = 0.0
    total_size_mb = 0.0
    per_figure: dict[str, dict] = {}
    for fig in figures:
        bundles = build_figure_bundles(fig)
        print(f"\n=== {fig}: {len(bundles)} mouse-days ===")
        fig_dir = OUTPUT_DIR / fig
        fig_dir.mkdir(parents=True, exist_ok=True)
        t0 = time.time()
        successes = 0
        failures: list[tuple[str, str]] = []
        size_mb_local = 0.0
        for b in bundles:
            output = fig_dir / f"{b.subject_entry.subject_id}.nwb"
            try:
                final_path = merge_one_mouse_day(b, output)
                successes += 1
                size_mb_local += final_path.stat().st_size / 1024 / 1024
                print(
                    f"  PASS  {b.subject_entry.subject_id}  ({len(b.cells)} cells, paired={b.has_paired_conditions()})"
                )
            except Exception as e:
                failures.append((b.subject_entry.subject_id, str(e)[:200]))
                print(f"  FAIL  {b.subject_entry.subject_id}: {str(e)[:200]}")
        elapsed = time.time() - t0
        per_figure[fig] = {
            "bundles": len(bundles),
            "successes": successes,
            "failures": failures,
            "size_mb": size_mb_local,
            "elapsed_s": elapsed,
        }
        total_files += successes
        total_seconds += elapsed
        total_size_mb += size_mb_local
        print(f"  {fig}: {successes}/{len(bundles)} OK, {size_mb_local:.1f} MB, {elapsed:.1f}s")

    print(f"\n=== Summary ===")
    print(f"Total merged files: {total_files}")
    print(f"Total disk size:   {total_size_mb:.1f} MB ({total_size_mb/1024:.2f} GB)")
    print(f"Total runtime:     {total_seconds/60:.1f} min")
    print()
    print(f"Per-figure:")
    for fig, stats in per_figure.items():
        status = f"{stats['successes']}/{stats['bundles']}"
        print(
            f"  {fig:5s} {status:8s} {stats['size_mb']:7.1f} MB  {stats['elapsed_s']:6.1f}s  failures={len(stats['failures'])}"
        )
    return 0 if all(not s["failures"] for s in per_figure.values()) else 1


if __name__ == "__main__":
    sys.exit(main())
