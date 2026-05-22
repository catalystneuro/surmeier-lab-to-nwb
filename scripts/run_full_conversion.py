"""Run every figure conversion script with --stub-test=False, outer tqdm progress.

Writes a per-script timing+pass/fail log to logs/full_conversion.log next to the
project root. Each script's own internal tqdm bars still print to stdout.

Usage:
    uv run scripts/run_full_conversion.py
"""

from __future__ import annotations

import subprocess
import sys
import time
from pathlib import Path

from tqdm import tqdm

REPO_ROOT = Path(__file__).resolve().parent.parent
LOG_PATH = REPO_ROOT / "logs" / "full_conversion.log"
SCRIPTS: list[tuple[str, str]] = [
    ("dendritic_excitability/figure_1_dendritic_excitability.py", "Fig1 dendritic"),
    ("dendritic_excitability/figure_3_dendritic_excitability.py", "Fig3 dendritic"),
    ("dendritic_excitability/figure_6_dendritic_excitability.py", "Fig6 dendritic"),
    ("dendritic_excitability/figure_7_dendritic_excitability.py", "Fig7 dendritic"),
    ("dendritic_excitability/figure_7_oxoM_dendritic_excitability.py", "Fig7 oxoM dendritic"),
    ("somatic_excitability/figure_1_somatic_excitability.py", "Fig1 somatic"),
    ("somatic_excitability/figure_3_somatic_excitability.py", "Fig3 somatic"),
    ("somatic_excitability/figure_6_somatic_excitability.py", "Fig6 somatic"),
    ("somatic_excitability/figure_7_somatic_excitability.py", "Fig7 somatic"),
    ("somatic_excitability/figure_8_somatic_excitability.py", "Fig8 somatic"),
    ("spine_density/figure_2_spine_density.py", "Fig2 spine"),
    ("spine_density/figure_4_spine_density.py", "Fig4 spine"),
    ("spine_density/figure_6_spine_density.py", "Fig6 spine"),
    ("spine_density/figure_7_spine_density.py", "Fig7 spine"),
    ("spine_density/figure_8_spine_density.py", "Fig8 spine"),
    ("confocal_spine_density/figure_4_confocal_spine_density_nikon.py", "Fig4 confocal Nikon"),
    ("confocal_spine_density/figure_4h_confocal_spine_density_olympus.py", "Fig4H Olympus"),
    ("aim_behavior/figure_7_behavioral_aim_experiments.py", "Fig7 AIM"),
    ("aim_behavior/figure_8_behavioral_aim_experiments.py", "Fig8 AIM"),
    ("videos/figure_7_behavioral_videos.py", "Fig7 videos"),
    ("videos/supplementary_figure_3_behavioral_videos.py", "SF3 videos"),
    ("acetylcholine_biosensor/figure_5_acetylcholine_biosensor.py", "Fig5 biosensor"),
    ("optical_stimulation/figure_2_optical_stimuli.py", "Fig2 optical"),
    ("optical_stimulation/figure_4_optical_stimuli.py", "Fig4 optical"),
]


def main() -> int:
    LOG_PATH.parent.mkdir(parents=True, exist_ok=True)
    overall_start = time.time()
    results: list[tuple[str, str, float, str]] = []

    with LOG_PATH.open("w") as log:
        log.write(f"# Full conversion started at {time.strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        log.flush()

        with tqdm(SCRIPTS, desc="overall", unit=" script", position=0) as bar:
            for rel, label in bar:
                bar.set_postfix_str(label)
                script_path = REPO_ROOT / "src" / "surmeier_lab_to_nwb" / "zhai2025" / "conversion_scripts" / rel
                t0 = time.time()
                proc = subprocess.run(
                    [sys.executable, str(script_path), "--stub-test=False"],
                    capture_output=True,
                    text=True,
                    cwd=REPO_ROOT,
                )
                dt = time.time() - t0
                status = "PASS" if proc.returncode == 0 else "FAIL"
                last_err = ""
                if proc.returncode != 0:
                    last_err = " | ".join(
                        line.strip() for line in proc.stderr.strip().splitlines()[-5:] if line.strip()
                    )[:500]
                results.append((label, status, dt, last_err))
                log.write(f"{status:5}  {dt:7.1f}s  {label}\n")
                if last_err:
                    log.write(f"       stderr: {last_err}\n")
                log.flush()

        overall_dt = time.time() - overall_start
        passed = sum(1 for _, s, _, _ in results if s == "PASS")
        summary = f"\n# Done in {overall_dt/60:.1f} min. {passed}/{len(results)} scripts passed.\n"
        log.write(summary)
        print(summary)

    return 0 if passed == len(results) else 1


if __name__ == "__main__":
    sys.exit(main())
