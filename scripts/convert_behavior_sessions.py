"""Write one NWB file per behaviour mouse and levodopa session (Figures 7 and 8, Supplementary Figure 3).

Outputs go to nwb_files/behavior/<figure>/<subject_id>_<YYYYMMDD>.nwb. The files link to copies of the lab's videos
without audio, made in nwb_files/behavior/<figure>/videos/ with the lab's folder layout.
"""

from __future__ import annotations

from pathlib import Path

from tqdm import tqdm

from surmeier_lab_to_nwb.zhai2025.behavior_sessions import (
    build_behavior_sessions,
    write_behavior_session,
)

OUTPUT_DIR = Path(__file__).resolve().parents[1] / "nwb_files" / "behavior"


def main() -> None:
    for figure in ("Figure 7", "Figure 8"):
        figure_dir = OUTPUT_DIR / figure.replace(" ", "_").lower()
        figure_dir.mkdir(parents=True, exist_ok=True)
        written = set()
        for session in tqdm(build_behavior_sessions(figure), desc=figure, unit=" session"):
            output_path = figure_dir / f"{session.subject.subject_id}_{session.aim.date.strftime('%Y%m%d')}.nwb"
            if output_path in written:
                raise ValueError(f"Two sessions would be written to {output_path}")
            written.add(output_path)
            write_behavior_session(session, output_path, video_folder=figure_dir / "videos")


if __name__ == "__main__":
    main()
