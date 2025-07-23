import json
import os
import re
from datetime import datetime
from pathlib import Path


def analyze_recording_dates_to_json(root_directory, output_file):
    all_sessions = []

    for dirpath, dirnames, filenames in os.walk(root_directory):
        for dirname in [d for d in dirnames if d.startswith("BOT_")]:
            session_folder = Path(os.path.join(dirpath, dirname))
            xml_file = session_folder / f"{session_folder.name}.xml"

            if xml_file.exists():
                try:
                    with open(xml_file, "r", encoding="utf-8") as f:
                        content = f.read()
                        match = re.search(r'<PVScan.*?date="(.*?)"', content)
                        if match:
                            date_str = match.group(1)
                            try:
                                dt_object = datetime.strptime(date_str, "%m/%d/%Y %I:%M:%S %p")
                                all_sessions.append(
                                    {
                                        "path": str(session_folder.relative_to(root_directory)),
                                        "session_start_time": dt_object,
                                    }
                                )
                            except ValueError:
                                print(f"Could not parse date string: {date_str} in file: {xml_file}")
                except Exception as e:
                    print(f"Error reading file {xml_file}: {e}")

    sorted_sessions = sorted(all_sessions, key=lambda x: x["session_start_time"])

    output_json = {}
    for session in sorted_sessions:
        path_parts = session["path"].split(os.sep)
        current_level = output_json
        for part in path_parts[:-1]:
            current_level = current_level.setdefault(part, {})

        session_name = path_parts[-1]
        current_level[session_name] = {
            "session_start_time": session["session_start_time"].strftime("%Y-%m-%d %H:%M:%S")
        }

    with open(output_file, "w") as f:
        json.dump(output_json, f, indent=4)


parent_directory = "/media/heberto/One Touch/Surmeier-CN-data-share/consolidated_data/LID_paper_Zhai_2025/Raw data for Figs/Figure 5_SF2"
output_file_path = (
    "/home/heberto/development/surmeier-lab-to-nwb/figure_scripts/acetylcholien_biosensor_session_times.json"
)
analyze_recording_dates_to_json(parent_directory, output_file_path)
