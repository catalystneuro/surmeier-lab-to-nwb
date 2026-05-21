"""Shared parsing helpers for Prairie View (Bruker) XML metadata files.

The two Prairie View interfaces in this package (line scan, brightness over time)
both parse the master `.xml` file that Prairie View writes alongside each recording.
This module centralises the parsing so the same logic is not reimplemented in each
interface.

See `obsidian_docs/source_formats/bruker/prairie_view_xml_format.md` for the XML
structure these helpers operate on. The roiextractors `BrukerTiffImagingExtractor`
parses the same XML with lxml; we use xmltodict here because the two custom
acquisition modes Surmeier uses (PVLinescanDefinition, PVBOTs) are not supported by
the upstream extractor, so reusing its instance-bound parser is not an option.
"""

from datetime import datetime
from pathlib import Path
from zoneinfo import ZoneInfo

import xmltodict

# Surmeier lab is in Evanston, IL (Central Time).
SURMEIER_LAB_TIMEZONE = ZoneInfo("America/Chicago")


def parse_prairie_view_xml(file_path: str | Path) -> dict:
    """Load and parse a Prairie View master XML file into a nested dict.

    Parameters
    ----------
    file_path
        Path to the master ``<recording>.xml`` file.

    Returns
    -------
    dict
        xmltodict representation of the document, with the ``PVScan`` element at the
        top level. Use the other helpers in this module to extract specific fields.
    """
    xml_path = Path(file_path)
    with open(xml_path, "r", encoding="utf-8") as f:
        return xmltodict.parse(f.read())


def get_session_start_time(xml_dict: dict, tz: ZoneInfo = SURMEIER_LAB_TIMEZONE) -> datetime:
    """Extract the session start time from a parsed Prairie View XML.

    The ``<PVScan date="...">`` attribute uses US-locale format
    ``M/D/YYYY H:MM:SS AM/PM`` and is recorded in the rig's local timezone.

    Parameters
    ----------
    xml_dict
        Parsed XML as returned by :func:`parse_prairie_view_xml`.
    tz
        Timezone for the returned datetime. Defaults to the Surmeier lab's
        Central Time zone.

    Returns
    -------
    datetime
        Timezone-aware session start time.
    """
    prairie_view_metadata = xml_dict["PVScan"]
    if "@date" not in prairie_view_metadata:
        raise ValueError("Could not find @date attribute in PVScan element")
    date_str = prairie_view_metadata["@date"]
    session_start_time = datetime.strptime(date_str, "%m/%d/%Y %I:%M:%S %p")
    return session_start_time.replace(tzinfo=tz)


def get_session_start_time_from_file(file_path: str | Path, tz: ZoneInfo = SURMEIER_LAB_TIMEZONE) -> datetime:
    """Convenience: parse the XML file and return the session start time."""
    return get_session_start_time(parse_prairie_view_xml(file_path), tz=tz)


def _pv_state_values(xml_dict: dict) -> list:
    """Return the top-level ``PVStateValue`` list from the parsed XML."""
    return xml_dict["PVScan"]["PVStateShard"]["PVStateValue"]


def get_pv_scalar(xml_dict: dict, key: str, default=None):
    """Look up a scalar ``PVStateValue`` by ``@key`` and return its ``@value``."""
    return next(
        (item["@value"] for item in _pv_state_values(xml_dict) if item["@key"] == key),
        default,
    )


def get_pv_indexed(xml_dict: dict, key: str) -> dict[str, str]:
    """Look up an indexed ``PVStateValue`` by ``@key`` and return ``{index: value}``.

    Used for keys like ``micronsPerPixel`` whose children are
    ``<IndexedValue index="XAxis" value="..."/>``.
    """
    state_value = next(
        (item for item in _pv_state_values(xml_dict) if item["@key"] == key),
        None,
    )
    if state_value is None or "IndexedValue" not in state_value:
        raise ValueError(f"Could not find indexed PVStateValue with key={key!r}")
    return {v["@index"]: v["@value"] for v in state_value["IndexedValue"]}
