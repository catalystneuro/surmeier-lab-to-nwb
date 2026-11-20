"""The dataset's metadata, read from metadata.yaml."""

from copy import deepcopy
from functools import cache
from pathlib import Path

from neuroconv.utils import load_dict_from_file

METADATA_PATH = Path(__file__).parent / "metadata.yaml"


@cache
def _load() -> dict:
    return load_dict_from_file(METADATA_PATH)


def shared_metadata() -> dict:
    """The NWBFile, Subject and ExternalResources metadata that every file shares."""
    return deepcopy(
        {key: value for key, value in _load().items() if key not in ("DeviceModels", "Devices", "Sessions")}
    )


def device_metadata(key: str) -> dict:
    """The name and description of a device that sessions of more than one kind add, by its metadata key."""
    return deepcopy(_load()["Devices"][key])


def device_model_metadata(name: str) -> dict:
    """The manufacturer, model number and description of a device model, by its name."""
    return deepcopy(_load()["DeviceModels"][name])


def session_metadata(key: str) -> dict:
    """What one kind of session adds to the shared metadata, by its metadata key; fields in braces are templates."""
    return deepcopy(_load()["Sessions"][key])
