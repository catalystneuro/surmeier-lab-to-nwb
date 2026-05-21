"""Shared DeviceModel factories for the Zhai 2025 conversion.

In pynwb 3.1+, the `Device.manufacturer`, `Device.model_number`, and
`Device.model_name` attributes are deprecated in favor of linking each Device
to a `DeviceModel` via the `Device.model` link. This module centralizes the
DeviceModel objects used across multiple interfaces and scripts so the model
metadata stays consistent and only lives in one place.

Each factory is idempotent: it returns the existing DeviceModel if one is
already present on the NWBFile, otherwise it creates and registers one.
"""

from pynwb import NWBFile
from pynwb.device import DeviceModel


def get_or_create_bruker_ultima_model(nwbfile: NWBFile) -> DeviceModel:
    """Return the Bruker Ultima In Vitro multiphoton microscope model.

    Shared across line-scan dendritic excitability (`ophys_interfaces`),
    Brightness-Over-Time biosensor imaging (`prairie_view_fluorescence_interface`),
    and two-photon spine-density acquisitions (`spine_density_utils`). The lab's
    paper identifies the same physical Bruker Ultima In Vitro rig mounted on an
    Olympus BX-51 chassis for all three.
    """
    name = "BrukerUltimaInVitroModel"
    if name in nwbfile.device_models:
        return nwbfile.device_models[name]
    model = DeviceModel(
        name=name,
        manufacturer="Bruker Corporation",
        model_number="Ultima In Vitro",
        description=(
            "Bruker Ultima In Vitro multiphoton microscope system (Bruker, Billerica, MA), "
            "mounted on an Olympus BX-51 upright microscope chassis for ex vivo two-photon "
            "laser-scanning imaging. Used in conjunction with patch-clamp electrophysiology "
            "and a Chameleon Ultra II (Coherent) tunable Ti:sapphire laser."
        ),
    )
    nwbfile.add_device_model(model)
    return model


def get_or_create_nikon_axr_model(nwbfile: NWBFile) -> DeviceModel:
    """Return the Nikon AXR confocal laser microscope model.

    Used for Figure 4 confocal spine density imaging of sparsely labeled iSPNs.
    The paper identifies the AXR specifically and credits the Center for Advanced
    Microscopy and Nikon Imaging Center at Northwestern University.
    """
    name = "NikonAXRModel"
    if name in nwbfile.device_models:
        return nwbfile.device_models[name]
    model = DeviceModel(
        name=name,
        manufacturer="Nikon",
        model_number="AXR",
        description=(
            "Nikon AXR confocal laser microscope from the Center for Advanced Microscopy "
            "and Nikon Imaging Center at Northwestern University. Used for high-resolution "
            "confocal imaging of dendritic spines in sparsely labeled iSPNs at 60x oil "
            "immersion (NA=1.49) with 0.09 um pixels and 0.125 um z-steps."
        ),
    )
    nwbfile.add_device_model(model)
    return model


def get_or_create_olympus_fv10i_model(nwbfile: NWBFile) -> DeviceModel:
    """Return the Olympus FV10i-DUC confocal laser scanning microscope model.

    Used for Figure 4H confocal spine density imaging. The paper identifies the
    system specifically as the FV10i-DUC.
    """
    name = "OlympusFV10iDUCModel"
    if name in nwbfile.device_models:
        return nwbfile.device_models[name]
    model = DeviceModel(
        name=name,
        manufacturer="Olympus",
        model_number="FV10i-DUC",
        description=(
            "Olympus FV10i-DUC confocal laser scanning microscope. Used for high-resolution "
            "confocal imaging of dendritic spines at 60x magnification (NA=1.35) with "
            "0.125 um z-steps and ~0.207 um pixels for spine density methodology validation."
        ),
    )
    nwbfile.add_device_model(model)
    return model
