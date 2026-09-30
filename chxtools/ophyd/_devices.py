from enum import Enum
from typing import ClassVar, TypedDict

from area_detectors.eiger import EigerFastTrigger
from area_detectors.timepix import TimePixDetector
from ophyd import Device


# Define the static structure
class DeviceDescriptor(TypedDict):
    object: Device
    name: str | None = None
    init_params: dict | None = None
    meta_structure: dict | None = None

class Devices(Enum):
    timepix: ClassVar[DeviceDescriptor] = {
        'object': TimePixDetector,
    }
    EigerFastTrigger: ClassVar[DeviceDescriptor] = {
        'object': EigerFastTrigger,
    }