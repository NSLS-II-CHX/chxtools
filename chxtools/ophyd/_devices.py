from enum import Enum
from typing import ClassVar, TypedDict

from area_detectors.eiger import EigerFastTrigger
from area_detectors.timepix import TimePixDetector
from beam_monitors import Elm, XBPM, XBpm, XBpmFeedbackEnable
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
    eigerfasttrigger: ClassVar[DeviceDescriptor] = {
        'object': EigerFastTrigger,
    }
    elm: ClassVar[DeviceDescriptor] = {
        'object': Elm
    }