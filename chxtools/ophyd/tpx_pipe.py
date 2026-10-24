from typing import Callable
import warnings

from ophyd import Device, Component, DeviceStatus, EpicsSignal, EpicsSignalRO, Signal
from ophyd.areadetector.paths import EpicsPathSignal
from ophyd.status import SubscriptionStatus
from pathlib import Path
from bisect import insort
from bluesky.protocols import Flyable

from .area_detectors.utils import str_to_wav, _asset_undefined_error

class TPX3Pipe(Device):
    sid = Component(EpicsSignal,"sid")
    scan = Component(EpicsSignal,"scan")
    path = Component(
        EpicsPathSignal,
        "path",
        string=True,
        kind="config",
        path_semantics="posix",
    )
    active = Component(EpicsSignal,"active")
    fire = Component(EpicsSignal,"fire")
    file_stream = Component(EpicsSignalRO,"file",string=True)
    files = Component(Signal, value=[])

    def __init__(self, *args,path_provider:Callable[[],str], **kwargs):
        self._detector = kwargs.pop('det',None)
        self.assets_path = path_provider
        super().__init__(*args, **kwargs)

    def stage(self):
        if det := self._detector:
            det.streaming = True
            det.stage()
            
        self._fileset = set()
        self._files = []
        self.sid.put(-1)
        self.path.set(self.assets_path()+"timepix_pipeline/")

        def file_monitor(obj,value=None,old_value=None,**kwargs):
            print(value)
            if value not in self._fileset:
                insort(self._files,value,
                              key=lambda x:int(Path(x).stem.split("_")[3]))
                self._fileset.add(value)
                self.files.put(self._files)

        
        self.sub_id = self.file_stream.subscribe(file_monitor,run=False)

    def await_revert(self,old_value=None,value=None,**kwargs):
        return (not value)

    def trigger(self):
        self._fileset = set()
        self._files = []
        self.files.put(self._files)
        if det := self._detector:
            det_stat = det.trigger()
        self.fire.put(True)
        status = SubscriptionStatus(self.fire,self.await_revert,run=False)
        return status
    
    def kickoff(self):
        self._fileset = set()
        self._files = []
        self.files.put(self._files)
        self.fire.put(True)
        return SubscriptionStatus(self.fire,lambda *args, **kwargs: not self.await_revert(*args,**kwargs),run=False)
        
    def complete(self):
        return SubscriptionStatus(self.fire,self.await_revert,run=False)
        

    def unstage(self):
        if det := self._detector:
            det.streaming = False
            det.stage()
        self.file_stream.unsubscribe(self.sub_id)
