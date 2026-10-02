from datetime import datetime

import bluesky.plan_stubs as bps
import numpy as np
from nslsii.ad33 import SingleTriggerV33
from ophyd import Component as Cpt
from ophyd import Device
from ophyd.areadetector.base import EpicsSignalWithRBV
from ophyd.areadetector.detectors import AreaDetector
from ophyd.areadetector.filestore_mixins import new_short_uid
from ophyd.areadetector.paths import EpicsPathSignal
from ophyd.areadetector.plugins import (
    ROIPlugin_V34,
    StatsPlugin_V34,
    TimeSeriesPlugin_V34,
)
from ophyd.signal import EpicsSignal, Signal


class PreciseDtypeSignal(Signal):
    def describe(self):
        ret = super().describe()
        ret[self.name].setdefault('dtype_str', str(np.asarray(self.get()).dtype))
        return ret

class Tpx3Files(Device):
    # The other PVs only take affect after this is set
    set_settings = Cpt(EpicsSignal, "WriteData", kind="omitted")

    raw_filepath = Cpt(
        EpicsPathSignal,
        "RawFilePath",
        string=True,
        kind="config",
        path_semantics="posix",
    )
    raw_file_template = Cpt(
        EpicsSignalWithRBV, "RawFileTemplate", string=True, kind="config"
    )
    raw_write_enable = Cpt(EpicsSignalWithRBV, "WriteRaw", string=True, kind="omitted", put_complete=True)

    img_filepath = Cpt(EpicsSignalWithRBV, "ImgFilePath", kind="config")
    img_file_template = Cpt(EpicsSignalWithRBV, "ImgFileTemplate", kind="config")
    img_write_enable = Cpt(EpicsSignalWithRBV, "WriteImg", kind="omitted")

    prv_filepath = Cpt(EpicsSignal, "PrvImgFilePath", kind="config")
    prv_file_template = Cpt(EpicsSignal, "PrvImgFileTemplate", kind="config")

    prv1_filepath = Cpt(EpicsSignal, "PrvImg1FilePath", kind="config")

    # HACK UNTIL WE GET HANDLERS SET UP
    raw_filepaths = Cpt(PreciseDtypeSignal, kind="normal")

    def __init__(self, *args, **kwargs):
        self.sequence_id_offset = 1
        # This is changed for when a datum is a slice
        # also used by ophyd
        self.filestore_spec = "TPX3_RAW"
        self.frame_num = None
        super().__init__(*args, **kwargs)
        self._datum_kwargs_map = {}  # store kwargs for each uid
        self._n = 0

        self.stage_sigs['img_write_enable'] = 'false'

    def stage(self):
        # TODO also do the images
        
        self._res_uid = res_uid = new_short_uid()
        write_path_template = "file:" + self.assets_path() + f"{self.parent.assets_name}/%Y/%m/%d/"
        self._write_path = write_path = datetime.now().strftime(write_path_template)
        self.raw_filepath.set(write_path).wait()

        # TODO check what the % formatting means to the server
        self.raw_file_template.set(f"{res_uid}_0").wait()

        # because we need to flush setting to actual server from IOC
        self.raw_write_enable.set(1).wait()
        self.set_settings.set(1).wait()
        self.set_settings.set(1).wait()

        # fill in first guess, this ill increment self._n but we really do not want to!
        self.update_file_template()
        # reset this back to 0!  
        self._n = 0
        
        super().stage()

    def update_file_template(self):
        # The server generates files with in a trigger that are formatted as
        # filepath/template{datetime only}_{j:d6}.tpx3
        # however the counting is only within the trigger and because we want to be able to predict the file names we can not
        # rely on the date formatting, this we need to set this template on every trigger

        # TODO check what the % formatting means to the server
        self.raw_file_template.set(f"{self._res_uid}_{self._n:05d}_").wait()
        # because we need to flush setting to actual server from IOC
        self.set_settings.set(1).wait()

        filenames = [f'{self._write_path}{self._res_uid}_{self._n:05d}_{j:06d}.tpx3' for j in range(self.parent.cam.num_images.get())]
        self._n += 1
        self.raw_filepaths.set(filenames).wait()

    def unstage(self) -> list[object]:
        # TODO reset these to their original values rather than junk
        self.raw_filepath.set('file:/media/nvme/raw/').wait()
        self.raw_file_template.set("garbage").wait()
        self.raw_write_enable.set(0).wait()
        # because we need to flush setting to actual server from IOC
        self.set_settings.set(1).wait()
        return super().unstage()


class Tpx3HDF(Device):
    # Use the HDF5 plugin directory creation feature to avoid having to
    # create directory creation logic for TPX FW explicitly.
    hdf5_file_path = Cpt(EpicsSignalWithRBV, "FilePath", kind="omitted")
    hdf5_create_directory = Cpt(EpicsSignalWithRBV, "CreateDirectory", kind="omitted")
    
    def stage(self):
         self.hdf5_create_directory.set(-4)
         write_path_template = self.assets_path() + f"{self.parent.assets_name}/%Y/%m/%d/"
         write_path = datetime.now().strftime(write_path_template)
         self.hdf5_file_path.put(write_path)


class TimePixDetector(SingleTriggerV33, AreaDetector):
    _default_configuration_attrs = None
    _default_read_attrs = None

    hdf_plugin = Cpt(Tpx3HDF, "HDF1:")

    files = Cpt(Tpx3Files, "cam1:")
    
    stats1 = Cpt(StatsPlugin_V34, "Stats1:")
    stats2 = Cpt(StatsPlugin_V34, "Stats2:")
    stats3 = Cpt(StatsPlugin_V34, "Stats3:")
    stats4 = Cpt(StatsPlugin_V34, "Stats4:")

    roi1 = Cpt(ROIPlugin_V34, "ROI1:")
    roi2 = Cpt(ROIPlugin_V34, "ROI2:")
    roi3 = Cpt(ROIPlugin_V34, "ROI3:")
    roi4 = Cpt(ROIPlugin_V34, "ROI4:")

    ts1 = Cpt(TimeSeriesPlugin_V34, "Stats1:TS:")
    ts2 = Cpt(TimeSeriesPlugin_V34, "Stats2:TS:")
    ts3 = Cpt(TimeSeriesPlugin_V34, "Stats3:TS:")
    ts4 = Cpt(TimeSeriesPlugin_V34, "Stats4:TS:")

    def __init__(self, *args, **kwargs):
        self.assets_name = kwargs.pop("assets_name", "timepix")
        super().__init__(*args, **kwargs)

        for j in range(1, 5):
            stat = getattr(self, f'stats{j}')
            stat.kind = 'normal'
            stat.total.kind = 'hinted'
            stat.ts_total.kind = 'normal'

    def trigger(self):
        self.files.update_file_template()
        for j in range(1, 5):
            getattr(self, f'ts{j}').ts_acquire.set(1).wait()    
        return super().trigger()
    
    def set_num_images(self, n):
        yield from bps.mv(
            self.cam.num_images, n, 
            self.ts1.ts_num_points, n,
            self.ts2.ts_num_points, n,
            self.ts3.ts_num_points, n,
            self.ts4.ts_num_points, n)
        
    def set_total_exposure(self, exp, *, max_exposure=10):
        num_frames = 1 + exp // max_exposure
        real_exp = exp / num_frames
        period = real_exp + .1

        yield from bps.mv(self.cam.acquire_period, period)
        yield from bps.mv(self.cam.acquire_time, real_exp)
        yield from self.set_num_images(num_frames)      
