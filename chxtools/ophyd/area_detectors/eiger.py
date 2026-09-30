import time as ttime  # tea time
from collections import OrderedDict
from datetime import datetime
from pathlib import PurePath
from types import SimpleNamespace

from nslsii.ad33 import StatsPluginV33
from ophyd import (
    AreaDetector,
    Device,
    EpicsSignal,
    EpicsSignalRO,
    ImagePlugin,
    ProcessPlugin,
    ROIPlugin,
    Signal,
    SingleTrigger,
    StatsPlugin,
)
from ophyd import Component as Cpt
from ophyd.areadetector.base import ADComponent, EpicsSignalWithRBV
from ophyd.areadetector.cam import AreaDetectorCam
from ophyd.areadetector.filestore_mixins import (
    FileStoreBase,
    new_short_uid,
)
from ophyd.device import Staged
from ophyd.status import StatusBase
from ophyd.utils import set_and_wait

from chxtools.ophyd.area_detectors.utils import name_dir_mapping


class EigerSimulatedFilePlugin(Device, FileStoreBase):
    sequence_id = ADComponent(EpicsSignalRO, 'SequenceId')
    file_path = ADComponent(EpicsSignalWithRBV, 'FilePath', string=True)
    create_dir_depth = ADComponent(EpicsSignalWithRBV, "CreateDirectory")
    file_write_name_pattern = ADComponent(EpicsSignalWithRBV, 'FWNamePattern',
                                          string=True)
    file_write_images_per_file = ADComponent(EpicsSignalWithRBV,
                                             'FWNImagesPerFile')
    current_run_start_uid = Cpt(Signal, value='', add_prefix=())
    enable = SimpleNamespace(get=lambda: True)

    def __init__(self, *args, **kwargs):
        self.sequence_id_offset = 1
        # This is changed for when a datum is a slice
        # also used by ophyd
        self.filestore_spec = "AD_EIGER2"
        self.frame_num = None
        super().__init__(*args, **kwargs)
        self._datum_kwargs_map = {}  # store kwargs for each uid

    def stage(self):
        res_uid = new_short_uid()
        self.create_dir_depth.set(-4).wait()
        write_path = datetime.now().strftime(self.write_path_template)
        #set_and_wait(self.file_path, write_path + '/')
        self.file_path.set(write_path + '/').wait()
        #set_and_wait(self.file_write_name_pattern, '{}_$id'.format(res_uid))
        self.file_write_name_pattern.set(f'{res_uid}_$id').wait()
        super().stage()
        fn = (PurePath(self.file_path.get()) / res_uid)
        ipf = int(self.file_write_images_per_file.get())
        # logger.debug("Inserting resource with filename %s", fn)
        self._fn = fn
        res_kwargs = {'images_per_file' : ipf}
        self._generate_resource(res_kwargs)

    def generate_datum(self, key, timestamp, datum_kwargs):
        # The detector keeps its own counter which is uses label HDF5
        # sub-files.  We access that counter via the sequence_id
        # signal and stash it in the datum.
        seq_id = int(self.sequence_id_offset) + int(self.sequence_id.get())  # det writes to the NEXT one
        datum_kwargs.update({'seq_id': seq_id})
        if self.frame_num is not None:
            datum_kwargs.update({'frame_num': self.frame_num})
        return super().generate_datum(key, timestamp, datum_kwargs)

    def describe(self,):
        ret = super().describe()
        if hasattr(self.parent.cam, 'bit_depth'):
            cur_bits = self.parent.cam.bit_depth.get()
            dtype_str_map = {8: '|u1', 16: '<u2', 32:'<u4'}
            ret[self.parent._image_name]['dtype_str'] = dtype_str_map[cur_bits]
        return ret


class EigerBase(AreaDetector):
    """
    Eiger, sans any triggering behavior.

    Use EigerSingleTrigger or EigerFastTrigger below.
    """
    num_triggers = ADComponent(EpicsSignalWithRBV, 'cam1:NumTriggers')
    file = Cpt(EigerSimulatedFilePlugin, suffix='cam1:',
               write_path_template='')
    beam_center_x = ADComponent(EpicsSignalWithRBV, 'cam1:BeamX')
    beam_center_y = ADComponent(EpicsSignalWithRBV, 'cam1:BeamY')
    wavelength = ADComponent(EpicsSignalWithRBV, 'cam1:Wavelength')
    det_distance = ADComponent(EpicsSignalWithRBV, 'cam1:DetDist')
    threshold_energy = ADComponent(EpicsSignalWithRBV, 'cam1:ThresholdEnergy')
    photon_energy = ADComponent(EpicsSignalWithRBV, 'cam1:PhotonEnergy')
    manual_trigger = ADComponent(EpicsSignalWithRBV, 'cam1:ManualTrigger')  # the checkbox
    special_trigger_button = ADComponent(EpicsSignal, 'cam1:Trigger')  # the button next to 'Start' and 'Stop'
    stream_enable = ADComponent(EpicsSignal,'cam1:StreamEnable')
    data_source = ADComponent(EpicsSignal,'cam1:DataSource') # data source like None, FileWriter, Stream
    image = Cpt(ImagePlugin, 'image1:')
    stats1 = Cpt(StatsPlugin, 'Stats1:')
    stats2 = Cpt(StatsPlugin, 'Stats2:')
    stats3 = Cpt(StatsPlugin, 'Stats3:')
    stats4 = Cpt(StatsPlugin, 'Stats4:')
    stats5 = Cpt(StatsPlugin, 'Stats5:')
    roi1 = Cpt(ROIPlugin, 'ROI1:')
    roi2 = Cpt(ROIPlugin, 'ROI2:')
    roi3 = Cpt(ROIPlugin, 'ROI3:')
    roi4 = Cpt(ROIPlugin, 'ROI4:')
    proc1 = Cpt(ProcessPlugin, 'Proc1:')

    shutter_mode = ADComponent(EpicsSignalWithRBV, 'cam1:ShutterMode')

    # hotfix: shadow non-existant PV
    size_link = None

    def stage(self, *args, **kwargs):
        # before parent
        ret = super().stage(*args, **kwargs)
        # after parent
        self.manual_trigger.set(1).wait()
        self.file.write_path_template = self.assets_path() + f'{name_dir_mapping[self.name]}/%Y/%m/%d/'
        self.file.reg_root = self.assets_path() + f'{name_dir_mapping[self.name]}'
        return ret

    def unstage(self):
        #set_and_wait(self.manual_trigger, 0)
        self.manual_trigger.set(0).wait()
        super().unstage()

    @property
    def hints(self):
        return {'fields': [self.stats1.total.name]}


class EigerDetectorCamV33(AreaDetectorCam):
    '''This is used to update the Eiger detector to AD33.
    '''
    firmware_version = Cpt(EpicsSignalRO, 'FirmwareVersion_RBV', kind='config')

    wait_for_plugins = Cpt(EpicsSignal, 'WaitForPlugins',
                           string=True, kind='config')

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.stage_sigs['wait_for_plugins'] = 'Yes'

    def ensure_nonblocking(self):
        self.stage_sigs['wait_for_plugins'] = 'Yes'
        for c in self.parent.component_names:
            cpt = getattr(self.parent, c)
            if cpt is self:
                continue
            if hasattr(cpt, 'ensure_nonblocking'):
                cpt.ensure_nonblocking()

class NewEigerDetectorCamV33(EigerDetectorCamV33):
    bit_depth = Cpt(EpicsSignalRO, 'BitDepthImage_RBV', kind='config')

class EigerBaseV33(EigerBase):
    cam = Cpt(EigerDetectorCamV33, 'cam1:')
    stats1 = Cpt(StatsPluginV33, 'Stats1:')
    stats2 = Cpt(StatsPluginV33, 'Stats2:')
    stats3 = Cpt(StatsPluginV33, 'Stats3:')
    stats4 = Cpt(StatsPluginV33, 'Stats4:')
    stats5 = Cpt(StatsPluginV33, 'Stats5:')


class EigerSingleTrigger(SingleTrigger, EigerBase):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.stage_sigs['cam.trigger_mode'] = 0
        self.stage_sigs['shutter_mode'] = 1  # 'EPICS PV'
        self.stage_sigs.update({'num_triggers': 1})

    def stage(self, *args, **kwargs):
        return super().stage(*args, **kwargs)

    def trigger(self, *args, **kwargs):
        status = super().trigger(*args, **kwargs)
        set_and_wait(self.special_trigger_button, 1)
        return status

    def read(self, *args, streaming=False, **kwargs):
        '''
            This is a test of using streaming read.
            Ideally, this should be handled by a new _stream_attrs property.
            For now, we just check for a streaming key in read and
            call super() if False, or read the one key we know we should read
            if True.

            Parameters
            ----------
            streaming : bool, optional
                whether to read streaming attrs or not
        '''
        #ret = super().read()
        #print("super read() : {}".format(ret))
        #return ret
        if streaming:
            key = self._image_name  # this comes from the SingleTrigger mixin
            read_dict = super().read()
            ret = OrderedDict({key: read_dict[key]})
            return ret
        else:
            ret = super().read(*args, **kwargs)
            return ret

    def describe(self, *args, streaming=False, **kwargs):
        '''
            This is a test of using streaming read.
            Ideally, this should be handled by a new _stream_attrs property.
            For now, we just check for a streaming key in read and
            call super() if False, or read the one key we know we should read
            if True.

            Parameters
            ----------
            streaming : bool, optional
                whether to read streaming attrs or not
        '''
        if streaming:
            key = self._image_name  # this comes from the SingleTrigger mixin
            read_dict = super().describe()
            ret = OrderedDict({key: read_dict[key]})
            return ret
        else:
            ret = super().describe(*args, **kwargs)
            return ret



class EigerSingleTrigger_AD37(SingleTrigger, EigerBaseV33):
    cam = Cpt(NewEigerDetectorCamV33, 'cam1:')
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.auto_trigger = True
        self.stage_sigs['cam.trigger_mode'] = 0
        self.stage_sigs['shutter_mode'] = 1  # 'EPICS PV'
        self.stage_sigs.update({'num_triggers': 1})

    def stage(self, *args, **kwargs):
        self.file.write_path_template = self.assets_path() + f'{name_dir_mapping[self.name]}/%Y/%m/%d/'
        self.file.reg_root = self.assets_path() + f'{name_dir_mapping[self.name]}'
        return super().stage(*args, **kwargs)

    def trigger(self, *args, **kwargs):
        status = super().trigger(*args, **kwargs)
        #set_and_wait(self.special_trigger_button, 1)
        if self.auto_trigger:
            self.special_trigger_button.set(1).wait()
        return status

    def read(self, *args, streaming=False, **kwargs):
        '''
            This is a test of using streaming read.
            Ideally, this should be handled by a new _stream_attrs property.
            For now, we just check for a streaming key in read and
            call super() if False, or read the one key we know we should read
            if True.

            Parameters
            ----------
            streaming : bool, optional
                whether to read streaming attrs or not
        '''
        #ret = super().read()
        #print("super read() : {}".format(ret))
        #return ret
        if streaming:
            key = self._image_name  # this comes from the SingleTrigger mixin
            read_dict = super().read()
            ret = OrderedDict({key: read_dict[key]})
            return ret
        else:
            ret = super().read(*args, **kwargs)
            return ret

    def describe(self, *args, streaming=False, **kwargs):
        '''
            This is a test of using streaming read.
            Ideally, this should be handled by a new _stream_attrs property.
            For now, we just check for a streaming key in read and
            call super() if False, or read the one key we know we should read
            if True.

            Parameters
            ----------
            streaming : bool, optional
                whether to read streaming attrs or not
        '''
        if streaming:
            key = self._image_name  # this comes from the SingleTrigger mixin
            read_dict = super().describe()
            ret = OrderedDict({key: read_dict[key]})
            return ret
        else:
            ret = super().describe(*args, **kwargs)
            return ret

class FastShutterTrigger(Device):
    """This represents the fast trigger *device*.

    See below, FastTriggerMixin, which defines the trigging logic.
    """
    auto_shutter_mode = Cpt(EpicsSignal, 'Mode-Sts', write_pv='Mode-Cmd')
    num_images = Cpt(EpicsSignal, 'NumImages-SP')
    exposure_time = Cpt(EpicsSignal, 'ExposureTime-SP')
    acquire_period = Cpt(EpicsSignal, 'AcquirePeriod-SP')
    acquire = Cpt(EpicsSignal, 'Acquire-Cmd', trigger_value=1)


class EigerFastTrigger(EigerBaseV33):
    tr = Cpt(FastShutterTrigger, '')

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.stage_sigs['cam.trigger_mode'] = 3  # 'External Enable' mode
        self.stage_sigs['shutter_mode'] = 0  # 'EPICS PV'
        self.stage_sigs['tr.auto_shutter_mode'] = 1  # 'Enable'

    def trigger(self):
        self.dispatch('image', ttime.time())
        return self.tr.trigger()

class EigerManualTrigger(SingleTrigger, EigerBase):
    '''
        Like Eiger Single Trigger but the triggering is done through the
        special trigger button.
    '''
    def __init__(self, *args, **kwargs):
        self._set_st = None
        super().__init__(*args, **kwargs)
        # in this case, we don't need this + 1 sequence id cludge
        # this is because we write datum after image is acquired
        self.file.sequence_id_offset = 0
        self.file.filestore_spec = "AD_EIGER_SLICE"
        self.file.frame_num = 0
        # set up order
        self.stage_sigs = OrderedDict()
        self.stage_sigs['cam.image_mode'] = 1
        self.stage_sigs['cam.trigger_mode'] = 0
        self.stage_sigs['shutter_mode'] = 1
        self.stage_sigs['manual_trigger'] = 1
        #self.stage_sigs['cam.acquire'] = 1
        self.stage_sigs['num_triggers'] = 10

        # monkey patch
        # override with special trigger button, not acquire
        #self._acquisition_signal = self.special_trigger_button

    def stage(self):
        self.file.frame_num = 0
        super().stage()
        # for some reason if doing this too fast in staging
        # this gets reset. so I do it here instead.
        # the bit gets unset when done but we should unset again in unstage
        ttime.sleep(1)
        self.cam.acquire.put(1)
        # need another sleep to ensure the sequence ID is updated
        ttime.sleep(1)

    def unstage(self):
        self.file.frame_num = 0
        super().unstage()
        ttime.time.sleep(.1)
        self.cam.acquire.put(0)


    def trigger(self):
        ''' custom trigger for Eiger Manual'''
        if self._staged != Staged.yes:
            raise RuntimeError("This detector is not ready to trigger."
                               "Call the stage() method before triggering.")

        if self._set_st is not None:
            raise RuntimeError(f'trying to set {self.name}'
                               ' while a set is in progress')

        st = StatusBase()
        # idea : look at the array counter and the trigger value
        # the logic here is that I want to check when the status is back to zero
        def counter_cb(value, timestamp, **kwargs):
            # whenevr it counts just move on
            #print("changed : {}".format(value))
            self._set_st = None
            self.cam.array_counter.clear_sub(counter_cb)
            st._finished()


        # first subscribe a callback
        self.cam.array_counter.subscribe(counter_cb, run=False)

        # then call trigger on the PV
        self.special_trigger_button.put(1, wait=False)
        self.dispatch(self._image_name, ttime.time())
        self.file.frame_num += 1

        return st

