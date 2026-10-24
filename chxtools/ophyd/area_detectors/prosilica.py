
from typing import Callable

from nslsii.ad33 import SingleTriggerV33, StatsPluginV33
from ophyd import Component as Cpt
from ophyd import (
    EpicsSignal,
    ImagePlugin,
    OverlayPlugin,
    PointGreyDetector,
    PointGreyDetectorCam,
    ProcessPlugin,
    ProsilicaDetector,
    ProsilicaDetectorCam,
    ROIPlugin,
    SingleTrigger,
    StatsPlugin,
    TransformPlugin,
)

from chxtools.ophyd.area_detectors.plugins import (
    TIFFPluginEnsuredOff,
    TIFFPluginWithFileStore,
)
from chxtools.ophyd.area_detectors.utils import name_dir_mapping


class ProsilicaDetectorCamV33(ProsilicaDetectorCam):
    '''This is used to update the Standard Prosilica to AD33. It adds the
process
    '''
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

class StandardProsilica(SingleTrigger, ProsilicaDetector):
    image = Cpt(ImagePlugin, 'image1:')
    stats1 = Cpt(StatsPlugin, 'Stats1:')
    stats2 = Cpt(StatsPlugin, 'Stats2:')
    stats3 = Cpt(StatsPlugin, 'Stats3:')
    stats4 = Cpt(StatsPlugin, 'Stats4:')
    stats5 = Cpt(StatsPlugin, 'Stats5:')
    trans1 = Cpt(TransformPlugin, 'Trans1:')
    roi1 = Cpt(ROIPlugin, 'ROI1:')
    roi2 = Cpt(ROIPlugin, 'ROI2:')
    roi3 = Cpt(ROIPlugin, 'ROI3:')
    roi4 = Cpt(ROIPlugin, 'ROI4:')
    proc1 = Cpt(ProcessPlugin, 'Proc1:')
    over1 = Cpt(OverlayPlugin, 'Over1:')

    # This class does not save TIFFs. We make it aware of the TIFF plugin
    # only so that it can ensure that the plugin is not auto-saving.
    tiff = Cpt(TIFFPluginEnsuredOff, suffix='TIFF1:')

    def __init__(self, *args,path_provider: Callable[[],str], **kwargs):
        self.assets_path = path_provider
        super().__init__(*args, **kwargs)

    @property
    def hints(self):
        return {'fields': [self.stats1.total.name]}

class StandardProsilicaV33(SingleTriggerV33, ProsilicaDetector):
    cam = Cpt(ProsilicaDetectorCamV33, 'cam1:')
    image = Cpt(ImagePlugin, 'image1:')
    stats1 = Cpt(StatsPluginV33, 'Stats1:')
    stats2 = Cpt(StatsPluginV33, 'Stats2:')
    stats3 = Cpt(StatsPluginV33, 'Stats3:')
    stats4 = Cpt(StatsPluginV33, 'Stats4:')
    stats5 = Cpt(StatsPluginV33, 'Stats5:')
    trans1 = Cpt(TransformPlugin, 'Trans1:')
    roi1 = Cpt(ROIPlugin, 'ROI1:')
    roi2 = Cpt(ROIPlugin, 'ROI2:')
    roi3 = Cpt(ROIPlugin, 'ROI3:')
    roi4 = Cpt(ROIPlugin, 'ROI4:')
    proc1 = Cpt(ProcessPlugin, 'Proc1:')
    over1 = Cpt(OverlayPlugin, 'Over1:')

    # This class does not save TIFFs. We make it aware of the TIFF plugin
    # only so that it can ensure that the plugin is not auto-saving.
    tiff = Cpt(TIFFPluginEnsuredOff, suffix='TIFF1:')

    def __init__(self, *args,path_provider: Callable[[],str], **kwargs):
        self.assets_path = path_provider
        super().__init__(*args, **kwargs)

    @property
    def hints(self):
        return {'fields': [self.stats1.total.name]}


class StandardProsilicaWithTIFF(StandardProsilica):
    tiff = Cpt(TIFFPluginWithFileStore,
               suffix='TIFF1:',
               write_path_template='')
    
    def stage(self, *args, **kwargs):
        self.tiff.write_path_template = self.assets_path() + f'{name_dir_mapping[self.name]}/%Y/%m/%d/'
        self.tiff.read_path_template = () + f'{name_dir_mapping[self.name]}/%Y/%m/%d/'
        self.tiff.reg_root = self.assets_path() + f'{name_dir_mapping[self.name]}'
        return super().stage(*args, **kwargs)

class StandardProsilicaWithTIFFV33(StandardProsilicaV33):
    tiff = Cpt(TIFFPluginWithFileStore,
               suffix='TIFF1:',
               write_path_template='')
    
    def stage(self, *args, **kwargs):
        self.tiff.write_path_template = self.assets_path() + f'{name_dir_mapping[self.name]}/%Y/%m/%d/'
        self.tiff.read_path_template = self.assets_path() + f'{name_dir_mapping[self.name]}/%Y/%m/%d/'
        self.tiff.reg_root = self.assets_path() + f'{name_dir_mapping[self.name]}'
        return super().stage(*args, **kwargs)


class PointGreyDetectorCamV33(PointGreyDetectorCam):
    '''This is used to update the Standard PointGrey to AD33. It adds the
process
    '''
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

class StandardPointGreyV33(SingleTriggerV33, PointGreyDetector):
    cam = Cpt(PointGreyDetectorCamV33, 'cam1:')
    image = Cpt(ImagePlugin, 'image1:')
    stats1 = Cpt(StatsPluginV33, 'Stats1:')
    stats2 = Cpt(StatsPluginV33, 'Stats2:')
    stats3 = Cpt(StatsPluginV33, 'Stats3:')
    stats4 = Cpt(StatsPluginV33, 'Stats4:')
    stats5 = Cpt(StatsPluginV33, 'Stats5:')
    trans1 = Cpt(TransformPlugin, 'Trans1:')
    roi1 = Cpt(ROIPlugin, 'ROI1:')
    roi2 = Cpt(ROIPlugin, 'ROI2:')
    roi3 = Cpt(ROIPlugin, 'ROI3:')
    roi4 = Cpt(ROIPlugin, 'ROI4:')
    proc1 = Cpt(ProcessPlugin, 'Proc1:')
    over1 = Cpt(OverlayPlugin, 'Over1:')

    # This class does not save TIFFs. We make it aware of the TIFF plugin
    # only so that it can ensure that the plugin is not auto-saving.
    tiff = Cpt(TIFFPluginEnsuredOff, suffix='TIFF1:')

    def __init__(self, *args,path_provider: Callable[[],str], **kwargs):
        self.assets_path = path_provider
        super().__init__(*args, **kwargs)

    @property
    def hints(self):
        return {'fields': [self.stats1.total.name]}


class StandardPointGreyWithTIFFV33(StandardPointGreyV33):
    tiff = Cpt(TIFFPluginWithFileStore,
               suffix='TIFF1:',
               write_path_template='')
    
    def stage(self, *args, **kwargs):
        self.tiff.write_path_template = self.assets_path() + f'{name_dir_mapping[self.name]}/%Y/%m/%d/'
        self.tiff.read_path_template = self.assets_path() + f'{name_dir_mapping[self.name]}/%Y/%m/%d/'
        self.tiff.reg_root = self.assets_path() + f'{name_dir_mapping[self.name]}'
        return super().stage(*args, **kwargs)

