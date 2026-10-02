from nslsii.ad33 import SingleTriggerV33, StatsPluginV33
from ophyd import Component as Cpt
from ophyd import (
    EpicsSignal,
    ImagePlugin,
    PilatusDetector,
    PilatusDetectorCam,
    ProcessPlugin,
    ROIPlugin,
)

from chxtools.ophyd.area_detectors.plugins import TIFFPluginWithFileStore


class PilatusDetectorCamV33(PilatusDetectorCam):
    """This is used to update the standard prosilica to AD33."""

    wait_for_plugins = Cpt(EpicsSignal, "WaitForPlugins", string=True, kind="config")

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.stage_sigs["wait_for_plugins"] = "Yes"

    def ensure_nonblocking(self):
        self.stage_sigs["wait_for_plugins"] = "Yes"
        for c in self.parent.component_names:
            cpt = getattr(self.parent, c)
            if cpt is self:
                continue
            if hasattr(cpt, "ensure_nonblocking"):
                cpt.ensure_nonblocking()


class PilatusV33(SingleTriggerV33, PilatusDetector):
    cam = Cpt(PilatusDetectorCamV33, "cam1:")
    image = Cpt(ImagePlugin, "image1:")
    stats1 = Cpt(StatsPluginV33, "Stats1:")
    stats2 = Cpt(StatsPluginV33, "Stats2:")
    stats3 = Cpt(StatsPluginV33, "Stats3:")
    stats4 = Cpt(StatsPluginV33, "Stats4:")
    stats5 = Cpt(StatsPluginV33, "Stats5:")
    roi1 = Cpt(ROIPlugin, "ROI1:")
    roi2 = Cpt(ROIPlugin, "ROI2:")
    roi3 = Cpt(ROIPlugin, "ROI3:")
    roi4 = Cpt(ROIPlugin, "ROI4:")
    proc1 = Cpt(ProcessPlugin, "Proc1:")

    tiff = Cpt(
        TIFFPluginWithFileStore,
        suffix="TIFF1:",
        write_path_template="/nsls2/chx/%Y/%m/%d/",
        root="/nsls2/chx",
    )


class Pilatus800V33(PilatusV33):

    tiff = Cpt(
        TIFFPluginWithFileStore,
        suffix="TIFF1:",
        read_path_template="",
        write_path_template="",
        root="",
    )