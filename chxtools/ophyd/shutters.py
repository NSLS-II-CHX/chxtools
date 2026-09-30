from ophyd import Component as Cpt
from ophyd import Device, EpicsSignal


class FourPVShutter(Device):
    """
    Ophyd structure object for 4 component shutters (2 button mode).
    
    used for fron-end (fe_sh) and front of enclosure (foe_sh) shutters
    """
    open_command = Cpt(EpicsSignal, 'Cmd:Opn-Cmd')
    open_status = Cpt(EpicsSignal, 'Pos-Sts')
    close_command = Cpt(EpicsSignal, 'Cmd:Cls-Cmd')
    close_status = Cpt(EpicsSignal, 'Pos-Sts')

    def open(self):
        self.open_command.put(1)

    def close(self):
        self.close_command.put(1)

    @property
    def is_open(self):
        self._check_sanity()
        return self.open_status.get()

    def _check_sanity(self):
        consistency = self.open_status.get() ^ self.close_status.get()
        assert consistency, "Shutter status is not self-consistent"


class TwoPVShutter(EpicsSignal):
    """
    Ophyd structure object for 2 component shutters (switch mode).

    used for fast shutter (fast_sh) and pre shutter (pre_sh)
    """
    def open(self):
        self.put(1)
    
    def close(self):
        self.put(0)