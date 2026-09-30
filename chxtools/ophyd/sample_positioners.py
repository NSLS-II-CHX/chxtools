from ophyd import Component as Cpt
from ophyd import Device, EpicsMotor, PseudoPositioner, PseudoSingle


class Diffractometer(Device):
    
    Del= Cpt( EpicsMotor, '-Ax:Del}Mtr')
    gam = Cpt(EpicsMotor, '-Ax:Gam}Mtr')
    om = Cpt(EpicsMotor, '-Ax:Om}Mtr')
    phi = Cpt(EpicsMotor, '-Ax:Ph}Mtr')
    xb = Cpt(EpicsMotor, '-Ax:XB}Mtr')
    yb = Cpt(EpicsMotor, '-Ax:YB}Mtr')
    chh = Cpt(EpicsMotor, '-Ax:ChH}Mtr')
    thh = Cpt(EpicsMotor, '-Ax:ThH}Mtr')
    phh = Cpt(EpicsMotor, '-Ax:PhH}Mtr')
    xh = Cpt(EpicsMotor, '-Ax:XH}Mtr')
    yh = Cpt(EpicsMotor, '-Ax:YH2}Mtr')
    zh = Cpt(EpicsMotor, '-Ax:ZH}Mtr')
    chv = Cpt(EpicsMotor, '-Ax:ChV}Mtr')
    #thv = Cpt(EpicsMotor, '-Ax:ThV}Mtr')
    phv = Cpt(EpicsMotor, '-Ax:ThV}Mtr')  #Xiao's requested change
    xv = Cpt(EpicsMotor, '-Ax:XV}Mtr')
    yv = Cpt(EpicsMotor, '-Ax:YV}Mtr')
    zv = Cpt(EpicsMotor, '-Ax:ZV}Mtr')
    xv2 = Cpt(EpicsMotor, '-Ax:XV2}Mtr')

    @property
    def hints(self):
        fields = []
        for name in self.component_names:
            motor = getattr(self, name)
            fields.extend(motor.hints['fields'])
        return {'fields': fields}


class Goniometer(Device):
    """
    Ophyd component for full space translation goniometer (3+theta).

    Primarily used for the final stage smaract goniometer "gonzo"
    """
    # rot = Cpt(EpicsMotor, '-Ax:X}Mtr') #spare:3
    # px = Cpt(EpicsMotor, '-Ax:X}Mtr') #spare: 12
    theta = Cpt(EpicsMotor,'3-Ax:X}Mtr', name='theta')
    x = Cpt(EpicsMotor,'12-Ax:X}Mtr', name='gonx')
    y = Cpt(EpicsMotor,'13-Ax:X}Mtr', name='gony')
    z = Cpt(EpicsMotor,'14-Ax:X}Mtr', name='gonz')


class SamplePositioner(PseudoPositioner):
    '''
    Maintains an offset between a master/slave set of positioners
    such that the slave movement is the negative of the master's relative
    motion (i.e. maintains constant, negative, relative offset).
    Assumes that the user has adjusted the axes to their initial positions.

    Example
    -------------
    psamp_x = SamplePositioner(prefix='', name='psamp_x', concurrent=True)
    '''
    physical_sample_holder = Cpt(EpicsMotor, 'XF:11IDB-ES{Dif-Ax:XH}Mtr')
    beamstop_x = Cpt(EpicsMotor, 'XF:11IDB-OP{BS:Samp-Ax:X}Mtr')
    sample_holder = Cpt(PseudoSingle, limits=(0, 0))

    def forward(self, pos):
        "pos is a self.PseudoPosition"
        delta = pos - pos.sample_holder.pos
        return self.RealPosition(physical_sample_holder=pos.sample_holder,
                                 beamstop_x=self.beamstop_x.position - delta)

    def inverse(self, pos):
        "pos is self.RealPosition"
        return self.PseudoPosition(sample_holder=pos.physical_sample_holder)