from ophyd import Component as Cpt
from ophyd import (
    Device,
    EpicsMotor,
    EpicsSignal,
    EpicsSignalRO,
    FormattedComponent,
    PVPositioner,
    PVPositionerPC,
)
from ophyd import DynamicDeviceComponent as DDC


class MotorCenterAndGap(Device):
    "Center and gap using Epics Motor records"
    xc = Cpt(EpicsMotor, '-Ax:XCtr}Mtr')
    yc = Cpt(EpicsMotor, '-Ax:YCtr}Mtr')
    xg = Cpt(EpicsMotor, '-Ax:XGap}Mtr')
    yg = Cpt(EpicsMotor, '-Ax:YGap}Mtr')

    @property
    def hints(self):
        fields = []
        for name in self.component_names:
            motor = getattr(self, name)
            fields.extend(motor.hints['fields'])
        return {'fields': fields}


class VirtualGap(PVPositioner):
    readback = Cpt(EpicsSignalRO, 't2.C')
    setpoint = Cpt(EpicsSignal, 'size')
    done = Cpt(EpicsSignalRO, 'DMOV')
    done_value = 1


class VirtualCenter(PVPositioner):
    readback = Cpt(EpicsSignalRO, 't2.D')
    setpoint = Cpt(EpicsSignal, 'center')
    done = Cpt(EpicsSignalRO, 'DMOV')
    done_value = 1


class VirtualMotorCenterAndGap(Device):
    "Center and gap with virtual motors"
    xc = Cpt(VirtualCenter, '-Ax:X}')
    yc = Cpt(VirtualCenter, '-Ax:Y}')
    xg = Cpt(VirtualGap, '-Ax:X}')
    yg = Cpt(VirtualGap, '-Ax:Y}')
    @property
    def hints(self):
        fields = []
        for name in self.component_names:
            motor = getattr(self, name)
            fields.extend(motor.hints['fields'])
        return {'fields': fields}    


class Blades(Device):
    top = Cpt(EpicsMotor, '-Ax:T}Mtr')
    bottom = Cpt(EpicsMotor, '-Ax:B}Mtr')
    outboard = Cpt(EpicsMotor, '-Ax:O}Mtr')
    inboard = Cpt(EpicsMotor, '-Ax:I}Mtr')


class MotorSlits(Blades, MotorCenterAndGap):
    "combine t b i o and xc yc xg yg"
    pass

class VirtualMotorSlits(Blades, VirtualMotorCenterAndGap):
    "combine t b i o and xc yc xg yg"
#    def __init__(self, *args, **kwargs):
#      super().__init__(*args, **kwargs)
#        self.xc.readback.name = self.name
#        self.yc.readback.name = self.name
 #       self.xg.readback.name = self.name
    pass


class XYMotor(Device):
    x = Cpt(EpicsMotor, '-Ax:X}Mtr')
    y = Cpt(EpicsMotor, '-Ax:Y}Mtr')


class XYThetaMotor(XYMotor):
    "used for GI mirror"
    th = Cpt(EpicsMotor, '-Ax:Th}Mtr')


class HorizontalDiffractionMirror(XYMotor):
    "x and y with pitch, which has different read and write PVs"
    #p = FormattedComponent(EpicsSignal, read_pv='{self.prefix}-Ax:P}}E-I', write_pv='{self.prefix}-Ax:P}}E-SP', add_prefix=('read_pv', 'write_pv', 'suffix'))
    p = FormattedComponent(EpicsSignal, read_pv='{self.prefix}-Ax:P}}Pos-I', write_pv='{self.prefix}-Ax:P}}PID-SP', add_prefix=('read_pv', 'write_pv', 'suffix'))
    # for some reason we cannot scan on E-SP. This is the actual piezo voltage (max 100) while our 'usual values' are converted to urad by some other laye rof logic in the ioc
    # the currrent SP is the input of the PID feedback loop. This requitred the feedback loop to be turned ON

class DCM(Device):
    en = Cpt(EpicsMotor, '-Ax:Energy}Mtr')
    b = Cpt(EpicsMotor, '-Ax:B}Mtr')
    r = Cpt(EpicsMotor, '-Ax:R}Mtr')
    x = Cpt(EpicsMotor, '-Ax:X}Mtr')
    fp = Cpt(EpicsMotor, '-Ax:FP}Mtr')
    p = Cpt(EpicsMotor, '-Ax:P}Mtr')

 

class SAXSBeamStop( Device):
    x = Cpt(     EpicsMotor, '-Ax:X}Mtr' )
    y1 = Cpt(     EpicsMotor, '-Ax:YFT}Mtr')
    x2 = Cpt(     EpicsMotor, '-Ax:XFB}Mtr')
    y2 = Cpt(     EpicsMotor, '-Ax:YFB}Mtr')
  
    @property
    def hints(self):
        fields = []
        for name in self.component_names:
            motor = getattr(self, name)
            fields.extend(motor.hints['fields'])
        return {'fields': fields}
 

class DMM(Device):
    # en = Cpt(EpicsMotor, '-Ax:Energy}Mtr')
    b = Cpt(EpicsMotor, '-Ax:B}Mtr')
    r = Cpt(EpicsMotor, '-Ax:R}Mtr')
    x = Cpt(EpicsMotor, '-Ax:X}Mtr')
    y = Cpt(EpicsMotor, '-Ax:Y}Mtr')
    fp = Cpt(EpicsMotor, '-Ax:FP}Mtr')


class Transfocator(Device):
    crl = DDC({'num%d' % i: (EpicsMotor, '%d-Ax:X}Mtr' % i, {})
               for i in range(1, 9)})
    x = Cpt(EpicsMotor, 'Ves-Ax:X}Mtr')
    y = Cpt(EpicsMotor, 'Ves-Ax:Y}Mtr')
    z = Cpt(EpicsMotor, 'Ves-Ax:Z}Mtr')
    ph = Cpt(EpicsMotor, 'Ves-Ax:Ph}Mtr')
    th = Cpt(EpicsMotor, 'Ves-Ax:Th}Mtr')


class Kinoform(Device):
    z = Cpt(EpicsMotor, '-Ax:ZB}Mtr')
    x = Cpt(EpicsMotor, '-Ax:XB}Mtr')
    y = Cpt(EpicsMotor, '-Ax:YB}Mtr')
    chi = Cpt(EpicsMotor, '-Ax:Ch}Mtr')
    theta = Cpt(EpicsMotor, '-Ax:Th}Mtr')
    phi = Cpt(EpicsMotor, '-Ax:Ph}Mtr')
    lx = Cpt(EpicsMotor, '-Ax:XT}Mtr')
    ly = Cpt(EpicsMotor, '-Ax:YT}Mtr')

class SmarPod_x(PVPositionerPC):
    readback = Cpt(EpicsSignalRO, '-Ax:2}Pos-I')
    setpoint = Cpt(EpicsSignal, '-Ax:2}Pos-SP')
    actuate = Cpt(EpicsSignal, '}Move-Cmd')
    actuate_value = 1

class SmarPod_y(PVPositionerPC):
    readback = Cpt(EpicsSignalRO, '-Ax:3}Pos-I')
    setpoint = Cpt(EpicsSignal, '-Ax:3}Pos-SP')
    actuate = Cpt(EpicsSignal, '}Move-Cmd')
    actuate_value = 1 

class SmarPod_z(PVPositionerPC):
    readback = Cpt(EpicsSignalRO, '-Ax:1}Pos-I')
    setpoint = Cpt(EpicsSignal, '-Ax:1}Pos-SP')
    actuate = Cpt(EpicsSignal, '}Move-Cmd')
    actuate_value = 1 

class SmarPod_rx(PVPositionerPC):
    readback = Cpt(EpicsSignalRO, '-Ax:2}Rot-I')
    setpoint = Cpt(EpicsSignal, '-Ax:2}Rot-SP')
    actuate = Cpt(EpicsSignal, '}Move-Cmd')
    actuate_value = 1

class SmarPod_ry(PVPositionerPC):
    readback = Cpt(EpicsSignalRO, '-Ax:3}Rot-I')
    setpoint = Cpt(EpicsSignal, '-Ax:3}Rot-SP')
    actuate = Cpt(EpicsSignal, '}Move-Cmd')
    actuate_value = 1

class SmarPod_rz(PVPositionerPC):
    readback = Cpt(EpicsSignalRO, '-Ax:1}Rot-I')
    setpoint = Cpt(EpicsSignal, '-Ax:1}Rot-SP')
    actuate = Cpt(EpicsSignal, '}Move-Cmd')
    actuate_value = 1

class SmarAct_single_axis(Device):
    x = Cpt(EpicsMotor, '{POS:1}Mtr')

class XBPM( Device):
   vt = Cpt( EpicsSignal, 'CtrlDAC:BLevel-SP' )

class amp_motor(Device):
    #x = EpicsMotor('XF:11IDB-OP{BS:Sam-Ax:X}Mtr')
    ampx = EpicsSignal('XF:11IDB-OP{BS:Samp-Ax:X}Mtr.VAL', name='ampx')
    ampy = EpicsSignal('XF:11IDB-OP{Stg:Samp-Ax:Phi}Mtr.VAL', name='ampy')
    ampz = EpicsSignal('XF:11IDB-OP{BS:Samp-Ax:Y}Mtr.VAL', name='ampz')
    #caput('XF:11IDB-ES{Det:Eig4M}cam1:NumImages',  fnum )
