from nslsii.detectors.xspress3 import Xspress3Channel, Xspress3Detector, XspressTrigger
from nslsii.detectors.zebra import (
    EpicsSignalWithRBV,
    ZebraAddresses,
    ZebraFrontOutput3,
    ZebraFrontOutput4,
    ZebraFrontOutput12,
    ZebraGate,
    ZebraPulse,
    ZebraRearOutput,
)
from ophyd import Cpt, Device, EpicsSignal, EpicsSignalRO
from ophyd import FormattedComponent as FC
from ophyd.areadetector.plugins import PluginBase


class ZebraINP(Device):
    use = FC(EpicsSignal, 
             '{self.prefix}_ENA:B{self._bindex}')
    source_addr = FC(EpicsSignalWithRBV, 
                     '{self.prefix}_INP{self.index}')
    source_str = FC(EpicsSignalRO, 
                    '{self.prefix}_INP{self.index}:STR', 
                    string=True)
    source_status = FC(EpicsSignalRO, 
                    '{self.prefix}_INP{self.index}:STA')
    invert = FC(EpicsSignal, 
                '{self.prefix}_INV:B{self._bindex}')

    def __init__(self, prefix, *, index, read_attrs=None, configuration_attrs=None,
                 **kwargs):
        if read_attrs is None:
            read_attrs = []
        if configuration_attrs is None:
            configuration_attrs = ['use', 'source_addr', 'source_str', 'invert']
        self.index = index
        self._bindex = index - 1
        super().__init__(prefix, read_attrs=read_attrs, configuration_attrs=configuration_attrs,
                         **kwargs)
    

class ZebraLogic(Device):

    inp1 = Cpt(ZebraINP, '', index=1)
    inp2 = Cpt(ZebraINP, '', index=2)
    inp3 = Cpt(ZebraINP, '', index=3)
    inp4 = Cpt(ZebraINP, '', index=4)
    out = Cpt(EpicsSignalRO, '_OUT')

    def __init__(self, *args, read_attrs=None, configuration_attrs=None,
                 **kwargs):
        if read_attrs is None:
            read_attrs = ['out']
        if configuration_attrs is None:
            configuration_attrs = [f'inp{j}' for j in range(1, 5)]

        super().__init__(*args, read_attrs=read_attrs, 
                          configuration_attrs=configuration_attrs, **kwargs)
        

class CHXXspress3Detector(XspressTrigger, Xspress3Detector):
    roi_data = Cpt(PluginBase, 'ROIDATA:')
    channel1 = Cpt(Xspress3Channel, 
                   'C1_', channel_num=1, 
                   read_attrs=['rois'])

class GateConfig(Device):
    trig_source = Cpt(EpicsSignalWithRBV, 'PC_GATE_SEL', string=True)
    start = Cpt(EpicsSignal, 'PC_GATE_START')
    width = Cpt(EpicsSignal, 'PC_GATE_WID')
    step = Cpt(EpicsSignal, 'PC_GATE_STEP')
    ngate = Cpt(EpicsSignal, 'PC_GATE_NGATE')
    status = Cpt(EpicsSignal, 'PC_GATE_OUT')
    

class PulseConfig(Device):
    trig_source = Cpt(EpicsSignalWithRBV, 'PC_PULSE_SEL', string=True)
    start = Cpt(EpicsSignal, 'PC_PULSE_START')
    width = Cpt(EpicsSignal, 'PC_PULSE_WID')
    step = Cpt(EpicsSignal, 'PC_PULSE_STEP')
    delay = Cpt(EpicsSignal, 'PC_PULSE_DLY')
    nmax = Cpt(EpicsSignal, 'PC_PULSE_MAX')
    status = Cpt(EpicsSignal, 'PC_PULSE_OUT')

class Zebra(Device):
    soft_input1 = Cpt(EpicsSignal, 'SOFT_IN:B0')
    soft_input2 = Cpt(EpicsSignal, 'SOFT_IN:B1')
    soft_input3 = Cpt(EpicsSignal, 'SOFT_IN:B2')
    soft_input4 = Cpt(EpicsSignal, 'SOFT_IN:B3')

    pulse1 = Cpt(ZebraPulse, 'PULSE1_', index=1)
    pulse2 = Cpt(ZebraPulse, 'PULSE2_', index=2)
    pulse3 = Cpt(ZebraPulse, 'PULSE3_', index=3)
    pulse4 = Cpt(ZebraPulse, 'PULSE4_', index=4)

    output1 = Cpt(ZebraFrontOutput12, 'OUT1_', index=1)
    output2 = Cpt(ZebraFrontOutput12, 'OUT2_', index=2)
    output3 = Cpt(ZebraFrontOutput3, 'OUT3_', index=3)
    output4 = Cpt(ZebraFrontOutput4, 'OUT4_', index=4)

    output5 = Cpt(ZebraRearOutput, 'OUT5_', index=5)
    output6 = Cpt(ZebraRearOutput, 'OUT6_', index=6)
    output7 = Cpt(ZebraRearOutput, 'OUT7_', index=7)
    output8 = Cpt(ZebraRearOutput, 'OUT8_', index=8)

    gate1 = Cpt(ZebraGate, 'GATE1_', index=1)
    gate2 = Cpt(ZebraGate, 'GATE2_', index=2)
    gate3 = Cpt(ZebraGate, 'GATE3_', index=3)
    gate4 = Cpt(ZebraGate, 'GATE4_', index=4)

    or1 = Cpt(ZebraLogic, 'OR1')
    or2 = Cpt(ZebraLogic, 'OR2')
    or3 = Cpt(ZebraLogic, 'OR3')
    or4 = Cpt(ZebraLogic, 'OR4')

    and1 = Cpt(ZebraLogic, 'AND1')
    and2 = Cpt(ZebraLogic, 'AND2')
    and3 = Cpt(ZebraLogic, 'AND3')
    and4 = Cpt(ZebraLogic, 'AND4')

    addresses = ZebraAddresses

    arm = Cpt(EpicsSignal, 'PC_ARM')
    arm_status = Cpt(EpicsSignal, 'PC_ARM_OUT')
    arm_trig = Cpt(EpicsSignalWithRBV, 'PC_ARM_SEL', string=True)
    
    gate_config = Cpt(GateConfig, '', read_attrs=[], configuration_attrs=['trig_source', 'start', 'width', 'step', 'ngate'])
    pulse_config = Cpt(PulseConfig, '', read_attrs=[], configuration_attrs=['trig_source', 'start', 'width', 'step', 'delay', 'nmax'])
    


class CHXZebra(Zebra):
    ...



class XspressZebra(Device):
    zebra = Cpt(CHXZebra, 'XF:11IDB-ES{Zebra}:', add_prefix=(),
                configuration_attrs=[f'{n}{j}' for n in ['pulse', 'gate', 'or', 'and']  
                                     for j in range(1,5)] + [f'output{i}' for i in range(1, 9)] + ['gate_config', 'pulse_config', 'arm_trig'],
                read_attrs=[])
    #xs = Cpt(CHXXspress3Detector, 'XSPRESS3-EXAMPLE:', add_prefix=())
    xs = Cpt(CHXXspress3Detector, 'XF:11IDB-ES{Xsp:1}:', add_prefix=())
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        

    def trigger(self):
        st = self.xs.trigger()
        # TODO make sure Xspress is actually ready before 
        # arming zebra
        self.zebra.arm.set(1)
        return st
    
    @property
    def hints(self):
        return {'fields': [  self.xs.channel1.rois.roi01.value.name ] } #,  self.xs.channel1.rois.roi01.value_sum.name,      ]}    
