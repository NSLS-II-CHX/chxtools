from ophyd import Component as Cpt
from ophyd import Device, EpicsSignal, EpicsSignalRO


class XBpmFeedbackEnable(Device):
    x = Cpt(EpicsSignal, 'Fdbk:AEn-SP')
    y = Cpt(EpicsSignal, 'Fdbk:BEn-SP')  

    def on(self):
        fb_x_val = self.x.set(1)
        fb_y_val = self.y.set(1) 
        print('Turning XBPM feedback on. \nYou are trusting that the control loop engaged.')
        print('\tTurn off: RE( bps.wait( xbpm.feedback_enable.off() ) )')
        return fb_x_val & fb_y_val

    def off(self):
        fb_x_val = self.x.set(0)
        fb_y_val =  self.y.set(0) 
        print('Turning XBPM feedback off.')
        print('\tTurn on: RE( bps.wait( xbpm.feedback_enable.on() ) )')
        return fb_x_val & fb_y_val


class XBpm(Device):
    x = Cpt(EpicsSignalRO, 'Pos:X-I')
    y = Cpt(EpicsSignalRO, 'Pos:Y-I')
    a = Cpt(EpicsSignalRO, 'Ampl:ACurrAvg-I')
    b = Cpt(EpicsSignalRO, 'Ampl:BCurrAvg-I')
    c = Cpt(EpicsSignalRO, 'Ampl:CCurrAvg-I')
    d = Cpt(EpicsSignalRO, 'Ampl:DCurrAvg-I')

class XBPM( Device):
   vt = Cpt( EpicsSignal, 'CtrlDAC:BLevel-SP' )

class Elm(Device):
	sum_x = Cpt(EpicsSignalRO, 'SumX:MeanValue_RBV')
	sum_y = Cpt(EpicsSignalRO, 'SumY:MeanValue_RBV')
	sum_all = Cpt(EpicsSignalRO, 'SumAll:MeanValue_RBV', kind='hinted')