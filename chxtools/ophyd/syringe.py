from ophyd import Component as Cpt
from ophyd import (
    Device,
    EpicsSignal,
    EpicsSignalRO,
)
from ophyd import DynamicDeviceComponent as DDC


class Syringe_Pump(Device):
    "Syringe pump controller"
    pump_series = lambda x, sig, r= [1,2]: {f'p{i}': (sig, f'{{Pmp:{i}}}{x}', {}) for i in r}

    vol_sp = DDC(pump_series('Val:Vol-SP',EpicsSignal))
    vol_rb = DDC(pump_series('Val:Vol-RB',EpicsSignalRO))
    
    rate_sp = DDC(pump_series('Val:Rate-SP',EpicsSignal))
    rate_rb = DDC(pump_series('Val:Rate-RB',EpicsSignalRO))
    
    dia_sp = DDC(pump_series('Val:Dia-SP',EpicsSignal))
    dia_rb = DDC(pump_series('Val:Dia-RB',EpicsSignalRO))
    
    dir_sp   =   DDC({f'p{i}': (EpicsSignal, '{Pmp:%s}Val:Dir-Sel'%i, {}) for i in range(1, 3)})
    dir_rb   =   DDC({f'p{i}': (EpicsSignal, '{Pmp:%s}Val:Dir-Sts'%i, {}) for i in range(1, 3)})
    
    vol_unit_sp = DDC({f'p{i}': (EpicsSignal, '{Pmp:%s}Unit:Vol-Sel'%i, {}) for i in range(1, 3)})
    rate_unit_sp = DDC({f'p{i}': (EpicsSignal, '{Pmp:%s}Unit:Rate-Sel'%i, {}) for i in range(1, 3)})      
    vol_unit_rb = DDC({f'p{i}': (EpicsSignalRO, '{Pmp:%s}Unit:Vol-Sel'%i, {}) for i in range(1, 3)})
    rate_unit_rb = DDC({f'p{i}': (EpicsSignalRO, '{Pmp:%s}Unit:Rate-Sel'%i, {}) for i in range(1, 3)})     
      
    
    Run   =    DDC({f'p{i}': (EpicsSignal, '{Pmp:%s}Cmd:Run-Cmd'%i, {}) for i in [1,2,'All'] } )  
    Stop   =   DDC({f'p{i}': (EpicsSignal, '{Pmp:%s}Cmd:Stop-Cmd'%i, {}) for i in [1,2,'All'] } ) 
    Purge  =   DDC({f'p{i}': (EpicsSignal, '{Pmp:%s}Cmd:Purge-Cmd'%i, {}) for i in [1,2,'All'] } ) 
          
    disI_vol_rb = DDC({f'p{i}': (EpicsSignal, '{Pmp:%s}Val:InfDisp-I'%i, {}) for i in [1,2] } ) 
    disI_unit_rb = DDC({f'p{i}': (EpicsSignal, '{Pmp:%s}Val:InfDisp-I.EGU'%i, {}) for i in [1,2] } )
    
    disW_vol_rb = DDC({f'p{i}': (EpicsSignal, '{Pmp:%s}Val:WdlDisp-I'%i, {}) for i in [1,2] } ) 
    disW_unit_rb = DDC({f'p{i}': (EpicsSignal, '{Pmp:%s}Val:WdlDisp-I.EGU'%i, {}) for i in [1,2] } )    
    
    Clr = DDC({f'p{i}': (EpicsSignal, '{Pmp:%s}Cmd:Clr-Cmd'%i, {}) for i in [1,2] } ) 
         
   
    
    p='All'
    run_all = Cpt(EpicsSignal, '{Pmp:%s}Cmd:Run-Cmd'%p)
    purge_all = Cpt(EpicsSignal, '{Pmp:%s}Cmd:Purge-Cmd'%p)
    stop_all = Cpt(EpicsSignal, '{Pmp:%s}Cmd:Stop-Cmd'%p)    
     
    def get_vol( self, pump =1 ):
        if pump==1:
            return self.vol_rb.p1.value
        else:
            return self.vol_rb.p2.value
            
    def get_vol_unit( self, pump =1 ):
        if pump==1:
            return self.vol_unit_rb.p1.value
        else:
            return self.vol_unit_rb.p2.value 
               
    def set_vol_unit( self, pump =1 ):
        if pump==1:
            return self.vol_unit_sp.p1.value
        else:
            return self.vol_unit_sp.p2.value  
            
    def get_rate_unit( self, pump =1 ):
        if pump==1:
            return self.rate_unit_rb.p1.value
        else:
            return self.rate_unit_rb.p2.value 
               
    def set_rate_unit( self, pump =1 ):
        if pump==1:
            return self.rate_unit_sp.p1.value
        else:
            return self.rate_unit_sp.p2.value    
                                           
            
    def set_vol( self, val, pump =1 ):
        if pump==1:
            return self.vol_sp.p1.set( val) 
        else:
            return self.vol_sp.p2.set( val) 
            
    def get_rate( self, pump =1 ):
        if pump==1:
            return self.rate_rb.p1.value
        else:
            return self.rate_rb.p2.value
            
    def set_rate( self, val, pump =1 ):
        if pump==1:
            return self.rate_sp.p1.set( val) 
        else:
            return self.rate_sp.p2.set( val) 
            
    def get_dia( self, pump =1 ):
        if pump==1:
            return self.dia_rb.p1.value
        else:
            return self.dia_rb.p2.value
            
    def set_dia( self, val, pump =1 ):
        if pump==1:
            return self.dia_sp.p1.set( val) 
        else:
            return self.dia_sp.p2.set( val)
            
    def get_dir( self, pump =1 ):
        if pump==1:
            return self.dir_rb.p1.value
        else:
            return self.dir_rb.p2.value
            
    def set_dir( self, val, pump =1 ):
        if pump==1:
            return self.dir_sp.p1.set( val) 
        else:
            return self.dir_sp.p2.set( val)  
            
    def run( self, pump=1 ):
         if pump==1:
            return self.Run.p1.set(1) 
         elif pump==2:
            return self.Run.p2.set(1)                 
         else:
            return self.Run.pAll.set(1)       # edit from sp to self TODO: check if safe fix     
            
    def stop( self, pump=1 ):
         if pump==1:
            return self.Stop.p1.set(1) 
         elif pump==2:
            return self.Stop.p2.set(1)                 
         else:
            return self.Stop.pAll.set(1)      
            
    def get_disvol( self,  pump=1, Dir=None, ): 
        '''dir: 0 for infusion, 1 for withdraw '''
        if Dir is None:
            Dir = self.get_dir(   pump  ) 
        if Dir ==0:
            if pump==1:
                return self.disI_vol_rb.p1.value
            else:
                return self.disI_vol_rb.p2.value  
        else:
            if pump==1:
                return self.disW_vol_rb.p1.value
            else:
                return self.disW_vol_rb.p2.value   
           
    def clr( self,  pump=1, Dir=None, ): 
        '''dir: 0 for infusion, 1 for withdraw '''
        if Dir is None:
            Dir = self.get_dir(   pump  ) 
        if pump==1:
            return self.Clr.p1.set( Dir )   
        else:
            return self.Clr.p2.set( Dir )   
                          
    
    @property
    def hints(self):
        fields = []
        for name in self.component_names:
            motor = getattr(self, name)
            fields.extend(motor.hints['fields'])
        return {'fields': fields}                   
