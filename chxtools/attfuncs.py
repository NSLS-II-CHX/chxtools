# -*- coding: utf-8 -*-
"""
Functions for CHX attenuator system:
physical setup of attenuator system -> att_setup
calculate optimal transmission -> calc_T
set optimal transmission -> set_T
by LW 12/08/2015
"""
import numpy as np
import chxtools.xfuncs as xf

from epics import *

##PVs
# XF:11IDB-BI{Attn:03}Cmd:In-Cmd
# XF:11IDB-BI{Attn:03}Cmd:Out-Cmd
#  XF:11IDB-BI{Attn:03}Pos-Sts
#  caget ('XF:11IDB-BI{Attn:03}Pos-Sts')  #to get the status
#   caput ('XF:11IDB-BI{Attn:03}Cmd:Out-Cmd',1)  #take out 
# caput ('XF:11IDB-BI{Attn:03}Cmd:In-Cmd',1)  #put in



####
#Wayne help create two pv
# for set value      XF:11IDB-BI{Attn}Val:Trans-SP
# for display value  XF:11IDB-BI{Attn}Val:Trans-I


#caput ('XF:11IDB-BI{Attn:03}Cmd:Out-Cmd',1)


def att_setup():
    """
    defining physical configuration of CHX attenuator system
    call: att_setup()
    returns:[[att_material list],[Cu thickness list],[Si thickness list]]
    dependencies: none
    by LW 12/10/2015
    """
    # slot9  slot8  slot7  slot6  slot5  slot4  slot3 slot2 slot1 <------ beam direction
    #  Si     Si     Si     Si     Si     Si     Si    Cu    Cu           absorber material
    # 25um  50um 100um  200um  400um  800um  1600um  120um  220um         absorber thickness
    Cu_th=[220,120];Si_th=[1600,800,400,200,100,50,25]
    att_mat=['Cu','Cu','Si','Si','Si','Si','Si','Si','Si']
    return [att_mat,Cu_th,Si_th]


def calc_T(T,E='auto', foil_mode='Si'):  
    """
    calc_T(T,E='auto',foil_mode='Si'): funtion to calculate CHX attenuator settings for a commanded transmission
    required arument:
    T: commanded tranmission (e.g. 1e-3)
    optional arguments:
    E='auto': get beamline energy from PV (currently: DCM) | E=7894: manual overwrite of energy parameter
    foil_mode='Si'|'mix': 'Si': calculating optimal configuration using only Si wafers (for coherence applications) | 'mix': using Si wafers and metal foils
    function returns dictionary of the form: {'tot_index':a,'Si_index':Si_ind,'Cu_index':Cu_ind}
    'tot_index': 'binary' np.array, where '1' means attenuator 'in', '0' means attenuator 'out'
    'Si_index': same as above, just for the Si wafer slots
    'Cu_index': same as above, just for the metal foil slots
    dependencies: imports numpy and xfuncs
    calls att_setup() to get physical configuration of attenuator system
    by LW 12/09/2015
    """      
    T=np.array(T)    
    if E is 'auto':
       #E=8000   # temporary: don't have channel access -> set E to 8000eV   
       E=caget('XF:11IDA-OP{Mono:DCM-Ax:Energy}Mtr.RBV')     ### get energy automatically with channel access
       print ('getting energy from global PV: E=',E,'eV (currently not implemented in test version (no channel access) -> 8000eV default)')   # future: add PV name for house keeping
       if E> 30000 or E< 2000:
           raise attfuncs_Exception("error: Input argument E has to be 2000<E<30000 [eV]")
    else:
        if E<= 30000 and E>=2000:
                E=np.array(E)
                print ('manual input: E= ',E,'eV')
        else:
            raise attfuncs_Exception("error: could not convert energy argument. Input argument E has to be 2000<E<30000 [eV] or 'auto'")
   
    if foil_mode is 'Si':
        print ('selected foil mode is "Si" -> only Si wafers will be used for attenuation')
    elif foil_mode is 'mix':
        print ('selected foil mode is "mix" -> both Si wafers and metal foils will be used for attenuation')
    else:
        raise attfuncs_Exception("error: foil_mode has to be either 'Si' or 'mix'")      
            
    # attenuator setup
    att_conf=att_setup()
    Si_th=np.array(att_conf[2])       
    Cu_th=np.array(att_conf[1])
    abs_mat=att_conf[0]
    
    if foil_mode is 'Si':
        Cu_th=Cu_th*0
    
    # calculate all available absorption possibilities
    abs_th=np.append(Cu_th,Si_th)
    f=str('{:0'+str(len(abs_th))+'b}')
    a=np.zeros(len(abs_th))
    T_tot=np.zeros(2**len(abs_th))
    sT=np.zeros(len(abs_th))
    for m in range(0, len(abs_th)):
            sT[m]=xf.get_T(abs_mat[m],E/1000,abs_th[m])
    for l in range(0, 2**len(abs_th)):
        k=f.format(l)
        a=np.zeros(len(abs_th))
        for h in range(0, len(abs_th)):
            a[h]=int(k[h])
        x=sT*a;x[x==0]=1
        T_tot[l]=np.product(x)
    
    # determine best attenuator configuration    
    diff=np.abs(T_tot-T);       
    best_T=T_tot[np.argmin(diff)]
    k=f.format(np.argmin(diff))
    for m in range(0, len(abs_th)):
        a[m]=int(k[m])
    Si_ind=a[len(Cu_th):len(a)]
    Cu_ind=a[0:len(Cu_th)]
    
    print ('*'*40)
    print (T, best_T)
    print ('*'*40)
    
    
    caput ('XF:11IDB-BI{Attn}Val:Trans-SP',  T)  #set point
    caput ('XF:11IDB-BI{Attn}Val:Trans-I', best_T)  #the best available
    
    # some output and return value
    print ('requested transmission: ','{:.2e}'.format(float(T)),'    closest match: ','{:.2e}'.format(best_T),'    difference: ','{:.2e}'.format(best_T-T))
    print ('Si wafer configuration: ',Si_th*Si_ind,'    Cu foil configuration: ',Cu_th*Cu_ind)
    
    return {'tot_index':a,'Si_index':Si_ind,'Cu_index':Cu_ind}    

def set_T(T,E='auto',foil_mode='Si'):
    """
    call: set_T(T,E='auto',foil_mode='Si')
    function calls calc_T(T,E='auto',foil_mode='Si')
    required arument:
    T: commanded tranmission (e.g. 1e-3)
    optional arguments:
    E='auto': get beamline energy from PV (currently: DCM) | E=7894: manual overwrite of energy parameter
    foil_mode='Si'|'mix': 'Si': calculating optimal configuration using only Si wafers (for coherence applications) | 'mix': using Si wafers and metal foils
    function sets attenuator configuration via channel access according to calculation by calc_T()
    function reads back configuration and cross-checks that attenuator configuration was successfully set    
    by LW 12/10/2015
    """
    setting=calc_T(T,E,foil_mode)
    
    # setting the attenuators according to calculation
    for m in range(0,len(setting['tot_index'])):
        if int(setting['tot_index'][m]):
            caput ('XF:11IDB-BI{Attn:%02d}Cmd:In-Cmd'%(m+1),1)  #put in
            print ('the %s--attenuator is put in.'%(m+1))
        else:
            caput ('XF:11IDB-BI{Attn:%02d}Cmd:Out-Cmd'%(m+1),1)  #take out               
            print ('the %s--attenuator is taken out.'%(m+1))
        


        
#just printing, adjust PV name, remove 'print' and remove quotes for activating channel access command
#    # cross check that attenuators were set correctly                  #### Uncomment when PVs are available
#        check_ind=np.zeros(len(setting['tot_index']))
#    for l in range(0,len(setting['tot_index'])):
#        check_ind[l]=caget("'XF:11IDB-ATT{slot:",str(int(m+1)),"}RBV-VALUE'")
##    if np.array_equal(np.array(setting['tot_index']),np.array(check_ind))==True:
#        print 'check: attenuators successfully set!'
#    else: raise attfuncs_Exception("ERROR: something went wrong when setting the attenuators")
##    NEED to set foil_mode, E_mode and E in CSS screen -> add when screen and PVs are available    
#    
def get_T(E='auto'):
    """
    get transmission from current CHX attenuator settings
    call: get_T(E='auto')
    E='auto': read beamline energy from PV (currently DCM) | e.g. E=7645: manual input of energy
    reports Si wafer and metal foil configuration as well as transmission
    return transmission
    dependencies: calls att_conf()
    imports xfuncs, numpy
    by LW 12/10/2015
    """
    # attenuator setup
    att_conf=att_setup()
    Si_th=np.array(att_conf[2])       
    Cu_th=np.array(att_conf[1])
    abs_mat=att_conf[0]
    if E is 'auto':
       #E=8000   # temporary: don't have channel access -> set E to 8000eV   
       E=caget('XF:11IDA-OP{Mono:DCM-Ax:Energy}Mtr.RBV')     ### get energy automatically with channel access
       print ('getting energy from global PV: E=',E,'eV (currently not implemented in test version (no channel access) -> 8000eV default)')   # future: add PV name for house keeping
       if E> 30000 or E< 2000:
           raise attfuncs_Exception("error: Input argument E has to be 2000<E<30000 [eV]")
    else:
        if E<= 30000 and E>=2000:
                E=np.array(E)
                print ('manual input: E= ',E,'eV')
        else:
            raise attfuncs_Exception("error: could not convert energy argument. Input argument E has to be 2000<E<30000 [eV] or 'auto'")
    abs_th=np.append(Cu_th,Si_th)
    read_conf=np.zeros(len(abs_th))
    for m in range(0,len(abs_th)):              ### uncomment when chennel access and PVs available
        read_conf[m] =  caget ('XF:11IDB-BI{Attn:%02d}Pos-Sts'%(m+1))  
         
        #read_conf[m]=caget("'XF:11IDB-ATT{slot:",str(int(m+1)),"}RBV-VALUE'")
#   read_conf=np.array([0,1,0,1,0,1,0,0,1])   ### just for testing!! comment/delete!
    
    
    
    sT=np.zeros(len(abs_th))
    for m in range(0, len(abs_th)):
            sT[m]=xf.get_T(abs_mat[m],E/1000.,abs_th[m])
    x=sT*read_conf;x[x==0]=1
    T_tot=np.product(x)
    Si_ind=read_conf[len(Cu_th):len(read_conf)]
    Cu_ind=read_conf[0:len(Cu_th)]
    print ('Found:')
    print ('Si wafer configuration: ',Si_th*Si_ind,'    Cu foil configuration: ',Cu_th*Cu_ind)
    print ('current Transmission at ',E, 'eV: T= ','{:.2e}'.format(T_tot))
    return T_tot
    

        
class attfuncs_Exception(Exception):
    pass
    """
    by LW 12/10/2015
    class to raise att_funcs specific exceptions
    """