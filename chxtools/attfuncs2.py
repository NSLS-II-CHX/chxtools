# -*- coding: utf-8 -*-
"""
Created on Mon Sep 25 15:05:26 2017

@author: lwiegart
"""

"""
Functions for CHX attenuators in monitor chamber:
physical setup of attenuator system -> att_setup
calculate optimal transmission -> calc_T
set optimal transmission -> set_T
by LW 09/25/2017
"""
import numpy as np
import xfuncs as xf
from epics import *

def att_setup():
    """
    defining physical configuration of CHX attenuators in monitor chamber
    call: att_setup()
    returns:[[position list],[Si thickness list]]
    position of attenuator is relative to center of 'empty' slot 
    dependencies: none
    by LW 09/25/2017
    """
    att_pos=[2.75,4.75,6.75,8.75,10.75,12.75,14.75,16.75,18.75]
    Si_th=[200,400,600,800,1000,1200,1400,1600,1800]
    return [att_pos,Si_th]
    
def calc_T(T,E='auto'):  
    """
    calc_T(T,E='auto'): funtion to calculate attenuator position for a commanded transmission (attenuator in monitor chamber)
    required arument:
    T: commanded tranmission (e.g. 1e-3)
    optional arguments:
    E='auto': get beamline energy from PV (currently: DCM) | E=7894: manual overwrite of energy parameter
    function returns attenuator position and transmission value
    calls att_setup() to get physical configuration of attenuators in monitor chamber
    by LW 09/25/2017
    """      
    T=np.array(T)    
    if E is 'auto':
       E=caget('XF:11IDA-OP{Mono:DCM-Ax:Energy}Mtr.RBV')     ### get energy automatically with channel access
       print('getting energy from global PV: E='+str(E)+'eV')
       if E> 30000 or E< 2000:
           raise attfuncs_Exception("error: Input argument E has to be 2000<E<30000 [eV]")
    else:
        if E<= 30000 and E>=2000:
                E=np.array(E)
        else:
            raise attfuncs_Exception("error: could not convert energy argument. Input argument E has to be 2000<E<30000 [eV] or 'auto'")           
    # attenuator setup
    att_conf=att_setup()
    Si_th=np.array(att_conf[1])       
    att_pos=att_conf[0]
    ava_T=xf.get_T('Si',E/1E3,Si_th)
    ratio_T=ava_T/T
    for i in range(len(ratio_T)):
        if ratio_T[i] <1:
            ratio_T[i]=1/ratio_T[i]
    ind=np.argmin(abs(ratio_T-1))
    bestT=ava_T[ind]
    th_Si=Si_th[ind]
    pos_att=att_pos[ind]
    if T==1:
        th_Si=0;pos_att=0;bestT=1
    elif T == 0:
        th_Si=max(Si_th);pos_att=max(np.array(att_pos));bestT=xf.get_T('Si',E/1E3,th_Si)  
    print('calculations for Si absorbers in monitor chamber:')
    print('requested transmission: '+str(T)+' at '+str(E/1E3)+'keV' )
    print('best match: '+str(bestT)+' using '+str(th_Si)+'um of Si -> T_request/T_available= '+str(T/bestT))
    return [pos_att,bestT]

def get_T(E='auto'):
    """
    call: get_T(E='auto')
    function calls att_setup()
    required arguments: none
    optional argument: E='auto': get beamline energy from PV (currently: DCM) | E=7894: manual overwrite of energy parameter
    function reads back configuration of Si absorbers in monitor chamber and returns transmission value
    pre-requisit: center position of 'empty' slot is correct
    by LW 09/25/2017
    """
    if E is 'auto':
       E=caget('XF:11IDA-OP{Mono:DCM-Ax:Energy}Mtr.RBV')     ### get energy automatically with channel access
       print('getting energy from global PV: E='+str(E)+'eV')
       if E> 30000 or E< 2000:
           raise attfuncs_Exception("error: Input argument E has to be 2000<E<30000 [eV]")
    else:
        if E<= 30000 and E>=2000:
                E=np.array(E)
        else:
            raise attfuncs_Exception("error: could not convert energy argument. Input argument E has to be 2000<E<30000 [eV] or 'auto'")   
    tol=0.5 #tolerance for accepting absorber position
    tol2=.8 #tolerance for acceptiong empty slot position
    empty_pos=-18.5                               # this needs to become a PV!!!
    current_pos=caget('XF:11IDB-OP{Mon:Foil-Ax:X}Mtr.RBV')
    # attenuator setup
    att_conf=att_setup()
    Si_th=np.array(att_conf[1])     
    att_pos=np.array(att_conf[0])
    if current_pos >=empty_pos+min(att_pos)-tol and current_pos <= empty_pos+max(att_pos)+tol:
        ind=np.argmin(abs(empty_pos+att_pos-current_pos))
        if abs(empty_pos-current_pos+att_pos[ind]) <= tol:
            currentT=xf.get_T('Si',E/1E3,Si_th[ind])
            print('current transmission at '+str(E/1E3)+'keV: '+str(currentT)+'  ('+str(Si_th[ind])+'um Si)')
    elif abs(current_pos-empty_pos) <= tol2:
        currentT=1
        print('empty slot is in the beam, transmission T=1')
    else:
        print('Neither emppty slot nor any Si absorber is aligned with the beam...transmission unknown')
        currentT=float('nan')
    return currentT
       
def set_T(T):
    """
    call: set_T(T)
    function calls calc_T(T)
    required arument:
    T: commanded tranmission (e.g. 1e-3)
    optional arguments:
    E='auto': get beamline energy from PV (currently: DCM) | E=7894: manual overwrite of energy parameter
    function moves Si absorbers in monitor chamber according to calculation by calc_T()
    function reads back configuration and cross-checks that absorber positions where successfully reached    
    by LW 09/25/2017
    """
    E=caget('XF:11IDA-OP{Mono:DCM-Ax:Energy}Mtr.RBV')     ### get energy automatically with channel access    
    tol=0.5
    empty_pos=-18.5                               # this needs to become a PV!!!
    [pos_att,bestT]=calc_T(T,E)
    target_pos=empty_pos+pos_att
    print('moving foil_x to target_pos')
    caput('XF:11IDB-OP{Mon:Foil-Ax:X}Mtr.VAL',target_pos)
    time.sleep(15)
    if abs(caget('XF:11IDB-OP{Mon:Foil-Ax:X}Mtr.RBV')-target_pos) <= tol:
        print('absorbers set successfully!')
    else: raise attfuncs_Exception("something went wrong...position for requested absorber not reached...")
    
def set_abs(abs_num=3):
    """
    call: set_abs(abs_num=3)
    function calls att_setup()
    required arument: abs_num -> number of Si absorber
    0: empry slot, 1: 200um, 2: 400um, ...9: 1800um
    function moves Si absorbers in monitor chamber to requested aborber number
    function returns current transmission value
    by LW 09/25/2017
    """
    tol=0.5
    empty_pos=-18.5
    # attenuator setup
    att_conf=att_setup()
    Si_th=np.array(att_conf[1])     
    att_pos=np.array(att_conf[0])
    if abs_num >=0 and abs_num <=9 and isinstance(abs_num,int):
        print('moving absorber #'+str(abs_num)+' in monitor chamber into the beam')
        if abs_num == 0:
            target_pos=empty_pos
            E=caget('XF:11IDA-OP{Mono:DCM-Ax:Energy}Mtr.RBV')
            curr_T=1
            curr_Sith=0
        else:
            if abs_num >0:
                target_pos=empty_pos+att_pos[abs_num-1]
                E=caget('XF:11IDA-OP{Mono:DCM-Ax:Energy}Mtr.RBV')     ### get energy automatically with channel access
                curr_T=xf.get_T('Si',E/1E3,Si_th[abs_num-1])
                curr_Sith=Si_th[abs_num-1]
        caput('XF:11IDB-OP{Mon:Foil-Ax:X}Mtr.VAL',target_pos)
        time.sleep(15)
        if abs(caget('XF:11IDB-OP{Mon:Foil-Ax:X}Mtr.VAL')-target_pos) <= tol:
            print('absorbers in monitor chamber set successfully!')
            print('transmission is '+str(curr_T)+'  with '+str(curr_Sith)+'um Si at '+str(E/1E3)+'keV')
        else: raise attfuncs_Exception("something went wrong...position for requested absorber not reached...")  
    else: raise attfuncs_Exception("error: absorber number needs to be an integer between 0 [empty slot] and 9 [9x200um Si]") 
        


class attfuncs_Exception(Exception):
    pass
    """
    by LW 12/10/2015
    class to raise att_funcs specific exceptions
    """
