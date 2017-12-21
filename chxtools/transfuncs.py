# -*- coding: utf-8 -*-
"""
Functions for CHX transfocator system:
physical setup of transfocator system -> trans_setup
calculate optimal settings for desired image position and energy -> calc_transsetup
calculate current image position -> get_ip
by LW 03/13/2016
"""
import numpy as np
from chxtools import xfuncs as xf
from scipy.optimize import fminbound
from epics import *
def trans_setup():
    """
    defining physical configuration of CHX transfocator system
    call: trans_setup()
    returns a dictionary: {'lens_material':lens_mat,'lens_number':lens_N,'lens_radius':lens_R,'trans_position':trans_pos}
    dependencies: none
    by LW 03/13/2016
    """
    # slot7  slot6  slot5  slot4  slot3 slot2 slot1   <------ beam direction (slot 8 is currently B-fiber only)
    #  Be     Be     Be     Be     Be     Be   Be            lens material
    #  1.5   1.5     0.5    0.5   0.5    0.5   0.5           lens radius [mm]
    #   1     1       5      8      4     2     1            number of lenses
    lens_R=[0.5,0.5,0.5,0.5,0.5,1.5,1.5]
    lens_mat=['Be','Be','Be','Be','Be','Be','Be']
    lens_N=[1,2,4,8,5,1,1]
    trans_pos=[35.2,35.8]
    return {'lens_material':lens_mat,'lens_number':lens_N,'lens_radius':lens_R,'trans_position':trans_pos}


def calc_transsetup(image_pos,E='auto',silent=False):
    """
    calc_transsetup(image_pos,E='auto'): funtion to calculate CHX transfocator settings
    required arument:
    image_pos: image position [m from source]
    optional arguments:
    E='auto': get beamline energy from PV (currently: DCM) | E=7894: manual overwrite of energy parameter
    silent='False': silence the output (e.g. to use function in a loop): False/True
    Note: 'virtual images' are discarded in the search
    Note: search is implemented to go through the lens stacks with large numbers of lenses first
    -> returned solution will be e.g. using lens stack #5 with 5 lenses, instead of #1 (1) and #3 (4)
    function returns dictionary of the form: {'lens_config':bin_index,'z_position':zmin}
    dependencies: imports numpy, fminbound (scipy.optimze) and xfuncs
    calls trans_setup() to get physical configuration of transfocator system
    by LW 03/14/2016
    """
    image_pos=np.array(image_pos)
    if E is 'auto':
       #E=8000   # temporary: don't have channel access -> set E to 8000eV
       E=caget('XF:11IDA-OP{Mono:DCM-Ax:Energy}Mtr.RBV')     ### get energy automatically with channel access
       print ("getting energy from global PV: E=",E,'eV')
       #(currently not implemented in test version (no channel access) -> 8000eV default)'   # future: add PV name for house keeping
       if E> 30000 or E< 2000:
           raise transfuncs_Exception("error: Input argument E has to be 2000<E<30000 [eV]")
    E=np.array(E)*1.


   # transfocator setup
    trans_conf=trans_setup()
    lens_R = np.array(trans_conf['lens_radius'])*1.
    zpos = np.array(trans_conf['trans_position'])*1.
    lens_N = np.array(trans_conf['lens_number'])
    lens_mat = np.array(trans_conf['lens_material'])


   # calculate all available focal lengths:
    f=str('{:0'+str(len(lens_R))+'b}')
    delta=np.zeros(len(lens_R))
    for w in range(0, len(lens_R)):
        delta[w]=np.real(xf.get_n(lens_mat[w],E/1000.))
    F=np.zeros(2**len(lens_R))
    zmin=np.zeros(2**len(lens_R))
    impos_diff=np.zeros(2**len(lens_R))
    q=np.zeros(2**len(lens_R))
    real_im=np.zeros(2**len(lens_R))
    for l in range(0, 2**len(lens_R)):
        k=f.format(l)
        a=np.zeros(len(lens_R))
        for h in range(0, len(lens_R)):
            a[h]=int(k[h])
        ln=lens_N*a
        F[l]=1/np.sum(2*delta*ln/(lens_R*1E-3))
        #print F[l]
        def image_func(x):
            return abs(x - image_pos + 1/(-1/x + 1/F[l]))

        zmin[l] = fminbound(image_func, zpos[0], zpos[1])
        q[l] = 1/(-1/zmin[l]+1/F[l])
        if q[l]<0:
            real_im[l]=0
            q[l]=float('NaN')
        impos_diff[l] = q[l]+zmin[l]-image_pos

    # looking for the best option
    index=np.nanargmin(abs(impos_diff))
    bi=f.format(index)
    bin_index=np.zeros(len(lens_R))
    for g in range(0, len(lens_R)):
            bin_index[g]=int(bi[g])
    F= F[index]
    zmin = zmin[index]
    q = q[index]
    impos_diff = impos_diff[index]
    if impos_diff <= 0: loc='upstream'
    if impos_diff >0: loc = 'downstream'
    #conf_lens_mat=lens_mat*bin_index
    conf_lensR=lens_R*bin_index
    conf_lensN=lens_N*bin_index
    if silent == False:
        print (' ')
        print ('optimized transfocator settings for E='+str(E)+' [eV]')
        print ('requested image position: '+str(image_pos)+' [m from source]   best match: '+str(image_pos+impos_diff)+' [m from source]')
        print ('difference in image location: '+str(impos_diff)+' [m], image '+loc+' of requested position')
        print ('requested demagnification: 1/'+str(zpos.mean()/(image_pos-zpos.mean()))+'  actual demagnification: 1/'+str(zmin/(image_pos-zmin)))
        print ('optimized  transfocator position: '+str(zmin)+' [m from source],   relative position: '+str(zmin-zpos.mean())+'[m]')
        print (' # slot7  slot6  slot5  slot4  slot3 slot2 slot1   <------ beam direction (slot 8 is currently B-fiber only)')
        print ('    '+str(conf_lensR[6])+'    '+str(conf_lensR[5])+'    '+str(conf_lensR[4])+'    '+str(conf_lensR[3])+'    '+str(conf_lensR[2])+'   '+str(conf_lensR[1])+'   '+str(conf_lensR[0])+'       lens radius [mm]')
        print ('    '+str(conf_lensN[6])+'    '+str(conf_lensN[5])+'    '+str(conf_lensN[4])+'    '+str(conf_lensN[3])+'    '+str(conf_lensN[2])+'   '+str(conf_lensN[1])+'   '+str(conf_lensN[0])+'       number of lenses')
        print (' ' )
    return {'lens_config':bin_index,'z_position':zmin}

def get_ip(conf,rel_z_pos,E='auto'):
    """
    get_ip(conf,real_z_pos,E='auto'): funtion to calculate the image position from current CHX transfocator configuration
    required aruments:
    conf: 'binary' array, reflecting the status of the lens stacks
    !!! omit slot #8 (B-fiber)
    e.g. [0,0,1,0,0,0,0] -> lens stack #3 in, all other out
    rel_z_pos: relative position of transfocator [m] ('0' -> 35.5m from source)
    optional arguments:
    E='auto': get beamline energy from PV (currently: DCM) | E=7894: manual overwrite of energy parameter
    Note: for 'virtual images' the image position is returned as 'NaN'
    function returns image position in [m] from source
    dependencies: imports numpy and xfuncs
    calls trans_setup() to get physical configuration of transfocator system
    by LW 03/15/2016
    """
    # transfocator setup
    trans_conf=trans_setup()
    lens_R = np.array(trans_conf['lens_radius'])*1.
    zpos = np.array(trans_conf['trans_position'])*1.
    lens_N = np.array(trans_conf['lens_number'])
    lens_mat = np.array(trans_conf['lens_material'])

   # checks:
    conf=np.array(conf)
    #print conf.min<0
    #print conf.max()>1
    # print len(conf) != len(lens_R)
    if len(conf) != len(lens_R):
        print ('thats interesting')
    if E is 'auto':
      #E=8000   # temporary: don't have channel access -> set E to 8000eV
      E=caget('XF:11IDA-OP{Mono:DCM-Ax:Energy}Mtr.RBV')     ### get energy automatically with channel access
      print ('getting energy from global PV: E=',E,'eV')# (currently not implemented in test version (no channel access) -> 8000eV default)'   # future: add PV name for house keeping
      if E> 30000 or E< 2000:
          raise transfuncs_Exception("error: Input argument E has to be 2000<E<30000 [eV]")
    E=np.array(E)*1.
    if conf.min()<0 or conf.max()>1 or len(conf) != len(lens_R):
        raise transfuncs_Exception("error: Input argument conf has to be an array or list of the form [1,0,1,...] and its lenght has to match the number of slots occupied by lenses")
    if rel_z_pos>zpos.max()-zpos.mean() or rel_z_pos<zpos.min()-zpos.mean():
        raise transfuncs_Exception("error: relative transfocator position has to be ["+str(zpos.min()-zpos.mean())+":"+str(zpos.max()-zpos.mean())+"] m")

    # calculate focusing:
    delta=np.zeros(len(lens_R))
    for w in range(0, len(lens_R)):
        delta[w]=np.real(xf.get_n(lens_mat[w],E/1000.))
    ln=lens_N*conf
    F=1/np.sum(2*delta*ln/(lens_R*1E-3))
    q=1/(-1/(rel_z_pos+zpos.mean()) + 1/F)
    imp=rel_z_pos+zpos.mean()+q
    print ('requested calculation for E= '+str(E)+'[eV], relative position of transfocator: '+str(rel_z_pos)+'[m]')
    print (' # slot7  slot6  slot5  slot4  slot3 slot2 slot1   <------ beam direction (slot 8 is currently B-fiber only)')
    print ('     '+str(conf[6])+'      '+str(conf[5])+'       '+str(conf[4])+'      '+str(conf[3])+'      '+str(conf[2])+'     '+str(conf[1])+'     '+str(conf[0])+'    (1: lens in| 0: lens out) ')
    if q<0:
        print ('image is imaginary...'+'(@'+str(imp)+' [m from source])')
    else: print ('image position: '+str(imp)+' [m from source]')
    return imp


class transfuncs_Exception(Exception):
    pass
    """
    by LW 03/13/2016
    class to raise transfuncs specific exceptions
    """
