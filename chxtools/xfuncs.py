# -*- coding: utf-8 -*-
"""
Created on Fri Apr 05 16:19:42 2013
package contains standard x-ray related functions. Needs to be linked to database!!!
@author: lwiegart
version 0.1: included help function in all macros [12/01/2014]
version 0.1.1: added get_gap & get_Es
version 0.1.2: added get_EBragg
version 0.1.3: fixed error messages in get_gap and get_Es
version 0.1.4: added lattice spacings for Si111 and Si220 @80K (from XOP) (02/26/2015)
version 0.2.0: fixed some problems wiht arrays as inputs, e.g. in get_gap  
               and changed the help files to use the standard python 'help' and docstring functions;
               added xfuncs_exeption class for handling xfuncs specific exeptions within python (03/23/2015)
               added get_pinflux: calculates photon flux from photo current from pin diode
version 0.2.1: (this file) changed ID selection for gap function to 'default_id'
"""

import pylab as pl
import numpy as np
import re

from pkg_resources import resource_filename as rs_fn
from pathlib import Path

# path to X-ray data files
# This is a Path object:
datapath = Path(rs_fn('chxtools', 'X-ray_database'))
xdatafiles = [str(f.relative_to(datapath)) for f in datapath.glob('*') if (Path(datapath) / Path(f)).is_file()]   
default_id = 'CHX_IVU20_05272017'

def get_Lambda(E,u='SI'):
    """
    calculates X-ray wavelength as a function of Energy [keV] in optional units.
    Syntax: getLambda(E,u), 
    where E=X-ray energy; optional: u= 'A','nm','um','cm','mm','m','SI' (='m'), default in the absence of u: 'SI'  
    
    """
    hPlank=6.62606876e-34;
    cvac=2.99792458e8;
    Qelectron=1.602176463e-19;
    scale=1
    #l=hPlank*cvac/(E*1000*Qelectron)
    l=hPlank*cvac/(E*1000*Qelectron);
    if u is 'A':
            scale=1e10;return l*scale # Angstroem
    elif u is 'nm':
            scale=1e9; return l*scale # nm
    elif u is 'um':
            scale=1e6; return l*scale # um
    elif u is 'mm':
            scale=1e3; return l*scale # mm
    elif u is 'cm':
            scale=1e2; return l*scale # cm
    elif u is 'm' or u is 'SI':
            scale=1; return l*scale
    else:
            print ('invalid option, type "get_Lambda(\'?\')" for available options and syntax')
  
   
def get_saxsQ(theta,E=8):
    """
    written by LW 13/03/2008.
    Calculates Q value for SAXS experiments [1/A], 
    input parameters: scattering angle theta [deg] (not theta/2!!!), 
    X-ray energy in keV (default: 8keV) function is vector compatible 
    calling sequence: get_saxsQ(theta,E)
    """    
    lam=get_Lambda(E,'A')
    return 4*np.pi/lam*np.sin(np.rad2deg(theta/2));

def get_ac(material,E=8):
    """
    by LW 10/03/2010
    function calculates the critical angle for total external reflection as a function of
    the material and the X-ray energy according to ac=sqrt(2*delta)
    index of refraction is a .dat file from http://henke.lbl.gov/optical_constants/getdb2.html
    (energy range: 2-30keV,delete the header % lines, name the file n_material.dat) % 
    calling sequence: ac=get_ac(material,E) where ac: critial angle in degrees, E [keV] (default: 8keV)
    type get_ac(\'materilal?\') to show list of supported materials"
    """
    
    #get list_of supported materials from data file directory:
    name=[]
    for i in range(len(xdatafiles)):
        m=re.search('(?<=n_)\w+', xdatafiles[i])
        if m is not None:
            name.append(m.group(0))             
    
    E=np.array(E)
    if material in name:
        loadn = Path(datapath) / Path('n_{}.dat'.format(material))
        n = pl.loadtxt(loadn, comments='%')
        if np.min(E)>=np.min(n[:,0]/1000) and np.max(E)<=np.max(n[:,0]/1000):
            d=np.interp(E*1000,n[:,0],n[:,1])
            return np.degrees(np.sqrt(2*d))
        else: print ('error: energy '+"%3.4f" %E +'[keV] out of range ('+"%3.4f" % np.min(n[:,0]/1000)+'=<E<='+"%3.4f" % np.max(n[:,0]/1000)+'keV)')
    elif material=='material?':
        print ('list of supported materials (based on data files in directory '+datapath+':')
        print (name   )    
    else: print ('error: non recognized material, please create index of refraction file first. Type "get_ac?" for instructions; type get_ac("material?") for list of supported materials') 
        
def get_n(material,E=8):
    """
    "by LW 07/04/2011 function get the index of refraction from stored data file,
    index of refraction is a .dat file from http://henke.lbl.gov/optical_constants/getdb2.html 
    (energy range: 2-30keV,delete the header lines, name the file n_material.dat) 
    calling sequence: n=get_n(material,E) where n is the complex refractive index detlta-i*beta, E: X-ray energy in keV"
    """
    #get list_of supported materials from data file directory:
    name=[]
    for i in range(len(xdatafiles)):
        m=re.search('(?<=n_)\w+', xdatafiles[i])
        if m is not None:
            name.append(m.group(0))             
    
    E=np.array(E)
    if material in name:
        loadn = Path(datapath) / Path('n_{}.dat'.format(material))
        n = pl.loadtxt(loadn, comments='%')
        if np.min(E)>=np.min(n[:,0]/1000) and np.max(E)<=np.max(n[:,0]/1000):
            d=np.interp(E*1000,n[:,0],n[:,1])
            b=np.interp(E*1000,n[:,0],n[:,2])
            return d-1j*b
        else: print ('error: energy '+"%3.4f" %E +'[keV] out of range ('+"%3.4f" % np.min(n[:,0]/1000)+'=<E<='+"%3.4f" % np.max(n[:,0]/1000)+'keV)')
    elif material=='material?':
        print ('list of supported materials (based on data files in directory '+datapath+':')
        print (name)
    else: print ('error: non recognized material, please create index of refraction file first. Type "get_n?" for instructions; type get_n("material?") for list of supported materials')

def get_mu(material,E=8):
    """
    by LW 07/04/2011
    function gets the attenuation length from stored data file, 
    attenuation length is a .dat file from http://henke.lbl.gov/optical_constants/getdb2.html
    (energy range: 2-30keV,delete the header lines or comment with '%', name the file n_material.dat)
    calling sequence: mu=get_mu(material,E) where mu [MICRONS!!!] is the 1/e attenuation length, E: X-ray energy in keV'
    """
    #get list_of supported materials from data file directory:
    name=[]
    for i in range(len(xdatafiles)):
        mm=re.search('(?<=mu_)\w+', xdatafiles[i])
        if mm is not None:
            name.append(mm.group(0))             
    
    E=np.array(E)
    if material in name:
        loadn = Path(datapath) / Path('mu_{}.dat'.format(material))
        m = pl.loadtxt(loadn, comments='%')
        if np.min(E)>=np.min(m[:,0]/1000) and np.max(E)<=np.max(m[:,0]/1000):
            mu=np.interp(E*1000,m[:,0],m[:,1])
            return mu
        else: print ('error: energy '+"%3.4f" %E +'[keV] out of range ('+"%3.4f" % np.min(m[:,0]/1000))+'=<E<='+"%3.4f" % np.max(m[:,0]/1000)+'keV)'
    elif material=='material?':
        print ('list of supported materials (based on data files in directory '+datapath+':')
        print (name)
    else: print ('error: non recognized material, please create index of refraction file first. Type get_mu("?") for instructions; type get_n("material?") for list of supported materials')
    
def get_T(material,E=8,l=1):
    """
    by LW 10/03/2010, 
    function calculates the transmission as a function of the material and the X-ray energy according to e^(-mul),
    where mu=4pi/lambda*beta 
    index of refraction is a .dat file from http://henke.lbl.gov/optical_constants/getdb2.html 
    (energy range: 2-30keV,delete the header lines, name the file n_material.dat) 
    calling sequence: T=get_T(material,E,l) 
    where T: transmission, material: E: X-ray energy in keV, l: thickness of the material [um],
    either E or l can be vectors; type get_T(\"material?\") for a list of supported materials
    """
    #get list_of supported materials from data file directory:
    name=[]
    for i in range(len(xdatafiles)):
        m=re.search('(?<=n_)\w+', xdatafiles[i])
        if m is not None:
            name.append(m.group(0))             
    
    E=np.array(E)
    l=np.array(l)
    #if len(E)==1 or len(l)==1:
    if E.size==1 or l.size==1:
        if material in name:
            loadn = Path(datapath) / Path('n_{}.dat'.format(material))
            n = pl.loadtxt(loadn, comments='%')
            if np.min(E)>=np.min(n[:,0]/1000) and np.max(E)<=np.max(n[:,0]/1000):
                b=np.interp(E*1000,n[:,0],n[:,2])
                mu=4*np.pi/get_Lambda(E,'um')*b;
                return np.exp(-mu*l);
            else: print ('error: energy '+"%3.4f" %E +'[keV] out of range ('+"%3.4f" % np.min(n[:,0]/1000)+'=<E<='+"%3.4f" % np.max(n[:,0]/1000)+'keV)')
        elif material=='material?':
            print ('list of supported materials (based on data files in directory '+datapath+':')
            print (name)
        else: print ('error: non recognized material, please create index of refraction file first. Type "get_T?" for instructions; type get_T("material?") for list of supported materials')
    else: print ('error: either energy or length must be a scalar, cannot scan both energy and length at the same time.')

def get_Bragg(reflection,E=8.):
    """
     by LW 17/03/2010
     function return the Bragg angle [deg.] of a given material and reflection at a given Energy. 
     Calling sequence: thetaB=get_Bragg(reflection,E),  thetaB(1)=Bragg angle[deg.] thetaB(2)=dhkl [A], thetaB(3)=I/Io [%].
     E: x-ray energy in keV (can be an array of energies),
     reflection: string, e.g. 'Si111'. Reflections implemented from http://database.iem.ac.ru/mincryst, T=25C or calculated from XOP, e.g for Si111&Si220 @80K
     type get_Bragg(\'reflections?\') for a list of currently availabel materials;
     """
    reflstr=['Si111cryo','Si220cryo','Si111', 'Si220', 'Si113', 'Si224', 'Si331', 'Si400','Ge111', 'Ge220', 'Ge113', 'Ge224', 'Ge331', 'Ge620', 'Ge531', 'Ge400', 'Ge115', 'Ge335','Ge440', 'Ge444', 'Ge333', 'C111', 'C220']
    dspace=np.array([3.13379852,1.91905183,3.13542,1.92004,1.63742,1.10854,1.24589,1.35767,3.26627,2.00018,1.70576,1.15480,1.29789,0.89451,0.95627,1.41434,1.08876,0.86274,1.00009,0.81657,1.08876,2.05929,1.26105])
    Irel=np.array([100,67.80,40.50,23.80,16.60,10.90,100,73.80,44.10,23.10,17.00,15.90,15.70,11.50,9.80,8.50,8.20,7.30,3.30,100,39.00])        
    if isinstance(reflection, str): # and all(isinstance(E, (int, long, float, complex)) for item in [E,E]): # <- bug in python: check for E is numeric works in standalone function, but not in this package => don't check
        E=np.array(E)
        lam=get_Lambda(E,'A')
        if reflection in reflstr:
            ind=reflstr.index(reflection)
            print (reflstr[ind] +': d_{hkl}=' + "%3.4f" %dspace[ind] +'A   I/I_o='+ "%3.4f" %Irel[ind])
            theta=np.degrees(np.arcsin(lam/2/dspace[ind]))
            ds=[];I=[]
            for l in range(0,np.size(theta)):
                ds.append(dspace[ind])
                I.append(Irel[ind])
            res=np.array([theta,np.array(ds),np.array(I)])
            return res.T
        elif reflection=='reflections?':
            print ('List of available reflections (T=25C):')
            print (reflstr )  
        else: print ('error: reflection not found. Type get_Bragg("reflections?") for list of available reflections.')
    else: print ('error: reflection has to be a string and E needs to be numeric. Type get_Bragg? for help')

def get_EBragg(reflection,theta_Bragg=12.0):
     """
     by LW 17/03/2010
     function returns the X-ray energy [keV] of a given crystal and Bragg reflection.
     Calling sequence: EBragg=get_EBragg(reflection,theta_Bragg (default: 12deg)),  thetaB(1)=Bragg angle[deg.] thetaB(2)=dhkl [A], thetaB(3)=I/Io [%]. 
     E: x-ray energy in keV (can be an array of energies), 
     reflection: string, e.g. 'Si111'. Reflections implemented from http://database.iem.ac.ru/mincryst, T=25C or calculated from XOP, e.g. for Si111 and Si220 @80K
     type get_Bragg(\'reflections?\') for a list of currently availabel materials
     """
     reflstr=['Si111cryo','Si220cryo','Si111', 'Si220', 'Si113', 'Si224', 'Si331', 'Si400','Ge111', 'Ge220', 'Ge113', 'Ge224', 'Ge331', 'Ge620', 'Ge531', 'Ge400', 'Ge115', 'Ge335','Ge440', 'Ge444', 'Ge333', 'C111', 'C220']
     dspace=np.array([3.13379852,1.91905183,3.13542,1.92004,1.63742,1.10854,1.24589,1.35767,3.26627,2.00018,1.70576,1.15480,1.29789,0.89451,0.95627,1.41434,1.08876,0.86274,1.00009,0.81657,1.08876,2.05929,1.26105])
     Irel=np.array([100,67.80,40.50,23.80,16.60,10.90,100,73.80,44.10,23.10,17.00,15.90,15.70,11.50,9.80,8.50,8.20,7.30,3.30,100,39.00])
     # constants needed:
     hPlank=6.62606876e-34
     cvac=2.99792458e8
     Qelectron=1.602176463e-19

     if isinstance(reflection,str): # and all(isinstance(E, (int, long, float, complex)) for item in [E,E]): # <- bug in python: check for E is numeric works in standalone function, but not in this package => don't check
         theta_Bragg=np.array(theta_Bragg)
         if reflection in reflstr:
             ind=reflstr.index(reflection)
             #print reflstr[ind] +': d_{hkl}=' + "%3.4f" %dspace[ind] +'A   I/I_o='+ "%3.4f" %Irel[ind]
             ds=[];I=[]
             for l in range(0,np.size(theta_Bragg)):
                 ds.append(dspace[ind])
                 I.append(Irel[ind])
             dspace=np.array(ds)
             lam=2*dspace*np.sin(theta_Bragg/180*np.pi)
             EthetaB=hPlank*cvac/(lam*Qelectron)*1e7;
             return EthetaB.T
         elif reflection=='reflections?':
             print ('List of available reflections (T=25C):')
             print (reflstr   )
         else: print ('error: reflection not found. Type get_EBragg("reflections?") for list of available reflections.')
     
     else: print ('error: reflection has to be a string and E needs to be numeric. Type get_Bragg? for help')

    
def get_Qpixel(center=[0,0],cord=[100,100],dis=1000,ps=75,E=8):
     """
     written by LW 13/03/2008. 
     calculates Q [1/A] value for SAXS experiments from detector position, 
     input parameters: position of the direct beam [x,y] in pixel numbers, pixel coordinate to calculte Q [x y], 
     sample detector distance in mm, pixel size [um], X-ray energy in keV.
     function is vector compatible for parameters dis OR ps OR E.
     calling sequence: get_Qpixel([xc,yc],[x,y],dis,ps,E)
     """
     center=np.array(center);cord=np.array(cord);dis=np.array(dis);ps=np.array(ps);E=np.array(E)      
     lam=get_Lambda(E,'A')
     xc=center[0]
     yc=center[1]
     x=cord[0]
     y=cord[1]
     ps=ps/1000 #convert pixel size to mm
     #real space distance on the CCD [mm]:
     u=np.sqrt(((x-xc)**2+(y-yc)**2))*ps
     #wavevector transfer [1/A]
     return 4*np.pi/lam*np.sin(np.arctan(u/dis/2))
         
def get_2theta(Qxy, E=8):
    """
    by LW 16/06/2005
    Returns the 2theta angle [deg] as a function of the input prameters Qxy [1/A] and Energy.
    Type get_2theta(Qxy [1/A], E [keV] (default: 8keV) ). This function is vector compatible.
    """
    Qxy=np.array(Qxy);E=np.array(E)
    lam=get_Lambda(E,'A')    
    return np.degrees(2*np.arcsin(Qxy/(4*np.pi)*lam))
    
def get_Gam(Qz,alpha_i=.12,E=8):
    hlpstr="by sof 01/12/2003   Returns the angle Gamma (out-of-plane angle in GIXD) defined in terms of Qz, the incident angle and the energy of the X-Ray beam. Type getGam(Qz [1/A], Incident Angle [deg] (default: 0.12deg) , E[keV] (default:8keV)), This function is vector compatible. Type get_Gam(\'?\') for help"
    if Qz=='?':
        print (hlpstr)
    else:
        Qz=np.array(Qz);alpha_i=np.array(alpha_i);E=np.array(E)
        lam=get_Lambda(E,'A')
        return np.degrees(lam*Qz/(2*np.pi)-np.sin(alpha_i))

def get_Qll(Ty,alpha_i=.12,d=1000,E=8):
    hlpstr="LW 26-01-2005 Function returns the parallel wavevector transfer in GI-XPCS geometry [cm^-1]. Type: get_Qll(Ty[mm],alpha_i[deg] (default: 0.12deg), d [mm] (default: 1000mm),E[keV] (default:8keV))=>qll[cm^-1]; type get_Qll(\'?\') for help"
    if Ty=='?':
        print ( hlpstr)
    else:
        Ty=np.array(Ty);alpha_i=np.array(alpha_i);E=np.array(E)
        lam=get_Lambda(E,'A')
        af=np.arctan((d*np.tan(np.deg2rad(alpha_i))+Ty)/d)
        return (2*np.pi/lam)*(np.cos(np.deg2rad(alpha_i))-np.cos(af))*1e8
    
def get_Qxy(theta,E=8):
    hlpstr="sof 26/11/2003   To return for in plane X-ray scattering the parallel wavevector transfer in 1/A. Type get_Qxy(FULL scattering angle [deg], E [keV] (default:8keV), the function returns the Q vector [1/A] in the surface plane. Angle measured in surface from specular direction! This function is vector compatible. Type get_QXY(\'?\') for help."
    if theta=='?':
        print (hlpstr)
    else:
        theta=np.array(theta);E=np.array(E)
        lam=get_Lambda(E,'A')
        return 4*np.pi*np.sin(np.deg2rad(theta)/2)/lam
    
def get_Qz(Gam,alpha_i=.12,E=8):
    hlpstr="function by sof 01/12/2003   Returns the Qz z-component of wavevector transfer defined in terms of the incident  and measured angles and the energy of the X-Ray beam.  Type get_Qz(Gam [deg], alpha_i [deg] (default: 0.12deg), E[keV] (default: 8keV)). This function is vector compatible. Type get_Qz(\'?\') for help."
    if Gam=='?':
        print (hlpstr)
    else:
        E=np.array(E);Gam=np.deg2rad(np.array(Gam));alpha_i=np.deg2rad(np.array(alpha_i))
        lam=get_Lambda(E,'A')
        return 2*np.pi*(np.sin(Gam)+np.sin(alpha_i))/lam

def ch2deg(chn,channel_per_degree=100,chOff=0):
    """
    LW 08-07-2005
    function converts the PSD channel (chn) into angles [deg],
    needed are the calibration parameters channels/deg (channel_per_degree)
    and the channel number of the direct beam (chOff).
    Calling sequence: ch2deg(chan,channel_per_degree (default: 100),chOff (default=0) -> corresponding angle [deg].
    Function is vector compatible.
    """
    chn=np.array(chn);channel_per_degree=np.array(channel_per_degree);chOff=np.array(chOff)
    return (chn-chOff)/channel_per_degree

def get_Vineyard(Qz,E=8,ra0=1):
    """
    by LW 20/07/2006
    function calculates the Vineyard (transmission) Function for a given (range) of Qz(s)
    and given X-ray energy and scattering length density ra0 of the subphase media.
    Type: getVineyard(Qz [A], E [keV] (default: 8keV), ra0 [SI] (default: 1 ).
    """
    lam=get_Lambda(E,'A')
    Qc=4*np.pi*np.sin(np.sqrt(np.pi*ra0)*lam/np.pi)/lam
    rx=Qz/Qc
    V=np.ones(np.size(Qz))
    for i in range(0,np.size(Qz)):
         if rx[i]<=1:
             V[i]=2*rx[i]
         else: V[i]=2*rx[i]/(rx[i]+np.sqrt(rx[i]**2-1))
    return V

def get_gap(E,harmonic=3,ID=default_id):
    """    
    by LW 12/03/2014, function calculates the undulator gap for a requested energy
    and harmonic based on magnetic measurement data in the database 
    type get_gap(E [kev], harmonic [integer] (default=3), id (default: defined by 'default_id')); 
    E can be an array of energies. Type get_gap\"ID?\") for a list of available magnetic datasets.
    """
    #get list_of available magnetic measurements from data file directory:
    name=[]
    for i in range(len(xdatafiles)):
        m=re.search('(?<=id_)\w+', xdatafiles[i])
        if m is not None:
            name.append(m.group(0))             

    if E=='ID?':
        print ('list of available magnetic measurements (based on data files in directory {}:'.format(str(datapath)))
        print (name)
    else:
        E=np.array(E)*1.0
        harm_check(harmonic)     
        if ID in name:
            loadn = Path(datapath) / Path('id_{}.dat'.format(ID))
            magdat = pl.loadtxt(loadn, comments='%')
            #harmonic=harmonic*1.0
            if np.min(E/harmonic)>=np.min(magdat[:,2]) and np.max(E/harmonic)<=np.max(magdat[:,2]):
                gap=np.interp(E/harmonic,magdat[:,2],magdat[:,0])
                return gap
            # this else should be a warning only and should return NaN    
            else: raise xfuncs_Exception ('error: energy '+"%3.4f" %E +'[keV] out of range for requested harmonic number, gap limit: ('+"%3.4f" % np.min(magdat[:,0])+'=<gap<='+"%3.4f" % np.max(magdat[:,0])+'mm), try using higher/lower harmonic number.')
        else: raise xfuncs_Exception('error: non recognized magnetic data. Type get_gap(\'ID?\') for a list of available mangetic datasets.')
       

def get_Es(gap,harmonic=[1,2,3,4,5],ID=default_id):
    """
    by LW 12/03/2014,
    function calculates the X-ray energies for a given undulator gap and set of harmonics
    based on magnetic measurement data in the database 
    type get_Es(gap [mm], harmonic [integer] (default=[1,2,3,4,5]), id (default: defined by 'default_id')); 
    harmonic can be a list of integers. Type get_Es(\"ID?\") for a list of available magnetic datasets
    """
    #get list_of available magnetic measurements from data file directory:
    name=[]
    for i in range(len(xdatafiles)):
        m=re.search('(?<=id_)\w+', xdatafiles[i])
        if m is not None:
            name.append(m.group(0))             
    if gap=='ID?':
        print ('list of available magnetic measurements (based on data files in directory '+datapath+':')
        print (name)
    else:
        for l in range(0, np.size(harmonic)):
            harm_check(harmonic)            
            if l==np.size(harmonic)-1:
                    gap=np.array(gap)
                    if ID in name:
                        loadn = Path(datapath) / Path('id_{}.dat'.format(ID))
                        magdat = pl.loadtxt(loadn, comments='%')
                        harmonic=np.array(harmonic)*1.0
                        if np.min(gap)>=np.min(magdat[:,0]) and np.max(gap)<=np.max(magdat[:,0]):
                            Es=np.interp(gap,magdat[:,0],magdat[:,2])
                            Eharmonics=np.array([harmonic,harmonic*Es])
                            return Eharmonics.T  
                        else: print ('error: gap '+"%3.4f" % gap +'[mm] out of range for gap ('+"%3.4f" % np.min(magdat[:,0])+'=<gap<='+"%3.4f" % np.max(magdat[:,0])+'mm), try higher/lower harmonic number.')
                    else: print ('error: non recognized magnetic data. Type get_gap(\'ID?\') for a list of available mangetic datasets.')

def get_pinflux(current,Energy,thickness=300):
    """
    by LW 03/26/2015
    function to calculate photon flux from pin-diode measurement
    uses scattering cross section for Si from NIST -> PhElAbsCross_si.dat in database
    assumes Silicon as material as other materials would be pretty exotic
    calling sequence: get_pinflux(current,Energy,thickness=300)
    current: current [A]
    Energy: X-ray energy [keV]
    thickness: Si diode thickness [um]
    """
    # some conversions and contstants:
    rho_Si=2.329   # density of Si in g/cm^3!!!
    epsilon=3.66   # energy for creation of electron-hole pair in Si [eV]    
    thickness=thickness/1.0E4  #conversion to cm
    # read datafile with scattering cross section:
    loadn=datapath+'PhElAbsCross_Si.dat'
    crossdat=pl.loadtxt(loadn,comments='%')
    # check for energy range:
    xmin=min(crossdat[:,0])*1E3
    xmax=max(crossdat[:,0])*1E3
    Energy=np.array(Energy)*1.0
    current=np.array(current)*1.0
    if np.max(Energy) > xmax or np.min(Energy) <xmin:
        raise xfuncs_Exception(['X-ray energy out of range for cross section data. ',xmin,'<=Energy<=,',xmax,' [keV]'])
    # calculate photo current PER 1E10 ph/s:
    PhCur=1E10*crossdat[:,0]*1e3*1.6022e-16*(1-np.exp(-crossdat[:,1]*thickness*rho_Si))/epsilon
    PhCurint=np.interp(Energy,crossdat[:,0]*1E3,PhCur)
    photon_flux=current/PhCurint*1E10
    print ('photo current for E= ',Energy,'keV: ',PhCurint*1E3,'mA/10^10ph/s')
    print ('flux for photo current ',current*1E3,'mA at E=',Energy,'keV: ',photon_flux,'ph/s')
    return photon_flux
            
####### Help functions
class xfuncs_Exception(Exception):
    pass
    """
    by LW 03/19/2015
    class to raise xfuncs specific exceptions
    """

def harm_check(harm_n):
    """
    by LW 03/19/2015
    simple function to check whether single value or list of values only contains
    integers >=1 (number of undulator harmonics)
    calling sequence: harm_check(harm_n), where harm_n is a single value or list of values
    function raises an exception is harm_n does contain values other than integers >=1
    """
    state=True
    
    if np.all(np.array(harm_n)>=1)!=True:
       raise xfuncs_Exception("Type Error: Harmonic numbers need to be >=1!") 
    
    try: 
        list(harm_n)         # b is the list to be checked
        for i in list(harm_n):
            #print isinstance(i,int)
            state=state*isinstance(i,int)
    except:      # check for b being a single integer
        try: 
                state=state*isinstance(harm_n, int)
                #print isinstance(harm_n, int)
        except: pass
    
    if state !=True:
        raise xfuncs_Exception("Type Error: Harmonic numbers need to be integers >=1!")
