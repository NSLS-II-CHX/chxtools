# -*- coding: utf-8 -*-
"""
by LW 03/20/2015
some wrapper functions for ophyd
"""
from __future__ import print_function


def det_select(det):
    """
    by LW 03/20/2015
    function to select detector for scanning and ct
    removes all existing detectors from ophyd function ascan.user_detectors
    and replaces them with the list of names passed in this function
    ONLY FOR SINGLE STRINGS (not lists) SO FAR! Only works if there is one default detector/trigger to start with
    => needs some work in the future
    calling sequence: det_select(det)
    """
    try:
        rm_det=ascan.user_detectors[0].name
        ascan.user_detectors.remove(session_mgr[rm_det])
        ascan.default_triggers = []
        ascan.default_detectors = []
        print('removed previous default detector: {}'.format(rm_det))
    except:
        print('list of detectors appeared to be empty...')

    ascan.user_detectors.append(det)        # add detector
    new_def_det=ascan.user_detectors[0].name

    print('\nnew default detector: {}'.format(new_def_det))


def cw_ascan(mot,xmin,xmax,npoints,acqt='default',pos_ret=True):
    """
    by LW 03/21/2015
    wrapper function for ophyd's acan
    added functionality:
        - default option: returns scanned axis back to the position where it was prior to starting the scan
        - option to change acquition times on detectors (default: use current setting from corresponding css screen)
        - stop camera and set image mode to single
        - after the scan: restarts camera and set image mode to where it was prior to starting the scan
    calling sequence: cw_ascan(det,xmin,xmax,npoints,acqt='default',pos_ret=True)
    WILL NOT WORK FOR A LIST OF DETECTORS!
    """
    # gather beamline information prior to starting the scan:
    ini_motpos = caget(mot.record+'.RBV')
    # current detector:
    acq_pv = session_mgr[ascan.user_detectors[0].name].pvname

    # set different exposure time for the scan, if requested:
    if acqt != 'default':
        try:
            ini_expt = caget(acq_pv[2]) # initial exposure time
            session_mgr[ascan.user_detectors[0].name].acquire_time = acqt
            print('successfully set exposure time to [s]: {}'.format(acqt))
        except:
            print('could not set exposure time to {}'.format(acqt))

    # execute the scan
    ascan(mot,xmin,xmax,npoints)

    # put beamline back into initial state
    if pos_ret:
        caput(mot.record+'.VAL',ini_motpos)
        print('returned axes to: {}'.format(ini_motpos))
    if acqt != 'default':
        try:
            session_mgr[ascan.user_detectors[0].name].acquire_time = ini_expt
            print('successfully reset exposure time to [s]: {}'.format(ini_expt))
        except:
            print('could not reset exposure time to {}'.format(ini_expt))

def cw_dscan(mot,mdx,pdx,npoints,acqt='default',pos_ret=True):
    """
    by LW 03/21/2015
    wrapper function for ophyd's ascan
    added functionality:
        - option to change acquition times on detectors (default: use current setting from corresponding css screen)
    calling sequence: cw_dscan(det,mdx,pdx,npoints,acqt='default',pos_ret=True)
    WILL NOT WORK FOR A LIST OF DETECTORS!
    """
    # current detector:
    acq_pv = session_mgr[ascan.user_detectors[0].name].pvname

    # set different exposure time for the scan, if requested:
    if acqt != 'default':
        try:
            ini_expt=caget(acq_pv[2]) # initial exposure time
            session_mgr[ascan.user_detectors[0].name].acquire_time = acqt
            print('successfully set exposure time to [s]: {}'.format(acqt))
        except:
            print('could not set exposure time to {}'.format(acqt))

    # execute the scan
    dscan(mot,mdx,pdx,npoints)
    #print('finished scan')

    if acqt != 'default':
        try:
            session_mgr[ascan.user_detectors[0].name].acquire_time = ini_expt
            print('successfully reset exposure time to [s]: {}'.format(ini_expt))
        except:
            print('could not reset exposure time to {}'.format(ini_expt))

def cw_CCDseries(folder,filename,detector,imnum='default',startn=1,acqt='default',acqperiod='default'):
    """
    by LW 04/06/2015
    wrapper function to take a time series with a CCD detector
    functionalities:
    - select/create folder via GUI (name for new folder in the GUI, new folder will be created)
    - select CCD
    - gather current CCD acquisition parameters and restore these after the timeseries
    - select acquisition time and acquiering period ('default': use current parameters from area detector
    - switch autosave on before the series and off when done
    - select start number for image series !!if file existed already, it will be overwritten !!
    - currently saving .tiff as default
    calling sequence: cw_CCDseries(folder,filename,detector,imnum='default',startn=1,acqt='default',acqperiod='default')
    """
    import time
    import Tkinter, tkFileDialog
    import os, stat, sys
    # get the detector name::
    detector.pvname.split('}')[0]+'}'
    #get folder interactively:

    if folder == 'ia':
      root=Tkinter.Tk()
      root.withdraw()
      directory=tkFileDialog.askdirectory()
      folder=directory
    # check whether target directory exists and create it, if it doesn't
    if not os.path.exists(folder):
      os.mkdir(folder)
      os.chmod(folder,436)      #make sure everybody can read an write
      os.chmod(folder,stat.S_IRWXO)
      print('successfully created new directory: {}'.format(folder))

    # put folder:
    r=caput(detector.pvname.split('}')[0]+'}TIFF1:FilePath', folder)      # have some problem with syntax here...
    if r == 1:
        print('changed saving directory to: {}'.format(folder))
    else:
        print('error: could not change directory for saving data files.')
    # put the filename:
    r=caput(detector.pvname.split('}')[0]+'}TIFF1:FileName', filename)      # have some problem with syntax here...
    if r == 1:
        print('filename for saving: {}'.format(filename))
    else:
        print('error: could not change file name for saving.')

    # put start number
    caput(detector.pvname.split('}')[0]+'}TIFF1:FileNumber',startn)

    #gather information about current camera settings
    acq_pv = session_mgr[ascan.user_detectors[0].name].pvname
    ini_acq = caget(acq_pv.split('}')[0]+'}cam1:Acquire')   # initial state: started or stopped
    ini_mode = caget(acq_pv.split('}')[0]+'}cam1:ImageMode')  # initial image mode: single [0], multiple [1], continous [2]
    ini_expt = caget(acq_pv.split('}')[0]+'}cam1:AcquireTime') # initial exposure time
    ini_acqperiod = caget(acq_pv.split('}')[0]+'}cam1:AcquirePeriod')  # initial acquiring period
    ini_imnum = caget(acq_pv.split('}')[0]+'}cam1:NumImages')  # initial image number

    if acqt != 'default':
        try:
            caput(acq_pv.split('}')[0]+'}cam1:AcquireTime',acqt)
            print('successfully set exposure time to [s]: {}'.format(acqt))
        except:
            print('could not set exposure time to {}'.format(acqt))
    # stop camara:
    try:
        caput(acq_pv.split('}')[0]+'}cam1:Acquire',0)
        print('successfully stopped camera')
    except:
        print('could not stop camera')
    # try to set image mode to multiple
    try:
        caput(acq_pv.split('}')[0]+'}cam1:ImageMode',1)
        print('successfully set ImageMode to "multiple"')
    except:
        print('could not set ImageMode to "multiple"')
    if acqperiod != 'default':
        try:
            caput(acq_pv.split('}')[0]+'}cam1:AcquirePeriod',acqperiod)
            print('successfully set acquiering period to: {}'.format(acqperiod))
        except:
            print('could not set aquireing period to {}'.format(acqperiod))

    # set number of images to be taken:
    if imnum != 'default':
        try:
            caput(acq_pv.split('}')[0]+'}cam1:NumImages', imnum)
            print('successfully set number of images to: {}'.format(imnum))
        except:
            print('could not set number of images to {}'.format(imnum))
    print('going to start the acquisition...')
    time.sleep(1)
    # start the acquisition!!!
    caput(acq_pv.split('}')[0]+'}cam1:ArrayCounter',0)  # set Images acquired (image counter) to 0

    caput(detector.pvname.split('}')[0]+'}TIFF1:AutoSave',1)    # start auto save
    caput(acq_pv.split('}')[0]+'}cam1:Acquire',1)               # start series acquisition

    counting=0
    current_file=caget(detector.pvname.split('}')[0]+'}TIFF1:FullFileName_RBV',as_string=True)
    while counting==0:
        time.sleep(.5)
        if caget(acq_pv.split('}')[0]+'}cam1:ArrayCounter_RBV') == imnum:
            counting=1
        if caget(detector.pvname.split('}')[0]+'}TIFF1:FullFileName_RBV',as_string=True) != current_file:
          current_file=caget(detector.pvname.split('}')[0]+'}TIFF1:FullFileName_RBV',as_string=True)
          print('file written: {}'.format(current_file))
        time.sleep(.1)

    print('going to stop the acquisition...')
    time.sleep(.5)
    caput(acq_pv.split('}')[0]+'}cam1:Acquire',0)               # stop series acquisition (redundent...should have stopped already
    caput(detector.pvname.split('}')[0]+'}TIFF1:AutoSave',0)    # stop auto save

    print('going to reset camera...')
    #time.sleep(5)
     # put camera back into initial state:
    if acqt != 'default':
        try:
            caput(acq_pv.split('}')[0]+'}cam1:AcquireTime',ini_expt)
            print('successfully reset exposure time to [s]: {}'.format(ini_expt))
        except:
            print('could not reset exposure time to {}'.format(ini_expt))
    try:
        caput(acq_pv.split('}')[0]+'}cam1:ImageMode', ini_mode)
        print('successfully reset ImageMode')
    except:
        print('could not reset ImageMode')
    try:
        caput(acq_pv.split('}')[0]+'}cam1:Acquire', ini_acq)
        print('successfully reset camera acqusition mode')
    except:
        print('could not reset camera acquisition mode')
    if acqperiod!='default':
        try:
            caput(acq_pv.split('}')[0]+'}cam1:AcquirePeriod', ini_acqperiod)
            print('successfully reset acquisition period to [s]: {}'.format(ini_acqperiod))
        except:
            print('could not reset acquisition period  to {}'.format(ini_acqperiod))
    if imnum!='default':
        try:
            caput(acq_pv.split('}')[0]+'}cam1:NumImages', ini_imnum)
            print('successfully reset image numbers to: {}'.format(ini_imnum))
        except:
            print('could not reset image numbers to {}'.format(ini_imnum))
    time.sleep(.5)
    try:
        caput(acq_pv.split('}')[0]+'}cam1:Acquire', ini_acq)      # restart camera if it was running before taking the series
        print('restarted camera')
    except:
        print('could not restart camera...')

