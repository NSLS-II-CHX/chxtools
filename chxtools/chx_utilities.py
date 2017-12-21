# -*- coding: utf-8 -*-
"""
Created on Wed Mar 25 14:02:59 2015
by LW March 2015
set of utility functions for beamline alingment and commissioning
v 0.0.1 (this version): might have created a typo in E-calibration!!!
                        added dcm_roll for calculating DCM Roll correction
"""

from pathlib import Path
from pkg_resources import resource_filename as rs_fn


############
##################
####
def E_calibration(file,Edge='Cu',xtal='Si111cryo',B_off=0):
    """
    by LW 3/25/2015
    function to read energy scan file and determine offset correction
    calling sequence: E_calibration(file,Edge='Cu',xtal='Si111cryo',B_off=0)
    file: path/filename of experimental data; 'ia' opens interactive dialog
    Edge: elment used for calibration
    xtal: monochromator crystal under calibration
    B_off (optional): apply offset to Bragg angle data
    currently there is no check on input parameters!
    """
    # read the data file
    import csv
    import numpy as np
    import matplotlib.pyplot as plt
    import xfuncs as xf
    import Tkinter, tkFileDialog

    if file=='ia':          # open file dialog
        root = Tkinter.Tk()
        root.withdraw()
        file_path = tkFileDialog.askopenfilename()
    else: file_path=file

    plt.close("all")
    Bragg=[]
    Gap=[]
    Intensity=[]
    Edge_data={'Cu': 8.979, 'Ti': 4.966}
    with open(file_path, 'rb') as csvfile:
        filereader = csv.reader(csvfile, delimiter=' ')
        filereader.next()   # skip header lines
        filereader.next()
        filereader.next()
        for row in filereader:              # read data
            try: Bragg.append(float(row[2]))
            except: print('could not convert: ',row[2])
            try: Gap.append(float(row[5]))
            except: print('could not convert: ',row[5])
            try: Intensity.append(float(row[7]))
            except: print('could not convert: ',row[8])
    B=np.array(Bragg)*-1.0+B_off
    G=np.array(Gap[0:len(B)])   # not currently used, but converted for future use
    Int=np.array(Intensity[0:len(B)])

    # normalize and remove background:
    Int=Int-min(Int)
    Int=Int/max(Int)

    plt.figure(1)
    plt.plot(B,Int,'ko-',label='experimental data')
    plt.plot([xf.get_Bragg(xtal,Edge_data[Edge])[0],xf.get_Bragg(xtal,Edge_data[Edge])[0]],[0,1],'r--',label='Edge for: '+Edge)
    plt.legend(loc='best')
    plt.xlabel(r'$\theta_B$ [deg.]')
    plt.ylabel('intensity')
    plt.title(['Energy Calibration using: ',file])
    plt.grid()

    plt.figure(2)
    Eexp=xf.get_EBragg(xtal,B)
    plt.plot(Eexp,Int,'ko-',label='experimental data')
    plt.plot([Edge_data[Edge],Edge_data[Edge]],[0,1],'r--',label='Edge for: '+Edge)
    plt.legend(loc='best')
    plt.xlabel('E [keV.]')
    plt.ylabel('intensity')
    plt.title(['Energy Calibration using: ',file])
    plt.grid()

    # calculate derivative and analyze:
    Bragg_Edge=xf.get_Bragg(xtal,Edge_data[Edge])[0]
    plt.figure(3)
    diffdat=np.diff(Int)
    plt.plot(B[0:len(diffdat)],diffdat,'ko-',label='diff experimental data')
    plt.plot([Bragg_Edge,Bragg_Edge],[min(diffdat),max(diffdat)],'r--',label='Edge for: '+Edge)
    plt.legend(loc='best')
    plt.xlabel(r'$\theta_B$ [deg.]')
    plt.ylabel('diff(int)')
    plt.title(['Energy Calibration using: ',file])
    plt.grid()

    plt.figure(4)
    plt.plot(xf.get_EBragg(xtal,B[0:len(diffdat)]),diffdat,'ko-',label='diff experimental data')
    plt.plot([Edge_data[Edge],Edge_data[Edge]],[min(diffdat),max(diffdat)],'r--',label='Edge for: '+Edge)
    plt.legend(loc='best')
    plt.xlabel('E [keV.]')
    plt.ylabel('diff(int)')
    plt.title(['Energy Calibration using: ',file])
    plt.grid()

    edge_index=np.argmax(diffdat)
    B_edge=xf.get_Bragg(xtal,Edge_data[Edge])[0]

    print('')
    print('Energy calibration for: ',file_path)
    print('Edge used for calibration: ',Edge)
    print('Crystal used for calibration: ',xtal)
    print('Bragg angle offset: ', B_edge-B[edge_index],'deg. (CHX coordinate system: ',-(B_edge-B[edge_index]),'deg.)')
    print('=> move Bragg to ',-B[edge_index],'deg. and set value to ',-Bragg_Edge,'deg.')
    print( 'Energy offset: ',Eexp[edge_index]-Edge_data[Edge],' keV')

def dcm_roll(Bragg,offset,distance,offmode='mm',pixsize=5.0):
    """
    by LW 03/27/2015
    function to calculate Roll correction on the DCM
    calling sequence: dcm_roll(Bragg,offset,distance,offmode='mm',pixsize=5.0)
    Bragg: set of Bragg angles
    offset: set of corresponding offsets
    offmode: units of offsets = mm or pixel (default:'mm')
    pixsize: pixel size for offset conversion to mm, if offsets are given in pixels
    default is 5um (pixsize is ignored, if offmoe is 'mm')
    distance: DCM center of 1st xtal to diagnostic/slit [mm]
    preset distances available: 'dcm_bpm',dcm_mbs', 'dcm_sample'
    """
    import numpy as np
    from scipy import optimize
    from matplotlib import pyplot as plt
    Bragg=np.array(Bragg)
    if offmode=='mm':
        offset=np.array(offset)
    elif offmode=='pixel':
        offset=np.array(offset)*pixsize/1000.0
    else: raise CHX_utilities_Exception('Eror: offmode must be either "mm" or "pixel"')
    if distance=='dcm_bpm':
        d=3000.0 # distance dcm-bpm in mm
    elif distance=='dcm_mbs':
        d=2697.6 #distance dcm-mbs in mm
    elif distance=='dcm_sample':
        d=16200 #distance dcm-sample in mm
    else:
        try:
            d=float(distance)
        except:
            raise CHX_utilities_Exception('Eror: distance must be a recognized string or numerical value')

    # data fitting
    fitfunc = lambda p, x: p[0]+2*d*p[1]*np.sin(x/180.*np.pi) # Target function
    errfunc = lambda p, x, y: fitfunc(p, Bragg) - y # Distance to the target function
    p0 = [np.mean(offset), -.5] # Initial guess for the parameters
    p1, success = optimize.leastsq(errfunc, p0[:], args=(Bragg, offset))

    # plotting the result:
    plt.close(1)
    plt.figure(1)
    B = np.linspace(Bragg.min(), Bragg.max(), 100)
    plt.plot(Bragg,offset,'ro',label='measured offset')
    plt.plot(B,fitfunc(p1,B),'k-',label=r'$x_o$+2*D*$\Delta$$\Phi$*sin($\theta_B$)')
    plt.legend(loc='best')
    plt.ylabel('beam offset [mm]')
    plt.xlabel('Bragg angle  [deg.]')
    print('x_0= ',p1[0],'mm')
    print('\Delta \Phi= ',p1[1]*180.0/np.pi,'deg')

def sum_image(filename, firstim=0, lastim=9999):
    """
    by LW 03/30/2015
    function to sum up tiff files and display cuts
    calling sequence: sum_image(filename, firstim=0, lastim=9999)
    filename can be 'ia' for chosing path and first image via GUI
    """
    from PIL import Image
    import Tkinter, tkFileDialog
    from matplotlib import pyplot as plt
    import numpy as np
    import time

    if filename=='ia':          # open file dialog
        root = Tkinter.Tk()
        root.withdraw()
        file_path = tkFileDialog.askopenfilename()
    else: file_path=filename
    print(file_path)

    # extract filename and first image number:
    a=file_path.split('/')
    fn=a[len(a)-1].split('.')[0]  # file name
    b=fn.split('.')
    #fe=b[len(b)-1]  # file ending, e.g .tiff
    fe=a[len(a)-1].split('.')[1]
    c=fn.split('_')
    fnum=c[len(c)-1].split('.')[0]    # file number
    nlen=len(fnum)  # length of number format
    fp=file_path.split(fn)[0]
    fn=fn.split('_'+fnum)[0]
    print('file path: ',fp)
    print('file name: ',fn)
    print('file ending: ',fe)
    print('file number: ',fnum)

    if filename=='ia':      # get first image number from filebame in interactive mode
        firstim=int(fnum)
        #lastim=30     # maximum number of images to search in autmatic mode

    imcount=firstim
    #print 'imcount: ',imcount
    #print 'firstim: ',firstim
    get_out=0
    image_data=0
    while imcount <= lastim:
        ic=str(imcount)
        while len(ic) < nlen:         # create file number
                ic='0'+ic
                #ic=str(ic)
                #print ic
        #assemble the image name to be loaded:
        file_to_open=fp+fn+'_'+ic+'.'+fe

        try:
            im = Image.open(file_to_open)
            #im.show()
            print('loaded file: ',file_to_open)
        except:
            print('waiting for file: '+file_to_open)
            #plt.close(1)    # use the waiting time to plot some intermediate result -> does not work, as python is waiting until the end before printing....
            #plt.figure(1)
            #imgplot= plt.imshow(image_data/(imcount-firstim+1),origin='lower')
            #imgplot.set_cmap('gray')
            time.sleep(5)
            try:
                im = Image.open(file_to_open)
                #im.show()
            except:get_out=1
        if get_out==1:
            return image_data
            plt.close(1)    # only plot at the end for now: too much overhead to plot every iteration
            plt.figure(1)
            raise CHX_utilities_Exception('no more images found, images summed up to: '+file_to_open)
        if imcount==firstim:
            # print 'made it here!'
            im=np.array(im)*1.0
            image_data=np.sum(im, axis=2)
            imgplot= plt.imshow(image_data,origin='lower')
            imgplot.set_cmap('gray')
        else: image_data=image_data+np.array(np.sum(im, axis=2))
        #plt.close(1)
        #plt.figure(1)
        #plt.imshow(image_data/((imcount-firstim)+1))
        imcount=imcount+1
    plt.close(1)    # only plot at the end for now: too much overhead to plot every iteration
    plt.figure(1)
    #finalimage=Image.fromarray(image_data/(imcount-firstim+1))
    imgplot= plt.imshow(image_data/(imcount-firstim+1),origin='lower')
    imgplot.set_cmap('gray')
    plt.title('summed :'+fp+fn+fnum+'-'+ic)
    #finalimage.show()
    return image_data/(imcount-firstim+1)

def gauss(x,*p):
    """
    LW 04/06/2015
    Gaussian help function of the form:
    yo+A*exp(-(x-xc).^2/(2*w^2))
    """
    import numpy as np
    #print p
    yo,A,xc,w=p
    return yo+A*np.exp(-(x-xc)**2/(2.*w**2))


###########################################
##New functions added by Yugang Zhang in June,2015

#####PVs
PV_dict = dict(
            #Acc_Vert='SR:C31-{AI}Aie11:Angle-y-Cal',
            T_DCM_Bragg = 'XF:11IDA-OP{Mono:DCM-Ax:Bragg}T-I',
            T_DCM_Pitch = 'XF:11IDA-OP{Mono:DCM-Ax:Pitch}T-I',
            T_DCM_X='XF:11IDA-OP{Mono:DCM-Ax:X}T-I',
            T_DCM_R = 'XF:11IDA-OP{Mono:DCM-Ax:Roll}T-I',
            T_DCM_WaterOut= 'XF:11IDA-OP{Mono:DCM}T:H2O_Out-I',
            T_DCM_plate1= 'XF:11IDA-OP{Mono:DCM}T:Plate1-I',
            T_Cryo = 'XF:11IDA-UT{Cryo:1}T:05-I',
            P_Cryo = 'XF:11IDA-UT{Cryo:1}P:03-I',
            #DCM_Pitch_Encoder='XF:11IDA-OP{Mono:DCM-Ax:P}Mtr.DRBV',
            E_DCM_Pitch_RBV='XF:11IDA-OP{Mono:DCM-Ax:P}Mtr.RBV',
            E_DCM_Pitch_OFF='XF:11IDA-OP{Mono:DCM-Ax:P}Mtr.OFF',
            E_DCM_Pitch_ERES='XF:11IDA-OP{Mono:DCM-Ax:P}Mtr.ERES',
            E_HDM_Pitch ='XF:11IDA-OP{Mir:HDM-Ax:P}Pos-I',
          )

PV_dict = dict(
            E_HDM_Pitch ='XF:11IDA-OP{Mir:HDM-Ax:P}Pos-I',
          )

label =  PV_dict.keys()
PV = [ PV_dict[lab] for lab in label    ]

PV_Ind = dict(
            Freq_Cryo = 'XF:11IDA-UT{Cryo:1-Pmp}Freq-I',
        )

labeli =  PV_Ind.keys()
PVi = [ PV_Ind[lab] for lab in labeli    ]


#############################
def determine_aspect(shape, extent):
    dx = (extent[1] - extent[0]) / float(shape[1])
    dy = (extent[3] - extent[2]) / float(shape[0])
    return dx / dy


def get_cuts(img, cut=1000,direction='horz'):
    if direction=='horz':
        cdat=img[cut,:];
    elif direction=='vert':
        cdat=img[:,cut];
    return cdat

def show_cut(img, cut=1000, direction='horz',width=10,title=None):
    #fig,ax=plt.subplots()
    img_=img.copy()

    if direction=='horz':
        img_[cut-width/2:cut+width/2,:] = 0
    elif direction=='vert':
        img_[:,cut-width/2:cut+width/2,] = 0
    if title==None:title = 'Img-&-Cut_%s'%cut

    show_img(img_,title=title)


def show_img(img, ax=None,save=False, vmin=None, vmax=None,cmap='spectral',fontsize=24,
             axis_on=True,title_on=True,xlabel=None,ylabel=None,aspect=None,
        title='img_', show=True, logs=False,outDir=None, sizex=9,sizey=9,
        ylim=None, xlim=None, xticks=True, yticks=True,extent=None):
             #rescale=None): #extent=None):

    """show a two-D image"""
    import matplotlib.pyplot as plt
    from matplotlib import cm
    import numpy as np

    if ax is None:fig, ax = plt.subplots(nrows=1)
    if vmin==None:vmin=img.min()
    if vmax==None:vmax=img.max()
    shape = img.shape
    dy,dx= shape

    #ax.set_xlim( 0,dx)
    #ax.set_ylim(0,dy)

    if xlim is not None:ax.set_xlim(xlim)
    if extent is not None:
        x1,x2,y1,y2=extent
        if ylim is  None:ylim=[y2,y1]
        aspect_=determine_aspect(shape, extent)
        #rsx,rsy = [ (x2-x1)/float(dx), (y2-y1)/float(dy)] #rescale x, and rescale y
    else:aspect_=None
    if ylim is  None:ylim=[0,dy]
    #print ylim
    ax.set_ylim([ylim[0], ylim[1]])
    if not logs:
        cax=ax.imshow( img,cmap=cmap, vmin=vmin,vmax=vmax,
            aspect = aspect_,extent=extent)#,interpolation='none')
    if logs:
        img= np.log(img)
        if vmin==None:vmin=img[nonzero(img)].min() #img.min()
        cax=ax.imshow( img,cmap=cmap,vmin=vmin, vmax=vmax,
            aspect = aspect_, extent=extent)#,interpolation='none')
    if aspect is not None:
        im = ax.get_images()
        x1,x2,y1,y2 =  im[0].get_extent()
        if ylim is not None:y1,y2=ylim
        if xlim is not None:x1.x2=xlim
        ax.set_aspect(abs((x2-x1)/(y2-y1))/aspect)
    if title_on:plt.title(  title,fontsize=fontsize )
    if not axis_on:plt.axis('off')
    if xlabel is not None:plt.xlabel(xlabel,fontsize=fontsize)
    else:plt.xlabel('')
    if ylabel is not None:plt.ylabel(ylabel,fontsize=fontsize)
    else:plt.ylabel('')

    if xticks is None:plt.xticks( []  )
    elif xticks is True:plt.setp( ax.get_xticklabels(),fontsize=fontsize*.8, visible=True)#
    else:plt.xticks( xticks,fontsize=fontsize *.8 )
    if yticks is None:plt.yticks( []  )
    elif yticks is True:plt.setp( ax.get_yticklabels(),fontsize=fontsize*.8, visible=True)#
    else:plt.yticks( yticks,fontsize=fontsize *.8 )
    if ax is None:cbar = fig.colorbar(cax, ticks=[vmin,vmax])

    #if ylim is not None:ax.set_ylim(ylim)
    #if xlim is not None:ax.set_xlim(xlim)

    if save:
        if outDir!=None:fp=outDir + title + '_.png'
        else:fp= title + '_.png'
        plt.savefig( fp )
    if show:plt.show()

def trans_data_to_pd(data, label=None,dtype='list'):
    #lists a [ list1, list2...] all the list have the same length
    from numpy import arange,array
    import pandas as pd,sys
    if dtype == 'list':
        data=array(data).T
    elif dtype == 'array':
        data=array(data)
    else:
        print("Wrong data type! Now only support 'list' and 'array' tpye")
    N,M=data.shape
    index =  arange( N )
    if label is None:label=['data%s'%i for i in range(M)]
    #print label
    df = pd.DataFrame( data, index=index, columns= label  )
    return df

def trans_tf_to_td(tf):
    #trans epoch time to datetime
    import pandas as pd
    import numpy as np
    import datetime
    td =np.array( [datetime.datetime.fromtimestamp( tf[i]  )
    for i in range( 0, len(tf)  )] )
    return td

def trans_td_to_tf(td):
    import time,numpy as np
    t = np.array(td)
    tf = np.array( [time.mktime( t[i].timetuple() )
                    for i in range( len(t))] )
    return tf

def plot_pv_values( dict_tv, time, keys=None,title=None, xlabel='Time', ylim_tv=None,
            data= None,  ylim_data=None, img=None,pixsize=0.79,
            aspect=5.0, ylim_img = None):

    import pandas.tools.plotting as ptg
    import matplotlib.dates as md
    #from numpy import arange,array
    import numpy as np
    import datetime
    import matplotlib.pyplot as plt

    if keys is None:
        keys=dict_tv.keys()
    M = len(keys)

    if img is not None:
        N=M+1
        N0=1
        N_img=N0-1
    else:
        N=M
        N0=0
        N_img=-1
    if data is not None:
        N +=1
        N0+=1
        N_data=N_img +1

    sharex=True
    #fig, axs = plt.subplots(N,sharex = True)
    fig,axs=ptg._subplots( N, sharex=sharex,sharey=False, layout=[N,1])
    #tf=time
    axs[0].set_xlim( np.min(time), np.max(time) )
    for n in range(N0,N):
        i=n-N0
        d = dict_tv[ keys[i] ]
        if ylim_tv is  None:
            ymean = d[ keys[i] ].mean()
            ymax,ymin = d[ keys[i] ].max(), d[ keys[i] ].min()
            width = min( [ ymax-ymean, ymean - ymin] )
            ymax_ = ymean + width
            ymin_ = ymean - width
            #print ymean, ymax_, ymin_
            #ylim=[ ymin_ - ymean*.2, ymax_ + ymean*.2]
            ylim = [ymean - width*5, ymean + width *5]
        else:ylim=ylim_tv[i]

        d.plot( x='tf', y = [ keys[i] ],subplots=False, ylim=ylim,
                marker='o',c='b', ls = '--', ax=axs[n])
        yt = axs[n].get_yticks()
        axs[n].set_yticklabels(yt,fontsize=24)
        axs[n].set_ylabel( keys[i] ,fontsize=24)
        ki=keys[i]
        if n!=N-1:
            m=axs[n].get_xticklabels()
            axs[n].set_xlabel('')
            for x in m:x.set_visible(False)

    #dd =np.array( [datetime.datetime.fromtimestamp( tf[i] )
                   #for i in range( 0, len(tf)  )] )
    if img is not None:
        dy,dx = img.shape
        shape=img.shape
        extent= [ np.min(time), np.max(time), dy* pixsize, 0 ]
        show_img(img, aspect =aspect, ax=axs[N_img],
                extent=extent,  title='', ylim=ylim_img,xticks=True)
        m=axs[N_img].get_xticklabels()
        for x in m:x.set_visible(False)
        #axs[N_img].set_xlim( np.min(time), np.max(time) )
        axs[N_img].set_ylabel('size, um',fontsize=24)

    if data is not None:
        data.plot(  x=time, subplots=False,  ylim=ylim_data,
                marker='s',c='g', ls = '--', ax=axs[N_data])  #x= tf,
        m=axs[N_data].get_xticklabels()
        for x in m:x.set_visible(False)
        axs[N_data].set_xlabel('')
        axs[N_data].set_ylabel( data.keys()[0] ,fontsize=24)
        yt = axs[N_data].get_yticks()
        axs[N_data].set_yticklabels(yt,fontsize=24)

    if title is not None:axs[0].set_title(title, fontsize=28)

    axs[n].set_xlim( np.min(time), np.max(time) )

    #print trans_tf_to_td( [np.min(time), np.max(time) ])
    xt = axs[n].get_xticks()
    xticks =  trans_tf_to_td(xt)
    axs[n].set_xlabel(xlabel,fontsize=24)
    axs[n].set_xticklabels([x.strftime('%m-%d %H:%M') for x in xticks],fontsize=24)

    axs[n].set_ylabel( ki ,fontsize=24)

    #axs[0].set_xlim( np.min(time), np.max(time) )


def make_wave_data2( x,y):
    import numpy as np
    x=np.array(x)
    y=np.array(y)
    X = np.zeros( len(x)*2 -1 )
    Y = np.zeros( len(x)*2 -1 )
    X[::2] = x
    X[1::2] = x[1:]
    Y[::2] = y
    Y[1::2] = y[:-1]
    return X,Y

def make_wave_data( x,dtype='x'):
    import numpy as np
    x=np.array(x)
    X = np.zeros( len(x)*2 -1, dtype = object )
    X[::2] = x
    if dtype=='x':
        X[1::2] = x[1:]
    else:
        X[1::2] = x[:-1]
    return X


def get_pvlist_values(PVs,labels):
    '''get current PV-list, PV-Label-list,
       Return a dict
       { label: value}'''
    dict_tv ={}
    N= len(PVs)
    for i in range(N):
        dict_tv[ labels[i]] = get_pv_value( PVs[i] )
    return dict_tv


def get_pv_value( PV ):
    '''get current PV-list, PV-Label-list,
       return PV-value'''
    import epics as ep
    return ep.caget( PV )


def get_archive_pvlist_values(PVs, labels, start_time, end_time,
               limit=100000,interpolation='raw'):
    '''get PV-list, PV-Label-list from archiver from start_time to end_time,
       Return:
       a dict { label:
       pandas.DataFrame with keys as time-str, time-float, label-values
              }
    '''

    from channelarchiver import Archiver
    archiver = Archiver('http://xf11id-ca.cs.nsls2.local/cgi-bin/ArchiveDataServer.cgi')
    archiver_acc = Archiver('http://ca.cs.nsls2.local:8888/cgi-bin/ArchiveDataServer.cgi')

    #archiver.scan_archives()

    dict_tv ={}
    N= len(PVs)
    for i in range(N):
        PV=PVs[i]
        label=labels[i]
        if i==0:scan_archives=True
        else:scan_archives=False

        scan_archives=True
        dict_tv[label]=get_archive_pv_value(PV, label, start_time, end_time,
        scan_archives=scan_archives, limit=limit,
                interpolation=interpolation,archiver=[archiver,archiver_acc] )

    return dict_tv


def get_archive_pv_value(PV, label, start_time, end_time,scan_archives=True,
                   limit=100000,interpolation='raw',archiver=None,):

    '''get PV, PV-Label from archiver from start_time to end_time

    Return:
       pandas.DataFrame with keys as time-str, time-float, label-values
    '''

    import pandas as pd
    import numpy as np
    import time

    if archiver  is None:
        from channelarchiver import Archiver
        archiver = [Archiver('http://xf11id-ca.cs.nsls2.local/cgi-bin/ArchiveDataServer.cgi')]
        #if scan_archives:
    #archiver.scan_archives()
    if label[:3] == 'Acc':arch=   archiver[1]
    else:arch=   archiver[0]
    data = arch.get(PV, start_time, end_time,scan_archives=scan_archives,
                limit=limit,interpolation=interpolation )

    td = np.array(data.times)
    v =np.array(data.values, dtype = float)
    td = make_wave_data( td,dtype='x')
    v = make_wave_data( v,dtype='y')

    tf = trans_td_to_tf(td)
    tv = np.array([ td, tf, v] ).T
    index = np.arange(len(tv) )
    data = tv
    df = pd.DataFrame( data, index=index, columns= ['td','tf',label] )

    print('seraching %s from: %s---to: %s'%(label, start_time, end_time))
    print('Found results:     from: %s  to :%s with  %s points.' %(
        df.td[0],df.td[ len(df.td)-1], len(df.td)))
    return df


def read_scan( sid,fill=True ):
    ''' read data from sid = scan_id'''

    from dataportal import (DataBroker as db,
                        StepScan as ss,
                        StepScan, DataBroker,
                        DataMuxer as dm)
    import datetime

    hdr = db[sid]
    ev=db.fetch_events(hdr,fill=fill)  #, fill=True)
    muxer = dm.from_events( ev )
    data = muxer.to_sparse_dataframe()
    dt = data.time
    #print dt
    #print dt[0]
    #print dt[len(dt)-1]
    #data = list( db.fetch_events(hdr))

    t1 = datetime.datetime.fromtimestamp(dt[0]).strftime('%Y-%m-%d %H:%M:%S')
    t2 = datetime.datetime.fromtimestamp(dt[len(dt)-1]).strftime('%Y-%m-%d %H:%M:%S')

    #t1 = dt[0].strftime('%Y-%m-%d %H:%M:%S')
    #t2 = dt[len(dt)-1].strftime('%Y-%m-%d %H:%M:%S')

    print('the first scan time is:   %s'%t1)
    print('the last scan time  is:   %s'%t2)
    start_time=t1
    end_time =t2
    return data, start_time, end_time


#w,t = get_waterfall( data,direction='vert', cuts=[1000],firstim=0, lastim=1200,detector="xray_cam_img_image_lightfield",)

#show_img(w[1000].T,aspect=1.2, ylim=[700,1400])

def get_waterfall( data,direction='horz', cuts=[1,2,3],firstim=0,
                lastim=1,detector="xray_cam_img_image_lightfield",):

    import numpy as np

    imcount=firstim
    #im_time=[]
    waterfall={}
    notime = lastim +1-firstim
    in_time = data.time
    imgs = data[detector]
    for n in range(0,notime):
        t= firstim + n
        #im_time.append( data[t].time )
        dat=imgs[n]
        #print t, data[t].time
        #dat = data[t].data[detector]
        for i in cuts:
            if t==0:waterfall[i] = []
            if direction=='horz':
                cdat=dat[i,:];
            elif direction=='vert':
                cdat=dat[:,i];
        waterfall[i].append( cdat )
    for i in cuts:waterfall[i]=np.array(waterfall[i])

    return waterfall, in_time


def get_img_waterfall( imgs,direction='horz', cuts=[1,2,3],firstim=0, lastim=1,
                       ):

    import numpy as np

    imcount=firstim
    #im_time=[]
    waterfall={}
    notime = lastim +1-firstim
    #in_time = data.time
    #imgs = data[detector]
    for n in range(0,notime):
        t= firstim + n
        #im_time.append( data[t].time )
        dat=imgs[n]
        #print t, data[t].time
        #dat = data[t].data[detector]
        for i in cuts:
            if t==0:waterfall[i] = []
            if direction=='horz':
                cdat=dat[i,:];
            elif direction=='vert':
                cdat=dat[:,i];
        waterfall[i].append( cdat )
    for i in cuts:waterfall[i]=np.array(waterfall[i])

    return waterfall


def line_focus(filename,imgs=None,times=None, direction='horz',
               cuts=[1,2,3], firstim=0, lastim=1,
               ):
    """
    LW 04/06/2015

    function to analyse line focus data for e.g. stability
    """
    from PIL import Image
    import Tkinter, tkFileDialog
    from matplotlib import pyplot as plt
    import numpy as np
    from scipy.optimize import curve_fit
    #import time
    import os.path, time
    import pandas as pd
    if filename=='ia':
        if filename=='ia':          # open file dialog
            root = Tkinter.Tk()
            root.withdraw()
            file_path = tkFileDialog.askopenfilename()
        else: file_path=filename
        print(file_path)

        # extract filename and first image number:
        a=file_path.split('/')
        fn=a[len(a)-1].split('.')[0]  # file name
        b=fn.split('.')
        #fe=b[len(b)-1]  # file ending, e.g .tiff
        fe=a[len(a)-1].split('.')[1]
        c=fn.split('_')
        fnum=c[len(c)-1].split('.')[0]    # file number
        nlen=len(fnum)  # length of number format
        fp=file_path.split(fn)[0]
        fn=fn.split('_'+fnum)[0]
        print('file path: ',fp)
        print('file name: ',fn)
        print('file ending: ',fe)
        print('file number: ',fnum)

        if filename=='ia':      # get first image number from filebame in interactive mode
            firstim=int(fnum)
            #lastim=30     # maximum number of images to search in autmatic mode

        imcount=firstim
        #print 'imcount: ',imcount
        #print 'firstim: ',firstim
        get_out=0
        image_data=0
        baseline=[]
        amp=[]
        center=[]
        width=[]
        im_time=[]
        while imcount <= lastim:
            ic=str(imcount)
            while len(ic) < nlen:         # create file number
                ic='0'+ic
                    #ic=str(ic)
                    #print ic
            #assemble the image name to be loaded:
            file_to_open=fp+fn+'_'+ic+'.'+fe

            try:
                im = Image.open(file_to_open);
                (mode, ino, dev, nlink, uid, gid, size, atime, mtime, ctime) = os.stat(file_to_open)
                im_time.append(time.ctime(mtime))
                #im.show()
                #im.show()
                print('loaded file: ',file_to_open)
            except:
                print('waiting for file: '+file_to_open)
                #plt.close(1)    # use the waiting time to plot some intermediate result -> does not work, as python is waiting until the end before printing....
                #plt.figure(1)
                #imgplot= plt.imshow(image_data/(imcount-firstim+1),origin='lower')
                #imgplot.set_cmap('gray')
                time.sleep(5)
                try:
                    im = Image.open(file_to_open);
                    (mode, ino, dev, nlink, uid, gid, size, atime, mtime, ctime) = os.stat(file_to_open)
                    im_time.append(time.ctime(mtime))
                    #im.show()
                except:get_out=1
            if get_out==1:
                return image_data ### THIS NEEDS TO BE CHANGED!!!
                plt.close(1)    # only plot at the end for now: too much overhead to plot every iteration
                plt.figure(1)
                raise CHX_utilities_Exception('no more images found, images summed up to: '+file_to_open)
            # do the analysis
            im=np.array(im)*1.0;
            #print 'shape of image: ',im.shape
            #plt.figure(1)
            #plt.imshow(im)
            try:
              dat=np.sum(im,axis=2); # sum RGB channels for a color image
            except:
              dat=im;
            for i in cuts:
                if direction=='horz':
                    cdat=dat[i,:];
                elif direction=='vert':
                    cdat=dat[:,i];
                    #print cdat;
                else: raise CHX_utilities_Exception('error: cut direction must be either "horz" or "vert".')
                # create initial guess:
                yo=np.mean(cdat)
                A=np.max(cdat)
                xc=np.argmax(cdat)*1.0
                w=8
                p0=[yo,A,xc,w]
                try:
                    coeff,var_matrix = curve_fit(gauss,np.array(xrange(len(cdat))),cdat,p0=p0);
                except: coeff=[0,0,0,0]
                baseline.append(coeff[0])
                amp.append(coeff[1])
                center.append(coeff[2])
                width.append(coeff[3])
            imcount=imcount+1

    if filename=='im':
        imcount=firstim
        im_time = times
        baseline=[]
        amp=[]
        center=[]
        width=[]
        notime = lastim +1-firstim

        for n in range(0,notime):
            t= firstim + n
            #im_time.append( data[t].time )
            im=imgs[t]

            try:
                dat=np.sum(im,axis=2); # sum RGB channels for a color image
            except:
                dat=im;

            for i in cuts:
                if direction=='horz':
                    cdat=dat[i,:];
                elif direction=='vert':
                    cdat=dat[:,i];
                        #print cdat;
                else: raise CHX_utilities_Exception('error: cut direction must be either "horz" or "vert".')
                # create initial guess:
                yo=np.mean(cdat)
                A=np.max(cdat)
                xc=np.argmax(cdat)*1.0
                w=8
                p0=[yo,A,xc,w]
                try:
                    coeff,var_matrix = curve_fit(gauss,
                                    np.array(xrange(len(cdat))),cdat,p0=p0);
                except: coeff=[0,0,0,0]
                baseline.append(coeff[0])
                amp.append(coeff[1])
                center.append(coeff[2])
                width.append(coeff[3])

    if filename=='wt':
        imcount=firstim
        im_time = times
        baseline=[]
        amp=[]
        center=[]
        width=[]
        notime = lastim +1-firstim

        for n in range(0,notime):
            t= firstim + n
            #im_time.append( data[t].time )
            im=imgs[ :,t]
            try:
                dat=np.sum(im,axis=2); # sum RGB channels for a color image
            except:
                dat=im;

            for i in cuts:
                cdat=im

                #else: raise CHX_utilities_Exception('error: cut direction must be either "horz" or "vert".')
                # create initial guess:
                yo=np.mean(cdat)
                A=np.max(cdat)
                xc=np.argmax(cdat)*1.0
                w=8
                p0=[yo,A,xc,w]
                try:
                    coeff,var_matrix = curve_fit(gauss,
                                    np.array( xrange(  len(cdat)  )),cdat,p0=p0);
                except: coeff=[0,0,0,0]
                baseline.append(coeff[0])
                amp.append(coeff[1])
                center.append(coeff[2])
                width.append(coeff[3])

    #df_cen = trans_data_to_pd([center],label=['center'],dtype='list')

    return np.array(baseline),np.array(amp),np.array(center),np.array(width),np.array(im_time);


def show_fit( data, para, func = gauss,title=None,ax=None):
    import matplotlib.pyplot as plt
    import numpy as np
    x,y = data
    x = np.array(x)
    y=np.array(y)
    #print 'here'
    x0 = np.linspace( x.min(),x.max(),2000)
    fit =  func(x0,*para)
    if ax is None:fig,ax=plt.subplots(  )
    if title is not None:ax.set_title( title )
    ax.plot( x,y,'bo')
    ax.plot(x0,fit,'r',ls='-')
    #plt.show()


def show_focus_line_fit( waterfall, para, lists=None,title=None):
    import matplotlib.pyplot as plt
    L, t = waterfall.shape
    if lists is None:
        lists=[]
        #lists.append( [t/4,t/2,3*t/4])
        lists.append( [t/5,2*t/5,3*t/5,4*t/5])
    #print lists
    fig,ax=plt.subplots( len(lists[0]) )
    #if title is not None:plt.title( title )
    #print title
    for i,l in enumerate(lists[0]):
        im=waterfall[:,l]
        p = para[i:i+1].values[0]
        #print i,l
        show_fit( [range(L), im], p,ax=ax[i], title =title+'_time@_%s'%l)


def test_fft():
    import matplotlib.pyplot as plt
    import numpy as np
    fig, ax = plt.subplots(2)
    t=np.linspace(0,20,1000)
    y= .5*np.sin(2*np.pi*5*t) + 2*np.sin(2*np.pi*t) + 1*np.sin(2*np.pi*14*t)
    ax[0].plot(t,y,'ro',ls='--')

    fx,fy=get_fft( t,y)
    ax[1].plot(fx,fy,'go',ls='--')
    plt.show()


def get_fft( t,y ):
    import scipy
    import scipy.fftpack
    import numpy as np

    L = len(t)
    tp = np.array(t)
    yp=np.array( y)
    dx = (t[L-1] - t[0])/float(L)
    tm= 2 * L * dx
    xs =  1/ tm
    #print ts,xs, tm
    #fx = xs * np.arange( L  )

    y= np.concatenate( ( y, np.zeros( L )))

    FFT = abs(scipy.fft(y))[:L]
    freqs = scipy.fftpack.fftfreq(y.size, dx)[:L]#
    #FFT    = np.fft.fft(y)
    #freqs = np.fft.fftfreq(len(y), dx )##
    # Find the peak in the coefficients
    idx = np.argmax(np.abs(FFT))
    freq = freqs[idx]
    freq_in_hertz = abs(freq ) #* dx)
    print('The maximum frequency is:  %s'%(freq_in_hertz ))##
    #freqs = trans_data_to_pd([freqs], label=['freqs'],dtype='list')##
    #20*scipy.log10(FFT)
    #FFT = trans_data_to_pd([  20*scipy.log10(FFT)  ], label=label,dtype='list')
    #FFT = trans_data_to_pd([  (FFT)  ], label=label,dtype='list')##
    return np.array(freqs), np.array( FFT   )


def plot_line_focus2(df_res,cuts,pix=0.79,title=None,
            times=None, ROI=None, FFT=True):
    from matplotlib import pyplot as plt
    import numpy as np
    import time,datetime
    import pandas.tools.plotting as ptg
    import matplotlib.dates as md
    import pandas as pd
    if ROI is not None:
        x1,x2=ROI
    else:
        x1,x2=0,len(times)

    times_ = times[x1:x2]
    df_res_ = df_res.loc[x1:x2-1]
    df_res_.index = range( 0, x2-x1 )

    N=len(times_)
    keys = df_res.keys()
    M = len(keys)
    td=times_
    tf =  trans_td_to_tf( td )

    if title==None:title='LF_'
    fig,axs=ptg._subplots( M, sharex=True,sharey=False,
                           layout=[M,1])
    axs[0].set_title(title)
    for n in range(M):
        df_res_.plot( x=td, y = [ keys[n] ],subplots=False,
                marker='o',c='b', ls = '--', ax=axs[n])

    #do analysis _center

    fig2,axs2=ptg._subplots( 2, sharex=True,sharey=True,
                           layout=[2,1])

    yc=df_res_.center
    #_df_res_ = df_res
    df_res_['cen_pos'] = pix*(yc - yc.mean())
    df_res_.plot( x=td, y = [ 'cen_pos' ],subplots=False,
                 title=title+'center position',
            marker='o',c='k', ls = '--', ax=axs2[0],label=str(cuts[0]))

    x=xrange(N)
    y= yc
    pol=np.polyfit(x,y,20) #20)
    polfit=np.polyval(pol,x)
    yfitted=y-polfit

    df_res_['cen_fitted'] = (yfitted-np.mean(yfitted))*pix
    #print dd.shape,df_res.cen_fitted.shape
    i=0
    mean = (  polfit- np.mean(polfit) )*pix
    std =  mean.std()
    #label=str(str(cuts[i]) +' PV: '+str(round(max(polfit)-min(polfit),2))+ 'um   rms: +/-'+str(  round(std0,2))+'um' )

    label=str(str(cuts[i])
        +' Center: '+str( round(max(polfit)-min(polfit),2) )+ 'um   rms: +/-'
        +str(  round(std,2))+'um' )

    axs2[0].plot( td, np.zeros( len(td) ) + std, c='red',ls='--',lw=4)
    axs2[0].plot( td, np.zeros( len(td) )- std,  c='red',ls='--',lw=4)

    axs2[0].text(  .3,.9,label,transform=axs2[0].transAxes)
    df_res_['cen_polyfit'] = (polfit-np.mean(polfit))*pix

    df_res_.plot( x=td, y = [ 'cen_polyfit' ],subplots=False,
                 title=title+'Drift of center position',
            marker='',c='b', ls = '-', lw=3, ax=axs2[0],label=label)

    #label=str( str(cuts[i]) +' PV: '+ str(round(max(yfitted)-min(yfitted),2)) +  'um   rms: +/-'+str(round(( (yfitted-np.mean(yfitted))*pix ).std(),2))+'um'  )

    mean = (  yfitted- np.mean(yfitted) )*pix
    std =  mean.std()
    label=str(str(cuts[i])
        +'PV: '+str( round(max(polfit)-min(polfit),2) )+ 'um   rms: +/-'
        +str(  round(std,2))+'um' )

    axs2[1].text(  .3,.9,label,transform=axs2[1].transAxes)
    axs2[1].plot( td, np.zeros( len(td) ) + std, c='red',ls='--',lw=4)
    axs2[1].plot( td, np.zeros( len(td) ) - std, c='red',ls='--',lw=4)

    df_res_.plot( x=td, y = [ 'cen_fitted' ],subplots=False,
                title=title+'drift corrected vibrations',
            marker='+',c='b', ls = '-', ax=axs2[1],label=label)

    axs2[1].set_xlabel('time')
    xfmt = md.DateFormatter('%m-%d %H:%M')
    axs2[1].xaxis.set_major_formatter(xfmt)

    axs[n].set_xlabel('time')
    xfmt = md.DateFormatter('%m-%d %H:%M')
    axs[n].xaxis.set_major_formatter(xfmt)

    if FFT:
        x,y0 =  get_fft( tf, df_res_['cen_pos'])
        x,y1 = get_fft( tf, df_res_['cen_fitted'])
        x,y2 = get_fft( tf, df_res_['cen_polyfit'])
        df_y = trans_data_to_pd(
            [x,y0,y1,y2],label=['freqs, Hz','fft_cen_pos',
                    'fft_cen_fitted','fft_cen_polyfit'
                        ],dtype='list')

        df_y.plot(subplots=True,x='freqs, Hz',logx=True, y=['fft_cen_pos',
                    'fft_cen_fitted','fft_cen_polyfit'], title=title + 'FFT')

    #return df_y


def plot_line_focus(res,cuts,pix=1,epochtime='false', times=None):
    """
    LW 04/06/2015
    function to plot results of line focus analysis
    """
    from matplotlib import pyplot as plt
    import numpy as np
    import time
    from dateutil import parser
    base=np.array(res[0])
    amp=np.array(res[1])
    center=np.array(res[2])
    width=np.array(res[3])
    epoch=[]
    plt.close()
    if times is None:
        for i in xrange(len(res[4][:])):
            dt =parser.parse(res[4][i])
            epoch.append(dt.strftime('%s'))
    else:epoch=times

    epoch=np.array(epoch)
    plt.figure(1)           #### plot and anlyze center position
    for i in xrange(len(cuts)):
        y=center[i::len(cuts)]
        #print y
        if epochtime=='false':
           x=xrange(len(y))
           x_str='image #'
        else:
          x=epoch
          x_str='epoch [s]'
        plt.plot(x,(y-np.mean(y))*pix,'+',label=str(cuts[i]))
    plt.legend(loc='best')
    plt.title('center position')
    plt.ylabel('x_0 - <x_0> [um]')
    plt.xlabel(x_str)

    plt.figure(5)
    plt.figure(6)           ##### plot and analyze center short time vibrations
    for i in xrange(len(cuts)):
        y=np.array(center[i::len(cuts)])
        #print y
        if epochtime=='false':
            x=xrange(len(y))
            x_str='image #'
        else:
          #x=epoch
          x=xrange(len(y))
          x_str='epoch [s]'

    pol=np.polyfit(x,y,20)
    polfit=np.polyval(pol,x)
    yfitted=y-polfit
    plt.figure(5)

    plt.plot(epoch,(yfitted-np.mean(yfitted))*pix,'+-',
    label=str(str(cuts[i])+' PV: '+str(round(max(yfitted)-min(yfitted),2))+
    'um   rms: +/-'+str(round(((yfitted-np.mean(yfitted))*pix).std(),2))+'um'))
    plt.figure(6)
    plt.plot(epoch,(polfit-np.mean(polfit))*pix,
        label=str(str(cuts[i])+' PV: '+str(round(max(polfit)-min(polfit),2))+
    'um   rms: +/-'+str(round(((polfit-np.mean(polfit))*pix).std(),2))+'um'))

    plt.figure(5)
    plt.legend(loc='best')
    plt.title('drift corrected vibrations')
    plt.ylabel('x_0 - <x_0> [um]')
    plt.xlabel(x_str)
    plt.figure(6)
    plt.legend(loc='best')
    plt.title('Drift of center position')
    plt.ylabel('x_0 - <x_0> [um]')
    plt.xlabel(x_str)

    plt.figure(6)           #### plot and analyze long time drifts

    plt.figure(2)   ############ plot and analyze fwhm
    for i in xrange(len(cuts)):
        y=width[i::len(cuts)]
        #print y
        y=2.35*y*pix
        plt.plot(x,y,'+',label=str(cuts[i]))
        print('minimum focus: ',y.min(),'um   lowest 10%: ',y[x<(y.min()+.1*y.min())].mean(),'+/-',y[x<(y.min()+.1*y.min())].std(),'um @pixel cut: ',cuts[i])
    plt.legend(loc='best')
    plt.title('FWHM')
    plt.ylabel('Focus FWHM [um]')
    plt.xlabel(x_str)

    plt.figure(3)
    for i in xrange(len(cuts)):
        y=amp[i::len(cuts)]
        #print y
        #x=xrange(len(y))
        plt.plot(x,y,'+',label=str(cuts[i]))
    plt.legend(loc='best')
    plt.title('Amplitude')
    plt.ylabel('Amplitude [arb.u.]')
    plt.xlabel(x_str)

    plt.figure(4)
    for i in xrange(len(cuts)):
        y=base[i::len(cuts)]
        #print y
        #x=xrange(len(y))
        plt.plot(x,y,'+',label=str(cuts[i]))
    plt.legend(loc='best')
    plt.title('Baseline')
    plt.ylabel('y_0 [arb.u.]')
    plt.xlabel(x_str)


def read_camonitor(filename,epoch='true'):
    import Tkinter, tkFileDialog
    import time
    import numpy as np
    import datetime
    if filename=='ia':          # open file dialog
        root = Tkinter.Tk()
        root.withdraw()
        file_path = tkFileDialog.askopenfilename()
    else: file_path=filename
    print('loading camonitor data from: ', file_path)
    with open(file_path,"r") as cafile:
     cadata=cafile.readlines()
    #datime=[]
    #datdate=[]
    value_camon=[]
    epoch_camon=[]
    for i in xrange(len(cadata)):
        line=cadata[i].split()
        value_camon.append(float(line[len(line)-1]))
        datdate=line[len(line)-3]
        datime=line[len(line)-2]
        #print datdate.split()[0]
        #print datime
        #print int(datdate.split('-')[0]),int(datdate.split('-')[1]),int(datdate.split('-')[2]),int(datime.split(':')[0]),int(datime.split(':')[1]),int(datime.split(':')[2].split('.')[0]),int(datime.split(':')[2].split('.')[1])
        epoch_camon.append(float(datetime.datetime(int(datdate.split('-')[0]),int(datdate.split('-')[1]),int(datdate.split('-')[2]),int(datime.split(':')[0]),int(datime.split(':')[1]),int(datime.split(':')[2].split('.')[0]),int(datime.split(':')[2].split('.')[1])).strftime('%s')))
    return np.array(epoch_camon),np.array(value_camon)


def knife_edge(filename, direction='horz', cuts=[1,2,3], firstim=0, lastim=1,pix=1):
    """
    LW 04/06/2015
    function to analyse knife edge image data for e.g. stability
    """
    from PIL import Image
    import Tkinter, tkFileDialog
    from matplotlib import pyplot as plt
    import numpy as np
    from scipy.optimize import curve_fit
    import time
    import scipy as sp

    if filename=='ia':          # open file dialog
        root = Tkinter.Tk()
        root.withdraw()
        file_path = tkFileDialog.askopenfilename()
    else: file_path=filename
    print(file_path)

    # extract filename and first image number:
    a=file_path.split('/')
    fn=a[len(a)-1].split('.')[0]  # file name
    b=fn.split('.')
    #fe=b[len(b)-1]  # file ending, e.g .tiff
    fe=a[len(a)-1].split('.')[1]
    c=fn.split('_')
    fnum=c[len(c)-1].split('.')[0]    # file number
    nlen=len(fnum)  # length of number format
    fp=file_path.split(fn)[0]
    fn=fn.split('_'+fnum)[0]
    print('file path: ',fp)
    print('file name: ',fn)
    print('file ending: ',fe)
    print('file number: ',fnum)

    if filename=='ia':      # get first image number from filebame in interactive mode
        firstim=int(fnum)
        #lastim=30     # maximum number of images to search in autmatic mode

    imcount=firstim
    #print 'imcount: ',imcount
    #print 'firstim: ',firstim
    get_out=0
    image_data=0
    center=[]
    while imcount <= lastim:
        ic=str(imcount)
        while len(ic) < nlen:         # create file number
                ic='0'+ic
                #ic=str(ic)
                #print ic
        #assemble the image name to be loaded:
        file_to_open=fp+fn+'_'+ic+'.'+fe

        try:
            im = Image.open(file_to_open);
            print('loaded file: ',file_to_open)
        except:
            print('waiting for file: '+file_to_open)
            #plt.close(1)    # use the waiting time to plot some intermediate result -> does not work, as python is waiting until the end before printing....
            #plt.figure(1)
            #imgplot= plt.imshow(image_data/(imcount-firstim+1),origin='lower')
            #imgplot.set_cmap('gray')
            time.sleep(5)
            try:
                im = Image.open(file_to_open);
                #im.show()
            except:get_out=1
        if get_out==1:
            return image_data ### THIS NEEDS TO BE CHANGED!!!
            plt.close(1)    # only plot at the end for now: too much overhead to plot every iteration
            plt.figure(1)
            raise CHX_utilities_Exception('no more images found, images summed up to: '+file_to_open)
        im=np.array(im)*1.0;
        try:
          dat=np.sum(im,axis=2); # sum RGB channels for a color image
        except:
          dat=im;
    if imcount==firstim:
      plt.close(4)
      plt.figure(4)
      plt.imshow(dat,origin='lower')
      for i in cuts:
        plt.plot([i,i],[0,2048],'r--')          #2048 is a fixed value for now, need to get dimension from data in the future
      plt.tight_layout()
      plt.xlabel('pixel')
      plt.ylabel('pixel')
      plt.title(file_to_open)
    for i in cuts:
            if direction=='vert':
                cdat=dat[i,:];
            elif direction=='horz':
                cdat=dat[:,i];
               #print cdat;
            else: raise CHX_utilities_Exception('error: cut direction must be either "horz" or "vert".')
            ######### try some fourier filtering before diff:
            fft=sp.fft(cdat)
            bp=fft[:];
            for i in range(len(bp)):
                if i>=100:bp[i]=0
            cdat=sp.ifft(bp);
        ###### try with normalization
        #devdat=np.diff(cdat/max(cdat))
        #######
            devdat=np.diff(cdat)
            ##########
            if imcount==firstim:
                plt.close(9)
                plt.figure(9)
                plt.plot(xrange(len(devdat)),devdat/max(devdat)*max(cdat))
                plt.plot(xrange(len(cdat)),cdat)

        #print devdat
        #np.argmax(devdat)
        #center.append(np.argmax(devdat))
        #### try helping with using the correct maximum (allows to use more Fourier components)
    center.append(np.argmax(devdat[900:1500]))
        ####
    imcount=imcount+1
    #plotting section:
    #plt.figure(3)
    #plt.plot(xrange(len(devdat)),devdat)
    #plt.figure(2)
    #plt.plot(xrange(len(cdat)),cdat,'o-')
    plt.close(51)
    plt.figure(51)
    for i in xrange(len(cuts)):
        y=center[i::len(cuts)]
        #print y
        x=xrange(len(y))
        dy=(y-np.mean(y))*pix;
        plt.plot(x,dy,'+',label=str(cuts[i]))
    plt.legend(loc='best')
    plt.title('knife edge position')
    plt.ylabel('x_0 - <x_0> [um]')
    plt.xlabel('image #')
    print('knife edge position: ',dy.mean(),'um +/- ',dy.std(),'um')
    plt.close(52)
    plt.figure(52)
    for i in xrange(len(cuts)):
        y=center[i::len(cuts)]
        #print y
        x=xrange(len(y))
        dy=y;
        plt.plot(x,dy,'+',label=str(cuts[i]))
    plt.legend(loc='best')
    plt.title('knife edge position')
    plt.ylabel('knife edge position [pixel]')
    plt.xlabel('image #')
    return center

#copied from 99-bluesky....didn't figure out how to import...
def detselect(detector_object, suffix="_stats1_total"):
    """Switch the active detector and set some internal state"""
    gs.DETS =[detector_object]
    gs.PLOT_Y = detector_object.name + suffix
    gs.TABLE_COLS = [gs.PLOT_Y]

def get_ID_calibration(gapstart, gapstop, xray_eye1=None, gapstep=.2, gapoff=0, sl=300):
    """
    by LW 04/20/2015
    function to automatically take a ID calibration curve_fit
    calling sequence: get_ID_calibration(gapstart,gapstop,gapstep=.2,gapoff=0,sl=300)
        gapstart: minimum gap used in calibration (if <5.2, value will be set to 5.2)
        gapstop: maximum gap used in calibration
        gapstep: size of steps between two gap points
        gapoff: offset applied to calculation gap vs. energy from xfuncs.get_Es()
        sl: sleep between two gap points (to avoid overheating the DCM Bragg motor)
    writes outputfile with fitted value for the center of the Bragg scan to:  '/home/xf11id/Repos/chxtools/chxtools/X-ray_database/
        changes 03/18/2016: made compatible with python V3 and latest version of bluesky (working on it!!!)
    """
    import numpy as np
    import xfuncs as xf
    #from dataportal import DataBroker as db, StepScan as ss, DataMuxer as dm
    import time
    from epics import caput, caget
    from matplotlib import pyplot as plt
    from scipy.optimize import curve_fit

    ### just testing passing bluesky objects (motors, detectors) into function
    gs.DETS =[xray_eye1]
    #detselect(xray_eye1)
    dscan(diff.xh,gapstart,gapstop,3)

    ### end of testing section  #################

    gaps=np.arange(gapstart,gapstop,gapstep)+gapoff
    print('ID calibration will contain the following gaps [mm]: ',gaps)
    if caget('XF:11IDA-OP{Mono:DCM-Ax:X}Pos-Sts') == 1:
        xtal='Si111cryo'
    elif caget('XF:11IDA-OP{Mono:DCM-Ax:X}Pos-Sts') == 2:
        xtal='Si220cryo'
    else: raise CHX_utilities_Exception('error: trying to do ID gap calibration with no crystal in the beam')
    print('using ',xtal,' for ID gap calibration')
    # create file for writing calibration data:
    fn='id_CHX_IVU20_'+str(time.strftime("%m"))+str(time.strftime("%d"))+str(time.strftime("%Y"))+'.dat'
    dat_file = Path(rs_fn('chxtools', 'X-ray_database')) / Path(fn)
    try:
        with dat_file.open('w') as outFile:
            outFile.write('% data from measurements {}\n'.format(time.strftime("%D")))
            outFile.write('% K column is a placeholder! \n')
            outFile.write('% ID gap [mm]     K      E_1 [keV] \n')
            print('successfully created output file: {}'.format(dat_file))
    except:
        raise CHX_utilities_Exception('error: could not create output file')

    ### do the scanning and data fitting, file writing,....
    center=[]
    E1=[]
    realgap=[]
    #detselect(xray_eye1)
    for i in gaps:
        if i>= 5.2:
            B_guess=-1.0*xf.get_Bragg(xtal,xf.get_Es(i-gapoff,5)[1])[0]
        else:
         i=5.2
         B_guess=-1.0*xf.get_Bragg(xtal,xf.get_Es(i,5)[1])[0]
    print('initial guess: Bragg= ',B_guess,' deg.   ID gap = ',i,' mm')
    if xf.get_Es(i,5)[1] < 9.5 and round(caget('XF:11IDA-OP{Mir:HDM-Ax:Y}Mtr.VAL'),1) != -7.5:
        caput('XF:11IDA-OP{Mir:HDM-Ax:Y}Mtr.VAL',-7.5)  # use HDM Si stripe
        time.sleep(20)
    elif xf.get_Es(i,5)[1] >= 9.5 and round(caget('XF:11IDA-OP{Mir:HDM-Ax:Y}Mtr.VAL'),1) != 7.5:
        caput('XF:11IDA-OP{Mir:HDM-Ax:Y}Mtr.VAL',7.5)   # use HDM Rh stripe
        time.sleep(20)
    if round(caget('XF:11IDA-BI{Foil:Bpm-Ax:Y}Mtr'),1) != 0.0:
        caput('XF:11IDA-BI{Foil:Bpm-Ax:Y}Mtr',0.0)
        time.sleep(30)
    else: pass
    print('moving DCM Bragg angle to: ',B_guess,' deg and ID gap to ',i,' mm')
    #dcm.b.timeout=1200	#make sure dcm motions don't timeout...
    #dcm.en.timeout=1200
    #mov(dcm.b,B_guess)
    #mov(ivu_gap,i)
    print('hurray, made it up to here!')
#  ascan(dcm.b,float(B_guess-.4),float(B_guess+.4),60)   # do the Bragg scan
    header = db[-1]  # retrive the data (first data point is often "wrong", so don't use
    data = get_table(header)
    B = data.dcm_b[2:]
    intdat = data.xray_eye1_stats1_total[2:]
    B=np.array(B)
    intdat=np.array(intdat)
    B=np.array(ss[-1].dcm_b)[2:]  # retrive the data (first data point is often "wrong", so don't use
    intdat=np.array(ss[-1].bpm_cam_stats_total1)[2:]
    A=np.max(intdat)          # initial parameter guess and fitting
    xc=B[np.argmax(intdat)]
    w=.2
    yo=mean(intdat)
    p0=[yo,A,xc,w]
    print('initial guess for fitting: ',p0)
    try:
        coeff,var_matrix = curve_fit(gauss,B,intdat,p0=p0)
        center.append(coeff[2])
        E1.append(xf.get_EBragg(xtal,-coeff[2])/5.0)
        realgap.append(caget('SR:C11-ID:G1{IVU20:1-LEnc}Gap'))
#   # append data file by i, 1 & xf.get_EBragg(xtal,-coeff[2]/5.0):
        with dat_file.open('a') as myfile:
            myfile.write(str(caget('SR:C11-ID:G1{IVU20:1-LEnc}Gap'))+'    1.0 '+str(float(xf.get_EBragg(xtal,-coeff[2])/5.0))+'\n')
        print('added data point: ',caget('SR:C11-ID:G1{IVU20:1-LEnc}Gap'),' ',1.0,'     ',str(float(xf.get_EBragg(xtal,-coeff[2])/5.0)))
    except: print('could not evaluate data point for ID gap = ',i,' mm...data point skipped!')
    time.sleep(sl)
    plt.close(234)
    plt.figure(234)
    plt.plot(E1,gaps,'ro-')
    plt.xlabel('E_1 [keV]')
    plt.ylabel('ID gap [mm]')
    plt.title('ID gap calibration in file: {}'.format(dat_file), size=12)
    plt.grid()


class CHX_utilities_Exception(Exception):
    pass
    """
    by LW 03/19/2015
    class to raise xfuncs specific exceptions
    """

