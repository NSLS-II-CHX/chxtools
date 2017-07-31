from epics import caput, caget
import numpy as np

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from datetime import datetime
import time as ti

pv_ca = 'XF:11IDB-BI{XBPM:02}FA-A'
pv_cb = 'XF:11IDB-BI{XBPM:02}FA-B'
pv_cc = 'XF:11IDB-BI{XBPM:02}FA-C'
pv_cd = 'XF:11IDB-BI{XBPM:02}FA-D'
pv_px = 'XF:11IDB-BI{XBPM:02}FA-X'
pv_py = 'XF:11IDB-BI{XBPM:02}FA-Y'
pv_sumi='XF:11IDB-BI{XBPM:02}FA-S'

meta_pvdict = dict(

pv_Px = 'XF:11IDB-BI{XBPM:02}Fdbk:AKp-SP',
pv_Ix = 'XF:11IDB-BI{XBPM:02}Fdbk:AKi-SP',
pv_Dx = 'XF:11IDB-BI{XBPM:02}Fdbk:AKd-SP',
pv_fdx= 'XF:11IDB-BI{XBPM:02}Fdbk:AEn-SP',
pv_sklx= 'XF:11IDB-BI{XBPM:02}Fdbk:ACtrlScaleFactor-SP',
pv_Kx = 'XF:11IDB-BI{XBPM:02}Pos:Kx-SP',

pv_Py = 'XF:11IDB-BI{XBPM:02}Fdbk:BKp-SP',  
pv_Iy = 'XF:11IDB-BI{XBPM:02}Fdbk:BKi-SP',
pv_Dy = 'XF:11IDB-BI{XBPM:02}Fdbk:BKd-SP',
pv_fdy= 'XF:11IDB-BI{XBPM:02}Fdbk:BEn-SP',
pv_skly= 'XF:11IDB-BI{XBPM:02}Fdbk:BCtrlScaleFactor-SP',
pv_Ky = 'XF:11IDB-BI{XBPM:02}Pos:Ky-SP',

pv_fdHz = 'XF:11IDB-BI{XBPM:02}Fdbk:delT-I',
    )
#'pv_Px', 'pv_Ix', 'pv_Dx', 'pv_fdx', 'pv_sklx', 'pv_Kx', 
#'pv_Py', 'pv_Iy', 'pv_Dy', 'pv_fdy', 'pv_skly', 'pv_Ky', 
#'pv_fdHz'

 


pv_num_sam = 'XF:11IDB-BI{XBPM:02}Trig:NumSamples-SP'
pv_trig = 'XF:11IDB-BI{XBPM:02}FaSoftTrig-SP'

data_path = '/XF11ID/analysis/Commissioning/BPM_Stability/Data/'
res_path = '/XF11ID/analysis/Commissioning/BPM_Stability/Result/'

#filename


def get_fft( t,y ): 
    import numpy as np
    import scipy.fftpack
    L = len(t)
    tp = np.array(t, dtype = float)
    yp=np.array( y, dtype = float)
    ts = (t[L-1] - t[0])/float(L)    
    tm= 2 * L * ts
    xs =  1/ tm 
    ps = np.abs( np.fft.fft(  yp ) )
    time_step =   ts    
    freqs = np.fft.fftfreq(yp.size, time_step)
    idx = np.argsort(freqs)
    F,Y = freqs[idx][int(L/2+1):], ps[idx][int(L/2+1):]
    ##  Find the peak in the coefficients
    idxm = np.argmax( Y )
    freq = F[idxm]
    freq_in_hertz = freq  
    print ('The maximum frequency is:  %s'%(freq_in_hertz))##

    return F,Y  


def plot_current(t, ca,cb,cc,cd,sumi, res_path, filename, metadata  ):
    fig = plt.figure(figsize=(8,12)) 
    plt.axis('off')
    plt.title(filename + '_' + str( round(1/metadata['pv_fdHz'],1) ) + ' Hz')
    #plt.title(filename)
    colms =  ['A', 'B','C','D', 'Sum']
    for n,i in   enumerate(  [ ca,cb,cc,cd,sumi  ]  ):
        y = i
        ax = fig.add_subplot( 5,1,n+1 ) 
        ax.plot(t, y, '--o', label="Col-- "+ str( colms[n] ))
        ax.set_xlabel("Time, (s)")
        ax.set_ylabel(str( colms[n] ))
        #ax.set_ylim(0, 6E2)
        #ax.set_xlim(0, 1)
        ax.legend( loc='best', fontsize = 16) 
    #try filename = filename.rstrip('.txt') 
    plt.savefig( res_path + filename  + '-current.png')    
    plt.show()



def plot_posxy_fft(t,posx,posy,res_path, filename, metadata, frequency_range=[0,500]  ):   
            
    fig = plt.figure(figsize=(12,10))
    plt.axis('off')
    plt.title(filename + '_' + str( round(1/metadata['pv_fdHz'],1) ) + ' Hz')
    ax = fig.add_subplot( 221 )
    #ax.set_title(filename )
    y = posx        
    ym = y.mean()

    std = (y-ym).std()
    ax.plot(   t, y - ym, '--o', label="X-<X>"   )
    #ax.plot(t, std+ np.zeros_like(t), '--o', label="%s +/-"%round(ym,2) + str( round(std,3) ) + ' um' ) 
    ax.plot(t, std+ np.zeros_like(t), '--o', label= 'rms-'+ str( round(std,3) ) + ' um' )     
    ax.plot(t, -std+ np.zeros_like(t), '--o')#, label="rms-- "+ str( colms[i] ))

    ax.set_xlabel("Time, (s)")
    ax.set_ylabel("X")
    ax.legend( loc='best', fontsize = 16)

    ax = fig.add_subplot( 223 )    
     
    y = posy
    ym=y.mean()

    std = (y-ym).std()
    ax.plot(t, y - ym, '--o', label="Y-<Y>")
   # ax.plot(t,  std + np.zeros_like(t), '--o', label="%s +/-"%round(ym,2) +  str( round(std,3)) + ' um' )
    ax.plot(t,  std + np.zeros_like(t), '--o', label="rms-" +  str( round(std,3)) + ' um' )
    ax.plot(t, -std+ np.zeros_like(t), '--o')#, label="rms-- "+ str( colms[i] ))


    ax.set_xlabel("Time, (s)")
    ax.set_ylabel("Y")    
    ax.legend( loc='best', fontsize = 16) 
    #try filename = filename.rstrip('.txt')
    
    ax = fig.add_subplot( 222 )
    #ax.set_title(filename )
    yt = posx
    freq,fft =  get_fft( t,yt )
    ax.plot(freq,fft, '--o', label="FFT-posX" )
    ax.set_xlabel("freq, (Hz)")
    ax.set_ylabel("fft_x")
    ax.set_xlim( frequency_range)
    ax.legend( loc='best', fontsize = 16)    

    pp = 'x'
    txt1 = 'P I D: %s  %s  %s, Scaling: %s'%( metadata['pv_P%s'%pp], metadata['pv_I%s'%pp], metadata['pv_D%s'%pp],metadata['pv_skl%s'%pp] )
    ax.text(x = .4, y= 0.7,  s= txt1, fontsize=14, transform=ax.transAxes) 
    txt2 = 'K%s: %s, Feedback: %s'%(pp, metadata['pv_K%s'%pp], metadata['pv_fd%s'%pp] )
    ax.text(x = .4, y= 0.6,  s= txt2, fontsize=14, transform=ax.transAxes) 
    
    ax = fig.add_subplot( 224 )  
    yt = posy
    freq,fft =  get_fft( t,yt )
    ax.plot(freq,fft, '--o', label="FFT-PosY")
    ax.set_xlabel("freq, (Hz)")
    ax.set_ylabel("fft_y")
    ax.set_xlim( frequency_range)
    ax.legend( loc='best', fontsize = 16)    
        
    pp = 'y'
    txt1 = 'P I D: %s  %s  %s, Scaling: %s'%( metadata['pv_P%s'%pp], metadata['pv_I%s'%pp], metadata['pv_D%s'%pp],metadata['pv_skl%s'%pp] )
    ax.text(x = .4, y= 0.7,  s= txt1, fontsize=14, transform=ax.transAxes) 
    txt2 = 'K%s: %s, Feedback: %s'%(pp, metadata['pv_K%s'%pp], metadata['pv_fd%s'%pp] )
    ax.text(x = .4, y= 0.6,  s= txt2, fontsize=14, transform=ax.transAxes) 
    
    fig.tight_layout()    
    plt.savefig( res_path + filename  + '-fft-time_posX-Y-rms.png')
    plt.show()    
    
    
def plot_posxy_rms(t,posx,posy,res_path, filename  ):
    fig = plt.figure(figsize=(8,6))
    plt.axis('off')
    ax = fig.add_subplot( 211 )
    ax.set_title(filename )
    y = posx        
    ym = y.mean()

    std = (y-ym).std()
    ax.plot(   t, y - ym, '--o', label="X-<X>"   )
    #ax.plot(t, std+ np.zeros_like(t), '--o', label="%s +/-"%round(ym,2) + str( round(std,3) ) + ' um' ) 
    ax.plot(t, std+ np.zeros_like(t), '--o', label= 'rms-'+ str( round(std,3) ) + ' um' )     
    ax.plot(t, -std+ np.zeros_like(t), '--o')#, label="rms-- "+ str( colms[i] ))

    ax.set_xlabel("Time, (s)")
    ax.set_ylabel("X")
    ax.legend( loc='best', fontsize = 16)

    ax = fig.add_subplot( 212 )    
     
    y = posy
    ym=y.mean()

    std = (y-ym).std()
    ax.plot(t, y - ym, '--o', label="Y-<Y>")
   # ax.plot(t,  std + np.zeros_like(t), '--o', label="%s +/-"%round(ym,2) +  str( round(std,3)) + ' um' )
    ax.plot(t,  std + np.zeros_like(t), '--o', label="rms-" +  str( round(std,3)) + ' um' )
    ax.plot(t, -std+ np.zeros_like(t), '--o')#, label="rms-- "+ str( colms[i] ))


    ax.set_xlabel("Time, (s)")
    ax.set_ylabel("Y")    


    ax.legend( loc='best', fontsize = 16) 
    #try filename = filename.rstrip('.txt')
    plt.savefig( res_path + filename  + '-time_posX-Y-rms.png')
    plt.show()



def plot_posxy(t,posx,posy,res_path, filename  ):
    fig = plt.figure(figsize=(8,6))
    plt.axis('off')
    ax = fig.add_subplot( 211 )
    ax.set_title(filename )
 
    y = posx
    ax.plot(t, y, '--o', label="Pos-X")
    ax.set_xlabel("Time, (s)")
    ax.set_ylabel("X")
    ax.legend( loc='best', fontsize = 16)
    
    ax = fig.add_subplot( 212 )    
 
    y = posy
    ax.plot(t, y, '--o', label="Pos-Y")
    ax.set_xlabel("Time, (s)")
    ax.set_ylabel("Y")    
    
    
    ax.legend( loc='best', fontsize = 16) 
    #try filename = filename.rstrip('.txt')
    plt.savefig( res_path + filename + '-time_posX-Y.png')
    plt.show()
    

    
    
    
def plot_fft_posxy(t,posx,posy,res_path, filename,frequency_range=[0,500]  ):
    fig = plt.figure(figsize=(8,6))
    plt.axis('off')    
    ax = fig.add_subplot( 211 )    
    ax.set_title(filename )
    yt = posx         
    freq,fft =  get_fft( t,yt )
    ax.plot(freq,fft, '--o', label="FFT-posX" )
    ax.set_xlabel("freq, (Hz)")
    ax.set_ylabel("fft_x")
    ax.set_xlim( frequency_range)
    ax.legend( loc='best', fontsize = 16)
    ax = fig.add_subplot( 212 ) 
    yt = posy
    freq,fft =  get_fft( t,yt )
    ax.plot(freq,fft, '--o', label="FFT-PosY")
    ax.set_xlabel("freq, (Hz)")
    ax.set_ylabel("fft_y")
    ax.set_xlim( frequency_range)
    ax.legend( loc='best', fontsize = 16) 
    #try filename = filename.rstrip('.txt')
    plt.savefig( res_path + filename.rstrip('.txt') + '-fft_posX-Y.png')
    plt.show()
    

def bpm_read( num_sample, filename=None,rate=10, show_trace=False ):
    '''rate: the data acq rate in unit K'''
    rate = rate*1000.  #10 KHz
    dt =datetime.now()
    times = '%s%02d%02d-%02d%02d' % (dt.year, dt.month, dt.day,dt.hour,dt.minute)
    
    if filename is None:
        filename = times
    else:
        filename = filename + '_' + times

    colms = ['A', 'B','C','D', 'X','Y','Sum']

    num_sample *= rate
    metadata = {}
    for k in list(meta_pvdict.keys()):
        metadata[k] = caget( meta_pvdict[k] )    
 
    if num_sample<13E4:
        caput( pv_num_sam, num_sample)
        caput( pv_trig, 0 )
        caput( pv_trig, 1)
        ti.sleep(num_sample/rate+3)
		# need to wait for the acquisition to finish before pulling data from buffer!
        for n, pv in enumerate([pv_ca, pv_cb, pv_cc, pv_cd, pv_px, pv_py, pv_sumi]):
            d = caget(pv)            
            if n==0: 
                try:
                    zero_point = np.where( d ==0 )[0][0]
                except:
                    zero_point = len( d ) 
                data = d[:zero_point]
            else:
                data = np.vstack( [  data,d[:zero_point]  ] )
                
        np.savetxt( data_path + filename +'.txt', data)        
        t = np.arange( 0, zero_point  ) * 1/rate  
        
        ca = data[0]        
        cb = data[1]
        cc = data[2]
        cd = data[3]
        posx = data[4]
        posy = data[5]
        sumi= data[6]
        #plot_posxy(t,posx,posy, save_path, filename )
        #plot_posxy_rms(t,posx,posy,res_path, filename  )
        #plot_fft_posxy(t,posx,posy,res_path, filename  ) 
        plot_posxy_fft(t,posx,posy,res_path, filename, metadata = metadata  )
        if show_trace:
            plot_current(t, ca,cb,cc,cd, sumi, res_path, filename, metadata = metadata  )
    else:
        print ('The sample number is too large,this number should < 13E4!')
                        

        
        
        
        
        
        
        
        
        
        
    
