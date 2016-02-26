#################################################################################### 
###some latest developed functions for chx beamline
##this is a dynamic collections
##Once functions are incorporated in skxray, those function will be deleted in this file. 
############################################################################



from __future__ import absolute_import, division, print_function

from databroker import DataBroker as db, get_images, get_table, get_events
from filestore.api import register_handler, deregister_handler
from filestore.retrieve import _h_registry, _HANDLER_CACHE

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import skxray.core.roi as roi
from datetime import datetime

import logging
import time
from math import isnan
 
    
import xray_vision
import xray_vision.mpl_plotting as mpl_plot  
from xray_vision.mpl_plotting import speckle
from xray_vision.mask.manual_mask import ManualMask

import skxray.core.roi as roi

import skxray.core.correlation as corr
import skxray.core.utils as utils



import h5py
from filestore.retrieve import HandlerBase
from eiger_io.pims_reader import EigerImages


###################################################    
##for radial average intensity, correct some bugs
###################################################  
def circular_average(image, calibrated_center, threshold=-1, nx=None,
                     pixel_size=None, mask=None):
    """Circular average of the the image data
    The circular average is also known as the radial integration
    Parameters
    ----------
    image : array
        Image to compute the average as a function of radius
    calibrated_center : tuple
        The center of the image in pixel units
        argument order should be (row, col)
    threshold : int, optional
        Ignore counts above `threshold`
    nx : int, optional
        Number of bins in R. Defaults to 100
    pixel_size : tuple, optional
        The size of a pixel (in a real unit, like mm).
        argument order should be (pixel_height, pixel_width)
    Returns
    -------
    bin_centers : array
        The center of each bin in R. shape is (nx, )
    ring_averages : array
        Radial average of the image. shape is (nx, ).
    """
    radial_val = utils.radial_grid(calibrated_center, image.shape,  pixel_size )
    
    if nx is None:
        ps =np.min( pixel_size )
        max_x = np.max(radial_val)/ps
        min_x = np.min(radial_val)/ps        
        nx = int(max_x - min_x)
    #print (nx)

    if mask is  None: mask =1

    bin_edges, sums, counts =  bin_1D(np.ravel(radial_val  * mask ),
                                           np.ravel(image  * mask), nx)
    th_mask = counts > threshold
    ring_averages = sums[th_mask] / counts[th_mask]

    bin_centers = utils.bin_edges_to_centers(bin_edges)[th_mask]

    return bin_centers, ring_averages


def bin_1D(x, y, nx=None, min_x=None, max_x=None):
    """
    Bin the values in y based on their x-coordinates
    Parameters
    ----------
    x : array
        position
    y : array
        intensity
    nx : integer, optional
        number of bins to use defaults to default bin value
    min_x : float, optional
        Left edge of first bin defaults to minimum value of x
    max_x : float, optional
        Right edge of last bin defaults to maximum value of x
    Returns
    -------
    edges : array
        edges of bins, length nx + 1
    val : array
        sum of values in each bin, length nx
    count : array
        The number of counts in each bin, length nx
    """

    # handle default values
    if min_x is None:
        min_x = np.min(x)
    if max_x is None:
        max_x = np.max(x)
    if nx is None:
        nx = int(max_x - min_x)
    
    # use a weighted histogram to get the bin sum
    bins = np.linspace(start=min_x, stop=max_x, num=nx+1, endpoint=True)
    val, _ = np.histogram(a=x, bins=bins, weights=y)
    # use an un-weighted histogram to get the counts
    count, _ = np.histogram(a=x, bins=bins)
    # return the three arrays
    return bins, val, count


###################################################    
##for  eiger detector files
###################################################  

EIGER_MD_DICT = {
    'y_pixel_size': 'entry/instrument/detector/y_pixel_size',
    'x_pixel_size': 'entry/instrument/detector/x_pixel_size',
    'detector_distance': 'entry/instrument/detector/detector_distance',
    'incident_wavelength': 'entry/instrument/beam/incident_wavelength',
    'frame_time': 'entry/instrument/detector/frame_time',
    'beam_center_x': 'entry/instrument/detector/beam_center_x',
    'beam_center_y': 'entry/instrument/detector/beam_center_y',
    'count_time': 'entry/instrument/detector/count_time',
    'pixel_mask': 'entry/instrument/detector/detectorSpecific/pixel_mask',
}

class FixedEigerImages(EigerImages):
    def __init__(self, path, metadata):
        super().__init__(path)
        self._metadata = metadata
    
    @property
    def md(self):
        return self._metadata
    
    @property
    def dtype(self):
        return self.pixel_type
    
    @property
    def shape(self):
        return self.frame_shape

class LazyEigerHandler(HandlerBase):
    specs = {'AD_EIGER'} | HandlerBase.specs
    def __init__(self, fpath, frame_per_point, mapping=None):
        # create pims handler
        self.vals_dict = EIGER_MD_DICT.copy()
        if mapping is not None:
            self.vals_dict.update(mapping)
        self._base_path = fpath
        self.fpp = frame_per_point

    def __call__(self, seq_id):
        import h5py
        master_path = '{}_{}_master.h5'.format(self._base_path, seq_id)
        md = {}
        print('hdf5 path = %s' % master_path)
        with h5py.File(master_path, 'r') as f:
            md = {k: f[v].value for k, v in self.vals_dict.items()}
        # the pixel mask from the eiger contains:
        # 1  -- gap
        # 2  -- dead
        # 4  -- under-responsive
        # 8  -- over-responsive
        # 16 -- noisy
        pixel_mask = md['pixel_mask']
        pixel_mask[pixel_mask>0] = 1
        pixel_mask[pixel_mask==0] = 2
        pixel_mask[pixel_mask==1] = 0
        pixel_mask[pixel_mask==2] = 1
        md['framerate'] = 1./md['frame_time']
        # TODO Return a multi-dimensional PIMS seq
        return FixedEigerImages(master_path, md)

deregister_handler('AD_EIGER')
_HANDLER_CACHE.clear()
register_handler('AD_EIGER', LazyEigerHandler)


def print_attrs(name, obj):
    print(name)
    for key, val in obj.attrs.items():
        print("    %s: %s" % (key, val))
        
        
        
###################################################    
##for  image process
###################################################  

class Reverse_Coordinate(object):
    def __init__(self, indexable, mask):
        self.indexable = indexable
        self.mask = mask
        self.shape = indexable.shape
        self.length= len(indexable)
    def __getitem__(self, key ):      
        if self.mask is not None:
            img =self.indexable[key] * self.mask  
        else:
            img = self.indexable[key]
            
        if len(img.shape) ==3:
            img_=img[:,::-1,:]
        if len(img.shape)==2:
            img_=img[::-1,:] 
        return img_
 

class RemoveHotSpots(object):
    def __init__(self, indexable, threshold= 1E7 ):
        self.indexable = indexable
        self.threshold = threshold        
        
        try:
            self.N = len( indexable  ) 
        except:
            self.N= indexable.length
    def _get_mask(self, Ns=None,Ne=None ):
        mask = np.ones_like(np.array(self.indexable[0]))
        if Ns is None:Ns=0
        if Ne is None:Ne=self.N
        #for key in range(self.N):
        for key in range( Ns,Ne ):    
            data = np.array( self.indexable[key])  #.copy()        
            badp = np.where(  data >= self.threshold )
            if len(badp[0])!=0:                
                mask[badp] = 0                            
        return mask
    def __getitem__(self, key):
        return self.indexable[key] * mask
    
 

class Masker(object):
    def __init__(self, indexable, mask):
        self.indexable = indexable
        self.mask = mask
        self.length = len( indexable)
    def __getitem__(self, key):        
        img =self.indexable[key] * self.mask        
        return img

  


def view_image(imgsr,i):
    #from ipywidgets import interact
    fig, ax = plt.subplots()
    ax.imshow(imgsr[i], interpolation='nearest', cmap='viridis',
                  origin='lower', norm= LogNorm(vmin=0.001, vmax=1e1 ) )
    ax.set_title("Browse the Image Stack")
    plt.show()
    
    
    
 
def view_image_movie(imgsr,sleeps=1, ims=0, ime = 1):    
    fig, ax = plt.subplots()  
    for i in range( ims, ime  ):
        ax.imshow(imgsr[i],  interpolation='nearest', cmap='viridis',
                  origin='lower', norm= LogNorm( vmin=0.001, vmax=1e1 ) )
        ax.set_title("images_%s"%i)
        time.sleep( sleeps )        
        plt.draw()
        if i!=ime-1:
            ax.cla()
        
        
        
def average_img( imgs, Ns=None,Ne = None ):
    ''' Do imgs average,
        Options:
        imgs: the image seriers
        Ns: the start image
        Ne: the last image
        e.g.,
        ave = average_img(imgs)'''
    import numpy as np 
    ave = np.zeros_like(imgs[0],dtype =float)
    #if Ns is None:Ns=0
    #if Ne is None:Ne=len(imgs)
    #if Ne>len(imgs):Ne=len(imgs)
    for i in range(Ns,Ne):
        ave += imgs[i]
    ave /= (Ne-Ns)
    return ave
        
    
    

###################################################    
##for  calculation time
################################################### 

def run_time(t0):
    '''Calculate running time of a program
    Parameters
    ----------
    t0: time_string, t0=time.time()
        The start time
    Returns
    -------
    Print the running time 
    
    One usage
    ---------
    t0=time.time()
    .....(the running code)
    run_time(t0)
    '''    
    
    elapsed_time = time.time() - t0
    print ('Total time: %.2f min' %(elapsed_time/60.)) 
    

###################################################    
##for  file operation
################################################### 
    
def cpopen(  filename=None, inDir=None,  ): 
    import _pickle as cPickle
    import os    
    if inDir!=None:filename=inDir + filename
    if os.path.isfile(filename):
        #fp=file(filename,'rb')   
        fp=open(filename,'rb')  
        data = cPickle.load(fp)
        fp.close()
        return data
    else:
        return None
    


##########################################
###Functions for GiSAXS
##########################################


def make_gisaxs_grid( qr_w= 10, qz_w = 12, dim_r =100,dim_z=120):
    y, x = np.indices( [dim_z,dim_r] )
    Nr = int(dim_r/qp_w)
    Nz = int(dim_z/qz_w)
    noqs = Nr*Nz
    
    ind = 1
    for i in range(0,Nr):
        for j in range(0,Nz):        
            y[ qr_w*i: qr_w*(i+1), qz_w*j:qz_w*(j+1)]=  ind
            ind += 1 
    return y 


###########################################
#for Q-map, convert pixel to Q
########################################### 


def get_incident_angles( inc_x0, inc_y0, refl_x0, refl_y0, pixelsize=[75,75], Lsd=5.0):
    ''' giving: incident beam center: bcenx,bceny
                reflected beam on detector: rcenx, rceny
                sample to detector distance: Lsd, in meters
                pixelsize: 75 um for Eiger4M detector
        get incident_angle (alphai), the title angle (phi)
    '''
    px,py = pixelsize
    phi = np.arctan2( (refl_x0 - inc_x0)*px *10**(-6), (refl_y0 - inc_y0)*py *10**(-6) )    
    alphai = np.arctan2( (refl_y0 -inc_y0)*py *10**(-6),  Lsd ) /2.     
    #thetai = np.arctan2(  (rcenx - bcenx)*px *10**(-6), Lsd   ) /2.  #??   
    
    return alphai,phi 
    
       
def get_reflected_angles(inc_x0, inc_y0, refl_x0, refl_y0, thetai=0.0,
                         pixelsize=[75,75], Lsd=5.0,dimx = 2070.,dimy=2167.):
    
    ''' giving: incident beam center: bcenx,bceny
                reflected beam on detector: rcenx, rceny
                sample to detector distance: Lsd, in meters                
                pixelsize: 75 um for Eiger4M detector
                detector image size: dimx = 2070,dimy=2167 for Eiger4M detector
        get  reflected angle alphaf (outplane)
             reflected angle thetaf (inplane )
    '''    
    
    alphai, phi =  get_incident_angles( inc_x0, inc_y0, refl_x0, refl_y0, pixelsize, Lsd)
    print ('The incident_angle (alphai) is: %s'%(alphai* 180/np.pi))
    px,py = pixelsize
    y, x = np.indices( [dimy,dimx] )    
    #alphaf = np.arctan2( (y-inc_y0)*py*10**(-6), Lsd )/2 - alphai 
    alphaf = np.arctan2( (y-inc_y0)*py*10**(-6), Lsd )  - alphai 
    thetaf = np.arctan2( (x-inc_x0)*px*10**(-6), Lsd )/2 - thetai   
    
    return alphaf,thetaf, alphai, phi    
    
    

def convert_gisaxs_pixel_to_q( inc_x0, inc_y0, refl_x0, refl_y0, 
                               pixelsize=[75,75], Lsd=5.0,dimx = 2070.,dimy=2167.,
                              thetai=0.0, lamda=1.0 ):
    
    ''' 
    
    giving: incident beam center: bcenx,bceny
                reflected beam on detector: rcenx, rceny
                sample to detector distance: Lsd, in meters                
                pixelsize: 75 um for Eiger4M detector
                detector image size: dimx = 2070,dimy=2167 for Eiger4M detector                
                wavelength: angstron               
                
        get: q_parallel (qp), q_direction_z (qz)
                
    '''         
    
    
    alphaf,thetaf,alphai, phi = get_reflected_angles( inc_x0, inc_y0, refl_x0, refl_y0, thetai, pixelsize, Lsd,dimx,dimy)
       
    pref = 2*np.pi/lamda
    
    qx = np.cos( alphaf)*np.cos( 2*thetaf) - np.cos( alphai )*np.cos( 2*thetai)  
    qy_ = np.cos( alphaf)*np.sin( 2*thetaf) - np.cos( alphai )*np.sin ( 2*thetai)    
    qz_ = np.sin(alphaf) + np.sin(alphai)   
    
    qy = qz_* np.sin( phi) + qy_*np.cos(phi) 
    qz = qz_* np.cos( phi) - qy_*np.sin(phi)   
    
    qr = np.sqrt( qx**2 + qy**2 ) 
    
    
    return qx*pref  , qy*pref  , qr*pref  , qz*pref  
        
    
    
    
def get_qedge( qstart,qend,qwidth,noqs,  ):
    ''' DOCUMENT make_qlist( )
    give qstart,qend,qwidth,noqs
    return a qedge by giving the noqs, qstart,qend,qwidth.
           a qcenter, which is center of each qedge 
    KEYWORD:  None    ''' 
    import numpy as np 
    qcenter = np.linspace(qstart,qend,noqs)
    #print ('the qcenter is:  %s'%qcenter )
    qedge=np.zeros(2*noqs) 
    qedge[::2]= (  qcenter- (qwidth/2)  ) #+1  #render  even value
    qedge[1::2]= ( qcenter+ qwidth/2) #render odd value
    return qedge, qcenter    
    
    
###########################################
#for plot Q-map 
###########################################     
    
def get_qmap_label( qmap, qedge ):
    import numpy as np
    '''give a qmap and qedge to bin the qmap into a label array'''
    edges = np.atleast_2d(np.asarray(qedge)).ravel()
    label_array = np.digitize(qmap.ravel(), edges, right=False)
    label_array = np.int_(label_array)
    label_array = (np.where(label_array % 2 != 0, label_array, 0) + 1) // 2
    label_array = label_array.reshape( qmap.shape )
    return label_array
        
    
    
def get_qzrmap(label_array_qz, label_array_qr, qz_center, qr_center   ):
    '''get   qzrmap  '''
    qzmax = label_array_qz.max()
    label_array_qr_ = np.zeros( label_array_qr.shape  )
    ind = np.where(label_array_qr!=0)
    label_array_qr_[ind ] =  label_array_qr[ ind ] + 1E4  #add some large number to qr
    label_array_qzr = label_array_qz * label_array_qr_  
    
    #convert label_array_qzr to [1,2,3,...]
    uqzr = np.unique( label_array_qzr )[1:]
    
    uqz = np.unique( label_array_qz )[1:]
    uqr = np.unique( label_array_qr )[1:]
    #print (uqzr)
    label_array_qzr_ = np.zeros_like( label_array_qzr )
    newl = np.arange( 1, len(uqzr)+1)
    
    qzc =list(qz_center) * len( uqr )
    qrc= [  [qr_center[i]]*len( uqz ) for i in range(len( uqr ))  ]
    
    for i, label in enumerate(uqzr):
        #print (i, label)
        label_array_qzr_.ravel()[ np.where(  label_array_qzr.ravel() == label)[0] ] = newl[i]    
    
    
    return np.int_(label_array_qzr_), np.array( qzc ), np.concatenate(np.array(qrc ))
    


def show_label_array_on_image(ax, image, label_array, cmap=None,norm=None, log_img=True,alpha=0.3,
                              imshow_cmap='gray', **kwargs):  #norm=LogNorm(), 
    """
    This will plot the required ROI's(labeled array) on the image
    Additional kwargs are passed through to `ax.imshow`.
    If `vmin` is in kwargs, it is clipped to minimum of 0.5.
    Parameters
    ----------
    ax : Axes
        The `Axes` object to add the artist too
    image : array
        The image array
    label_array : array
        Expected to be an unsigned integer array.  0 is background,
        positive integers label region of interest
    cmap : str or colormap, optional
        Color map to use for plotting the label_array, defaults to 'None'
    imshow_cmap : str or colormap, optional
        Color map to use for plotting the image, defaults to 'gray'
    norm : str, optional
        Normalize scale data, defaults to 'Lognorm()'
    Returns
    -------
    im : AxesImage
        The artist added to the axes
    im_label : AxesImage
        The artist added to the axes
    """
    ax.set_aspect('equal')
    if log_img:
        im = ax.imshow(image, cmap=imshow_cmap, interpolation='none',norm=LogNorm(norm),**kwargs)  #norm=norm,
    else:
        im = ax.imshow(image, cmap=imshow_cmap, interpolation='none',norm=norm,**kwargs)  #norm=norm,
        
    im_label = mpl_plot.show_label_array(ax, label_array, cmap=cmap, norm=norm, alpha=alpha,
                                **kwargs)  # norm=norm,
    
    
    return im, im_label    
    
    
 

    
def get_qr_tick_label( qr, label_array_qr, inc_x0):
    ''' 
    Dec 16, 2015, Y.G.@CHX
    get zticks,zticks_label 
    
    Parameters:
         
        qr:  2-D array, qr of a gisaxs image (data)
        label_array_qr: a labelled array of qr map, get by:
                        label_array_qr = get_qmap_label( qr, qz_edge)
    Options:
                        
        inc_x0: x-center of incident beam
    Return:
        rticks: list, r-tick positions in unit of pixel
        rticks_label: list, r-tick positions in unit of real space   
    
    Examples:
        rticks,rticks_label = get_qr_tick_label( qr, label_array_qr)

    '''
        
    rticks =[]
    rticks_label = []
    num =  len( np.unique( label_array_qr ) )
    for i in range( 1, num   ):
        ind =  np.where( label_array_qr==i )[1]
        tick = round( qr[label_array_qr==i].mean(),2)
        if ind[0] < inc_x0 and ind[-1]>inc_x0:
             
            mean1 = int( (ind[np.where(ind < inc_x0)[0]]).mean() )
            mean2 = int( (ind[np.where(ind > inc_x0)[0]]).mean() )
            rticks.append( mean1)
            rticks.append(mean2)
             
            rticks_label.append( tick )
            rticks_label.append( tick )
             
        else:    
            mean = int( ind.mean() )
            rticks.append(mean)
            rticks_label.append( tick )
    #print (rticks)
        
    return np.array(rticks), np.array(rticks_label)
    
def get_qz_tick_label( qz, label_array_qz):  
    ''' 
    Dec 16, 2015, Y.G.@CHX
    get zticks,zticks_label 
    
    Parameters:
         
        qz:  2-D array, qz of a gisaxs image (data)
        label_array_qz: a labelled array of qz map, get by:
                        label_array_qz = get_qmap_label( qz, qz_edge)
 
    Return:
        zticks: list, z-tick positions in unit of pixel
        zticks_label: list, z-tick positions in unit of real space   
    
    Examples:
        zticks,zticks_label = get_qz_tick_label( qz, label_array_qz)

    '''
        
    num =  len( np.unique( label_array_qz ) )
    zticks = np.array( [ int( np.where( label_array_qz==i )[0].mean() ) for i in range( 1,num ) ])
    zticks_label = np.array( [ round( qz[label_array_qz==i].mean(),2) for i in range( 1, num ) ])
    return  zticks,zticks_label 

def show_qz(qz):
    ''' 
    plot qz mape

    '''
        
        
    fig, ax = plt.subplots()
    im=ax.imshow(qz, origin='lower' ,cmap='viridis',vmin=qz.min(),vmax= qz.max() )
    fig.colorbar(im)
    ax.set_title( 'Q-z')
    plt.show()
    
def show_qr(qr):
    ''' 
    plot qr mape

    '''
    fig, ax = plt.subplots()    
    im=ax.imshow(qr, origin='lower' ,cmap='viridis',vmin=qr.min(),vmax= qr.max() )
    fig.colorbar(im)
    ax.set_title( 'Q-r')
    plt.show()    

def show_alphaf(alphaf,):
    ''' 
     plot alphaf mape

    '''
        
    fig, ax = plt.subplots()
    im=ax.imshow(alphaf*180/np.pi, origin='lower' ,cmap='viridis',vmin=-1,vmax= 1.5 )
    #im=ax.imshow(alphaf, origin='lower' ,cmap='viridis',norm= LogNorm(vmin=0.0001,vmax=2.00))    
    fig.colorbar(im)
    ax.set_title( 'alphaf')
    plt.show()    
    

def show_qzr_map(  qr, qz, inc_x0, data=None, Nzline=10,Nrline=10  ):
    
    ''' 
    Dec 16, 2015, Y.G.@CHX
    plot a qzr map of a gisaxs image (data) 
    
    Parameters:
        qr:  2-D array, qr of a gisaxs image (data)
        qz:  2-D array, qz of a gisaxs image (data)
        inc_x0:  the incident beam center x 
         
    Options:
        data: 2-D array, a gisaxs image, if None, =qr+qz
        Nzline: int, z-line number
        Nrline: int, r-line number
        
        
    Return:
        zticks: list, z-tick positions in unit of pixel
        zticks_label: list, z-tick positions in unit of real space
        rticks: list, r-tick positions in unit of pixel
        rticks_label: list, r-tick positions in unit of real space
   
    
    Examples:
        
        ticks = show_qzr_map(  qr, qz, inc_x0, data = None, Nzline=10, Nrline= 10   )
        ticks = show_qzr_map(  qr,qz, inc_x0, data = avg_imgmr, Nzline=10,  Nrline=10   )
    '''
        
    
    import matplotlib.pyplot as plt    
    import copy
    import matplotlib.cm as mcm
    
    
    
    cmap='viridis'
    _cmap = copy.copy((mcm.get_cmap(cmap)))
    _cmap.set_under('w', 0)
    
    
    
    qr_start, qr_end, qr_num = qr.min(),qr.max(), Nzline
    qz_start, qz_end, qz_num = qz.min(),qz.max(), Nrline 
    qr_edge, qr_center = get_qedge(qr_start , qr_end, ( qr_end- qr_start)/(qr_num+100), qr_num )
    qz_edge, qz_center = get_qedge( qz_start,   qz_end,   (qz_end - qz_start)/(qz_num+100 ) ,  qz_num )

    label_array_qz = get_qmap_label( qz, qz_edge)
    label_array_qr = get_qmap_label( qr, qr_edge)
 
    labels_qz, indices_qz = roi.extract_label_indices( label_array_qz  )
    labels_qr, indices_qr = roi.extract_label_indices( label_array_qr  )
    num_qz = len(np.unique( labels_qz ))
    num_qr = len(np.unique( labels_qr ))
     


    fig, ax = plt.subplots()
    if data is None:
        data=qr+qz        
        im = ax.imshow(data, cmap='viridis',origin='lower') 
    else:
        im = ax.imshow(data, cmap='viridis',origin='lower',  norm= LogNorm(vmin=0.001, vmax=1e1)) 

    imr=ax.imshow(label_array_qr, origin='lower' ,cmap='viridis', vmin=0.5,vmax= None  )#,interpolation='nearest',) 
    imz=ax.imshow(label_array_qz, origin='lower' ,cmap='viridis', vmin=0.5,vmax= None )#,interpolation='nearest',) 

    caxr = fig.add_axes([0.81, 0.1, 0.03, .8])  #x,y, width, heigth
     
    cba = fig.colorbar(im, cax=caxr    )      
    ax.set_xlabel(r'$q_r$', fontsize=18)
    ax.set_ylabel(r'$q_z$',fontsize=18)
 
    zticks,zticks_label  = get_qz_tick_label(qz,label_array_qz)
    #rticks,rticks_label  = get_qr_tick_label(label_array_qr,inc_x0)

    rticks,rticks_label = zip(*sorted(  zip( *get_qr_tick_label( qr, label_array_qr, inc_x0) ))  )

    stride = int(len(zticks)/7)
    ax.set_yticks( zticks[::stride] )
    yticks =  zticks_label[::stride] 
    ax.set_yticklabels(yticks, fontsize=9)

    stride = int(len(rticks)/7)
    ax.set_xticks( rticks[::stride] )
    xticks =  rticks_label[::stride]
    ax.set_xticklabels(xticks, fontsize=9)


    ax.set_title( 'Q-zr_Map', y=1.03,fontsize=18)
    plt.show()    
    return  zticks,zticks_label,rticks,rticks_label



 

def show_qzr_roi( data, rois, inc_x0, ticks, alpha=0.3):    
        
    ''' 
    Dec 16, 2015, Y.G.@CHX
    plot a qzr map of a gisaxs image with rois( a label array) 
    
    Parameters:
        data: 2-D array, a gisaxs image 
        rois:  2-D array, a label array         
        inc_x0:  the incident beam center x 
        ticks:  zticks, zticks_label, rticks, rticks_label = ticks
            zticks: list, z-tick positions in unit of pixel
            zticks_label: list, z-tick positions in unit of real space
            rticks: list, r-tick positions in unit of pixel
            rticks_label: list, r-tick positions in unit of real space
         
    Options:
        alpha:  transparency of the label array on top of data
        
    Return:
         a plot of a qzr map of a gisaxs image with rois( a label array)
   
    
    Examples:        
        show_qzr_roi( avg_imgr, box_maskr, inc_x0, ticks)
         
    '''
      
        
        
    #import matplotlib.pyplot as plt    
    #import copy
    #import matplotlib.cm as mcm
    
    #cmap='viridis'
    #_cmap = copy.copy((mcm.get_cmap(cmap)))
    #_cmap.set_under('w', 0)
    zticks, zticks_label, rticks, rticks_label = ticks
    avg_imgr, box_maskr = data, rois
    num_qzr = len(np.unique( box_maskr)) -1
    fig, ax = plt.subplots(figsize=(8,8))
    ax.set_title("ROI--Labeled Array on Data")
    im,im_label = show_label_array_on_image(ax, avg_imgr, box_maskr, imshow_cmap='viridis',
                            cmap='Paired', alpha=alpha,
                             vmin=0.01, vmax=30. ,  origin="lower")


    for i in range( 1, num_qzr+1 ):
        ind =  np.where( box_maskr == i)[1]
        indz =  np.where( box_maskr == i)[0]
        c = '%i'%i
        y_val = int( indz.mean() )

        #print (ind[0], ind[1])
        if ind[0] < inc_x0 and ind[-1]>inc_x0:             
            x_val1 = int( (ind[np.where(ind < inc_x0)[0]]).mean() )
            x_val2 = int( (ind[np.where(ind > inc_x0)[0]]).mean() )
            ax.text(x_val1, y_val, c, va='center', ha='center')
            ax.text(x_val2, y_val, c, va='center', ha='center')

        else:
            x_val = int( ind.mean() )
            #print (xval, y)
            ax.text(x_val, y_val, c, va='center', ha='center')

        #print (x_val1,x_val2)

    stride = int(len(zticks)/7)
    ax.set_yticks( zticks[::stride] )
    yticks =  zticks_label[::stride] 
    ax.set_yticklabels(yticks, fontsize=9)

    stride = int(len(rticks)/7)
    ax.set_xticks( rticks[::stride] )
    xticks =  rticks_label[::stride]
    ax.set_xticklabels(xticks, fontsize=9)

    #caxr = fig.add_axes([0.95, 0.1, 0.03, .8])  #x,y, width, heigth     
    #cba = fig.colorbar(im_label, cax=caxr    )  

    #fig.colorbar(im_label, shrink =.85)
    fig.colorbar(im, shrink =.82)
    ax.set_xlabel(r'$q_r$', fontsize=22)
    ax.set_ylabel(r'$q_z$',fontsize=22)

    plt.show()


    
########################
# get one-d of I(q) as a function of qr for different qz
##################### 
    

    
    
def get_1d_qr(  data, Qr,Qz, qr, qz, inc_x0,  mask=None, show_roi=True,  ticks=None, alpha=0.3 ): 
    '''
       plot one-d of I(q) as a function of qr for different qz
       data: a dataframe
       Qr: info for qr, = qr_start , qr_end, qr_width, qr_num
       Qz: info for qz, = qz_start,   qz_end,  qz_width , qz_num
       qr: qr-map
       qz: qz-map
       inc_x0: x-center of incident beam
       mask: a mask for qr-1d integration
       show_roi: boolean, if ture, show the interest ROI
       ticks: ticks for the plot, = zticks, zticks_label, rticks, rticks_label
       alpha: transparency of ROI
       Return: qr_1d, a dict, with keys as qz number
                      qr_1d[key]: 
               Plot 1D cureve as a function of Qr for each Qz  
               
    Examples:
        #to make two-qz, from 0.018 to 0.046, width as 0.008,
        qz_width = 0.008
        qz_start = 0.018 + qz_width/2
        qz_end = 0.046  -  qz_width/2
        qz_num= 2


        #to make one-qr, from 0.02 to 0.1, and the width is 0.1-0.012
        qr_width =  0.1-0.02
        qr_start =    0.02 + qr_width  /2
        qr_end =  0.01 -  qr_width  /2
        qr_num = 1

        Qr = [qr_start , qr_end, qr_width, qr_num]
        Qz=  [qz_start,   qz_end,  qz_width , qz_num ]
        new_mask[ :, 1020:1045] =0
        ticks = show_qzr_map(  qr,qz, inc_x0, data = avg_imgmr, Nzline=10,  Nrline=10   )
        qx, qy, qr, qz = convert_gisaxs_pixel_to_q( inc_x0, inc_y0,refl_x0,refl_y0, lamda=lamda, Lsd=Lsd )
        
        qr_1d = get_1d_qr( avg_imgr, Qr, Qz, qr, qz, inc_x0,  new_mask,  True, ticks, .8)


    '''               
 
    
    qr_start , qr_end, qr_width, qr_num =Qr
    qz_start,   qz_end,  qz_width , qz_num =Qz
    qr_edge, qr_center = get_qedge(qr_start , qr_end, qr_width, qr_num )    
    qz_edge, qz_center = get_qedge( qz_start,   qz_end,  qz_width , qz_num ) 
     
    print ('The qr_edge is:  %s\nThe qr_center is:  %s'%(qr_edge, qr_center))
    print ('The qz_edge is:  %s\nThe qz_center is:  %s'%(qz_edge, qz_center))    
    label_array_qr = get_qmap_label( qr, qr_edge)

    if show_roi:
        label_array_qz0 = get_qmap_label( qz , qz_edge)
        label_array_qzr0,qzc0,qrc0 = get_qzrmap(label_array_qz0, label_array_qr,qz_center, qr_center  )  
        
        if mask is not None:label_array_qzr0 *= mask
        #data_ = data*label_array_qzr0           
        show_qzr_roi( data,label_array_qzr0, inc_x0, ticks, alpha)

    fig, ax = plt.subplots()
    qr_1d ={}
    for i,qzc_ in enumerate(qz_center):
        
        #print (i,qzc_)
        label_array_qz = get_qmap_label( qz, qz_edge[i*2:2*i+2])
        #print (qzc_, qz_edge[i*2:2*i+2])
        label_array_qzr,qzc,qrc = get_qzrmap(label_array_qz, label_array_qr,qz_center, qr_center  )
        #print (np.unique(label_array_qzr ))    
        if mask is not None:label_array_qzr *=   mask
        roi_pixel_num = np.sum( label_array_qzr, axis=0)
        qr_ = qr  *label_array_qzr
        data_ = data*label_array_qzr    
        qr_ave = np.sum( qr_, axis=0)/roi_pixel_num
        data_ave = np.sum( data_, axis=0)/roi_pixel_num     
        qr_1d[i]= [qr_ave, data_ave]
        ax.plot( qr_ave, data_ave,  '--o', label= 'qz= %f'%qzc_)
        
        
    #ax.set_xlabel( r'$q_r$', fontsize=15)
    ax.set_xlabel('$q $'r'($\AA^{-1}$)', fontsize=18)
    ax.set_ylabel('$Intensity (a.u.)$', fontsize=18)
    ax.set_yscale('log')
    #ax.set_xscale('log')
    ax.set_xlim(   qr.max(),qr.min()  )
    ax.legend(loc='best')
    return qr_1d
   
def interp_zeros(  data ): 
    from scipy.interpolate import interp1d
    gf = data.ravel() 
    indice, = gf.nonzero() 
    start, stop = indice[0], indice[-1]+1 
    dx,dy = data.shape 
    x=np.arange( dx*dy ) 
    f = interp1d(x[indice], gf[indice]) 
    gf[start:stop] = f(x[start:stop]) 
    return gf.reshape([dx,dy])     
    
#GiSAXS End
###############################







#import numpy as np

#from . import utils as core
#from . import roi



#####################
##old one-time
#####

from lmfit import minimize, Model, Parameters

logger = logging.getLogger(__name__)


def multi_tau_auto_corr(num_levels, num_bufs, labels, images):
    ##comments, please add start_image, end_image, the default as None
    
    
    from skxray.core import roi
    from skxray.core  import utils as core
    
    """
    This function computes one-time correlations.
    It uses a scheme to achieve long-time correlations inexpensively
    by downsampling the data, iteratively combining successive frames.
    The longest lag time computed is num_levels * num_bufs.
    Parameters
    ----------
    num_levels : int
        how many generations of downsampling to perform, i.e.,
        the depth of the binomial tree of averaged frames
    num_bufs : int, must be even
        maximum lag step to compute in each generation of
        downsampling
    labels : array
        labeled array of the same shape as the image stack;
        each ROI is represented by a distinct label (i.e., integer)
    images : iterable of 2D arrays
        dimensions are: (rr, cc)
    Returns
    -------
    g2 : array
        matrix of normalized intensity-intensity autocorrelation
        shape (num_levels, number of labels(ROI))
    lag_steps : array
        delay or lag steps for the multiple tau analysis
        shape num_levels
    Notes
    -----
    The normalized intensity-intensity time-autocorrelation function
    is defined as
    :math ::
        g_2(q, t') = \frac{<I(q, t)I(q, t + t')> }{<I(q, t)>^2}
    ; t' > 0
    Here, I(q, t) refers to the scattering strength at the momentum
    transfer vector q in reciprocal space at time t, and the brackets
    <...> refer to averages over time t. The quantity t' denotes the
    delay time
    This implementation is based on code in the language Yorick
    by Mark Sutton, based on published work. [1]_
    References
    ----------
    .. [1] D. Lumma, L. B. Lurio, S. G. J. Mochrie and M. Sutton,
        "Area detector based photon correlation in the regime of
        short data batches: Data reduction for dynamic x-ray
        scattering," Rev. Sci. Instrum., vol 70, p 3274-3289, 2000.
    """
    # In order to calculate correlations for `num_bufs`, images must be
    # kept for up to the maximum lag step. These are stored in the array
    # buffer. This algorithm only keeps number of buffers and delays but
    # several levels of delays number of levels are kept in buf. Each
    # level has twice the delay times of the next lower one. To save
    # needless copying, of cyclic storage of images in buf is used.

    if num_bufs % 2 != 0:
        raise ValueError("number of channels(number of buffers) in "
                         "multiple-taus (must be even)")

    if hasattr(images, 'frame_shape'):
        # Give a user-friendly error if we can detect the shape from pims.
        if labels.shape != images.frame_shape:
            raise ValueError("Shape of the image stack should be equal to"
                             " shape of the labels array")

    # get the pixels in each label
    label_mask, pixel_list = roi.extract_label_indices(labels)

    num_rois = np.max(label_mask)

    # number of pixels per ROI
    num_pixels = np.bincount(label_mask, minlength=(num_rois+1))
    num_pixels = num_pixels[1:]

    if np.any(num_pixels == 0):
        raise ValueError("Number of pixels of the required roi's"
                         " cannot be zero, "
                         "num_pixels = {0}".format(num_pixels))

    # G holds the un normalized auto-correlation result. We
    # accumulate computations into G as the algorithm proceeds.
    G = np.zeros(((num_levels + 1)*num_bufs/2, num_rois),
                 dtype=np.float64)

    # matrix of past intensity normalizations
    past_intensity_norm = np.zeros(((num_levels + 1)*num_bufs/2, num_rois),
                                   dtype=np.float64)

    # matrix of future intensity normalizations
    future_intensity_norm = np.zeros(((num_levels + 1)*num_bufs/2, num_rois),
                                     dtype=np.float64)

    # Ring buffer, a buffer with periodic boundary conditions.
    # Images must be keep for up to maximum delay in buf.
    buf = np.zeros((num_levels, num_bufs, np.sum(num_pixels)),
                   dtype=np.float64)

    # to track processing each level
    track_level = np.zeros(num_levels)

    # to increment buffer
    cur = np.ones(num_levels, dtype=np.int64)

    # to track how many images processed in each level
    img_per_level = np.zeros(num_levels, dtype=np.int64)

    start_time = time.time()  # used to log the computation time (optionally)

    for n, img in enumerate(images):

        cur[0] = (1 + cur[0]) % num_bufs  # increment buffer

        # Put the image into the ring buffer.
        buf[0, cur[0] - 1] = (np.ravel(img))[pixel_list]

        # Compute the correlations between the first level
        # (undownsampled) frames. This modifies G,
        # past_intensity_norm, future_intensity_norm,
        # and img_per_level in place!
        _process(buf, G, past_intensity_norm,
                 future_intensity_norm, label_mask,
                 num_bufs, num_pixels, img_per_level,
                 level=0, buf_no=cur[0] - 1)

        # check whether the number of levels is one, otherwise
        # continue processing the next level
        processing = num_levels > 1

        # Compute the correlations for all higher levels.
        level = 1
        while processing:
            if not track_level[level]:
                track_level[level] = 1
                processing = False
            else:
                prev = 1 + (cur[level - 1] - 2) % num_bufs
                cur[level] = 1 + cur[level] % num_bufs

                buf[level, cur[level] - 1] = (buf[level - 1, prev - 1] +
                                              buf[level - 1,
                                                  cur[level - 1] - 1])/2

                # make the track_level zero once that level is processed
                track_level[level] = 0

                # call the _process function for each multi-tau level
                # for multi-tau levels greater than one
                # Again, this is modifying things in place. See comment
                # on previous call above.
                _process(buf, G, past_intensity_norm,
                         future_intensity_norm, label_mask,
                         num_bufs, num_pixels, img_per_level,
                         level=level, buf_no=cur[level]-1,)
                level += 1

                # Checking whether there is next level for processing
                processing = level < num_levels

    # ending time for the process
    end_time = time.time()

    logger.info("Processing time for {0} images took {1} seconds."
                "".format(n, (end_time - start_time)))

    # the normalization factor
    if len(np.where(past_intensity_norm == 0)[0]) != 0:
        g_max = np.where(past_intensity_norm == 0)[0][0]
    else:
        g_max = past_intensity_norm.shape[0]

    # g2 is normalized G
    g2 = (G[:g_max] / (past_intensity_norm[:g_max] *
                       future_intensity_norm[:g_max]))

    # Convert from num_levels, num_bufs to lag frames.
    tot_channels, lag_steps = core.multi_tau_lags(num_levels, num_bufs)
    lag_steps = lag_steps[:g_max]

    return g2, lag_steps


def _process(buf, G, past_intensity_norm, future_intensity_norm,
             label_mask, num_bufs, num_pixels, img_per_level, level, buf_no):
    """
    Internal helper function. This modifies inputs in place.
    This helper function calculates G, past_intensity_norm and
    future_intensity_norm at each level, symmetric normalization is used.
    Parameters
    ----------
    buf : array
        image data array to use for correlation
    G : array
        matrix of auto-correlation function without
        normalizations
    past_intensity_norm : array
        matrix of past intensity normalizations
    future_intensity_norm : array
        matrix of future intensity normalizations
    label_mask : array
        labels of the required region of interests(roi's)
    num_bufs : int, even
        number of buffers(channels)
    num_pixels : array
        number of pixels in certain roi's
        roi's, dimensions are : [number of roi's]X1
    img_per_level : array
        to track how many images processed in each level
    level : int
        the current multi-tau level
    buf_no : int
        the current buffer number
    Notes
    -----
    :math ::
        G   = <I(\tau)I(\tau + delay)>
    :math ::
        past_intensity_norm = <I(\tau)>
    :math ::
        future_intensity_norm = <I(\tau + delay)>
    """
    img_per_level[level] += 1

    # in multi-tau correlation other than first level all other levels
    #  have to do the half of the correlation
    if level == 0:
        i_min = 0
    else:
        i_min = num_bufs//2

    for i in range(i_min, min(img_per_level[level], num_bufs)):
        t_index = level*num_bufs/2 + i

        delay_no = (buf_no - i) % num_bufs

        past_img = buf[level, delay_no]
        future_img = buf[level, buf_no]

        #  get the matrix of auto-correlation function without normalizations
        tmp_binned = (np.bincount(label_mask,
                                  weights=past_img*future_img)[1:])
        G[t_index] += ((tmp_binned / num_pixels - G[t_index]) /
                       (img_per_level[level] - i))

        # get the matrix of past intensity normalizations
        pi_binned = (np.bincount(label_mask,
                                 weights=past_img)[1:])
        past_intensity_norm[t_index] += ((pi_binned/num_pixels
                                         - past_intensity_norm[t_index]) /
                                         (img_per_level[level] - i))

        # get the matrix of future intensity normalizations
        fi_binned = (np.bincount(label_mask,
                                 weights=future_img)[1:])
        future_intensity_norm[t_index] += ((fi_binned/num_pixels
                                           - future_intensity_norm[t_index]) /
                                           (img_per_level[level] - i))

    return None  # modifies arguments in place!


 

 