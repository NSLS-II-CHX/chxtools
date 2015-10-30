import matplotlib.pyplot as plt
from numpy import indices, where, array, int32, nonzero,log,arange, sqrt
#from GetEdf import get_edf,sumfile
from matplotlib.colors import LogNorm
from modest_image import imshow  #for huge array show
from scipy.interpolate import interp1d
import scipy.misc, scipy.ndimage



def plot(x,y, title='img_',xlabel=None,ylabel=None,
         save=False,show=True,outDir=None):
    plt.plot( x, y, 'bo', ls='--');
    if title!=None:plt.title(  title )
    if xlabel!=None:plt.xlabel(xlabel,fontsize=18)
    if ylabel!=None:plt.ylabel(ylabel,fontsize=18)
    if save:
        if outDir!=None:fp=outDir + title + '_.png'
        else:fp= title + '_.png'
        plt.savefig( fp )        
    if show:plt.show()    


def plot_panels(x, y, subtit=None, title='img_',xlabel=None,ylabel=None,
         save=False,show=True,outDir=None):
    nx,ny = y.shape
    sx= int( round (sqrt(ny)) )
    if ny%sx==0:sy=ny/sx
    else:sy=ny/sx+1
    fig = plt.figure()
    fig.set_size_inches(20,10.5)
    plt.title(title,fontsize=24)
    for sn in range(0, ny):        
        ax = fig.add_subplot(sx,sy,sn+1 )
        if subtit !=None:
            ax.text(.5,.9, '%s'%subtit[sn],fontsize=16,
             horizontalalignment='center',transform=ax.transAxes)
            #sub = [ 'q=%.5f'%s for s in qradi*qperpixel]
        #plt.title('q= %.3f'%qradi[sn],fontsize=16)
        #x=self.dly[1:m+1]        
        y_=y[:,sn]
        plt.plot(x, y_, 'o',linewidth=3, ls='-', color='b',markersize=8)            
        if xlabel!=None:plt.xlabel(xlabel,fontsize=18)
        if ylabel!=None:plt.ylabel(ylabel,fontsize=18)
    if save:
        if outDir!=None:fp=outDir + title + '_.png'
        else:fp= title + '_.png'
        plt.savefig( fp )        
    if show:plt.show()    

        
def show_img(img, save=False, vmin=None, vmax=None,cmap='winter',fontsize=24,
             axis_on=True,title_on=True,xlabel=None,ylabel=None,aspect=1.0,
        title='img_', show=True, logs=False,outDir=None, sizex=9,sizey=9 ):
    """show a two-D image"""
    #show img
    ax = plt.gca()
    F = plt.gcf()
    #print F.get_size_inches()
    F.set_size_inches(sizex,sizey)
    #plt.figure(figsize=(12, 12))
    if vmin==None:vmin=img.min()
    if vmax==None:vmax=img.max()
    if not logs:
        artist = imshow(ax, img,cmap=plt.get_cmap(cmap), vmin=vmin,vmax=vmax )
        
    if logs:
        img= log(img)        
        if vmin==None:vmin=img[nonzero(img)].min() #img.min()
        artist = imshow(ax, img,cmap=plt.get_cmap(cmap),vmin=vmin, vmax=vmax,)
    
    dx,dy=img.shape
    #print dx,dy
    ax.set_ylim(0,dx-1)
    ax.set_xlim(0,dy-1)
    #ax.set_aspect('auto')
    if aspect!=1.0:
        im = ax.get_images()
        extent =  im[0].get_extent()
        ax.set_aspect(abs((extent[1]-extent[0])/(extent[3]-extent[2]))/aspect)
    if title_on:plt.title(  title,fontsize=fontsize )
    #else:plt.title( None)
    if not axis_on:plt.axis('off')
    if xlabel!=None:plt.xlabel(xlabel,fontsize=18)
    if ylabel!=None:plt.ylabel(ylabel,fontsize=18)
    plt.setp(ax.get_xticklabels(),  fontsize=fontsize)
    plt.setp(ax.get_yticklabels(),  fontsize=fontsize)
    if save:
        if outDir!=None:fp=outDir + title + '_.png'
        else:fp= title + '_.png'
        plt.savefig( fp )        
    if show:plt.show()

    

def show_img2(img, save=False, title='img_', show=True, logs=False,outDir=None,aspect=1.0 ):
    """show a two-D image"""
    #show img
    fig, ax = plt.subplots(nrows=1)
    if not logs:
        ax.imshow(img, cmap=plt.get_cmap("Spectral"), interpolation='nearest')
    if logs:
        vmax=img.max()
        vmin=img[nonzero(img)].min() #img.min()
        #vmin=img.min()
        ax.imshow(img, cmap=plt.get_cmap("Spectral"),interpolation='nearest',
                  norm=LogNorm(vmin=vmin, vmax=vmax) )    
    #vmax=img.max()
    #vmin= np.min(img[np.nonzero(img)])
    #plt.imshow(img, cmap=plt.get_cmap("Spectral"), vmin=vmin, vmax=vmax,interpolation='spline36');plt.show()    
    #plt.xlabel('t1',fontsize=12)
    #plt.ylabel('t2',fontsize=12) 
    dx,dy=img.shape 
    ax.set_xlim(0,dx-1)
    ax.set_ylim(0,dy-1) 

    if aspect!=1.0:
        im = ax.get_images()
        extent =  im[0].get_extent()
        ax.set_aspect(abs((extent[1]-extent[0])/(extent[3]-extent[2]))/aspect)
   
    
    
    plt.title(  title )       
    if save:
        if outDir!=None:fp=outDir + title + '_.png'
        else:fp= title + '_.png'
        plt.savefig( fp )        
    if show:plt.show()

def show_mask( pixellist, qind,dimx=150,dimy=150, reverse=False):
    """give the pixellist, qind to show the area defined by pixellist and qind """
    y,x = indices( [ dimx,dimy] )    
    img= (   y*dimy + x ).flatten()
    #x[:]=0
    #img = x.flatten()
    img[:]= -1
    img[pixellist]= qind
    img=img.reshape( [ dimx,dimy] )
    #img=img.T
    if reverse:img = img[::-1]
    show_img( img )
    #return img

def show_img_with_mask(img_fname,pixellist, qind,vmin=None,vmax=None,
                       img_format='EDF',
                       DK=None,show_mask=True,logs=True,reverse=False):
    """give the edf filename, pixellist, qind to show the edf file
        with the interested area defined by pixellist and qind """

    if img_format=='EDF':
        img_fname+='.edf';
        if DK!=None:img= get_edf( img_fname ) -  DK
        else:img=get_edf( img_fname )
    elif img_format=='TIFF':
        img_fname+='.tiff';
        img = scipy.misc.imread(img_fname,flatten=1)
    else:img= cpopen( n= n,prefix= 'data_', inDir=DATA_DIR)

    dimx,dimy= img.shape
    if show_mask:
        img = img.flatten()
        
        img = array( img, dtype = int32)
        #img.dtype='int32'
        img[pixellist]= qind
        img=img.reshape( [ dimx,dimy] )
    show_img( img, logs=logs,vmin=vmin,vmax=vmax )
    #return img  


def reverse_pixellist( pixellist, qind,dimx=150,dimy=150):
    """reverse for python reading pixellist, qind
        generated by a yorick program"""
    y,x = indices( [ dimx,dimy] )
    x[:]=0
    img = x.flatten()
    img[pixellist]= qind
    img=img.reshape( [ dimx,dimy] )
    #img=img.T
    img = img[::-1]
    img = img.flatten()
    pixellist = where( img!=0 )[0]
    qind = img[pixellist]    
    return  pixellist, qind

def interp_zeros(  data ):
    gf = data.ravel()
    indice, = gf.nonzero()
    start, stop = indice[0], indice[-1]+1
    dx,dy = data.shape
    x=arange( dx*dy )
    f = interp1d(x[indice], gf[indice])
    gf[start:stop] = f(x[start:stop])
    return gf.reshape([dx,dy])

    
#def update_progress(progress):
   # print '\r[{0}] {1}%'.format('#'*(progress/10), progress)













    
