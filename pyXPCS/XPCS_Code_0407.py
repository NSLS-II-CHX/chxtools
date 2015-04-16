######################################################################################
##Revised Based on Yorick_Multitau_Code Obtained from Dr. Andrei Fluerasu###############
####################   Coded by Dr. Yugang Zhang #################################
####################   631-885-4714   ############################################
####################   Email: yuzhang@bnl.gov     ####################################
########################################################################################
#################### or: yuzhangnew@icloud.com      ##################################
####################  At Brookhaven National Lab   #########################################
############################## Dec 20, 2015, Ver0##########################################
##############################Jan 26, 2015 ,Ver1, Add array_method################################
#######Feb 01,Ver2, Add read_pixellist, qind, qradi from a txt file############
########Feb 18,Ver3, Add two-time correlation function from large array ########
######Feb 22,Ver4, Add two-time correlation function by image_read method#######
#######Feb 28,Ver5, corrected time-index of two-time correlation function by image_read method##########
##########March 1,Ver6 add the method get one-time from the diag of two-time##################
##########March 2,Ver7 do cleaning##############################################
##########March 7,Ver8 add frame_inten,fram_pixel_waterfall,two_time_show2##################################
#######################################################################################
######################################################################################
######################################################################################
from numpy import arange, array,zeros,linspace,round,hypot, ones,sqrt, int_, hstack, vstack
from numpy import digitize,where,copy,indices,histogram,bincount,tril,diag
from numpy import savetxt, log, column_stack,intersect1d, save, load
from math import sin,cos,pi
#from scipy.interpolate import interp1d
#import matplotlib.pyplot as plt
import time
import os
import sys
#from GetEdf import get_edf,sumfile
#from create_and_read_data import cpopen,cpdump
#from create_and_read_data import make_img_array,make_img_pixel_array
#from load_data import  read_ring_data, read_txt
#from img_process import show_img,show_mask,show_edf_with_mask, interp_zeros
from img_process import *
from Init_for_Bfiber import * # the setup file

T=True
F=False

class xpcs( object):
    def __init__(self, PAR=PAR):
        """ DOCUMENT __init__( PAR )
        the initilization of the XPCS class
        KEYWORD:  PAR, gives the parameters, including nobuf,nolev,begframe,
        noframes, if dark_img, it will use DK as dark_img """
            
        global nobuf,nolev,begframe,noframes
        self.version='version_11'
        #print PAR
        if PAR['dark_img']:
            if DK==None:print 'Get dark_img from data.';self.backg()  #to get a background
            else:print 'Get dark_img from a file.'
        else:print 'No Dark Images are used!'
        nobuf=PAR['nobuf']
        nolev=PAR['nolev']
        begframe=PAR['begframe'] 
        noframes=PAR['noframes'] 
        self.noframes=noframes
        self.begframe = begframe
    def clear(self,):
        ''' DOCUMENT clear( )
        clear all the global variables
        KEYWORD:  None   '''
        global buf,G,IAP,IAF,num,cts,cur,g2,gmax,sigg2
        buf=0;G=0;IAP=0;
        IAF=0;num=0;cts=0;
        cur=0;g2=0;gmax=0;sigg2=0
        ttx=0;Ndel=0;Npix=0;
        g2a=0;g2a2=0;g12=0;g12a=0;g12L=0;g12s=0;
        g12x=0;g12y=0;g12z=0
        g12r=0;
        
        
    def backg(self,):
        ''' DOCUMENT backg( )
        return DK: the average of a series of dark images.
        KEYWORD:  None   '''  
        global DK
        DK = sumfile( FILENAME, first_dark, first_dark+number_darks-1, avgflag=1, ext='.edf')

    def delays(self,time=1,nolevs=None,nobufs=None, correct=True): 
        ''' DOCUMENT delays(time=)
        return array of delays.
        KEYWORD:  time: scale delays by time ( should be time between frames)
        '''
        global nolev,nobuf,dly, dict_dly
        if nolevs==None:nolevs=nolev
        if nobufs==None:nobufs=nobuf
        if nobufs%2!=0:print "nobuf must be even!!!"    
        dly=zeros( (nolevs+1)*nobufs/2 +1  )        
        dict_dly ={}
        for i in range( 1,nolevs+1):
            if i==1:imin= 1
            else:imin=nobufs/2+1
            ptr=(i-1)*nobufs/2+ arange(imin,nobufs+1)
            dly[ptr]= arange( imin, nobufs+1) *2**(i-1)            
            dict_dly[i] = dly[ptr-1]            
        dly*=time
        #self.dly=dly*timeperframe
        if correct:self.dly=dly[1:]*timeperframe
        else:self.dly=dly[:-1]*timeperframe
        
    def get_min_dly(self):
        ''' DOCUMENT get_min_dly( )
        return the min value of the dly_array.
        KEYWORD:  None
        '''
                
        global min_ind_dly  
        if len(where(dly-noframes >=0 )[0])==0:min_ind_dly = len(dly)
        else:min_ind_dly =  min( where(dly-noframes >0)[0] ) -1
        #print min_ind_dly
        
    def make_qlist(self):
        ''' DOCUMENT make_qlist( )
        return a qlist by giving the noqs, qstart,qend,qwidth.
        KEYWORD:  None
        ''' 
        global qwidth,noqs,qend,qstart,qnum, qlist, qradi
        global qlist,qradi
        qradi = linspace(qstart,qend,noqs)
        qlist=zeros(2*noqs)
        qlist[::2]= round(qradi-qwidth/2)  #render  even value
        qlist[1::2]= round(qradi+(1+qwidth)/2) #render odd value

        qlist[::2]= int_(qradi-qwidth/2)  #render  even value
        qlist[1::2]= int_(qradi+(1+qwidth)/2) #render odd value
        if qlist_!=None:qlist=qlist_       

    def calqlist(self, qmask=None ,  shape='circle',
                 asym=False,theta=0, center=None):
        ''' DOCUMENT calqlist( qmask=,shape=, )
        calculate the equvilent pixel with a shape,
        return pixellist,nopixels,qind,nopr        
        KEYWORD:  qmask, a mask file;
            shape='circle', give a circle shaped qlist
            shape='column', give a column shaped qlist
            shape='row', give a row shaped qlist
            shape='slope', give a line-shape qlist  by
            the rotation of a row line with center and tan(theta) degree 
            "to be finished yet for the theta"
            
        '''
        global pixellist,nopixels,qind,nopr
        self.make_qlist()
        y, x = indices( [dimx,dimy] )
        if shape=='circle':
            y_= y- ybar +1;x_=x-xbar+1        
            r= int_( hypot(x_, y_)    + 0.5  )#why add 0.5?
        elif shape=='column': 
            r= x
        elif shape=='row':
            r=y
        elif shape=='slope':
            if center==None:cx,cy=[(dimy/2.),(dimx/2.)]
            else:cx,cy=center[1]/1.,center[0]/1.
            a=  cx/(cx+theta*cy) 
            b=theta*cx/(cx+theta*cy)
            print a,b,cx,cy
            r = x*a + y*b
            r=array(r,dtype=int32)
            #print r
        r= r.flatten() 
        #print qlist
        noqrs = len(qlist)    
        qind = digitize(r, qlist)        
        if qmask==None:
            w_= where( (qind)%2 )# qind should be odd;print 'Yes'
            w=w_[0]
        else:
            a=where( (qind)%2 )[0]            
            b=where(  mask.flatten() )[0]            
            #print a.shape,b.shape
            w= intersect1d(a,b)
        
        nopixels=len(w)
        qind=qind[w]/2
        #pixellist= (   y*imgwidth +x ).flatten() [w]
        pixellist= (   y*dimy +x ).flatten() [w]
        nopr,bins=histogram( qind, bins= range( len(qradi) +1 ))       

        
    ###########################################################################
    ########for frame_inten######################################
    ##################################################################
    def fram_inten( self, frame_start, frame_end):
        global fint
        nofra = frame_end+1-frame_start
        fint= zeros( [nofra, noqs])
        start_time = time.time()        
        for i in range(0, nofra ):
            n= frame_start + i
            fp = FILENAME + '%04d'%n 
            
            if img_format=='EDF':fp+='.edf';img= get_edf( fp ) -  DK
            elif img_format=='TIFF':fp+='.tiff';img = scipy.misc.imread(fp,flatten=1)
            else:img= cpopen( n= n,prefix= 'data_', inDir=DATA_DIR)
            
            imgf=img.flatten()[pixellist]
            img=[]            
            fint[i]=  (histogram(qind, bins=noqs, weights= imgf))[0]
            if  int(i %( (frame_end+1-frame_start) /10)) ==0:                
                    sys.stdout.write("#")
                    sys.stdout.flush()
        elapsed_time = time.time() - start_time
        print 'Total time: %.2f min' %(elapsed_time/60.)
        #return fint 
        
    ###########################################################################
    ########for fram_pixel_waterfall######################################
    ##################################################################
    ##fram_pixel_waterfall
        
    def fram_pixel_waterfall( self, frame_start, frame_end, q=None):
        global pintq,pint_dict        
        start_time = time.time()
        nofra = frame_end+1-frame_start        
        if q==None:
            pint_dict={}
            for i in range(noqs):
                pint_dict[i]=zeros( [nofra, nopr[i]])
        else:
            pintq = zeros( [nofra, nopr[q]])
            #print pintq.shape
        for i in range(0, nofra ):
            n= frame_start + i
            fp = FILENAME + '%04d'%n

            if img_format=='EDF':fp+='.edf';img= get_edf( fp ) -  DK
            elif img_format=='TIFF':fp+='.tiff';img = scipy.misc.imread(fp,flatten=1)
            else:img= cpopen( n= n,prefix= 'data_', inDir=DATA_DIR)        
                      
            imgf=img.flatten()[pixellist]
            img=[]            
            if q!=None:                
                pintq[i]= imgf[ qind==q]
                #print pintq[i].shape
            else:
                for q_ in range(noqs):
                    pint_dict[q_][i] = imgf[ qind==q_] 
            if  int(n %( (frame_end+1-frame_start) /10)) ==0:                
                    sys.stdout.write("#")
                    sys.stdout.flush()
        elapsed_time = time.time() - start_time
        print 'Total time: %.2f min' %(elapsed_time/60.)

        
        
    ###########################################################################
    ########for one_time correlation function using image-reading method
    ##################################################################
    def process(self,lev,bufno, n=None):             
        num[lev]+=1  
        if lev==0:imin=0
        else:imin=nobuf/2        
        for i in range(imin, min(num[lev],nobuf) ):
            ptr=lev*nobuf/2+i    
            delayno=(bufno-i)%nobuf 
            IP=buf[lev,delayno]
            IF=buf[lev,bufno]
            G[ptr]+= (  (histogram(qind, bins=noqs, weights= IF*IP))[0]/nopr-G[ptr] )/ (num[lev]-i)
            IAP[ptr]+= (  (histogram(qind, bins=noqs, weights= IP))[0]/nopr-IAP[ptr] )/ (num[lev]-i)
            IAF[ptr]+= (  (histogram(qind, bins=noqs, weights= IF))[0]/nopr-IAF[ptr] )/ (num[lev]-i)
    def insertimg(self, n, norm=None, print_=False, brute=False):        
        cur[0]=1+cur[0]%nobuf 
        fp = FILENAME + '%04d'%n #+'.edf'

        if img_format=='EDF':fp+='.edf';img= get_edf( fp ) -  DK
        elif img_format=='TIFF':fp+='.tiff';img = scipy.misc.imread(fp,flatten=1)
        else:img= cpopen( n= n,prefix= 'data_', inDir=DATA_DIR)
          
        if print_:print 'The insert image %s is %s' %(n,fp)
        buf[0, cur[0]-1 ]=img.flatten()[pixellist]        
        img=[] #//save space    
        self.process(lev=0, bufno=cur[0]-1, n=n )    
        processing=1
        if not brute:
            lev=1
            while processing: 
                if cts[lev]:
                    prev=  1+ (cur[lev-1]-1-1+nobuf)%nobuf
                    cur[lev]=  1+ cur[lev]%nobuf
                    buf[lev,cur[lev]-1] = ( buf[lev-1,prev-1] + buf[lev-1,cur[lev-1]-1] ) /2
                    cts[lev]=0                 
                    self.process(lev= lev, bufno= cur[lev]-1 , n=n)        
                    lev+=1
                    if lev<nolev:processing = 1
                    else:processing = 0                                
                else:
                    cts[lev]=1      #// set flag to process next time
                    processing=0    #// can stop until more images are accumulated
    def autocor( self, print_=False, save_=True, brute=False,filename=None):
        global buf,G,IAP,IAF,num,cts,cur,g2,gmax,sigg2
        global Ndel,Npix
        
        start_time = time.time()
        #initialize all arrays
        buf=zeros([nolev,nobuf,nopixels])  #// matrix of buffers
        cts=zeros(nolev)
        cur=ones(nolev) * nobuf        
        G=zeros( [(nolev+1)*nobuf/2,noqs])
        IAP=zeros( [(nolev+1)*nobuf/2,noqs])
        IAF=zeros( [(nolev+1)*nobuf/2,noqs])
        num= array(zeros(  nolev ),dtype='int')        
        ttx=0        
        for n in range(1,noframes +1 ):
            self.insertimg(begframe+n-1, print_=print_,brute=brute)            
            if  n %(noframes/10) ==0:
                sys.stdout.write("#")
                sys.stdout.flush()
        elapsed_time = time.time() - start_time
        print 'Total time: %.2f min' %(elapsed_time/60.)
        #print G.shape    
        if len(where(IAP==0)[0])!=0:gmax = where(IAP==0)[0][0]        
        else:gmax=IAP.shape[0]
        #g2=G/(IAP*IAF)
        g2=(G[:gmax]/(IAP[:gmax]*IAF[:gmax]))       

        if save_:
            if filename==None:filename='g2_-%s-%s_ImgReadMethod_'%(
            begframe,begframe+noframes-1)
            save( RES_DIR + filename+FOUT, g2)
            print 'the %s was stored in %s'%(filename,RES_DIR)            
        return g2,elapsed_time/60.




    ###########################################################################
    ########for two_time correlation function using image-reading method
    ##################################################################
    def process_two_time(self,lev,bufno,n):               
        num[lev]+=1  
        if lev==0:imin=0
        else:imin=nobuf/2        
        for i in range(imin, min(num[lev],nobuf) ):
            ptr=lev*nobuf/2+i    
            delayno=(bufno-i)%nobuf #//cyclic buffers            
            IP=buf[lev,delayno]
            IF=buf[lev,bufno]
            I_t12 =  (histogram(qind, bins=noqs, weights= IF*IP))[0]
            I_t1  =  (histogram(qind, bins=noqs, weights= IP))[0]
            I_t2  =  (histogram(qind, bins=noqs, weights= IF))[0]
            tind1 = (n-1);tind2=(n -dly[ptr] -1)
            if not isinstance( n, int ):                
                nshift = 2**(lev-1)                
                for i in range( -nshift+1, nshift +1 ):
                    #print tind1+i
                    g12[ int(tind1 + i), int(tind2 + i) ] =I_t12/( I_t1 * I_t2) * nopr
            else:
                #print tind1
                g12[ tind1, tind2 ]  =   I_t12/( I_t1 * I_t2) * nopr
                
    def insertimg_twotime(self, n, norm=None, print_=False):        
        cur[0]=1+cur[0]%nobuf  # increment buffer        
        fp = FILENAME + '%04d'%n

        if img_format=='EDF':fp+='.edf';img= get_edf( fp ) -  DK
        elif img_format=='TIFF':fp+='.tiff';img = scipy.misc.imread(fp,flatten=1)
        else:img= cpopen( n= n,prefix= 'data_', inDir=DATA_DIR)            
                 
        if print_:print 'The insert image %s is %s' %(n,fp)
        buf[0, cur[0]-1 ]=img.flatten()[pixellist]        
        img=[] #//save space 
        countl[0] = 1+ countl[0]
        current_img_time = n - begframe +1
        self.process_two_time(lev=0, bufno=cur[0]-1,n=current_img_time )     
        time_ind[0].append(  current_img_time   )
        processing=1
        lev=1
        while processing:
            if cts[lev]:
                prev=  1+ (cur[lev-1]-1-1+nobuf)%nobuf
                cur[lev]=  1+ cur[lev]%nobuf
                countl[lev] = 1+ countl[lev]                                
                buf[lev,cur[lev]-1] = ( buf[lev-1,prev-1] + buf[lev-1,cur[lev-1]-1] ) /2
                cts[lev]=0                
                t1_idx=   (countl[lev]-1) *2
                current_img_time = ((time_ind[lev-1])[t1_idx ] +  (time_ind[lev-1])[t1_idx +1 ] )/2. 
                time_ind[lev].append(  current_img_time      )                
                self.process_two_time(lev= lev, bufno= cur[lev]-1,n=current_img_time )                 
                lev+=1
                #//Since this level finished, test if there is a next level for processing
                if lev<nolev:processing = 1
                else:processing = 0                                
            else:
                cts[lev]=1      #// set flag to process next time
                processing=0    #// can stop until more images are accumulated


    def autocor_two_time(self,print_=False,save_=True,filename=None):
        global buf,num,cts,cur,g12, countl        
        global Ndel,Npix
        global time_ind  #generate a time-frame for each level
        global g12x, g12y, g12z #for interpolate
        start_time = time.time()
        buf=zeros([nolev,nobuf,nopixels])  #// matrix of buffers, for store img
        cts=zeros(nolev)
        cur=ones(nolev) * nobuf
        countl = array(zeros(  nolev ),dtype='int')        
        g12 =  zeros( [ noframes,noframes, noqs] ) 
        g12x=[]
        g12y=[]
        g12z=[]        
        num= array(zeros(  nolev ),dtype='int')        
        time_ind ={key: [] for key in range(nolev)}         
        ttx=0        
        for n in range(1,noframes +1 ):   ##do the work here
            self.insertimg_twotime(begframe+n-1, print_=print_)
            if  n %(noframes/10) ==0:
                sys.stdout.write("#")
                sys.stdout.flush()                
        for q in range(noqs):            
            x0 =  g12[:,:,q]
            g12[:,:,q] = tril(x0) +  tril(x0).T - diag(diag(x0))            
        elapsed_time = time.time() - start_time
        print 'Total time: %.2f min' %(elapsed_time/60.)
        if save_:
            if filename==None:
                filename =  'g12_-%s-%s_ImgReadMethod_'%(
            begframe,begframe+noframes-1)+FOUT
                
            save(  RES_DIR + filename+FOUT, g12)
            print 'the %s was stored in %s'%(filename,RES_DIR)
        return g12, (elapsed_time/60.)




    ###########################################################################
    ########for dislay the time correlation function 
    ##################################################################
            
    def show_two_time_correlation( self, data,  interp=False,
            vmin=None,vmax=None,filename='g12', title_on= True,
                    show=True,show_tau_lines= False):
        '''the data should be a N*N array'''        
        if show_tau_lines:
            linS = [];linE=[]
            linS.append( zip( [0]*len(dly[:-1]), dly[:-1] ))
            linE.append( zip(  noframes-1- dly[:-1], [noframes-1]*len(dly[:-1])))
        if interp:data=interp_zeros(data)            
        if show_tau_lines:
            ax = plt.gca()                
            imshow(ax, data, vmin=vmin, vmax=vmax)
            print 'show_tau'
            for i, [ps,pe] in enumerate(zip(linS[0],linE[0])):                    
                if i<nobuf:lev=0
                else:lev = (i-nobuf)/(nobuf/2) +1                
                levn = lev/(nolev-1.0)
                red = 1. - levn
                #blue = levn
                green=levn
                linewidth=(2+lev*.4)*1.0
                linewidth= 2
                ax.plot( [ps[0], pe[0]],[ps[1],pe[1]],
                        linewidth=linewidth, color=(red , green,0)  )
                ax.set_xlim(0,noframes-1)
                ax.set_ylim(0,noframes-1)
                
            fp=RES_DIR + filename+ '.png'
            plt.savefig( fp )  
            if show:plt.show()
        else:
            show_img(data, save=True, show=show, title_on=title_on,
                     title=filename, outDir=RES_DIR) 

    def show_multi_two_time_correlation( self, data=None,  interp=False,
        vmin=None, vmax=None, filename='g12', show=True,show_tau_lines= False):
        if data==None:data=g12            
        if len(data.shape)==3:dx,dy, dq = data.shape
        else:dq=1
        for q in range(dq):
            print q
            if dq>1:y=data[:,:,q]
            else:y=data
            xp.show_two_time_correlation(
                data = y, interp=interp, vmin=vmin, vmax=vmax,
            filename=filename + '_q%s'%q, show=show,show_tau_lines=show_tau_lines)
            y=0
            
    def showg2(self,data='g2',show=True,save_=True,filename='g2' ):

        sx= int( round (sqrt(noqs)) )
        if noqs%sx==0:sy=noqs/sx
        else:sy=noqs/sx+1
        if data=='g2':data=g2 
        m,n=data.shape        
        fig = plt.figure()
        fig.set_size_inches(20,10.5)
        title=filename
        plt.title(title,fontsize=24)
        plt.axis('off')
        for sn in range(0, n):
            #ax = fig.add_subplot('%i'%sx,'%i'%sy,'%i'%(sn+1) )
            ax = fig.add_subplot(sx,sy,sn+1 )            
            ax.text(.5,.9, 'q= %.5f A-1'%(qradi[sn]*qperpixel),fontsize=16,  horizontalalignment='center',
                transform=ax.transAxes)
            #plt.title('q= %.3f'%qradi[sn],fontsize=16)
            
            x=self.dly[:m]
            #x=self.dly[1:m+1]
            
            y=data[:,sn]
            plt.plot(x, y, 'o',linewidth=3, ls='-', color='b',markersize=8)            
            #ax.set_xscale('log')
            plt.xscale('log')
            #print max(y) 
            plt.ylim([    min(y) , max(y[1:]) ])
            plt.xlim([   min(x)+0*1e-6, max(x)])
            #plt.ylim([  1, 2])
            plt.xlabel('time delay',fontsize=12)
            plt.ylabel('g2',fontsize=12)            
        if save_:
            plt.savefig( RES_DIR + filename +'.png' )            
            #cpdump(data,FOUT+'_g2',RES_DIR)
            #cpdump(sigg2,FOUT+'_sigg2',RES_DIR)
            
        if show:plt.show()


############################################
#to get the waterfall plot
       
xp=xpcs();
xp.calqlist( shape='column');
xp.fram_pixel_waterfall(1,1400)
show_img(pint_dict[0][:,:],aspect=3)

#to show a image with mask
#show_img_with_mask(FILENAME +'0001',pixellist, qind,img_format='TIFF',show_mask=T,logs=F)









