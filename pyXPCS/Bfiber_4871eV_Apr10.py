from setupQ import qpix
#from create_and_read_data import read_txt, read_yorick_imgtxt
from numpy import array,load


img_format = 'TIFF'  #'EDF'
PAR= dict(dark_img=False)
mask = None
DK=0
qlist_=None

inDir = '/nfs/xf11id/pyXPCS/'
DATA_DIR='/nfs/xf11id/data/2015/4/10/Bfiber_4871/'
RES_DIR='/nfs/xf11id/data/2015/4/10/'

FILENAME = DATA_DIR+"Bfiber__"

FOUT = "Bfiber_4871eV_"
FITS_FILE = "Bfiber_4871eV_"

title = "Bfiber_4871_waterfall"


EDF=False


####for Bfiber
first_dark = 0
number_darks = 0
lastfile2open=1500 
first_file = 1
number_files = 1400


dimx=2050;
dimy=2448;

qstart = 1400
qstart = 1310
qend = 1900
qwidth = 50



qstart = 1000
qend = 2000
qwidth = 1
noqs = 1
xbar = 143.
ybar =  256-133.
dpix=0.0007
dead_time = 29.8
exposuretime  = 0.2
timeperframe = dead_time +  exposuretime

qperpixel =qpix(dpix, lambda_ = 1.546, Ldet = 2230.)



nobuf= 8  #//must be even!
nolev=  6

begframe = first_file
noframes = number_files

ccdbox = [1,1,dimx,dimy]
imgwidth = dimx

PAR['nobuf']=nobuf
PAR['nolev']=nolev
PAR['begframe'] = begframe
PAR['noframes'] = noframes

mask= load(inDir + 'mask.npy')

#######to read mask generated by Yorick
#mask = read_yorick_imgtxt( inDir + 'mask.txt', dimx,dimy)
#DK = read_yorick_imgtxt( inDir + 'DK.txt', dimx,dimy)
#qlist_= array( [13,17,18,22,23,28,29,34,37,39,48,51,52,57,58,63] )





























