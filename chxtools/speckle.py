"""
X-ray speckle visibility spectroscopy(XSVS) - Dynamic information of
the speckle patterns are obtained by analyzing the speckle statistics
and calculating the speckle contrast in single scattering patterns.
This module will provide XSVS analysis tools
"""

from __future__ import (absolute_import, division, print_function)
import six
import numpy as np
import time

from skxray.core import roi
from skxray.core.utils import bin_edges_to_centers, geometric_series

import logging
logger = logging.getLogger(__name__)

import sys

def xsvs(image_sets, label_array, number_of_img, timebin_num=2, time_bin=None,
         max_cts=None, bad_images = None, threshold=None):   
    """
    This function will provide the probability density of detecting photons
    for different integration times.
    The experimental probability density P(K) of detecting photons K is
    obtained by histogramming the speckle counts over an ensemble of
    equivalent pixels and over a number of speckle patterns recorded
    with the same integration time T under the same condition.
    Parameters
    ----------
    image_sets : array
        sets of images
    label_array : array
        labeled array; 0 is background.
        Each ROI is represented by a distinct label (i.e., integer).
    number_of_img : int
        number of images (how far to go with integration times when finding
        the time_bin, using skxray.utils.geometric function)
    timebin_num : int, optional
        integration time; default is 2
    max_cts : int, optional
       the brightest pixel in any ROI in any image in the image set.
       defaults to using skxray.core.roi.roi_max_counts to determine
       the brightest pixel in any of the ROIs
       
       
    bad_images: array, optional
        the bad images number list, the XSVS will not analyze the binning image groups which involve any bad images
    threshold: float, optional
        If one image involves a pixel with intensity above threshold, such image will be considered as a bad image.
    
    
    Returns
    -------
    prob_k_all : array
        probability density of detecting photons
    prob_k_std_dev : array
        standard deviation of probability density of detecting photons
    Notes
    -----
    These implementation is based on following references
    References: text [1]_, text [2]_
    .. [1] L. Li, P. Kwasniewski, D. Oris, L Wiegart, L. Cristofolini,
       C. Carona and A. Fluerasu , "Photon statistics and speckle visibility
       spectroscopy with partially coherent x-rays" J. Synchrotron Rad.,
       vol 21, p 1288-1295, 2014.
    .. [2] R. Bandyopadhyay, A. S. Gittings, S. S. Suh, P.K. Dixon and
       D.J. Durian "Speckle-visibilty Spectroscopy: A tool to study
       time-varying dynamics" Rev. Sci. Instrum. vol 76, p  093110, 2005.
    There is an example in https://github.com/scikit-xray/scikit-xray-examples
    It will demonstrate the use of these functions in this module for
    experimental data.
    """
    if max_cts is None:
        max_cts = roi.roi_max_counts(image_sets, label_array)

    # find the label's and pixel indices for ROI's
    labels, indices = roi.extract_label_indices(label_array)
    nopixels = len(indices)
    # number of ROI's
    u_labels = list(np.unique(labels))
    num_roi = len(u_labels)

    # create integration times
    if time_bin is None:time_bin = geometric_series(timebin_num, number_of_img)

    # number of times in the time bin
    num_times = len(time_bin)

    # number of pixels per ROI
    num_pixels = np.bincount(labels, minlength=(num_roi+1))[1:]

    # probability density of detecting photons
    prob_k_all = np.zeros([num_times, num_roi], dtype=np.object)

    # square of probability density of detecting photons
    prob_k_pow_all = np.zeros_like(prob_k_all)

    # standard deviation of probability density of detecting photons
    prob_k_std_dev = np.zeros_like(prob_k_all)

    # get the bin edges for each time bin for each ROI
    bin_edges = np.zeros(prob_k_all.shape[0], dtype=prob_k_all.dtype)
    for i in range(num_times):
        bin_edges[i] = np.arange(max_cts*2**i)

    start_time = time.time()  # used to log the computation time (optionally)

    for i, images in enumerate(image_sets):
        # Ring buffer, a buffer with periodic boundary conditions.
        # Images must be keep for up to maximum delay in buf.
        #buf = np.zeros([num_times, timebin_num], dtype=np.object)  # matrix of buffers
        
        buf =  np.ma.zeros([num_times, timebin_num, nopixels]) 
        buf.mask = True   

        # to track processing each time level
        track_level = np.zeros( num_times )
        track_bad_level = np.zeros( num_times )
        # to increment buffer
        cur = np.full(num_times, timebin_num)

        # to track how many images processed in each level
        img_per_level = np.zeros(num_times, dtype=np.int64)

        prob_k = np.zeros_like(prob_k_all)
        prob_k_pow = np.zeros_like(prob_k_all)
 
        try:
            noframes= len(images)
        except:
            noframes= images.length
            
        
        #Num= { key: [0]* len(  dict_dly[key] ) for key in list(dict_dly.keys())  }
        
        for n, img in enumerate(images):
            cur[0] = 1 + cur[0]% timebin_num
            # read each frame
            # Put the image into the ring buffer.
            
            img_ = (np.ravel(img))[indices]   
        
            if threshold is not None:
                if img_.max() >= threshold:
                    print ('bad image: %s here!'%n  )
                    img_ =  np.ma.zeros( len(img_) )
                    img_.mask = True    
                
            if bad_images is not None:        
                if n in bad_images:
                    print ('bad image: %s here!'%n)
                    img_ =  np.ma.zeros( len(img_) )
                    img_.mask = True         
        
            buf[0, cur[0] - 1] = img_

            _process(num_roi, 0, cur[0] - 1, buf, img_per_level, labels,
                     max_cts, bin_edges[0], prob_k, prob_k_pow,track_bad_level)
            
            #print (0, img_per_level)

            # check whether the number of levels is one, otherwise
            # continue processing the next level
            level = 1
            processing=1   
            #print ('track_level: %s'%track_level)
            #while level < num_times:
                #if not track_level[level]:
                    #track_level[level] = 1
            while processing:
                if track_level[level]:
                    prev = 1 + (cur[level - 1] - 2) % timebin_num
                    cur[level] = 1 + cur[level] % timebin_num
                    
                    bufa = buf[level-1,prev-1]
                    bufb=  buf[level-1,cur[level-1]-1] 
                
                    if (bufa.data==0).all():
                        buf[level,cur[level]-1] =  bufa
                    elif (bufb.data==0).all():
                        buf[level,cur[level]-1] = bufb 
                    else:
                        buf[level, cur[level]-1] =  bufa + bufb
                    
                    #print (level, cur[level]-1)
                    
                    
                    track_level[level] = 0

                    _process(num_roi, level, cur[level]-1, buf, img_per_level,
                             labels, max_cts, bin_edges[level], prob_k,
                             prob_k_pow,track_bad_level)
                    level += 1
                    if level < num_times:processing = 1
                    else:processing = 0
                    
                else:
                    track_level[level] = 1
                    processing = 0
                #print ('track_level: %s'%track_level)
            
            if  noframes>=10 and n %( int(noframes/10) ) ==0:
                sys.stdout.write("#")
                sys.stdout.flush() 
            

            prob_k_all += (prob_k - prob_k_all)/(i + 1)
            prob_k_pow_all += (prob_k_pow - prob_k_pow_all)/(i + 1)

    prob_k_std_dev = np.power((prob_k_pow_all -
                               np.power(prob_k_all, 2)), .5)
    
 
    for i in range(num_times):
        if   isinstance(prob_k_all[i,0], float ):
            for j in range( len(u_labels)):
                prob_k_all[i,j] = np.array(  [0] * (len(bin_edges[i]) -1 ) )
                prob_k_std_dev[i,j] = np.array(  [0] * (len(bin_edges[i]) -1 ) )
        
        
    logger.info("Processing time for XSVS took %s seconds."
                "", (time.time() - start_time))
    elapsed_time = time.time() - start_time
    #print (Num)
    print ('Total time: %.2f min' %(elapsed_time/60.)) 
    
    #print (img_per_level - track_bad_level)
    #print (buf)
    
    return prob_k_all, prob_k_std_dev


def _process(num_roi, level, buf_no, buf, img_per_level, labels,
             max_cts, bin_edges, prob_k, prob_k_pow,track_bad_level):
    """
    Internal helper function. This modifies inputs in place.
    This helper function calculate probability of detecting photons for
    each integration time.
    .. warning :: This function mutates the input values.
    Parameters
    ----------
    num_roi : int
        number of ROI's
    level : int
        current time level(integration time)
    buf_no : int
        current buffer number
    buf : array
        image data array to use for XSVS
    img_per_level : int
        to track how many images processed in each level
    labels : array
        labels of the required region of interests(ROI's)
    max_cts: int
        maximum pixel count
    bin_edges : array
        bin edges for each integration times and each ROI
    prob_k : array
        probability density of detecting photons
    prob_k_pow : array
        squares of probability density of detecting photons
    """
    img_per_level[level] += 1
    data = buf[level, buf_no]
    if ( data.data ==0).all():  
        track_bad_level[level] += 1   
 
    #print (img_per_level,track_bad_level)
    
    u_labels = list(np.unique(labels))
    
    if not (data.data  ==0).all():
        for j, label in enumerate(u_labels):        
            roi_data = data[labels == label]
            spe_hist, bin_edges = np.histogram(roi_data, bins=bin_edges,
                                           density=True)
            spe_hist = np.nan_to_num(spe_hist)
            prob_k[level, j] += (spe_hist -
                             prob_k[level, j])/( img_per_level[level] - track_bad_level[level] )

            prob_k_pow[level, j] += (np.power(spe_hist, 2) -
                                 prob_k_pow[level, j])/(img_per_level[level] - track_bad_level[level])


def normalize_bin_edges(num_times, num_rois, mean_roi, max_cts):
    """
    This will provide the normalized bin edges and bin centers for each
    integration time.
    Parameters
    ----------
    num_times : int
        number of integration times for XSVS
    num_rois : int
        number of ROI's
    mean_roi : array
        mean intensity of each ROI
        shape (number of ROI's)
    max_cts : int
        maximum pixel counts
    Returns
    -------
    norm_bin_edges : array
        normalized speckle count bin edges
         shape (num_times, num_rois)
    norm_bin_centers :array
        normalized speckle count bin centers
        shape (num_times, num_rois)
    """
    norm_bin_edges = np.zeros((num_times, num_rois), dtype=object)
    norm_bin_centers = np.zeros_like(norm_bin_edges)
    for i in range(num_times):
        for j in range(num_rois):
            norm_bin_edges[i, j] = np.arange(max_cts*2**i)/(mean_roi[j]*2**i)
            norm_bin_centers[i, j] = bin_edges_to_centers(norm_bin_edges[i, j])

    return norm_bin_edges, norm_bin_centers


def get_bin_edges(num_times, num_rois, mean_roi, max_cts):
    """
    This will provide the normalized bin edges and bin centers for each
    integration time.
    Parameters
    ----------
    num_times : int
        number of integration times for XSVS
    num_rois : int
        number of ROI's
    mean_roi : array
        mean intensity of each ROI
        shape (number of ROI's)
    max_cts : int
        maximum pixel counts
    Returns
    -------
    norm_bin_edges : array
        normalized speckle count bin edges
         shape (num_times, num_rois)
    norm_bin_centers :array
        normalized speckle count bin centers
        shape (num_times, num_rois)
    """
    norm_bin_edges = np.zeros((num_times, num_rois), dtype=object)
    norm_bin_centers = np.zeros_like(norm_bin_edges)
    
    bin_edges = np.zeros((num_times, num_rois), dtype=object)
    bin_centers = np.zeros_like(bin_edges)    
    
    for i in range(num_times):
        for j in range(num_rois):
            bin_edges[i, j] = np.arange(max_cts*2**i)
            bin_centers[i, j] = bin_edges_to_centers(bin_edges[i, j])
            norm_bin_edges[i, j] = bin_edges[i, j]/(mean_roi[j]*2**i)
            norm_bin_centers[i, j] = bin_edges_to_centers(norm_bin_edges[i, j])

    return bin_edges, bin_centers, norm_bin_edges, norm_bin_centers



#################
##for fit
###################

from scipy import stats
from scipy.special import gamma, gammaln






def gammaDist(x,params):
    '''Gamma distribution function
    M,K = params, where K is  average photon counts <x>,
    M is the number of coherent modes,
    In case of high intensity, the beam behavors like wave and
    the probability density of photon, P(x), satify this gamma function.
    '''
    
    K,M = params
    K = float(K)
    M = float(M)
    coeff = np.exp(M*np.log(M) + (M-1)*np.log(x) - gammaln(M) - M*np.log(K))
    Gd = coeff*np.exp(-M*x/K)
    return Gd

def gamma_dist(bin_values, K, M):
    """
    Gamma distribution function
    Parameters
    ----------
    bin_values : array
        scattering intensities
    K : int
        average number of photons
    M : int
        number of coherent modes
    Returns
    -------
    gamma_dist : array
        Gamma distribution
    Notes
    -----
    These implementations are based on the references under
    nbinom_distribution() function Notes

    : math ::
        P(K) =(\frac{M}{<K>})^M \frac{K^(M-1)}{\Gamma(M)}\exp(-M\frac{K}{<K>})
    """

    #gamma_dist = (stats.gamma(M, 0., K/M)).pdf(bin_values)
    x= bin_values
    coeff = np.exp(M*np.log(M) + (M-1)*np.log(x) - gammaln(M) - M*np.log(K))
    gamma_dist = coeff*np.exp(-M*x/K)
    return gamma_dist

 


def nbinom_dist(bin_values, K, M):
    """
    Negative Binomial (Poisson-Gamma) distribution function
    Parameters
    ----------
    bin_values : array
        scattering bin values
    K : int
        number of photons
    M : int
        number of coherent modes
    Returns
    -------
    nbinom : array
        Negative Binomial (Poisson-Gamma) distribution function
    Notes
    -----
    The negative-binomial distribution function
    :math ::
        P(K) = \frac{\\Gamma(K + M)} {\\Gamma(K + 1) ||Gamma(M)}(\frac {M} {M + <K>})^M (\frac {<K>}{M + <K>})^K

    These implementation is based on following references

    References: text [1]_
    .. [1] L. Li, P. Kwasniewski, D. Oris, L Wiegart, L. Cristofolini,
       C. Carona and A. Fluerasu , "Photon statistics and speckle visibility
       spectroscopy with partially coherent x-rays" J. Synchrotron Rad.,
       vol 21, p 1288-1295, 2014.

    """
    co_eff = np.exp(gammaln(bin_values + M) -
                    gammaln(bin_values + 1) - gammaln(M))

    nbinom = co_eff * np.power(M / (K + M), M) * np.power(K / (M + K), bin_values)
    
    return nbinom



#########poisson
def poisson(x,K):    
    '''Poisson distribution function.
    K is  average photon counts    
    In case of low intensity, the beam behavors like particle and
    the probability density of photon, P(x), satify this poisson function.
    '''
    K = float(K)    
    Pk = np.exp(-K)*power(K,x)/gamma(x+1)
    return Pk



def poisson_dist(bin_values, K):
    """
    Poisson Distribution
    Parameters
    ---------
    K : int
        average counts of photons
    bin_values : array
        scattering bin values
    Returns
    -------
    poisson_dist : array
       Poisson Distribution
    Notes
    -----
    These implementations are based on the references under
    nbinom_distribution() function Notes
    :math ::
        P(K) = \frac{<K>^K}{K!}\exp(-<K>)
    """
    #poisson_dist = stats.poisson.pmf(K, bin_values)
    K = float(K)
    poisson_dist = np.exp(-K) * np.power(K, bin_values)/gamma(bin_values + 1)
    return poisson_dist


def diff_mot_con_factor(times, relaxation_rate,
                        contrast_factor, cf_baseline=0):
    """
    This will provide the speckle contrast factor of samples undergoing
    a diffusive motion.

    Parameters
    ----------
    times : array
        integration times

    relaxation_rate : float
        relaxation rate

    contrast_factor : float
        contrast factor

    cf_baseline : float, optional
        the baseline for the contrast factor

    Return
    ------
    diff_contrast_factor : array
        speckle contrast factor for samples undergoing a diffusive motion

    Notes
    -----
    integration times more information - geometric_series function in
    skxray.core.utils module

    These implementations are based on the references under
    negative_binom_distribution() function Notes

    """
    co_eff = (np.exp(-2*relaxation_rate*times) - 1 +
              2*relaxation_rate*times)/(2*(relaxation_rate*times)**2)

    return contrast_factor*co_eff + cf_baseline




def get_roi(data, threshold=1e-3):
    roi = np.where(data>threshold)
    if len(roi[0]) > len(data)-3:
        roi = (np.array(roi[0][:-3]),)                    
    elif len(roi[0]) < 3:
        roi = np.where(data>=0)
    return roi[0]





