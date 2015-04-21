from skxray.correlation import multi_tau_auto_corr


class MTAutocorr(object):
    """
    Example
    -------
    >>> header = db[-1] 
    >>> images = SubtractedImages(headers, 'lightfield1', 'darkfield1', flatfield=some_array, dead_pixels=dead_array)
    >>> multi_tau_autocorr = MTAutocorr(**header)
    >>> whatever = make_labellist_pixellist(q_start, q_end, q_width)
    >>> xpcs.microns_per_pixel = 5
    >>> xpcs.q_list = [...]
    >>> xpcs(images, whatever, num_buffers, num_levels)
    >>> print xpcs.q_list
    >>> xpcs.q_list = [...]
    >>> xpcs(frames, num_buffers, num_levels)
    """

    def __init__(self, detector_distance, pixel_size, wavelength,
                 direct_beam_location, **kwargs):
    # ^ some of this is for creating the mask/labels

    def __call__(images, dict_of_indices, num_buffers, num_levels):
        shape = images.image_shape
        labels = indices_to_labels(dict_of_indices, shape)
        result = mutli_tau_auto_corr(num_levels, num_buffers, labels, images)
        return g2, lag_times
