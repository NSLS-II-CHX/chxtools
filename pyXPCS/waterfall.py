import logging
import time
import six
import numpy as np


logger = logging.getLogger(__name__)


def labels_to_indices(labeled_image):
    """
    Convert an image of labels to a dictionary mapping labels to indices.
    
    Parameters
    ----------
    labeled_image : array

    Returns
    -------
    indices_dict : dict
        mapping labels onto list of indices with that label
        e.g., {1: [1, 5, 6], 2: [2, 3, 4, 8]}
    """
    result = {}
    shape = labeled_image.shape
    for label in np.unique(labeled_image):
        if label == 0:
            continue
        indices = np.where(labeled_image == label)
        result[label] = np.ravel_multi_index(indices, shape)
    return result


def indices_to_labels(dict_of_indices, shape):
    """
    Convert a dictionary mapping labels to indices to an image of labels.

    Parameters
    ----------
    dict_of_indices : dict
        mapping labels onto 1D list of integers counting pixels on the detector
        (This is a "counting" position along the unraveled array, not a 2D
        Cartesian coordinate.)
        along the "unraveled" 1D array of image data
        e.g., {1: [0, 1, 2], 2: [3, 4, 5]}
    shape: tuple
        dimensions of labeled array

    Returns
    -------
    labeled_image : array
    """
    labeled_image = np.empty(shape)
    for label, indices in six.iteritems(dict_of_indices):
        unraveled_indices = np.unravel_index(indices, shape)
        labeled_image[unraveled_indices] = label
    return labeled_image


def waterfall(frames, dict_of_indices):
    """
    Generate a waterfall array from each labeled collection of pixels."

    Parameters
    ----------
    frames : iterable of arrays
       This can be a 3D array, a list of 2D arrays, or any
       iterable of 2D array-like objects.
    dict_of_indices : dict
        mapping labels onto 1D list of integers counting pixels on the detector
        (This is a "counting" position along the unraveled array, not a 2D
        Cartesian coordinate.)
        along the "unraveled" 1D array of image data
        e.g., {1: [0, 1, 2], 2: [3, 4, 5]}

    Returns
    -------
    result : dict
       mapping

    Example
    -------
    # Frames can be any iterable of arrays.
    >>> frames = [my_opening_function(filename) for filename in filenames]
    # Or use pims, which handles memory management nicely:
    >>> frames = pims.ImageSequence(...)
    >>> label_dict = {1: [0, 1, 2, 5, 6, 7], 2: ...}
    >>> waterfall(frames, label_dict)  # all frames are used
    >>> waterfall(frames[100:200])  # frames 100-200 are used
    """
    start_time = time.time()
    labels = dict_of_indices.keys()
    result = {label: np.empty((len(frames), len(dict_of_indices[label])))
              for label in labels}
    for i, image in enumerate(frames):
        logger.info("image %d of %d" % (i, len(frames)))
        for label, indices in six.iteritems(dict_of_indices):
            roi = image.ravel()[indices]
            result[label][i] = roi

    elapsed_time = time.time() - start_time
    logger.info('Total time: %.2f min' , (elapsed_time/60.))

    return result
