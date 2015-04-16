from ..waterfall import waterfall, labels_to_indices, indices_to_labels
import numpy as np
from numpy.testing.utils import assert_array_equal
from nose.tools import assert_equal

def assert_dicts_of_arrays_equal(actual, expected):
    for key in expected:
        assert_array_equal(actual[key], expected[key])

def test_basic_usage():
    waterfall(np.ones((5, 5)), {1: [1,2,3]})

def test_labels_indices_roundtrip():
    dict_of_indices = {1: [0, 2, 4], 2: [1, 3, 5]}
    labeled_image = np.array([[1, 2, 1], [2, 1, 2]])
    shape = labeled_image.shape

    # image > dict
    actual = labels_to_indices(labeled_image)
    assert_dicts_of_arrays_equal(actual, dict_of_indices)

    # dict > image
    actual = indices_to_labels(dict_of_indices, shape)
    assert_array_equal(actual, labeled_image)

    # image > dict > image
    roundtrip = indices_to_labels(labels_to_indices(labeled_image), shape)
    assert_array_equal(labeled_image, roundtrip)

    # dict > image > dict
    roundtrip = labels_to_indices(indices_to_labels(dict_of_indices, shape))
    assert_dicts_of_arrays_equal(dict_of_indices, roundtrip)
