import nose
from skxray.diff_roi_rings import roi_rings
from ..pims_readers import XrayImageSequence
from ..XPCS_Code_0407 import xpcs


def test_legacy_pipeline():
    a = xpcs()  # configured with globally-defined PAR
    a.calqlist()  # compute labels
    a.fram_inten(1, 10)  # load images
    g2, elapsed_minutes = a.autocor(save_=False)

def test_new_pipeline():
    filename = os.path.join(os.path.dirname(__file__),
                            'data', 'Bfiber_s0p4_overweekend_')
    images = XrayImageSequence(filename)

    ### TODO: WHAT???

    qstart = 1000  # mm
    qend = 2000  # mm
    qwidth = 1  # mm
    noqs = 1
    wavelength = 1.546  # Angstroms
    Ldet = 2230  # mm? detector size?
    dpix=0.0007  
    dead_time = 29.995
    exposuretime  = 0.005
    timeperframe = dead_time +  exposuretime

    qperpixel =qpix(dpix, lambda_ = 1.546, Ldet = 2230.)
    #####

    center = (143, 256 - 133)
    shape = images.image_shape
    label_array = roi_rings(shape, center, noqs, qstart, qwidth)
    mta = MTAutocorr(detector_distance, pixel_size, wavelength=wavelength,
                     direct_beam_location=center)
    actual = mta(frames, label_array, num_buffers=8, num_levels=6)
    actual_g2, actual_lag_times = actual
