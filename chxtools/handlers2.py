import numpy as np
import h5py
from pims import FramesSequence, Frame
from databroker import db


class EigerImages2(FramesSequence):
    def __init__(self, master_filepath, images_per_file, *, md=None):
        self._md = md
        self.master_filepath = master_filepath
        self.images_per_file = images_per_file
        self._handle = h5py.File(master_filepath, 'r')
        try:
            self._entry = self._handle['entry']['data']  # Eiger firmware v1.3.0 and onwards
        except KeyError:
            self._entry = self._handle['entry']          # Older firmwares

    @property
    def md(self):
        return self._md

    @property
    def valid_keys(self):
        valid_keys = []
        for key in sorted(self._entry.keys()):
            try:
                self._entry[key]
            except KeyError:
                pass  # This is a link that leads nowhere.
            else:
                valid_keys.append(key)
        return valid_keys

    def get_frame(self, i):
        dataset = self._entry['data_{:06d}'.format(1 + (i // self.images_per_file))]
        img = dataset[i % self.images_per_file]
        return Frame(img, frame_no=i)

    def __len__(self):
        return sum(self._entry[k].shape[0] for k in self.valid_keys)

    @property
    def frame_shape(self):
        return self[0].shape

    @property
    def pixel_type(self):
        return self[0].dtype

    @property
    def md(self):
        return self._metadata

    @property
    def dtype(self):
        return self.pixel_type

    @property
    def shape(self):
        return self.frame_shape

    def close(self):
        self._handle.close()


class EigerHandler2:
    EIGER_MD_LAYOUT = {
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
    specs = {'AD_EIGER'}
    def __init__(self, fpath, images_per_file):
        # create pims handler
        self._base_path = fpath
        self._images_per_file = images_per_file

    def __call__(self, seq_id):
        master_path = '{}_{}_master.h5'.format(self._base_path, seq_id)
        with h5py.File(master_path, 'r') as f:
            md = {k: f[v].value for k, v in self.EIGER_MD_LAYOUT.items()}
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
        # TODO Return a multi-dimensional PIMS seq.
        return EigerImages2(master_path, self._images_per_file, md=md)

db.fs.register_handler('AD_EIGER2', EigerHandler2)
