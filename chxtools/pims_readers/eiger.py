"""
This reader opens images taken using an Eiger detector at NSLS-II.
It expects a "master file" and a "data file" in the same directory.

It requires h5py and the hdf5-lz4 library from Dectris. The latter
has been packaged for conda and is available on the NSLS-II internal
conda channels and externally from:

conda install -c danielballan hdf5-lz4
"""
import re
import os
import h5py
from pims import FramesSequence, Frame


class EigerImages(FramesSequence):

    pattern = re.compile('(.*)master.*')

    def __init__(self, master_filepath):
        # The 'master' file points to data in other files.
        # Construct a list of those filepaths and check that they exist.
        self.master_filepath = master_filepath
        m = self.pattern.match(os.path.basename(master_filepath))
        if m is None:
            raise ValueError("This reader expects filenames containing "
                             "the word 'master'. If the file was renamed, "
                             "revert to the original name given by the "
                             "detector.")
        prefix = m.group(1)
        with h5py.File(master_filepath) as f:
            self.keys = [k for k in f['entry'].keys() if k.startswith('data')]
            lengths = [f['entry'][key].shape[0] for key in self.keys]
        for k in self.keys:
            filename = prefix + k
            filepath = os.path.join(os.path.dirname(master_filename), filename)
            if not os.path.isfile(filepath):
                raise IOError("Cannot locate expected data file: {0}".format(
                        filepath))
        # Table of Contents return a tuple:
        # self._toc[5] -> [which file, which element in that file]
        self._toc = np.concatenate(
                [zip(i*np.ones(length), np.arange(length))
                for i, length in enumerate(lengths)])

    def get_frame(self, i):
        key_number, elem_number = self._toc[i]
        key = self.keys[key_number]
        with h5py.File(self.master_filepath) as f:
            img = f[key][elem_number]
            return Frame(img, frame_no=i)

    def __len__(self):
        return len(self._toc)

    @property
    def frame_shape(self):
        return self[0].shape

    @property
    def pixel_type(self):
        return self[0].dtype
