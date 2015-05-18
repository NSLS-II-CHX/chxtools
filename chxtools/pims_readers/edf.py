import fabio
from pims import ImageSequence


class EDFImages(ImageSequence):

    def imread(self, filepath, plugin=None):
        return fabio.edfimage(filepath).data
