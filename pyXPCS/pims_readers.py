from pims import FramesSequence
from ..broker import DataBroker
from filestore.api import retrieve
from skimage import img_as_float
import fabio


class Images(FramesSequence):
    def __init__(self, headers, name, flat_field, dead_pixels, 
                 process_func=None, dtype=None, as_grey=False):
        """
        Load images from a detector for given Header(s).

        Parameters
        ----------
        headers : Header or list of Headers
        name : str
            alias (data key) of a detector
        flat_field : array
        dead_pixels : array
        process_func: callable, optional
            function to be applied to each image
        dtype : numpy.dtype or str, optional
            data type to cast each image as
        as_grey : boolean, optional
            False by default
            quick-and-dirty way to ensure images are reduced to greyscale
            To take more control over how conversion is performed,
            use process_func above.

        TODO: Add beamstop mask!

        Example
        -------
        >>> header = DataBroker[-1]
        >>> images = Images(header, 'my_detector_lightfield')
        >>> for image in images:
                # do something
        """
        self._dtype = dtype
        events = DataBroker.fetch_events(headers, fill=False)
        self._datum_uids = [event.data[name] for event in events]

        self._validate_process_func(process_func)
        self._as_grey(as_grey, process_func)

    def get_frame(self, i):
        img = retrieve(self._datum_uids[i])
        if self._dtype is not None and img.dtype != self._dtype:
            img = img.astype(self._dtype)
        return Frame(self.process_func(img), frame_no=i)


class SubtractedImages(FramesSequence):
    def __init__(self, headers, lightfield_name, darkfield_name,
                 flat_field, dead_pixels,
                 process_func=None, dtype=None, as_grey=False):
        """
        Load images from a detector for given Header(s). Subtract
        dark images from each corresponding light image automatically.

        Parameters
        ----------
        headers : Header or list of Headers
        lightfield_name : str
            alias (data key) of lightfield images
        darkfield_name : str
            alias (data key) of darkfield images
        flat_field : array
        dead_pixels : array
        process_func: callable, optional
            function to be applied to each image
        dtype : numpy.dtype or str, optional
            data type to cast each image as
        as_grey : boolean, optional
            False by default
            quick-and-dirty way to ensure images are reduced to greyscale
            To take more control over how conversion is performed,
            use process_func above.

        Example
        -------
        >>> header = DataBroker[-1]
        >>> images = SubtractedImages(header, 'my_lightfield', 'my_darkfield')
        >>> for image in images:
                # do something
        """
        self.light = Images(
                headers, lightfield_name, process_func, dtype, as_grey)
        self.dark = Images(
                headers, darkfield_name, process_func, dtype, as_grey)

    def get_frame(self, i):
        # Convert to float to avoid out-of-bounds wrap-around errors,
        # as in 10-11 = 255.
        return img_as_float(self.light[i]) - img_as_float(self.dark[i])

class XrayImageSequence(ImageSequence):
    def __init__(self, filepath, flat_field=None, dead_pixels=None,
                 process_func=None, dtype=None, as_grey=False):
        super(XrayImageSequence, self).__init__(
                filepath, process_func=process_func, as_grey=as_grey)

    def get_frame(self, i):
        frame = super(XrayImageSequence, self).get_frame(i)
        result = img_as_float(frame)
        if flat_field is not None:
            result -= img_as_float(flat_field)
        # TODO implement dead_pixels and possibly beam_stop
        return result


class XraySubtractedImageSequence(FramesSequence):
    def __init__(self, light_filepath, dark_filepath, flat_field=None, dead_pixels=None,
                 process_func=None, dtype=None, as_grey=False):
        self.light = XrayImageSequence(light_filepath, flat_field, dead_pixels,
                                       process_func, dtype, as_grey)
        self.dark = XrayImageSequence(dark_filepath, flat_field, dead_pixels,
                                      process_func, dtype, as_grey)

    def get_frame(self, i):
        # Convert to float to avoid out-of-bounds wrap-around errors,
        # as in 10-11 = 255.
        return img_as_float(self.light[i]) - img_as_float(self.dark[i])

class EDFSequence(XrayImageSequence):

    def imread(self, filename, **kwargs):
        return fabio.edfimage(filename).data

class SubtractedEDFSequence(SubtractedXrayImageSequence):

    def imread(self, filename, **kwargs):
        return fabio.edfimage(filename).data
