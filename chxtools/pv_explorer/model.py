from atom.api import *
from matplotlib.figure import Figure
from matplotlib.axes import Axes
import numpy as np
import datetime

def _get_from_channel_archiver(pv_name, t0=0, t1=10):
    x = np.linspace(t0, t1, 1000)
    y = np.sin(x) * 10
    y += np.random.randn(len(x))
    return x, y


class Model(Atom):
    pv1 = Str()
    pv2 = Str()
    pv3 = Str()
    pv4 = Str()
    _fig = Typed(Figure)
    _data = Dict()
    autolim_axes = Bool(True)
    _axes = Dict()

    t0 = Float()
    t1 = Float()

    dt0 = Typed(datetime.datetime)

    def __init__(self):
        with self.suppress_notifications():
            # plotting initialization
            self.dt0 = datetime.datetime.utcnow()
            self._fig = Figure(figsize=(1, 1))
            self._fig.set_tight_layout(True)
            for name, position in zip(['pv1', 'pv2', 'pv3', 'pv4'],
                                      range(1, 5)):
                self._axes[name] = self._fig.add_subplot(4, 1, position)

    @observe('dt0')
    def dt0_changed(self, changed):
        print(changed)

    @observe('pv1', 'pv2', 'pv3', 'pv4')
    def get_pv1(self, changed):
        print(changed)
        # get the data from the channel archiver
        pv_name = changed['value']
        axes = self._axes[changed['name']]
        self._update_data(pv_name, axes)

    def _update_data(self, pv_name, axes):
        x, y = _get_from_channel_archiver(pv_name, self.t0, self.t1)
        # self._data[pv_name] = (x, y)
        axes.cla()
        axes.plot(x, y)
        self.reformat_view()

    @observe('t0', 't1')
    def change_time(self, changed):
        for k, axes in self._axes.items():
            pv_name = getattr(self, k)
            self._update_data(pv_name, axes)

    def reformat_view(self, *args, **kwargs):
        """
        Recompute the limits, rescale the view, reformat the legend and redraw
        the canvas
        """
        # ignore the args and kwargs. They are here so that any function can be
        # connected to this one
        try:
            for k, ax in self._axes.items():
                ax.relim(visible_only=True)
                ax.autoscale_view(tight=True)
            self._fig.canvas.draw()
        except AttributeError:
            pass
