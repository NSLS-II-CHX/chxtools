import numpy as np
import subprocess
from dataportal import DataBroker, DataMuxer
from dataportal.broker import EventQueue
import matplotlib.pyplot as plt
import time as ttime
import sys
from ophyd.userapi.scan_api import estimate

def new_queue(header, queue=None):
    if queue is None:
        queue = EventQueue(header)
        return header, queue
    hdr = DataBroker[-1]
    if header.scan_id != hdr.scan_id:
        print("New header found: Scan id = %s. uid = %s" % 
              (hdr.scan_id, hdr.run_start_uid))
        sys.stdout.flush()
        queue = EventQueue(hdr)
        return hdr, queue
    return header, queue

vlines = {'center_of_mass': {'color': 'red'},
          'cen': {'color': 'red', 'ls': '--'},}
hlines = {'avgy': {'color': 'blue',  'ls': '-'}, 
          'ymin': {'color': 'black', 'ls': '--'}, 
          'ymax': {'color': 'black', 'ls': '--'}, }
points = {'cen': {'color': 'red', 'marker': 'o'},
          'fwmh_left': {'color': 'red', 'marker': '<'}, 
          'fwhm_right': {'color': 'red', 'marker': '>'}}
          
def plot1d(y, x=None, scans=None, live=True, sleep_time=1):
    """Plot live data and on-the-fly peak stats estimator

    Parameters
    ----------
    y : str
        The name of the y value to plot
    x : str, optional
        The name of the value to plot on the x axis. If None, defaults 
        to the sequence number of the event (Note that this probably works, 
        but I'm not sure as it has not been tested!)
    scans : list, optional
        List of other scan indices to plot. uses db[] syntax, so any valid
        entry to [] will work
    live : bool, optional
        Grab new data and plot it as it comes off. Defaults to True.
    sleep_time : float, optional
        Time to sleep between data updates. Defaults to 1 sec
    """
    if scans is None:
        scans = []
    lines1 = {}
    fig, (ax1, ax2) = plt.subplots(nrows=2, figsize=(15,10), sharex=True)
    fig.show()
    for scan_id in scans:
        hdr = DataBroker[scan_id]
        events = DataBroker.fetch_events(hdr)
        dm = DataMuxer.from_events(events)
        df = dm.to_sparse_dataframe()
        if x is None:
            old_x = np.asarray(df.index)
        else:
            old_x = np.asarray(df[x])
        old_y = np.asarray(df[y])
        lines1[scan_id], = ax1.plot(old_x, old_y, 'o', ms=15, label=scan_id)
    if x is None:
        ax1.set_xlabel('scan point index')
        ax2.set_xlabel('scan point index')
    else:
        ax1.set_xlabel(x)
        ax2.set_xlabel(x)
    ax1.set_ylabel(y)
    ax2.set_ylabel(y)
    ax1.set_title('data stream')
    ax2.set_title('peak estimator')

    if live:
        hdr = DataBroker[-1]
        scan_id = hdr.scan_id
        while scan_id in lines1:
            ttime.sleep(.5)
            hdr = DataBroker[-1]
            scan_id = hdr.scan_id
        lines1[scan_id], = ax1.plot([], [], 'o', ms=15, label=scan_id)
        queue = None
        prev_stats = None
        while True:
            # loop until killed
            hdr, queue = new_queue(hdr, queue)
            scan_id = hdr.scan_id
            queue.update()
            new_events = queue.get()
            try:
                old_x, old_y = lines1[scan_id].get_data()
                old_x = list(old_x)
                old_y = list(old_y)
            except KeyError:
                lines1[scan_id], = ax1.plot([], [], 'o', ms=15, label=scan_id)
                old_x, old_y = [], []
            if x is None:
                new_x = [event.seq_num for ev in new_events]
            else:
                new_x = [ev['data'][x] for ev in new_events]
            new_y = [ev['data'][y] for ev in new_events]
            new_x = old_x + new_x
            new_y = old_y + new_y
            lines1[scan_id].set_data(new_x, new_y)
            ax1.relim(visible_only=True)
            ax1.legend(loc=0).draggable()
            
            # now deal with axis 2
            try:
                stats = estimate(np.asarray(new_x), np.asarray(new_y))
            except ValueError:
                stats = prev_stats
    #        print(stats)
            if stats != prev_stats:
                ax2.cla()
                ax2.plot(new_x, new_y, 'o', ms=15, label=scan_id)
                ax2.set_title('peak estimator')
                for stat, vals in stats.items():
                    if stat in points:
                        # sometimes 'cen' comes back as one or two values. This
                        # try/except block is a way to do the right thing when 
                        # this happens
                        try:
                            vals[0]
                            ax2.scatter(vals[0], vals[1], label=stat, **points[stat])
                        except IndexError:
                            ax2.axvline(vals, label=stat, **vlines[stat])
                    elif stat in hlines:
                        # draw a horizontal line
                        ax2.axhline(vals, label=stat, **hlines[stat])
                    elif stat in vlines:
                        # draw a vertical line
                        ax2.axvline(vals, label=stat, **vlines[stat])
                prev_stats = stats
                ax2.relim(visible_only=True)
                ax2.legend(loc=0).draggable()
            fig.canvas.draw()
            fig.canvas.flush_events()
            ttime.sleep(sleep_time)

