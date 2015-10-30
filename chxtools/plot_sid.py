#from chxtools.plot import plot_scan
#plot_scan([-1,-2], x='ivu_gap',y='xray_cam_stats_total1')

def plot_scan( sid = -1, x=None, y=None ):
    '''plot scan_ids,
        Options:
        sid: the scan id, a number or a list
        x: the x-axis, a string
        y: the y-axis, a string
    '''
    from databroker import DataBroker as db,  get_table #get_events, get_images,
    import matplotlib.pyplot as plt

    if not isinstance(sid,list):sid=[sid]
    if x is None:
        x='time'
    if y is None:
        y='time'
    fig,ax=plt.subplots()
    for s in sid:
        dat = get_table(  db[s] )
        #print ('here')
        if x not in dat.keys():
            print('Wrong x input!')
            print('The available X includes:  %s'%list(dat.keys()))
            break
        if y not in dat.keys():
            print('Wrong y input!')
            print('The available Y includes:  %s'%list(dat.keys()))

            break
        #datx=dat[x]
        #daty=dat[y]
        #print (x,y)

        dat.plot( x=x,y=y,ax=ax, label='sid: %s'%s )
        ax.set_ylabel( y )
