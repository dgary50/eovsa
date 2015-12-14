#   2015-May-29  DG
#      Converted from using datime() to using Time() based on astropy

import stateframe, os, datetime, struct
from util import Time
import numpy as np
import pylab as plt

def get_val_from_stflog(filename,loc):
    arr = False
    if len(loc) == 3:
        print 'Cannot plot multiple values--will plot only index 0'
        arr = True
    with open(filename,'rb') as f:
        data = f.read(100)
        recsize = struct.unpack_from('i',data,16)[0]
        tstamp = stateframe.extract(data,['d',0])
        t = Time(tstamp,format='lv')
        f.seek(0)
        data = f.read(recsize)
        i = 0
        nrec = os.stat(filename).st_size/recsize
        tm = np.zeros(nrec,'float')
        val = np.zeros(nrec,'float')
        while len(data) == recsize:
            # Read timestamps
            tm[i] = stateframe.extract(data,['d',0])
            if arr:
                val[i] = stateframe.extract(data,loc)[0]
            else:
                val[i] = stateframe.extract(data,loc)
            data = f.read(recsize)
            i += 1
    # Convert timestamps to Time() array
    times = Time(tm,format='lv')
    bad = np.where(val == 0.0)[0]
    val[bad] = np.NaN
    tm[bad] = np.NaN
    fig = plt.figure()
    sub_plot = fig.add_subplot(111)
    sub_plot.plot_date(times.plot_date,val,'+')
    fmt = plt.DateFormatter('%H:%M')
    sub_plot.xaxis.set_major_formatter(fmt)
    plt.title(t.iso[0:10])
    return times,val
