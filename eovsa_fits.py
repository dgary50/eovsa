from astropy.io import fits
import glob
from util import Time
import numpy as np
import matplotlib.pylab as plt
from matplotlib.dates import DateFormatter

def eovsa_fits2plot(files):
    '''
    '''
    for file in files:
        spec = fits.getdata(file,ext=0)
        freq = fits.getdata(file,ext=1)
        ut = fits.getdata(file,ext=2)
        fghz = freq['sfreq']
        time = Time(ut['mjd']+ut['time']/86400000.,format='mjd').plot_date
        if file != files[0]:
            specs = np.concatenate((specs,spec),1)
            times = np.concatenate((times,time))
        else:
            specs = spec
            times = time
            date = Time(time[0],format='plot_date').iso[:10]
            header = fits.open(file)[0].header
            if header['TYPE'] == 0:
                typstr = 'Undefined'
            elif header['TYPE'] == 1:
                typstr = 'Total Power'
            elif header['TYPE'] == 2:
                typstr = 'Cross Power'
            
    f, ax = plt.subplots(1,1,figsize=(14,5))
    X = np.sort(specs.flatten())   # Sorted, flattened array
    # Set any time gaps to 0
    tdif = times[1:] - times[:-1]
    bad, = np.where(tdif > 120./86400)  # Time gaps > 2 minutes
    specs[:,bad] = 0
    # Overwrite bottom of frequency gap with 0
    fdif = fghz[1:] - fghz[:-1]
    bad, = np.where(fdif > 0.1)
    specs[bad,:] = 0
    vmax = X[int(len(X)*0.85)]  # Clip at 10% of points
    im = ax.pcolormesh(times,fghz,specs,vmax=vmax,vmin=0)
    #    plt.colorbar(im,ax=ax,label='Amplitude [arb. units]')
    ax.xaxis_date()
    ax.xaxis.set_major_formatter(DateFormatter("%H:%M"))
    ax.set_ylim(fghz[0], fghz[-1])
    ax.set_xlabel('Time [UT]')
    ax.set_ylabel('Frequency [GHz]')
    ax.set_title('EOVSA '+typstr+' Dynamic Spectrum for '+date)
    #    f.autofmt_xdate(bottom=0.15)
