#
# GOES SXR Plotting
#
# These routines are used to read and plot the daily/hourly GOES SXR flux.
# The routine get_goes() reads from the NOAA online data in json format.
# The routine goes_std_plots() calls that several times to read all of
# the relevant data needed for the standard 3-day and 6-hour plots and
# generates the plots for the EOVSA status page.
#
# History
#  2020-03-21  DG
#    Initial complete version
#  2020-05-13  DG
#    Removed calls to EOVSA local routines
#  2022-06-24  DG
#    Suddenly the GOES Hi and Lo arrays began having different numbers of
#    times in the json file, so this case has now been taken care of by
#    returning only common times.
#

if __name__ == '__main__':
    import matplotlib
    matplotlib.use("Agg")
    
import urllib2
import json
from astropy.time import Time
import matplotlib.pylab as plt
import numpy as np
from matplotlib.dates import DateFormatter
from util import common_val_idx

def smooth(x,window_len=11,window='hanning'):
    """smooth the data using a window with requested size.

    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.

    input:
        x: the input signal
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal

    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)

    see also:

    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter

    TODO: the window parameter could be the window itself if an array instead of a string
    NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.
    """

    if x.ndim != 1:
        raise ValueError, "smooth only accepts 1 dimension arrays."

    if x.size < window_len:
        raise ValueError, "Input vector needs to be bigger than window size."


    if window_len<3:
        return x


    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"


    s = np.r_[x[window_len-1:0:-1], x, x[-2:-window_len-1:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w = np.ones(window_len,'d')
    else:
        w = eval('np.'+window+'(window_len)')

    y = np.convolve(w/w.sum(), s, mode='valid')
    return y

def get_goes(url='https://services.swpc.noaa.gov/json/goes/primary/xrays-7-day.json'):
    ''' Read a single one of the standard GOES files (determined by the value of the
        url string.)
        
        The routine reads the json file pointed to by the url, and replaces any zeros with nan.
        
        Returns arrays of the GOES low energy (1-8 A) flux, the GOES high-energy (0.5-4 A) flux,
        and a Time object that is the array of UT times.
    '''
    try:
        f = urllib2.urlopen(url)
    except:
        print 'http url error for',url
        return [], [], []
    txt = f.readline()
    goes = json.loads(txt)
    goeshi = []
    goeslo = []
    goestime = []
    goeslotime = []
    for i in goes:
        if i['energy'] == '0.05-0.4nm':
            goeshi.append(i['flux'])
            goestime.append(i['time_tag'])
        else:
            goeslo.append(i['flux'])
            goeslotime.append(i['time_tag'])
    if len(goestime) == 0:
        return [], [], []
    # Convert to arrays
    goestime = np.array(goestime)
    goeslotime = np.array(goeslotime)
    goeslo = np.array(goeslo)
    goeshi = np.array(goeshi)
    # Reduce to common times in case they are different
    idx1, idx2 = common_val_idx(goeslotime, goestime)
    goestime = goestime[idx1]
    goeslo = goeslo[idx1]
    goeshi = goeshi[idx2]
    # Set zeros to nans
    goeslo[np.where(goeslo == 0.0)] = np.nan
    goeshi[np.where(goeshi == 0.0)] = np.nan
    return goeslo, goeshi, Time(goestime)

def goes_std_plots(outpath='./'):
    ''' Generates the urls of the standard files of GOES data, reads them in sequence,
        and makes well-formed plots of the data.  No arguments.  Nothing is returned.
        
        Inputs:
           outpath    Path to folder in which two plots (as .png files) are created,
                         with names goes3-day.png and goes6-hour.png. Default location
                         is the current folder.
    '''
    types = ['3-day','6-hour']
    sources = ['primary','secondary']
    classes = ['A','B','C','M','X']
    for type in types:
        f, ax = plt.subplots(1,1)
        ax.xaxis_date()
        if type == '6-hour':
            ax.xaxis.set_major_formatter(DateFormatter("%H:%M"))
        else:
            ax.xaxis.set_major_formatter(DateFormatter("%d-%H"))
        for source in sources:
            url = 'https://services.swpc.noaa.gov/json/goes/'+source+'/xrays-'+type+'.json'
            lo, hi, t = get_goes(url)
            if len(t) > 0:
                if type == '3-day':
                    hi = smooth(hi,20,'blackman')[10:-9]
                    lo = smooth(lo,20,'blackman')[10:-9]
                ax.plot_date(t.plot_date,lo,'-',label='  1-8 A '+source)
                ax.plot_date(t.plot_date,hi,'-',label='0.5-4 A '+source)
        if type == '3-day':
            # Set start time according to current date
            today = np.floor(Time.now().mjd)
            tstart = Time(today-2, format='mjd')
            tend = Time(today+1, format='mjd')
            plt.xticks(tstart.plot_date + np.arange(10)/3.)
            ax.plot((tstart.plot_date+1)*np.ones(2),[1e-9,1e-2],'k',linewidth=1)
            ax.plot((tstart.plot_date+2)*np.ones(2),[1e-9,1e-2],'k',linewidth=1)
            ax.set_xlabel('Start Date '+tstart.iso[:10]+' [DD-HH]',fontsize=14)
            ax.set_title('GOES SXR 3-Day Plot',fontsize=16)
        else:
            # Set start time according to current hour
            now = int(Time.now().mjd * 24)/24.
            tstart = Time(now-5/24., format='mjd')
            tend = Time(now+1/24., format='mjd')                
            plt.xticks(tstart.plot_date + np.arange(7)/24.)
            ax.set_xlabel(tstart.iso[:10]+' [HH:MM]',fontsize=14)
            ax.set_title('GOES SXR 6-Hour Plot',fontsize=16)
        f.autofmt_xdate(rotation=0, ha='center')
        ax.set_ylim(1e-9,1e-2)
        for i in range(5):
            ax.text(1.01,1/7.*i + 1.4/7,classes[i],transform=ax.transAxes)
        ax.text(0.01,0.97,'Created '+Time.now().iso[:19],fontsize=6,transform=ax.transAxes)
        ax.set_ylabel(r'SXR Flux [W m$\mathregular{^{-2}}$]',fontsize=14)
        ax.set_xlim(tstart.plot_date,tend.plot_date)
        ax.set_yscale('log')
        ax.grid(linestyle='--')
        ax.legend()
        f.set_tight_layout(True)
        plt.savefig(outpath+'goes'+type+'.png',bbox_inches='tight')

if __name__ == '__main__':
    outpath = '/common/webplots/flaremon/'
    goes_std_plots(outpath)
