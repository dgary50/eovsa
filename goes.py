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
#

if __name__ == '__main__':
    import matplotlib
    matplotlib.use("Agg")
    
import urllib2
import json
from util import Time
import matplotlib.pylab as plt
import numpy as np
from matplotlib.dates import DateFormatter
from autocorrect_tp import smooth

def get_goes(url='https://services.swpc.noaa.gov/json/goes/primary/xrays-3-day.json'):
    ''' Read a single one of the standard GOES files (determined by the value of the
        url string.)
        
        The routine reads the json file pointed to by the url, and replaces any zeros with nan.
        
        Returns arrays of the GOES low energy (1-8 A) flux, the GOES high-energy (0.5-4 A) flux,
        and a Time object that is the array of UT times.
    '''
    f = urllib2.urlopen(url)
    txt = f.readline()
    goes = json.loads(txt)
    goeshi = []
    goeslo = []
    goestime = []
    for i in goes:
        if i['energy'] == '0.05-0.4nm':
            goeshi.append(i['flux'])
            goestime.append(i['time_tag'])
        else:
            goeslo.append(i['flux'])
    goeslo = np.array(goeslo)
    goeslo[np.where(goeslo == 0.0)] = np.nan
    goeshi = np.array(goeshi)
    goeshi[np.where(goeshi == 0.0)] = np.nan
    if len(goestime) == 0:
        return [], [], []
    return goeslo, goeshi, Time(goestime)

def goes_std_plots():
    ''' Generates the urls of the standard files of GOES data, reads them in sequence,
        and makes well-formed plots of the data.  No arguments.  Nothing is returned.
        
        Generates two plots in /common/webplots/flaremon/ as .png files.
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
        plt.savefig('/common/webplots/flaremon/goes'+type+'.png',bbox_inches='tight')

if __name__ == '__main__':
    goes_std_plots()
