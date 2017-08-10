# Purpose: Plot elevation vs time for VLA calibrators and optionally tp sources
#          with C Band flux greater than minflux (in Jy). Thick lines indicate
#          time range visible by the 27-m (due to the +-55 deg hour angle range)        
# History:
#  2016-08-02  BC
#    Hacked from JV's plot_bright_src.py
#  2016-11-12  DG
#    Cleaned up code a bit, changed routine name to whatup()
#  2017-03-07  DG
#    27-m antenna is currently restricted to +/- 50 HA, so change thick 
#    lines accordingly
#  2017-03-23  DG
#    Numerous changes to clean up the plot.  Eliminates horizontal lines
#    bridging gaps, changes legend for readability, and plots three different
#    line thicknesses to distinguish equatorial mount coverage (thickest),
#    az_el mount coverage (intermediate), and source outside of either coverage
#    (thinnest).
#  2017-07-13  DG
#    27-m cable tension is fixed, so go back to +/- 55-degree HA limits.

import os, util
import eovsa_cat
from eovsa_visibility import *
from pylab import *
import copy
from readvla import readvlacaldb
import numpy as np
import matplotlib.cm as cm
from astropy.time import TimeDelta
import matplotlib.dates as mdates

import tooltip as tt

def deg(rad):
    return (rad * 180./pi + 180) % 360 - 180

def whatup(minflux=7., dur=24., t=None, showtp=False):
    ''' Find out what sources are up at what times.  Produces a plot with some
        mouse interactivity (click once to activate, click a second time to deactivate)
        
        Optional arguments:
           minflux: low-flux cutoff for the plot
           dur: duration in hours from the start time
           t: Time() object giving start time of the plot, on None for current time.
           showtp: If True, show total power sources. Default False
    '''
    os.chdir(os.path.expanduser('~')+'/Dropbox/PythonCode/Current')

    if showtp:
        srclist0 = ['Sun','Moon','CAS-A','CYG-A','TAU-A','VIR-A']
    else:
        srclist0 = ['Sun']
    ntp = len(srclist0)
    morelist = []
    moreflux = []
    # additional srclist from VLA calibrator catalogue
    cal = readvlacaldb()
    # find all bright sources with flux > 10 Jy at C band
    for i, c in enumerate(cal):
        if 'C' in c.band and c.flux[c.band.index('C')] > minflux:
            morelist.append(c.name)
            moreflux.append(c.flux[c.band.index('C')])
    # sort the VLA cal src list by flux
    morelist=[l for (f,l) in sorted(zip(moreflux,morelist))]
    moreflux.sort()
    moreflux.reverse()
    morelist.reverse()

    srclist = srclist0 + morelist

    if t is None:
        t = util.Time.now()

    try:
        float(dur)
    except ValueError:
        print 'Bad value for dur. Use default dur=8.'
        dur = 8.

    # Add times in 1-min steps up to duration dur (hours)
    ts = t + TimeDelta(arange(0.,dur,1./60.)/24.,format='jd')

    aa = eovsa_cat.eovsa_array_with_cat()

    nt = len(ts)
    ra = zeros((len(srclist),nt))
    ha = zeros((len(srclist),nt))
    dec = zeros((len(srclist),nt))
    alt = zeros((len(srclist),nt))
    for i in range(nt):
        aa.set_jultime(ts[i].jd)    
        lst = aa.sidereal_time()

        for j,srcname in enumerate(srclist):
            src = aa.cat[srcname]
            src.compute(aa)
            ra[j,i] = src.ra
            ha[j,i] = lst - src.ra
            dec[j,i] = src.dec
            alt[j,i] = src.alt

    ra_deg = deg(ra)
    ha_deg = deg(ha)
    dec_deg = deg(dec)
    alt = deg(alt)
    eq_alt = copy.deepcopy(alt)
    az_alt = copy.deepcopy(alt)
    colors = cm.rainbow(np.linspace(0, 1, len(srclist)))
    f=plt.figure(figsize=(10,6))
    ax=plt.subplot(111)

    for j in range(len(srclist)):
        #idx=where(abs(ha_deg[j]) < 55.)[0]
        bad, = where(abs(ha_deg[j]) > 55.)   # Back to nominal +/- 55 HA
        eq_alt[j,bad] = nan
        #ax.plot_date(ts.plot_date[idx],alt[j,idx],'-',linewidth=6,color=colors[j])
        ax.plot_date(ts.plot_date,eq_alt[j],'-',linewidth=6,color=colors[j])#,alpha=0.5)

    for j in range(len(srclist)):
        # Plot thin lines
        ax.plot_date(ts.plot_date,alt[j],'-',color=colors[j])
    for j in range(len(srclist)):
        bad, = where(abs(alt[j]) < 10.)   # Eliminate points < 10-degrees altitude
        az_alt[j,bad] = nan
        # Plot width = 2 lines
        ax.plot_date(ts.plot_date,az_alt[j],'-',linewidth=2,color=colors[j])
    for j in range(len(srclist)):
        # Replot width = 6 lines
        ax.plot_date(ts.plot_date,eq_alt[j],'-',linewidth=6,color=colors[j])#,alpha=0.5)

    srclegend=[]
    for i,s in enumerate(srclist):
        if i < ntp:
            srclegend.append(s)
        else:
            srclegend.append(s + ' ({0:.1f} Jy)'.format(moreflux[i-ntp]))

    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d %H:%M'))
    ax.set_ylim([0,90])
    def format_coord(x, y, ts, srclist, alt):
        col = np.argmin(np.absolute(ts.plot_date - x))
        nsrc, nt = alt.shape
        if col >= 0 and col < len(ts):
            tm = ts.plot_date[col]
            ft = mdates.DateFormatter('%Y-%m-%d %H:%M')
            timstr = str(ft(tm))
            # find out source name at the given time
            alt_t = alt[:, col]
            ind = np.argmin(np.absolute(alt_t - y))
            srcname = srclist[ind]
            ra0 = np.degrees(ra[ind, col])/15.
            ha0 = ha_deg[ind, col]/15.
            dec0 = dec_deg[ind, col]
            if ind < ntp:
                return '{2} (RA:{3:.1f}h, HA:{4:.1f}h, DEC:{5:.1f}d)'.format(timstr, y, srcname, ra0, ha0, dec0)
            else:
                flux=moreflux[ind-ntp]
            return '{2} ({6:.1f}Jy RA:{3:.1f}h, HA:{4:.1f}h, DEC:{5:.1f}d)'.format(timstr, y, srcname, ra0, ha0, dec0, flux)
        else:
            return 'x = {0:.1f}, y = {1:.1f}'.format(x,y)
#    ax.format_coord = format_coord

    plt.xticks(rotation=40)
    box = ax.get_position()
    ax.set_position([box.x0-0.05, box.y0+0.15, box.width * 0.82, box.height * 0.85])
    ax.legend(srclegend,loc='center left', bbox_to_anchor=(1., 0.5),
                      fancybox=True, shadow=True)
    ax.set_xlabel('UTC Time')
    ax.set_ylabel('Altitude (deg)')
    ax.set_title('Bright sources visible by EOVSA 27-m Antenna')
    txt = tt.tooltip(f, ax, format_coord, ts, srclist, alt)
