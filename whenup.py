# Purpose: Contains a bunch of routines for finding and plotting when
#          calibration sources and the Sun are up.  Also contains the
#          key routine make_sched(), which creates a solar schedule
#          for a given date.
# History:
#  2018-09-16  DG
#    First wrote whenup(), based on whatup.py
#  2018-09-18  DG
#    Completed the routines whenup(), sunup(), plot_sun(), and make_sched()
#  2019-01-18  DG
#    A bug occurred on some dates due to source not being within 0.1 degree of
#    10 degrees altitude.  Introduced a for loop and test in both whenup() and
#    make_sched() to use a wider window.
#  2019-05-20  DG
#    Added a 1-min PHASECAL with sequence solar.fsq just before the SKYCALTEST.
#  2020-10-01  DG
#    Added a check for day numbers 259-287, when 3C273 (1229+020) is too close
#    to the Sun.  On those dates, replace any 1229+020 lines with 1331+305.
#  2020-10-29  DG
#    Comment out a bunch of lines due to 27-m not working--lines are commented
#    out using #**.  Note that two lines were added that have to be removed,
#    and they also have #** in a comment on those lines.
#  2020-10-31  DG
#    Looks like the 27-m is back, so I reverted the code back to the original.
#  2021-09-09  DG
#    The schedule broke today because it is day 251, which is a transition between
#    three and two calibrations during the day.  I changed the iday range from >251
#    to >=251, and now it seems to work.
#  2021-09-22  DG
#    My replacement of 1229+020 by 1331+305 for day numbers 259-287 was broken, 
#    because 1331+305 was not up yet!  I now add an appropriate number of minutes 
#    to the times of those lines.
#  2022-03-14  DG
#    Changes to reduce the length of PHASECALs to from 25-30 minutes to 20 minutes.
#  2022-05-14  DG
#    Added a remove_cal() function to remove all 27-m calibration lines from
#    the solar schedule created by make_sched().

import os
from util import Time
import eovsa_cat
from eovsa_visibility import *
import numpy as np
from astropy.time import TimeDelta

def deg(rad):
    return (rad * 180./np.pi + 180) % 360 - 180

def deg2(rad):
    return rad * 180./np.pi

def whenup(date=None,verbose=False):
    ''' Find out the times when preferred sources are up for a given date. 
        Prints the result of source rise and set times for preferred sources:
        Sun, 0319+415, 1229+020, 1331+305, 2136+006, 2253+161
        
        Optional arguments:
           date: Time() object giving start date (time of day is ignored) 
                 or None for current date.
           verbose: If True, actually print the table of times to the screen
           
        Returns: dictionary of time objects with keys 'taz', 'teq', 'tgap', each
                 being a two-element Time() object, where first time is
                 start (rise) and second is end (set), and 'lines', which is a
                 printable table of times.
    '''
    srclist = ['Sun     ', '0319+415', '1229+020', '1331+305', '2136+006', '2253+161']

    # Use the 24-h day specified by date (i.e. drop the time of day)
    if date is None:
        # Use today's date
        mjd = int(Time.now().mjd)
    else:
        mjd = int(date.mjd)
    t = Time(mjd,format='mjd')

    # Add times in 1-min steps up to duration dur (hours)
    ts = t + TimeDelta(np.arange(0.,24.,1./60.)/24.,format='jd')

    aa = eovsa_cat.eovsa_array_with_cat()

    nt = len(ts)
    ra = np.zeros((len(srclist),nt))
    ha = np.zeros((len(srclist),nt))
    dec = np.zeros((len(srclist),nt))
    alt = np.zeros((len(srclist),nt))
    az = np.zeros((len(srclist),nt))
    taz = []
    teq = []
    tgap = []
    lines = []
    for i in range(nt):
        aa.set_jultime(ts[i].jd)    
        lst = aa.sidereal_time()

        for j,srcname in enumerate(srclist):
            src = aa.cat[srcname.split()[0]]
            src.compute(aa)
            ra[j,i] = src.ra
            ha[j,i] = lst - src.ra
            dec[j,i] = src.dec
            alt[j,i] = src.alt
            az[j,i] = src.az

    ra_deg = deg(ra)
    ha_deg = deg(ha)
    dec_deg = deg(dec)
    alt = deg(alt)
    az_deg = deg2(az)
    lines.append('Source     Alt-Az    Equatorial     Gap')
    lines.append(' Name    Rise   Set  Rise   Set  Start  End')
    lines.append('-------- ----- ----- ----- ----- ----- -----')
    for j in range(len(srclist)):
        iset_eq = np.where(np.abs(ha_deg[j] - 55.0) < 0.15)[0][0]
        irise_eq = np.where(np.abs(ha_deg[j] + 55.0) < 0.15)[0][0]
        trise_eq = ts[irise_eq]
        tset_eq = ts[iset_eq]
        if iset_eq < irise_eq:
            if j == 0:
                tset_eq += TimeDelta(1,format='jd')
            else:
                tset_eq += TimeDelta(1436./60./24.,format='jd')
        for iwindow in np.arange(0.1,0.2,0.01):
            try:
                i1 = np.where(np.abs(alt[j] - 10.0) < iwindow)[0][0]
                i2 = np.where(np.abs(alt[j] - 10.0) < iwindow)[0][1]
                if i2 == i1+1:
                    i2 = np.where(np.abs(alt[j] - 10.0) < iwindow)[0][2]
                break
            except:
                print 'Window',iwindow,'did not work.  Trying again'
        if alt[j,i1] < alt[j,i1+1]:
            irise_az = i1
            iset_az = i2
        else:
            irise_az = i2
            iset_az = i1
        trise_az = ts[irise_az]
        tset_az = ts[iset_az]
        if iset_az < irise_az:
            if j == 0:
                tset_az += TimeDelta(1,format='jd')
            else:
                tset_az += TimeDelta(1436./60./24.,format='jd')
        if j == 1:
            up, = np.where(alt[j] > 10.)
            j1 = np.where(np.abs(az_deg[j,up] - 35.0) < 1.0)[0][0]
            j2 = np.where(np.abs(az_deg[j,up] - 325.0) < 1.0)[0][0]
            jset = up[j1]
            jrise = up[j2]
            tset_gap = ts[jset]
            trise_gap = ts[jrise]
            if jrise < jset:
                trise_gap += TimeDelta(1436./60./24.,format='jd')
            taz.append(Time([trise_az.iso,tset_az.iso]))
            teq.append(Time([trise_eq.iso,tset_eq.iso]))
            tgap.append(Time([tset_gap.iso,trise_gap.iso]))
            lines.append('{0} {1} {2} {3} {4} {5} {6}'.format(srclist[j], trise_az.iso[11:16], tset_az.iso[11:16], trise_eq.iso[11:16], tset_eq.iso[11:16], tset_gap.iso[11:16], trise_gap.iso[11:16]))
        else:
            taz.append(Time([trise_az.iso,tset_az.iso]))
            teq.append(Time([trise_eq.iso,tset_eq.iso]))
            tgap.append((None,None))
            lines.append('{0} {1} {2} {3} {4} {5} {6}'.format(srclist[j], trise_az.iso[11:16], tset_az.iso[11:16], trise_eq.iso[11:16], tset_eq.iso[11:16], ' --- ',' --- '))
    if verbose:
        for line in lines:
            print line
    return {'source':srclist,'taz':taz, 'teq':teq, 'tgap':tgap, 'lines':lines}

def sunup(daterange):
    ''' Find out the times when the Sun is up for a given date range. 
           
        Returns: dictionary of time objects with keys 'taz', 'teq', each
                 being a two-element Time() object, where first time is
                 start (rise) and second is end (set)
    '''
    # Use the 24-h day specified by daterange (i.e. drop the time of day)
    mjd1 = int(daterange[0].mjd)
    mjd2 = int(daterange[1].mjd)
    aa = eovsa_cat.eovsa_array_with_cat()
    taz_rise = []
    teq_rise = []
    taz_set = []
    teq_set = []
    mjd_list = []
    for mjd in range(mjd1,mjd2+1):
        t = Time(mjd,format='mjd')
        # Add times in 1-min steps up to duration dur (hours)
        ts = t + TimeDelta(np.arange(0.,24.,1./60.)/24.,format='jd')

        nt = len(ts)
        ha = np.zeros(nt)
        alt = np.zeros(nt)
        az = np.zeros(nt)
        lines = []
        for i in range(nt):
            aa.set_jultime(ts[i].jd)    
            lst = aa.sidereal_time()

            src = aa.cat['Sun']
            src.compute(aa)
            ha[i] = lst - src.ra
            alt[i] = src.alt
            az[i] = src.az

        ha_deg = deg(ha)
        alt = deg(alt)
        az_deg = deg2(az)
        iset_eq = np.where(np.abs(ha_deg - 55.0) < 0.15)[0][0]
        irise_eq = np.where(np.abs(ha_deg + 55.0) < 0.15)[0][0]
        trise_eq = ts[irise_eq]
        tset_eq = ts[iset_eq]
        if iset_eq < irise_eq:
            tset_eq += TimeDelta(1,format='jd')
        for iwindow in np.arange(0.1,0.2,0.01):
            try:
                idx, = np.where(np.abs(alt - 10.0) < iwindow)  #****Was 0.1
                i1 = idx[0]
                i2 = idx[1]
                if i2 == i1+1:
                    if len(idx) > 2:
                        i2 = idx[2]
                    else:
                        i2 = 1439  # Transition is at end
                break
            except:
                print 'Window',iwindow,'did not work.  Trying again'
        if alt[i1] < alt[i1+1]:
            irise_az = i1
            iset_az = i2
        else:
            irise_az = i2
            iset_az = i1
        trise_az = ts[irise_az]
        tset_az = ts[iset_az]
        if iset_az < irise_az:
            tset_az += TimeDelta(1,format='jd')
        taz_rise.append(trise_az.mjd)
        taz_set.append(tset_az.mjd)
        teq_rise.append(trise_eq.mjd)
        teq_set.append(tset_eq.mjd)
        mjd_list.append(mjd)
    return {'date':Time(mjd_list,format='mjd'),
            'taz_rise':Time(taz_rise,format='mjd'),'teq_rise':Time(teq_rise,format='mjd'),
             'taz_set':Time(taz_set, format='mjd'), 'teq_set':Time(teq_set, format='mjd')}


def plot_sun(sun):
    ''' Plot the output of sunup() in a nice format for visualizing the
        relevant times.
        
        Returns:
           ax    The axis of the plot, for potential use by make_sched for overplotting
    '''
    import matplotlib.pylab as plt
    t = np.array([0.0, 0.141, 0.1639, 0.3049, 0.3805, 0.6847, 0.7278, 0.7625, 0.816, 1.1188])
    mjd0 = sun['date'][0].mjd
    mjd1 = sun['date'][-1].mjd
    y = np.array([mjd0, mjd0, mjd1, mjd1])
    out = whenup(sun['date'][0])
    t0 = out['teq'][1][0].mjd % 1
    f, ax = plt.subplots(1,1)
    ax.plot(sun['taz_rise'].mjd-sun['date'].mjd,sun['date'].mjd,'--',color='C1')
    ax.plot(sun['taz_set'].mjd-sun['date'].mjd,sun['date'].mjd,'--',color='C1')
    ax.plot(sun['taz_rise'].mjd-sun['date'].mjd-0.0583,sun['date'].mjd,':',color='C1')
    ax.plot(sun['taz_set'].mjd-sun['date'].mjd+0.0583,sun['date'].mjd,':',color='C1')
    ax.plot(sun['teq_rise'].mjd-sun['date'].mjd,sun['date'].mjd,'-',color='C0')
    ax.plot(sun['teq_set'].mjd-sun['date'].mjd,sun['date'].mjd,'-',color='C0')
    colors = ['C3','C3','C4','C5','C6','C7']
    for k,i in enumerate([0,2,4,5,7,8]):
        x = np.array([t0 + t[i], t0+t[i+1], (mjd0-mjd1)*0.0027378+t0+t[i+1],(mjd0-mjd1)*0.0027378+t0+t[i]])
        ax.fill(x,y,color=colors[k],alpha=0.25)
        ax.fill(x+1,y,color=colors[k],alpha=0.25)
        ax.fill(x-1,y,color=colors[k],alpha=0.25)
    ax.plot(np.ones(2)*0.7708,[mjd0,mjd1],color='black') # 18:30
    ax.plot(np.ones(2)*0.8958,[mjd0,mjd1],color='black') # 21:30
    ax.plot(np.ones(2)*0.6319,[mjd0,mjd1],color='gray') # 15:10
    ax.plot(np.ones(2)*0.6563,[mjd0,mjd1],color='gray') # 15:45
    ax.plot(np.ones(2)*0.7153,[mjd0,mjd1],color='gray') # 17:10
    ax.plot(np.ones(2)*0.7396,[mjd0,mjd1],color='gray') # 17:45
    ax.plot(np.ones(2)*0.8264,[mjd0,mjd1],color='gray') # 19:50
    ax.plot(np.ones(2)*0.8507,[mjd0,mjd1],color='gray') # 20:25
    ax.plot(np.ones(2)*0.9236,[mjd0,mjd1],color='gray') # 22:10
    ax.plot(np.ones(2)*0.9479,[mjd0,mjd1],color='gray') # 22:45
    ax.plot(np.ones(2)*1.0069,[mjd0,mjd1],color='gray') # 00:10
    ax.plot(np.ones(2)*1.0313,[mjd0,mjd1],color='gray') # 00:45
    ax.set_xlim(0.0,1.5)
    ax.set_title(sun['date'][0].iso[:10]+' to '+sun['date'][-1].iso[:10])
    ax.set_ylabel('MJD')
    ax.set_xlabel('Time [in fraction of a day]')
    return ax

def make_sched(sun=None, t=None, ax=None, verbose=False):
    ''' Main routine to create a daily schedule, given at date as a Time() object.
        If called with no arguments, a schedule for the current day is returned.
        Optional arguments:
           sun    A "sun" dictionary returned by sunup(), with the restriction
                    that the daterange used in the sunup() call must contain
                    the date in t.  This is a time-saving measure in case many 
                    schedules are being generated for a given year, since the 
                    internal call to sunup() is skipped.  If omitted, sunup() is
                    called for the given date.
           t      A Time() object giving the date for which the schedule is
                    desired.  If omitted, the current date is used.
           ax     An axis object for an existing plot.  This is mainly a debugging
                    tool.  If ax = plot_sun() is called, then make_sched(ax=ax) 
                    will cause blue lines to be overplotted for the given date, 
                    showing the time ranges chosen for the calibrators.
          verbose If True, will also print the schedule lines to the screen.
          
        Returns:
          lines   A list of text lines representing the schedule.
    '''
    # From give time, get mjd0 = MJD of Jan 1 for that year
    if t is None:
        t = Time.now()
    imjd = int(t.mjd)
    year = t.iso[:4]
    mjd0 = int(Time(year+'-01-01').mjd)
    # If no sun dictionary is give, create a 2-day one
    if sun is None:
        sun = sunup(Time([t.mjd,t.mjd+1],format='mjd'))
    # Calibration durations, minutes
    refdur = 84.
    caldur = 20.  #35.
    # Get calibrator reference time.  Rather than calculating exact timing for
    # calibrators for each day, we just get it at one reference time and then
    # calculate it for additional days by subtracting 0.0027387*dday
    out = whenup(sun['date'][0])
    t0 = out['teq'][1][0].mjd % 1   # Reference time for all calibrators
    # Calibrator windows, as fraction of a day
    ts = np.array([0.0000, 0.1639, 0.3805, 0.6847, 0.7625, 0.8160])+t0
    te = np.array([0.1410, 0.3049, 0.6847, 0.7278, 0.8160, 1.1188])+t0
    # Day numbers of source transitions (this will shift for leap-years,
    # but it is not critical, so that is ignored.
    pcal2trans = np.array([[0,20,45,83,251,366],[0,75,83,251,291,309,330,366]])
    pcal3trans = np.array([[0,83,182,199,251,366],[0,83,111,128,197,251,366],[0,83,132,251,366]])
    refcaltrans = np.array([[0,121,210,242,305,366],[0,31,37,84,210,366]])
    # Nominal phasecal scan start times
    pcal2s = np.array([(17. + 10./60.)/24.,(22. + 10./60.)/24.])  # 17:10 and 22:10 UT
    pcal3s = np.array([(15. + 10./60.)/24.,(19. + 50./60.)/24.,(24. + 10./60.)/24.]) # 15:10, 19:50 and 24:10 UT
    # Source designations for each range
    pcal2srcs  = [[3,4,5,-1,2],[5,0,-1,2,3,4,5]] 
    pcal3srcs  = [[-1,5,0,1,-1],[-1,5,0,1,2,-1],[-1,1,2,-1]]
    refcalsrcs = [[2,5,0,1,2],[5,0,1,2,5]]
    calnames = ['0319+415','0319+415', '1229+020', '1331+305', '2136+006', '2253+161']
#    for imjd = range(mjd0,mjd1+1):
    iday = int(imjd - mjd0)
    nday = imjd - int(sun['date'][0].mjd)   # Number of days into sun dictionary corresponding to this day
    # For each day, identify the sources
    #   Morning refcal
    refcal1 = refcalsrcs[0][np.where(iday - np.array(refcaltrans[0]) >= 0.0)[0][-1]]
    refset = (te[refcal1]-nday*0.0027378 + 1.0) % 1
    sunrise = sun['taz_rise'][nday].mjd % 1
    rc1end = min(refset,sunrise)
    rc1start = rc1end - refdur/1440.
    lines = []
    lines.append('{:} {:} {:}'.format(Time(imjd + rc1start,format='mjd').iso[:19],'ACQUIRE',calnames[refcal1]))
    lines.append('{:} {:}'.format(Time(imjd + rc1start + 1./1440.,format='mjd').iso[:19],'LOSELECT'))
    lines.append('{:} {:} {:} {:}'.format(Time(imjd + rc1start + 4./1440.,format='mjd').iso[:19],'PHASECAL_LO',calnames[refcal1],'pcal_lo.fsq'))
    lines.append('{:} {:}'.format(Time(imjd + rc1start + 24./1440.,format='mjd').iso[:19],'HISELECT'))
    lines.append('{:} {:} {:} {:}'.format(Time(imjd + rc1start + 25./1440.,format='mjd').iso[:19],'PHASECAL',calnames[refcal1],'pcal_hi-all.fsq'))
    if verbose:
        print Time(imjd + rc1start,format='mjd').iso[:19],'ACQUIRE',calnames[refcal1]
        print Time(imjd + rc1start + 1./1440.,format='mjd').iso[:19],'LOSELECT'
        print Time(imjd + rc1start + 4./1440.,format='mjd').iso[:19],'PHASECAL_LO',calnames[refcal1],'pcal_lo.fsq'
        print Time(imjd + rc1start + 24./1440.,format='mjd').iso[:19],'HISELECT'
        print Time(imjd + rc1start + 25./1440.,format='mjd').iso[:19],'PHASECAL',calnames[refcal1],'pcal_hi-all.fsq'
    if rc1end != sunrise:
        lines.append('{:} {:}'.format(Time(imjd + rc1start + 85./1440.,format='mjd').iso[:19],'STOW'))
        if verbose: print Time(imjd + rc1start + 85./1440.,format='mjd').iso[:19],'STOW'
    lines.append('{:} {:}'.format(Time(imjd + sunrise,format='mjd').iso[:19],'SUN'))
    if verbose: print Time(imjd + sunrise,format='mjd').iso[:19],'SUN'
    if ax:
        ax.plot([rc1start,rc1start+refdur/1440.],[imjd,imjd],color='C0',alpha=0.25)
    #   Phasecals
    if iday < 83 or iday >= 251:
        #import pdb; pdb.set_trace()
        pcal1 = pcal2srcs[0][np.where(iday - np.array(pcal2trans[0]) >= 0.0)[0][-1]]
        pc1rise = (ts[pcal1]-nday*0.0027378 + 1.0) % 1
        pc1set = (te[pcal1]-nday*0.0027378 + 1.0) % 1
        if pc1set < 0.5: pc1set += 1.0
        #print Time(pc1rise+imjd,format='mjd').iso, pc1set
        pc1start = max(pc1rise,pcal2s[0])
        #print Time(pc1start+imjd,format='mjd').iso
        pc1start = min(pc1start,pc1set-caldur/1440.)
        #print Time(pc1start+imjd,format='mjd').iso,Time(pc1set+imjd,format='mjd').iso
        #print pc1rise, pc1set, Time(imjd + pc1rise,format='mjd').iso,Time(imjd + pc1set,format='mjd').iso,
        pcal2 = pcal2srcs[1][np.where(iday - np.array(pcal2trans[1]) >= 0.0)[0][-1]]
        pc2rise = (ts[pcal2]-nday*0.0027378 + 1.0) % 1
        pc2set = (te[pcal2]-nday*0.0027378 + 1.0) % 1
        if pc2set < 0.5: pc2set += 1.0
        pc2start = max(pc2rise,pcal2s[1])
        pc2start = min(pc2start,pc2set-caldur/1440.)
        lines.append('{:} {:} {:}'.format(Time(imjd + pc1start,format='mjd').iso[:19],'ACQUIRE',calnames[pcal1]))
        lines.append('{:} {:} {:} {:}'.format(Time(imjd + pc1start + 4./1440.,format='mjd').iso[:19],'PHASECAL',calnames[pcal1],'pcal_hi-all.fsq'))
        lines.append('{:} {:} {:} {:}'.format(Time(imjd + pc1start + 20./1440.,format='mjd').iso[:19],'PHASECAL',calnames[pcal1],'solar.fsq'))
        lines.append('{:} {:} {:}'.format(Time(imjd + pc1start + 21./1440.,format='mjd').iso[:19],'SKYCALTEST',calnames[pcal1]))
        lines.append('{:} {:}'.format(Time(imjd + pc1start + 25./1440.,format='mjd').iso[:19],'SUN'))
#        lines.append('{:} {:} {:} {:}'.format(Time(imjd + pc1start + 30./1440.,format='mjd').iso[:19],'PHASECAL',calnames[pcal1],'solar.fsq'))
#        lines.append('{:} {:} {:}'.format(Time(imjd + pc1start + 31./1440.,format='mjd').iso[:19],'SKYCALTEST',calnames[pcal1]))
#        lines.append('{:} {:}'.format(Time(imjd + pc1start + 35./1440.,format='mjd').iso[:19],'SUN'))
        lines.append('{:} {:}'.format(Time(imjd + (18.*60. + 30.)/1440.,format='mjd').iso[:19],'SOLPNTCAL solar.fsq solpnt.trj'))
        lines.append('{:} {:}'.format(Time(imjd + (18.*60. + 35.)/1440.,format='mjd').iso[:19],'SUN'))
        lines.append('{:} {:}'.format(Time(imjd + (20.*60. + 00.)/1440.,format='mjd').iso[:19],'GAINCALTEST'))
        lines.append('{:} {:}'.format(Time(imjd + (20.*60. + 03.)/1440.,format='mjd').iso[:19],'SUN'))
        lines.append('{:} {:}'.format(Time(imjd + (21.*60. + 30.)/1440.,format='mjd').iso[:19],'SOLPNTCAL solar.fsq solpnt.trj'))
        lines.append('{:} {:}'.format(Time(imjd + (21.*60. + 35.)/1440.,format='mjd').iso[:19],'SUN'))
        lines.append('{:} {:} {:}'.format(Time(imjd + pc2start,format='mjd').iso[:19],'ACQUIRE',calnames[pcal2]))
        lines.append('{:} {:} {:} {:}'.format(Time(imjd + pc2start + 4./1440.,format='mjd').iso[:19],'PHASECAL',calnames[pcal2],'pcal_hi-all.fsq'))
        lines.append('{:} {:}'.format(Time(imjd + pc2start + 20./1440.,format='mjd').iso[:19],'SUN'))
#        lines.append('{:} {:}'.format(Time(imjd + pc2start + 35./1440.,format='mjd').iso[:19],'SUN'))
        if verbose:
            print Time(imjd + pc1start,format='mjd').iso[:19],'ACQUIRE',calnames[pcal1]
            print Time(imjd + pc1start + 4./1440.,format='mjd').iso[:19],'PHASECAL',calnames[pcal1],'pcal_hi-all.fsq'
            print Time(imjd + pc1start + 20./1440.,format='mjd').iso[:19],'PHASECAL',calnames[pcal1],'solar.fsq'
            print Time(imjd + pc1start + 21./1440.,format='mjd').iso[:19],'SKYCALTEST',calnames[pcal1]
            print Time(imjd + pc1start + 25./1440.,format='mjd').iso[:19],'SUN'
#            print Time(imjd + pc1start + 30./1440.,format='mjd').iso[:19],'PHASECAL',calnames[pcal1],'solar.fsq'
#            print Time(imjd + pc1start + 31./1440.,format='mjd').iso[:19],'SKYCALTEST',calnames[pcal1]
#            print Time(imjd + pc1start + 35./1440.,format='mjd').iso[:19],'SUN'
            print Time(imjd + (18.*60. + 30.)/1440.,format='mjd').iso[:19],'SOLPNTCAL solar.fsq solpnt.trj'
            print Time(imjd + (18.*60. + 35.)/1440.,format='mjd').iso[:19],'SUN'
            print Time(imjd + (20.*60. + 00.)/1440.,format='mjd').iso[:19],'GAINCALTEST'          
            print Time(imjd + (20.*60. + 03.)/1440.,format='mjd').iso[:19],'SUN'
            print Time(imjd + (21.*60. + 30.)/1440.,format='mjd').iso[:19],'SOLPNTCAL solar.fsq solpnt.trj'
            print Time(imjd + (21.*60. + 35.)/1440.,format='mjd').iso[:19],'SUN'
            print Time(imjd + pc2start,format='mjd').iso[:19],'ACQUIRE',calnames[pcal2]
            print Time(imjd + pc2start + 4./1440.,format='mjd').iso[:19],'PHASECAL',calnames[pcal2],'pcal_hi-all.fsq'
            print Time(imjd + pc2start + 20./1440.,format='mjd').iso[:19],'SUN'
#            print Time(imjd + pc2start + 35./1440.,format='mjd').iso[:19],'SUN'
        if ax:
            ax.plot([pc1start,pc1start+caldur/1440.],[imjd,imjd],color='C0',alpha=0.25)
            ax.plot([pc2start,pc2start+caldur/1440.],[imjd,imjd],color='C0',alpha=0.25)
    else:
        pcal1 = pcal3srcs[0][np.where(iday - np.array(pcal3trans[0]) >= 0.0)[0][-1]]
        pc1rise = (ts[pcal1]-nday*0.0027378 + 1.0) % 1
        pc1set = (te[pcal1]-nday*0.0027378 + 1.0) % 1
        pc1start = max(pc1rise,pcal3s[0])
        pc1start = min(pc1start,pc1set-caldur/1440.)
        pcal2 = pcal3srcs[1][np.where(iday - np.array(pcal3trans[1]) >= 0.0)[0][-1]]
        pc2rise = (ts[pcal2]-nday*0.0027378 + 1.0) % 1
        pc2set = (te[pcal2]-nday*0.0027378 + 1.0) % 1
        if pc2set < 0.5: pc2set += 1.0
        pc2start = max(pc2rise,pcal3s[1])
        pc2start = min(pc2start,pc2set-caldur/1440.)
        pcal3 = pcal3srcs[2][np.where(iday - np.array(pcal3trans[2]) >= 0.0)[0][-1]]
        pc3rise = (ts[pcal3]-nday*0.0027378 + 1.0) % 1
        if pc3rise < 0.5: pc3rise += 1.0
        pc3set = (te[pcal3]-nday*0.0027378 + 1.0) % 1
        if pc3set < 0.5: pc3set += 1.0
        pc3start = max(pc3rise,pcal3s[2])
        pc3start = min(pc3start,pc3set-caldur/1440.)

        lines.append('{:} {:} {:}'.format(Time(imjd + pc1start,format='mjd').iso[:19],'ACQUIRE',calnames[pcal1]))
        lines.append('{:} {:} {:} {:}'.format(Time(imjd + pc1start + 4./1440.,format='mjd').iso[:19],'PHASECAL',calnames[pcal1],'pcal_hi-all.fsq'))
        lines.append('{:} {:} {:} {:}'.format(Time(imjd + pc1start + 20./1440.,format='mjd').iso[:19],'PHASECAL',calnames[pcal1],'solar.fsq'))
        lines.append('{:} {:} {:}'.format(Time(imjd + pc1start + 21./1440.,format='mjd').iso[:19],'SKYCALTEST',calnames[pcal1]))
        lines.append('{:} {:}'.format(Time(imjd + pc1start + 25./1440.,format='mjd').iso[:19],'SUN'))
#        lines.append('{:} {:} {:} {:}'.format(Time(imjd + pc1start + 30./1440.,format='mjd').iso[:19],'PHASECAL',calnames[pcal1],'solar.fsq'))
#        lines.append('{:} {:} {:}'.format(Time(imjd + pc1start + 31./1440.,format='mjd').iso[:19],'SKYCALTEST',calnames[pcal1]))
#        lines.append('{:} {:}'.format(Time(imjd + pc1start + 35./1440.,format='mjd').iso[:19],'SUN'))
        lines.append('{:} {:}'.format(Time(imjd + (17.*60. + 00.)/1440.,format='mjd').iso[:19],'GAINCALTEST'))
        lines.append('{:} {:}'.format(Time(imjd + (17.*60. + 03.)/1440.,format='mjd').iso[:19],'SUN'))
        lines.append('{:} {:}'.format(Time(imjd + (18.*60. + 30.)/1440.,format='mjd').iso[:19],'SOLPNTCAL solar.fsq solpnt.trj'))
        lines.append('{:} {:}'.format(Time(imjd + (18.*60. + 35.)/1440.,format='mjd').iso[:19],'SUN'))
        lines.append('{:} {:} {:}'.format(Time(imjd + pc2start,format='mjd').iso[:19],'ACQUIRE',calnames[pcal2]))
        lines.append('{:} {:} {:} {:}'.format(Time(imjd + pc2start + 4./1440.,format='mjd').iso[:19],'PHASECAL',calnames[pcal2],'pcal_hi-all.fsq'))
        lines.append('{:} {:}'.format(Time(imjd + pc2start + 20./1440.,format='mjd').iso[:19],'SUN'))
#        lines.append('{:} {:}'.format(Time(imjd + pc2start + 35./1440.,format='mjd').iso[:19],'SUN'))
        lines.append('{:} {:}'.format(Time(imjd + (21.*60. + 30.)/1440.,format='mjd').iso[:19],'SOLPNTCAL solar.fsq solpnt.trj'))
        lines.append('{:} {:}'.format(Time(imjd + (21.*60. + 35.)/1440.,format='mjd').iso[:19],'SUN'))
        lines.append('{:} {:} {:}'.format(Time(imjd + pc3start,format='mjd').iso[:19],'ACQUIRE',calnames[pcal3]))
        lines.append('{:} {:} {:} {:}'.format(Time(imjd + pc3start + 4./1440.,format='mjd').iso[:19],'PHASECAL',calnames[pcal3],'pcal_hi-all.fsq'))
        lines.append('{:} {:}'.format(Time(imjd + pc3start + 20./1440.,format='mjd').iso[:19],'SUN'))
#        lines.append('{:} {:}'.format(Time(imjd + pc3start + 35./1440.,format='mjd').iso[:19],'SUN'))
        if verbose:
            print Time(imjd + pc1start,format='mjd').iso[:19],'ACQUIRE',calnames[pcal1]
            print Time(imjd + pc1start + 4./1440.,format='mjd').iso[:19],'PHASECAL',calnames[pcal1],'pcal_hi-all.fsq'
            print Time(imjd + pc1start + 20./1440.,format='mjd').iso[:19],'PHASECAL',calnames[pcal1],'solar.fsq'
            print Time(imjd + pc1start + 21./1440.,format='mjd').iso[:19],'SKYCALTEST',calnames[pcal1]
            print Time(imjd + pc1start + 25./1440.,format='mjd').iso[:19],'SUN'
#            print Time(imjd + pc1start + 30./1440.,format='mjd').iso[:19],'PHASECAL',calnames[pcal1],'solar.fsq'
#            print Time(imjd + pc1start + 31./1440.,format='mjd').iso[:19],'SKYCALTEST',calnames[pcal1]
#            print Time(imjd + pc1start + 35./1440.,format='mjd').iso[:19],'SUN'
            print Time(imjd + (17.*60. + 00.)/1440.,format='mjd').iso[:19],'GAINCALTEST'          
            print Time(imjd + (17.*60. + 03.)/1440.,format='mjd').iso[:19],'SUN'
            print Time(imjd + (18.*60. + 30.)/1440.,format='mjd').iso[:19],'SOLPNTCAL solar.fsq solpnt.trj'
            print Time(imjd + (18.*60. + 35.)/1440.,format='mjd').iso[:19],'SUN'
            print Time(imjd + pc2start,format='mjd').iso[:19],'ACQUIRE',calnames[pcal2]
            print Time(imjd + pc2start + 4./1440.,format='mjd').iso[:19],'PHASECAL',calnames[pcal2],'pcal_hi-all.fsq'
            print Time(imjd + pc2start + 20./1440.,format='mjd').iso[:19],'SUN'
#            print Time(imjd + pc2start + 35./1440.,format='mjd').iso[:19],'SUN'
            print Time(imjd + (21.*60. + 30.)/1440.,format='mjd').iso[:19],'SOLPNTCAL solar.fsq solpnt.trj'
            print Time(imjd + (21.*60. + 35.)/1440.,format='mjd').iso[:19],'SUN'
            print Time(imjd + pc3start,format='mjd').iso[:19],'ACQUIRE',calnames[pcal3]
            print Time(imjd + pc3start + 4./1440.,format='mjd').iso[:19],'PHASECAL',calnames[pcal3],'pcal_hi-all.fsq'
            print Time(imjd + pc3start + 20./1440.,format='mjd').iso[:19],'SUN'
#            print Time(imjd + pc3start + 35./1440.,format='mjd').iso[:19],'SUN'
        if ax:
            ax.plot([pc1start,pc1start+caldur/1440.],[imjd,imjd],color='C0',alpha=0.25)
            ax.plot([pc2start,pc2start+caldur/1440.],[imjd,imjd],color='C0',alpha=0.25)
            ax.plot([pc3start,pc3start+caldur/1440.],[imjd,imjd],color='C0',alpha=0.25)
    #   Evening refcal
    refcal2 = refcalsrcs[1][np.where(iday - np.array(refcaltrans[1]) >= 0.0)[0][-1]]
    refrise = (ts[refcal2]-nday*0.0027378 + 1.0) % 1
    if refrise < 0.5: refrise += 1.0
    sunset = sun['taz_set'][nday].mjd % 1
    if sunset < 0.5: sunset += 1.0
    rc2start = max(refrise,sunset)
#        rc2end = rc2start + refdur/1440.
    if refrise > sunset:
        lines.append('{:} {:}'.format(Time(imjd + sunset,format='mjd').iso[:19],'STOW'))
        if verbose: print Time(imjd + sunset,format='mjd').iso[:19],'STOW'
    lines.append('{:} {:} {:}'.format(Time(imjd + rc2start,format='mjd').iso[:19],'ACQUIRE',calnames[refcal2]))
    lines.append('{:} {:}'.format(Time(imjd + rc2start + 1./1440.,format='mjd').iso[:19],'LOSELECT'))
    lines.append('{:} {:} {:} {:}'.format(Time(imjd + rc2start + 4./1440.,format='mjd').iso[:19],'PHASECAL_LO',calnames[refcal2],'pcal_lo.fsq'))
    lines.append('{:} {:}'.format(Time(imjd + rc2start + 24./1440.,format='mjd').iso[:19],'HISELECT'))
    lines.append('{:} {:} {:} {:}'.format(Time(imjd + rc2start + 25./1440.,format='mjd').iso[:19],'PHASECAL',calnames[refcal2],'pcal_hi-all.fsq'))
    lines.append('{:} {:}'.format(Time(imjd + rc2start + 85./1440.,format='mjd').iso[:19],'STOW'))
    if verbose:
        print Time(imjd + rc2start,format='mjd').iso[:19],'ACQUIRE',calnames[refcal2]
        print Time(imjd + rc2start + 1./1440.,format='mjd').iso[:19],'LOSELECT'
        print Time(imjd + rc2start + 4./1440.,format='mjd').iso[:19],'PHASECAL_LO',calnames[refcal2],'pcal_lo.fsq'
        print Time(imjd + rc2start + 24./1440.,format='mjd').iso[:19],'HISELECT'
        print Time(imjd + rc2start + 25./1440.,format='mjd').iso[:19],'PHASECAL',calnames[refcal2],'pcal_hi-all.fsq'
        print Time(imjd + rc2start + 85./1440.,format='mjd').iso[:19],'STOW'
    if ax:
        ax.plot([rc2start,rc2start+refdur/1440.],[imjd,imjd],color='C0',alpha=0.25)
    lines = chk_sched(lines)
    return lines

def chk_sched(lines):
    # Post-creation check on the schedule to address the fact that
    # 3C273 is too close to the Sun from 9/15 to 10/13 each year (day-of-year 259-287).
    doy = int(Time(lines[0][:10]).yday[5:8])
    next = -1  # Impossible line number
    if doy >= 259 and doy <= 287:
        # This is a date range when the 27-m antenna should not point at source 1229+020
        # Simply replace 1229+020 with 1331+305 (3C286).  That is a much weaker source, but it
        # may be okay for a PHASECAL
        if doy < 272:
            dt = (3000. - (doy - 259)*240.)/86400.  # Starts at 50 min and decreases by 4 min each day
        else:
            dt = 0.
        for i in range(len(lines)):
            if lines[i].find('1229+020') != -1:
                # This line has to be adjusted by dt later and change source to '1331+305'
                lines[i] = lines[i].replace('1229+020','1331+305')
                lines[i] = (Time(lines[i][:19])+dt).iso[:19]+lines[i][19:]
                next = i+1  # Set to next line number, which also has to increment by dt
            elif next == i:
                lines[i] = (Time(lines[i][:19])+dt).iso[:19]+lines[i][19:]
                next = -1
    return lines

def remove_cal(lines):
    from util import Time
    keepidx = []
    sunidx = []
    rmidx = []
    for i,line in enumerate(lines):
        if line.find('SUN') > 0:
            keepidx.append(i)
            keepidx.append(i+1)
        if line.find('SKYCAL') > 0:
            keepidx.append(i-1)
            keepidx.append(i)
    keepixd = np.array(keepidx)
    keeplines = np.array(lines)
    keeplines = keeplines[keepidx]
    for i,line in enumerate(keeplines):
        if line.find('SUN') > 0:
            sunidx.append(i)
        if line.find('PHASECAL') > 0:
            # Subtract 1 min from time in this line and use that time in previous (ACQUIRE) line
            mjd = Time(line[:19]).mjd - 60./86400
            keeplines[i-1] = Time(mjd,format='mjd').iso[:19] + keeplines[i-1][19:]
    for i in sunidx[1:-1]:
        if keeplines[i+1].find('ACQUIRE') > 1:
            rmidx.append(i+1)
            rmidx.append(i+2)
    outlines = []
    for i,line in enumerate(keeplines):
        if not i in rmidx:
            outlines.append(line)
    outlines[-1] = outlines[-1][:19]+' STOW'
    return outlines
    
