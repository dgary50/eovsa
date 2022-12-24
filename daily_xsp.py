# Name:
#   DAILY_XSP
# Create daily "all_day" dynamic spectrum plots of EOVSA solar 
# cross-correlation data
#
#  2017-08-29  DG
#    First written, based on routine in read_idb.py
#  2019-02-26  DG
#    Added a second URL for GOES data, in case the first is unreachable.
#    Also added a FITS argument to command line, for optional output of FITS file.
#    FITS-ONLY means make a FITS file and stop.  FITS by itself means make a FITS file
#    in addition to a plot.
#  2020-05-10  DG
#    Updated cal_qual() to use util.get_idbdir() to find IDB root path.
#  2020-05-11  DG
#    Further update to make this work on the DPP.
#  2020-05-17  DG
#    The GOES data continues to be problematic because the Goddard respositories
#    are not kept up to date.  I added a call to by new goes.py routine get_goes(),
#    which grabs the latest 7 days of GOES data from NOAA.  If that fails or the
#    requested date is more than 7 days from the current date, then it falls back
#    to the original code and tries to get the data from Goddard.
#  2020-05-23  DG
#    A classic blunder, using median() instead of nanmedian() on GOES data with nans!
#    Now fixed.
#  2021-07-11  DG
#    Changes to cal_qual() to default to creating web plots, to display both TP and XP
#    calibration, and to fix the limitation that plots only covered a max of 600 s.
#    Also added some line plots (freq indexes 100,300) to make it more quantitative.
#    An important change is to make cal_qual() part of the daily plot generation.
#  2021-07-29  DG
#    Clip line plots for frequency indexes 100 and 300 so as not to exceed the number
#    number of frequencies (451).
#  2022-02-10  DG
#    In cal_qual(), writing to the /common/webplots/flaremon/daily/ folder failed
#    if run by someone other than user.  In that case, the plot is now created in /tmp/.
#  2022-03-05  DG
#    Small change to call get_projects() instead of flare_monitor(), along with the
#    added nosql=True argument so that this works without the SQL server.  Only downside
#    is that ACQUIRE states are not displayed (they are in SQL but not fdb files).
#  2022-03-10  DG
#    A couple of other changes due to loss of SQL, marked with comments.
#  2022-06-23  DG
#    Suddenly getting GOES data returns unequal-length arrays for time and data,
#    so truncate output from get_goes() to equal lengths.
#

if __name__ == "__main__":
    import matplotlib
    matplotlib.use("Agg")
    import os
    try:
        user_paths = os.environ['PYTHONPATH'].split(os.pathsep)
    except:
        user_paths = []
    print user_paths


import read_idb as ri
from util import Time
import numpy as np
import glob
from goes import get_goes

def get_goes_data(t=None,sat_num=None):
    ''' Reads GOES data from https://umbra.nascom.nasa.gov/ repository, for date
        and satellite number provided.  If sat_num is None, data for all available 
        satellites are downloaded, with some sanity check used to decide the best.
        If the Time() object t is None, data for the day before the current date 
        are read (since there is a delay of 1 day in availability of the data).
        
        Returns:
           goes_t    GOES time array in plot_date format
           goes_data GOES 1-8 A lightcurve
        '''
    # Can short-circuit the entire code below this block by using my goes.get_goes() routine
    lo, hi, goes_t = get_goes()
    if len(goes_t) != 0:
        # Got the data, now isolate the requested day
        good, = np.where(np.floor(goes_t.mjd) == np.floor(t.mjd))
        if len(good) != 0:
            lent = len(goes_t.plot_date) 
            lend = len(lo)
            if lent != lend:
                lmin = min(lent,lend)
                goes_t.plot_date = goes_t.plot_date[:lmin]
                lo = lo[:lmin]
            return goes_t.plot_date,lo
                
    from sunpy.util.config import get_and_create_download_dir
    import shutil
    from astropy.io import fits
    import urllib2
    if t is None:
        t = Time(Time.now().mjd - 1,format='mjd')
    yr = t.iso[:4]
    datstr = t.iso[:10].replace('-','')
    try:
        if sat_num is None:
            try:
                f = urllib2.urlopen('https://umbra.nascom.nasa.gov/goes/fits/'+yr, timeout=3)
            except:
                f = urllib2.urlopen('https://hesperia.gsfc.nasa.gov/goes/'+yr,timeout=3)
            lines = f.readlines()
            sat_num = []
            for line in lines:
                idx = line.find(datstr)
                if idx != -1:
                    sat_num.append(line[idx-2:idx])
        if type(sat_num) is int:
            sat_num = [str(sat_num)]
        filenames = []
        for sat in sat_num:
            filename = 'go'+sat+datstr+'.fits'
            try:
                url = 'https://umbra.nascom.nasa.gov/goes/fits/'+yr+'/'+filename
                f = urllib2.urlopen(url, timeout=3)
            except:
                url = 'https://hesperia.gsfc.nasa.gov/goes/'+yr+'/'+filename
                f = urllib2.urlopen(url, timeout=3)
            with open(get_and_create_download_dir()+'/'+filename,'wb') as g:
                shutil.copyfileobj(f,g)
            filenames.append(get_and_create_download_dir()+'/'+filename)
        pmerit = 0
        for file in filenames:
            gfits = fits.open(file)
            data = gfits[2].data['FLUX'][0][:,0]
            good, = np.where(data > 1.e-8)
            tsecs = gfits[2].data['TIME'][0]
            merit = len(good)
            date_elements = gfits[0].header['DATE-OBS'].split('/')
            if merit > pmerit:
                print 'File:',file,'is best'
                pmerit = merit
                goes_data = data
                goes_t = Time(date_elements[2]+'-'+date_elements[1]+'-'+date_elements[0]).plot_date + tsecs/86400.
        try:
            return goes_t, goes_data
        except:
            print 'No good GOES data for',datstr
            return None, None
    except:
        print 'GOES site unreachable?'
        return None, None
        
def allday_udb(t=None, doplot=True, goes_plot=True, savfig=False, savfits=False, gain_corr=True):
    if savfits:
        import xspfits #jmm, 2018-01-05
    # Plots (and returns) UDB data for an entire day
    from flare_monitor import get_projects
    if t is None:
        t = Time.now()
    # Cannot get a GOES plot unless doplot is True
    if goes_plot: doplot = True
    date = t.iso[:10].replace('-','')
    # Look also at the following day, up to 9 UT
    date2 = Time(t.mjd + 1,format='mjd').iso[:10].replace('-','')
    year = date[:4]
    files = glob.glob('/data1/eovsa/fits/UDB/'+year+'/UDB'+date+'*')
    files.sort()
    files2 = glob.glob('/data1/eovsa/fits/UDB/'+year+'/UDB'+date2+'0*')
    files2.sort()
    files = np.concatenate((np.array(files),np.array(files2)))
    # Eliminate files starting before 10 UT on date (but not on date2)
    for i,file in enumerate(files):
        if file[-6] != '0':
            break
    try:
        files = files[i:]
    except:
        print 'No files found in /data1/eovsa/fits/UDB/ for',date
        return {}
    #import pdb; pdb.set_trace()
    out = ri.read_idb(files,src='Sun')
    if out.keys() == []:
        print 'Read error, or no Sun data in',files
        return {}
    if gain_corr:
        import gaincal2 as gc
        out = gc.apply_gain_corr(out)
    trange = Time(out['time'][[0,-1]], format = 'jd')
    fghz = out['fghz']
    pdata = np.nansum(np.nansum(np.abs(out['x'][0:11,:]),1),0)  # Spectrogram to plot
    if savfits:
        print "***************** PDATA OUTPUT *********"
        print pdata.shape
        xspfits.daily_xsp_writefits(out, pdata)
    if doplot:
        import matplotlib.pylab as plt
        from matplotlib.dates import DateFormatter
        f, ax = plt.subplots(1,1,figsize=(14,5))
        X = np.sort(pdata.flatten())   # Sorted, flattened array
        # Set any time gaps to nan
        tdif = out['time'][1:] - out['time'][:-1]
        bad, = np.where(tdif > 120./86400)  # Time gaps > 2 minutes
        pdata[:,bad] = 0
        vmax = X[int(len(X)*0.85)]  # Clip at 15% of points
        im = ax.pcolormesh(Time(out['time'],format='jd').plot_date,out['fghz'],pdata,vmax=vmax)
        plt.colorbar(im,ax=ax,label='Amplitude [arb. units]')
        ax.xaxis_date()
        ax.xaxis.set_major_formatter(DateFormatter("%H:%M"))
        ax.set_ylim(fghz[0], fghz[-1])
        ax.set_xlabel('Time [UT]')
        ax.set_ylabel('Frequency [GHz]')
        ax.set_title('EOVSA 1-min Data for '+t.iso[:10])
        f.autofmt_xdate(bottom=0.15)

        if goes_plot:
            #from sunpy import lightcurve
            #from sunpy.time import TimeRange
            # Initially assign GOES times as None
            goes_t = None
            goes_t2 = None
            # Get GOES data for overplotting
            #goes_tr = TimeRange(trange.iso)
            goes_label = [' A',' B',' C',' M',' X']
            # The GOES label is placed to start 20 min into the day
            goes_label_time = Time(out['time'][[0]], format = 'jd').plot_date + 0.014
            rightaxis_label_time = trange[1].plot_date

            # Retrieve GOES data for the day, but this only goes to end of UT day
            goes_t, goes_data = get_goes_data(trange[0])
            if int(trange[1].mjd) != int(trange[0].mjd):
                goes_t2, goes_data2 = get_goes_data(trange[1])
            if goes_t is None and goes_t2 is None:
                ax.text (goes_label_time, 12, 'GOES soft x-ray data missing', color = 'yellow')
            else:
                if not goes_t is None:
                    goes_data = 2* (np.log10(goes_data + 1.e-9)) + 26
                    ax.plot_date(goes_t, goes_data,'-',color='yellow')
                    ytext = np.nanmedian(goes_data) - 1
                else:
                    ytext = None
                if not goes_t2 is None:
                    goes_data2 = 2* (np.log10(goes_data2 + 1.e-9)) + 26
                    ax.plot_date(goes_t2, goes_data2,'-',color='yellow')
                    ytext2 = np.nanmedian(goes_data2) - 1
                    if ytext:
                        ytext = (ytext+ytext2)/2
                    else:
                        ytext = ytext2
                ax.text (goes_label_time, ytext, 'GOES soft x-ray data', color = 'yellow')
            # try:
                # goes = lightcurve.GOESLightCurve.create(goes_tr)
                # if len(np.where(goes.data['xrsb'] != 0.0)[0]) < 100:
                    # # Looks like the GOES data are all zero, so just skip it
                    # ax.text (goes_label_time, 12, 'GOES soft x-ray data missing', color = 'yellow')
                # else:
                    # goes.data['xrsb'] = 2* (np.log10(goes.data['xrsb'] + 1.e-9)) + 26
                    # ytext = np.median(goes.data['xrsb']) - 1
                    # ax.text (goes_label_time, ytext, 'GOES soft x-ray data', color = 'yellow')
                    # goes.data['xrsb'].plot(color = 'yellow')
            # except:
                # # Looks like the GOES data do not exist, so just skip it
                # ax.text (goes_label_time, 12, 'GOES soft x-ray data missing', color = 'yellow')
            for k,i in enumerate([10,12,14,16,18]):
                ax.text(rightaxis_label_time, i-0.4, goes_label[k], fontsize = '12')
                ax.plot_date(rightaxis_label_time + np.array([-0.005,0.0]),[i,i],'-',color='yellow')
            # try:
                # # If the day goes past 0 UT, get GOES data for the next UT day
                # if int(trange[1].mjd) != int(trange[0].mjd):
                    # goes_tr2 = TimeRange([trange[1].iso[:10], trange[1].iso])
                    # goesday2 = lightcurve.GOESLightCurve.create(goes_tr2)
                    # if len(np.where(goesday2.data['xrsb'] != 0.0)[0]) < 100:
                        # pass
                    # else:
                        # goesday2.data['xrsb'] = 2* (np.log10(goesday2.data['xrsb'] + 1.e-9)) + 26
                        # goesday2.data['xrsb'].plot(color = 'yellow')
            # except:
                # # Looks like the GOES data do not exist, so just skip it
                # pass
        pstart = Time(t.iso[:10]+' 13:30').plot_date
        prange = [pstart,pstart+13./24]
        ax.set_xlim(prange)

        projdict = get_projects(t, nosql=True)   # Hopefully temporary call that is independent of SQL server
        if projdict == {}:
            print 'No annotation can be added to plot for',t.iso[:10]
        else:
            defcolor = '#ff7f0e'
            nscans = len(projdict['Timestamp'])
            SOS = Time(projdict['Timestamp'],format='lv').plot_date
            EOS = Time(projdict['EOS'],format='lv').plot_date
            yran = np.array(ax.get_ylim())
            for i in range(nscans):
                uti = SOS[i]*np.array([1.,1.])
                #if uti[0] >= trange[0].plot_date:
                ax.plot_date(uti,yran,'g',lw=0.5)
                if projdict['Project'][i] == 'NormalObserving' or projdict['Project'][i] == 'Normal Observing':
                    ax.text(uti[0],yran[1]*0.935,'SUN',fontsize=8, color = defcolor, clip_on=True)
                elif projdict['Project'][i] == 'None':
                    ax.text(uti[0],yran[1]*0.975,'IDLE',fontsize=8, color = defcolor, clip_on=True)
                elif projdict['Project'][i][:4] == 'GAIN':
                    ax.text(uti[0],yran[1]*0.955,'GCAL',fontsize=8, color = defcolor, clip_on=True)
                elif projdict['Project'][i] == 'SOLPNTCAL':
                    ax.text(uti[0],yran[1]*0.955,'TPCAL',fontsize=8, color = defcolor, clip_on=True)
                elif projdict['Project'][i] == 'PHASECAL':
                    ax.text(uti[0],yran[1]*0.955,'PCAL',fontsize=8, color = defcolor, clip_on=True)
                else:
                    ax.text(uti[0],yran[1]*0.975,projdict['Project'][i],fontsize=8, color = defcolor, clip_on=True)
            known = ['GAIN','PHAS','SOLP']  # known calibration types (first 4 letters)
            for i in range(nscans):
                uti = EOS[i]*np.array([1.,1.])
                ax.plot_date(uti,yran,'r--',lw=0.5)
                uti = np.array([SOS[i],EOS[i]])
                if projdict['Project'][i] == 'NormalObserving':
                    ax.plot_date(uti,yran[1]*np.array([0.93,0.93]),ls='-',marker='None',color='#aaffaa',lw=2,solid_capstyle='butt')
                elif projdict['Project'][i][:4] in known:
                    ax.plot_date(uti,yran[1]*np.array([0.95,0.95]),ls='-',marker='None',color='#aaaaff',lw=2,solid_capstyle='butt')
                else:
                    ax.plot_date(uti,yran[1]*np.array([0.97,0.97]),ls='-',marker='None',color='#ffaaaa',lw=2,solid_capstyle='butt')

            if savfig:
                plt.savefig('/common/webplots/flaremon/daily/'+date[:4]+'/XSP'+date+'.png',bbox_inches='tight')
    return out

def cal_qual(t=None, savfig=True):
    ''' Check the quality of the total power and gain calibrations for a given date
    '''
    import cal_header as ch
    from stateframe import extract
    import dump_tsys as dt
    import pipeline_cal as pc
    import matplotlib.pylab as plt
    import rstn
    from util import get_idbdir
    import socket
    
    if t is None:
        t = Time.now()
    mjd = t.mjd
    # First check whether the total power calibration is current
    caltype = 10
    xml, buf = ch.read_cal(caltype,t=t)
    tp_mjd = Time(extract(buf,xml['SQL_timestamp']),format='lv').mjd
    if mjd - tp_mjd > 0.5:
        print 'CAL_QUAL: Warning, TP Calibration not (yet) available for this date.'
    # Find GCAL scan for this date
    fdb = dt.rd_fdb(Time(mjd,format='mjd'))
    gcidx, = np.where(fdb['PROJECTID'] == 'GAINCALTEST')
    if len(gcidx) == 1:
        datadir = get_idbdir(t) + fdb['FILE'][gcidx][0][3:11]+'/'
        # List of GCAL files
        gcalfile = [datadir+i for i in fdb['FILE'][gcidx]]
    else:
        print 'CAL_QUAL: Warning, no GAINCALTEST scan for this date.  Will try using the GAINCALTEST from previous day.'
        fdb = dt.rd_fdb(Time(mjd-1,format='mjd'))
        gcidx, = np.where(fdb['PROJECTID'] == 'GAINCALTEST')
        if len(gcidx) == 1:
            datadir = get_idbdir(t)
            # Add date path if on pipeline
            # if datadir.find('eovsa') != -1: datadir += fdb['FILE'][gcidx][0][3:11]+'/'
            host = socket.gethostname()
            if host == 'pipeline': datadir += fdb['FILE'][gcidx][0][3:11]+'/'
            # List of GCAL files
            gcalfile = [datadir+i for i in fdb['FILE'][gcidx]]
        else:
            print 'CAL_QUAL: Error, no GAINCALTEST scan for previous day.'
            return
    # Find SOLPNTCAL scan for this date
    fdb = dt.rd_fdb(Time(mjd,format='mjd'))
    gcidx, = np.where(fdb['PROJECTID'] == 'SOLPNTCAL')
    if len(gcidx) > 0:
        datadir = get_idbdir(t)
        # Add date path if on pipeline
        # if datadir.find('eovsa') != -1: datadir += fdb['FILE'][gcidx][0][3:11]+'/'
        host = socket.gethostname()
        if host == 'pipeline': datadir += fdb['FILE'][gcidx][0][3:11]+'/'
        # List of SOLPNTCAL files
        solpntfile = [datadir+i for i in fdb['FILE'][gcidx]]
    else:
        print 'CAL_QUAL: Error, no SOLPNTCAL scan(s) for this date.'
        return
    files = gcalfile+solpntfile
    outnames = []
    for file in files:
        outnames.append(pc.udb_corr(file, calibrate=True, attncal=True, desat=True))
    out = ri.read_idb(outnames, srcchk=False)
    nt = len(out['time'])
    nf = len(out['fghz'])
    tpfac = 500./nf

    frq, flux = rstn.rd_rstnflux(t)
    s = rstn.rstn2ant(frq, flux, out['fghz']*1000., t)
    fluximg = s.repeat(nt).reshape(nf,nt)
    f, ax = plt.subplots(4,7)
    f.set_size_inches(16,7,forward=True)
    f.tight_layout(rect=[0.0,0.0,1,0.95])
    ax.shape = (2, 14)
    for i in range(13):
        for j in range(2):
            ax[j,i].imshow(out['p'][i,j],aspect='auto',origin='lower',vmax=np.max(s),vmin=0)
            ax[j,i].plot(np.clip(out['p'][i,j,int(nf/3.)]/tpfac,0,nf),linewidth=1)
            ax[j,i].plot(np.clip(out['p'][i,j,int(2*nf/3.)]/tpfac,0,nf),linewidth=1)
            ax[j,i].set_title('Ant '+str(i+1)+[' X Pol',' Y Pol'][j],fontsize=10)
    for j in range(2): 
        ax[j,13].imshow(fluximg,aspect='auto',origin='lower',vmax=np.max(s),vmin=0)
        ax[j,13].set_title('RSTN Flux',fontsize=10)
    for i in range(13):
        for j in range(2):
            ax[j,i].plot(np.clip(fluximg[int(nf/3.)]/tpfac,0,nf),'--',linewidth=1,color='C0')
            ax[j,i].plot(np.clip(fluximg[int(2*nf/3.)]/tpfac,0, nf),'--',linewidth=1,color='C1')

    f.suptitle('Total Power Calibration Quality for '+t.iso[:10])
    date = t.iso[:10].replace('-','')
    if savfig:
        try:
            plt.savefig('/common/webplots/flaremon/daily/'+date[:4]+'/QUAL_'+date+'TP.png')
        except:
            plt.savefig('/tmp/'+date[:4]+'/QUAL_'+date+'TP.png')
            print 'The .png file could not be created in the /common/webplots/flaremon/daily/ folder.'
            print 'A copy was created in /tmp/.'

    f, ax = plt.subplots(4,7)
    f.set_size_inches(16,7,forward=True)
    f.tight_layout(rect=[0.0,0.0,1,0.95])
    ax.shape = (2, 14)
    for i in range(13):
        for j in range(2):
            ax[j,i].imshow(np.real(out['a'][i,j]),aspect='auto',origin='lower',vmax=np.max(s),vmin=0)
            ax[j,i].plot(np.clip(np.real(out['a'][i,j,int(nf/3.)]/tpfac),0,nf),linewidth=1)
            ax[j,i].plot(np.clip(np.real(out['a'][i,j,int(2*nf/3.)]/tpfac),0,nf),linewidth=1)
            ax[j,i].set_title('Ant '+str(i+1)+[' X Pol',' Y Pol'][j],fontsize=10)
    for j in range(2): 
        ax[j,13].imshow(fluximg,aspect='auto',origin='lower',vmax=np.max(s),vmin=0)
        ax[j,13].set_title('RSTN Flux',fontsize=10)
    for i in range(13):
        for j in range(2):
            ax[j,i].plot(np.clip(fluximg[int(nf/3.)]/tpfac,0,nf),'--',linewidth=1,color='C0')
            ax[j,i].plot(np.clip(fluximg[int(2*nf/3.)]/tpfac,0,nf),'--',linewidth=1,color='C1')
    f.suptitle('Cross-Power Calibration Quality for '+t.iso[:10])
    date = t.iso[:10].replace('-','')
    if savfig:
        try:
            plt.savefig('/common/webplots/flaremon/daily/'+date[:4]+'/QUAL_'+date+'XP.png')
        except:
            plt.savefig('/tmp/'+date[:4]+'/QUAL_'+date+'XP.png')
            print 'The .png file could not be created in the /common/webplots/flaremon/daily/ folder.'
            print 'A copy was created in /tmp/.'
    

if __name__ == "__main__":
    ''' For non-interactive use, use a backend that does not require a display
        Usage python /common/python/current/daily_xsp.py <date>
        
        If optional argument date is given (as YYYY-MM-DD), data are processed for 
           that date only.
        If the date is omitted, data are processed for the previous two UT days
          (yesterday and day-before-yesterday)
    '''
    import glob, sys
    t = None
    t2 = None
    # Default parameters
    savfits = False
    savfig = True
    goes_plot = True
    doplot=True
    # ************ This line added due to loss of SQL **************
    #gain_corr = False
    gain_corr = True
    argin = ''
    if len(sys.argv) >= 2:
        try:
            t = Time(sys.argv[1])
            print t.iso
            if len(sys.argv) == 3:
                argin = sys.argv[2].upper()
        except:
            argin = sys.argv[1].upper()  # Any following arguments are ignored
    if argin == 'FITS-ONLY':
        # Asking for FITS-ONLY means skip all other features.
        savfits = True
        savfig = False
        goes_plot = False
        doplot=False
    elif argin == 'FITS':
        # Asking for FITS means also do all other features.
        savfits = True
    elif argin != '':
        print 'Cannot interpret',argin,'as valid time, or string FITS or FITS-ONLY.'
        exit()
    if t is None:
        t = Time.now()   # Get today's date
        t2 = Time(t.mjd-2,format='mjd')   # Set t2 to day-before-yesterday
        t = Time(t.mjd-1,format='mjd')    # Set t to yesterday
    print t.iso[:19],': ',
    blah = allday_udb(t=t, doplot=doplot, goes_plot=goes_plot, savfig=savfig, savfits=savfits, gain_corr=gain_corr)   # Process time t
    if goes_plot and not t2 is None:
        # Do this second date only if goes_plot is True
        blah = allday_udb(t=t2, savfig=True, gain_corr=gain_corr)   # Process time t2
    # ************ This line commented out due to loss of SQL **************
    # cal_qual(Time(t.iso[:10]))  # Make daily plot of calibration quality