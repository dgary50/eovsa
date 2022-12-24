'''
   Module for plotting the median of the front-end RF detector voltages
   from the stateframe SQL database, as a crude flare monitor'''
#
# History:
#   2014-Dec-20  DG
#      First written.
#   2014-Dec-21  DG
#      Added annotation and information about source.
#   2014-Dec-22  DG
#      Cleaned up error handling so that procedure will work as cron job.
#   2014-Dec-24  DG
#      Added printing of date, for cron log file.
#   2014-Dec-26  DG
#      Fix bug when there are no gaps in scan_off_times.  Also set xlim to 15-24 h
#   2015-Feb-20  DG
#      Add v38/v39 stateframe table update.
#   2015-Mar-10  DG
#      Add v39/v42 stateframe table update.
#   2015-Mar-29  DG
#      Add v42/v45 stateframe table update.
#   2015-Mar-31  DG
#      Add v45/v46 stateframe table update.
#   2015-Apr-02  DG
#      Finally made this version independent!  Added calls to new routine
#      dbutil.find_table_version().
#   2015-May-29  DG
#      Converted from using datime() to using Time() based on astropy.
#   2015-Jul-04  DG
#      Added xdata_display() routine and plotting of cross-correlation
#      amplitude in files named /common/webplots/flaremon/XSP*.png
#   2015-Aug-30  DG
#      Added flaremeter() routine to calculate medians of cross-correlation
#      amplitudes across all baselines, polarizations, and frequencies. Added
#      code to plot this median information on flare_monitor plot.  Also,
#      extend timerange of plots to 13 UT on current day to 02 UT on next day.
#      A new series of binary files are created containing the flaremeter
#      information for each day, currently in /common/webplots/flaremon/flaremeter/.
#   2015-Sep-06  DG
#      Added code in xdata_display and __main__ to read fdb files into next 
#      UT day, so that scans that extend past 24 UT are fully plotted.
#   2015-Sep-07  DG
#      Attempt to fix bug in extending date past 24 UT
#   2016-Jun-30  DG
#      Update to work with 16-ant-correlator data and new routine read_idb()
#      Also does correct scaling of x-corr level in case of extraneous inf
#   2016-Jul-15  DG
#      Add sk_flag to xdata display.
#   2016-Aug-04  DG
#      After update of numpy, my medians no longer worked.  Changed to nanmedian.
#   2017-Mar-20  DG
#      Changes to get this working for 300 MHz correlator.  Skip SK flagging, skip
#      call to flaremeter() [just takes too long!], and skip get_history() call.
#   2017-Apr-06  DG
#      Changes to put identifying text on plot for more types of calibration.
#   2017-Aug-11  DG
#      Changes to rationalize the handling of the scan start-stop information in
#      flare_monitor(), and tweaks to the line plot.
#   2017-Sep-08  DG
#      Fixed small bug that caused crash on error return from xdata_display()
#   2021-Nov-13  DG
#      Added RT_flare_monitor(), which reads a file RT_latest.txt and creates
#      a summary file RT_<datstr>.txt in the same /data1/RT folder on the dpp.
#   2022-Feb-09  DG
#      Changed the flare_monitor() code to keep data if within 0.1 units of
#      previous measurement instead of stringent 0.01 units.
#   2022-Mar-05  DG
#      Due to a (hopefully) temporary outage of the SQL server, I wrote a
#      new version of get_projects() called get_projects_nosql() that
#      does almost the same thing using the FDB files.  Calling get_projects()
#      with nosql=True also works.  The __main__ routine has some changes
#      to avoid reading from SQL.
#   2022-Jul-26  DG
#      Fix error in projdict['EOS'] returned from get_projects_nosql().
#   2022-Sep-05  DG
#      For some reason, rd_fdb is returning an extra blank file now.
#      I added a check for that in get_projects_nosql().
#
import numpy as np
from util import Time, get_idbdir

def RT_flare_monitor():
    ''' Read "real-time" data file once and obtain the median over antenna
        and frequency, appending the result to /data1/RT/RT_<date>.txt.
        
        Returns ut times in plot_date format and median values.
    '''
    tp = np.zeros((100,16))
    amp = np.zeros((100,16))
    f = open('/data1/RT/RT_latest.txt','r')
    lines = f.readlines()
    f.close()
    mjd = Time.now().mjd - 0.5
    datstr = Time(mjd,format='mjd').iso[:10]
    tstr = lines[1].split(':')[1][1:7]
    tstr = tstr[:2]+':'+tstr[2:4]+':'+tstr[4:]
    print datstr+' '+tstr+' ('+Time.now().iso[:19]+')'
    t = Time(datstr+' '+tstr)

    try:
        # Read last line of output and compare with input. Skip if the same.
        f = open('/data1/RT/RT_'+datstr+'.txt','r')
        lastline = f.readlines()[-1]
        f.close()
        lastt = Time(lastline[:19])
        if lastt == t:
            # Times are the same, so do not write to output file
            return
    except:
        pass

    # Read entire file and form median values of TP and AMP        
    for i,line in enumerate(lines[10::10]):
        tp[i] = np.array(line.split()).astype(float)
    for i,line in enumerate(lines[11::10]):
        amp[i] = np.array(line.split()).astype(float)
    tpmed = np.median(tp[:,:13])
    ampmed = np.median(amp[:,:13])

    # Open output file for appending, write median values, and close.
    f = open('/data1/RT/RT_'+datstr+'.txt','a')
    f.write(t.iso+' {:8.4f} {:8.4f}'.format(tpmed,ampmed)+'\n')
    f.close()

def get_projects(t, nosql=False):
    ''' Read all projects from SQL for the current date and return a summary
        as a dictionary with keys Timestamp, Project, and EOS (another timestamp)
    '''
    if nosql == True:
        return get_projects_nosql(t)
    import dbutil
    # timerange is 12 UT to 12 UT on next day, relative to the day in Time() object t
    trange = Time([int(t.mjd) + 12./24,int(t.mjd) + 36./24],format='mjd')
    tstart, tend = trange.lv.astype('str')
    cursor = dbutil.get_cursor()
    mjd = t.mjd
    # Get the project IDs for scans during the period
    verstrh = dbutil.find_table_version(cursor,trange[0].lv,True)
    if verstrh is None:
        print 'No scan_header table found for given time.'
        return {}
    query = 'select Timestamp,Project from hV'+verstrh+'_vD1 where Timestamp between '+tstart+' and '+tend+' order by Timestamp'
    projdict, msg = dbutil.do_query(cursor, query)
    if msg != 'Success':
        print msg
        return {}
    elif len(projdict) == 0:
        # No Project ID found, so return data and empty projdict dictionary
        print 'SQL Query was valid, but no Project data were found.'
        return {}
    projdict['Timestamp'] = projdict['Timestamp'].astype('float')  # Convert timestamps from string to float
    for i in range(len(projdict['Project'])): projdict['Project'][i] = projdict['Project'][i].replace('\x00','')
    projdict.update({'EOS':projdict['Timestamp'][1:]})
    projdict.update({'Timestamp':projdict['Timestamp'][:-1]})
    projdict.update({'Project':projdict['Project'][:-1]})
    cursor.close()
    return projdict

def get_projects_nosql(t):
    ''' Read all projects from FDB file for the current date and return a summary
        as a dictionary with keys Timestamp, Project, and EOS (another timestamp)
    '''
    import dump_tsys as dt
    # timerange is 12 UT to 12 UT on next day, relative to the day in Time() object t
    trange = Time([int(t.mjd) + 12./24,int(t.mjd) + 36./24],format='mjd')
    tstart = t.iso[2:10].replace('-','')+'120000'
    t2 = Time(int(t.mjd)+1, format='mjd')
    tend = t2.iso[2:10].replace('-','')+'120000'
    fdb = dt.rd_fdb(t)
    fdb2 = dt.rd_fdb(t2)
    if fdb == {}:
        # No FDB file found, so return empty project dictionary
        print 'No Project data [FDB file] found for the given date.'
        return {}
    if fdb2 == {}:
        pass
    else:
        #  Concatenate the two dicts into one
        if fdb['FILE'][-1] == '':
            fdb = dict([(k, np.concatenate((fdb.get(k,[])[:-1],fdb2.get(k,[])))) for k in set(fdb)|set(fdb2)])
        else:
            fdb = dict([(k, np.concatenate((fdb.get(k,[]),fdb2.get(k,[])))) for k in set(fdb)|set(fdb2)])  
    # Get "good" indexes for times between 12 UT on date and 12 UT on next date
    gidx, = np.where(np.logical_and(fdb['SCANID']>tstart,fdb['SCANID']<tend))        
    scanid,idx = np.unique(fdb['SCANID'][gidx],return_index=True)
    sidx = gidx[idx]   # Indexes into fdb for the start of each scan
    eidx = np.concatenate((sidx[1:]-1,np.array([gidx[-1]])))   # Indexes into fdb for the end of each scan
    # Get the project IDs for scans during the period
    projdict = {'Timestamp':fdb['ST_TS'][sidx].astype(float),
                'Project':fdb['PROJECTID'][sidx],
                'EOS':fdb['EN_TS'][eidx].astype(float)}
    return projdict

def flare_monitor(t):
    ''' Get all front-end power-detector voltages for the given day
        from the stateframe SQL database, and obtain the median of them, 
        to use as a flare monitor.
        
        Returns ut times in plot_date format and median voltages.
    '''
    import dbutil
    # timerange is 12 UT to 12 UT on next day, relative to the day in Time() object t
    trange = Time([int(t.mjd) + 12./24,int(t.mjd) + 36./24],format='mjd')
    tstart, tend = trange.lv.astype('str')
    cursor = dbutil.get_cursor()
    mjd = t.mjd
    verstr = dbutil.find_table_version(cursor,tstart)
    if verstr is None:
        print 'No stateframe table found for given time.'
        return tstart, [], {}
    query = 'select Timestamp,Ante_Fron_FEM_HPol_Voltage,Ante_Fron_FEM_VPol_Voltage from fV'+verstr+'_vD15 where timestamp between '+tstart+' and '+tend+' order by timestamp'
    data, msg = dbutil.do_query(cursor, query)
    if msg != 'Success':
        print msg
        return tstart, [], {}
    for k,v in data.items():
        data[k].shape = (len(data[k])/15,15)
    hv = []
    try:
        ut = Time(data['Timestamp'][:,0].astype('float'),format='lv').plot_date 
    except:
        print 'Error for time',t.iso
        print 'Query:',query,' returned msg:',msg
        print 'Keys:', data.keys()
        print data['Timestamp'][0,0]
    hfac = np.median(data['Ante_Fron_FEM_HPol_Voltage'].astype('float'),0)
    vfac = np.median(data['Ante_Fron_FEM_VPol_Voltage'].astype('float'),0)
    for i in range(4):
        if hfac[i] > 0:
            hv.append(data['Ante_Fron_FEM_HPol_Voltage'][:,i]/hfac[i])
        if vfac[i] > 0:
            hv.append(data['Ante_Fron_FEM_VPol_Voltage'][:,i]/vfac[i])
    #import pdb; pdb.set_trace()
    flm = np.median(np.array(hv),0)
    good = np.where(abs(flm[1:]-flm[:-1])<0.1)[0]

    projdict = get_projects(t)
    return ut[good],flm[good],projdict

def xdata_display(t,ax=None):
    ''' Given the time as a Time object, search the FDB file for files
        associated with the scan for that time and create a dynamic spectrogram
        on the axis specified by ax, or on a new plot if no ax. If the requested
        time is more than 10 minutes after the last file of that scan, returns
        None to indicate no plot.
        
        Skip SK flagging [2017-Mar-20 DG]
    '''
    import time, os
    import dump_tsys
    #import get_X_data2 as gd
    import read_idb as ri
    import spectrogram_fit as sp

    fdb = dump_tsys.rd_fdb(t)
    # Get files from next day, in case scan extends past current day
    t1 = Time(t.mjd + 1,format='mjd')
    fdb1 = dump_tsys.rd_fdb(t1)
    # Concatenate the two days (if the second day exists)
    if fdb1 != {}:
        for key in fdb.keys():
            fdb[key] = np.concatenate((fdb[key],fdb1[key]))
            
    # Find unique scan IDs
    scans, idx = np.unique(fdb['SCANID'],return_index=True)

    # Limit to scans in 'NormalObserving' mode
    good, = np.where(fdb['PROJECTID'][idx] == 'NormalObserving')
    if len(good) > 0:
        scans = scans[good]
    else:
        print 'No NormalObserving scans found.'
        return None, None, None, None

    # Find scanID that starts earlier than, but closest to, the current time
    for i,scan in enumerate(scans):
        print scan
        dt = t - Time(time.strftime('%Y-%m-%d %H:%M:%S',time.strptime(scan,'%y%m%d%H%M%S')))
        if dt.sec > 0.:
            iout = i
    scan = scans[iout]

    # Find files for this scan
    fidx, = np.where(fdb['SCANID'] == scan)
    tlevel = None
    bflag = None
    if len(fidx) > 0:
        files = fdb['FILE'][fidx]
        # Find out how old last file of this scan is, and proceed only if less than 20 minutes
        # earlier than the time given in t.
        try:
            dt = t - Time(time.strftime('%Y-%m-%d %H:%M:%S',time.strptime(files[-1],'IDB%Y%m%d%H%M%S')))
        except:
            dt = 10000.  # Forces skip of plot creation
            print 'Unexpected FDB file format.'
            scan = None
        if dt.sec < 1200.:
            # This is a currently active scan, so create the figure
            path = '/data1/IDB/'
            if not os.path.isdir(path+files[0]):
                # Look in /dppdata1
                datstr = t.iso[:10].replace('-','')
                # path = '/data1/eovsa/fits/IDB/'+datstr+'/'
                path = get_idbdir(t)+datstr+'/'
                if not os.path.isdir(path+files[0]):
                    print 'No files found for this scan ID',scan
                    scan = None
                    times = None
                    return scan, tlevel, bflag, times
            filelist = files
            files = []
            for i,file in enumerate(filelist):
                files.append(path+file)
            # data, uvw, fghz, times = gd.get_X_data(files)
            out = ri.read_idb(files)
            #out = ri.flag_sk(out)  # Skip flagging for sk
            fghz = out['fghz']
            times = Time(out['time'],format='jd')
            data = out['x']
            if ax is not None:
                datstr = times[0].iso[:10]
                ax.set_xlabel('Time [UT on '+datstr+']')
                ax.set_ylabel('Frequency [GHz]')
                ax.set_title('EOVSA Summed Cross-Correlation Amplitude for '+datstr)
            pdata = np.nansum(np.nansum(np.abs(data[0:11,:]),1),0)  # Spectrogram to plot
            X = np.sort(pdata.flatten())   # Sorted, flattened array
            dmax = X[int(len(X)*0.95)]  # Clip at 5% of points
            sp.plot_spectrogram(fghz, times, pdata, 
                                ax=ax, logsample=None, xdata=True, cbar=True, dmax=dmax)
            #tlevel, bflag = flaremeter(data)
        else:
            print 'Time',dt.sec,'is > 1200 s after last file of last NormalObserving scan.  No plot created.'
            scan = None
            times = None
    else:
        print 'No files found for this scan ID',scan
        scan = None
    return scan, tlevel, bflag, times

def flaremeter(data):
    ''' Obtain median of data across baselines, polarizations, and frequencies to create a
        time series indicated whether a flare has occurred.  Values returned will be close
        to unity if no flare.  Returns:
            tlevel:      Array of levels at each time, nominally near unity
            bflag:       Array of flags indicating nominal background (where True) or
                            elevated background (where False) indicating possible flare
    '''
    nbl,npol,nf,nt = data.shape
    tlevel = np.zeros(nt,'float')
    background = np.sqrt(np.abs(data[:,0,:,:])**2 + np.abs(data[:,1,:,:])**2)
    init_bg = np.nanmedian(background,2)  # Initially take background as median over entire time range
    bflag = np.ones(nt,'bool')   # flags indicating "good" background times (not in flare)
    for i in range(nt):
        good, = np.where(bflag[:i] == True)   # List of indexes of good background times up to current time
        ngood = len(good)                  # Truncate list of indexes to last 100 elements (or fewer)
        if ngood > 100:
            good = good[ngood-100:]
            # Calculate median over good background times
            bg = np.nanmedian(background[:,:,good],2)
        else:
            # If there haven't been 100 times with good backgrounds yet, just use the initial one.
            # This is supposed to avoid startup transients.
            bg = init_bg
        # Generate levels for each baseline and frequency for this time
        level = np.sqrt(abs(data[:,0,:,i])**2 + abs(data[:,1,:,i])**2)/bg
        # Take median over baseline and frequency to give a single number for this time
        tlevel[i] = np.nanmedian(level)
        if tlevel[i] > 1.05:
            # If the level of the current time is higher than 1.05, do not include this time in future backgrounds
            bflag[i] = False
    return tlevel, bflag

def cleanup(bflag):
    ''' Cleans up the background flag array to remove rapid fluctuations
        and provide better in-flare designations.
    '''
    return bflag
    
def get_history(times, tlevel, bflag):
    ''' Given newly determined tlevel and bflag, see if a file already
        exists for this date and append or replace with new information, 
        if so, otherwise create a new file.
        
        File created is a binary file of records, 9 bytes per record: 
           float time, float tlevel, bool bflag
        Returns data for entire day (contents of any existing file plus
        the new data)
    '''
    import glob
    import dump_tsys
    import struct
    
    datstr = times[0].iso[:10].replace('-','')
    filename = '/common/webplots/flaremon/flaremeter/FLM'+datstr+'.dat'
    if len(glob.glob(filename)) == 1:
        # Filename exists, so read entire file at once
        f = open(filename,'rb')
        buf = f.read()
        f.close()
        nrec = len(buf)/13  # 13 bytes per record: double time, float level, bool flag
        t = np.zeros(nrec,'double')
        l = np.zeros(nrec,'float')
        b = np.zeros(nrec,'bool')
        for i in range(nrec):
            t[i],l[i],b[i] = struct.unpack('dfB',buf[i*13:(i+1)*13])
        # Since unique also sorts, and takes the first instance, it should be enough to
        # concatenate times with t and get unique indexes
        times_lv = np.concatenate(((times.lv+0.5).astype('int'), (t+0.5).astype('int')))
        tlevel = np.concatenate((tlevel, l))
        bflag = np.concatenate((bflag, b))
        blah, idx = np.unique(times_lv,return_index=True)
        times = Time(times_lv[idx],format='lv')
        tlevel = tlevel[idx]
        bflag = bflag[idx]
        
    # Open same filename for writing (overwrites contents if file exists)
    f = open(filename,'wb')
    for i in range(len(times)):
        f.write(struct.pack('dfB',*(times[i].lv,tlevel[i],bflag[i])))
    f.close()
    return times, tlevel, bflag

def rd_RT(filename=None):
    import os
    from util import Time
    if filename is None:
        datstr = Time.now().iso[:10]
        filename = '/data1/RT/RT_'+datstr+'.txt'
    if not os.path.exists(filename):
        print 'File',filename,'not found.'
        return {}
    f = open(filename,'r')
    lines = f.readlines()
    f.close()
    t = []
    amp = []
    tp = []
    from util import Time
    for line in lines:
        t.append(Time(line[:19]).plot_date)
        tp.append(float(line[24:31]))
        amp.append(float(line[32:-3]))
    t = np.array(t)
    dt = t[1:] - t[:-1]
    jmp, = np.where(dt < -0.75)
    if len(jmp) == 1:
        # This is a day jump, so add one to later times
        t[jmp[0]+1:] += 1
    
    return {'time':t,'tp':tp,'amp':amp}

    
if __name__ == "__main__":
    ''' For non-interactive use, use a backend that does not require a display
        Usage python /common/python/current/flare_monitor.py "2014-12-20" <skip>
        If optional argument skip is given, the time-consuming creation of the
        xdata spectrum (XSP file) is skipped.
    '''
    import glob, shutil
    import matplotlib, sys, util
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    known = ['GAIN','PHAS','SOLP']  # known calibration types (first 4 letters)
    t = Time.now()
    skip = False
    if len(sys.argv) >= 2:
        try:
            t = Time(sys.argv[1])
        except:
            if sys.argv[1].upper() == 'RT':
                RT_flare_monitor()
            else:
                print 'Cannot interpret',sys.argv[1],'as a valid date/time string.'
            exit()
        if len(sys.argv) == 3:
            skip = True
    print t.iso[:19],': ',
    # if (t.mjd % 1) < 3./24:
        # # Special case of being run at or before 3 AM (UT), so change to late "yesterday" to finish out
        # # the previous UT day
        # imjd = int(t.mjd)
        # t = Time(float(imjd-0.001),format='mjd')

    tlevel = None
    if not skip:
        # Check if cross-correlation plot already exists
        f, ax = plt.subplots(1,1)
        f.set_size_inches(14,5)
        scanid, tlevel, bflag, times = xdata_display(t,ax)
        if times is None:
            pass
        else:
            plt.savefig('/common/webplots/flaremon/XSP20'+scanid+'.png',bbox_inches='tight')
            plt.close(f)
            print 'Plot written to /common/webplots/flaremon/XSP20'+scanid+'.png'
            bflag = cleanup(bflag)
        # See if a file for this date already exists, and if so, read it and 
        # append or replace with the newly determined levels
        #times, tlevel, bflag = get_history(times, tlevel, bflag)

    # ************ This block commented out due to loss of SQL **************
    projdict = get_projects_nosql(t)
    ut = [Time(projdict['Timestamp'],format='lv').plot_date[0]]*2
    out = rd_RT()
    # ut, fl, projdict = flare_monitor(t)
    # if fl == []:
        # print 'Error retrieving data for',t.iso[:10],'from SQL database.'
        # exit()
    f, ax = plt.subplots(1,1)
    f.set_size_inches(10,3)
    if len(out) == 3:
        ax.plot_date(out['time'],out['amp'],'-',drawstyle='steps')
    # plt.plot_date(ut,fl,'b')
    # if tlevel:
        # plt.plot_date(times.plot_date,tlevel,'r,')
        ax.set_xlabel('Time [UT]')
        ax.set_ylabel('RF Detector [arb. units]')
        ax.set_title('EOVSA Flare Monitor for '+t.iso[:10])
#    ymax = 1.4
    # if np.max(fl) > ymax: ymax = np.max(fl)
    # # Get level max, ignoring nan and inf
    # #lmax = np.max(tlevel[np.isfinite(tlevel)])
    # #if lmax > ymax: ymax = lmax
#    ax.set_ylim(0.8,ymax)
    ax.set_xlim(int(ut[0])+13/24.,int(ut[0])+26/24.)  # Time plot ranges from 13 UT to 02 UT
    if projdict == {}:
        print 'No annotation can be added to plot for',t.iso[:10]
    else:
        nscans = len(projdict['Project'])
        SOS = Time(projdict['Timestamp'],format='lv').plot_date
        EOS = Time(projdict['EOS'],format='lv').plot_date
        yran = np.array(ax.get_ylim())
        for i in range(nscans):
            uti = SOS[i]*np.array([1.,1.])
            plt.plot_date(uti,yran,'g',lw=0.5)
            if projdict['Project'][i] == 'NormalObserving' or projdict['Project'][i] == 'Normal Observing':
                ax.text(uti[0],yran[1]*0.935,'SUN',fontsize=8, clip_on=True)
            elif projdict['Project'][i] == 'None':
                ax.text(uti[0],yran[1]*0.975,'IDLE',fontsize=8, clip_on=True)
            elif projdict['Project'][i][:4] == 'GAIN':
                ax.text(uti[0],yran[1]*0.955,'GCAL',fontsize=8, clip_on=True)
            elif projdict['Project'][i] == 'SOLPNTCAL':
                ax.text(uti[0],yran[1]*0.955,'TPCAL',fontsize=8, clip_on=True)
            elif projdict['Project'][i] == 'PHASECAL':
                ax.text(uti[0],yran[1]*0.955,'PCAL',fontsize=8, clip_on=True)
            else:
                ax.text(uti[0],yran[1]*0.975,projdict['Project'][i],fontsize=8)
        if len(projdict['EOS']) == nscans:
            for i in range(nscans):
                uti = EOS[i]*np.array([1.,1.])
                plt.plot_date(uti,yran,'r--',lw=0.5)
                uti = np.array([SOS[i],EOS[i]])
                if projdict['Project'][i] == 'NormalObserving':
                    plt.plot_date(uti,yran[1]*np.array([0.93,0.93]),ls='-',marker='None',color='#aaffaa',lw=2,solid_capstyle='butt')
                elif projdict['Project'][i][:4] in known:
                    plt.plot_date(uti,yran[1]*np.array([0.95,0.95]),ls='-',marker='None',color='#aaaaff',lw=2,solid_capstyle='butt')
                else:
                    plt.plot_date(uti,yran[1]*np.array([0.97,0.97]),ls='-',marker='None',color='#ffaaaa',lw=2,solid_capstyle='butt')
    datstr = t.iso[:10].replace('-','')
    plt.savefig('/common/webplots/flaremon/FLM'+datstr+'.png',bbox_inches='tight')
    plt.close(f)
    print 'Plot written to /common/webplots/flaremon/FLM'+datstr+'.png'
    # Copy the most recent two files to fixed names so that the web page can find them.
    flist = np.sort(glob.glob('/common/webplots/flaremon/XSP20*'))
    shutil.copy(flist[-1],'/common/webplots/flaremon/XSP_latest.png')
    shutil.copy(flist[-2],'/common/webplots/flaremon/XSP_later.png')

