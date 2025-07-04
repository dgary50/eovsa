''' Routine to read and analyze REFCAL observations to generate REFCAL calibration SQL database.
'''
#
#  2017-05-13 BC
#    Began writing the code
#  2017-06-26 DG
#    Added phacal_anal() routine to plot phase difference wrt refcal, and return
#    the fitted phase slopes and offsets
#  2017-06-27 DG
#    Extended graph_results() to save figures and print string for pasting
#    into the wiki table (if savefigs argument is set to True).  Also updated
#    phacal_anal() to default to not including a phase offset in the fit, but
#    will do so if fitoffsets=True.
#  2017-06-27 BC
#    Added key 'src' into the output dictionary of refcal_anal() to return the 
#    source name
#  2017-06-28 DG
#    Added doplot keyword to refcal_anal(), default is True, to allow the summary 
#    plot of the results to be optional.
#  2017-07-03  DG
#    Added return of HA and DEC keys in rd_refcal(), as keys 'has' and 'decs'.
#    My plan is to do baseline analysis of multi-source data, and these are
#    needed for that.  Also added routine unrot_refcal(), which applies the
#    correction for feed rotation and returns another refcal in the same format.
#    Right now, it will only work for data taken since 2017-07-01.  Also added
#    partially complete routine fit_blerror(), which will use the data from
#    an all-night multi-source, multi-frequency calibration to determine
#    baseline errors.  Currently works for only Bx and By.
#  2017-09-02  DG
#    Skip scans with 0 times, in rd_refcal().  Also read and correct for parallactic
#    angle only for times in a scan, in unrot_refcal().  This improves the time taken 
#    to process a long timerange with many scans.
#  2017-11-10  DG
#    Add hour to savefig output file names, so that multiple times on a given
#    date can be saved.
#  2018-01-08  DG
#    Update unrot_refcal() to read from SQL and apply new xi_rot term.
#  2018-03-05  NK
#    Update refcal_anal() to calculate band 4 phase from lohi scan if lohi=True,
#    and estimate band 4 phase using the most recent lohi=True result if otherwise.
#  2025-05-22  DG
#    NB: I believe this code is superseded by equivalent code in calwidget.py, so
#    it has not been updated to work with 16 antennas.  It doesn't even work for
#    52 bands...
#
import read_idb as ri
from util import Time, ant_str2list, lobe, nearest_val_idx
import numpy as np
import matplotlib.pyplot as plt
from glob import glob
import matplotlib.dates as mdates
import numpy.ma as ma
import os
# import pcapture2 as p
import pdb

# from IPython import embed

# import chan_util_bc as cu

bl2ord = ri.bl2ord


def findfiles(trange, projid='PHASECAL', srcid=None):
    '''identify refcal files
    ***Optional Keywords***
    projid: String--PROJECTID in UFBD records. Default is PHASECAL
    srcid: String--SOURCEID in UFBD records. Can be a string or a list
    '''
    from util import nearest_val_idx
    import struct, time, glob, sys, socket
    import dump_tsys
    fpath = '/data1/eovsa/fits/UDB/' + trange[0].iso[:4] + '/'
    t1 = trange[0].to_datetime()
    t2 = trange[1].to_datetime()
    daydelta = (t2.date() - t1.date()).days
    tnow = Time.now()
    if t1.date() != t2.date():
        # End day is different than start day, so read and concatenate two fdb files
        ufdb = dump_tsys.rd_ufdb(trange[0])
        for ll in xrange(daydelta):
            ufdb2 = dump_tsys.rd_ufdb(Time(trange[0].mjd + ll + 1, format='mjd'))
            if ufdb2:
                for key in ufdb.keys():
                    ufdb.update({key: np.append(ufdb[key], ufdb2[key])})
    else:
        # Both start and end times are on the same day
        ufdb = dump_tsys.rd_ufdb(trange[0])

    if srcid:
        if type(srcid) is str:
            srcid = [srcid]
        sidx_ = np.array([])
        for sid in srcid:
            sidx, = np.where((ufdb['PROJECTID'] == projid) & (ufdb['SOURCEID'] == sid))
            sidx_ = np.append(sidx_, sidx)
        scanidx = np.sort(sidx_).astype('int')
    else:
        scanidx, = np.where(ufdb['PROJECTID'] == projid)
    # List of scan start times
    tslist = Time(ufdb['ST_TS'][scanidx].astype(float).astype(int), format='lv')
    # List of PHASECAL scan end times
    telist = Time(ufdb['EN_TS'][scanidx].astype(float).astype(int), format='lv')

    k = 0  # Number of scans within timerange
    m = 0  # Pointer to first scan within timerange
    flist = []
    status = []
    tstlist = []
    tedlist = []
    srclist = []
    for i in range(len(tslist)):
        if tslist[i].jd >= trange[0].jd and tslist[i].jd <= trange[1].jd:
            flist.append(fpath + ufdb['FILE'][scanidx[i]].astype('str'))
            tstlist.append(tslist[i])
            tedlist.append(telist[i])
            srclist.append(ufdb['SOURCEID'][scanidx[i]])
            k += 1

    if k == 0:
        print 'No scans found within given time range for ' + projid
        return None
    else:
        print 'Found', k, 'scans in timerange.'

    return {'scanlist': flist, 'srclist': srclist, 'tstlist': tstlist, 'tedlist': tedlist}


def rd_refcal(trange, projid='PHASECAL', srcid=None, quackint=180., navg=3):
    '''take a time range from the Time object, e.g., trange=Time(['2017-04-08T05:00','2017-04-08T15:30']),
       a projectid, and source id, return visibility data for all baselines correlated with ant 14.
       ***Optional keywords***
       projid: string -- predefined PROJECTID when setting up the observations. Default is 'PHASECAL'
       srcid: string -- if provided, then only use the specified source. E.g., '1229+020' is often used for
              reference calibration. Default is to use all scans.
       quackint: interval in seconds to skip at the beginning of each scan
       navg: number of data points to average
       ***Output dictionary***

    '''
    import struct, time, glob, sys, socket
    import dbutil as db
    sclist = findfiles(trange, projid, srcid)
    scanlist = sclist['scanlist']
    srclist = sclist['srclist']

    # read scans one by one
    nscan = len(scanlist)
    vis = []
    bandnames = []
    times = []
    fghzs = []
    has = []
    decs = []
    good = []
    for n, scan in enumerate(scanlist):
        print 'Reading scan: ' + scan
        out = ri.read_idb([scan], navg=navg, quackint=quackint)
        nt = len(out['time'])
        if nt == 0:
            # If there are no times in the scan, skip this scan entirely
            continue
        good.append(n)
        times.append(out['time'])
        has.append(out['ha'])
        decs.append(out['dec'])
        bds, sidx = np.unique(out['band'], return_index=True)
        nbd = len(bds)
        eidx = np.append(sidx[1:], len(out['band']))
        vs = np.zeros((15, 4, 34, nt), dtype=complex)
        fghz = np.zeros(34)
        # average over channels within each band
        for b, bd in enumerate(bds):
            fghz[bd - 1] = np.nanmean(out['fghz'][sidx[b]:eidx[b]])
            for a in range(13):
                for pl in range(4):
                    vs[a, pl, bd - 1] = np.mean(out['x'][bl2ord[a, 13], pl, sidx[b]:eidx[b]], axis=0)
        vis.append(vs)
        fghzs.append(fghz)
        bandnames.append(bds)
    # Keep only good scans
    gscanlist = []
    gsrclist = []
    gtstlist = []
    gtedlist = []
    for i in good:
        gscanlist.append(sclist['scanlist'][i])
        gsrclist.append(sclist['srclist'][i])
        gtstlist.append(sclist['tstlist'][i])
        gtedlist.append(sclist['tedlist'][i])
    return {'scanlist': gscanlist, 'srclist': gsrclist, 'tstlist': gtstlist, 'tedlist': gtedlist, 'vis': vis, 'bandnames': bandnames, 'fghzs': fghzs,
            'times': times, 'has': has, 'decs': decs}


def unrot_refcal(refcal_in):
    ''' Apply feed-rotation correction to data read with rd_refcal(), returning updated data in
        the same format for further processing.
    '''
    import dbutil as db
    import copy
    import chan_util_bc as cu
    import cal_header as ch
    from stateframe import extract
    refcal = copy.deepcopy(refcal_in)
    xml, buf = ch.read_cal(11, Time(refcal['times'][0][0], format='jd'))
    dph = extract(buf, xml['XYphase'])
    xi_rot = extract(buf, xml['Xi_Rot'])
    freq = extract(buf, xml['FGHz'])
    freq = freq[np.where(freq != 0)]
    band = []
    for f in freq:
        band.append(cu.freq2bdname(f))
    bds, sidx = np.unique(band, return_index=True)
    nbd = len(bds)
    eidx = np.append(sidx[1:], len(band))
    dxy = np.zeros((14, 34), dtype=np.float)
    xi = np.zeros(34, dtype=np.float)
    fghz = np.zeros(34)
    # average dph and xi_rot frequencies within each band, to convert to 34-band representation
    for b, bd in enumerate(bds):
        fghz[bd - 1] = np.nanmean(freq[sidx[b]:eidx[b]])
        xi[bd - 1] = np.nanmean(xi_rot[sidx[b]:eidx[b]])
        for a in range(14):
            dxy[a, bd - 1] = np.angle(np.sum(np.exp(1j * dph[a, sidx[b]:eidx[b]])))
    nscans = len(refcal['scanlist'])
    for i in range(nscans):
        # Read parallactic angles for this scan
        trange = Time([refcal['tstlist'][i].iso, refcal['tedlist'][i].iso])
        times, chi = db.get_chi(trange)
        tchi = times.jd
        t = refcal['times'][i]
        if len(t) > 0:
            vis = copy.deepcopy(refcal['vis'][i])
            idx = nearest_val_idx(t, tchi)
            pa = chi[idx]  # Parallactic angle for the times of this refcal.
            pa[:, [8, 9, 10, 12]] = 0.0
            nt = len(idx)  # Number of times in this refcal
            # Apply X-Y delay phase correction
            for a in range(13):
                a1 = lobe(dxy[a] - dxy[13])
                a2 = -dxy[13] - xi
                a3 = dxy[a] - xi + np.pi
                for j in range(nt):
                    vis[a, 1, :, j] *= np.exp(1j * a1)
                    vis[a, 2, :, j] *= np.exp(1j * a2)
                    vis[a, 3, :, j] *= np.exp(1j * a3)
            for j in range(nt):
                for a in range(13):
                    refcal['vis'][i][a, 0, :, j] = vis[a, 0, :, j] * np.cos(pa[j, a]) + vis[a, 3, :, j] * np.sin(pa[j, a])
                    refcal['vis'][i][a, 2, :, j] = vis[a, 2, :, j] * np.cos(pa[j, a]) + vis[a, 1, :, j] * np.sin(pa[j, a])
                    refcal['vis'][i][a, 3, :, j] = vis[a, 3, :, j] * np.cos(pa[j, a]) - vis[a, 0, :, j] * np.sin(pa[j, a])
                    refcal['vis'][i][a, 1, :, j] = vis[a, 1, :, j] * np.cos(pa[j, a]) - vis[a, 2, :, j] * np.sin(pa[j, a])
    return refcal


def graph(out, refcal=None, ant_str='ant1-13', bandplt=[5, 11, 17, 23], scanidx=None, pol=0, tformat='%H:%M'):
    '''Produce a figure showing phases and amplitudes of selected antennas, bands, and polarization.
       Optionally takes in time averaged 'refcal' (can be from refcal_anal() or sql database) 
       to show the selected time range and plot the averaged phase and amplitudes'''
    # takes input from rd_refcal (out)
    date_format = mdates.DateFormatter(tformat)
    # make a color plot showing the phase
    # ri.summary_plot_pcal(out)
    if scanidx:
        scanlist = [out['scanlist'][i] for i in scanidx]
        tstlist = [out['tstlist'][i] for i in scanidx]
        tedlist = [out['tedlist'][i] for i in scanidx]
        bandnames = [out['bandnames'][i] for i in scanidx]
        vis = [out['vis'][i] for i in scanidx]
        times = [out['times'][i] for i in scanidx]
    else:
        scanlist = out['scanlist']
        tstlist = out['tstlist']
        tedlist = out['tedlist']
        bandnames = out['bandnames']
        vis = out['vis']
        times = out['times']

    nscan = len(scanlist)
    ant_list = ant_str2list(ant_str)
    nant = len(ant_list)
    bds0 = bandplt
    nband = len(bds0)
    f1, ax1 = plt.subplots(nant, nband, figsize=(12, 8))
    f2, ax2 = plt.subplots(nant, nband, figsize=(12, 8))
    for n, scan in enumerate(scanlist):
        ts = Time(times[n], format='jd').plot_date
        try:
            # make plots of phase and amp vs. time on all 34 bands
            for ant in ant_list:
                for b, bd in enumerate(bds0):
                    ph_deg = np.angle(vis[n][ant, pol, bd - 1], deg=True)
                    amp = np.abs(vis[n][ant, pol, bd - 1])
                    ax1[ant, b].plot_date(ts, ph_deg, '.', markersize=5)
                    ax1[ant, b].xaxis.set_major_formatter(date_format)
                    ax2[ant, b].plot_date(ts, amp, '.', markersize=5)
                    ax2[ant, b].xaxis.set_major_formatter(date_format)
                    ax1[ant, b].set_ylim([-180., 180.])
                    if ant == 0:
                        ax1[ant, b].set_title('Band ' + str(bd))
                        ax2[ant, b].set_title('Band ' + str(bd))
                    if b == 0:
                        ax1[ant, b].text(-0.4, 0.5, 'Ant ' + str(ant + 1), ha='center', va='center', transform=ax1[ant, b].transAxes, fontsize=10)
                        ax2[ant, b].text(-0.4, 0.5, 'Ant ' + str(ant + 1), ha='center', va='center', transform=ax2[ant, b].transAxes, fontsize=10)
                    else:
                        ax1[ant, b].set_yticks([])
                        ax2[ant, b].set_yticks([])
        except:
            print 'Failure in plotting scan: ', scan
            continue
    if refcal:
        for ant in ant_list:
            for b, bd in enumerate(bds0):
                phavg = np.degrees(refcal['pha'][ant, pol, bd - 1])
                ampavg = refcal['amp'][ant, pol, bd - 1]
                flg = refcal['flag'][ant, pol, bd - 1]
                if flg == 0:
                    if 't_bg' in refcal.keys():
                        bt = Time(refcal['t_bg'], format='jd').plot_date
                    else:
                        bt = Time(times[0][0], format='jd').plot_date
                    if 't_ed' in refcal.keys():
                        et = Time(refcal['t_ed'], format='jd').plot_date
                    else:
                        et = Time(times[-1][-1], format='jd').plot_date
                    ax1[ant, b].plot_date([bt, et], [phavg, phavg], '-r')
                    ax1[ant, b].plot_date([bt, bt], [-200, 200], '--k')
                    ax1[ant, b].plot_date([et, et], [-200, 200], '--k')
                    ax2[ant, b].plot_date([bt, et], [ampavg, ampavg], '-r')
                    ax2[ant, b].plot_date([bt, bt], [0, np.max(amp)], '--k')
                    ax2[ant, b].plot_date([et, et], [0, np.max(amp)], '--k')
                    ax1[ant, b].plot_date(refcal['timestamp'].plot_date, phavg, 'o')
                    ax2[ant, b].plot_date(refcal['timestamp'].plot_date, ampavg, 'o')


def refcal_anal(out, timerange=None, scanidx=None, minsnr=0.7, bandplt=[5, 11, 17, 23], doplot=True, lohi=False):
    '''Analyze the visibility data from rd_refcal and return time averaged visibility values and flags.
       ***Optional Keywords***
       timerange: time range to obtain the average. E.g., timerange=Time(['2017-04-08T05:00','2017-04-08T07:00'])
       scanidx: index numbers of scans to select. Useful when other (undesired) types of observations exist. 
               !!!! The selected scans should have exactly the same frequency/band setup !!!!
       minsnr: minimum signal to noise to consider. Data with smaller SNRs will be flagged (as 1 in the flag array)
       bandplt: bands to show in the figure
       doplot: if True, display plots of results (default)
       ***Outputs***
       refcal: complex array of shape (15, 2, 34) (nant, npol, nband) as the result of the reference calibration
       flag: int array of shape (15, 2, 34). 0 is unflagged and 1 is flagged.
       timestamp: midpoint of the time range used for averaging to obtain the refcal values 
    '''
    if scanidx:
        scanlist = [out['scanlist'][i] for i in scanidx]
        srclist = [out['srclist'][i] for i in scanidx]
        tstlist = [out['tstlist'][i] for i in scanidx]
        tedlist = [out['tedlist'][i] for i in scanidx]
        bandnames = [out['bandnames'][i] for i in scanidx]
        vis = [out['vis'][i] for i in scanidx]
        times_ = [out['times'][i] for i in scanidx]
        fghzs = [out['fghzs'][i] for i in scanidx]
    else:
        scanlist = out['scanlist']
        srclist = out['srclist']
        tstlist = out['tstlist']
        tedlist = out['tedlist']
        bandnames = out['bandnames']
        vis = out['vis']
        times_ = out['times']
        fghzs = out['fghzs']

    scanidx = range(len(scanlist))
    if len(set(np.array(srclist)[scanidx])) > 1:
        prompt = ''
        while not (prompt.lower() in ['y', 'n']):
            prompt = raw_input('Multiple sources are selected. Are you sure to continue? [y/n]')
        if prompt.lower() == 'n':
            print 'Abort...'
            return None

    fghz = fghzs[0]
    times = np.concatenate(times_)
    vis = np.concatenate(vis, axis=3)
    # only keep the first 2 polarizations
    vis = vis[:, :2]
    if timerange:
        tidx, = np.where((times > timerange[0].jd) & (times < timerange[1].jd))
        if len(tidx) == 0:
            print 'no records within the selected timerange. Abort...'
        for i in range(len(scanlist)):
            sidx, = np.where((times_[i] > timerange[0].jd) & (times_[i] < timerange[1].jd))
            if len(sidx) > 0:
                src = srclist[i]
                break
        vis = vis[:, :, :, tidx]
        timeavg = times[tidx]
    else:
        timeavg = times
        src = srclist[0]
    # vismean = np.nanmean(np.angle(vis),axis=3)
    vis_ = np.zeros(vis.shape[:3], dtype=complex)
    flag = np.zeros(vis.shape[:3], dtype=int)
    sigma = np.zeros(vis.shape[:3]) + 1e10
    # compute standard deviation of the visibilities
    sigma_ = np.nanstd(vis, axis=3)
    # sigma = np.nanstd(np.abs(vis),axis=3)
    # mask out records with phases > 3 sigma
    for bd in range(34):
        # print 'band: ',bd
        for ant in range(13):
            # print 'ant: ',ant
            for pol in range(2):
                # print 'pol: ',pol
                amp = np.abs(vis[ant, pol, bd])
                amp_median = np.median(np.abs(vis[ant, pol, bd]))
                ind, = np.where(np.abs(amp - amp_median) < sigma_[ant, pol, bd])
                snr = amp_median / sigma_[ant, pol, bd]
                # pdb.set_trace()
                if snr < minsnr or np.isnan(snr):
                    flag[ant, pol, bd] = 1
                else:
                    if len(ind) > len(timeavg) / 2:
                        vis_[ant, pol, bd] = np.nanmean(vis[ant, pol, bd, ind])
                        sigma[ant, pol, bd] = np.nanstd(vis[ant, pol, bd, ind])
                    else:
                        vis_[ant, pol, bd] = np.nanmean(vis[ant, pol, bd])
                        flag[ant, pol, bd] = 1
                        # print '# of valid datapoints: ',len(ind)
        # count how many datapoints are flagged in a given band
        nflag = np.count_nonzero(flag[:13, :, bd])
        print '{0:d} of 26 measurements are flagged due to SNR < {1:.1f} in Band {2:d}'.format(nflag, minsnr, bd + 1)
    # flag zero values (e.g., antennas or bands not observed)
    zeroind = np.where(np.abs(vis_ == 0))
    flag[zeroind] = 1
    # timestamps
    timestamp = Time(np.mean(timeavg), format='jd')
    timestamp_gcal = Time((tstlist[0].jd + tedlist[0].jd) / 2., format='jd')
    
    #Calculate band 4 phases
    refcal_lohi = sql2refcal(timestamp, lohi=True)  #obtain from SQL database
    dph = np.zeros((15,2,30))
    for i in range(13):
         for j in range(2):
             dph[i,j,:] = np.unwrap(lobe(np.angle(vis_)[i,j,4:] - refcal_lohi['pha'][i,j,4:]))
                 
    dphfitxxyy = np.zeros((13,2,2))
    w = np.ones(31)
    w[0] += 99
    for i in range(13):
         for j in range(2):
             dphfitxxyy[i,j,:] = np.polyfit(np.append([0.],fghz[4:]),np.append([0.],[dph[i,j]]),1,w=w)
    
    fghz[3] = refcal_lohi['fghz'][3]
    dph4 = np.zeros((13,2))
    for i in range(13):
         for j in range(2):
             dph4[i,j] = np.polyval(dphfitxxyy[i,j,:],fghz[3])        
    dph4 = dph4
    for i in range(13):  #Insert band 4 phase in high frequency receiver
        for j in range(2):
            vis_[i,j,3] = np.complex(np.cos(refcal_lohi['pha'][i,j,3] + dph4[i,j]), np.sin(refcal_lohi['pha'][i,j,3] + dph4[i,j]))  #Amp = 1.0
    
    if doplot:
        visavg = {'pha': np.angle(vis_), 'amp': np.abs(vis_), 'timestamp': timestamp, 't_bg': timeavg[0], 't_ed': timeavg[-1], 'flag': flag}
        graph(out, visavg, scanidx=scanidx, bandplt=bandplt)
        graph(out, visavg, scanidx=scanidx, bandplt=bandplt, pol=1)
        f2, ax2 = plt.subplots(2, 13, figsize=(12, 5))
        if lohi == 0:
            f2.suptitle('Orange = Estimated band 4 phases based on lohi scans on ' + str(refcal_lohi['timestamp'].iso))
        plt.title('source: {}'.format(srclist[scanidx[0]]))
        f3, ax3 = plt.subplots(2, 13, figsize=(12, 5))
        plt.title('source: {}'.format(srclist[scanidx[0]]))
        allbands = np.arange(34) + 1
        for ant in range(13):
            for pol in range(2):
                ind, = np.where(flag[ant, pol, :] == 0)
                ax2[pol, ant].plot(allbands[ind], np.unwrap(visavg['pha'][ant, pol, ind]), '.', markersize=5)
                if lohi == 0:
                    ax2[pol, ant].plot(4., np.unwrap(visavg['pha'])[ant, pol, 3], '.', markersize=5)  #Add band 4
                ax2[pol, ant].set_ylim([-20, 20])
                ax2[pol, ant].set_xlim([1, 34])
                ax3[pol, ant].plot(allbands[ind], visavg['amp'][ant, pol, ind], '.', markersize=5)
                ax3[pol, ant].set_xlim([1, 34])
                ax3[pol, ant].set_ylim([0, 1.])
                if ant == 0:
                    ax2[pol, ant].set_ylabel('Phase (radian)')
                    ax3[pol, ant].set_ylabel('Amplitude')
                else:
                    ax2[pol, ant].set_yticks([])
                    ax3[pol, ant].set_yticks([])
                if pol == 1 and ant == 0:
                    ax2[pol, ant].set_xlabel('Band #')
                    ax3[pol, ant].set_xlabel('Band #')
                if pol == 0:
                    ax2[pol, ant].set_title('Ant ' + str(ant + 1))
                    ax2[pol, ant].set_xticks([])
                    ax3[pol, ant].set_title('Ant ' + str(ant + 1))
                    ax3[pol, ant].set_xticks([])
                    
        if lohi == 0:
            f, ax = plt.subplots(2,13,figsize=(12,5))  #Plot delay curve used for band 4 phase estimation
            f.suptitle('Phase differences with respect to ' + str(refcal_lohi['timestamp'].iso))
            for i in range(13):
                 ax[0,i].set_title('Ant ' + str(i+1))
                 for j in range(2):
                     ax[j,i].plot(fghz[4:],dph[i,j,:],'.')
                     ax[j,i].plot(np.append([0.],fghz[4:]),np.polyval(dphfitxxyy[i,j,:],np.append([0.],fghz[4:])))
                     ax[j,i].plot(fghz[3],np.polyval(dphfitxxyy[i,j,:],fghz[3]),'.')
                     ax[j,i].set_xlim([0.,18.])
                     ax[j,i].set_ylim([-np.pi,np.pi])
                     ax[0,0].set_ylabel('XX')
                     ax[1,0].set_ylabel('YY')
                     ax[1,0].set_xlabel('fghz')

    if lohi:
        fghz = out['fghzs'][0]
        fghz[3] = out['fghzs'][1][3]
        bandnames = np.append(4,out['bandnames'][0])
        com_idx = [4,5,6,7,8,9,10,11]  #because common_val_idx somehow failed
        phlo = np.angle(np.sum(out['vis'][1],3))[:,:2,:]
        phhi = np.angle(vis_)
        
        #Caliculate band 4 phase
        dphlohi = np.unwrap(lobe(phlo[:,:,com_idx] - phhi[:,:,com_idx]))  #dph over bands 5-12
        dphfitlohi = np.zeros((13,2,3))
        for i in range(13):  #Poly-fit dph curve
             for j in range(2):
                 dphfitlohi[i,j,:] = np.polyfit(out['fghzs'][0][com_idx],dphlohi[i,j,:],2)
        dph4 = np.zeros((13,2))  #Extrapolate the fit to band 4
        for i in range(13):
             for j in range(2):
                 dph4[i,j] = np.polyval(dphfitlohi[i,j,:],out['fghzs'][1][3])
        for i in range(13):  #Insert band 4 phase in high frequency receiver
            for j in range(2):
                vis_[i,j,3] = np.complex(np.cos(phlo[i,j,3] - dph4[i,j]), np.sin(phlo[i,j,3] - dph4[i,j]))  #Amp = 1.0
        phhi_new = np.angle(vis_)
        
        #Plot band 4 phase
        f, ax = plt.subplots(2,13,figsize=(12,5))
        f.suptitle('Calculated band 4 phases (orange)')
        plt.title('source: {}'.format(srclist[scanidx[0]]))
        for ant in range(13):
             ax[0,ant].set_title('Ant ' + str(ant+1))
             for pol in range(2):
                 ax[pol,ant].plot(bandnames[1:], np.unwrap(phhi_new)[ant,pol,4:],'.')
                 ax[pol,ant].plot(bandnames[0], np.unwrap(phhi_new)[ant,pol,3],'.')                                  
                 ax[pol,ant].set_ylim([-20,20])
                 ax[pol,ant].set_xlim([1,34])
                 if pol == 0: ax[pol,ant].set_xticks([])
                 if ant >= 1: ax[pol,ant].set_yticks([])
        ax[1,0].set_xlabel('Band #')
        ax[0,0].set_ylabel('XX Phase (radian)')
        ax[1,0].set_ylabel('YY Phase (radian)')
        
        return {'src': src, 'vis': vis_, 'pha': np.angle(vis_), 'amp': np.abs(vis_), 'fghz': fghz, 'flag': flag,
                'sigma': sigma, 'timestamp': timestamp, 't_gcal': timestamp_gcal, 't_bg': Time(timeavg[0], format='jd'),
                't_ed': Time(timeavg[-1], format='jd')}

    return {'src': src, 'vis': vis_, 'pha': np.angle(vis_), 'amp': np.abs(vis_), 'fghz': fghz, 'flag': flag, 'sigma': sigma, 'timestamp': timestamp,
            't_gcal': timestamp_gcal, 't_bg': Time(timeavg[0], format='jd'), 't_ed': Time(timeavg[-1], format='jd')}


def graph_results(refcal, unwrap=True, savefigs=False):
    '''Provide an output from sql2refcal() (single refcal result) or sql2refcalX() (list of refcal results).
       Plot phase and amp vs. bands
       
       If savefigs is True, creates plot files and prints the entry needed for pasting into the wiki page.
    '''
    f1, ax1 = plt.subplots(2, 13, figsize=(12, 5))
    f2, ax2 = plt.subplots(2, 13, figsize=(12, 5))
    allbands = np.arange(34) + 1
    bad = ''
    for ant in range(13):
        for pol in range(2):
            if type(refcal) is dict:
                refcal = [refcal]
            for ref in refcal:
                ind, = np.where(ref['flag'][ant, pol, :] == 0)
                if unwrap:
                    pha = np.unwrap(ref['pha'][ant, pol, ind])
                    ax1[pol, ant].set_ylim([-20, 20])
                else:
                    pha = ref['pha'][ant, pol, ind]
                    ax1[pol, ant].set_ylim([-4, 4])
                ax1[pol, ant].plot(allbands[ind], pha, '.', markersize=5)
                ax1[pol, ant].set_xlim([1, 34])
                ax2[pol, ant].plot(allbands[ind], ref['amp'][ant, pol, ind], '.', markersize=5)
                ax2[pol, ant].set_xlim([1, 34])
                ax2[pol, ant].set_ylim([0, 1.])
                if len(ind) == 0:
                    ax1[pol, ant].text(17, 0, 'No Cal', ha='center')
                    ax2[pol, ant].text(17, 0.5, 'No Cal', ha='center')
                    if pol == 0:
                        bad += ' Ant ' + str(ant + 1)

                if ant == 0:
                    ax1[pol, ant].set_ylabel('Phase (radian)')
                    ax2[pol, ant].set_ylabel('Amplitude')
                else:
                    ax1[pol, ant].set_yticks([])
                    ax2[pol, ant].set_yticks([])
                if pol == 1 and ant == 0:
                    ax1[pol, ant].set_xlabel('Band #')
                    ax2[pol, ant].set_xlabel('Band #')
                if pol == 0:
                    ax1[pol, ant].set_title('Ant ' + str(ant + 1))
                    ax1[pol, ant].set_xticks([])
                    ax2[pol, ant].set_title('Ant ' + str(ant + 1))
                    ax2[pol, ant].set_xticks([])
    if savefigs:
        dstr = ref['timestamp'].iso[:10]
        tstr = ref['timestamp'].iso[11:19]
        file1 = dstr.replace('-', '') + '_' + tstr[:2] + '_refcal_pha.png'
        f1.savefig('/common/webplots/refcal/' + file1)
        file2 = dstr.replace('-', '') + '_' + tstr[:2] + '_refcal_amp.png'
        f2.savefig('/common/webplots/refcal/' + file2)
        maxlen = 0
        idx, = np.where(ref['flag'][0, 0] == 0)
        for ant in range(13):
            for pol in range(2):
                ind, = np.where(ref['flag'][ant, pol] == 0)
                if len(ind) > maxlen:
                    idx = ind
                    maxlen = len(idx)
        bstr = str(idx[0] + 1) + '~' + str(idx[-1] + 1)
        if bad != '':
            badstr = 'No calibration for:' + bad
        else:
            badstr = ''
        print '|', dstr.replace('-', '/'), '||', tstr, '|| ', ref[
            'src'], ' || || 0 ||  ||', bstr, '|| [http://ovsa.njit.edu/refcal/' + file1, 'Phase] || [http://ovsa.njit.edu/refcal/' + file2, 'Amp] ||', badstr
        print '|-'


def mbdfunc0(fghz, mbd):
    # fghz: frequency in GHz
    # ph0 = 0: phase offset identically set to zero (not fitted)
    # mbd: multi-band delay associated with the phase_phacal - phase_refcal in ns
    return 2. * np.pi * fghz * mbd


def mbdfunc1(fghz, ph0, mbd):
    # fghz: frequency in GHz
    # ph0: phase offset in radians
    # mbd: multi-band delay associated with the phase_phacal - phase_refcal in ns
    return ph0 + 2. * np.pi * fghz * mbd


def phase_diff(phacal, refcal=None, strictness=0.5, fitoffsets=False, verbose=False):
    '''Fit the phase difference between a phase calibration (or another refcal)
       and the reference calibration.  Returns the phase slopes and, optionally,
       the offsets, along with the relevant times, as a dictionary.
       strictness from 0 to 1 gives the increasing strictness of quality control over the phase_diff solution.
       
       Returns:
       pfit    A python dictionary with the following keys:
           't_pha' : Central time of the phase calibration data
           't_ref' : Central time of the reference calibration data
           'poff'  : Phase offsets (zero if fitoffsets=False) for each antenna, polarization
           'pslope': Phase slope for each antenna, polarization
           'flag'  : Indication of bad data for each antenna, polarization (0 = good, 1 = bad)
           'phacal': Actual instance of phase calibration dictionary, containing frequency-dependent
                       amplitudes, phases, etc.
    '''
    from scipy.optimize import curve_fit

    t_pha = phacal['timestamp']
    if refcal is None:
        refcal = sql2refcal(t_pha)
    t_ref = refcal['timestamp']
    src = phacal['src']
    dpha = phacal['pha'] - refcal['pha']
    flag_pha = phacal['flag']
    flag_ref = refcal['flag']
    nants = 15

    poff = [[], []]
    pslope = [[], []]
    prms = [[], []]
    flag = [[], []]
    for ant in range(nants):
        if verbose: print 'ant: ', ant
        for pol in range(2):
            if ant < 13:
                if verbose: print 'pol: ', pol
                ind, = np.where((flag_pha[ant, pol] == 0) & (flag_ref[ant, pol] == 0))
                dpha_unw = np.unwrap(dpha[ant, pol, ind])
                fghz = phacal['fghz'][ind]
                if len(fghz) > 3:
                    # Ensure that offset is close to zero, modulo 2*pi
                    if dpha_unw[0] > np.pi: dpha_unw -= 2 * np.pi
                    if dpha_unw[0] < -np.pi: dpha_unw += 2 * np.pi
                sig_pha = phacal['sigma'][ant, pol, ind]
                sig_ref = refcal['sigma'][ant, pol, ind]
                amp_pha = phacal['amp'][ant, pol, ind]
                amp_ref = refcal['amp'][ant, pol, ind]
                sigma = ((sig_pha / amp_pha) ** 2. + (sig_ref / amp_ref) ** 2.) ** 0.5
                # Do a linear fit on the residual phases
                if len(fghz) > 3:
                    if fitoffsets:
                        popt, pcov = curve_fit(mbdfunc1, fghz, dpha_unw, p0=[0., 0.], sigma=sigma, absolute_sigma=False)
                        poff[pol].append(popt[0])
                        pslope[pol].append(popt[1])
                        flag[pol].append(0)
                        residuals = mbdfunc1(fghz, *popt) - dpha_unw
                        if verbose: print 'Phase offset (deg):', np.degrees(popt[0])
                        if verbose: print 'MBD (ns):', popt[1]
                    else:
                        popt, pcov = curve_fit(mbdfunc0, fghz, dpha_unw, p0=[0.], sigma=sigma, absolute_sigma=False)
                        poff[pol].append(0.0)
                        pslope[pol].append(popt[0])
                        flag[pol].append(0)
                        residuals = mbdfunc0(fghz, *popt) - dpha_unw
                        if verbose:
                            print 'Phase offset (deg): 0.0 (not fit)'
                        if verbose:
                            print 'MBD (ns):', popt[0]
                    rms = np.sqrt(np.dot(residuals, residuals) / len(residuals))
                    prms[pol].append(rms)
                    if rms > 1.0:
                        hist, bins = np.histogram(np.abs(residuals), bins=len(residuals))
                        bins = ma.masked_greater(bins[1:], 1.0)
                        hist = ma.masked_array(hist, mask=bins.mask)
                        if np.sum(hist, dtype=np.float) / np.sum(hist.data, dtype=np.float) < strictness:
                            flag[pol][-1] = 1
                else:
                    poff[pol].append(0.0)
                    pslope[pol].append(0.0)
                    flag[pol].append(1)
                    prms[pol].append(0.0)
            else:
                poff[pol].append(0.0)
                pslope[pol].append(0.0)
                flag[pol].append(1)
                prms[pol].append(0.0)

    antstr = lambda ant: ' Ant ={:6s}                 '.format(ant)
    titlestr = ' '.join(['{:12s}{:2s}'.format('rms', 'F')] * 2)
    sepstr = '{0}|{0}'.format('-' * 14)
    caltbstr = lambda rms, flg: '{:10.5f}  {}  {:10.5f}  {} '.format(rms[0], flg[0], rms[1], flg[1])

    ncols = 4
    nrows = np.int(np.ceil(nants / 4.0))
    print('------------------------------------------------------------------------------------------------------------------------')
    print('PHASECAL quality assessment (rms in degree) ----- Field = {0:10s}, Time = {1}~{2}'.format(src, phacal['t_bg'].iso[:-4],
                                                                                                     phacal['t_ed'].iso[:-4]))
    print('------------------------------------------------------------------------------------------------------------------------')
    for row in range(nrows):
        if row < nrows - 1:
            ants = range(ncols * row, ncols * (row + 1))
        else:
            ants = range(ncols * row, nants)
        print '|'.join(map(antstr, ['eo{:02d}'.format(ll + 1) for ll in ants])) + '|'
        print '|'.join([titlestr] * len(ants)) + '|'
        print '|'.join([sepstr] * ncols) + '|'
        rms = []
        flg = []
        for ant in ants:
            rms.append([prms[pol][ant] for pol in range(2)])
            fg = []
            for pol in range(2):
                if flag[pol][ant] == 1:
                    fg.append(flag[pol][ant])
                else:
                    fg.append(' ')
            flg.append(fg)
        # rms = [[prms[pol][ant] for pol in range(2)] for ant in ants]
        print ' '.join(map(caltbstr, rms, flg))

    return {'t_pha': t_pha, 't_ref': t_ref, 'poff': np.transpose(poff), 'pslope': np.transpose(pslope), 'prms': np.transpose(prms),
            'flag': np.transpose(flag), 'phacal': phacal}


def graph_pdiff(p_diff, refcal, strictness=0.5, verbose=False, plot_rms=False):
    '''Produce a figure showing phase difference from phase_diff().
    strictness from 0 to 1 gives the increasing strictness of quality control over the phase_diff solution.'''
    t_pha = p_diff['t_pha']
    t_ref = p_diff['t_ref']
    phacal = p_diff['phacal']
    poff = p_diff['poff']
    pslope = p_diff['pslope']
    dpha = phacal['pha'] - refcal['pha']
    flag_pha = phacal['flag']
    flag_ref = refcal['flag']

    if plot_rms:
        pflag = p_diff['flag']
        f, ax = plt.subplots(4, 13, figsize=(13, 5))
        f.suptitle('Phase Diff residuals vs. Freq from ' + t_pha.iso + ', relative to ' + t_ref.iso)
        for ant in range(15):
            if verbose: print 'ant: ', ant
            for pol in range(2):
                if ant < 13:
                    if verbose: print 'pol: ', pol
                    ind, = np.where((flag_pha[ant, pol] == 0) & (flag_ref[ant, pol] == 0))
                    dpha_unw = np.unwrap(dpha[ant, pol, ind])
                    fghz = phacal['fghz'][ind]
                    if len(fghz) > 3:
                        # Ensure that offset is close to zero, modulo 2*pi
                        if dpha_unw[0] > np.pi: dpha_unw -= 2 * np.pi
                        if dpha_unw[0] < -np.pi: dpha_unw += 2 * np.pi
                    if len(fghz) > 3:
                        residuals = dpha_unw - mbdfunc1(fghz, poff[ant, pol], pslope[ant, pol])
                        # rms = np.sqrt(np.dot(residuals, residuals) / len(residuals))
                        ax[pol, ant].plot(fghz, residuals, '.k')
                        ax[pol, ant].set_ylim([-10, 10])
                        ax[pol, ant].set_xlim([0, 18])
                        ax[pol, ant].axhline(0.0, ls='--', color='r')
                        ax[pol, ant].text(9, 8, 'Ant ' + str(ant + 1), ha='center')
                        nres = len(residuals)
                        if pflag[ant, pol] == 1:
                            ax[pol + 2, ant].hist(np.abs(residuals), bins=nres, cumulative=True, range=[0, 5], color='r')
                        else:
                            ax[pol + 2, ant].hist(np.abs(residuals), bins=nres, cumulative=True, range=[0, 5], color='royalblue')
                        ax[pol + 2, ant].axvline(1, ls=':', c='orange')
                        ax[pol + 2, ant].axhline(strictness * nres, ls=':', c='orange')
                    else:
                        ax[pol, ant].text(9, 8, 'Ant ' + str(ant + 1), ha='center')
                        ax[pol, ant].text(9, 0, 'No Cal', ha='center')
                    ax[pol, ant].text(9, -8, 'Flag {}'.format(pflag[ant, pol]), ha='center')
                    ax[pol + 2, ant].text(0.1, 0.2, 'Flag {}'.format(pflag[ant, pol]), ha='left', transform=ax[pol + 2, ant].transAxes)
                    ax[pol + 2, ant].set_xlim([0, 5])
        ax[0, 0].set_ylabel('Phase Diff residuals [rad]')
        ax[2, 0].set_ylabel('Cumulative Histogram')
        for i in range(13):
            ax[1, i].set_xlabel('f [GHz]')
            if i != 0:
                for j in range(2): ax[j, i].set_yticklabels([])
    else:
        f, ax = plt.subplots(2, 13, figsize=(13, 5))
        f.suptitle('Phase Diff vs. Freq from ' + t_pha.iso + ', relative to ' + t_ref.iso)
        for ant in range(15):
            if verbose: print 'ant: ', ant
            for pol in range(2):
                if ant < 13:
                    if verbose: print 'pol: ', pol
                    ind, = np.where((flag_pha[ant, pol] == 0) & (flag_ref[ant, pol] == 0))
                    dpha_unw = np.unwrap(dpha[ant, pol, ind])
                    fghz = phacal['fghz'][ind]
                    if len(fghz) > 3:
                        # Ensure that offset is close to zero, modulo 2*pi
                        if dpha_unw[0] > np.pi: dpha_unw -= 2 * np.pi
                        if dpha_unw[0] < -np.pi: dpha_unw += 2 * np.pi
                    if len(fghz) > 3:
                        ax[pol, ant].plot(fghz, dpha_unw, '.k')
                        ax[pol, ant].set_ylim([-10, 10])
                        ax[pol, ant].set_xlim([0, 18])
                        ax[pol, ant].plot(fghz, mbdfunc1(fghz, poff[ant, pol], pslope[ant, pol]), 'r--')
                        ax[pol, ant].text(9, 8, 'Ant ' + str(ant + 1), ha='center')
                    else:
                        ax[pol, ant].text(9, 8, 'Ant ' + str(ant + 1), ha='center')
                        ax[pol, ant].text(9, 0, 'No Cal', ha='center')
        for j in range(2): ax[j, 0].set_ylabel('Phase Diff [rad]')
        for i in range(13):
            ax[1, i].set_xlabel('f [GHz]')
            if i != 0:
                for j in range(2): ax[j, i].set_yticklabels([])


def sql2refcal(t, lohi=False):
    '''Supply a timestamp in Time format, return the closest refcal data'''
    import cal_header as ch
    import stateframe as stf
    if lohi:
        caltype = 12
    else:
        caltype = 8
    xml, buf = ch.read_cal(caltype, t=t)
    refcal = stf.extract(buf, xml['Refcal_Real']) + stf.extract(buf, xml['Refcal_Imag']) * 1j
    flag = stf.extract(buf, xml['Refcal_Flag'])
    fghz = stf.extract(buf, xml['Fghz'])
    sigma = stf.extract(buf, xml['Refcal_Sigma'])
    timestamp = Time(stf.extract(buf, xml['Timestamp']), format='lv')
    tbg = Time(stf.extract(buf, xml['T_beg']), format='lv')
    ted = Time(stf.extract(buf, xml['T_end']), format='lv')
    pha = np.angle(refcal)
    amp = np.absolute(refcal)
    return {'pha': pha, 'amp': amp, 'flag': flag, 'fghz': fghz, 'sigma': sigma, 'timestamp': timestamp, 't_bg': tbg, 't_ed': ted}


def sql2refcalX(trange, lohi=False, *args, **kwargs):
    '''same as sql2refcal. trange can be either a timestamp or a timerange.'''
    import cal_header as ch
    import stateframe as stf
    if lohi:
        caltype = 12
    else:
        caltype = 8
    xml, bufs = ch.read_calX(caltype, t=trange, *args, **kwargs)
    if isinstance(bufs, list):
        refcals = []
        for i, buf in enumerate(bufs):
            try:
                ref = stf.extract(buf, xml['Refcal_Real']) + stf.extract(buf, xml['Refcal_Imag']) * 1j
                flag = stf.extract(buf, xml['Refcal_Flag'])
                fghz = stf.extract(buf, xml['Fghz'])
                sigma = stf.extract(buf, xml['Refcal_Sigma'])
                timestamp = Time(stf.extract(buf, xml['Timestamp']), format='lv')
                tbg = Time(stf.extract(buf, xml['T_beg']), format='lv')
                ted = Time(stf.extract(buf, xml['T_end']), format='lv')
                pha = np.angle(ref)
                amp = np.absolute(ref)
                refcals.append({'pha': pha, 'amp': amp, 'flag': flag, 'fghz': fghz, 'sigma': sigma, 'timestamp': timestamp, 't_bg': tbg, 't_ed': ted})
            except:
                print 'failed to load record {} ---> {}'.format(i + 1, Time(stf.extract(buf, xml['Timestamp']), format='lv').iso)
        return refcals
    elif isinstance(bufs, str):
        refcal = stf.extract(bufs, xml['Refcal_Real']) + stf.extract(bufs, xml['Refcal_Imag']) * 1j
        flag = stf.extract(bufs, xml['Refcal_Flag'])
        fghz = stf.extract(bufs, xml['Fghz'])
        sigma = stf.extract(bufs, xml['Refcal_Sigma'])
        timestamp = Time(stf.extract(bufs, xml['Timestamp']), format='lv')
        tbg = Time(stf.extract(bufs, xml['T_beg']), format='lv')
        ted = Time(stf.extract(bufs, xml['T_end']), format='lv')
        pha = np.angle(refcal)
        amp = np.absolute(refcal)
        return {'pha': pha, 'amp': amp, 'flag': flag, 'fghz': fghz, 'sigma': sigma, 'timestamp': timestamp, 't_bg': tbg, 't_ed': ted}


def sql2phacalX(trange, *args, **kwargs):
    '''Supply a timestamp in Time format, return the closest phacal data.
        If a time range is provided, return records within the time range.'''
    import cal_header as ch
    import stateframe as stf
    xml, bufs = ch.read_calX(9, t=trange, *args, **kwargs)
    if isinstance(bufs, list):
        phacals = []
        for i, buf in enumerate(bufs):
            try:
                phacal_flag = stf.extract(buf, xml['Phacal_Flag'])
                fghz = stf.extract(buf, xml['Fghz'])
                sigma = stf.extract(buf, xml['Phacal_Sigma'])
                timestamp = Time(stf.extract(buf, xml['Timestamp']), format='lv')
                tbg = Time(stf.extract(buf, xml['T_beg']), format='lv')
                ted = Time(stf.extract(buf, xml['T_end']), format='lv')
                pha = stf.extract(buf, xml['Phacal_Pha'])
                amp = stf.extract(buf, xml['Phacal_Amp'])
                tmp = stf.extract(buf, xml['MBD'])
                poff, pslope = tmp[:, :, 0], tmp[:, :, 1]
                flag = stf.extract(buf, xml['Flag'])[:, :, 0]
                t_ref = Time(stf.extract(buf, xml['T_refcal']), format='lv')
                phacals.append({'pslope': pslope, 't_pha': timestamp, 'flag': flag, 'poff': poff, 't_ref': t_ref,
                                'phacal': {'pha': pha, 'amp': amp, 'flag': phacal_flag, 'fghz': fghz, 'sigma': sigma, 'timestamp': timestamp,
                                           't_bg': tbg, 't_ed': ted}})
            except:
                print 'failed to load record {} ---> {}'.format(i + 1, Time(stf.extract(buf, xml['Timestamp']), format='lv').iso)
        return phacals
    elif isinstance(bufs, str):
        phacal_flag = stf.extract(bufs, xml['Phacal_Flag'])
        fghz = stf.extract(bufs, xml['Fghz'])
        sigma = stf.extract(bufs, xml['Phacal_Sigma'])
        timestamp = Time(stf.extract(bufs, xml['Timestamp']), format='lv')
        tbg = Time(stf.extract(bufs, xml['T_beg']), format='lv')
        ted = Time(stf.extract(bufs, xml['T_end']), format='lv')
        pha = stf.extract(bufs, xml['Phacal_Pha'])
        amp = stf.extract(bufs, xml['Phacal_Amp'])
        tmp = stf.extract(bufs, xml['MBD'])
        poff, pslope = tmp[:, :, 0], tmp[:, :, 1]
        flag = stf.extract(bufs, xml['Flag'])[:, :, 0]
        t_ref = Time(stf.extract(bufs, xml['T_refcal']), format='lv')
        return {'pslope': pslope, 't_pha': timestamp, 'flag': flag, 'poff': poff, 't_ref': t_ref,
                'phacal': {'pha': pha, 'amp': amp, 'flag': phacal_flag, 'fghz': fghz, 'sigma': sigma, 'timestamp': timestamp, 't_bg': tbg,
                           't_ed': ted}}


def fit_blerror(out):
    ''' Determines baseline errors on baselines with Ant 14, by fitting the
        phase variations due to baseline error dependences on HA and Dec.
        (Currently only fits Bx and By errors, Bz to be added later)
        
        Required inputs:
        out   a list of calibrator observation results returned by rd_refcal(),
                corrected for feed rotation by unrot_refcal().
        
        Output:
        dbx, dby   Arrays of errors, in cm, for all baselines with Ant 14, on
                     both polarizations and in all bands (1-34).  Size is 
                     [nant, npol, nband] = [13, 2, 34]
    '''
    from scipy.optimize import curve_fit

    def bxyfunc(ha, poff, dbx, dby):
        # ha: hour angle
        # poff: constant phase offset
        # dbx: baseline x error [ns] * cos(dec) * fghz
        # dby: baseline y error [ns] * cos(dec) * fghz
        return (2. * np.pi) * (dbx * np.cos(ha) - dby * np.sin(ha)) + poff

    gdbands = np.array([4, 5, 6, 7, 8, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25])
    halist = []  # List used later for Bz
    declist = []  # List used later for Bz
    srcs = np.unique(out['srclist'])
    # First fit for Bx and By errors for all sources with longer than 4-hour HA range 
    for src in srcs:
        idxlist, = np.where(np.array(out['srclist']) == src)
        ha = []
        ph = []
        dec = out['decs'][idxlist[0]]
        declist.append(dec)
        fghz = out['fghzs'][idxlist[0]]
        for i in idxlist:
            ha.append(out['has'][i])
            ph.append(np.angle(out['vis'][i]))
        ha = np.concatenate(ha)
        halist.append(ha)
        print 'Result for source', src, ':'
        if ha[-1] - ha[0] < np.pi / 3.:
            print '***HA range is too short to fit for dBx, dBy'
            print '***Must be at least 4 hours.  Will skip this source.'
        else:
            ph = np.concatenate(ph, 3)
            nant, npol, nf, nt = ph.shape
            dbx = np.zeros((13, 2, nf), np.float)
            dby = np.zeros((13, 2, nf), np.float)
            for a in range(13):
                for pol in range(2):
                    for f in gdbands:
                        popt, pcov = curve_fit(bxyfunc, ha, np.unwrap(ph[a, pol, f]), p0=[0., 0., 0.], sigma=None, absolute_sigma=False)
                        dbx[a, pol, f] = popt[1] / (np.cos(dec) * fghz[f])  # Convert to ns
                        dby[a, pol, f] = popt[2] / (np.cos(dec) * fghz[f])  # Convert to ns
            dBx = np.median(np.mean(dbx[:, :, gdbands], 1), 1)
            xstd = np.std(np.mean(dbx[:, :, gdbands], 1), 1)
            dBy = np.median(np.mean(dby[:, :, gdbands], 1), 1)
            ystd = np.std(np.mean(dby[:, :, gdbands], 1), 1)

            print '          dBx [m]       dBy [m]'
            for i in range(13):
                print 'Ant {:2d}'.format(i + 1), '{:6.3f}+/-{:6.4f} {:6.3f}+/-{:6.4f}'.format(dBx[0] - dBx[i], xstd[i], dBy[0] - dBy[i], ystd[i])
            print 'Ant 14', '{:6.3f}+/-{:6.4f} {:6.3f}+/-{:6.4f}'.format(dBx[0], xstd[0], dBy[0], ystd[0])
        print ' '
        # Now fit for Bz errors (makes assumption that Bx and By errors are already small
        # nhamin = len(halist[0])
        # minsrc = 0
        # for i,src in enumerate(srcs):
        # idxlist, = np.where(np.array(out['srclist']) == src)
        # # Find which source has the smalled number of measurements,
        # # and find the indexes in all sources nearest to that source's HAs
        # if nhamin > len(halist[i]):
        # nhamin = len(halist[i])
        # minsrc = i
        # for i,src in enumerate(srcs):
        # idxs.append(nearest_val_idx(halist[minsrc],halist[k])
    return dbx, dby
