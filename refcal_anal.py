''' Routine to read and analyze REFCAL observations to generate REFCAL calibration SQL database.
'''
#
#  2017-05-13 BC
#    Began writing the code
import read_idb as ri
from util import Time
import numpy as np
import matplotlib.pyplot as plt
from glob import glob
import matplotlib.dates as mdates
import os
import pcapture2 as p
import pdb
import chan_util_bc as cu

bl2ord = p.bl_list()


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
    t1 = str(trange[0].mjd)
    t2 = str(trange[1].mjd)
    daydelta = int(trange[1].mjd - trange[0].mjd)
    tnow = Time.now()
    if t1[:5] != t2[:5]:
        # End day is different than start day, so read and concatenate two fdb files
        ufdb1 = dump_tsys.rd_ufdb(trange[0])
        ufdb = ufdb1
        for ll in xrange(daydelta):
            ufdb2 = dump_tsys.rd_ufdb(Time(trange[0].mjd + ll, format='mjd'))
            if ufdb2:
                for key in ufdb.keys():
                    ufdb.update({key: np.append(ufdb[key], ufdb2[key])})
    else:
        # Both start and end times are on the same day
        ufdb = dump_tsys.rd_ufdb(trange[0])

    if srcid:
        if type(srcid) is str:
            srcid=[srcid]
        sidx_=np.array([])
        for sid in srcid:
            sidx, = np.where((ufdb['PROJECTID'] == projid) & (ufdb['SOURCEID']  == sid))
            sidx_=np.append(sidx_,sidx)
        scanidx=np.sort(sidx_).astype('int')
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
        if tslist[i].jd >= trange[0].jd and telist[i].jd <= trange[1].jd:
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


def rd_refcal(trange, projid='PHASECAL', srcid=None, quackint=180.,navg=3):
    '''take a time range from the Time object, e.g., trange=Time(['2017-04-08T05:00','2017-04-08T15:30']),
       a projectid, and source id, return visibility data for all baselines correlated with ant 14.
       ***Optional keywords***
       projid: string -- predefined PROJECTID when setting up the observations. Default is 'PHASECAL'
       srcid: string -- if provided, then only use the specified source. E.g., '1229+020' is often used for
              reference calibration
       quackint: interval in seconds to skip at the beginning of each scan
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
    for n, scan in enumerate(scanlist):
        print 'Reading scan: '+scan
        out = ri.read_idb([scan], navg=navg, quackint=quackint)
        times.append(out['time'])
        nt = len(out['time'])
        if nt == 0:
            continue 
        bds, sidx = np.unique(out['band'], return_index=True)
        nbd = len(bds)
        eidx = np.append(sidx[1:], len(out['band']))
        vs = np.zeros((15, 2, 34, nt), dtype=complex)
        fghz = np.zeros(34)
        # average over channels within each band
        for b, bd in enumerate(bds):
            fghz[bd-1]=np.nanmean(out['fghz'][sidx[b]:eidx[b]])
            for a in range(13):
                for pl in range(2):
                    vs[a, pl, bd - 1] = np.mean(out['x'][bl2ord[a, 13], pl, sidx[b]:eidx[b]], axis=0)
        vis.append(vs)
        fghzs.append(fghz)
        bandnames.append(bds)
    return {'scanlist': scanlist, 'srclist': srclist, 'tstlist': sclist['tstlist'], 'tedlist': sclist['tedlist'],
            'vis': vis, 'bandnames': bandnames, 'fghzs':fghzs, 'times': times}

def graph(out, refcal=None, ant_str='ant1-13', bandplt=[5, 11, 17, 23], scanidx=None, pol=0,
          tformat='%H:%M'):
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
    ant_list = p.ant_str2list(ant_str)
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
                        ax1[ant, b].text(-0.4, 0.5, 'Ant ' + str(ant + 1), ha='center', va='center',
                                         transform=ax1[ant, b].transAxes, fontsize=10)
                        ax2[ant, b].text(-0.4, 0.5, 'Ant ' + str(ant + 1), ha='center', va='center',
                                         transform=ax2[ant, b].transAxes, fontsize=10)
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
                    if 't_bg'  in refcal.keys():
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

def refcal_anal(out, timerange=None, scanidx=None, minsnr=0.7, bandplt=[5, 11, 17, 23]):
    '''Analyze the visibility data from rd_refcal and return time averaged visibility values and flags.
       ***Optional Keywords***
       timerange: time range to obtain the average. E.g., timerange=Time(['2017-04-08T05:00','2017-04-08T07:00'])
       scanidx: index numbers of scans to select. Useful when other (undesired) types of observations exist. 
               !!!! The selected scans should have exactly the same frequency/band setup !!!!
       minsnr: minimum signal to noise to consider. Data with smaller SNRs will be flagged (as 1 in the flag array)
       bandplt: bands to show in the figure
       ***Outputs***
       refcal: complex array of shape (15, 2, 34) (nant, npol, nband) as the result of the reference calibration
       flag: int array of shape (15, 2, 34). 0 is unflagged and 1 is flagged.
       timestamp: midpoint of the time range used for averaging to obtain the refcal values 
    '''
    if scanidx:
        scanlist = [out['scanlist'][i] for i in scanidx]
        tstlist = [out['tstlist'][i] for i in scanidx]
        tedlist = [out['tedlist'][i] for i in scanidx]
        bandnames = [out['bandnames'][i] for i in scanidx]
        vis = [out['vis'][i] for i in scanidx]
        times = [out['times'][i] for i in scanidx]
        fghzs = [out['fghzs'][i] for i in scanidx]
    else:
        scanlist = out['scanlist']
        tstlist = out['tstlist']
        tedlist = out['tedlist']
        bandnames = out['bandnames']
        vis = out['vis']
        times = out['times']
        fghzs = out['fghzs']

    fghz=fghzs[0]
    times = np.concatenate(times)
    vis = np.concatenate(vis, axis=3)
    if timerange:
        tidx, = np.where((times > timerange[0].jd) & (times < timerange[1].jd))
        vis = vis[:, :, :, tidx]
        timeavg = times[tidx]
    else:
        timeavg = times
    # vismean = np.nanmean(np.angle(vis),axis=3)
    vis_ = np.zeros(vis.shape[:3], dtype=complex)
    flag = np.zeros(vis.shape[:3], dtype=int)
    sigma = np.zeros(vis.shape[:3])+1e10
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
    visavg = {'pha': np.angle(vis_), 'amp': np.abs(vis_), 'timestamp': timestamp, 
              't_bg':timeavg[0], 't_ed':timeavg[-1], 'flag': flag}
    graph(out, visavg, scanidx=scanidx, bandplt=bandplt)
    f2, ax2 = plt.subplots(2, 13, figsize=(12, 5))
    f3, ax3 = plt.subplots(2, 13, figsize=(12, 5))
    allbands = np.arange(34) + 1
    for ant in range(13):
        for pol in range(2):
            ind, = np.where(flag[ant, pol, :] == 0)
            ax2[pol, ant].plot(allbands[ind], np.unwrap(visavg['pha'][ant, pol, ind]), '.', markersize=5)
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
    return {'vis': vis_, 'pha':np.angle(vis_), 'amp':np.abs(vis_), 'fghz':fghz, 'flag': flag, 
            'sigma':sigma, 'timestamp': timestamp, 't_gcal': timestamp_gcal,
            't_bg': Time(timeavg[0], format='jd'), 't_ed': Time(timeavg[-1], format='jd')}

def graph_results(refcal,unwrap=True):
    '''Provide an output from sql2refcal() (single refcal result) or sql2refcalX() (list of refcal results).
       Plot phase and amp vs. bands'''
    f1, ax1 = plt.subplots(2, 13, figsize=(12, 5))
    f2, ax2 = plt.subplots(2, 13, figsize=(12, 5))
    allbands = np.arange(34) + 1
    for ant in range(13):
        for pol in range(2):
            if type(refcal) is dict:
                refcal=[refcal]
            for ref in refcal:
                ind, = np.where(ref['flag'][ant, pol, :] == 0)
                if unwrap:
                    pha=np.unwrap(ref['pha'][ant, pol, ind])
                    ax1[pol, ant].set_ylim([-20, 20])
                else:
                    pha=ref['pha'][ant, pol, ind]
                    ax1[pol, ant].set_ylim([-4, 4])
                ax1[pol, ant].plot(allbands[ind], pha, '.', markersize=5)
                ax1[pol, ant].set_xlim([1, 34])
                ax2[pol, ant].plot(allbands[ind], ref['amp'][ant, pol, ind], '.', markersize=5)
                ax2[pol, ant].set_xlim([1, 34])
                ax2[pol, ant].set_ylim([0, 1.])
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

def sql2refcal(t):
    '''Supply a timestamp in Time format, return the closest refcal data'''
    import cal_header as ch
    import stateframe as stf
    xml, buf = ch.read_cal(8, t=t)
    refcal = stf.extract(buf, xml['Refcal_Real']) + stf.extract(buf, xml['Refcal_Imag']) * 1j
    flag = stf.extract(buf, xml['Refcal_Flag'])
    fghz = stf.extract(buf, xml['Fghz'])
    sigma = stf.extract(buf, xml['Refcal_Sigma'])
    timestamp = Time(stf.extract(buf, xml['Timestamp']), format='lv')
    pha = np.angle(refcal)
    amp = np.absolute(refcal)
    return {'pha': pha, 'amp': amp, 'flag':flag, 'fghz':fghz, 'sigma':sigma, 'timestamp': timestamp}

def sql2refcalX(trange):
    import eovsapy.cal_header as ch
    import stateframe as stf
    xml, bufs = ch.read_calX(8,t=trange)
    refcals=[]
    for i,buf in enumerate(bufs):
        try:
            ref = stf.extract(buf, xml['Refcal_Real']) + stf.extract(buf, xml['Refcal_Imag']) * 1j
            flag = stf.extract(buf, xml['Refcal_Flag'])
            timestamp = Time(stf.extract(buf, xml['Timestamp']), format='lv')
            pha = np.angle(ref)
            amp = np.absolute(ref)
            refcals.append({'pha': pha, 'amp': amp, 'flag':flag, 'timestamp': timestamp})
        except:
            continue
    return refcals

