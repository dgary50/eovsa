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
    '''identify refcal files'''
    from util import nearest_val_idx
    import struct, time, glob, sys, socket
    import dump_tsys
    fpath = '/data1/eovsa/fits/UDB/'+trange[0].iso[:4]+'/'
    t1 = str(trange[0].mjd)
    t2 = str(trange[1].mjd)
    tnow = Time.now()
    if t1[:5] != t2[:5]:
        # End day is different than start day, so read and concatenate two fdb files
        ufdb = {}
        ufdb1 = dump_tsys.rd_ufdb(trange[0])
        ufdb2 = dump_tsys.rd_ufdb(trange[1])
        for key in ufdb1.keys():
            ufdb.update({key:np.append(ufdb1[key],ufdb2[key])})
    else:
        # Both start and end times are on the same day
        ufdb = dump_tsys.rd_ufdb(trange[0])

    if srcid:
        scanidx, = np.where((ufdb['PROJECTID'] == projid) & (ufdb['SOURCEID'] == srcid)) 
    else:
        scanidx, = np.where(ufdb['PROJECTID'] == projid)
    # List of scan start times
    tslist = Time(ufdb['ST_TS'][scanidx].astype(float).astype(int),format='lv')
    # List of PHASECAL scan end times
    telist = Time(ufdb['EN_TS'][scanidx].astype(float).astype(int),format='lv')

    k = 0         # Number of scans within timerange
    m = 0         # Pointer to first scan within timerange
    flist = []
    status = []
    tstlist = []
    srclist = []
    for i in range(len(tslist)):
        if tslist[i].jd >= trange[0].jd and telist[i].jd <= trange[1].jd:
            flist.append(fpath+ufdb['FILE'][scanidx[i]].astype('str'))
            tstlist.append(tslist[i])
            srclist.append(ufdb['SOURCEID'][scanidx[i]])
            k += 1
        
    if k == 0: 
        print 'No scans found within given time range for '+projid
        return None
    else: 
        print 'Found',k,'scans in timerange.'

    return {'scanlist':flist,'srclist':srclist,'tstlist':tstlist}

def rd_refcal(trange, projid='PHASECAL', srcid=None, quackint=180.):
    '''take a time range from the Time object, e.g., trange=Time(['2017-04-08T05:00','2017-04-08T15:30']),
       a projectid, and source id, return visibility data for all baselines correlated with ant 14.
       **Optional keywords***
       projid: string -- predefined PROJECTID when setting up the observations. Default is 'PHASECAL'
       srcid: string -- if provided, then only use the specified source. E.g., '1229+020' is often used for
              reference calibration
       quackint: interval in seconds to skip at the beginning of each scan
    '''
    import struct, time, glob, sys, socket
    import dbutil as db
    sclist=findfiles(trange,projid,srcid)
    scanlist = sclist['scanlist']
    srclist = sclist['srclist']
    starttimelist = sclist['tstlist']

    # read scans one by one
    nscan=len(scanlist)
    vis = []
    bands = []
    times = []
    for n, scan in enumerate(scanlist):
        out = ri.read_idb([scan],navg=3,quackint=quackint)
        times.append(out['time'])
        nt = len(out['time'])
        bds,sidx=np.unique(out['band'],return_index=True)
        nbd = len(bds)
        eidx = np.append(sidx[1:],len(out['band'])) 
        vs = np.zeros((34,15,2,nt),dtype=complex)
        # average over channels within each band
        for b,bd in enumerate(bds):
            for a in range(13): 
                for pl in range(2):
                    vs[bd-1,a,pl] = np.mean(out['x'][bl2ord[a,13],pl,sidx[b]:eidx[b]],axis=0)
        vis.append(vs)
        bands.append(bds)
    return {'scanlist':scanlist,'srclist':srclist,'tstlist':starttimelist,'vis':vis, 'bands':bands, 'times':times} 

def graph(out, visavg=None, ant_str='ant1-13', bandplt=[4,10,16,22], scanidx=None, pol=0):
    '''Produce a figure showing phases and amplitudes of selected antennas, bands, and polarization.
       Optionally takes in a time averaged value 'visavg' to show the selected time range and plot the
       averaged phase and amplitudes'''
    # takes input from rd_refcal (out)
    date_format = mdates.DateFormatter('%H')
    # make a color plot showing the phase
    # ri.summary_plot_pcal(out)
    if scanidx:
        scanlist=[out['scanlist'][i] for i in scanidx]
        bands=[out['bands'][i] for i in scanidx]
        vis=[out['vis'][i] for i in scanidx]
        times=[out['times'][i] for i in scanidx]
    else:
        scanlist=out['scanlist']
        bands=out['bands']
        vis=out['vis']
        times=out['times']

    nscan=len(scanlist)
    ant_list = p.ant_str2list(ant_str)
    nant = len(ant_list)
    bds0=bandplt
    nband = len(bds0)
    f1, ax1 = plt.subplots(nant,nband,figsize=(12,8))
    f2, ax2 = plt.subplots(nant,nband,figsize=(12,8))
    for n, scan in enumerate(scanlist):
        ts = Time(times[n],format='jd').plot_date

        # make plots of phase and amp vs. time on all 34 bands
        for ant in ant_list:
            for b,bd in enumerate(bds0):
                ph_deg=np.angle(vis[n][bd-1,ant,pol],deg=True)
                amp=np.abs(vis[n][bd-1,ant,pol])
                ax1[ant,b].plot_date(ts,ph_deg,'.',markersize=5)
                ax1[ant,b].xaxis.set_major_formatter(date_format)
                ax2[ant,b].plot_date(ts,amp,'.',markersize=5)
                ax2[ant,b].xaxis.set_major_formatter(date_format)
                if visavg:
                    phavg=np.angle(visavg['vis'][bd-1,ant,pol],deg=True)
                    ampavg=np.abs(visavg['vis'][bd-1,ant,pol])
                    flg=visavg['flag'][bd-1,ant,pol]
                    if flg == 0:
                        bt=Time(visavg['time'][0],format='jd').plot_date
                        et=Time(visavg['time'][-1],format='jd').plot_date
                        ax1[ant,b].plot_date([bt,et],[phavg,phavg],'-r')
                        ax1[ant,b].plot_date([bt,bt],[-200,200],'--b')
                        ax1[ant,b].plot_date([et,et],[-200,200],'--b')
                        ax2[ant,b].plot_date([bt,et],[ampavg,ampavg],'-r')
                        ax2[ant,b].plot_date([bt,bt],[0,np.max(amp)],'--k')
                        ax2[ant,b].plot_date([et,et],[0,np.max(amp)],'--k')
                ax1[ant,b].set_ylim([-180.,180.])
                if ant == 0:
                    ax1[ant,b].set_title('Band '+str(bd))
                    ax2[ant,b].set_title('Band '+str(bd))
                if b == 0:
                    ax1[ant,b].text(-0.4,0.5,'Ant '+str(ant+1),ha='center',va='center',transform=ax1[ant,b].transAxes,fontsize=10)
                    ax2[ant,b].text(-0.4,0.5,'Ant '+str(ant+1),ha='center',va='center',transform=ax2[ant,b].transAxes,fontsize=10)
                else:
                    ax1[ant,b].set_yticks([])
                    ax2[ant,b].set_yticks([])

def refcal_anal(out, timerange=None, scanidx=None, minsnr=0.7, bandplt=[4,10,16,22]):
    '''Analyze the visibility data from rd_refcal and return time averaged visibility values and flags.
       ***Optional Keywords***
       timerange: time range to obtain the average. E.g., timerange=Time(['2017-04-08T05:00','2017-04-08T07:00'])
       scanidx: index numbers of scans to select. Useful when other (undesired) types of observations exist.
       minsnr: minimum signal to noise to consider. Data with smaller SNRs will be flagged (as 1 in the flag array)
       bandplt: bands to show in the figure
    '''
    times = np.concatenate(out['times'])
    vis = np.concatenate(out['vis'],axis=3)
    if timerange:
        tidx, = np.where((times > timerange[0].jd) & (times < timerange[1].jd))
        vis = vis[:,:,:,tidx]
        timeavg = times[tidx]
    else:
        timeavg = times
    #vismean = np.nanmean(np.angle(vis),axis=3)
    vis_ = np.zeros(vis.shape[:3],dtype=complex)
    flag = np.zeros(vis.shape[:3],dtype=int)
    # compute standard deviation of the phases
    sigma = np.nanstd(vis,axis=3)
    #sigma = np.nanstd(np.abs(vis),axis=3)
    # mask out records with phases > 3 sigma
    for bd in range(34):
        #print 'band: ',bd
        for ant in range(13):
            #print 'ant: ',ant
            for pol in range(2):
                #print 'pol: ',pol
                #ph=np.unwrap(np.angle(vis[bd,ant,pol]))
                #ph=np.angle(vis[bd,ant,pol])
                #ph_mean=np.angle(np.nanmean(vis[bd,ant,pol]))
                amp=np.abs(vis[bd,ant,pol])
                amp_median=np.median(np.abs(vis[bd,ant,pol]))
                #ind,=np.where(np.abs(np.angle(vis[bd,ant,pol])-np.angle(vismean[bd,ant,pol])) % (2.*np.pi) < sigma_ph[bd,ant,pol])
                #ind,=np.where(np.abs(ph-ph_mean) < np.radians(100.))
                ind, = np.where(np.abs(amp-amp_median) < sigma[bd,ant,pol])
                snr = amp_median/sigma[bd,ant,pol]
                #pdb.set_trace()
                if snr < minsnr or np.isnan(snr):
                    flag[bd,ant,pol]=1
                else:
                    if len(ind) > len(timeavg)/2:
                        vis_[bd,ant,pol]=np.nanmean(vis[bd,ant,pol,ind])
                    else:
                        vis_[bd,ant,pol]=np.nanmean(vis[bd,ant,pol])
                        flag[bd,ant,pol]=1
                #print '# of valid datapoints: ',len(ind)
                #pdb.set_trace()
                #if len(ind) > len(timeavg)/2: #at least half of the data points are good
                #    vis_[bd,ant,pol]=np.nanmean(vis[bd,ant,pol,ind])
                #else:
                #    flag[bd,ant,pol]=1
        # count how many datapoints are flagged in a given band
        nflag = np.count_nonzero(flag[bd,:13])
        print '{0:d} of 26 measurements are flagged due to SNR < {1:.1f} in Band {2:d}'.format(nflag,minsnr,bd+1)
    #flag zero values (e.g., antennas or bands not observed)
    zeroind = np.where(np.abs(vis_ == 0))
    flag[zeroind]=1
    visavg = {'pha':np.angle(vis_),'amp':np.abs(vis_), 'vis':vis_, 'time':timeavg, 'flag':flag}
    graph(out,visavg,scanidx=scanidx,bandplt=bandplt)
    f2, ax2 = plt.subplots(2,13,figsize=(12,5)) 
    allbands=np.arange(34)+1
    for ant in range(13):
        for pol in range(2):
            ind,=np.where(flag[:,ant,pol] == 0)
            ax2[pol,ant].plot(allbands[ind],np.unwrap(visavg['pha'][ind,ant,pol]),'.',markersize=5)
            ax2[pol,ant].set_ylim([-20,20])
            ax2[pol,ant].set_xlim([1,34])
            if ant == 0:
                ax2[pol,ant].set_ylabel('Phase (radian)')
            else:
                ax2[pol,ant].set_yticks([])
            if pol == 1 and ant ==0:
                ax2[pol,ant].set_xlabel('Band #')
            if pol == 0:
                ax2[pol,ant].set_title('Ant '+str(ant+1))
                ax2[pol,ant].set_xticks([])
    refcal=vis_
    timestamp=Time(np.mean(timeavg),format='jd')
    return refcal, flag, timestamp


        

        

            
             

                
