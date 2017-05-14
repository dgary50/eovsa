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

# identify REFCAL files
def findfiles(trange, projid='PHASECAL', srcid=None):
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
    # takes input from findfiles (sclist)
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
        out = ri.read_idb([scan],quackint=quackint)
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

def graph(out, visavg=None, ant_str='ant1-13', scanidx=None,pol=0):
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
    bds0=[4,6,8,10,12,14]
    nband = len(bds0)
    f, ax = plt.subplots(nant,nband,figsize=(12,8))
    for n, scan in enumerate(scanlist):
        ts = Time(times[n],format='jd').plot_date
        #if n == 0:
        #    bds0 = bands[0] 
        #    nband = len(bds0)

        # make a line plot of phase vs. time on all 34 bands
        for ant in ant_list:
            for b,bd in enumerate(bds0):
                ax[ant,b].plot_date(ts,np.degrees(np.angle(vis[n][bd-1,ant,pol])),'.',markersize=1)
                ax[ant,b].xaxis.set_major_formatter(date_format)
                if visavg:
                    phavg=np.degrees(np.angle(visavg['vis'][bd-1,ant,pol]))
                    bt=Time(visavg['time'][0],format='jd').plot_date
                    et=Time(visavg['time'][-1],format='jd').plot_date
                    ax[ant,b].plot_date([bt,et],[phavg,phavg],'-r')
                    ax[ant,b].plot_date([bt,bt],[-200,200],'--b')
                    ax[ant,b].plot_date([et,et],[-200,200],'--b')
                ax[ant,b].set_ylim([-180.,180.])
                if ant == 0:
                    ax[ant,b].set_title('Band '+str(bd))
                if b == 0:
                    ax[ant,b].text(-0.6,0.5,'Ant '+str(ant+1),ha='center',va='center',transform=ax[ant,b].transAxes,fontsize=10)
                else:
                    ax[ant,b].set_yticks([])
    return f,ax

def refcal_anal(out,timerange=None,scanidx=None):
    times = np.concatenate(out['times'])
    vis = np.concatenate(out['vis'],axis=3)
    if timerange:
        tidx, = np.where((times > timerange[0].jd) & (times < timerange[1].jd))
        vis = vis[:,:,:,tidx]
        timeavg = times[tidx]
    vis_ = np.mean(vis,axis=3)
    flag = np.zeros(vis_.shape)
    visavg = {'pha':np.angle(vis_),'amp':np.abs(vis_), 'vis':vis_, 'time':timeavg, 'flag':flag}
    f1, ax1 = graph(out,visavg,scanidx=scanidx)
    f2, ax2 = plt.subplots(2,13,figsize=(12,5)) 
    for ant in range(13):
        for pol in range(2):
            ax2[pol,ant].plot(np.arange(34)+1,np.unwrap(visavg['pha'][:,ant,pol]),'.',markersize=5)
            ax2[pol,ant].set_ylim([-20,20])
            if ant == 0:
                ax2[pol,ant].set_ylabel('Phase (radian)')
            else:
                ax2[pol,ant].set_yticks([])
            if pol == 1 and ant ==0:
                ax2[pol,ant].set_xlabel('Band #')
            if pol == 0:
                ax2[pol,ant].set_title('Ant '+str(ant+1))
                ax2[pol,ant].set_xticks([])

    return visavg


        

        

            
             

                
