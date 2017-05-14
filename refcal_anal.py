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

def rd_refcal(sclist):
    # takes input from findfiles (sclist)
    import struct, time, glob, sys, socket
    import dbutil as db
    scanlist = sclist['scanlist']
    srclist = sclist['srclist']
    starttimelist = sclist['tstlist']

    # read scans one by one
    nscan=len(scanlist)
    vis = []
    bands = []
    times = []
    for n, scan in enumerate(scanlist):
        out = ri.read_idb([scan],quackint=180.)
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

def graph(out,ant_str='ant1-13',pol=0):
    # takes input from rd_refcal (out)
    date_format = mdates.DateFormatter('%H')
    # make a color plot showing the phase
    # ri.summary_plot_pcal(out)
    scanlist=out['scanlist']
    nscan=len(scanlist)
    bands=out['bands']
    vis=out['vis']
    ant_list = p.ant_str2list(ant_str)
    nant = len(ant_list)
    for n, scan in enumerate(scanlist):
        ts = Time(out['times'][n],format='jd').plot_date
        if n == 0:
            bds0 = bands[0] 
            nband = len(bds0)
            f, ax = plt.subplots(nant,nband,figsize=(12,8))

        # make a line plot of phase vs. time on all 34 bands
        for ant in ant_list:
            for b,bd in enumerate(bds0):
                ax[ant,b].plot_date(ts,np.degrees(np.angle(vis[n][bd-1,ant,pol])),'.',markersize=1)
                ax[ant,b].xaxis.set_major_formatter(date_format)
                ax[ant,b].set_ylim([-180.,180.])
                if ant == 0:
                    ax[ant,b].set_title('Band '+str(bd))
                if b == 0:
                    ax[ant,b].text(-1.,0.5,'Ant '+str(ant+1),ha='center',va='center',transform=ax[ant,b].transAxes,fontsize=10)
                else:
                    ax[ant,b].set_yticks([])

def refcal_anal(out,timerange=None):
    times = out['times']
    if timerange:
        print 'do something'
        

            
             

                
