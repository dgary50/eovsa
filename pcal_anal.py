#
# History:
#  2017-Mar-05 YC
#    Initially written by Yi Chai
#  2017-Mar-30 DG
#    Extensive rewrites to streamline code and fix some bugs
#  2017-Jul-11 DG
#    Changed time plot frequency to middle frequency of range, and label 
#    frequency on plot
#  2018-Jan-02 DG
#    Added mv_pcal_files().  Also fixed a bug in graph(), to avoid
#    a crash when a bad/short IDB file is analyzed. 
#  2018-Jun-09 DG
#    Change plot to plot phases of antennas other than ant1 relative to ant1
#  2019-May-04 DG
#    Added mv_ptg_files().
#  2019-Jun-21 DG
#    Fixed a bug in filenames when running on pipeline.
#  2020-Jan-26 DG
#    mv_pcal_files() only worked for years 201*, so now works for 20*
#  2020-May-29 SY
#    update calpntanal() to use util.get_idbdir() to find IDB root path.
#  2021-Aug-02 DG
#    Fix a bug in findfile() introduced when FDB files were lost from the DPP.
#    Now on pipeline the IFDB files are used, which do not have a key ST_SEC.
#  2022-Mar-08 DG
#    Temporarily commented out check for Windscram due to loss of SQL
#

import numpy as np
from util import Time, lobe, fname2mjd,get_idbdir
ten_minutes = 600./86400.
one_minute = 60./86400.

def mv_pcal_files():
    ''' Moves (renames) files in the /common/webplots/phasecal folder
        into new folders according to date.  Leaves the last 20 .npz
        and associated files in the main folder.
    '''
    import glob, os
    from time import sleep
    npzfiles = glob.glob('/common/webplots/phasecal/20????????*')
    npzfiles.sort()
    #datstr = ''
    if len(npzfiles) > 20:
        for file in npzfiles[:-20]:
            #datstr_prev = datstr
            datstr = file[26:34]
            #if datstr != datstr_prev:
            directory = file[:34]+'/'
            if not os.path.exists(directory):
                #print 'mkdir',directory 
                os.makedirs(directory)
                sleep(0.1)
            #print 'mv',file,directory+os.path.basename(file) 
            os.rename(file,directory+os.path.basename(file))
            files = glob.glob('/common/webplots/phasecal/pc?'+datstr+'*')
            files.sort()    
            for f in files:
                #print 'mv',f,directory+os.path.basename(f) 
                os.rename(f,directory+os.path.basename(f))

def mv_ptg_files():
    ''' Moves (renames) files in the /common/webplots/PTG folder
        into new folders according to date.  Leaves the last 20 PTG
        files in the main folder.
    '''
    import glob, os
    from time import sleep
    files = glob.glob('/common/webplots/PTG/P*.png')
    files.sort()
    #datstr = ''
    if len(files) > 20:
        for file in files[:-20]:
            #datstr_prev = datstr
            datstr = file[26:30]
            #if datstr != datstr_prev:
            directory = file[:21]+file[24:30]+'/'
            if not os.path.exists(directory):
                #print 'mkdir',directory 
                os.makedirs(directory)
                sleep(0.1)
            #print 'mv',file,directory+os.path.basename(file) 
            os.rename(file,directory+os.path.basename(file))

def findfile(trange):

    from util import nearest_val_idx
    import struct, time, glob, sys, socket
    import dump_tsys

    host = socket.gethostname()
    if host == 'dpp':
        fpath = '/data1/IDB/'
    else:
        fpath = get_idbdir(trange[0])
    t1 = str(trange[0].mjd)
    t2 = str(trange[1].mjd)
    tnow = Time.now()

    if t1[:5] != t2[:5]:
        # End day is different than start day, so read and concatenate two fdb files
        fdb = {}
        fdb1 = dump_tsys.rd_fdb(trange[0])
        fdb2 = dump_tsys.rd_fdb(trange[1])
        for key in fdb1.keys():
            fdb.update({key:np.append(fdb1[key],fdb2[key])})
    else:
        # Both start and end times are on the same day
        fdb = dump_tsys.rd_fdb(trange[0])

    scanidx, = np.where(fdb['PROJECTID'] == 'PHASECAL')
    scans,sidx = np.unique(fdb['SCANID'][scanidx],return_index=True)
    eidx = np.append(sidx[1:],len(scanidx)) - 1
    # List of PHASECAL scan start times
    tslist = Time(fdb['ST_TS'][scanidx[sidx]].astype(float).astype(int),format='lv')
    # List of PHASECAL scan end times
    telist = Time(fdb['EN_TS'][scanidx[eidx]].astype(float).astype(int),format='lv')
    # Remove any bad values (i.e. those with ST_SEC = 0)
    try:
        good, = np.where(fdb['ST_SEC'][scanidx[sidx]] != '0')
        tslist = tslist[good]
        telist = telist[good]
    except KeyError:
        # Key 'ST_SEC' not found so just continue (this happens on pipeline when IFDB file is used)
        pass
        
    k = 0         # Number of scans within timerange
    m = 0         # Pointer to first scan within timerange
    flist = []
    status = []
    tstlist = []
    for i in range(len(tslist)):
        if tslist[i].jd >= trange[0].jd and telist[i].jd <= trange[1].jd:
            # Time is in range, so add it
            k += 1
        else:
            # Time is too early, so skip it
            m += 1
        
    if k == 0: 
        print 'No phase calibration data within given time range'
        return None
    else: 
        print 'Found',k,'scans in timerange.'
        for i in range(k):
            f1 = fdb['FILE'][np.where(fdb['SCANID'] == scans[m+i])].astype('str')
            # if fpath == '/data1/eovsa/fits/IDB/':
            #     f2 = [fpath + f[3:11] + '/' + f for f in f1]
            # else:
            #     f2 = [fpath + f for f in f1]
            if host == 'dpp':
                f2 = [fpath + f for f in f1]
            else:
                f2 = [fpath + f[3:11] + '/' + f for f in f1]
            flist.append(f2)
            tstlist.append(tslist[m+i])
            ted = telist[m+i]
            # Mark all files done except possibly the last
            fstatus = ['done']*len(f1)
            # Check if last file end time is less than 10 min ago
            if (tnow.jd - ted.jd) < (600./86400):
                # Current time is less than 10 min after this scan
                fstatus[-1] = 'undone'
            status.append(fstatus)

    return {'scanlist':flist,'status':status,'tstlist':tstlist}

def graph(f,navg=None,path=None):

    import matplotlib.pyplot as plt
    from matplotlib.ticker import FormatStrFormatter
    import struct, time, glob, sys, socket
    import read_idb as ri
    import dbutil as db

    if navg is None:
        navg = 60

    if path is None:
        path = ''

    out = ri.read_idb(f,navg=navg)
    if len(out['fghz']) == 0:
        # This file is no good, so skip it
        return
    fig, ax = plt.subplots(4,13,sharex=True, sharey=True)
    trange = Time([fname2mjd(f[0]),fname2mjd(f[-1]) + ten_minutes],format='mjd')
    # ************ This block commented out due to loss of SQL **************
    # times, wscram, avgwind = db.a14_wscram(trange)
    # nwind = len(wscram)
    # nbad = np.sum(wscram)
    nbad = 0    # Skip Windscram check
    if nbad != 0:
        warn = ' --> Windscram! ('+str(nbad)+' of '+str(nwind)+')'
        color = '#d62728'   # Plot points with "warning" Red color
    else:
        warn = ''
        color = '#1f77b4'   # Plot points with "normal" Blue color
    fig.set_size_inches(18,6)
    nf = len(out['fghz'])
    fstr = str(out['fghz'][nf/2]*1000)[:5]+' MHz '
    for k in range(13):
        for j in range(4):
            ax[j,k].cla()
            ax[j,k].plot(out['ha'],np.angle(out['x'][ri.bl2ord[k,13],j,nf/2]),'.',color=color)
            ax[j,k].set_ylim(-4, 4)
            ax[j,k].xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
            if k in range(1,13): ax[j,k].yaxis.set_visible(False)
            if j in range(3): ax[j,k].xaxis.set_visible(False)
            if j == 0: ax[0,k].title.set_text('antenna %d' %(k+1))
    fig.suptitle(out['source']+'  '+Time(out['time'][0],format='jd').iso[:19]+' UT  '+fstr+warn)
    ax[0,0].set_ylabel('XX Phase')
    ax[1,0].set_ylabel('YY Phase')
    ax[2,0].set_ylabel('XY Phase')
    ax[3,0].set_ylabel('YX Phase')
    fig.text(0.5, 0.04, 'Hour Angle', ha = 'center')
    t = Time(out['time'][0],format='jd').iso[:19].replace('-','').replace(':','').replace(' ','')
    s = out['source']
    ofile = path + t[:14] +'_'+ s +'.npz'
    np.savez(open(ofile,'wb'),out = out)
    plt.savefig(path + 'pcT'+t+'_'+ s +'.png',bbox_inches='tight')
    plt.close(fig)

    ph = np.angle(np.sum(out['x'],3))
    fig, ax = plt.subplots(4,13)
    fig.set_size_inches(18,6)
    for k in range(13):
        for j in range(4):
            ax[j,k].cla()
            if k == 0:
                ax[j,k].plot(out['fghz'],ph[ri.bl2ord[k,13],j],'.',color=color)
            else:
                ax[j,k].plot(out['fghz'],lobe(ph[ri.bl2ord[k,13],j]-ph[ri.bl2ord[0,13],j]),'.',color=color)                
            ax[j,k].set_ylim(-4, 4)
            ax[j,k].xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
            if k in range(1,13): ax[j,k].yaxis.set_visible(False)
            if j in range(3): ax[j,k].xaxis.set_visible(False)
            if j == 0: ax[0,k].title.set_text('antenna %d' %(k+1))
    fig.suptitle(out['source']+' '+Time(out['time'][0],format='jd').iso[:19]+' UT'+warn)
    ax[0,0].set_ylabel('XX Phase')
    ax[1,0].set_ylabel('YY Phase')
    ax[2,0].set_ylabel('XY Phase')
    ax[3,0].set_ylabel('YX Phase')
    fig.text(0.5, 0.04, 'Frequency[GHz]', ha = 'center')
    t = Time(out['time'][0],format='jd').iso[:19].replace('-','').replace(':','').replace(' ','')
    plt.savefig(path + 'pcF'+t+'_'+ s +'.png',bbox_inches='tight')
    plt.close(fig)

    
def pcal_anal(trange,path=None):

    import os
    import os.path
    import socket
    import glob

    if path is None:
        path = ''

    out = findfile(trange)
    if out is None:
        return
    host = socket.gethostname()
    filelist = out['scanlist']
    statuslist = out['status']
    starttimelist = out['tstlist']
    nscans = len(filelist)
    print 'Found',nscans,'scans to process.'
    for i in range(nscans):
        good, = np.where(np.array(statuslist[i]) == 'done')
        flist = np.array(filelist[i])[good].tolist()   # List of "done" files
        first_file = filelist[i][0]
        last_file = filelist[i][-1]
        mjd = fname2mjd(last_file)
        tdif = Time.now().mjd - mjd
        if len(good) == len(filelist[i]) and tdif > ten_minutes:
            # All files in this scan are marked "done", so process the scan only if the plots do not already exist
            tmark = fname2mjd(first_file)
            tmarkp = tmark+one_minute
            tmarkn = tmark-one_minute
            tmark = Time(tmark,format='mjd').iso.replace('-','').replace(':','').replace(' ','')[:12]
            tmarkp = Time(tmarkp,format='mjd').iso.replace('-','').replace(':','').replace(' ','')[:12]
            tmarkn = Time(tmarkn,format='mjd').iso.replace('-','').replace(':','').replace(' ','')[:12]
            f1 = glob.glob(path + 'pcT*'+tmark+'*.png')
            f2 = glob.glob(path + 'pcT*'+tmarkp+'*.png')
            f3 = glob.glob(path + 'pcT*'+tmarkn+'*.png')
            if f1 == [] and f2 == [] and f3 == []:
                #print 'No files:',tmarkn,tmark,tmarkp,'found.'
                print 'Processing completed scan',i+1
                graph(filelist[i],path=path)
            else:
                print 'Scan processing already complete.  Skipping scan',i+1
        elif len(good) == len(filelist[i]) and tdif < ten_minutes:
            # All files in this scan are marked "done", but it has been less than 10 min, so process the scan
            print 'Processing completed scan',i+1
            graph(flist,path=path)        
        elif len(good) == len(filelist[i])-1:
            # This scan is still active, so process all files up to this point.
            print 'Processing active scan',i+1
            if flist != []:
                graph(flist,path=path)
            else:
                print 'No files to process (yet).'

            
if __name__ == '__main__':
    import matplotlib
    matplotlib.use('Agg')
    path='/common/webplots/phasecal/'
    t1 = Time.now().jd-0.25
    t2 = Time.now().jd
    trange = Time([t1,t2],format='jd')
    print trange.iso
    pcal_anal(trange,path=path)
