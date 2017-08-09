''' Main routine to read the auto- and cross-correlation data from IDB files.
'''
#
#  2016-03-30  DG
#    Began writing the code, adapted from get_X_data2.py
#  2016-04-01  DG
#    Added summary_plot()
#  2016-05-12  DG
#    Truncate arrays in readXdata() when the file ends early
#  2016-05-14  DG
#    Several changes.  Changed filter keyword in readXdata() 
#    to default to False, since we are missing a lot of 
#    frequencies.  Instead any filtering will be done in 
#    read_idb() after arrays are concatenated.  Also added
#    time-averaging in read_idb().  Also added flag_sk(),
#    which eliminates bad data due to RFI (sort of)
#  2016-06-23 BC
#    Changed the default search directory of the EOVSA IDB files
#    in get_trange_files() from '/dppdata1/IDB' to /data1/eovsa/fits/IDB'
#  2016-06-30  DG
#    Change to allow input of a file list to read_idb() in place of
#    trange (used by the realtime pipeline).
#  2016-07-23  DG
#    Widen range of allowed SK to 0.7 - 1.5.
#  2016-11-19  DG
#    Fix glitches in uvw when averaging over time, by replacing
#    missing values with nan before averaging.  Also added srcchk boolean
#    to read_idb(), to skip the name check if False.  Also added the
#    option to specify a source name, and changed behavior so that
#    if a different source name is found, it just skips the file
#    instead of bailing out.
#  2016-12-28  DG
#    Fixed error in wrapping of HA.  Added new "production" routine
#    unrot_miriad(), which does the correction of raw IDB files for
#    differential feed rotation and writes out a new file.  This
#    routine probably needs to move to another module, but it is
#    still undergoing debugging.
#  2017-Jan-05  BC, DG
#    Added "reverse" parameter to unrot() so that simulation can go
#    either unrotating (reverse = False, default), or forward rotating
#    (reverse = True).  Also remove hard-coded 500 in xsampler and 
#    ysampler reading, now that IDB data no longer have non-existent
#    frequencies.
#  2017-Jan-11  DG
#    Added effect of multi-band delay to both unrot_miriad() [Jones
#    matrix approach] and unrot_miriad2() [Mueller matrix approach].
#  2017-Jan-19  DG
#    Change unrot_miriad() to use difference in feed angle on each
#    baseline rather than the angles themselves.  Also implements a
#    multi-band delay on each antenna.
#  2017-Jan-25  DG
#    Further changes to unrot_miriad(), and removed unrot_miriad2().
#    Also added code to handle extraneous nulls ('\x00') at end of
#    string vars written by aipy.
#  2017-Jan-27  DG
#    Fixed some bugs associated with tp_only.
#  2017-Feb-01  DG
#    Added read_udb() routine.
#  2017-Apr-04  DG
#    Added read_npz() routine to read and concatenate multiple NPZ files.
#  2017-Apr-14  DG
#    Very occasional 0-filled record in IDB file was throwing 
#    off summed times when navg was set in read_idb().  Now readXdata()
#    detects this and skips that record (a 1-s data gap).
#  2017-Apr-27  DG
#    Added allday_udb() routine, to read all UDB files for a given
#    day and optionally make a nice overview plot of the data.
#  2017-May-02  DG
#    Added saving of allday_udb() figure if requested by savfig keyword.
#  2017-May-13  BC
#    Added summary_plot_pcal() to plot only baselines correlating with Ant 14
#    Added "quackint" parameter in read_idb() to get rid of first quackint seconds of each file 
#    Truncate out['uvw'] to have the same shape of the time axis as out['time']
#  2017-May-14  BC
#    Added a key (out['band']) in the output of readXdata() to indicate the band name of a specific frequency
#    in fghz
#  2017-Jun-12  DG
#    Fix allday_udb() to find files on following date up to 09 UT
#  2017-Jul-12
#    Fixed long-standing bug in get_trange_files(), which no longer worked on DPP.
#  2017-Jul-13  DG
#    Major update to allday_udb() to overplot GOES data (if available) and scan type info.
#  2017-Jul-15  DG
#    Simpler calculation for band names from frequency, instead of calling freq2bdname(),
#    which anyway does not work...
#

import aipy
import os
from util import Time, nearest_val_idx
import glob
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.dates import DateFormatter
import spectrogram_fit as sp
import pcapture2 as p
import eovsa_lst as el
import copy
import chan_util_bc as cu

bl2ord = p.bl_list()

def read_udb(filename):
    ''' This routine reads the data from a UDB file.
    '''

    # Open uv file for reading
    uv = aipy.miriad.UV(filename)
    nf = len(uv['sfreq'])
    nt = uv['ntimes']

    freq = uv['sfreq']
    npol = uv['npol']
    nants = uv['nants']
    nbl = nants*(nants-1)/2
    outa = np.zeros((nants,npol,nf,nt),dtype=np.complex64)  # Auto-correlations
    outx = np.zeros((nbl,npol,nf,nt),dtype=np.complex64)  # Cross-correlations
    #outp = np.zeros((nants,2,nf,600),dtype=np.float)
    #outp2 = np.zeros((nants,2,nf,600),dtype=np.float)
    #outm = np.zeros((nants,2,nf,600),dtype=np.int)
    uvwarray = np.zeros((nbl,nt,3),dtype=np.float)
    timearray = []
    lstarray = []
    l = -1
    tprev = 0
    tsav = 0
    # Use antennalist if available
    ants = uv['antlist']
    while ants[-1] == '\x00': ants = ants[:-1]
    antlist = map(int, ants.split())

    src = uv['source']
    while src[-1] == '\x00': src = src[:-1]

    for preamble, data in uv.all():
        uvw, t, (i0,j0) = preamble
        i = antlist.index(i0+1)
        j = antlist.index(j0+1)
        if i > j:
            # Reverse order of indices
            j = antlist.index(i0+1)
            i = antlist.index(j0+1)
        # Assumes uv['pol'] is one of -5, -6, -7, -8
        k = -5 - uv['pol']
        if t != tprev:
            # New time 
            l += 1
            if l == nt:
                break
            tprev = t
            timearray.append(t)
            #xdata = uv['xsampler'].reshape(nf_orig,nants,3)
            #ydata = uv['ysampler'].reshape(nf_orig,nants,3)
            #outp[:,0,:,l] = np.swapaxes(xdata[:,:,0],0,1)
            #outp[:,1,:,l] = np.swapaxes(ydata[:,:,0],0,1)
            #outp2[:,0,:,l] = np.swapaxes(xdata[:,:,1],0,1)
            #outp2[:,1,:,l] = np.swapaxes(ydata[:,:,1],0,1)
            #outm[:,0,:,l] = np.swapaxes(xdata[:,:,2],0,1)
            #outm[:,1,:,l] = np.swapaxes(ydata[:,:,2],0,1)

        if i0 == j0:
            # This is an auto-correlation
            outa[i0,k,:,l] = data
            if k < 2 and np.sum(data != np.real(data)) > 0:
                print preamble,uv['pol'], 'has imaginary data!'
        else:
            outx[bl2ord[i,j],k,:,l] = data
            if k == 3: uvwarray[bl2ord[i,j],l] = uvw

    # Truncate in case of early end of data
    #nt = len(timearray)
    #outp = outp[:,:,:,:nt]
    #outp2 = outp2[:,:,:,:nt]
    #outm = outm[:,:,:,:nt]
    #if not tp_only:
    #    outa = outa[:,:,:,:nt]
    #    outx = outx[:,:,:,:nt]

    if len(lstarray) != 0:
        pass
    else:
        tarray = Time(timearray,format='jd')
        for t in tarray:
            lstarray.append(el.eovsa_lst(t))
    ha = np.array(lstarray) - uv['ra']
    ha[np.where(ha > np.pi)] -= 2*np.pi
    ha[np.where(ha < -np.pi)] += 2*np.pi
    out = {'a':outa, 'x':outx, 'uvw':uvwarray, 'fghz':freq, 'time':np.array(timearray),'source':src,'ha':ha,'ra':uv['ra'],'dec':uv['dec']}#,'p':outp,'p2':outp2,'m':outm
    return out
    
def readXdata(filename, filter=False, tp_only=False, src=None):
    ''' This routine reads the data from a single IDBfile.
        
        Optiona Keywords:
        filter   boolean--if True, returns only non-zero frequencies 
                    if False (default), returns uniform set of 500 frequencies
        tp_only  boolean--if True, returns only TP information
                    if False (default), returns everything (including 
                    auto & cross correlations)
    '''

    # Open uv file for reading
    uv = aipy.miriad.UV(filename)
    nf_orig = len(uv['sfreq'])
    good_idx = np.arange(nf_orig)
    if filter:
        good_idx = []
        # Read a bunch of records to get number of good frequencies, i.e. those with at least 
        # some non-zero data.  Read 20 records for baseline 1-2, XX pol
        uv.select('antennae',0,2,include=True)
        uv.select('polarization',-5,-5,include=True)
        for i in range(20):
            preamble, data = uv.read()
            idx, = data.nonzero()
            if len(idx) > len(good_idx):
                good_idx = copy.copy(idx)
        uv.select('clear',0,0)
        uv.rewind()

    if 'source' in uv.vartable:
        source = uv['source']
        while source[-1] == '\x00': source = source[:-1]
        if src is None:
            # If no source name is given, return the source from the file and keep going
            src = source
        elif src != source:
            # If a specific source name is given, and it does not match the file, stop and return None
            return source
        else:
            # If a source is given, and it matches the file, keep going
            pass
    else:
        if src:
            # If a specific source name is given, and there is no source in the file, stop and return
            pass#return '<no "source" var!>'
    nf = len(good_idx)
    freq = uv['sfreq'][good_idx]
    npol = uv['npol']
    nants = uv['nants']
    nbl = nants*(nants-1)/2
    if not tp_only:
        outa = np.zeros((nants,npol,nf,600),dtype=np.complex64)  # Auto-correlations
        outx = np.zeros((nbl,npol,nf,600),dtype=np.complex64)  # Cross-correlations
    outp = np.zeros((nants,2,nf,600),dtype=np.float)
    outp2 = np.zeros((nants,2,nf,600),dtype=np.float)
    outm = np.zeros((nants,2,nf,600),dtype=np.int)
    uvwarray = np.zeros((nbl,600,3),dtype=np.float)
    timearray = []
    lstarray = []
    l = -1
    tprev = 0
    tsav = 0
    # Use antennalist if available
    if 'antlist' in uv.vartable:
        ants = uv['antlist']
        while ants[-1] == '\x00': ants = ants[:-1]
        antlist = map(int, ants.split())
    else:
        antlist = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]

    for preamble, data in uv.all():
        uvw, t, (i0,j0) = preamble
        i = antlist.index(i0+1)
        j = antlist.index(j0+1)
        if i > j:
            # Reverse order of indices
            j = antlist.index(i0+1)
            i = antlist.index(j0+1)
        # Assumes uv['pol'] is one of -5, -6, -7, -8
        k = -5 - uv['pol']
        if filter:
            if len(data.nonzero()[0]) == nf:
                if t != tprev:
#                    print preamble
                    if t == 2440587.5:
                        # Time is 1970-01-01, which means a zero-filled record, so skip
                        # the entire thing.
                        continue
                    # New time 
                    l += 1
                    if l == 600:
                        break
                    tprev = t
                    timearray.append(t)
                    try:
                        lstarray.append(uv['lst'])
                    except:
                        pass
                    xdata = uv['xsampler'].reshape(nf_orig,nants,3)
                    ydata = uv['ysampler'].reshape(nf_orig,nants,3)
                    outp[:,0,:,l] = np.swapaxes(xdata[good_idx,:,0],0,1)
                    outp[:,1,:,l] = np.swapaxes(ydata[good_idx,:,0],0,1)
                    outp2[:,0,:,l] = np.swapaxes(xdata[good_idx,:,1],0,1)
                    outp2[:,1,:,l] = np.swapaxes(ydata[good_idx,:,1],0,1)
                    outm[:,0,:,l] = np.swapaxes(xdata[good_idx,:,2],0,1)
                    outm[:,1,:,l] = np.swapaxes(ydata[good_idx,:,2],0,1)
    
                if tp_only:
                    outa = None
                    outx = None
                else:
                    if i0 == j0:
                        # This is an auto-correlation
                        outa[i0,k,:,l] = data[data.nonzero()]
                    else:
                        outx[bl2ord[i,j],k,:,l] = data[data.nonzero()]
                        if k == 3: uvwarray[bl2ord[i,j],l] = uvw[data.nonzero()]
        else:
            if t != tprev:
                # New time 
                if t == 2440587.5:
                    # Time is 1970-01-01, which means a zero-filled record, so skip
                    # the entire thing.
                    continue
                l += 1
                if l == 600:
                    break
                tprev = t
                timearray.append(t)
                xdata = uv['xsampler'].reshape(nf_orig,nants,3)
                ydata = uv['ysampler'].reshape(nf_orig,nants,3)
                outp[:,0,:,l] = np.swapaxes(xdata[:,:,0],0,1)
                outp[:,1,:,l] = np.swapaxes(ydata[:,:,0],0,1)
                outp2[:,0,:,l] = np.swapaxes(xdata[:,:,1],0,1)
                outp2[:,1,:,l] = np.swapaxes(ydata[:,:,1],0,1)
                outm[:,0,:,l] = np.swapaxes(xdata[:,:,2],0,1)
                outm[:,1,:,l] = np.swapaxes(ydata[:,:,2],0,1)

            if tp_only:
                outa = None
                outx = None
            else:
                if i0 == j0:
                    # This is an auto-correlation
                    outa[i0,k,:,l] = data
                    if k < 2 and np.sum(data != np.real(data)) > 0:
                        print preamble,uv['pol'], 'has imaginary data!'
                else:
                    outx[bl2ord[i,j],k,:,l] = data
                    if k == 3: uvwarray[bl2ord[i,j],l] = uvw

    # Truncate in case of early end of data
    nt = len(timearray)
    outp = outp[:,:,:,:nt]
    outp2 = outp2[:,:,:,:nt]
    outm = outm[:,:,:,:nt]
    uvwarray = uvwarray[:,:nt]
    if not tp_only:
        outa = outa[:,:,:,:nt]
        outx = outx[:,:,:,:nt]

    if len(lstarray) != 0:
        pass
    else:
        tarray = Time(timearray,format='jd')
        for t in tarray:
            lstarray.append(el.eovsa_lst(t))
    ha = np.array(lstarray) - uv['ra']
    ha[np.where(ha > np.pi)] -= 2*np.pi
    ha[np.where(ha < -np.pi)] += 2*np.pi
    # Find out band name for each frequency
    bd = (freq*2 - 1).astype(np.int)
    #bd=[]
    #for f in freq:
    #    bd.append(cu.freq2bdname(f))
    out = {'a':outa, 'x':outx, 'uvw':uvwarray, 'fghz':freq, 'band':np.array(bd),'time':np.array(timearray),'source':src,'p':outp,'p2':outp2,'m':outm,'ha':ha,'ra':uv['ra'],'dec':uv['dec']}
    return out

def readXdatmp(filename):
    # This temporary routine reads the data from a single IDBfile where the
    # polarization was written incorrectly. (-1,-2,0,0) instead of (-5,-6,-7,-8)
    
    # the IDB array sorts basline pairs into correct order for an output
    # array that just has 18 baselines (2 sets of 4 antennas correlated separately
    # now just need to convert i,j into slot in 136-slot array: 15*16/2 corrs +16 auto
    # this array corresponds to the first index (auto corr) for each antenna
    # 16*(iant-1)-(iant-1)(iant-2)/2
    ibl = np.array( [ 0,16,31,45,58,70,81,91,100,108,115,121,126,130,133,135 ])

    # Open uv file for reading
    uv = aipy.miriad.UV(filename)
    good_idx = []
    # Read a bunch of records to get number of good frequencies, i.e. those with at least 
    # some non-zero data.  Read 20 records for baseline 1-2, XX pol
    uv.select('antennae',0,2,include=True)
    uv.select('polarization',-1,-1,include=True)
    for i in range(20):
        preamble, data = uv.read()
        idx, = data.nonzero()
        if len(idx) > len(good_idx):
            good_idx = copy.copy(idx)
    if 'source' in uv.vartable:
        src = uv['source']
    uv.select('clear',0,0)
    uv.rewind()
    nf = len(good_idx)
    freq = uv['sfreq'][good_idx]
    npol = uv['npol']
    nants = uv['nants']
    nbl = nants*(nants-1)/2
    outa = np.zeros((nants,npol,nf,600),dtype=np.complex64)  # Auto-correlations
    outx = np.zeros((nbl,npol,nf,600),dtype=np.complex64)  # Cross-correlations
    outp = np.zeros((nants,2,nf,600),dtype=np.float)
    outp2 = np.zeros((nants,2,nf,600),dtype=np.float)
    outm = np.zeros((nants,2,nf,600),dtype=np.int)
    uvwarray = []
    timearray = []
    l = -1
    tprev = 0
    # Use antennalist if available
    if 'antlist' in uv.vartable:
        ants = uv['antlist']
        antlist = map(int, ants.split())
    else:
        antlist = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]
    kp = 0
    for preamble, data in uv.all():
        uvw, t, (i0,j0) = preamble
        i = antlist.index(i0+1)
        j = antlist.index(j0+1)
        if i > j:
            # Reverse order of indices
            j = antlist.index(i0+1)
            i = antlist.index(j0+1)
        # Assumes uv['pol'] is -1, -2, 0, 0
        if i == 0 and j == 0:
            if uv['pol'] == -1:
                k = 0
            elif uv['pol'] == -2:
                k = 1
            elif uv['pol'] == 0:
                if kp == 1: 
                    k = 2
                    kp = 0
                else:
                    k = 3
                    kp = 1
        if k == 3:
            uvwarray.append(uvw)
        if len(data.nonzero()[0]) == nf:
            if t != tprev:
                # New time 
                l += 1
                if l == 600:
                    break
                tprev = t
                timearray.append(t)
                xdata = uv['xsampler'].reshape(500,nants,3)
                ydata = uv['ysampler'].reshape(500,nants,3)
                outp[:,0,:,l] = np.swapaxes(xdata[good_idx,:,0],0,1)
                outp[:,1,:,l] = np.swapaxes(ydata[good_idx,:,0],0,1)
                outp2[:,0,:,l] = np.swapaxes(xdata[good_idx,:,1],0,1)
                outp2[:,1,:,l] = np.swapaxes(ydata[good_idx,:,1],0,1)
                outm[:,0,:,l] = np.swapaxes(xdata[good_idx,:,2],0,1)
                outm[:,1,:,l] = np.swapaxes(ydata[good_idx,:,2],0,1)

            if i0 == j0:
                # This is an auto-correlation
                outa[i0,k,:,l] = data[data.nonzero()]
            else:
                outx[bl2ord[i,j],k,:,l] = data[data.nonzero()]
    out = {'a':outa, 'x':outx, 'uvw':np.array(uvwarray), 'fghz':freq, 'time':np.array(timearray),'source':src,'p':outp,'p2':outp2,'m':outm}
    return out
    
def summary_plot(out,ant_str='ant1-13',ptype='phase',pol='XX-YY'):
    ''' Makes a summary amplitude or phase plot for all baselines from ants in ant_str
        in out dictionary.
    '''
    import matplotlib.pyplot as plt
    
    ant_list = p.ant_str2list(ant_str)
    nant = len(ant_list)
    if ptype != 'amp' and ptype != 'phase':
        print "Invalid plot type.  Must be 'amp' or 'phase'."
        return
    poloff = 0
    if pol != 'XX-YY':
        poloff = 2
    f, ax = plt.subplots(nant,nant)
    f.subplots_adjust(hspace=0,wspace=0)
    for axrow in ax:
        for a in axrow:
            a.xaxis.set_visible(False)
            a.yaxis.set_visible(False)
    for i in range(nant-1):
        ai = ant_list[i]
        for j in range(i+1,nant):
            aj = ant_list[j]
            if ptype == 'phase':
                ax[i,j].imshow(np.angle(out['x'][bl2ord[ai,aj],0+poloff]))
                ax[j,i].imshow(np.angle(out['x'][bl2ord[ai,aj],1+poloff]))
            elif ptype == 'amp':
                ax[i,j].imshow(np.abs(out['x'][bl2ord[ai,aj],0+poloff]))
                ax[j,i].imshow(np.abs(out['x'][bl2ord[ai,aj],1+poloff]))
    for i in range(nant):
        ai = ant_list[i]
        ax[i,i].text(0.5,0.5,str(ai+1),ha='center',va='center',transform=ax[i,i].transAxes,fontsize=14)

def summary_plot_pcal(out,ant_str='ant1-14',ptype='phase',pol='XX-YY'):
    ''' Makes a summary amplitude or phase plot for all baselines from ants in ant_str
        in out dictionary.
    '''
    import matplotlib.pyplot as plt
    ha=out['ha']
    fghz=out['fghz']
    ant_list = p.ant_str2list(ant_str)
    nant = len(ant_list)
    if ptype != 'amp' and ptype != 'phase':
        print "Invalid plot type.  Must be 'amp' or 'phase'."
        return
    poloff = 0
    if pol != 'XX-YY':
        poloff = 2
    f, ax = plt.subplots(nant-1,2,figsize=(5,8))
    #f.subplots_adjust(hspace=0,wspace=0)
    for axrow in ax:
        for a in axrow:
            a.xaxis.set_visible(False)
            a.yaxis.set_visible(False)
    for i in range(nant-1):
        ai = ant_list[i]
        if ptype == 'phase':
            #ax[i,0].pcolormesh(ha,fghz,np.angle(out['x'][bl2ord[ai,13],0+poloff]))
            #ax[i,1].pcolormesh(ha,fghz,np.angle(out['x'][bl2ord[ai,13],1+poloff]))
            ax[i,0].imshow(np.angle(out['x'][bl2ord[ai,13],0+poloff]))
            ax[i,1].imshow(np.angle(out['x'][bl2ord[ai,13],1+poloff]))
        elif ptype == 'amp':
            #ax[i,0].pcolormesh(ha,fghz,np.abs(out['x'][bl2ord[ai,13],0+poloff]))
            #ax[i,1].pcolormesh(ha,fghz,np.abs(out['x'][bl2ord[ai,13],1+poloff]))
            ax[i,0].imshow(np.abs(out['x'][bl2ord[ai,13],0+poloff]))
            ax[i,1].imshow(np.abs(out['x'][bl2ord[ai,13],1+poloff]))
        ax[i,0].text(-0.1,0.5,str(ai+1),ha='center',va='center',transform=ax[i,0].transAxes,fontsize=14)
        polstr=pol.split('-')
        for j in range(2):
            ax[0,j].text(0.5,1.3,polstr[j],ha='center',va='center',transform=ax[0,j].transAxes,fontsize=14)

def allday_udb(t=None, doplot=True, goes_plot=True, savfig=False, gain_corr=False):
    # Plots (and returns) UDB data for an entire day
    from sunpy import lightcurve
    from sunpy.time import TimeRange
    from flare_monitor import flare_monitor
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
    out = read_idb(files,src='Sun')
    if gain_corr:
        import gaincal2 as gc
        out = gc.apply_gain_corr(out)
    trange = Time(out['time'][[0,-1]], format = 'jd')
    fghz = out['fghz']
    if doplot:
        f, ax = plt.subplots(1,1,figsize=(14,5))
        pdata = np.sum(np.sum(np.abs(out['x'][0:11,:]),1),0)  # Spectrogram to plot
        X = np.sort(pdata.flatten())   # Sorted, flattened array
        # Set any time gaps to nan
        tdif = out['time'][1:] - out['time'][:-1]
        bad, = np.where(tdif > 120./86400)  # Time gaps > 2 minutes
        pdata[:,bad] = 0
        vmax = X[int(len(X)*0.95)]  # Clip at 5% of points
        im = ax.pcolormesh(Time(out['time'],format='jd').plot_date,out['fghz'],pdata,vmax=vmax)
        plt.colorbar(im,ax=ax,label='Amplitude [arb. units]')
        ax.xaxis_date()
        ax.xaxis.set_major_formatter(DateFormatter("%H:%M"))
        ax.set_ylim(fghz[0], fghz[-1])
        ax.set_xlabel('Time [UT]')
        ax.set_ylabel('Frequency [GHz]')
        ax.set_title('EOVSA 1-min Data for '+t.iso[:10])

        if goes_plot:
            # Get GOES data for overplotting
            goes_tr = TimeRange(trange.iso)
            goes_label = ['A','B','C','M','X']
            # The GOES label is placed to start 20 min into the day
            goes_label_time = Time(out['time'][[0]], format = 'jd').plot_date + 0.014
            rightaxis_label_time = trange[1].plot_date

            # Retrieve GOES data for the day, but this only goes to end of UT day
            try:
                goes = lightcurve.GOESLightCurve.create(goes_tr)
                goes.data['xrsb'] = 2* (np.log10(goes.data['xrsb'])) + 26
                ytext = np.median(goes.data['xrsb']) - 1
                ax.text (goes_label_time, ytext, 'GOES soft x-ray data', color = 'yellow')
                goes.data['xrsb'].plot(color = 'yellow')
            except:
                # Looks like the GOES data do not exist, so just skip it
                pass
            for k,i in enumerate([10,12,14,16,18]):
                ax.text(rightaxis_label_time, i-0.4, goes_label[k], fontsize = '12')
                ax.plot_date(rightaxis_label_time + np.array([-0.005,0.0]),[i,i],'-',color='yellow')
            try:
                # If the day goes past 0 UT, get GOES data for the next UT day
                if int(trange[1].mjd) != int(trange[0].mjd):
                    goes_tr2 = TimeRange([trange[1].iso[:10], trange[1].iso])
                    goesday2 = lightcurve.GOESLightCurve.create(goes_tr2)
                    goesday2.data['xrsb'] = 2* (np.log10(goesday2.data['xrsb'])) + 26
                    goesday2.data['xrsb'].plot(color = 'yellow')
            except:
                # Looks like the GOES data do not exist, so just skip it
                pass
        ax.set_xlim(trange.plot_date)

        ut, fl, projdict = flare_monitor(t)
        if fl == []:
            print 'Error retrieving data for',t.iso[:10],'from SQL database.'
            return
        if projdict == {}:
            print 'No annotation can be added to plot for',t.iso[:10]
        else:
            nscans = len(projdict['Project'])
            SOS = Time(projdict['Timestamp'],format='lv').plot_date
            EOS = Time(projdict['EOS'],format='lv').plot_date
            yran = np.array(ax.get_ylim())
            for i in range(nscans):
                uti = SOS[i]*np.array([1.,1.])
                #if uti[0] >= trange[0].plot_date:
                ax.plot_date(uti,yran,'g',lw=0.5)
                if projdict['Project'][i] == 'NormalObserving' or projdict['Project'][i] == 'Normal Observing':
                    ax.text(uti[0],yran[1]*0.935,'SUN',fontsize=8, color = 'white')
                elif projdict['Project'][i] == 'None':
                    ax.text(uti[0],yran[1]*0.975,'IDLE',fontsize=8, color = 'white')
                elif projdict['Project'][i][:4] == 'GAIN':
                    ax.text(uti[0],yran[1]*0.955,'GCAL',fontsize=8, color = 'white')
                elif projdict['Project'][i] == 'SOLPNTCAL':
                    ax.text(uti[0],yran[1]*0.955,'TPCAL',fontsize=8, color = 'white')
                elif projdict['Project'][i] == 'PHASECAL':
                    ax.text(uti[0],yran[1]*0.955,'PCAL',fontsize=8, color = 'white')
                else:
                    ax.text(uti[0],yran[1]*0.975,projdict['Project'][i],fontsize=8, color = 'white')
            if len(projdict['EOS']) == nscans:
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
                plt.savefig('/common/webplots/flaremon/daily/XSP'+date+'.png',bbox_inches='tight')
    return out

# def allday_udb(t=None, doplot=True, savfig=False):
    # # Plots (and returns) UDB data for an entire day
    # if t is None:
        # t = Time.now()
    # date = t.iso[:10].replace('-','')
    # # Look also at the following day, up to 9 UT
    # date2 = Time(t.mjd + 1,format='mjd').iso[:10].replace('-','')
    # year = date[:4]
    # files = glob.glob('/data1/eovsa/fits/UDB/'+year+'/UDB'+date+'*')
    # files.sort()
    # files2 = glob.glob('/data1/eovsa/fits/UDB/'+year+'/UDB'+date2+'0*')
    # files2.sort()
    # files = np.concatenate((np.array(files),np.array(files2)))
    # # Eliminate files starting before 10 UT on date (but not on date2)
    # for i,file in enumerate(files):
        # if file[-6] != '0':
            # break
    # try:
        # files = files[i:]
    # except:
        # print 'No files found in /data1/eovsa/fits/UDB/ for',date
        # return {}
    # out = read_idb(files,src='Sun')
    # if doplot:
        # f, ax = plt.subplots(1,1)
        # f.set_size_inches(14,5)
        # pdata = np.sum(np.sum(np.abs(out['x'][0:11,:]),1),0)  # Spectrogram to plot
        # X = np.sort(pdata.flatten())   # Sorted, flattened array
        # # Set any time gaps to nan
        # tdif = out['time'][1:] - out['time'][:-1]
        # bad, = np.where(tdif > 120./86400)  # Time gaps > 2 minutes
        # pdata[:,bad] = 0
        # vmax = X[int(len(X)*0.95)]  # Clip at 5% of points
        # ax.pcolormesh(Time(out['time'],format='jd').plot_date,out['fghz'],pdata,vmax=vmax)
        # ax.xaxis_date()
        # ax.xaxis.set_major_formatter(DateFormatter("%H:%M"))
        # ax.set_xlabel('Time [UT]')
        # ax.set_ylabel('Frequency [GHz]')
        # ax.set_title('EOVSA 1-min Data for '+t.iso[:10])
        # if savfig:
            # plt.savefig('/common/webplots/flaremon/XSP_later.png',bbox_inches='tight')
    # return out
        
def get_IDBfiles(showthelast=10):
    #This will return the most recent IDB files saved to the 
    #  directory /data1/IDB/. They will be returned in a list
    #  with the format '/data1/IDB/IDByyyymmddhhmmss'
    #  We can adjust the number of files shown with 'showthelast' optional variable.
    #  It will not show the file being written currently
    s = showthelast
    IDBfiles = os.listdir('/dppdata1/IDB/')
    IDBfiles.sort()
    IDBfiles = IDBfiles[-(s+2):-2]
    for i in range(len(IDBfiles)):
        IDBfiles[i] = '/dppdata1/IDB/'+IDBfiles[i] 
    return IDBfiles
    
def read_idb(trange,navg=None,quackint=0.,filter=True,srcchk=True,src=None,tp_only=False):
    ''' This finds the IDB files within a given time range and concatenates 
        the times into a single dictionary.  If trange is not a Time() object,
        assume that it is the list of files to read.
        
        Keywords:
          src      string--if not None, files in trange will be skipped unless
                    their source name matches this string.
        Optional Keywords:
          filter   boolean--if True (default), returns only non-zero frequencies 
                    if False, returns uniform set of 500 frequencies, with gaps
          srcchk   boolean--if True (default), stops reading files when source 
                    name is different from initial source name.  If set to
                    False, only the source name in the first file is returned
          tp_only  boolean--if True, returns only TP information
                    if False (default), returns everything (including 
                    auto & cross correlations)
          quackint  float--first time range (in seconds) to skip in the beginning of
                    each file. Default is 0., or no quack.
    '''
    if type(trange) == Time:
        files = get_trange_files(trange)
    else:
        # If input type is not Time, assume that it is the list of files to read
        files = trange

    datalist = []
    for file in files:
        #This will skip any files that give us errors.
        #  The names of the bad or unreadable files will
        #  be printed.
        try:
            out = readXdata(file,tp_only=tp_only,src=src)
            if type(out) is str:
                print 'Source name:',out,'does not match requested name:',src+'.  Will skip',file
            else:
                if srcchk and src is None:
                    # This is the first file, and we care about the source, so set source name
                    src = out['source']
                if navg:
                    # Perform time average over navg seconds. Note that this does not do the
                    # right thing over time gaps (yet)        
                    # First set any time-frequency bins with M value 0 to nan
                    badidx = np.where(out['m'][0,0] == 0)
                    out['p'][:,:,badidx[0],badidx[1]] = np.nan
                    out['p2'][:,:,badidx[0],badidx[1]] = np.nan
                    if not tp_only:
                        out['a'][:,:,badidx[0],badidx[1]] = np.nan
                        out['x'][:,:,badidx[0],badidx[1]] = np.nan
                        out['uvw'][:,badidx[1],:] = np.nan
                    # Truncate arrays so that times are evenly divisible by navg
                    nt = len(out['time'])
                    nout = nt/navg
                    ntnew = nout*navg
                    out['m'] = out['m'][:,:,:,:ntnew]
                    out['p'] = out['p'][:,:,:,:ntnew]
                    out['p2'] = out['p2'][:,:,:,:ntnew]
                    if not tp_only:
                        out['a'] = out['a'][:,:,:,:ntnew]
                        out['x'] = out['x'][:,:,:,:ntnew]
                    out['time'] = out['time'][:ntnew]
                    out['ha'] = out['ha'][:ntnew]
                    out['uvw'] = out['uvw'][:,:ntnew,:]
                    # Recast shape 
                    out['m'].shape = out['m'].shape[0:3]+(nout,navg)
                    out['p'].shape = out['p'].shape[0:3]+(nout,navg)
                    out['p2'].shape = out['p2'].shape[0:3]+(nout,navg)
                    if not tp_only:
                        out['a'].shape = out['a'].shape[0:3]+(nout,navg)
                        out['x'].shape = out['x'].shape[0:3]+(nout,navg)
                    out['time'].shape = (nout,navg)
                    out['ha'].shape = (nout,navg)
                    out['uvw'].shape = (out['uvw'].shape[0],nout,navg,3)
                    # Perform the average (mean) power and add to out dictionary
                    out.update({'meanp':np.nanmean(out['p'],4)})
                    # Perform sum over m, p and p2, to preserve SK
                    out['m'] = np.nansum(out['m'],4)
                    out['p'] = np.nansum(out['p'],4)
                    out['p2'] = np.nansum(out['p2'],4)
                    if not tp_only:
                        # Perform the average (mean), over non-nan values
                        out['a'] = np.nanmean(out['a'],4)
                        out['x'] = np.nanmean(out['x'],4)
                    out['time'] = np.mean(out['time'],1)  # Weighted average time
                    out['uvw'] = np.nanmean(out['uvw'],2)
                    ha = np.mean(np.unwrap(out['ha']),1)  # Weighted average "unwrapped" ha
                    # Wrap it again...
                    ha[np.where(ha > np.pi)] -= 2*np.pi
                    ha[np.where(ha < -np.pi)] += 2*np.pi
                    out['ha'] = ha

                datalist.append(out)
        except:
            print 'The problematic file is:',file
            
    # Have to concatenate outa, outx, uvw, time, and ha arrays
    outa = []
    outx = []
    outp = []
    outp2 = []
    outm = []
    uvw = []
    time = []
    ha = []
    # This could be a lot of data, so handle P data first to determine
    # frequencies to eliminate
    for out in datalist:
        outp.append(out['p'])
    out['p'] = np.concatenate(outp,3)
    if filter:
        # Eliminate frequencies where there is no nonzero value
        # sums power over every dimension except freq.
        goodidx, = np.sum(np.sum(np.sum(out['p'],3),1),0).nonzero()
        # Eliminate frequencies where there is no nonzero value
        out['p'] = out['p'][:,:,goodidx]
    for out in datalist:
        if filter:
            # Eliminate frequencies where there is no nonzero value
            # before concatenation, to reduce memory load
            out['p2'] = out['p2'][:,:,goodidx]
            out['m'] = out['m'][:,:,goodidx]
            if not tp_only:
                out['a'] = out['a'][:,:,goodidx]
                out['x'] = out['x'][:,:,goodidx]
            out['fghz'] = out['fghz'][goodidx]
            out['band'] = out['band'][goodidx]
        if quackint > 0.:
            # time interval between data points in seconds
            dt = np.nanmedian(np.diff(out['time']))*86400.             
            nt = np.rint(quackint/dt).astype('int')
            # check if nt is too large
            if nt < len(out['time']):
                out['a']=out['a'][:,:,:,nt:]
                out['x']=out['x'][:,:,:,nt:]
                out['p2']=out['p2'][:,:,:,nt:]
                out['m']=out['m'][:,:,:,nt:]
                out['uvw']=out['uvw'][:,nt:]
                out['time']=out['time'][nt:]
                out['ha']=out['ha'][nt:]
        outa.append(out['a'])
        outx.append(out['x'])
        outp2.append(out['p2'])
        outm.append(out['m'])
        uvw.append(out['uvw'])
        time.append(out['time'])
        ha.append(out['ha'])
    if tp_only:
        out['a'] = outa
        out['x'] = outx
    else:
        out['a'] = np.concatenate(outa,3)
        out['x'] = np.concatenate(outx,3)
    out['p2'] = np.concatenate(outp2,3)
    out['m'] = np.concatenate(outm,3)
    out['uvw'] = np.concatenate(uvw,1)
    out['time'] = np.concatenate(time)
    out['ha'] = np.concatenate(ha)
    return out
    
def read_npz(files):
    ''' This reads already-processed data from Numpy compressed files
        in the given file-list, and concatenates the times into a single 
        dictionary.  The result is the same as read_idb() on several files.
        
        files     A list of npz files.
    '''
    # Have to concatenate outa, outx, uvw, time, and ha arrays
    outa = []
    outx = []
    outp = []
    outp2 = []
    outm = []
    uvw = []
    time = []
    ha = []
    for file in files:
        f = open(file,'rb')
        data = np.load(f)
        if file == files[0]:
            out = data[data.keys()[0]].item()
        outp.append(data[data.keys()[0]].item()['p'])
        outa.append(data[data.keys()[0]].item()['a'])
        outx.append(data[data.keys()[0]].item()['x'])
        outp2.append(data[data.keys()[0]].item()['p2'])
        outm.append(data[data.keys()[0]].item()['m'])
        uvw.append(data[data.keys()[0]].item()['uvw'])
        time.append(data[data.keys()[0]].item()['time'])
        ha.append(data[data.keys()[0]].item()['ha'])
        f.close()
    out['p'] = np.concatenate(outp,3)
    out['a'] = np.concatenate(outa,3)
    out['x'] = np.concatenate(outx,3)
    out['p2'] = np.concatenate(outp2,3)
    out['m'] = np.concatenate(outm,3)
    out['uvw'] = np.concatenate(uvw,1)
    out['time'] = np.concatenate(time)
    out['ha'] = np.concatenate(ha)
    return out
    
def flag_sk(out):
    m = out['m']
    sk = (m+1.)/(m-1.)*(m*out['p2']/(out['p']**2) - 1)
    u_lim = 1.5
    l_lim = 0.7
    sk_flag = np.logical_or(sk > u_lim,sk < l_lim)
    nant,npol,nf,nt = m.shape
    for i in range(nant-1):
        for j in range(i+1,nant):
            flags = np.logical_or(sk_flag[i],sk_flag[j])
            out['x'][bl2ord[i,j],flags] = np.nan
    for i in range(nant):
        if 'meanp' in out:
            out['meanp'][i,sk_flag[i]] = np.nan
        out['a'][i,sk_flag[i]] = np.nan
    return out

def fname2mjd(filename):
    fstem = filename.split('/')[-1]
    fstr = fstem[3:7]+'-'+fstem[7:9]+'-'+fstem[9:11]+' '+fstem[11:13]+':'+fstem[13:15]+':'+fstem[15:17]
    t = Time(fstr)
    return t.mjd

def get_trange_files(trange):
    #Given a timerange, this routine will take all relevant IDBfiles from
    #  that time range, put them in a list, and return that list.
    #  This function is used in get_X_data(data).
    fstr = trange[0].iso
    # look for environmental variable EOVSADB first
    datadir=os.getenv('EOVSADB')
    if not datadir:
        # go to default directory on pipeline
        datadir='/data1/eovsa/fits/IDB/'

    folder=datadir+fstr.replace('-','').split()[0]
    try:
        os.listdir(folder)
    except:
        try:
            folder = '/data1/IDB'
            os.listdir(folder)
        except:
            print 'Something wrong with the definition of EOVSA data directory.'
            print 'Best to define a EOVSADB variable in .cshrc (c-shell) or .bashrc (bash)'
            return

    files = glob.glob(folder+'/IDB'+fstr.replace('-','').split()[0]+'*')
    files.sort()
    mjd1, mjd2 = trange.mjd.astype('int')
    if mjd2 != mjd1:
        if (mjd2 - 1) != mjd1:
            usage('Second date must differ from first by at most 1 day')
        else:
            fstr2 = trange[1].iso
            files2 = glob.glob(folder+'/IDB'+fstr2.replace('-','').split()[0]+'*')
            files2.sort()
            files += files2

    def fname2mjd(filename):
        fstem = filename.split('/')[-1]
        fstr = fstem[3:7]+'-'+fstem[7:9]+'-'+fstem[9:11]+' '+fstem[11:13]+':'+fstem[13:15]+':'+fstem[15:17]
        t = Time(fstr)
        return t.mjd

    filelist = []
    for filename in files:
        mjd = fname2mjd(filename)
        if mjd >= trange[0].mjd and mjd < trange[1].mjd:
            filelist.append(filename)
    return filelist
    
def show_selfcalibrated(index, trange, plot='multipanel', antennas=0):
    #for index, choose a time when the flare starting, with a large slope. 
    #the options for plot are 'multipanel' , 'saturation' , and 'uncalibrated'
    #antennas controls which pair the saturation plot shows. 
  
    data = get_trange_files(trange)
    IDBdata, uvw, freq, times = get_X_data(data)
    
    s = sp.Spectrogram(trange)
    s.fidx = [0,IDBdata.shape[2]]
    tsys, std = s.get_median_data()
    
    pcal = np.angle(IDBdata[:,:, :, index])
    acal = abs(IDBdata[:,:, :, index])
    calout = copy.copy(IDBdata)
    # Calibrate for time 'index', just before the initial peak of the flare.
    for i in range(IDBdata.shape[3]):
        calout[:,:,:,i] = calout[:,:,:,i]*(np.cos(pcal)-1j*np.sin(pcal))/acal
    # Normalize to the total power spectrum at the same time, with reference
    # to the shortest baseline (preserves the relative amplitudes on various
    # baselines and polarizations.
    norm = abs(calout[:,:, :, index])
    for i in range(12):
        for j in range(2):
            norm[i,j,:] = tsys[:,index]*abs(calout[i,j,:,index]) / abs(calout[0,0,:,index])
    for i in range(IDBdata.shape[3]):
        calout[:,:,:,i] = calout[:,:,:,i]*norm       
    if plot == 'multipanel':
        # Multi-panel Plot
        f, ax = plt.subplots(4,5)
        sbl = ['1-2','1-3','1-4','2-3','2-4','3-4','5-7','5-8','7-8']
        for i,ibl in enumerate([0,1,2,3,4,5,7,8,11]):
            if (i > 4):
                ax[2,i % 5].imshow(abs(calout[ibl,0,50:,:]))
                ax[2,i % 5].text(100,10,sbl[i]+' Amp',color='white')
                ax[3,i % 5].imshow(np.angle(calout[ibl,0,50:,:]))
                ax[3,i % 5].text(100,10,sbl[i]+' Phase')
            else:
                ax[0,i % 5].imshow(abs(calout[ibl,0,50:,:]))
                ax[0,i % 5].text(100,10,sbl[i]+' Amp',color='white')
                ax[1,i % 5].imshow(np.angle(calout[ibl,0,50:,:]))
                ax[1,i % 5].text(100,10,sbl[i]+' Phase')
            ax[2,4].imshow(tsys[50:,:])
            ax[2,4].text(100,10,'Total Power',color='white')
            plt.subplots_adjust(left=0.02, bottom=0.03, right=0.99, top=0.97, wspace=0.20, hspace=0.20)
    else:
        if plot == 'saturation':
            # Saturation Plot
            ants_ = 'Ants ' + str(antennas+1) + '-' + str(antennas+2)
            plt.figure()
            plt.plot(tsys[IDBdata.shape[2]-110,:],abs(calout[antennas, 0, IDBdata.shape[2]-110,:]),'.',label=str(freq[IDBdata.shape[2]-110])[:5]+' GHz')
            plt.plot(tsys[IDBdata.shape[2]-60,:],abs(calout[antennas, 0, IDBdata.shape[2]-60,:]),'.',label=str(freq[IDBdata.shape[2]-60])[:5]+' GHz')
            plt.plot(tsys[IDBdata.shape[2]-10,:],abs(calout[antennas, 0, IDBdata.shape[2]-10,:]),'.',label=str(freq[IDBdata.shape[2]-10])[:5]+' GHz')
            plt.xlabel('Total Power [sfu]')
            plt.ylabel('Correlated Power (' + ants_ + ') [sfu]')
            plt.legend(loc='lower right')    
        else:
            if plot == 'uncalibrated':                          
                #"uncalibrated" 
                f, ax = plt.subplots(4,5)
                sbl = ['1-2','1-3','1-4','2-3','2-4','3-4','5-7','5-8','7-8']
                for i,ibl in enumerate([0,1,2,3,4,5,7,8,11]):
                    if (i > 4):
                        ax[2,i % 5].imshow(abs(IDBdata[ibl,0,50:,:]))
                        ax[2,i % 5].text(100,10,sbl[i]+' Amp',color='white')
                        ax[3,i % 5].imshow(np.angle(IDBdata[ibl,0,50:,:]))
                        ax[3,i % 5].text(100,10,sbl[i]+' Phase')
                    else:
                        ax[0,i % 5].imshow(abs(IDBdata[ibl,0,50:,:]))
                        ax[0,i % 5].text(100,10,sbl[i]+' Amp',color='white')
                        ax[1,i % 5].imshow(np.angle(IDBdata[ibl,0,50:,:]))
                        ax[1,i % 5].text(100,10,sbl[i]+' Phase')
                ax[2,4].imshow(tsys[50:,:])
                ax[2,4].text(100,10,'Total Power',color='white')
                plt.subplots_adjust(left=0.02, bottom=0.03, right=0.99, top=0.97, wspace=0.20, hspace=0.20)
            else:
                print 'please choose valid plot type'

def unrot(data,params=[1.,0.,0.],reverse=False):
    ''' params gives non-ideal behavior of second ant, as a, d, chi0,
        where a is gain factor Y/X, d is cross-talk as fraction of X,
        and chi0 is feed rotation in degrees.
    '''
    import feed_rot_simulation as frs
    import dbutil as db
    nant = 14
    trange = Time([data['time'][0],data['time'][-1]],format='jd')
    times, chi = db.get_chi(trange)
    nt, = times.shape
    nf, = data['fghz'].shape
    # Find only "good" chi values, i.e. those not too close
    # to zero, since missing values are zero.
    good = []
    for i in range(nant):
        good.append(np.where(np.abs(chi[:,i]) > 0.001)[0])
        if len(good[i]) == 0:
            # This antenna must be offline, so go ahead and use all times.
            # They will be ignored anyway
            good[i] = np.arange(nt)
    # Set chi = 0 for equatorial antennas
    eq = np.array([8,9,10,12,13])
    chi = chi.astype(np.float)
    chi[:,eq] = 0.0
    chigood = np.zeros((nant,len(data['time'])),np.float)
    # Use only times corresponding to good values of chi, to find
    # nearest times to the data for each antenna
    for i in range(nant):
        idx = nearest_val_idx(data['time'],times[good[i]].jd)
        chigood[i] = chi[good[i][idx],i]
    outdata = np.zeros_like(data['x'])
    pol = np.array([0,2,3,1])  # Reorder polarizations to XX, XY, YX, YY
    for i in range(0,nant-1):
        for j in range(i+1,nant):
            for ifghz in range(nf):
                indict = {'data':data['x'][bl2ord[i,j],pol,ifghz],'chi1':chigood[i]+params[2],'chi2':chigood[j],'a1':1,'a2':params[0],'d1':0,'d2':params[1],'unrot': (not reverse)}
                outdata[bl2ord[i,j],:,ifghz] = frs.rot_sim(indict)
    return outdata
    
def unrot_miriad(filename,timeit=False,post=''):
    ''' Given a Miriad UV filename, read the parallactic angle information 
        from the stateframe database and unrotate the correlated data on
        each baseline.
    '''
    if timeit:
        import time
        ts = time.time()
    # Open existing Miriad file
    try:
        uvi = aipy.miriad.UV(filename)
    except:
        print 'Error opening Miriad file',filename
        return
    # Create new file in current directory (this will have to be changed to
    # something more rational after testing...)
    try:
        outfile = os.path.basename(filename)+post
        uvo = aipy.miriad.UV(outfile, status='new',corrmode='j')
        uvo.init_from_uv(uvi)
    except:
        print 'Could not open new Miriad file in current directory',outfile
        return
    # Do some sanity checks on this Miriad file
    npol = uvi['npol']
    if npol != 4:
        print 'Expecting four polarizations in this Miriad file, but got',npol
        return
    # Do test read of file to find out how many good frequencies:
    preamble, data = uvi.read()
    ifgood, = data.nonzero()
    uvi.rewind()
    nant = uvi['nants']
    nbl = nant*(nant-1)/2
    nbla = nbl + nant
    # Read the stateframe database and get the parallactic angle array
    # Clean up by removing any near-zero values
    import dbutil as db
    # Convert filename string to Julian Date
    jd = fname2mjd(filename) + 2400000.5
    # Set a timerange based on filename, ending 10 minutes after file start
    trange = Time([jd,jd+10./60./24.],format='jd')
    # Get parallactic angle for this timerange from SQL database
    times, chi = db.get_chi(trange)
    jd = times.jd
    nt, = times.shape
    fghz = uvi['sfreq']
    nf = len(fghz)
    nfg = len(ifgood)
    # Find only "good" chi values, i.e. those not too close
    # to zero, to eliminate missing values.
    good = []
    for i in range(nant):
        good.append(np.where(np.abs(chi[:,i]) > 0.001)[0])
        if len(good[i]) == 0:
            # This antenna must be offline, so go ahead and use all times.
            # They will be ignored anyway
            good[i] = np.arange(nt)
    # Set chi = 0 for equatorial antennas (and for non-existent ants 15 & 16)
    eq = np.array([8,9,10,12,13,14,15])
    chiant = chi.astype(np.float)
    chiant[:,eq] = 0.0
    # Next open the information on non-ideal feed behavior (these are constants for each
    # feed).  For now, just create an array with nominal values in place. This will become
    # a table, probably in the SQL database.  The values for each antenna are:
    #    a = Y feed amplitude, relative to X  [default is 1.0]
    #    d = proportion of cross-talk of Y relative to X (assumed symmetric) [default is 0.0]
    #    chi_0 = angle [degrees] of X and Y feeds relative to that of Ant 14 [default is 0]
    #    t = Y feed delay [ns], relative to X [default is 0.0]
    # Jones matrices for each antenna are the product of:
    #    A = [ 1      d     ]       R = [ cos(chi-chi_0)   sin(chi-chi_0]
    #        [-d  a*exp(iwt)]           [-sin(chi-chi_0)   cos(chi-chi_0]
    antdata = np.array([[1., 0., 0., 0.134],   # Ant 1
                        [1., 0., 0., 0.592],   # Ant 2
                        [1., 0., 0., 0.449],   # Ant 3
                        [1., 0., 0., 0.443],   # Ant 4
                        [1., 0., 0., 0.318],   # Ant 5
                        [1., 0., 0., -0.675],   # Ant 6
                        [1., 0., 0., 0.000],   # Ant 7
                        [1., 0., 0., 0.438],   # Ant 8
                        [1., 0., 0., 0.0],   # Ant 9
                        [1., 0., 0., 0.0],   # Ant 10
                        [1., 0., 0., 0.0],   # Ant 11
                        [1., 0., 0., -0.121],   # Ant 12
                        [1., 0., 0., 0.0],   # Ant 13
                        [1., 0., 0., 0.0],   # Ant 14
                        [1., 0., 0., 0.],   # Ant 15
                        [1., 0., 0., 0.]])   # Ant 16
    antdata[:,3] = 0.0  # ********Temporary override of delays, for test of effect of delay on results.
    # Create an array of matrices for each baseline, frequency and time
    M = np.zeros((nbl,4,4,nfg,nt),np.complex)
    nact = 14  # Actual number of antennas
    for i in range(0,nact-1):
        for j in range(i,nact):
            ai, di, chi_0i, ti = antdata[i]
            argi = 2*np.pi*ti*fghz
            aj, dj, chi_0j, tj = antdata[j]
            argj = 2*np.pi*tj*fghz
            # Calculate the feed angle difference for this baseline (for all times)
            Dchi = chiant[:,j] - chi_0j - (chiant[:,i] - chi_0i)
            # R has shape (2,2,nt)
            R = np.array([[np.cos(Dchi), -np.sin(Dchi)],[np.sin(Dchi), np.cos(Dchi)]])
            ai *= np.cos(argi)+1j*np.sin(argi)
            aj *= np.cos(argj)+1j*np.sin(argj)
            for k in range(nfg):
                Ai = np.array([[1,di],[-di,ai[k]]])
                # Jones matrix for first antenna of pair, with rotation and delay applied, shape (2,2,nt)
                JA = np.dot(Ai,R)
                # Jones matrix for second antenna of pair, with delay applied, shape (2,2)
                JB = np.array([[1,dj],[-dj,aj[k]]])
                for n in range(nt):
                     M[bl2ord[i,j],:,:,k,n] = np.kron(JA[:,:,n],np.conj(JB))

    if timeit:
        print 'Matrix complete after',time.time() - ts,'s'
    # Read a block of data for the current time, which consists of NBLA data samples
    # for each of 4 polarizations
    tprev = 0  # Indicates a time change, so that time lookup only happens then
    # Lists of premable-data pairs, one for each polarization
    XX = []
    YY = []
    XY = []
    YX = []
    cnt_in = -1
    cnt_out = 0
    for pd in uvi.all():
        cnt_in += 1
        if pd[0][1] != tprev:
            # Got a new time, so process current data and write it out
            if tprev != 0:
                # Find time index in M nearest to this time
                idx, = nearest_val_idx([tprev],jd)
                for ibl in range(len(XX)):
                    # Step through baselines, rotating the data according to the matrix
                    i0, j0 =  XX[ibl][0][2]
                    if i0 == j0:
                        # This is an auto-correlation, so nothing to do
                        pass
                    elif i0 >= nact or j0 >= nact:
                        # Skip over baselines with non-existent antennas
                        pass
                    else:
                        # Loop over non-zero frequencies, applying rotation matrix to each
                        for j,k in enumerate(XX[ibl][1].nonzero()[0]):
                            # Data vector for this baseline and frequency
                            data = np.array([XX[ibl][1][k],XY[ibl][1][k],YX[ibl][1][k],YY[ibl][1][k]])
                            try:
                                iM = np.linalg.inv(M[bl2ord[i0,j0],:,:,j,idx])
                                XX[ibl][1][k], XY[ibl][1][k], YX[ibl][1][k], YY[ibl][1][k] = np.dot(iM,data)
                                #if i0 == 0 and j0 == 13 and k == XX[ibl][1].nonzero()[0][0]:
                                #    print_calc(iM,data,tprev)
                            except:
                                # Case of singular matrix, generally because it is all zero
                                break
                uvo['pol'] = -5
                for ibl in range(nbla):
                    cnt_out += 1
                    uvo.write(XX[ibl][0],XX[ibl][1])
                uvo['pol'] = -6
                for ibl in range(nbla):
                    cnt_out += 1
                    uvo.write(YY[ibl][0],YY[ibl][1])
                uvo['pol'] = -7
                for ibl in range(nbla):
                    cnt_out += 1
                    uvo.write(XY[ibl][0],XY[ibl][1])
                uvo['pol'] = -8
                for ibl in range(nbla):
                    cnt_out += 1
                    uvo.write(YX[ibl][0],YX[ibl][1])
                #print '.',
                #if (cnt_out % (136*4*50)) == 0:
                #    print ' '
                XX = []
                YY = []
                XY = []
                YX = []
                # Things to do on getting a new time--update variables:
                uvo['ut'] = uvi['ut']
                if 'lst' in uvi.vartable:
                    uvo['lst'] = uvi['lst']
                else:
                    uvo['lst'] = 0.0
                uvo['xsampler'] = uvi['xsampler']
                uvo['ut'] = uvi['ut']
                uvo['ysampler'] = uvi['ysampler']
                uvo['delay'] = uvi['delay']
            tprev = pd[0][1]
            # Find nearest index in parallactic angles array, for this new time
            idx = nearest_val_idx([tprev],jd)
        # Separate data into the four polarizations
        if uvi['pol'] == -5:
            XX.append(pd)
        elif uvi['pol'] == -6:
            YY.append(pd)
        elif uvi['pol'] == -7:
            XY.append(pd)
        elif uvi['pol'] == -8:
            YX.append(pd)
        else:
            print 'Error-unrecognized polarization.'
            return
    # Don't forget to write out the last time sample
    # Find time index in M nearest to this time
    idx, = nearest_val_idx([tprev],jd)
    for ibl in range(len(XX)):
        # Step through baselines, rotating the data according to the matrix
        i0, j0 =  XX[ibl][0][2]
        if i0 == j0:
            # This is an auto-correlation, so nothing to do
            pass
        elif i0 >= nact or j0 >= nact:
            # Skip over baselines with non-existent antennas
            pass
        else:
            # Loop over non-zero frequencies, applying rotation matrix to each
            for j,k in enumerate(XX[ibl][1].nonzero()[0]):
                # Data vector for this baseline and frequency
                data = np.array([XX[ibl][1][k],XY[ibl][1][k],YX[ibl][1][k],YY[ibl][1][k]])
                try:
                    iM = np.linalg.inv(M[bl2ord[i0,j0],:,:,j,idx])
                    XX[ibl][1][k], XY[ibl][1][k], YX[ibl][1][k], YY[ibl][1][k] = np.dot(iM,data)
                    #if i0 == 0 and j0 == 13 and k == XX[ibl][1].nonzero()[0][0]:
                    #    print_calc(iM,data,tprev)
                except:
                    # Case of singular matrix, generally because it is all zero
                    break
    uvo['pol'] = -5
    for ibl in range(nbla):
        cnt_out += 1
        uvo.write(XX[ibl][0],XX[ibl][1])
    uvo['pol'] = -6
    for ibl in range(nbla):
        cnt_out += 1
        uvo.write(YY[ibl][0],YY[ibl][1])
    uvo['pol'] = -7
    for ibl in range(nbla):
        cnt_out += 1
        uvo.write(XY[ibl][0],XY[ibl][1])
    uvo['pol'] = -8
    for ibl in range(nbla):
        cnt_out += 1
        uvo.write(YX[ibl][0],YX[ibl][1])
    if timeit:
        print 'Done after',time.time() - ts,'s'
    return
  
def print_calc(M,data,t):
    res = np.dot(M,data)
    dat = data
    print Time(t,format='jd').iso[:19], '[',('{:6.3f} '*4).format(*M[0]),'] [ {:6.3f} ]   [ {:6.3f} ]'.format(dat[0],res[0])
    print '                    [',('{:6.3f} '*4).format(*M[1]),'] [ {:6.3f} ] = [ {:6.3f} ]'.format(dat[1],res[1])
    print '                    [',('{:6.3f} '*4).format(*M[2]),'] [ {:6.3f} ]   [ {:6.3f} ]'.format(dat[2],res[2])
    print '                    [',('{:6.3f} '*4).format(*M[3]),'] [ {:6.3f} ]   [ {:6.3f} ]'.format(dat[3],res[3])

def funrot_miriad(filename,timeit=False,post=''):
    ''' Given a Miriad UV filename, read the parallactic angle information 
        from the stateframe database and unrotate the correlated data on
        each baseline.  This is supposed to be the "fast" version that
        does this without doing delays, hence does not need to be done
        separately for each frequency.
    '''
    if timeit:
        import time
        ts = time.time()
    # Open existing Miriad file
    try:
        uvi = aipy.miriad.UV(filename)
    except:
        print 'Error opening Miriad file',filename
        return
    # Create new file in current directory (this will have to be changed to
    # something more rational after testing...)
    try:
        outfile = os.path.basename(filename)+post
        uvo = aipy.miriad.UV(outfile, status='new',corrmode='j')
        uvo.init_from_uv(uvi)
    except:
        print 'Could not open new Miriad file in current directory',outfile
        return
    # Do some sanity checks on this Miriad file
    npol = uvi['npol']
    if npol != 4:
        print 'Expecting four polarizations in this Miriad file, but got',npol
        return
    # Do test read of file to find out how many good frequencies:
    preamble, data = uvi.read()
    ifgood, = data.nonzero()
    uvi.rewind()
    nant = uvi['nants']
    nbl = nant*(nant-1)/2
    nbla = nbl + nant
    # Read the stateframe database and get the parallactic angle array
    # Clean up by removing any near-zero values
    import dbutil as db
    # Convert filename string to Julian Date
    jd = fname2mjd(filename) + 2400000.5
    # Set a timerange based on filename, ending 10 minutes after file start
    trange = Time([jd,jd+10./60./24.],format='jd')
    # Get parallactic angle for this timerange from SQL database
    times, chi = db.get_chi(trange)
    jd = times.jd
    nt, = times.shape
    fghz = uvi['sfreq']
    nf = len(fghz)
    nfg = len(ifgood)
    # Find only "good" chi values, i.e. those not too close
    # to zero, to eliminate missing values.
    good = []
    for i in range(nant):
        good.append(np.where(np.abs(chi[:,i]) > 0.001)[0])
        if len(good[i]) == 0:
            # This antenna must be offline, so go ahead and use all times.
            # They will be ignored anyway
            good[i] = np.arange(nt)
    # Set chi = 0 for equatorial antennas (and for non-existent ants 15 & 16)
    eq = np.array([8,9,10,12,13,14,15])
    chiant = chi.astype(np.float)
    #chiant = -chi.astype(np.float)  # Test reversing sign of chi...
    chiant[:,eq] = 0.0
    # Next open the information on non-ideal feed behavior (these are constants for each
    # feed).  For now, just create an array with nominal values in place. This will become
    # a table, probably in the SQL database.  The values for each antenna are:
    #    a = Y feed amplitude, relative to X  [default is 1.0]
    #    d = proportion of cross-talk of Y relative to X (assumed symmetric) [default is 0.0]
    #    chi_0 = angle [degrees] of X and Y feeds relative to that of Ant 14 [default is 0]
    # Jones matrices for each antenna are the product of:
    #    A = [ 1      d     ]       R = [ cos(chi-chi_0)   sin(chi-chi_0]
    #        [-d      a     ]           [-sin(chi-chi_0)   cos(chi-chi_0]
    antdata = np.array([[1., 0., 0.],   # Ant 1
                        [1., 0., 0.],   # Ant 2
                        [1., 0., 0.],   # Ant 3
                        [1., 0., 0.],   # Ant 4
                        [1., 0., 0.],   # Ant 5
                        [1., 0., 0.],   # Ant 6
                        [1., 0., 0.],   # Ant 7
                        [1., 0., 0.],   # Ant 8
                        [1., 0., 0.],   # Ant 9
                        [1., 0., 0.],   # Ant 10
                        [1., 0., 0.],   # Ant 11
                        [1., 0., 0.],   # Ant 12
                        [1., 0., 0.],   # Ant 13
                        [1., 0., 0.],   # Ant 14
                        [1., 0., 0.],   # Ant 15
                        [1., 0., 0.]])   # Ant 16
    # Create an array of matrices for each baseline, frequency and time
    iM = np.zeros((nbl,4,4,nt),np.complex)
    nact = 14  # Actual number of antennas
    for i in range(0,nact-1):
        for j in range(i,nact):
            ai, di, chi_0i = antdata[i]
            aj, dj, chi_0j = antdata[j]
            # Calculate the feed angle difference for this baseline (for all times)
            Dchi = chiant[:,j] - chi_0j - (chiant[:,i] - chi_0i)
            # R has shape (2,2,nt)
            R = np.array([[np.cos(Dchi), -np.sin(Dchi)],[np.sin(Dchi), np.cos(Dchi)]])
            Ai = np.array([[1,di],[-di,ai]])
            # Jones matrix for first antenna of pair, with rotation and delay applied, shape (2,2,nt)
            JA = np.dot(Ai,R)
            # Jones matrix for second antenna of pair, with delay applied, shape (2,2)
            JB = np.array([[1,dj],[-dj,aj]])
            for n in range(nt):
                try:
                    iM[bl2ord[i,j],:,:,n] = np.linalg.inv(np.kron(JA[:,:,n],np.conj(JB)))
                except:
                    iM[bl2ord[i,j],:,:,n] = np.array([[1,0],[0,1]])

    if timeit:
        print 'Matrix complete after',time.time() - ts,'s'
    # Read a block of data for the current time, which consists of NBLA data samples
    # for each of 4 polarizations
    tprev = 0  # Indicates a time change, so that time lookup only happens then
    # Lists of premable-data pairs, one for each polarization
    XX = []
    YY = []
    XY = []
    YX = []
    cnt_in = -1
    cnt_out = 0
    for pd in uvi.all():
        cnt_in += 1
        if pd[0][1] != tprev:
            # Got a new time, so process current data and write it out
            if tprev != 0:
                # Find time index in M nearest to this time
                idx, = nearest_val_idx([tprev],jd)
                for ibl in range(len(XX)):
                    # Step through baselines, rotating the data according to the matrix
                    i0, j0 =  XX[ibl][0][2]
                    if i0 == j0:
                        # This is an auto-correlation, so nothing to do
                        pass
                    elif i0 >= nact or j0 >= nact:
                        # Skip over baselines with non-existent antennas
                        pass
                    else:
                        # try:
                            # iM = np.linalg.inv(M[bl2ord[i0,j0],:,:,idx])
                        # except:
                            # # In case of singular matrix, just use the identity matrix
                            # iM = np.array([[1,0],[0,1]])
                        # Loop over non-zero frequencies, applying rotation matrix to each
                        for k in XX[ibl][1].nonzero()[0]:
                            # Data vector for this baseline and frequency
                            data = np.array([XX[ibl][1][k],XY[ibl][1][k],YX[ibl][1][k],YY[ibl][1][k]])
                            XX[ibl][1][k], XY[ibl][1][k], YX[ibl][1][k], YY[ibl][1][k] = np.dot(iM[bl2ord[i0,j0],:,:,idx],data)
                            #if i0 == 0 and j0 == 13 and k == XX[ibl][1].nonzero()[0][0]:
                            #    print_calc(iM,data,tprev)
                uvo['pol'] = -5
                for ibl in range(nbla):
                    cnt_out += 1
                    uvo.write(XX[ibl][0],XX[ibl][1])
                uvo['pol'] = -6
                for ibl in range(nbla):
                    cnt_out += 1
                    uvo.write(YY[ibl][0],YY[ibl][1])
                uvo['pol'] = -7
                for ibl in range(nbla):
                    cnt_out += 1
                    uvo.write(XY[ibl][0],XY[ibl][1])
                uvo['pol'] = -8
                for ibl in range(nbla):
                    cnt_out += 1
                    uvo.write(YX[ibl][0],YX[ibl][1])
                #print '.',
                #if (cnt_out % (136*4*50)) == 0:
                #    print ' '
                XX = []
                YY = []
                XY = []
                YX = []
                # Things to do on getting a new time--update variables:
                uvo['ut'] = uvi['ut']
                if 'lst' in uvi.vartable:
                    uvo['lst'] = uvi['lst']
                else:
                    uvo['lst'] = 0.0
                uvo['xsampler'] = uvi['xsampler']
                uvo['ut'] = uvi['ut']
                uvo['ysampler'] = uvi['ysampler']
                uvo['delay'] = uvi['delay']
            tprev = pd[0][1]
            # Find nearest index in parallactic angles array, for this new time
            idx = nearest_val_idx([tprev],jd)
        # Separate data into the four polarizations
        if uvi['pol'] == -5:
            XX.append(pd)
        elif uvi['pol'] == -6:
            YY.append(pd)
        elif uvi['pol'] == -7:
            XY.append(pd)
        elif uvi['pol'] == -8:
            YX.append(pd)
        else:
            print 'Error-unrecognized polarization.'
            return
    # Don't forget to write out the last time sample
    # Find time index in M nearest to this time
    idx, = nearest_val_idx([tprev],jd)
    for ibl in range(len(XX)):
        # Step through baselines, rotating the data according to the matrix
        i0, j0 =  XX[ibl][0][2]
        if i0 == j0:
            # This is an auto-correlation, so nothing to do
            pass
        elif i0 >= nact or j0 >= nact:
            # Skip over baselines with non-existent antennas
            pass
        else:
            # try:
                # iM = np.linalg.inv(M[bl2ord[i0,j0],:,:,idx])
            # except:
                # # In case of singular matrix, just use the identity matrix
                # iM = np.array([[1,0],[0,1]])
            # Loop over non-zero frequencies, applying rotation matrix to each
            for k in XX[ibl][1].nonzero()[0]:
                # Data vector for this baseline and frequency
                data = np.array([XX[ibl][1][k],XY[ibl][1][k],YX[ibl][1][k],YY[ibl][1][k]])
                XX[ibl][1][k], XY[ibl][1][k], YX[ibl][1][k], YY[ibl][1][k] = np.dot(iM[bl2ord[i0,j0],:,:,idx],data)
    uvo['pol'] = -5
    for ibl in range(nbla):
        cnt_out += 1
        uvo.write(XX[ibl][0],XX[ibl][1])
    uvo['pol'] = -6
    for ibl in range(nbla):
        cnt_out += 1
        uvo.write(YY[ibl][0],YY[ibl][1])
    uvo['pol'] = -7
    for ibl in range(nbla):
        cnt_out += 1
        uvo.write(XY[ibl][0],XY[ibl][1])
    uvo['pol'] = -8
    for ibl in range(nbla):
        cnt_out += 1
        uvo.write(YX[ibl][0],YX[ibl][1])
    if timeit:
        print 'Done after',time.time() - ts,'s'
    return

def applycal_miriad(filename,tcal=None,timeit=False):
    ''' Given a Miriad UV filename, and a Time() object corresponding to
        a reference calibration time, read the gain state of the refcal
        and the gain state for the period of the file, and apply the gain
        differences on each baseline, and write out a new file.
    '''
    import gaincal2 as gc
    one_minute = 60./86400.  # 1 minute, expressed in days
    if timeit:
        import time
        ts = time.time()
    # Get the gain state of the refcal
    trange = Time([tcal.iso,Time(tcal.lv+61,format='lv').iso])
    rc_gs =  gc.get_gain_state(trange)  # refcal gain state for 60 s
    # Get median of refcal gain state (which should be constant anyway)
    rc_gs['h1'] = np.median(rc_gs['h1'],1)
    rc_gs['h2'] = np.median(rc_gs['h2'],1)
    rc_gs['v1'] = np.median(rc_gs['v1'],1)
    rc_gs['v2'] = np.median(rc_gs['v2'],1)
    tsun = fname2mjd(filename)
    # Read ample timerange to ensure coverage
    trange = Time([tsun-2./86400.,tsun+602./86400.],format='mjd')
    sc_gs = gc.get_gain_state(trange)   # solar gain state for timerange of file
    # Open existing Miriad file
    try:
        uvi = aipy.miriad.UV(filename)
    except:
        print 'Error opening Miriad file',filename
        return
    # Create giant array of gains, translated to baselines and frequencies
    nt = len(sc_gs['times'])
    jd = Time(sc_gs['times'],format='lv').jd
    fghz = uvi['sfreq']
    nf = len(fghz)
    blist = (fghz*2 - 1).astype(int) - 1
    nact = 13  # Actual number of antennas to use
    antgain = np.zeros((nact,2,34,nt),float)   # Antenna-based gains vs. band
    gain = np.zeros((120,4,nf,nt),float)     # Baseline-based gains vs. frequency
    for i in range(nact):
        for j in range(34):
            antgain[i,0,j] = sc_gs['h1'][i] + sc_gs['h2'][i] - rc_gs['h1'][i] - rc_gs['h2'][i] + sc_gs['dcmattn'][i,0,j] - rc_gs['dcmattn'][i,0,j]
            antgain[i,1,j] = sc_gs['v1'][i] + sc_gs['v2'][i] - rc_gs['v1'][i] - rc_gs['v2'][i] + sc_gs['dcmattn'][i,1,j] - rc_gs['dcmattn'][i,1,j]
    for i in range(nact-1):
        for j in range(i+1,nact):
             gain[bl2ord[i,j],0] = 10**((antgain[i,0,blist] + antgain[j,0,blist])/20.)
             gain[bl2ord[i,j],1] = 10**((antgain[i,1,blist] + antgain[j,1,blist])/20.)
             gain[bl2ord[i,j],2] = 10**((antgain[i,0,blist] + antgain[j,1,blist])/20.)
             gain[bl2ord[i,j],3] = 10**((antgain[i,1,blist] + antgain[j,0,blist])/20.)
             
    # Create new file in current directory (this will have to be changed to
    # something more rational after testing...)
    try:
        outfile = os.path.basename(filename)
        uvo = aipy.miriad.UV(outfile, status='new',corrmode='j')
        uvo.init_from_uv(uvi)
    except:
        print 'Could not open new Miriad file in current directory',outfile
        return
    # Do some sanity checks on this Miriad file
    npol = uvi['npol']
    if npol != 4:
        print 'Expecting four polarizations in this Miriad file, but got',npol
        return
    # Do test read of file to find out how many good frequencies:
    preamble, data = uvi.read()
    ifgood, = data.nonzero()
    uvi.rewind()
    nant = uvi['nants']
    nbl = nant*(nant-1)/2
    nbla = nbl + nant

    if timeit:
        print 'Gain states have been calculated.',time.time() - ts,'s'
    # Read a block of data for the current time, which consists of NBLA data samples
    # for each of 4 polarizations
    tprev = 0  # Indicates a time change, so that time lookup only happens then
    # Lists of premable-data pairs, one for each polarization
    XX = []
    YY = []
    XY = []
    YX = []
    cnt_in = -1
    cnt_out = 0
    for pd in uvi.all():
        cnt_in += 1
        if pd[0][1] != tprev:
            # Got a new time, so process current data and write it out
            if tprev != 0:
                # Find time index in gain array nearest to this time
                idx, = nearest_val_idx([tprev],jd)
                for ibl in range(len(XX)):
                    # Step through baselines, rotating the data according to the matrix
                    i0, j0 =  XX[ibl][0][2]
                    if i0 == j0:
                        # This is an auto-correlation, so nothing to do
                        pass
                    elif i0 >= nact or j0 >= nact:
                        # Skip over baselines with non-existent antennas
                        pass
                    else:
                        # Loop over non-zero frequencies, applying gain factor to each
                        k = XX[ibl][1].nonzero()[0]
                        # Data vector for this baseline and frequency
                        XX[ibl][1][k] *= gain[bl2ord[i0,j0],0,k,idx]
                        YY[ibl][1][k] *= gain[bl2ord[i0,j0],1,k,idx]
                        XY[ibl][1][k] *= gain[bl2ord[i0,j0],2,k,idx]
                        YX[ibl][1][k] *= gain[bl2ord[i0,j0],3,k,idx]
                uvo['pol'] = -5
                for ibl in range(nbla):
                    cnt_out += 1
                    uvo.write(XX[ibl][0],XX[ibl][1])
                uvo['pol'] = -6
                for ibl in range(nbla):
                    cnt_out += 1
                    uvo.write(YY[ibl][0],YY[ibl][1])
                uvo['pol'] = -7
                for ibl in range(nbla):
                    cnt_out += 1
                    uvo.write(XY[ibl][0],XY[ibl][1])
                uvo['pol'] = -8
                for ibl in range(nbla):
                    cnt_out += 1
                    uvo.write(YX[ibl][0],YX[ibl][1])
                #print '.',
                #if (cnt_out % (136*4*50)) == 0:
                #    print ' '
                XX = []
                YY = []
                XY = []
                YX = []
                # Things to do on getting a new time--update variables:
                uvo['ut'] = uvi['ut']
                if 'lst' in uvi.vartable:
                    uvo['lst'] = uvi['lst']
                else:
                    uvo['lst'] = 0.0
                uvo['xsampler'] = uvi['xsampler']
                uvo['ut'] = uvi['ut']
                uvo['ysampler'] = uvi['ysampler']
                uvo['delay'] = uvi['delay']
            tprev = pd[0][1]
            # Find nearest index in parallactic angles array, for this new time
            idx = nearest_val_idx([tprev],jd)
        # Separate data into the four polarizations
        if uvi['pol'] == -5:
            XX.append(pd)
        elif uvi['pol'] == -6:
            YY.append(pd)
        elif uvi['pol'] == -7:
            XY.append(pd)
        elif uvi['pol'] == -8:
            YX.append(pd)
        else:
            print 'Error-unrecognized polarization.'
            return
    # Don't forget to write out the last time sample
    # Find time index in gain array nearest to this time
    idx, = nearest_val_idx([tprev],jd)
    for ibl in range(len(XX)):
        # Step through baselines, rotating the data according to the matrix
        i0, j0 =  XX[ibl][0][2]
        if i0 == j0:
            # This is an auto-correlation, so nothing to do
            pass
        elif i0 >= nact or j0 >= nact:
            # Skip over baselines with non-existent antennas
            pass
        else:
            # Loop over non-zero frequencies, applying gain factor to each
            k = XX[ibl][1].nonzero()[0]
            # Data vector for this baseline and frequency
            XX[ibl][1][k] *= gain[bl2ord[i0,j0],0,k,idx]
            YY[ibl][1][k] *= gain[bl2ord[i0,j0],1,k,idx]
            XY[ibl][1][k] *= gain[bl2ord[i0,j0],2,k,idx]
            YX[ibl][1][k] *= gain[bl2ord[i0,j0],3,k,idx]
    uvo['pol'] = -5
    for ibl in range(nbla):
        cnt_out += 1
        uvo.write(XX[ibl][0],XX[ibl][1])
    uvo['pol'] = -6
    for ibl in range(nbla):
        cnt_out += 1
        uvo.write(YY[ibl][0],YY[ibl][1])
    uvo['pol'] = -7
    for ibl in range(nbla):
        cnt_out += 1
        uvo.write(XY[ibl][0],XY[ibl][1])
    uvo['pol'] = -8
    for ibl in range(nbla):
        cnt_out += 1
        uvo.write(YX[ibl][0],YX[ibl][1])
    if timeit:
        print 'Done after',time.time() - ts,'s'
    return
