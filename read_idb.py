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
#  2017-Aug-09  DG
#    Change allday_udb() output plot to have fixed timerange, 13:30 - 02:30 UT.  Also
#    removed the import of spectrogram_fit and pcapture2, which load a lot of stuff for
#    very little purpose.  The pcapture2 routines bl2ord and ant_str2list were moved 
#    to util.py.  I just commented out the silly selfcal routine that is never used.
#  2017-Aug-11  DG
#    The allday_udb() plot was crashing due to GOES-15 data being all zero, so add
#    handling of that situation.  Also worked on the plot labeling a bit.
#  2018-Mar-18  DG
#    Attempt to make read_idb skip incompatible files (not matching shape of first file).
#  2019-Feb-20  DG
#    Fix long-standing error mode, where an empty datalist in read_idb() would cause a 
#    crash.  Now returns an empty dictionary.
#  2019-Jun-22  DG
#    Fix bandname list to work with either 34 or 52 band data, depending on date.  Dates
#    prior to 2019-Feb-22 use 34 bands.  Also, added allow_pickle=True in np.load() call
#    in read_npz().  Looks like a new security "feature" of Python.
#  2020-Jan-25  DG
#    Suppress warnings about imaginary total power data in a file, after the first one.
#  2020-05-10  DG
#    Updated cal_qual() to use util.get_idbdir() to find IDB root path.
#  2020-05-11  DG
#    Further update to make this work on the DPP.
#  2020-05-29  SY
#    Fix bug with datadir.find('eovsa') in case of using /nas3/IDB/
#  2020-12-13  DG
#    Important change to honor flags in the IDB files by setting flagged data to NaN.
#  2020-12-16  DG
#    Commented out irrelevant (and potentially dangerous, because out of date) read_udb() 
#    routine.  The standard read_idb() works fine for UDB files.
#  2021-01-10  DG
#    Important changes to add desat parameter to readXdata (and read_idb), to call the
#    new routine autocorr_desat() to correct for correlator saturation.  Default is False
#    right now, but will likely be changed to True after some tests.
#  2021-01-12  DG
#    I found a better way to determine the correction factor using the SK m values to
#    adjust for variations in power level due to different numbers of subchannels.
#    Also added a test for date that allows this to work for older and newer data.
#  2021-06-01  DG
#    It seems the IDB files occasionally have time glitches where the time jumps back,
#    I think due to resetting the network every 15 minutes.  Now readXdata() just skips
#    such records.
#  2022-01-30  DG
#    Found a bug in read_idb(), where the frequency list was from the last file read,
#    even if it did not match the first file.
#  2022-03-15  DG
#    Added an nmax parameter to read_idb() and readXdata() to override previous limitation
#    of reading only 600 times from a file.  With the new 20-ms files there can be 30000
#    records in a 10-min file!
#  2024-03-29  DG
#    Change unrot() to agree with the one in pipeline_cal.py.
#  2025-05-22  DG
#    Changed to work for 16 antennas.
#

import aipy
import os
from util import Time, nearest_val_idx, bl2ord, ant_str2list, common_val_idx, lobe, get_idbdir, freq2bdname
import glob
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.dates import DateFormatter
#import spectrogram_fit as sp
#import pcapture2 as p
import eovsa_lst as el
import copy
#import chan_util_bc as cu
#import chan_util_52 as cu52

#bl2ord = p.bl_list()

# def read_udb(filename):
    # ''' This routine reads the data from a UDB file.
    # '''

    # # Open uv file for reading
    # uv = aipy.miriad.UV(filename)
    # nf = len(uv['sfreq'])
    # nt = uv['ntimes']

    # print_warning = True     # Print a warning about imaginary total power data, if found.
    # freq = uv['sfreq']
    # npol = uv['npol']
    # nants = uv['nants']
    # nbl = nants*(nants-1)/2
    # outa = np.zeros((nants,npol,nf,nt),dtype=np.complex64)  # Auto-correlations
    # outx = np.zeros((nbl,npol,nf,nt),dtype=np.complex64)  # Cross-correlations
    # #outp = np.zeros((nants,2,nf,600),dtype=np.float)
    # #outp2 = np.zeros((nants,2,nf,600),dtype=np.float)
    # #outm = np.zeros((nants,2,nf,600),dtype=np.int)
    # uvwarray = np.zeros((nbl,nt,3),dtype=np.float)
    # timearray = []
    # lstarray = []
    # l = -1
    # tprev = 0
    # tsav = 0
    # # Use antennalist if available
    # ants = uv['antlist']
    # while ants[-1] == '\x00': ants = ants[:-1]
    # antlist = map(int, ants.split())

    # src = uv['source']
    # while src[-1] == '\x00': src = src[:-1]

    # for preamble, data in uv.all():
        # uvw, t, (i0,j0) = preamble
        # i = antlist.index(i0+1)
        # j = antlist.index(j0+1)
        # if i > j:
            # # Reverse order of indices
            # j = antlist.index(i0+1)
            # i = antlist.index(j0+1)
        # # Assumes uv['pol'] is one of -5, -6, -7, -8
        # k = -5 - uv['pol']
        # if t != tprev:
            # # New time 
            # l += 1
            # if l == nt:
                # break
            # tprev = t
            # timearray.append(t)
            # #xdata = uv['xsampler'].reshape(nf_orig,nants,3)
            # #ydata = uv['ysampler'].reshape(nf_orig,nants,3)
            # #outp[:,0,:,l] = np.swapaxes(xdata[:,:,0],0,1)
            # #outp[:,1,:,l] = np.swapaxes(ydata[:,:,0],0,1)
            # #outp2[:,0,:,l] = np.swapaxes(xdata[:,:,1],0,1)
            # #outp2[:,1,:,l] = np.swapaxes(ydata[:,:,1],0,1)
            # #outm[:,0,:,l] = np.swapaxes(xdata[:,:,2],0,1)
            # #outm[:,1,:,l] = np.swapaxes(ydata[:,:,2],0,1)

        # if i0 == j0:
            # # This is an auto-correlation
            # outa[i0,k,:,l] = data
            # if print_warning and k < 2 and np.sum(data != np.real(data)) > 0:
                # print preamble,uv['pol'], 'has imaginary data! Additional warnings suppressed.'
                # print_warning = False
        # else:
            # outx[bl2ord[i,j],k,:,l] = data
            # if k == 3: uvwarray[bl2ord[i,j],l] = uvw

    # # Truncate in case of early end of data
    # #nt = len(timearray)
    # #outp = outp[:,:,:,:nt]
    # #outp2 = outp2[:,:,:,:nt]
    # #outm = outm[:,:,:,:nt]
    # #if not tp_only:
    # #    outa = outa[:,:,:,:nt]
    # #    outx = outx[:,:,:,:nt]

    # if len(lstarray) != 0:
        # pass
    # else:
        # tarray = Time(timearray,format='jd')
        # for t in tarray:
            # lstarray.append(el.eovsa_lst(t))
    # ha = np.array(lstarray) - uv['ra']
    # ha[np.where(ha > np.pi)] -= 2*np.pi
    # ha[np.where(ha < -np.pi)] += 2*np.pi
    # out = {'a':outa, 'x':outx, 'uvw':uvwarray, 'fghz':freq, 'time':np.array(timearray),'source':src,'ha':ha,'ra':uv['ra'],'dec':uv['dec']}#,'p':outp,'p2':outp2,'m':outm
    # return out
    
def readXdata(filename, filter=False, tp_only=False, src=None, desat=False, nmax=600):
    ''' This routine reads the data from a single IDBfile.
        
        Optional Keywords:
        filter   boolean--if True, returns only non-zero frequencies 
                    if False (default), returns uniform set of 500 frequencies
        tp_only  boolean--if True, returns only TP information
                    if False (default), returns everything (including 
                    auto & cross correlations)
        nmax     max number of times to read from the file.  Defaults to 600,
                    which is the number of records in a 10-min file of 1-s
                    records.  This affects memory usage if it is much larger 
                    than the number of times in the file, but with the new 
                    20-ms files there can be 30000 recs in a 10-min file.
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
    print_warning = True    # Print a warning about imaginary total power data, if found.
    nf = len(good_idx)
    freq = uv['sfreq'][good_idx]
    npol = uv['npol']
    nants = uv['nants']
    nbl = nants*(nants-1)/2
    if not tp_only:
        outa = np.zeros((nants,npol,nf,nmax),dtype=np.complex64)  # Auto-correlations
        outx = np.zeros((nbl,npol,nf,nmax),dtype=np.complex64)  # Cross-correlations
    outp = np.zeros((nants,2,nf,nmax),dtype=np.float)
    outp2 = np.zeros((nants,2,nf,nmax),dtype=np.float)
    outm = np.zeros((nants,2,nf,nmax),dtype=np.int)
    uvwarray = np.zeros((nbl,nmax,3),dtype=np.float)
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
                    if t < tprev:
                        # Some kind of glitch, so skip the whole thing
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
                if t < tprev:
                    # Some kind of glitch, so skip the whole thing
                    continue
                l += 1
                if l == nmax:
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
                    outa[i0,k,:,l] = data.filled(fill_value=np.nan+np.nan*1j)
                    if print_warning and k < 2 and np.sum(data != np.real(data)) > 0:
                        print preamble,uv['pol'], 'has imaginary data! Additional warnings suppressed.'
                        print_warning = False
                else:
                    outx[bl2ord[i,j],k,:,l] = data.filled(fill_value=np.nan+np.nan*1j)
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
    bd = freq2bdname(freq, Time(timearray[0],format='jd'))
    out = {'a':outa, 'x':outx, 'uvw':uvwarray, 'fghz':freq, 'band':bd,'time':np.array(timearray),'source':src,'p':outp,'p2':outp2,'m':outm,'ha':ha,'ra':uv['ra'],'dec':uv['dec']}
    if desat:
        out = autocorr_desat(out)
    return out

def autocorr_desat(out):
    ''' Corrects for correlator saturation effects.  Applies a correction to 
        auto- and cross-correlation amplitudes based on total power amplitudes.
        
        Calculates the function eta = (x + d - c)/[(a*erf((x-c)/b) + d], where x = log(P),
        a,b,c,d = [1.22552, 1.37369, 2.94536, 2.14838], and erf() is the error function.
        However, eta is set to 1 for A < 50.
        
        Applies the function to autocorrelations A_i and cross-correlations xi_ij to obtain
        A'_i = A**eta_i and xi'_ij = xi_ij**[(eta_i + eta_j)/2].
    '''
    from scipy.special import erf    
    def eta_f(x, A):
        ''' This is the desaturation fuction for data taken with equalizer coefficient 8.0,
            which is data prior to 5/16/2021.
        '''
        a,b,c,d = [1.22552, 1.37369, 2.94536, 2.14838]  # Parameters define the invariant saturation curve
        eta = (x + d - c)/(a*erf((x-c)/b) + d)
        bad = np.where(A < 50)
        eta[bad] = 1.0
        return eta
        
    def eta_2(x, A):
        ''' This is the desaturation fuction for data taken with equalizer coefficient 2.0,
            which is data after to 5/16/2021.
        '''
        a1,b1,c1,d1 = [0.88025122, 1.0221639 , 4.39845723, 2.38911615]
        a2,b2,c2,d2 = [2.28517281, 2.64619331, 4.38657476, 2.37753165]
        hi = np.where(A > 300)
        low = np.where(A <= 300)
        eta = np.ones_like(x)
        eta[hi] = (x[hi] + d1 - c1)/(a1*erf((x[hi]-c1)/b1) + d1)
        eta[low] = (x[low] + d2 - c2)/(a2*erf((x[low]-c2)/b2) + d2)
        return eta
        
    nant = 16
    # Determine required "m" value for standardized power level.  The power changes
    # depending on number of channels averaged, etc., and the SK m value keeps track
    # of all of that.
    mjd = Time(out['time'][0],format='jd').mjd
    if mjd > 58536:
        m0 = 745472.  # Standard for most recent data (325 MHz bandwidth)
    else:
        m0 = 721536.  # Standard for earlier data (ca. 2017)
    nf, = out['fghz'].shape
    nt, = out['time'].shape
    npol = 4
    # Copy power for each polarization so modification below does not affect
    # actual data
    Px = out['p'][:,0]*m0/out['m'][:,0]
    Py = out['p'][:,1]*m0/out['m'][:,0]
#    n0 = len(np.where(out['band'] == np.max(out['band']))[0])
#    # Modify measured power to correct for variable science-channel bandwidth
#    for i in range(nf):
#        bnd = out['band'][i]
#        n = np.float(len(np.where(out['band'] == bnd)[0]))
#        Px[:,i] *= n/n0
#        Py[:,i] *= n/n0
    x = np.log10(Px)
    # Calculate correction for X pol (returns size [nant, nf, nt])
    if mjd < 59350:  # 2021-05-16
        eta_x = eta_f(x, abs(out['a'][:,0]))
    else:
        eta_x = eta_2(x, abs(out['a'][:,0]))

    x = np.log10(Py)
    # Calculate correction for Y pol (returns size [nant, nf, nt])
    if mjd < 59350:  # 2021-05-16
        eta_y = eta_f(x, abs(out['a'][:,1]))
    else:
        eta_y = eta_2(x, abs(out['a'][:,1]))
    eta = np.zeros_like(out['x'])
    # Loop over baselines for cross-correlation and calculate correction for 
    # different polarization states.
    for i in range(nant-1):
        for j in range(i+1,nant):
            eta[bl2ord[i,j],0] = eta_x[i]+eta_x[j]
            eta[bl2ord[i,j],1] = eta_y[i]+eta_y[j]
            eta[bl2ord[i,j],2] = eta_x[i]+eta_y[j]
            eta[bl2ord[i,j],3] = eta_x[j]+eta_y[i]
    eta = eta/2.
    # Correction is to log of values, so apply by raising amplitudes to eta power,
    # but preserve the phase
    amp = abs(out['x'])
    pha = np.angle(out['x'])
    amp = amp**eta
    out['x'] = amp*np.exp(1j*pha)
    # Repeat for auto-correlation by looping over antennas
    # This "re-uses" eta to overwrite first nant values
    for i in range(nant):
        eta[i,0] = eta_x[i]
        eta[i,1] = eta_y[i]
        eta[i,2] = (eta_x[i] + eta_y[i])/2.
        eta[i,3] = (eta_x[i] + eta_y[i])/2.
    amp = abs(out['a'])
    pha = np.angle(out['a'])
    amp = amp**eta[:nant]
    out['a'] = amp*np.exp(1j*pha)
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
    
def summary_plot(out,ant_str='ant1-15',ptype='phase',pol='XX-YY'):
    ''' Makes a summary amplitude or phase plot for all baselines from ants in ant_str
        in out dictionary.
    '''
    import matplotlib.pyplot as plt
    
    ant_list = ant_str2list(ant_str)
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

def summary_plot_pcal(out,ant_str='ant1-16',ptype='phase',pol='XX-YY'):
    ''' Makes a summary amplitude or phase plot for all baselines from ants in ant_str
        in out dictionary.
    '''
    import matplotlib.pyplot as plt
    ha=out['ha']
    fghz=out['fghz']
    ant_list = ant_str2list(ant_str)
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
            
def read_idb(trange,navg=None, nmax=600, quackint=0.,filter=True,srcchk=True,src=None,tp_only=False, desat=False):
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
            out = readXdata(file,tp_only=tp_only,src=src, desat=desat, nmax=nmax)
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
            
    if len(datalist) == 0:
        return {}
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
    # Keep track of files whose shape matches the first file
    shape1 = datalist[0]['p'].shape[:-1]
    match = []
    for i,out in enumerate(datalist):
        shape2 = out['p'].shape[:-1]
        if shape1 == shape2:
            match.append(True)
            outp.append(out['p'])
        else:
            match.append(False)
            print 'Scan/file',i+1,'skipped. Array shape',shape2,'does not match shape',shape1,'of first scan/file'
    out['p'] = np.concatenate(outp,3)
    if filter:
        # Eliminate frequencies where there is no nonzero value
        # sums power over every dimension except freq.
        goodidx, = np.sum(np.sum(np.sum(out['p'],3),1),0).nonzero()
        # Eliminate frequencies where there is no nonzero value
        out['p'] = out['p'][:,:,goodidx]
    for i,out in enumerate(datalist):
        if match[i]:
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
    out['fghz'] = datalist[0]['fghz']
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
        data = np.load(f, allow_pickle=True)
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

#def fname2mjd(filename):
#    fstem = filename.split('/')[-1]
#    fstr = fstem[3:7]+'-'+fstem[7:9]+'-'+fstem[9:11]+' '+fstem[11:13]+':'+fstem[13:15]+':'+fstem[15:17]
#    t = Time(fstr)
#    return t.mjd

def get_trange_files(trange):
    #Given a timerange, this routine will take all relevant IDBfiles from
    #  that time range, put them in a list, and return that list.
    #  This function is used in get_X_data(data).
    from util import get_idbdir, fname2mjd
    
    fstr = trange[0].iso
    # Get path to root of IDB data
    datadir = get_idbdir(trange[0])

    # Add date path if on pipeline
    import socket
    host = socket.gethostname()
    if host == 'pipeline': datadir += fstr.replace('-','').split()[0]+'/'
    # if datadir.find('eovsa') != -1: datadir += fstr.replace('-','').split()[0]+'/'
    folder=datadir
    try:
        os.listdir(folder)
    except:
        print 'Something wrong with the definition of path to root of IDB files.'
        print 'See util.get_idbdir() for details.'
        return

    files = glob.glob(folder+'IDB'+fstr.replace('-','').split()[0]+'*')
    files.sort()
    mjd1, mjd2 = trange.mjd.astype('int')
    if mjd2 != mjd1:
        if (mjd2 - 1) != mjd1:
            usage('Second date must differ from first by at most 1 day')
        else:
            fstr2 = trange[1].iso
            files2 = glob.glob(folder+'IDB'+fstr2.replace('-','').split()[0]+'*')
            files2.sort()
            files += files2

#    def fname2mjd(filename):
#        fstem = filename.split('/')[-1]
#        fstr = fstem[3:7]+'-'+fstem[7:9]+'-'+fstem[9:11]+' '+fstem[11:13]+':'+fstem[13:15]+':'+fstem[15:17]
#        t = Time(fstr)
#        return t.mjd

    filelist = []
    for filename in files:
        mjd = fname2mjd(filename)
        if mjd >= trange[0].mjd and mjd < trange[1].mjd:
            filelist.append(filename)
    return filelist
    

def unrot(data, azeldict=None):
    ''' Apply the correction to differential feed rotation to data, and return
        the corrected data.  This also applies flags to data whose antennas are
        not tracking.

        Inputs:
          data     A dictionary returned by read_idb.py's readXdata().
          azeldict The dictionary returned from get_sql_info(), or if None, the appropriate
                     get_sql_info() call is done internally.

        Output:
          cdata    A dictionary with the phase-corrected data.  Only the key
                     x is updated.
    '''
    import copy
    from pipeline_cal import get_sql_info
    import cal_header as ch
    from stateframe import extract

    trange = Time(data['time'][[0, -1]], format='jd')

    if azeldict is None:
        azeldict = get_sql_info(trange)
    chi = azeldict['ParallacticAngle'] * np.pi / 180.  # (nt, nant)
    # Correct parallactic angle for equatorial mounts, relative to Ant A
    if trange[0] < Time('2025-05-22'):
        chi[:, [8, 9, 10, 12, 13]] = 0  # Currently 0, but can be measured and updated
        nant = 15
    else:
        chi[:, 15] = 0   # Only Ant A is equatorial after 2025-05-22
        nant = 16

    # Which antennas are tracking
    track = azeldict['TrackFlag']  # True if tracking

    # Ensure that nearest valid parallactic angle is used for times in the data
    good = np.where(azeldict['ActualAzimuth'] != 0)
    tidx = []  # List of arrays of indexes for each antenna
    gd = []
    for i in range(nant):
        gd.append(good[0][np.where(good[1] == i)])
        tidx.append(nearest_val_idx(data['time'], azeldict['Time'][gd[i]].jd))

    # Read X-Y Delay phase from SQL database and get common frequencies
    xml, buf = ch.read_cal(11, t=trange[0])
    fghz = extract(buf, xml['FGHz'])
    good, = np.where(fghz != 0.)
    fghz = fghz[good]
    dph = extract(buf, xml['XYphase'])
    dph = dph[:, good]
    xi_rot = extract(buf, xml['Xi_Rot'])
    xi_rot = xi_rot[good]
    fidx1, fidx2 = common_val_idx(data['fghz'], fghz, precision=4)
    missing = np.setdiff1d(np.arange(len(data['fghz'])), fidx1)

    nbl, npol, nf, nt = data['x'].shape
    nf = len(fidx1)
    # Correct data for X-Y delay phase
    for i in range(nant-1):
        for j in range(i + 1, nant):
            k = bl2ord[i, j]
            if j == 13:                  # xi_rot was applied for all antennas, but this
                xi = xi_rot[fidx2]       # is wrong.  Now it is only done for ant14.
            else:
                xi = 0.0                 # xi_rot for other antennas is just zero.
            a1 = lobe(dph[i, fidx2] - dph[j, fidx2])
            a2 = -dph[j, fidx2] - xi
            a3 = dph[i, fidx2] - xi + np.pi
            data['x'][k, 1, fidx1] *= np.repeat(np.exp(1j * a1), nt).reshape(nf, nt)
            data['x'][k, 2, fidx1] *= np.repeat(np.exp(1j * a2), nt).reshape(nf, nt)
            data['x'][k, 3, fidx1] *= np.repeat(np.exp(1j * a3), nt).reshape(nf, nt)

    # Correct data for differential feed rotation
    cdata = copy.deepcopy(data)
    for n in range(nt):
        for i in range(nant-1):
            for j in range(i + 1, nant):
                k = bl2ord[i, j]
                ti = tidx[i][n]
                tj = tidx[j][n]
                if track[ti, i] and track[tj, j]:
                    # If something goes wrong with chi difference calculation, just default to chi = 0
                    try:
                        dchi = chi[gd[i][ti], i] - chi[gd[j][tj], j]
#                        dchi = -chi[gd[i][ti], i] + chi[gd[j][tj], j]  # Try opposite sign of rotation *****
                    except:
                        dchi = 0.0
                    cchi = np.cos(dchi)
                    schi = np.sin(dchi)
                    cdata['x'][k, 0, :, n] = data['x'][k, 0, :, n] * cchi + data['x'][k, 3, :, n] * schi
                    cdata['x'][k, 2, :, n] = data['x'][k, 2, :, n] * cchi + data['x'][k, 1, :, n] * schi
                    cdata['x'][k, 3, :, n] = data['x'][k, 3, :, n] * cchi - data['x'][k, 0, :, n] * schi
                    cdata['x'][k, 1, :, n] = data['x'][k, 1, :, n] * cchi - data['x'][k, 2, :, n] * schi
                else:
                    cdata['x'][k, :, :, n] = np.nan

    # Set flags for any missing frequencies (hopefully this also works when "missing" is np.array([]))
    # cdata['x'][missing] = np.ma.masked
    return cdata

