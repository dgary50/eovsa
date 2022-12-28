'''Reads and averages UDB, IDB files, hacked from
read_idb_.py. Replaces dump_tsys_ext.py. New version feb 2 2017, This
creates UDB files with the same variables as in IDB files (xsampler,
ysampler) rather than (xtyys, ytsys)'''
#jmm, 2017-02-03
#
# dg, 2017-08-30 -- Changes to allow these routines to work on IDB or UDB
#                   reprocessing via calls to readXdata() followed by udb_write()
# dg, 2021-01-10 -- Added autocorr_desat() routine to correct for correlator
#                   saturation that was discovered recently.  This is a big
#                   change, with more to come as we attempt to fix this.
# dg, 2021-01-12 -- I found a better way to determine the correction factor using 
#                   the SK m values to adjust for variations in power level due 
#                   to different numbers of subchannels. Also added a test for
#                   date that allows this to work for older and newer data.
# dg, 2021-03-18 -- Long-standing bug where the uv['ut'] (header) array variable 
#                   apparently has glitches.  Now the test for a bad record is 
#                   done using the variable t in each individual visibility record
#                   (which is how read_idb.py has always done it).
# dg, 2021-05-18 -- The correlator equalizer coefficient was changed from 8.0
#                   to 2.0, necessitating a change in the saturation correction
#                   factor in autocorr_desat().  This is applied to all data
#                   after 2021-05-16, when the change was made.

#needed for file creation
import time, os
import aipy
from astropy.io import fits
from util import Time
import util

import numpy as np
#data is a masked array
import numpy.ma as ma
#pcapture2 gives us a baseline array
import pcapture2 as p 
#eovsa_lst gives LST if it isn't present in the file
import eovsa_lst as el
#copy is used for filter option in idb_read
import copy
#to strip non-printable characters from antenna list
def strip_non_printable(string_in):
    ''' Only allow certain printable ascii characters to get through
    to fits write routine'''
#Keep values of 32 to 127
    stripped = (c for c in string_in if 31 < ord(c) < 127)
    return ''.join(stripped)
#End strip_non_printable

def avXdata(x, nsec=60):
    '''Averages UDB data over nsec seconds. '''
    # The input should be output from read_udb.readXdata
    one_day = 24.0*3600.0
    t = x['time']
    ntimes = len(t)
    dt = t[1::]-t[0:ntimes-1]
    dtsec = int(round(np.median(dt)*one_day))
    if dtsec >= nsec:
        print 'avXdata: Averaging time is too short, returning'
        return x
    #endif
    #time in seconds from the start
    #need a loop for some python reason
    tsec = t-t[0]
    for j in range(ntimes):
        tsec[j] = int(round(tsec[j]*one_day))
    #endfor
    #we'll create time bin edges here
    dtsec_all = int(np.max(tsec))
    nnew = 2+dtsec_all/nsec
    tsec_new = np.arange(0, nnew*nsec, nsec)
    #one less time than the bin edges
    ntnew = len(tsec_new)-1
    #define the output
    xx_shape = np.shape(x['x'])
    nf = xx_shape[0]
    nblc = xx_shape[1]
    npol = xx_shape[2]
    ntimes1 = xx_shape[3]
    #nf = np.size(x['x'][:,0,0,0])
    #nblc = np.size(x['x'][0,:,0,0])
    #npol = np.size(x['x'][0,0,:,0])
    outx0 = np.zeros((nf, nblc, npol, ntnew),dtype=np.complex64)
    omask = np.zeros((nf, nblc, npol, ntnew),dtype=np.int32)
    outx = ma.masked_array(outx0, mask = omask)
    nantsnf3 = np.size(x['px'][:,0])
    outpx = np.zeros((nantsnf3, ntnew), dtype=np.float32)
    outpy = np.zeros((nantsnf3, ntnew), dtype=np.float32)
    uvwarray = np.zeros((3, nblc, ntnew), dtype=np.float)
    lstarray = np.zeros(ntnew, dtype=np.float)
    nants = np.size(x['delay'][:,0])
    delayarray = np.zeros((nants, ntnew), dtype=np.float)
    utarray = np.zeros(ntnew, dtype=np.float)
    nsjarray = np.zeros(ntnew, dtype=np.int32)
    #step through and average, jc is the time bin that we're working
    #on, njc is the number of shorter intervals inside of this interval
    jc = 0
    njc = 0
    for j in range(ntnew):
        tj = tsec_new[j]
        tj1 = tsec_new[j+1]
        ssj = []
        ssj0 = []
        for i in range(ntimes):
            if tsec[i] == tj:
                ssj0.append(i)
            if tsec[i] >= tj and tsec[i] < tj1:
                ssj.append(i)
            #endif
        #endfor
#        print j, ssj0
#        print ssj
        if len(ssj) > 0:
            nsj = float(len(ssj))
            nsjarray[j] = nsj
            xxj = x['x'][:,:,:,ssj]
            outx[:,:,:,j]=ma.average(xxj, axis=3)
            pxj = x['px'][:,ssj]
            outpx[:,j]=np.sum(pxj, axis=1)
            pyj = x['py'][:,ssj]
            outpy[:,j]=np.sum(pyj, axis=1)
        #endif
    #endfor
    #time arrays, delays and uvw are interpolated: to interval center times
    tsec_mid = (tsec_new[1::]+tsec_new[0:ntnew])*0.5
    lstarray = np.interp(tsec_mid, tsec, x['lst'])
    utarray = np.interp(tsec_mid, tsec, x['ut'])
    for k in range(nants):
        delayk = x['delay'][k,:]
        delayarray[k, :] = np.interp(tsec_mid, tsec, delayk)
    #endfor
    for k in range(3):
        for i in range(nblc):
            uvwki = x['uvw'][k,i,:]
            uvwarray[k,i,:] = np.interp(tsec_mid, tsec, uvwki)
        #endfor
    #endfor
    tnew = t[0]+tsec_mid/one_day
    out = {'x':outx,'uvw':uvwarray,'time':tnew,'px':outpx,'py':outpy,'i0':x['i0'],
           'j0':x['j0'],'lst':lstarray,'pol':x['pol'],'delay':delayarray,'ut':utarray, 
           'file0':x['file0'], 'nsamples':nsjarray,'fghz':x['fghz']}
    return out
#END of avXdata

def udbfile_write(y, ufile_in, ufilename):
    '''Read in a UDB dataset average in time and write out the file. Y is
    the output from avXdata or readXdata, ufile_in is the input
    filename (needed for source, scan, etc...), ufilename is the
    output filename.

    '''

    if len(y) == 0:
        print 'udbfile_write: No data input'
        return []
    #endif
    if len(ufile_in) == 0:
        print 'udbfile_write: No file input'
        return []
    #endif
    if len(ufilename) == 0:
        print 'udbfile_write: No output file'
        return []
    #endif

    # Ready to output
    # Open the file and use that to replicate the NRV
    # (non-record-variable) variables
    print ufile_in
    uv = aipy.miriad.UV(ufile_in)
    if 'source' in uv.vartable: 
        src = uv['source']
    else: 
        src = 'None' #fix for bad files produced for GAINCALC 2018-04-06 to 2018-04-09, jmm
    #endif
    scanid = uv['scanid']
    nants = uv['nants']
    # The assumption here is that all the variables are going to be
    # there since it's already been processed
    uvout = aipy.miriad.UV(ufilename, 'new')

    uvout.add_var('name', 'a')
    uvout['name'] = strip_non_printable(ufilename)

    #handle source separately, jmm, 2018-04-09
    uvout.add_var('source', 'a')
    uvout['source'] = strip_non_printable(src)

    nrv_varlist_string = ['telescop', 'project', 'operator', 'version', 
                          'scanid', 'proj', 'antlist', 'obstype']
    for j in range(len(nrv_varlist_string)):
        uvout.add_var(nrv_varlist_string[j], 'a')
        uvout[nrv_varlist_string[j]] = strip_non_printable(uv[nrv_varlist_string[j]])
    #endfor

    nrv_varlist_int = ['nants', 'npol']
    for j in range(len(nrv_varlist_int)):
        uvout.add_var(nrv_varlist_int[j], 'i')
        uvout[nrv_varlist_int[j]] = uv[nrv_varlist_int[j]]
    #endfor

    nrv_varlist_rl = ['vsource', 'veldop', 'epoch']
    for j in range(len(nrv_varlist_rl)):
        uvout.add_var(nrv_varlist_rl[j], 'r')
        uvout[nrv_varlist_rl[j]] = uv[nrv_varlist_rl[j]]
    #endfor

    nrv_varlist_rl8 = ['freq', 'restfreq', 'antpos', 'ra', 'dec', 'obsra', 'obsdec']
    for j in range(len(nrv_varlist_rl8)):
        uvout.add_var(nrv_varlist_rl8[j], 'd')
        uvout[nrv_varlist_rl8[j]] = uv[nrv_varlist_rl8[j]]
    #endfor

    #sfreq, sdf, nchan, nspect don't change
    sfreq_in = uv['sfreq']
    na = len(sfreq_in)

    #add these vars
    uvout.add_var('nspect', 'i')
    uvout['nspect'] = uv['nspect']
    uvout.add_var('sfreq', 'd')
    uvout['sfreq'] = sfreq_in
    uvout.add_var('sdf', 'd')
    uvout['sdf'] = uv['sdf']

    #spectral windows
    nschan = np.ones(na, dtype=np.int32)
    ischan = np.arange(na, dtype=np.int32)+1
    uvout.add_var('nschan', 'i')
    uvout['nschan'] = nschan
    uvout.add_var('ischan', 'i')
    uvout['ischan'] = ischan

    #add a variable for ntimes
    uvout.add_var('ntimes', 'i')
    ntimes = len(y['time'])
    uvout['ntimes'] = ntimes
    # If input has nsamples key for the number of 1-s samples in each interval, 
    # add a variable for that, and one for integration time
    if 'nsamples' in y.keys():
        #Add a variable for delta_time in seconds
        dtsec = np.float(np.max(y['nsamples']))  # This will typically be 60.0
        uvout.add_var('inttime', 'r')
        uvout['inttime'] = dtsec
        uvout.add_var('nsamples', 'i')
        # Array of number of 1-s samples in each interval
        uvout['nsamples'] = y['nsamples']

    #define the record variables here
    uvout.add_var('ut', 'd')
    uvout.add_var('lst', 'd')
    uvout.add_var('xsampler', 'r')
    uvout.add_var('ysampler', 'r')
    uvout.add_var('delay', 'd')
    uvout.add_var('pol', 'i')

    #Need version info here
    version = "3.0"

    #Loop through times and add the other variables
    yy_shape = np.shape(y['x'])
    nf = yy_shape[0]
    nblc = yy_shape[1]
    npol = yy_shape[2]
    ntimes1 = yy_shape[3]
    for j in range(ntimes):
        tj = y['time'][j]
        utj = y['ut'][j]
        lstj = y['lst'][j]
        #odd things happen to xsampler, ysampler
        pxj = np.zeros(3*nants*nf, dtype = np.float32)
        pyj = np.zeros(3*nants*nf, dtype = np.float32)
        for k in range(3*nants*nf):
            pxj[k] = y['px'][k,j]
            pyj[k] = y['py'][k,j]
        #endfor
        dj = y['delay'][:,j]
        #xsampler
        uvout['ut'] = utj
        uvout['lst'] = lstj
        uvout['xsampler'] = pxj
        #ysampler
        uvout['ut'] = utj
        uvout['lst'] = lstj
        uvout['ysampler'] = pyj
        #delays
        uvout['delay'] = dj
        #for each polarization:
        for k in range(npol):
            uvout['pol'] = y['pol'][k]
            uvout['ut'] = utj
            uvout['lst'] = lstj
            #for each baseline
            for i in range(nblc):
                uvwij = y['uvw'][:,i,j]
                i0i = y['i0'][i]
                j0i = y['j0'][i]
                #this may work
                preamble = uvwij, tj, (i0i, j0i)
                data = y['x'][:,i,k,j]
                uvout.write(preamble, data)
            #endfor (baseline)
        #endfor (polarization)
    #endfor (time)

    del(uv) #done
    return ufilename

#End of udbfile_write

def readXdata(filename, filter=False, desat=False):
    '''This routine reads the data from a single IDB or UDB file.
       Optional Keywords: filter boolean--if True, returns only
       non-zero frequencies if False (default), returns all
       frequencies. This differs from Dale's version in that it
       includes all correlations, drops the tp_only option, and the
       outputs that are not in the UDB files.
       
       Added desat keyword so that saturation can be applied or not as desired.
    '''

    # Open uv file for reading
    print 'Processing: ', filename
    try:
        uv = aipy.miriad.UV(filename)
    except:
        print "UDB_UTIL.READXDATA: Bad File at initialzation: "+filename
        return []
    #endexcept
    # Read all to get a time array
    utcount = 0
    ut = 0.0
    try:
        for preamble, data in uv.all():
            # Look for time change
            if preamble[1] != ut:
                ut = preamble[1]
                utcount = utcount+1
            #endif
        #endfor
        uv.rewind()
    except:
        utcount = 0
        print "UDB_UTIL.READXDATA: Bad File: "+filename
    #endexcept
    if utcount == 0:
        print 'Returning: '
        return uv #done to fool the program to avoid segmentation fault -  core dump, jmm, 2019-08-08
    #endif
    nf_orig = len(uv['sfreq'])
    good_idx = np.arange(nf_orig)
    if filter:
        good_idx = []
        # Read a bunch of records to get number of good frequencies,
        # i.e. those with at least some non-zero data.  Read 20
        # records for baseline 1-2, XX pol
        uv.select('antennae',0,2,include=True)
        uv.select('polarization',-5,-5,include=True)
        for i in range(20):
            preamble, data = uv.read()
            idx, = data.nonzero()
            if len(idx) > len(good_idx):
                good_idx = copy.copy(idx)
        uv.select('clear',0,0)
        uv.rewind()
    #endif

    #set up outputs
    nf = len(good_idx)
#    print 'NF: ', nf
    freq = uv['sfreq'][good_idx]
    npol = uv['npol']
    polarr = np.array([-5, -6, -7, -8])
    nants = uv['nants']
    if 'nsamples' in uv.vartable:
        nsamples = uv['nsamples']
    else:
        nsamples = None
    nbl = nants*(nants-1)/2
    nblc = nbl+nants
    # all-correlations, add a mask for the output vis array
    outx0 = np.zeros((nf, nblc, npol, utcount),dtype=np.complex64)
    omask = np.zeros((nf, nblc, npol, utcount),dtype=np.int32)
    outx = ma.masked_array(outx0, mask = omask)
    i0array = np.zeros((nblc, utcount), dtype = np.int32)
    j0array = np.zeros((nblc, utcount), dtype = np.int32)
    outpx = np.zeros((3*nf*nants, utcount), dtype=np.float)
    outpy = np.zeros((3*nf*nants, utcount), dtype=np.float)
    uvwarray = np.zeros((3, nblc, utcount), dtype=np.float)
    delayarray = np.zeros((nants, utcount), dtype=np.float)
    #lists for time arrays
    utarray = []
    timearray = []
    lstarray = []
    l = -1
    tprev = 0
    tsav = 0
    # Use antennalist if available
    if 'antlist' in uv.vartable:
        ants = strip_non_printable(uv['antlist'])
        antlist = map(int, ants.split())
    else:
        antlist = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]
    #endelse
    #keep autocorrelations in the array
    bl2ord = p.bl_list()
    for ij in range(len(antlist)):
        bl2ord[ij, ij] = nbl+ij
    #endfor
    for preamble, data in uv.all():
        uvw, t, (i0,j0) = preamble
        i = antlist.index(i0+1)
        j = antlist.index(j0+1)
        if i > j:
            # Reverse order of indices
            j = antlist.index(i0+1)
            i = antlist.index(j0+1)
        #endif
        # Assumes uv['pol'] is one of -5, -6, -7, -8
        k = -5 - uv['pol']
        if filter:
            if len(data.nonzero()[0]) == nf:
                if t != tprev:
                    # New time 
                    if t == 2440587.5:
                        # Time is 1970-01-01, which means a zero-filled record, so skip
                        # the entire thing.
                        continue
                    l += 1
                    if l == utcount:
                        break
                    #endif
                    tprev = t
                    timearray.append(t)
                    utarray.append(uv['ut'])
                    try:
                        lstarray.append(uv['lst'])
                    except:
                        pass
                    #endexcept
                    xdata0 = uv['xsampler'].reshape(nf_orig,nants,3)
                    xdata = xdata0[good_idx, :, :]
                    outpx[:,l] = xdata.reshape(nf*nants*3)
                    ydata0 = uv['ysampler'].reshape(nf_orig,nants,3)
                    ydata = ydata0[good_idx, :, :]
                    outpy[:,l] = ydata.reshape(nf*nants*3)
                    delayarray[:,l] = uv['delay']
                #endif
                outx[:, bl2ord[i,j],k,l] = data[data.nonzero()]
                if k == 3:
                    uvwarray[:, bl2ord[i,j],l] = uvw
                    i0array[bl2ord[i,j],l] = i0
                    j0array[bl2ord[i,j],l] = j0
                #endif
            #endif, for len(data.nonzero()) == nf and uv['ut'] > 0
        else:
#            if uv['ut'] > 0:
            if t != tprev:
                # New time 
                if t == 2440587.5:
                    # Time is 1970-01-01, which means a zero-filled record, so skip
                    # the entire thing.
                    continue
                l += 1
                if l == utcount:
                    break
                #endif
                tprev = t
                timearray.append(t)
                utarray.append(uv['ut'])
                try:
                    lstarray.append(uv['lst'])
                except:
                    pass
                    #endexcept
                xdata = uv['xsampler']
                outpx[:,l] = xdata
                ydata = uv['ysampler']
                outpy[:,l] = ydata
                delayarray[:,l] = uv['delay']
            #endif
            outx[:,bl2ord[i,j],k,l] = data
            if k == 3:
                uvwarray[:,bl2ord[i,j],l] = uvw
                i0array[bl2ord[i,j],l] = i0
                j0array[bl2ord[i,j],l] = j0
            #endif 
#            #endif (uv['ut'] > 0)
        #endelse (not filter)
    #endfor
    # Truncate in case of early end of data, return if there is no good data
    nt = len(timearray)
    if nt == 0:
        out = []
    else:
        outpx = outpx[:,:nt]
        outpy = outpy[:,:nt]
        outx = outx[:,:,:,:nt]
        uvwarray = uvwarray[:, :, :nt]
        delayarray = delayarray[:, :nt]
        if len(lstarray) != 0:
            pass
        else:
            tarray = Time(timearray,format='jd')
            for t in tarray:
                lstarray.append(el.eovsa_lst(t))
            #endfor
        #endelse
        #i0 and j0 should always be the same
        i0array = i0array[:,0]
        j0array = j0array[:,0]
        bd = util.freq2bdname(freq, Time(timearray[0],format='jd'))
        #timearray, lstarray and utarray are lists
        out = {'x':outx,'uvw':uvwarray,'time':np.array(timearray),'px':outpx,'py':outpy,
               'i0':i0array,'j0':j0array,'lst':np.array(lstarray),'pol':polarr,
               'delay':delayarray,'ut':np.array(utarray),'file0':filename,'fghz':freq,'band':bd}
        if nsamples is None:
            pass
        else:
            out.update({'nsamples':nsamples})
    #endelse
    if desat:
        out = autocorr_desat(out)
    return out
#end of readXdata

def autocorr_desat(out):
    ''' Corrects for correlator saturation effects.  Applies a correction to 
        auto- and cross-correlation amplitudes based on total power amplitudes.
        
        Calculates the function eta = (x + d - c)/[(a*erf((x-c)/b) + d], where x = log(P),
        a,b,c,d = [1.22552, 1.37369, 2.94536, 2.14838], and erf() is the error function.
        However, eta is set to 1 for A < 50.
        
        Applies the function to autocorrelations A_i and cross-correlations xi_ij to obtain
        A'_i = A**eta_i and xi'_ij = xi_ij**[(eta_i + eta_j)/2].
        
        Since 2021-05-16, the correlator equalizer coefficient changed from 8.0 to 2.0,
        necessitating a different fit function.  This one combines two fits, one for
        auto-correlation amplitudes < 300 and another for > 300.
    '''
    from scipy.special import erf
    def eta_f(x, A):
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
        
    # Determine required "m" value for standardized power level.  The power changes
    # depending on number of channels averaged, etc., and the SK m value keeps track
    # of all of that.
    mjd = Time(out['time'][0],format='jd').mjd
    if mjd > 58536:
        m0 = 745472.  # Standard for most recent data (325 MHz bandwidth)
    else:
        m0 = 721536.  # Standard for earlier data (ca. 2017)
    bl2ord = util.bl2ord
    nant = 16
    nf, = out['fghz'].shape
    nt, = out['time'].shape
    npol, = out['pol'].shape
    Px = copy.deepcopy(out['px'].reshape(nf, nant, 3, nt))
    Py = copy.deepcopy(out['py'].reshape(nf, nant, 3, nt))
#    n0 = len(np.where(out['band'] == np.max(out['band']))[0])
#    # Modify measured power to correct for variable science-channel bandwidth
#    for i in range(nf):
#        bnd = out['band'][i]
#        n = np.float(len(np.where(out['band'] == bnd)[0]))
#        Px[i] *= n/n0
#        Py[i] *= n/n0
    Px = Px[:,:,0]*m0/Px[:,:,2]
    Py = Py[:,:,0]*m0/Py[:,:,2]
    
    x = np.log10(Px)
    # Calculate correction for X pol (returns size [nf, nant, nt])
    if mjd < 59350:   # 2021-05-16
        eta_x = eta_f(x, abs(out['x'][:,np.arange(120,136),0]))
    else:
        eta_x = eta_2(x, abs(out['x'][:,np.arange(120,136),0]))
    x = np.log10(Py)
    # Calculate correction for Y pol (returns size [nf, nant, nt])
    if mjd < 59350:   # 2021-05-16
        eta_y = eta_f(x, abs(out['x'][:,np.arange(120,136),1]))
    else:
        eta_y = eta_2(x, abs(out['x'][:,np.arange(120,136),1]))
    eta = np.zeros_like(out['x'])
    # Loop over baselines and calculate correction for different polarization states.
    for i in range(nant):
        for j in range(i,nant):
            eta[:,bl2ord[i,j],0] = eta_x[:,i]+eta_x[:,j]
            eta[:,bl2ord[i,j],1] = eta_y[:,i]+eta_y[:,j]
            eta[:,bl2ord[i,j],2] = eta_x[:,i]+eta_y[:,j]
            eta[:,bl2ord[i,j],3] = eta_x[:,j]+eta_y[:,i]
    eta = eta/2.
    #  Correction is to log of values, so apply by raising to eta power.
    amp = abs(out['x'])
    pha = np.angle(out['x'])
    amp = amp**eta
    out['x'] = amp*np.exp(1j*pha)
    return out
    
def concatXdata(x0, x):
    ''' Concatenates readXdata outputs'''
    #check for ok variable
    try:
        x0
    except:
        print 'udb_util.concatXdata: No initial input'
        return []
    #endexcept
    #check for ok variable
    try:
        x
    except:
        print 'udb_util.concatXdata: No concat input'
        return []
    #endexcept

    #Sometimes, the frequencies do not match -- typically this means
    #that the first x has crappy data, at least that is true in Jan
    #2017, jmm. So keep x and ditch x0
    xx_shape0 = np.shape(x0['x'])
    nf0 = xx_shape0[0]
    xx_shape = np.shape(x['x'])
    nf = xx_shape[0]
    if nf != nf0:
        print 'Frequency mismatch -- throwing out the first Xdata'
        return x
    #endif

    #vis array, is masked
    outx = ma.concatenate((x0['x'], x['x']), axis = 3)
    #uvw array
    uvwarray = np.concatenate((x0['uvw'], x['uvw']), axis = 2)
    #baseline to nt arrays
    i0array = x0['i0']
    j0array = x0['j0']
    #polarizations
    polarr = x0['pol']
    #power (sampler) arrays
    outpx = np.concatenate((x0['px'], x['px']), axis = 1)
    outpy = np.concatenate((x0['py'], x['py']), axis = 1)
    #delays
    delayarray = np.concatenate((x0['delay'], x['delay']), axis = 1)
    #times
    timearray = np.concatenate((x0['time'], x['time']))
    lstarray = np.concatenate((x0['lst'], x['lst']))
    utarray = np.concatenate((x0['ut'], x['ut']))
    #done
    out = {'x':outx,'uvw':uvwarray,'time':timearray,'px':outpx,'py':outpy,
           'i0':i0array,'j0':j0array,'lst':lstarray,'pol':polarr,'delay':delayarray,
           'ut':utarray,'file0':x0['file0'],'fghz':x0['fghz']}
    return out
#end of concatXdata

def valid_miriad_dataset(filelist0):
    '''Returns True or False for valid or invalid Miriad datasets,
    checks for existnce of the directory, and then for flags, header,
    vartable, and visdata. Also returns names of valid datasets, and
    invalid ones'''

    if len(filelist0) == 0:
        print 'valid_miriad_file: No files input'
        return False
    #endif

    #need a list input, otherwise all sorts of things are messed up
    if not isinstance(filelist0, list):
        filelist = [filelist0]
    else:
        filelist = filelist0
    #endelse

    n = len(filelist)
    otp = []
    ok_filelist = []
    bad_filelist = []
    for j in range(n):
        filename = filelist[j]
        tempvar = True
        if (os.path.isdir(filename) == False):
            tempvar = False
        if(os.path.isfile(filename+'/flags') == False or os.path.getsize(filename+'/flags') == 0):
            tempvar = False
        if(os.path.isfile(filename+'/header') == False or os.path.getsize(filename+'/header') == 0):
            tempvar = False
        if(os.path.isfile(filename+'/vartable') == False or os.path.getsize(filename+'/vartable') == 0):
            tempvar = False
        if(os.path.isfile(filename+'/visdata') == False or os.path.getsize(filename+'/visdata') == 0):
            tempvar = False
        #end ifs
        otp.append(tempvar)
        if tempvar == True:
            ok_filelist.append(filelist[j])
        #end if
        if tempvar == False:
            bad_filelist.append(filelist[j])
        #end if
    #endfor
    return otp, ok_filelist, bad_filelist
#End of valid_miriad_dataset

def udbfile_create(filelist, ufilename, nsec=60):
    '''Given a list of IDB filenames, create the appropriate UDB file, by
    averaging over energy bands, but keep 1 second time resolution'''
    print 'UDBFILE_CREATE: UFILENAME: ', ufilename

    if len(filelist) == 0:
        print 'udbfile_create: No files input'
        return [], []
    #endif

    # Be sure that files exist, and has all of the appropriate elements
    filelist_test, ok_filelist, bad_filelist = valid_miriad_dataset(filelist)
    if len(ok_filelist) == 0:
        print 'udbfile_create: No valid files input'
        return [], []
    #endif

    #For each file, read in the data, then concatenate and average
    bad_filename = []
    ufile_out = []
    fc = 0
    for filename in ok_filelist:
        xj = readXdata(filename, desat=True)
        #print 'Out of readXdata'

        #test for bad file, due to miriad bug, a bad file returns the
        #aipy.miriad.UV object itself, but a 0 filename will return []
        #print type(xj)
        if isinstance(xj, aipy.miriad.UV) == False:
            if len(xj) > 0:
                print 'concat :'+filename
                fc = fc+1
                if fc == 1:
                    x = xj
                    print x['x'].shape
                else:
                    x = concatXdata(x, xj)
                    print x['x'].shape
                #endelse
            else:
                print 'file skipped: ', filename
                return [], filename
            #endeles
        else:
            #as soon as you are here, aipy is toast. Pass out the
            #object, and let UDB_PROCESS handle the repercussions
            #by getting the bad IFDB entry out of the database.
            #Jmm, 2019-08-08
            print 'file skipped: ', filename
            return xj, filename
        #endeles
        #endif, error check, 2019-08-08, jmm
    #endfor
    #Now do the time average
    if fc == 0:
        print 'UDB_UTIL: No good data?'
        return ufilename, bad_filename
    #endif
    #average data here
    y = avXdata(x,nsec=nsec)
    print y['x'].shape
    #Now write the file
    ufile_out = udbfile_write(y, ok_filelist[0], ufilename)
    print 'UDBFILE_CREATE: UFILE_OUT: ', ufile_out

    return ufile_out, bad_filename
#End of udbfile_create

def xpx_comp(x):

    ''' Compares autocorrelations with Power calculations '''
    try:
        x
    except:
        print 'udb_util.xpx_comp: No input:'
        return []
    #endexcept
    
    xfactor = 64 #could be 16??
    px = x['px']
    py = x['py']

    nt = len(x['time'])
    nf = len(x['x'][:,0,0,0])

    px.shape = (nf, 16, 3, nt)
    py.shape = (nf, 16, 3, nt)

    M = px[:,0,2,0]/1792.0
    
    xcfrac = np.zeros((nf, 16, nt), dtype = np.float)
    ycfrac = np.zeros((nf, 16, nt), dtype = np.float)
    i = 0
    i1 = 1
    for kk in range(16):
        k = kk+120
        for n in range(nt):
            xcfrac[:,kk,n]  = x['x'][:,k,i,n]*M[:]/(xfactor*px[:,kk,0,n])
            ycfrac[:,kk,n] = x['x'][:,k,i1,n]*M[:]/(xfactor*py[:,kk,0,n])
        #endfor n
    #endfor kk

    out = {'xcfrac':xcfrac, 'ycfrac':ycfrac}
    return out
#End of xpx_comp
