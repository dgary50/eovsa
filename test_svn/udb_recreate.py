# !/usr/bin/env python 
#recreates a Miriad UDB file
import time, os
import numpy as np

def udbfile_recreate(ufile_in, ufilename):
    '''Read in a UDB dataset, correct the variables nschan, ischan,
    and write out the file'''

    import aipy
    if len(ufile_in) == 0:
        print 'udbfile_recreate: No file input'
        return []
    #endif

    # Replicate rd_miriad_tsys_file here: Open first file and use that
    # to replicate the NRV (non-record-variable) variables
    uv = aipy.miriad.UV(ufile_in)
    src = uv['source']
    scanid = uv['scanid']
    nants = uv['nants']
    # The assumption here is that all the variables are going to be
    # there since it's already been processed

    uvout = aipy.miriad.UV(ufilename, 'new')

    nrv_varlist_string = ['name', 'telescop', 'project', 'operator', 'version', 'source', 'scanid', 'proj', 'antlist', 'obstype']
    for j in range(len(nrv_varlist_string)):
        uvout.add_var(nrv_varlist_string[j], 'a')
        uvout[nrv_varlist_string[j]] = uv[nrv_varlist_string[j]]
    #endfor

    nrv_varlist_int = ['nants', 'npol']
    for j in range(len(nrv_varlist_int)):
        uvout.add_var(nrv_varlist_int[j], 'i')
        uvout[nrv_varlist_int[j]] = uv[nrv_varlist_int[j]]
    #endfor

    nrv_varlist_rl = ['vsource', 'veldop', 'inttime', 'epoch']
    for j in range(len(nrv_varlist_rl)):
        uvout.add_var(nrv_varlist_rl[j], 'r')
        uvout[nrv_varlist_rl[j]] = uv[nrv_varlist_rl[j]]
    #endfor

    nrv_varlist_rl8 = ['freq', 'restfreq', 'antpos', 'ra', 'dec', 'obsra', 'obsdec']
    for j in range(len(nrv_varlist_rl8)):
        uvout.add_var(nrv_varlist_rl8[j], 'd')
        uvout[nrv_varlist_rl8[j]] = uv[nrv_varlist_rl8[j]]
    #endfor

    #sfreq, sdf, nchan, nspect will not change
    sfreq_in = uv['sfreq']
    na = len(sfreq_in)
    print "shape of sfreq: ", np.shape(sfreq_in)

    #add these vars
    uvout.add_var('nspect', 'i')
    uvout['nspect'] = uv['nspect']
    uvout.add_var('sfreq', 'd')
    uvout['sfreq'] = sfreq_in
    uvout.add_var('sdf', 'd')
    uvout['sdf'] = uv['sdf']
    uvout.add_var('sfedg', 'd')
    uvout['sfedg'] = uv['sfedg']

    #spectral windows
    nschan = np.ones(na, dtype=np.int32)
    ischan = np.arange(1, na, dtype=np.int32)
    uvout.add_var('nschan', 'i')
    uvout['nschan'] = nschan
    uvout.add_var('ischan', 'i')
    uvout['ischan'] = ischan

    #define the record variables here
    uvout.add_var('ut', 'd')
    uvout.add_var('lst', 'd')
    uvout.add_var('xtsys', 'r')
    uvout.add_var('ytsys', 'r')
    uvout.add_var('delay', 'd')
    uvout.add_var('pol', 'i')

    #Need version info here
    if 'lst' in uv.vartable:
        version = "3.0"
    #endif

    # Loop and add the other variables
    init = False
    init_pol = False
    pol = -71 #dummy
    ut = 0.0 #not necessarily the same ut variable as in the preamble
    xcount = 0
    utcount = 0

    for preamble, data in uv.all():
        # Look for time change
        if preamble[1] != ut or init == False:
            ut = preamble[1]
            init = True
#xtsys
            uvout['ut'] = uv['ut']
            if version == "3.0":
                uvout['lst'] = uv['lst']
            else:
                uvout['lst'] = 0.0
            #endelse
            uvout['xtsys'] = uv['xtsys']
#ytsys
            uvout['ut'] = uv['ut']
            if version == "3.0":
                uvout['lst'] = uv['lst']
            else:
                uvout['lst'] = 0.0
            #endelse
            uvout['ytsys'] = uv['ytsys']
            uvout['delay'] = uv['delay']
            utcount = utcount+1
        #endif
#output polarization if it has changed:
        if uv['pol'] != pol or init_pol == False:
            init_pol = True
            pol = uv['pol']
            uvout['pol'] = pol
            uvout['ut'] = uv['ut']
            if version == "3.0":
                uvout['lst'] = uv['lst']
            else:
                uvout['lst'] = 0.0
            #endelse
        #endif
#Data
        uvout.write(preamble, data)
        xcount=xcount+1
    #endfor

    del(uv) #done
    print xcount, utcount
    return ufilename

#End of udbfile_recreate

def udbfile_comp(ufile1, ufile2):
    '''Compares two udb files'''

    import aipy
    if len(ufile1) == 0 or len(ufile2) == 0:
        print 'udbfile_comp: Need two files input'
        return
    #endif

    #Open files and compare vairables
    uv1 = aipy.miriad.UV(ufile1)
    uv2 = aipy.miriad.UV(ufile2)

    nrv_varlist_string = ['name', 'telescop', 'project', 'operator', 'version', 'source', 'scanid', 'proj', 'antlist', 'obstype']
    for j in range(len(nrv_varlist_string)):
        x1 = uv1[nrv_varlist_string[j]]
        x2 = uv2[nrv_varlist_string[j]]
        if x1 == x2:
            print nrv_varlist_string[j]+' matches'
        else:
            print nrv_varlist_string[j]+' does not match'
            print len(x1), x1
            print len(x2), x2
        #endelse
    #endfor
    nrv_varlist_int = ['nants', 'npol']
    for j in range(len(nrv_varlist_int)):
        x1 = uv1[nrv_varlist_int[j]]
        x2 = uv2[nrv_varlist_int[j]]
        if x1 == x2:
            print nrv_varlist_int[j]+' matches'
        else:
            print nrv_varlist_int[j]+' does not match'
            print x1
            print x2
        #endelse
    #endfor
    nrv_varlist_rl = ['vsource', 'veldop', 'inttime', 'epoch']
    for j in range(len(nrv_varlist_rl)):
        x1 = uv1[nrv_varlist_rl[j]]
        x2 = uv2[nrv_varlist_rl[j]]
        if x1 == x2:
            print nrv_varlist_rl[j]+' matches'
        else:
            print nrv_varlist_rl[j]+' does not match'
            print x1
            print x2
        #endelse
    #endfor
    nrv_varlist_rl8 = ['freq', 'restfreq', 'antpos', 'ra', 'dec', 'obsra', 'obsdec']
    for j in range(len(nrv_varlist_rl8)):
        x1 = uv1[nrv_varlist_rl8[j]]
        x2 = uv2[nrv_varlist_rl8[j]]
        if isinstance(x1, np.ndarray) == True:
            if np.array_equal(x1, x2) == True:
                print nrv_varlist_rl8[j]+' matches'
            else:
                print nrv_varlist_rl8[j]+' does not match'
                print x1
                print x2
            #endelse
        else:
            if x1 == x2:
                print nrv_varlist_rl8[j]+' matches'
            else:
                print nrv_varlist_rl8[j]+' does not match'
                print x1
                print x2
            #endelse
        #endelse
    #endfor
    nrv_varlist_ext = ['nspect', 'sfreq', 'sdf', 'sfedg', 'nschan', 'ischan']
    for j in range(len(nrv_varlist_ext)):
        x1 = uv1[nrv_varlist_ext[j]]
        x2 = uv2[nrv_varlist_ext[j]]
        if isinstance(x1, np.ndarray) == True:
            if np.array_equal(x1, x2) == True:
                print nrv_varlist_ext[j]+' matches'
            else:
                print nrv_varlist_ext[j]+' does not match'
                print x1
                print x2
            #endelse
        else:
            if x1 == x2:
                print nrv_varlist_ext[j]+' matches'
            else:
                print nrv_varlist_ext[j]+' does not match'
                print x1
                print x2
            #endelse
        #endelse
    #endfor
    count = 0
    rv_varlist = ['ut', 'xtsys', 'ytsys', 'delay', 'pol']
#    for preamble, data in uv1.all():
    for jkl in range(140):
        count = count+1
        preamble1, data1 = uv1.read()
        preamble2, data2 = uv2.read()
        pol1 = uv1['pol']
        pol2 = uv2['pol']
        print count, pol1, pol2
        for j in range(len(rv_varlist)):
            x1 = uv1[rv_varlist[j]]
            x2 = uv2[rv_varlist[j]]
            if isinstance(x1, np.ndarray) == True:
                if np.array_equal(x1, x2) == False:
                    print rv_varlist[j]+' does not match'
                    print count
                #endif
            else:
                if x1 != x2:
                    print rv_varlist[j]+' does not match'
                    print count
                    print x1, x2
                #endif
            #endelse
        #endfor
        
        if np.array_equal(preamble1[0], preamble2[0]) == False:
            print 'Preamble[0] does not match'
            print count
        #endif
        if preamble1[1] != preamble2[1]:
            print 'Preamble[1] does not match'
            print count
        #endif
        if preamble1[2][0] != preamble2[2][0]:
            print 'Preamble[2][0] does not match'
            print count
        #endif
        if preamble1[2][1] != preamble2[2][1]:
            print 'Preamble[2][1] does not match'
            print count
        #endif


        if np.array_equal(data1, data2) == False:
            print 'Data does not match'
            print count
        #endif
    #endfor

#End of udbfile_comp
