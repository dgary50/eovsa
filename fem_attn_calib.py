# !/usr/bin/env python 
# Hacked from Jim's dump_tsys_ext.py (pipeline:test_svn/python/, version 2016_may_25)
# 2016-08-01 BC: added rd_fem_attn() to read the FEM attenuation records fron the stateframe
#                added cal_fem_gain() to convert FEM attenuation records to gain change, 
#                       by taking a fem_attn_inc parameter from appropriate calibration
#                added calc_fem_attn_inc() to calculate the additional FEM attenuation calibrations
#                modified udbfile_create() to incorporate FEM attenuation corrections 
import os
import numpy as np
from util import Time
import dump_tsys_ext
import dbutil as db
from util import Time
import read_idb as ri
import pdb

def fem_attn_anal(idb_calib='/dppdata1/IDB/IDB20160731231934/',doplot=False, wrt2sql=False):
    import cal_header as ch
    '''Calculate additional corrections to the FEM attenuators at each bit change (0, 1, 2, 4, 8, 16)dB;
       Values are based on the measurement IDB20160731231934, the sequence is specified in 
       helios:Dropbox/PythonCode/Current/FEATTNTEST2.ctl:
       1. FEMAUTO-OFF, 2. FEMATTN 15 (both 31 dB, to get the bkg), 
       3. Change the 1st H and V attn to 0, 1, 2, 4, 8, 16 dB every 30s, while keeping the 2nd to be 8 dB,
       4. Change the 2nd H and V attn to 0, 1, 2, 4, 8, 16 dB every 30s, while keeping the 1st to be 8 dB

       return value: 
            fem_attn_inc (nant, npol, # of FEM attn, 5): 
                        additional corrections in dB w.r.t. the nominal values when bit changes'''
    import matplotlib.pyplot as plt
    out=ri.read_idb([idb_calib])
    nant, npol, nf, nt=out['p'].shape
    if doplot:
        # show the auto-correlation, use it to find time indices for each attn state
        f, ax = plt.subplots(5,3)
        for i in range(15):
            ax[i / 3, i % 3].imshow(out['p'][i,0])
    # define time idx ranges
    # each state lasted 30 s
    tidxs=[25+i*30 for i in range(13)] #begin idx for avg
    tidxe=[idx+20 for idx in tidxs] #end idx for avg
    bkg=np.mean(out['p'][:,:,:,tidxs[0]:tidxe[0]],axis=3)
    # measurements for 12 attn states
    p_1=np.zeros([nant,npol,nf,6]) #power values for 6 states of 1st FEM attn change
    p_2=np.zeros([nant,npol,nf,6]) #power values for 6 states of 2nd FEM attn change
    rp_1=np.zeros([nant,npol,nf,6]) #ratio of the power regarding to the reference state
    rp_2=np.zeros([nant,npol,nf,6])
    rdb_1=np.zeros([nant,npol,nf,6]) #power ratio converted to dB
    rdb_2=np.zeros([nant,npol,nf,6])
    # nominal dB values
    attns1=np.array([0.,1.,2.,4.,8.,16.])
    attns2=np.array([0.,1.,2.,4.,8.,16.])
    # reference state
    attn_idx_ref=0 # 0 dB state
    attns1-=attns1[attn_idx_ref]
    attns2-=attns2[attn_idx_ref]
    for i in range(12):
        if i < 6:
            p_1[:,:,:,i]=np.mean(out['p'][:,:,:,tidxs[i+1]:tidxe[i+1]],axis=3)-bkg
        else:
            p_2[:,:,:,i-6]=np.mean(out['p'][:,:,:,tidxs[i+1]:tidxe[i+1]],axis=3)-bkg

    for i in range(6):
        rp_1[:,:,:,i]=p_1[:,:,:,i]/p_1[:,:,:,attn_idx_ref]
        rp_2[:,:,:,i]=p_2[:,:,:,i]/p_2[:,:,:,attn_idx_ref]
    rdb_1=-10.*np.log10(rp_1)
    rdb_2=-10.*np.log10(rp_2)
    ddb_1=rdb_1-attns1 # additional dB correction wrt the nominal values, FEM attn 1
    ddb_2=rdb_2-attns2 # additional dB correction wrt the nominal values, FEM attn 2
    # plot the additional dB correction 
    if doplot:
        f1, ax1 = plt.subplots(6,5)
        for i in range(15):
            ax1[i / 5, i % 5].imshow(ddb_1[i,0],vmin=-2,vmax=2)
            ax1[i / 5+3, i % 5].imshow(ddb_1[i,1],vmin=-2,vmax=2)
        f2, ax2 = plt.subplots(6,5)
        for i in range(15):
            ax2[i / 5, i % 5].imshow(ddb_2[i,0],vmin=-2,vmax=2)
            ax2[i / 5+3, i % 5].imshow(ddb_2[i,1],vmin=-2,vmax=2)

    # generate corrections for the nominal values
    chran=[0,90] #range of frequency channels to average
    fem_attn_bitv=np.zeros([16,npol,2,5],dtype=np.complex)
    fem_attn_bitv[:,:,0,:]=np.nan_to_num(np.mean(rdb_1[:,:,:,1:],axis=2))+0j
    fem_attn_bitv[:,:,1,:]=np.nan_to_num(np.mean(rdb_2[:,:,:,1:],axis=2))+0j
    if wrt2sql:
        ch.fem_attn_val2sql(fem_attn_bitv,t=Time(out['time'][0],format='jd'))
    return fem_attn_bitv

def fem_attn_update(fem_attn, t=None, rdfromsql=True):
    '''Given a record of the frontend attenuation levels from the stateframe, recalculate the corrected attenuations levels.
       fem_attn_in: recorded attn levels in a 10-min duration
       fem_attn_bitv: complex corrections to be applied to the data. Read from the stateframe or provided as a (16, 2, 2, 5) array'''
    import cal_header as ch
    import stateframe as stf
    if rdfromsql:
        xml, buf = ch.read_cal(7,t)
        fem_attn_bitv=np.nan_to_num(stf.extract(buf, xml['FEM_Attn_Real'])) + np.nan_to_num(stf.extract(buf, xml['FEM_Attn_Imag'])) * 1j
     
    h1=fem_attn['h1']
    h2=fem_attn['h2']
    v1=fem_attn['v1']
    v2=fem_attn['v2']
    attn=np.concatenate((np.concatenate((h1[...,None],v1[...,None]),axis=2)[...,None],
                         np.concatenate((h2[...,None],v2[...,None]),axis=2)[...,None]),axis=3)
    # Start with 0 attenuation as reference
    fem_attn_out=attn*0
    # Calculate resulting attenuation based on bit attn values (1,2,4,8,16)
    for i in range(5):
        fem_attn_out = fem_attn_out + (np.bitwise_and(attn,2**i)>>i)*fem_attn_bitv[...,i]
    #fem_gain=10.**(-fem_gain_db/10.)
    return fem_attn_out

def udbfile_create(filelist, ufilename, verbose=False):
    '''Given a list of IDB filenames, create the appropriate UDB file, by
    averaging over energy bands, but keep 1 second time resolution. 
    FEM attn values and corrections are taken from the stateframe'''

    import aipy
    if len(filelist) == 0:
        print 'udbfile_create: No files input'
        return []
    #endif

    # Be sure that files exist, and has all of the appropriate elements
    filelist_test, ok_filelist, bad_filelist = dump_tsys_ext.valid_miriad_dataset(filelist)
    if len(ok_filelist) == 0:
        print 'udbfile_create: No valid files input'
        return []
    #endif

    # Replicate rd_miriad_tsys_file here: Open first file and use that
    # to replicate the NRV (non-record-variable) variables
    uv = aipy.miriad.UV(ok_filelist[0])
    src = uv['source']
    scanid = uv['scanid']
    nants = uv['nants']
    # The assumption here is that all the variables are going to be
    # there since the valid_miriad_dataset was invoked

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

    #sfreq, sdf, nchan, nspect will change
    #navg = 10
    navg = 1
    sfreq_in = uv['sfreq']
    nchan_in = len(sfreq_in)
    nch_avg = nchan_in/navg #will crash if this is not an integer

    #arange gives me an array
    a = np.arange(0, nchan_in, navg)
    b = a+navg-1
    #strip out zero frequencies
    ppp = np.where(sfreq_in[a] > 0)
    a = a[ppp]
    b = b[ppp]
    na = len(a)
    indexmax = np.amax(np.where(sfreq_in > 0))
    if b[na-1] > indexmax:
        b[na-1] = indexmax
    #end if
    #Now you have start end end subscripts for the frequency bands
    sfedg = sfreq_in[a]
    sfedg = np.append(sfedg, sfreq_in[b[na-1]])
    #bin midpoints
    sfreq_out = 0.5*(sfreq_in[a]+sfreq_in[b])
    sdf_out = sfreq_in[b]-sfreq_in[a]
    #add these vars
    uvout.add_var('nspect', 'i')
    uvout['nspect'] = na
    uvout.add_var('sfreq', 'd')
    uvout['sfreq'] = sfreq_out
    uvout.add_var('sdf', 'd')
    uvout['sdf'] = sdf_out
    uvout.add_var('sfedg', 'd')
    uvout['sfedg'] = sfedg

    #spectral windows
    nschan = np.zeros(na, dtype=np.int)
    ischan = np.zeros(na, dtype=np.int)
    for j in range(na):
        nschan[j] = 1
        ischan[j] = j+1
    #endfor
    uvout.add_var('nschan', 'i')
    uvout['nschan'] = nschan
    uvout.add_var('ischan', 'i')
    uvout['ischan'] = ischan

    #define the record variables here
    uvout.add_var('ut', 'd')
    uvout.add_var('xtsys', 'r')
    uvout.add_var('ytsys', 'r')
    uvout.add_var('delay', 'd')
    uvout.add_var('pol', 'i')
    #version test
    if 'xtsys' in uv.vartable:
        version = "1.0"
    else:
        version = "2.0"
    #endelse

    # Loop over filenames, and add the other variables
    init = False
    init_pol = False
    pol = -71 #dummy
    ut = 0.0 #not necessarily the same ut variable as in the preamble
    xcount = 0
    utcount = 0
    for filename in ok_filelist:
        uv = aipy.miriad.UV(filename)

# generate FEM attn gain corrections
        fem_attn=rd_fem_attn(uv)
        fem_attn_out=fem_attn_update(fem_attn)
        fem_gain=10.**(-fem_attn_out/10.)
        timejd=Time(fem_attn['timestamp'].astype('int'),format='lv').jd
        if uv['source'] != src or uv['scanid'] != scanid:
            print 'Source name:',uv['source'],'is different from initial source name:',src
            print 'Or scanid:',uv['scanid'],'is different from initial source name:',scanid
            print 'Will stop processing files.'
            break
        #endif

        for preamble, data in uv.all():
# look up for gain correction in fem_gain
            uvw, t, (ant1, ant2) = preamble
            polstr=aipy.miriad.pol2str[uv['pol']]
            tidx=np.abs(timejd-t).argmin()-1
            if np.abs(timejd-t).min() < 1./24./3600. and (ant1 < 15) and (ant2 < 15):
                if polstr[0] == 'x':
                    p1=0
                if polstr[0] == 'y':
                    p1=1
                gain1=fem_gain[tidx,ant1,p1,0]*fem_gain[tidx,ant1,p1,1]
                if polstr[1] == 'x':
                    p2=0
                if polstr[1] == 'y':
                    p2=1
                gain2=fem_gain[tidx,ant2,p2,0]*fem_gain[tidx,ant2,p2,1]
                #### additional background should be considered before converting dB changes to actual gain changes 
                ### should really be (data - bkg) *= 1./((gain1*gain2)**0.5) + bkg
                data *= 1./((gain1*gain2)**0.5)
                #if tidx % 100 == 0 and ant1 == ant2:
                #    print 'gain change at t, ant1, ant2, pol:', tidx, ant1, ant2, polstr, (gain1*gain2)**0.5

            # Look for time change, correct for power and power square
            if preamble[1] != ut or init == False:
                if verbose:
                    print 'processed ',utcount,' times\r',
                ut = preamble[1]
                init = True
#power, power^2, and m
                if version == "1.0":
                    xts = uv['xtsys']
                    xts.shape = (nchan_in, nants)
                else:
                    xts = uv['xsampler'] #nfreq,nants,3; we only keep nfreq,nants
                    xts.shape = (nchan_in, nants, 3)
                    xts = xts[:,:,0] # the 1st element in the 3rd dimension is power
#power FEM attn gain correction for xx
                pxgain=fem_gain[tidx,:,0,0]*fem_gain[tidx,:,0,1]
                xts = xts/np.abs(pxgain)
                #### additional background should be considered if it is not small enough 
                # obtain the background (if there is a measurement)
                #if utcount == 75:
                #    xbkg=np.copy(xts)
                #if utcount < 75:
                #    xts[:,:15] *= 1./np.abs(pxgain)
                #else:
                #    xts[:,:15] = (xts[:,:15]-xbkg[:,:15])*1./np.abs(pxgain) + xbkg[:,:15]

#reshape and contract along axis
                xts.shape = (nch_avg, navg, nants)
                xts_flag = np.zeros_like(xts) #flag to account for nonzero values
                ok = np.where(xts != 0)
                xts_flag[ok] = 1.0
                xts_new = np.sum(xts, axis=1, dtype=np.float32)
                xts_flag_new = np.sum(xts_flag, axis=1, dtype=np.float32)
                ok = np.where(xts_flag_new != 0)
                xts_new[ok] = xts_new[ok]/xts_flag_new[ok]
#only keep the first na (# of spectral windows) values
                xts_new = xts_new[:na]
                xts_new.shape = (na*nants)
                uvout['ut'] = uv['ut']
                uvout['xtsys'] = xts_new
#ytsys
                if version == "1.0":
                    yts = uv['ytsys']
                    yts.shape = (nchan_in, nants)
                else:
                    yts = uv['ysampler'] #nfreq,nants,3; we only keep nfreq,nants
                    yts.shape = (nchan_in, nants, 3)
                    yts = yts[:,:,0] # the 1st element in the 3rd dimension is power

#power FEM attn gain correction for yy
                pygain=fem_gain[tidx,:,1,0]*fem_gain[tidx,:,1,1]
                yts = yts/np.abs(pygain)
                #### additional background should be considered if it is not small enough 
                # obtain the background (if there is a measurement)
                #if utcount == 75:
                #    ybkg=np.copy(yts)
                #if utcount < 75:
                #    yts[:,:15] *= 1./np.abs(pygain)
                #else:
                #    yts[:,:15] = (yts[:,:15]-ybkg[:,:15])*1./np.abs(pygain) + ybkg[:,:15]

#reshape and contract along axis
                yts.shape = (nch_avg, navg, nants)
                yts_flag = np.zeros_like(yts) #flag to account for nonzero values
                ok = np.where(yts != 0)
                yts_flag[ok] = 1.0
                yts_new = np.sum(yts, axis=1, dtype=np.float32)
                yts_flag_new = np.sum(yts_flag, axis=1, dtype=np.float32)
                ok = np.where(yts_flag_new != 0)
                yts_new[ok] = yts_new[ok]/yts_flag_new[ok]

                yts_new = yts_new[:na, :]
                yts_new.shape = (na*nants)
                uvout['ut'] = uv['ut']
                uvout['ytsys'] = yts_new
                uvout['delay'] = uv['delay']
                utcount = utcount+1


#output polarization if it has changed:
            if uv['pol'] != pol or init_pol == False:
                init_pol = True
                pol = uv['pol']
                uvout['pol'] = pol
                uvout['ut'] = uv['ut']

#auto- and cross-correlation data
            dflag = np.zeros(len(data))
            ok_data = np.where(np.absolute(data) > 0)
            dflag[ok_data] = 1.0
#shape and contract
            data.shape = (nch_avg, navg)
            dflag.shape = (nch_avg, navg)
            dataout = np.sum(data, axis=1, dtype=np.complex64)
            dflagout = np.sum(dflag, axis=1)
            okij = np.where(dflagout != 0)
            dataout[okij] = dataout[okij]/dflagout[okij]
#only keep first na values
            dataout = dataout[:na]
            uvout.write(preamble, dataout)
            xcount=xcount+1

    del(uv) #done
    #print xcount, utcount
    return ufilename
