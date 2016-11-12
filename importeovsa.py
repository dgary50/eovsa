
''' Main routine to convert the EOVSA IDB files to CASA measurement set.
'''
import os
import sys
import gc
import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as constants
import time
import glob
import casacore.tables.table as tb
import aipy
import chan_util_bc as cu
import read_idb as ri
from dbutil import Time

def bl_list2(nant=16):
    ''' Returns a two-dimensional array bl2ord that will translate
        a pair of antenna indexes (antenna number - 1) to the ordinal
        number of the baseline in the 'x' key.  Note bl2ord(i,j) = bl2ord(j,i),
        and bl2ord(i,i) = -1.
    '''
    bl2ord = np.ones((nant,nant),dtype='int')*(-1)
    k = 0
    for i in range(nant):
        for j in range(i,nant):
            bl2ord[i,j] = k
            # bl2ord[j,i] = k
            k+=1
    return bl2ord

def get_band_edge(nband=None):
    # Input the frequencies from UV, returen the indices frequency edges of all bands    
    if not nband:
        nband=34
    idx_start_freq = [0]
    ntmp=0
    for i in range(1,nband+1):
        ntmp+=len(cu.start_freq(i))
        idx_start_freq.append(ntmp)
    return np.asarray(idx_start_freq)

def idb2ms(trange, modelms=None, outpath=None, nocreatms=False, nowritems=False):
    ''' This is the main routine to convert the EOVSA IDB files to CASA measurement set.
        Parameters
        ----------
        trange:  Finding the IDB files within the given time range. If trange is not a Time() object,
        assume that it is the list of files to read.    
    '''



    if type(trange) == Time:
        filelist = ri.get_trange_files(trange)
    else:
        # If input type is not Time, assume that it is the list of files to read
        filelist = trange

    if not modelms:
        modelms='/home/user/sjyu/20160531/ms/sun/SUN/SUN_20160531T142234-10m.1s.ms'
    try:
        for f in filelist:
            os.path.exists(f)
    except ValueError:
        print("Some files in filelist are invalid. Aborting...")
    if not outpath:
        # use current directory
        outpath='./'

    nocreatms = True
    nowritems = False

    # if len(filelist) == 0:
    #     print 'eovsa_idb2ms: No files input'
    #     return None

    # # Be sure that files exist, and has all of the appropriate elements
    # filelist_test, ok_filelist, bad_filelist = valid_miriad_dataset(filelist)
    # if len(ok_filelist) == 0:
    #     print 'eovsa_idb2ms: No valid files input'
    #     return None

    for filename in filelist:     
        uv = aipy.miriad.UV(filename)
        #if uv['source'].lower() == 'sun':
        #    outpath=outpath+'sun/'
        #    if not os.path.exists(outpath):
        #        os.mkdir(outpath)
        #else:
        #    outpath=outpath+'calibrator/'    
        #    if not os.path.exists(outpath):
        #        os.mkdir(outpath)
        #uv.rewind()

        start_time = 0 # The start and stop times are referenced to ref_time_jd in second
        end_time=600
        delta_time = 1
        time_steps = (end_time-start_time)/delta_time
        time0 = time.time()

        if 'antlist' in uv.vartable:
            ants = uv['antlist']
            antlist = map(int, ants.split())
        else:
            antlist = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]    

        good_idx = np.arange(len(uv['sfreq']))

        ref_time_jd = uv['time']
        ref_time_mjd = (ref_time_jd -  2400000.5)*24.*3600.+0.5*delta_time
        nf = len(good_idx)
        freq = uv['sfreq'][good_idx]
        npol = uv['npol']
        nants = uv['nants']
        sdf=uv['sdf']
        project= uv['proj']
        source_id = uv['source']
        ra,dec = uv['ra'], uv['dec']
        scan_id = uv['scanid']
        nbl = nants*(nants-1)/2
        bl2ord = bl_list2(nants)
        npairs = nbl+nants
        flag = np.ones((npol,nf,time_steps,npairs), dtype=bool)
        out = np.zeros((npol,nf,time_steps,npairs),dtype=np.complex64)  # Cross-correlations
        uvwarray = np.zeros((3,time_steps,npairs),dtype=np.float)
        bandedge = get_band_edge()
        nband = len(bandedge)-1

        uv.rewind()
        l = -1
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
            l += 1
            out[k,:,l/(npairs*npol),bl2ord[i0,j0]] = data.data
            flag[k,:,l/(npairs*npol),bl2ord[i0,j0]] = data.mask
            if i!=j:
                if k == 3: 
                    uvwarray[:,l/(npairs*npol),bl2ord[i0,j0]] = uvw*constants.speed_of_light/1e9

        nrows = time_steps*npairs
        out=out.reshape(npol,nf,nrows)
        flag=flag.reshape(npol,nf,nrows)
        uvwarray = uvwarray.reshape(3,nrows).swapaxes(0,1)
        uvwarray = np.tile(uvwarray,(nband,1))
        sigma = np.ones((nrows,4),dtype=np.float)
        sigma = np.tile(sigma,(nband,1))

        print 'IDB File {0} is readed in --- {1:10.2f} seconds ---'.format(filename,(time.time() - time0))


        # nocreatms=True
        # nowritems=True
        # nocreatms=False
        # nowritems=False
        # msname = list(filename.split('/')[-1])
        # msname.insert(11,'T')
        # msname = outpath+source_id.upper()+'_'+''.join(msname[3:])+'-10m.1s.ms'
        msname=outpath+filename.split('/')[-1]+'.ms'



        if not nocreatms:
            print 'Empty MS {0} created in --- {1:10.2f} seconds ---'.format(msname,(time.time() - time0))
        else:
            os.system("rm -fr %s"%msname)
            os.system("cp -r "+" %s"%modelms+" %s"%msname)     
            print 'Standard MS is copied to {0} in --- {1:10.2f} seconds ---'.format(msname,(time.time() - time0))        



        if not nowritems:
            print '----------------------------------------'
            print "Updating the main table of" '%s'%msname
            print '----------------------------------------'
            tabl=tb(msname,readonly=False)
            for l,bdedge in enumerate(bandedge[:-1]):
                time1 = time.time()
                nchannels = (bandedge[l+1]-bandedge[l])
                for row in range(nrows):
                    tabl.putcell('DATA',(row+l*nrows),out[:,bandedge[l]:bandedge[l+1],row].swapaxes(0,1))
                    tabl.putcell('FLAG',(row+l*nrows),flag[:,bandedge[l]:bandedge[l+1],row].swapaxes(0,1))
                print '---spw {0:02d} is updated in --- {1:10.2f} seconds ---'.format((l+1),time.time() - time1)
            tabl.putcol('UVW',uvwarray)
            tabl.putcol('SIGMA',sigma) 
            tabl.putcol('WEIGHT',1.0/sigma**2) 
            timearr = np.arange((time_steps),dtype=np.float)
            timearr=timearr.reshape(1,time_steps,1)
            timearr=np.tile(timearr,(nband,1,npairs))
            timearr=timearr.reshape(nband*npairs*time_steps)+ref_time_mjd
            tabl.putcol('TIME',timearr)
            tabl.putcol('TIME_CENTROID',timearr)
            tabl.close()

            print '----------------------------------------'
            print "Updating the OBSERVATION table of" '%s'%msname
            print '----------------------------------------'
            tabl=tb(msname+'/OBSERVATION',readonly=False)
            tabl.putcol('TIME_RANGE',np.asarray([ref_time_mjd-0.5*delta_time,ref_time_mjd+end_time-0.5*delta_time]).reshape(2,1).swapaxes(0,1))
            tabl.putcol('OBSERVER',['EOVSA team'])
            tabl.close()

            print '----------------------------------------'
            print "Updating the POINTING table of" '%s'%msname
            print '----------------------------------------'
            tabl=tb(msname+'/POINTING',readonly=False)
            timearr=np.arange((time_steps),dtype=np.float).reshape(1,time_steps,1)
            timearr=np.tile(timearr,(nband,1,nants))
            timearr=timearr.reshape(nband*time_steps*nants)+ref_time_mjd        
            tabl.putcol('TIME',timearr)
            tabl.putcol('TIME_ORIGIN',timearr-0.5*delta_time)
            direction=tabl.getcol('DIRECTION')
            direction[:,0,0]=ra
            direction[:,0,1]=dec
            tabl.putcol('DIRECTION',direction)
            target=tabl.getcol('TARGET')
            target[:,0,0]=ra
            target[:,0,1]=dec        
            tabl.putcol('TARGET',target)
            tabl.close()

            print '----------------------------------------'
            print "Updating the SOURCE table of" '%s'%msname
            print '----------------------------------------'
            tabl=tb(msname+'/SOURCE',readonly=False)
            radec=tabl.getcol('DIRECTION')
            radecshape=radec.shape
            radecflatten=radec.flatten()
            radecflatten[0],radecflatten[1]=ra,dec
            tabl.putcol('DIRECTION',radecflatten.reshape(radecshape))        
            name=[source_id]
            tabl.putcol('NAME',name)
            tabl.close()  


            print '----------------------------------------'
            print "Updating the DATA_DESCRIPTION table of" '%s'%msname
            print '----------------------------------------'
            tabl=tb(msname+'/DATA_DESCRIPTION/',readonly=False)  
            pol_id=tabl.getcol('POLARIZATION_ID')
            pol_id*=0
            tabl.putcol('POLARIZATION_ID',pol_id)
            tabl.close()

            if not nocreatms:
                print '----------------------------------------'
                print "Updating the POLARIZATION table of" '%s'%msname
                print '----------------------------------------'        
                tabl=tb(msname+'/POLARIZATION/',readonly=False)  
                tabl.removerows(rownrs=np.arange(1,nband,dtype=int))
                tabl.close()         

            print '----------------------------------------'
            print "Updating the FIELD table of" '%s'%msname
            print '----------------------------------------'        
            tabl=tb(msname+'/FIELD/',readonly=False) 
            radec=tabl.getcol('DELAY_DIR')
            radecshape=radec.shape
            radecflatten=radec.flatten()
            radecflatten[0],radecflatten[1]=ra,dec
            tabl.putcol('DELAY_DIR',radecflatten.reshape(radecshape))

            radec=tabl.getcol('PHASE_DIR')
            radecshape=radec.shape
            radecflatten=radec.flatten()
            radecflatten[0],radecflatten[1]=ra,dec
            tabl.putcol('PHASE_DIR',radecflatten.reshape(radecshape))

            radec=tabl.getcol('REFERENCE_DIR')
            radecshape=radec.shape
            radecflatten=radec.flatten()
            radecflatten[0],radecflatten[1]=ra,dec
            tabl.putcol('REFERENCE_DIR',radecflatten.reshape(radecshape))

            name=np.array([source_id],dtype='|S{0}'.format(len(source_id)+1))
            tabl.putcol('NAME',name)     
            tabl.close()   

        print 'finished in --- {0:10.2f} seconds ---'.format(time.time() - time0)      


            
