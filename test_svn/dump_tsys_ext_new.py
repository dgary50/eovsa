# !/usr/bin/env python 
# Hacked from Dale's dump_tsys.py
# Needs an sdf variable
import time, os
import numpy as np
from util import Time
from astropy.io import fits

def strip_non_printable(string_in):
    ''' Only allow certain printable ascii characters to get through
    to fits write routine'''

#Keep values of 32 to 127
    stripped = (c for c in string_in if 31 < ord(c) < 127)
    return ''.join(stripped)
#End


def rd_miriad_tsys_file(filelist):
    ''' Read total power data (TSYS) directly from Miriad files
        Major change to standardize output to ut_mdj, sfreq, and order
        of indices of tsys as (npol, nant, nf, ntimes)
        Added an sdf output
        Added check for file existence, 2015-10-24, jmm
    '''
    import aipy

    # Find files corresponding to times

    if len(filelist) == 0:
        print 'rd_miriad_tsys_file: No files input'
        return None

    # Be sure that files exist, and has all of the appropriate elements
    filelist_test, ok_filelist, bad_filelist = valid_miriad_dataset(filelist)
    if len(ok_filelist) == 0:
        print 'rd_miriad_tsys_file: No valid files input'
        return None

    # Open first file and check that it has correct form
    uv = aipy.miriad.UV(ok_filelist[0])
    uvok = True
    if 'source' in uv.vartable:
        src = uv['source']
    else: uvok = False
    if 'scanid' in uv.vartable:
        scanid = uv['scanid']
    else: 
        scanid = 'Unknown'
    if 'proj' in uv.vartable:
        proj = uv['proj']
    else: 
        proj = 'Unknown'
    if 'sfreq' in uv.vartable:
        sfreq = uv['sfreq']
        if len(sfreq) == 0:
            uvok = False
    else: uvok = False
    if 'sdf' in uv.vartable:
        sdf = uv['sdf']
        if len(sdf) == 0:
            uvok = False
    else: uvok = False
    if 'nants' in uv.vartable:
        nants = uv['nants']
    else: 
        nants = 16
        uvok = False
    if 'antlist' in uv.vartable:
        antennalist = uv['antlist']
    else:
        antennalist = ' '.join(map(str, range(1, nants)))
    if 'ut' in uv.vartable:
        pass
    else: uvok = False
    if 'xtsys' in uv.vartable or 'xsampler' in uv.vartable:
        if 'xtsys' in uv.vartable:
            version = "1.0"
        if 'xsampler' in uv.vartable:
            version = "2.0"
    else: uvok = False
    if 'ytsys' in uv.vartable or 'ysampler' in uv.vartable:
        pass
    else: uvok = False
    if 'lst' in uv.vartable:
        version = "3.0"
    if not uvok:
        print 'Miriad file has bad format'
        return None

    if version == "3.0":
        lst = []
    #endif
    utd = []
    xtsys = []
    ytsys = []
    # Loop over filenames
    for filename in ok_filelist:
        uv = aipy.miriad.UV(filename)
        if uv['source'] != src:
            print 'Source name:',uv['source'],'is different from initial source name:',src
            print 'Will stop reading files.'
            break
        #endif
        # uv.select('antennae',0,1,include=True)
        # Read first record of data
        preamble, data = uv.read()
        ut = preamble[1]
        ut1 = ut - 2400000.5
        utd.append(ut1)
        if version == "1.0":
            xtsys.append(uv['xtsys'])
            ytsys.append(uv['ytsys'])
        else:
            xtsys.append(uv['xsampler'])
            ytsys.append(uv['ysampler'])
        #endif
        if version == "3.0":
            lst.append(uv['lst'])
        #endif
        for preamble, data in uv.all():
            # Look for time change
            if preamble[1] != ut:
                # Time has changed, so read new xtsys and ytsys
                ut = preamble[1]
                if version == "1.0":
                    xtsys.append(uv['xtsys'])
                    ytsys.append(uv['ytsys'])
                else:
                    xtsys.append(uv['xsampler'])
                    ytsys.append(uv['ysampler'])
                #endif
                ut1 = ut - 2400000.5
                utd.append(ut1)
                if version == "3.0":
                    lst.append(uv['lst'])
                #endif
            #endif
        #endfor
    #endfor

    utd = np.array(utd)
    xtsys = np.array(xtsys)
    ytsys = np.array(ytsys)
#    print len(xtsys)
    if version == "3.0":
        lst = np.array(lst)
    #endif

#    print version
    if version == "1.0":
        xtsys.shape = (len(utd),len(sfreq),nants)
        ytsys.shape = (len(utd),len(sfreq),nants)
    else:
        xtsys.shape = (len(utd),len(sfreq),nants,3)
        ytsys.shape = (len(utd),len(sfreq),nants,3)
    #endif
    if version == "1.0":
        tsys = np.array((xtsys,ytsys))    #order is npol, nt, nf, nants
        tsys = np.swapaxes(tsys,1,3)  # Order is now npol, nants, nf, nt
        tsys = np.swapaxes(tsys,0,1)  # Order is now nants, npol, nf, nt, as desired
        good, = np.where(tsys.sum(0).sum(0).sum(1) != 0.0)
        if len(good) > 0:
            tsys = tsys[:,:,good,:]
            sfreq = sfreq[good]
            sdf = sdf[good]
            return {'source':src, 'scanid':scanid, 'proj':proj, 'sfreq':sfreq, 'sdf':sdf, 'ut_mjd':utd, 'tsys':tsys, 'file0':ok_filelist[0], 'antennalist':antennalist, 'nants':nants, 'version':version}
        else:
            print 'RD_MIRIAD_TSYS_FILE: No Good DATA: '
            return None
        #endelse
    else:
        tsys = np.array((xtsys,ytsys))    #order is npol, nt, nf, nants, 3
        tpwr = tsys[:,:,:,:,0]        #order is npol, nt, nf, nants
        tpwr = np.swapaxes(tpwr,1,3)  # Order is now npol, nants, nf, nt
        tpwr = np.swapaxes(tpwr,0,1)  # Order is now nants, npol, nf, nt, as desired
        tpwr2 = tsys[:,:,:,:,1]
        tpwr2 = np.swapaxes(tpwr2,1,3)
        tpwr2 = np.swapaxes(tpwr2,0,1)
        nsamp = tsys[:,:,:,:,2]
        nsamp = np.swapaxes(nsamp,1,3)
        nsamp = np.swapaxes(nsamp,0,1)
        good, = np.where(tpwr.sum(0).sum(0).sum(1) != 0.0)
        if len(good) > 0:
            tpwr = tpwr[:,:,good,:]
            tpwr2 = tpwr2[:,:,good,:]
            nsamp = nsamp[:,:,good,:]
            sfreq = sfreq[good]
            sdf = sdf[good]
            if version == "3.0":
                return {'source':src, 'scanid':scanid, 'proj':proj, 'sfreq':sfreq, 'sdf':sdf, 'ut_mjd':utd, 'tpwr':tpwr, 'tpwr2':tpwr2, 'nsamp':nsamp, 'file0':ok_filelist[0], 'antennalist':antennalist, 'nants':nants, 'version':version, 'lst':lst}
            else:
                return {'source':src, 'scanid':scanid, 'proj':proj, 'sfreq':sfreq, 'sdf':sdf, 'ut_mjd':utd, 'tpwr':tpwr, 'tpwr2':tpwr2, 'nsamp':nsamp, 'file0':ok_filelist[0], 'antennalist':antennalist, 'nants':nants, 'version':version}
            #endelse
        else:
            print 'RD_MIRIAD_TSYS_FILE: No Good DATA: '
            return None
        #endelse
    #endif
#end of rd_miriad_tsys_file

def tsys_writetofits(xdat, calflag):
    ''' This takes the dictionary xdat of rd_miriad_tsys_file and creates
    a FITS file. Unlike the UDB fits files this puts the tsys data as the
    primary output because I cannot figure out how to get a
    multidimensional array into a column... '''
    file_out = ''

    if xdat == None or len(xdat) == 0:
        print 'tsys_writetofits: No data input'
        return file_out
#   end

#IDB fits files
    idbfitsdir = '/data1/eovsa/fits/fullres/'

#create a filename
    file0 = xdat['file0'].split("/")
    print file0
    file0 = file0[3]
    yr = file0[3:7]
    yr2 = file0[5:7]
    mn = file0[7:9]
    dy = file0[9:11]
    hh = file0[11:13]
    mm = file0[13:15]
    ss = file0[15:]
    file_out = 'eovsa_1-18GHz_sp_fullres_'+yr+mn+dy+'_'+hh+mm+ss+'.fts'
#add directory
    outdir = idbfitsdir+yr+mn+dy
    if os.path.isdir(outdir) == False:
        print "tsys_writetofits: creating "+outdir
        os.mkdir(outdir)
    #end if
    file_out = outdir+'/'+file_out

    date_obs = yr+'-'+mn+'-'+dy+'T'+hh+':'+mm+':'+ss+'.000'
#Convert to unix time, and add seconds to get date_end
#http://astropy.readthedocs.org/en/latest/time
    ut = xdat['ut_mjd']
    t0 = Time(date_obs)
    dt = int(86400*(max(ut)-min(ut)))
    t1 = Time(t0.unix+dt, format = 'unix')
    date_end = t1.isot

# Create the primary header
    version = xdat['version']
    if version == "1.0":
        tsys = xdat['tsys']
    else:
        tpwr = xdat['tpwr']
    #endif
    hdu = fits.PrimaryHDU(tpwr)

# Set up the extensions: sfreq, sdf, ut
    sfreq = xdat['sfreq']
    col1 = fits.Column(name='sfreq', format='E', array = sfreq)
    cols1 = fits.ColDefs([col1])
    tbhdu1 = fits.BinTableHDU.from_columns(cols1)
    tbhdu1.name = 'SFREQ'
    sdf = xdat['sdf']
    col2 = fits.Column(name='sdf', format='E', array = sdf)
    cols2 = fits.ColDefs([col2])
    tbhdu2 = fits.BinTableHDU.from_columns(cols2)
    tbhdu2.name = 'SDF'
#    ut = xdat['ut']
# Split up mjd into days and msec, really intuitive syntax, thanks python
    ut_int = ut.astype(np.int32)
    ut_msec = 1000.0*86400.0*(ut-ut_int)
    ut_ms1 = ut_msec.astype(np.int32)

# J is the format code for a 32 bit integer, who would have thought
# http://astropy.readthedocs.org/en/latest/io/fits/usage/table.html
    col3 = fits.Column(name='mjd', format='J', array = ut_int)
    col4 = fits.Column(name='time', format='J', array = ut_ms1)

    cols3 = fits.ColDefs([col3, col4])
    tbhdu3 = fits.BinTableHDU.from_columns(cols3)
    tbhdu3.name = 'UT'

#create an HDUList object to put in header information
    hdulist = fits.HDUList([hdu, tbhdu1, tbhdu2, tbhdu3])

#primary header
    prihdr = hdulist[0].header
#Header information
    temp_out = file_out.split("/")
    temp_out = temp_out[len(temp_out)-1]
    prihdr.set('FILENAME', temp_out)
    prihdr.set('ORIGIN', 'NJIT', 'Institute where file was written')
    prihdr.set('TELESCOP', 'EOVSA', 'Expanded Owens Valley Solar Array')
    prihdr.set('OBJ_ID', xdat['source'], 'Object ID')
    prihdr.set('SCAN_ID', xdat['scanid'], 'Scan ID for this dataset')
    print xdat['proj']
    proj_tmp = xdat['proj']
    proj_tmp = strip_non_printable(proj_tmp)
    print proj_tmp
    prihdr.set('PROJECT_', proj_tmp, 'EOVSA Project ID')
    prihdr.set('ID', int(yr2+mn+dy+hh+mm+ss), 'Catalog ID, yymmddhhmm')
    prihdr.set('TYPE', 1, 'Spectrum')
    prihdr.set('DATE_OBS', date_obs, 'Start date/time of observation')
    prihdr.set('DATE_END', date_end, 'End date/time of observation')
    prihdr.set('FREQMIN', min(sfreq), 'Min freq in observation (GHz)')
    prihdr.set('FREQMAX', max(sfreq), 'Max freq in observation (GHz)')
    prihdr.set('XCEN', 0.0, 'Antenna pointing in arcsec from Sun centre')
    prihdr.set('YCEN', 0.0, 'Antenna pointing in arcsec from Sun centre')
    prihdr.set('POLARIZA', 'XX, YY', 'Polarizations present')
    prihdr.set('RESOLUTI', 0.0, 'Resolution value')
    prihdr.set('NANTS', xdat['nants'], 'Number of Antennae')
    prihdr.set('ANTENNA', xdat['antennalist'], 'Used antennae')
    prihdr.set('VERSION', version, 'SW version')
    if(calflag == True):
        prihdr.set('CAL_FLAG', 1, 'Calibration Flag: 1 for calibrated, 0 for not')
    else:
        prihdr.set('CAL_FLAG', 0, 'Calibration Flag: 1 for calibrated, 0 for not')
    #endif
# Write the file
    hdulist.writeto(file_out, clobber=True)

    return file_out
#END of tsys_writetofits

def tsys_writeudbfits(xdat, calflag):
    ''' This takes the dictionary xdat of rd_miriad_tsys_file and
    creates a UDB FITS file. Unlike the original UDB fits files this
    puts the tsys data as the primary output because I cannot figure
    out how to get a multidimensional array into a column...'''
    file_out = ''

    if xdat == None or len(xdat) == 0:
        print 'tsys_writeudbfits: No data input'
        return file_out
#   end

#UDB fits files
    udbfitsdir = '/data1/eovsa/fits/'

#create a filename
    file0 = xdat['file0'].split("/")
    lf0 = len(file0)
    file0 = file0[lf0-1]
    print file0
    yr = file0[3:7]
    yr2 = file0[5:7]
    mn = file0[7:9]
    dy = file0[9:11]
    hh = file0[11:13]
    mm = file0[13:15]
    ss = file0[15:]
    file_out = 'eovsa_1-18GHz_sp_'+yr+mn+dy+'_'+hh+mm+ss+'.fts'
#add directory
    outdir = udbfitsdir+yr+mn+dy
    if os.path.isdir(outdir) == False:
        print "tsys_writeudbfits: creating "+outdir
        os.mkdir(outdir)
    #end if
    file_out = outdir+'/'+file_out

    date_obs = yr+'-'+mn+'-'+dy+'T'+hh+':'+mm+':'+ss+'.000'
#Convert to unix time, and add seconds to get date_end
#http://astropy.readthedocs.org/en/latest/time
    ut = xdat['ut_mjd']
    print "date_obs: ", date_obs
    t0 = Time(date_obs)
    dt = int(86400*(max(ut)-min(ut)))
    t1 = Time(t0.unix+dt, format = 'unix')
    date_end = t1.isot

# Create the primary header
    version = xdat['version']
    if version == "1.0":
        tpwr = xdat['tsys']
    else:
        tpwr = xdat['tpwr']
    #endif
    hdu = fits.PrimaryHDU(tpwr)

# Set up the extensions: sfreq, sdf, ut
    sfreq = xdat['sfreq']
    col1 = fits.Column(name='sfreq', format='E', array = sfreq)
    cols1 = fits.ColDefs([col1])
    tbhdu1 = fits.BinTableHDU.from_columns(cols1)
#    tbhdu1.update_ext_name('SFREQ')
    tbhdu1.name = 'SFREQ'
    sdf = xdat['sdf']
    col2 = fits.Column(name='sdf', format='E', array = sdf)
    cols2 = fits.ColDefs([col2])
    tbhdu2 = fits.BinTableHDU.from_columns(cols2)
    tbhdu2.name = 'SDF'
#    ut = xdat['ut']
# Split up mjd into days and msec, really intuitive syntax, thanks python
    ut_int = ut.astype(np.int32)
    ut_msec = 1000.0*86400.0*(ut-ut_int)
    ut_ms1 = ut_msec.astype(np.int32)

# J is the format code for a 32 bit integer, who would have thought
# http://astropy.readthedocs.org/en/latest/io/fits/usage/table.html
    col3 = fits.Column(name='mjd', format='J', array = ut_int)
    col4 = fits.Column(name='time', format='J', array = ut_ms1)

    cols3 = fits.ColDefs([col3, col4])
    tbhdu3 = fits.BinTableHDU.from_columns(cols3)
    tbhdu3.name = 'UT'

#create an HDUList object to put in header information
    hdulist = fits.HDUList([hdu, tbhdu1, tbhdu2, tbhdu3])

#primary header
    prihdr = hdulist[0].header
#Header information, strip last character from strings
    obj_id = xdat['source']
    obj_id = obj_id[:len(obj_id)-1]
    scan_id = xdat['scanid']
    scan_id = scan_id[:len(scan_id)-1]
    proj_id = xdat['proj']
    proj_id = proj_id[:len(proj_id)-1]
    proj_id = strip_non_printable(proj_id)
    ant_list = xdat['antennalist']
    ant_list = ant_list[:len(ant_list)-1]

    temp_out = file_out.split("/")
    temp_out = temp_out[len(temp_out)-1]
    prihdr.set('FILENAME', temp_out)
    prihdr.set('ORIGIN', 'NJIT', 'Institute where file was written')
    prihdr.set('TELESCOP', 'EOVSA', 'Expanded Owens Valley Solar Array')
    prihdr.set('OBJ_ID', obj_id, 'Object ID')
    prihdr.set('SCAN_ID', scan_id, 'Scan ID for this dataset')
    prihdr.set('PROJECT_', proj_id, 'EOVSA Project ID')
    prihdr.set('ID', int(yr2+mn+dy+hh+mm+ss), 'Catalog ID, yymmddhhmm')
    prihdr.set('TYPE', 1, 'Spectrum')
    prihdr.set('DATE_OBS', date_obs, 'Start date/time of observation')
    prihdr.set('DATE_END', date_end, 'End date/time of observation')
    prihdr.set('FREQMIN', min(sfreq), 'Min freq in observation (GHz)')
    prihdr.set('FREQMAX', max(sfreq), 'Max freq in observation (GHz)')
    prihdr.set('XCEN', 0.0, 'Antenna pointing in arcsec from Sun centre')
    prihdr.set('YCEN', 0.0, 'Antenna pointing in arcsec from Sun centre')
    prihdr.set('POLARIZA', 'XX, YY', 'Polarizations present')
    prihdr.set('RESOLUTI', 0.0, 'Resolution value')
    prihdr.set('NANTS', xdat['nants'], 'Number of Antennae')
    prihdr.set('ANTENNA', ant_list, 'Used antennae')
    prihdr.set('VERSION', version, 'SW version')
    if(calflag == True):
        prihdr.set('CAL_FLAG', 1, 'Calibration Flag: 1 for calibrated, 0 for not')
    else:
        prihdr.set('CAL_FLAG', 0, 'Calibration Flag: 1 for calibrated, 0 for not')
    #endif
# Write the file
    hdulist.writeto(file_out, clobber=True)

    return file_out
#END of tsys_writeudbfits

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
        if os.path.isdir(filename) == False or os.path.isfile(filename+'/flags') == False or os.path.isfile(filename+'/header') == False or os.path.isfile(filename+'/vartable') == False or os.path.isfile(filename+'/visdata') == False:
            tempvar = False
        #end if
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

def udbfile_create(filelist, ufilename):
    '''Given a list of IDB filenames, create the appropriate UDB file, by
    averaging over energy bands, but keep 1 second time resolution'''

    import aipy
    if len(filelist) == 0:
        print 'udbfile_create: No files input'
        return []
    #endif

    # Be sure that files exist, and has all of the appropriate elements
    filelist_test, ok_filelist, bad_filelist = valid_miriad_dataset(filelist)
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
    navg = 10
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
    nschan = np.ones(na, dtype=np.int)
    ischan = np.arange(1, na, dtype=np.int)
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
    #version test
    if 'xtsys' in uv.vartable:
        version = "1.0"
    else:
        version = "2.0"
    #endelse
    if 'lst' in uv.vartable:
        version = "3.0"
    #endif
    # Loop over filenames, and add the other variables
    init = False
    init_pol = False
    pol = -71 #dummy
    ut = 0.0 #not necessarily the same ut variable as in the preamble
    xcount = 0
    utcount = 0
    for filename in ok_filelist:
        uv = aipy.miriad.UV(filename)
        if uv['source'] != src or uv['scanid'] != scanid:
            print 'Source name:',uv['source'],'is different from initial source name:',src
            print 'Or scanid:',uv['scanid'],'is different from initial source name:',scanid
            print 'Will stop processing files.'
            break
        #endif

        for preamble, data in uv.all():
            # Look for time change
            if preamble[1] != ut or init == False:
                ut = preamble[1]
                init = True
#xtsys
                if version == "1.0":
                    xts = uv['xtsys']
                    xts.shape = (nchan_in, nants)
                else:
                    xts = uv['xsampler'] #nfreq,nants,3; we only keep nfreq,nants
                    xts.shape = (nchan_in, nants, 3)
                    xts = xts[:,:,0] #xts is now nfreq, nants
                #endelse
#reshape and contract along axis
                xts.shape = (nch_avg, navg, nants)
                xts_flag = np.zeros_like(xts) #flag to account for nonzero values
                ok = np.where(xts != 0)
                xts_flag[ok] = 1.0
                xts_new = np.sum(xts, axis=1, dtype=np.float32)
                xts_flag_new = np.sum(xts_flag, axis=1, dtype=np.float32)
                ok = np.where(xts_flag_new != 0)
                xts_new[ok] = xts_new[ok]/xts_flag_new[ok]
#only keep the first na values
                xts_new = xts_new[:na, :]
                xts_new.shape = (na*nants)
                uvout['ut'] = uv['ut']
                if version == "3.0":
                    uvout['lst'] = uv['lst']
                else:
                    uvout['lst'] = 0.0
                #endelse
                uvout['xtsys'] = xts_new
#ytsys
                if version == "1.0":
                    xts = uv['ytsys']
                    xts.shape = (nchan_in, nants)
                else:
                    xts = uv['ysampler'] #nfreq,nants,3; we only keep nfreq,nants
                    xts.shape = (nchan_in, nants, 3)
                    xts = xts[:,:,0] #xts is now nfreq, nants
                #endelse
#reshape and contract along axis
                xts.shape = (nch_avg, navg, nants)
                xts_flag = np.zeros_like(xts) #flag to account for nonzero values
                ok = np.where(xts != 0)
                xts_flag[ok] = 1.0
                xts_new = np.sum(xts, axis=1, dtype=np.float32)
                xts_flag_new = np.sum(xts_flag, axis=1, dtype=np.float32)
                ok = np.where(xts_flag_new != 0)
                xts_new[ok] = xts_new[ok]/xts_flag_new[ok]

                xts_new = xts_new[:na, :]
                xts_new.shape = (na*nants)
                uvout['ut'] = uv['ut']
                if version == "3.0":
                    uvout['lst'] = uv['lst']
                else:
                    uvout['lst'] = 0.0
                #endelse
                uvout['ytsys'] = xts_new
                uvout['delay'] = uv['delay']
                utcount = utcount+1
#                print utcount
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
        #endfor
    #endfor

    del(uv) #done
    print xcount, utcount
    return ufilename

#End of udbfile_create

