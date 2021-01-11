#!/usr/bin/env python 
# Hacked from Dale's dump_tsys.py
# Needs an sdf variable
import time, os
import numpy as np
from util import Time
from astropy.io import fits

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
    filelist_test, ok_filelist = valid_miriad_dataset(filelist)
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
    else: uvok = False
    if 'sdf' in uv.vartable:
        sdf = uv['sdf']
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
    if not uvok:
        print 'Miriad file has bad format'
        return None

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
            #endif
        #endfor
    #endfor

    utd = np.array(utd)
    xtsys = np.array(xtsys)
    ytsys = np.array(ytsys)
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
        tsys = tsys[:,:,good,:]
        sfreq = sfreq[good]
        sdf = sdf[good]
        return {'source':src, 'scanid':scanid, 'proj':proj, 'sfreq':sfreq, 'sdf':sdf, 'ut_mjd':utd, 'tsys':tsys, 'file0':filelist[0], 'antennalist':antennalist, 'nants':nants}
    else:
        tsys = np.array((xtsys,ytsys))# order is npol, nt, nf, nants, 3
        tsys = np.swapaxes(tsys,1,3)  # Order is now npol, nants, nf, nt, 3
        tsys = np.swapaxes(tsys,0,1)  # Order is now nants, npol, nf, nt, 3
        good, = np.where(tsys.sum(0).sum(0).sum(1) != 0.0)
        tsys = tsys[:,:,good,:,:]
        sfreq = sfreq[good]
        sdf = sdf[good]
        return {'source':src, 'scanid':scanid, 'proj':proj, 'sfreq':sfreq, 'sdf':sdf, 'ut_mjd':utd, 'tpwr':tsys, 'file0':filelist[0], 'antennalist':antennalist, 'nants':nants}
    #endif
#end of read_miriad_tsys_tile

def tsys_writetofits(xdat, calflag):
    ''' This takes the dictionary xdat of rd_miriad_tsys_file and creates
    a FITS file. Unlike the UDB fits files this puts the tsys data as the
    primary output because I cannot figure out how to get a
    multidimensional array into a column... '''
    file_out = ''

    if xdat == None or len(xdat) == 0:
        print 'tsys_writetofits: No data input'
        return file_out
    #endif

    if 'tsys' in xdat:
        version = "1.0"
    else:
        version = "2.0"
    #endif

#IDB fits files
    idbfitsdir = '/data1/eovsa/fits/fullres/'

#create a filename
    file0 = xdat['file0'].split("/")
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
    if version == "1.0":
        tsys = xdat['tsys']
        hdu = fits.PrimaryHDU(tsys)
    else:
        tpwr = xdat['tpwr'][:,:,:,:,0]
        hdu = fits.PrimaryHDU(tpwr)
    #endif

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
    prihdr.set('PROJECT_', xdat['proj'], 'EOVSA Project ID')
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
    prihdr.set('NANTS', len(tsys), 'Number of Antennae')
    prihdr.set('ANTENNA', xdat['antennalist'], 'Used antennae (e.g. 1 2 3 9 5 6 14 8)')
    prihdr.set('VERSION', version, 'Output file version')
    if(calflag == True):
        prihdr.set('CAL_FLAG', 1, 'Calibration Flag: 1 for calibrated, 0 for uncali')
    else:
        prihdr.set('CAL_FLAG', 0, 'Calibration Flag: 1 for calibrated, 0 for uncali')
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
    #endif

    if 'tsys' in xdat:
        version = "1.0"
    else:
        version = "2.0"
    #endif

#IDB fits files
    udbfitsdir = '/data1/eovsa/fits/'

#create a filename
    file0 = xdat['file0'].split("/")
    file0 = file0[3]
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
    t0 = Time(date_obs)
    dt = int(86400*(max(ut)-min(ut)))
    t1 = Time(t0.unix+dt, format = 'unix')
    date_end = t1.isot

# Create the primary header
    if version == "1.0":
        tsys = xdat['tsys']
        hdu = fits.PrimaryHDU(tsys)
    else:
        tpwr = xdat['tpwr'][:,:,:,:,0]
        hdu = fits.PrimaryHDU(tpwr)
    #endif

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
#Header information
    temp_out = file_out.split("/")
    temp_out = temp_out[len(temp_out)-1]
    prihdr.set('FILENAME', temp_out)
    prihdr.set('ORIGIN', 'NJIT', 'Institute where file was written')
    prihdr.set('TELESCOP', 'EOVSA', 'Expanded Owens Valley Solar Array')
    prihdr.set('OBJ_ID', xdat['source'], 'Object ID')
    prihdr.set('SCAN_ID', xdat['scanid'], 'Scan ID for this dataset')
    prihdr.set('PROJECT_', xdat['proj'], 'EOVSA Project ID')
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
    prihdr.set('ANTENNA', xdat['antennalist'], 'Used antennae (e.g. 1 2 3 9 5 6 14 8)')
    prihdr.set('VERSION', version, 'Output file version')
    if(calflag == True):
        prihdr.set('CAL_FLAG', 1, 'Calibration Flag: 1 for calibrated, 0 for uncali')
    else:
        prihdr.set('CAL_FLAG', 0, 'Calibration Flag: 1 for calibrated, 0 for uncali')
    #endif
# Write the file
    hdulist.writeto(file_out, clobber=True)

    return file_out
#END of tsys_writeudbfits

def valid_miriad_dataset(filelist0):
    '''Returns True or False for valid or invalid Miriad datasets,
    checks for existnce of the directory, and then for flags, header,
    vartable, and visdata. Also returns names of valid and invalid datasets'''

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
        else:
            bad_filelist.append(filelist[j])
        #end if
    #endfor
    return otp, ok_filelist, bad_filelist

#End of valid_miriad_dataset
