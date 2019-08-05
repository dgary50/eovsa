#
# FITS writing routines, originally intended for XSP files, but now
# expanded to other types of dynamic spectrum files.
#
# History
#  2019-08-05  DG
#    Updates to tp_writefits() to extend to other types of files.
#    Also changed deprecated "clobber" keyword to "overwrite"


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

def daily_xsp_writefits(xdat, pdata):
    '''This takes the dictionary output from a read_idb, and pdata from
    daily_xsp, and creates a FITS file. The data is the primary output.
    '''

    file_out = ''
    if xdat == None or len(xdat) == 0:
        print 'xsp_writefits: No data input'
        return file_out
#   end

#UDB allday fits files
    xspfitsdir = '/data1/eovsa/fits/XSP/'
    if os.path.isdir(xspfitsdir) == False:
        print "daily_xsp_writefits: creating "+xspfitsdir
        os.mkdir(xspfitsdir)
    #end if
#create a filename, just use the start time
    t = xdat['time']
    print t[0]
    t0 = Time(t[0], format='jd') #The format is to tell the Time object about the input time
    t01 = t0.isot
    yr = t01[0:4]
    mn = t01[5:7]
    dy = t01[8:10]
    file_out = 'XSP'+yr+mn+dy+'.fts'
#add directory
    outdir = xspfitsdir+'/'+yr+'/'
    if os.path.isdir(outdir) == False:
        print "daily_xsp_writefits: creating "+outdir
        os.mkdir(outdir)
    #end if
    file_out = outdir+'/'+file_out

    date_obs = t01
#Convert to unix time, and add seconds to get date_end
#http://astropy.readthedocs.org/en/latest/time
    print "date_obs: ", date_obs
    t0 = Time(date_obs)
    dt = int(86400*(max(t)-min(t)))
    t1 = Time(t0.unix+dt, format = 'unix')
    date_end = t1.isot
    print "date_end: ", date_end

# Create the primary header
    tpwr = pdata
    hdu = fits.PrimaryHDU(tpwr)

# Set up the extensions: sfreq, ut
    sfreq = xdat['fghz']
    col1 = fits.Column(name='sfreq', format='E', array = sfreq)
    cols1 = fits.ColDefs([col1])
    tbhdu1 = fits.BinTableHDU.from_columns(cols1)
    tbhdu1.name = 'SFREQ'

# Split up mjd into days and msec, really intuitive syntax, thanks python
    ut = Time(t, format='jd') #Format is for the input time
    ut = ut.mjd
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
    hdulist = fits.HDUList([hdu, tbhdu1, tbhdu3])

#primary header
    prihdr = hdulist[0].header
#Header information, strip last character from strings
    obj_id = strip_non_printable(xdat['source'])
    temp_out = file_out.split("/")
    temp_out = temp_out[len(temp_out)-1]
    prihdr.set('FILENAME', temp_out)
    prihdr.set('ORIGIN', 'NJIT', 'Institute where file was written')
    prihdr.set('TELESCOP', 'EOVSA', 'Expanded Owens Valley Solar Array')
    prihdr.set('OBJ_ID', obj_id, 'Object ID')
    prihdr.set('TYPE', 1, 'Spectrum')
    prihdr.set('DATE_OBS', date_obs, 'Start date/time of observation')
    prihdr.set('DATE_END', date_end, 'End date/time of observation')
    prihdr.set('FREQMIN', min(sfreq), 'Min freq in observation (GHz)')
    prihdr.set('FREQMAX', max(sfreq), 'Max freq in observation (GHz)')
    prihdr.set('XCEN', 0.0, 'Antenna pointing in arcsec from Sun centre')
    prihdr.set('YCEN', 0.0, 'Antenna pointing in arcsec from Sun centre')
    prihdr.set('POLARIZA', 'XX, YY', 'Polarizations present')
    prihdr.set('RESOLUTI', 0.0, 'Resolution value')
# Write the file
    hdulist.writeto(file_out, overwrite=True)

    return file_out

def tp_writefits(out, med, filestem='', outpath='/data1/eovsa/fits/flares/'):
    '''This takes the dictionary output from read_idb, corrected with
       autocorrect_tp.py, and the background-subtracted median data 
       from it, and creates a FITS file.  Output is the filename.
    '''

    import os
    file_out = ''
    if out == None or len(out) == 0:
        print 'tp_writefits: No data input'
        return file_out
#   end

#create a filename, just use the start time
    t = out['time']
    print t[0]
    t0 = Time(t[0], format='jd') #The format is to tell the Time object about the input time
    t01 = t0.iso
    yr = t01[0:4]
    mm = t01[5:7]
    dy = t01[8:10]
    hr = t01[11:13]
    mn = t01[14:16]
    file_out = 'EOVSA_'+filestem+yr+mm+dy+hr+mn+'.fts'
#flare fits files
    if os.path.isdir(outpath) == False:
        print "tp_writefits: creating "+outpath
        os.mkdir(outpath)
#add yr directory
    outdir = outpath+'/'+yr+'/'
    if os.path.isdir(outdir) == False:
        print "daily_xsp_writefits: creating "+outdir
        os.mkdir(outdir)
#add mn directory
    outdir = outpath+'/'+yr+'/'+mm+'/'
    if os.path.isdir(outdir) == False:
        print "daily_xsp_writefits: creating "+outdir
        os.mkdir(outdir)
#add dy directory
    outdir = outpath+'/'+yr+'/'+mm+'/'+dy+'/'
    if os.path.isdir(outdir) == False:
        print "daily_xsp_writefits: creating "+outdir
        os.mkdir(outdir)
    file_out = outdir+'/'+file_out

    date_obs = t01
#Convert to unix time, and add seconds to get date_end
#http://astropy.readthedocs.org/en/latest/time
    print "date_obs: ", date_obs
    t0 = Time(date_obs)
    dt = int(86400*(max(t)-min(t)))
    t1 = Time(t0.unix+dt, format = 'unix')
    date_end = t1.iso
    print "date_end: ", date_end

# Create the primary header
    tpwr = med
    hdu = fits.PrimaryHDU(tpwr)

# Set up the extensions: sfreq, ut
    sfreq = out['fghz']
    col1 = fits.Column(name='sfreq', format='E', array = sfreq)
    cols1 = fits.ColDefs([col1])
    tbhdu1 = fits.BinTableHDU.from_columns(cols1)
    tbhdu1.name = 'SFREQ'

# Split up mjd into days and msec, really intuitive syntax, thanks python
    ut = Time(t, format='jd') #Format is for the input time
    ut = ut.mjd
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
    hdulist = fits.HDUList([hdu, tbhdu1, tbhdu3])

#primary header
    prihdr = hdulist[0].header
#Header information, strip last character from strings
    obj_id = strip_non_printable(out['source'])
    temp_out = file_out.split("/")
    temp_out = temp_out[len(temp_out)-1]
    prihdr.set('FILENAME', temp_out)
    prihdr.set('ORIGIN', 'NJIT', 'Institute where file was written')
    prihdr.set('TELESCOP', 'EOVSA', 'Expanded Owens Valley Solar Array')
    prihdr.set('OBJ_ID', obj_id, 'Object ID')
    if filestem == 'TP_':
        prihdr.set('TYPE', 1, 'Total Power Dynamic Spectrum')
    elif filestem == 'X_':
        prihdr.set('TYPE', 2, 'Cross Power Dynamic Spectrum')
    else:
        prihdr.set('TYPE', 0, 'Spectrum Type Undefined')
    prihdr.set('DATE_OBS', date_obs, 'Start date/time of observation')
    prihdr.set('DATE_END', date_end, 'End date/time of observation')
    prihdr.set('FREQMIN', min(sfreq), 'Min freq in observation (GHz)')
    prihdr.set('FREQMAX', max(sfreq), 'Max freq in observation (GHz)')
    prihdr.set('XCEN', 0.0, 'Antenna pointing in arcsec from Sun centre')
    prihdr.set('YCEN', 0.0, 'Antenna pointing in arcsec from Sun centre')
    prihdr.set('POLARIZA', 'I', 'Polarizations present')
    prihdr.set('RESOLUTI', 0.0, 'Resolution value')
# Write the file
    hdulist.writeto(file_out, overwrite=True)

    return file_out
    