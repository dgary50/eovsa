#
# READBSC
#
# Routines to read and manipulate the Bright Star Catalog.  These sources
# are used for optical pointing measurements, not radio observations.
# 
# 2013-Jun-23  DEG
#   Converted from original to use PyEphem routines
# 2015-May-29  DG
#   Converted from using datime() to using Time() based on astropy
# 2015-Oct-14  DG
#   Fixed a couple of bugs, including need to login with user/password
# 2015-Oct-18  DG
#   Added "mount" flag to startracktable() to indicate if this is an
#   azel mount (default, mount='azel'), or equatorial mount (mount=anything else)
# 2015-Oct-20  DG
#   Strange.  For some reason, eovsa_ha() did not work because the sources were
#   not "computed" for the EOVSA array.  Somehow, my earlier changes to some
#   routines must have broken this.  Now I just do the compute in starobs2dxeldel().
#   I realize that I also needed an hadec switch in that routine, to write out
#   HA and Dec coordinates for equatorial mounts (see eq_mountcal() in mountcal.py).
# 2016-Jan-09  DG
#   Added analyze_stars() routine, to automatically upload image files to
#   Astrometry.net, get the solutions, and create the starsolutions file.
#   Really cool!
#
# Must be run from Dropbox/PythonCode/Current directory

from numpy import array, zeros, ones, arange, where, argsort, sort, pi
from util import *
import ephem
import datetime as dt
#from Coordinates import *
#from Astrometry import *
from eovsa_lst import *
from eovsa_array import *
from ftplib import FTP

def readbsc(filename=None):
    ''' Read entire Bright Star Catalog and return a list of PyEphem sources
    '''
    if filename is None:
        filename = '../Star Pointing/BrightStarCatalog.txt'

    try:
        f = open(filename,'r')
    except:
        print 'readbsc: Could not open file',filename
        return None

    # Read two header lines (and ignore them)
    line = f.readline()
    line = f.readline()
    srcs = []

    # Loop over lines in file
    for line in f.readlines():
        srcs.append(ephem.readdb('HR'+line[1:5]+',f,'+line[6:17]+','+line[18:30]+','+line[53:57]+',2000'))
    f.close()
    # Do an initial compute so that name, ra, dec, etc. are accessible.  Will override with compute for
    # observer later.
    for src in srcs:
        src.compute()

    return srcs
    
def selectbsc(t, srcs, magrange):
    ''' Given a list of star sources, select from the list based on hour angle, 
        elevation, and magnitude range provided in magrange as [maglow, maghi]
    '''
    # Pare down the list by selecting stars within +3 and -9 hours of HA.
    # Pare down further by requiring stars within -30 and +45 degrees in Dec
    # Finally, find stars with magnitudes between 3.5 and 5
    halow = Dec_Angle('-09:00:00','hms').radians
    hahi  = Dec_Angle('03:00:00','hms').radians
    declow = Dec_Angle('-30:00:00','dms').radians
    dechi  = Dec_Angle('+45:00:00','dms').radians
    maglow = magrange[0]
    maghi  = magrange[1]
    # Make a list of good sources (those in the ranges specified)
    sel = []
    # Loop over sources
    for i,src in enumerate(srcs):
        # Make sure ha is in range -pi to pi
        ha = eovsa_ha(src,t)
        dec = src.dec
        mag = float(src.mag)
        good = (ha > halow) & (ha < hahi) & (dec > declow) & (dec < dechi) \
            & (mag > maglow) & (mag < maghi)
        if good:
            sel.append(i)

    # These selected stars should be reachable by the antennas, so now check them
    # for angular separation.
    nstars = len(sel)   # Number of stars so far selected
    widesep = zeros((nstars,nstars),bool)   # Boolean nstars x nstars array 
    for i in range(nstars-1):
        for j in range(i+1,nstars):
           # True where pairs of stars are > 20 deg apart
           widesep[i,j] = ephem.separation(srcs[sel[i]],srcs[sel[j]]) > (21.5*pi/180)
           widesep[j,i] = widesep[i,j]

    # Go through asep array row by row and mark stars for deletion with sep < 20 degrees
    idx = ones((nstars),bool)
    for i in range(nstars):
        if idx[i]:
            # Only look at rows for "True" columns
            x = widesep[:,i]
            x[:i+1] = True  # Do not delete stars in lower part of array
            idx = idx & x

    # This should be a list of remaining good stars
    ids = array(sel)[idx]
    print len(ids),'stars selected for date/time starting at:',t.iso
    print 'Number      RA         Dec      Magnitude'   
    fmt = '{0:<4} {1:>12} {2:>12} {3:>6}'
    for i in ids:
        dec = srcs[i].a_dec
        decstr = str(srcs[i].a_dec).zfill(10)
        if dec < 0:
            decstr = '-'+str(srcs[i].a_dec)[1:].zfill(10)
        print fmt.format(srcs[i].name[2:],str(srcs[i].a_ra),decstr,srcs[i].mag)

    return ids
    
def getbscnames(num=None, filename=None):
    ''' Looks up the common star name for each HR number in the list given in array num
    '''
    if filename is None:
        filename = '../Star Pointing/bsc5.dat'
        try:
            f = open(filename,'r')
        except:
            print 'getbscnames: Could not open file', filename
            return ''
        
    if num is None:
        print 'getbscnames: Must specify an ordered list of star numbers'
        return ''

    line = '9999'
    name = []
    for i in range(len(num)):
        while int(line[0:4]) != num[i]:
           line = f.readline()
        # This should be the line that matches the star number i
        n = line[4:7]
        greek = line[7:10]
        cnstl = line[10:14]
        alt = line[14:24]
        if greek == '   ':
            # No greek letter
            if n == '   ':
                # No number, so use alt
                name.append(alt)      # e.g. 'BD-10 6177' or 'CD-4015285'
            else:
                # Number, but no greek letter
                name.append(n+cnstl)  # e.g. ' 80 Peg' or '108 Aqr'
        else:
            # Greek letter
            name.append(greek+cnstl)    # e.g. 'Phi Peg' or 'Gam1And'

    f.close()
    
    return array(name)

def startracktable(t, names, srcs, ids, npts=25, mount='azel'):
    ''' Generate an RA_Dec track table for observing this list of stars
    '''
    nstars = len(ids)
    min2radians = RA_Angle('00:01:00','hms').radians
    min2mjd = 60.0001/86400.   # Slightly greater to avoid annoying times like 59.999.
    # Convert radians to antenna controller "user units" of 1/10000th of a degree
    r2u = 1800000./pi

    ovsa = ephem.Observer()
    ovsa.date = t.mjd - 15019.5  # Converts MJD to ephem date
    ovsa.lon = '-118.286953'
    ovsa.lat = '37.233170'
    ovsa.elevation = 1200
    ovsa.pressure = 0  # Eliminates refraction due to atmosphere
    
 
    filename = '../Star Pointing/startracktable.radec'
    try:
        f = open(filename,'w')
    except:
        print 'startracktable: Could not open output file',filename
        return

    # Generate star table file name
    datstr = t.iso[:10]
    outfile = '../Star Pointing/startable-'+datstr+'.txt'
    o = open(outfile,'w')
    
    ha_start = zeros(nstars)
    # Loop over sources
    for i in range(nstars):
        src = srcs[ids[i]]
        # First calculate HA
        ha_start[i] = eovsa_ha(src,t)

    # We now have nstars coordinates of date, with RA converted to an
    # initial HA.  Now we have to increment over time in steps of ha_minutes
    # for a specified number of pointings and identify the stars we want to observe based
    # on their Az, El visibility.

    # Find first star of interest (the one with the greatest hour angle)
    idx = argsort(ha_start)

    # Determine sky limits accounting for distance star can travel during
    # observation
    azlo = RA_Angle('45:00:00','dms').radians
    azhi = RA_Angle('325:00:00','dms').radians
    ello = Dec_Angle('15:00:00','dms').radians
    elhi = Dec_Angle('80:00:00','dms').radians
    # HA/Dec limits for equatorial mount antennas
    halo = Angle('-55:00:00','dms').radians
    hahi = Angle('55:00:00','dms').radians
    declo = Dec_Angle('-23:00:00','dms').radians
    dechi = Dec_Angle('44:00:00','dms').radians

    written = False   # Flag to say a star was written
    j = idx[nstars-1]  # Pointer to current star of interest, for cycling among the stars
    nchecked = 0
    ha_minutes = 0
    # Print to screen
    print 'Num   Name        RA(J2000)   Dec(J2000)    RA(Date)     Dec(Date)    Az(deg)      El(deg)    Mag    Time'
    print '==== ========== ============ ============ ============ ============ ============ =========== ===== ==========='
    # And write to outfile
    o.write('Num   Name        RA(J2000)   Dec(J2000)    RA(Date)     Dec(Date)    Az(deg)      El(deg)    Mag    Time   \n')
    o.write('==== ========== ============ ============ ============ ============ ============ =========== ===== ===========\n')
    dha_minutes = 4   # Length of time to stay on each star
    firstline = True
    # Entire duration is npts pointings, or npts*dha_minutes
    while ha_minutes < dha_minutes*npts:
        # print 'HA_Minutes is',ha_minutes
        for i in range(j,nstars):
            nchecked += 1
            newt = Time(t.mjd+ha_minutes*min2mjd,format='mjd')
            ovsa.date = newt.mjd - 15019.5  # Converts MJD to ephem date
            src = srcs[ids[i]]
            src.compute(ovsa)
            ha = eovsa_ha(src,newt)  # updates source to new time
            dec = src.dec
            # Get Az, El
            az, el = src.az, src.alt 
            azgood = (az > azlo) & (az < azhi)
            elgood = (el > ello) & (el < elhi)
            # Case of equatorial mount--just use same "good" variables but now
            # corresponding to HA and Dec
            if mount != 'azel':
                azgood = (ha > halo) & (ha < hahi)
                elgood = (dec > declo) & (dec < dechi)
            if azgood & elgood:
                # This star is okay, so create entries for it and
                # jump out of loop
                rad = int(src.ra*r2u)     
                decd = int(src.dec*r2u)
                mjd = t.mjd + ha_minutes*min2mjd # Time for this entry
                mjdint = int(mjd)   # Integer part of day
                ms = int(round((mjd-mjdint)*86400*1000))   # Time of day as ms
                if firstline:
                    line = str(rad)+' '+str(decd)+' '+str(mjdint)+' '+str(0)+'\n'  # Start first line at 0 UT (because why not?)
                else:
                    line = str(rad)+' '+str(decd)+' '+str(mjdint)+' '+str(ms)+'\n'
                f.write(line)
                mjd = t.mjd + (ha_minutes + 0.010/60.)*min2mjd # Add 10 msec
                mjdint = int(mjd)   # Integer part of day
                ms = int(round((mjd-mjdint)*86400*1000))   # Time of day as ms
                if firstline:
                    line = str(rad)+' '+str(decd)+' '+str(mjdint)+' '+str(10)+'\n'
                    firstline = False
                else:
                    line = str(rad)+' '+str(decd)+' '+str(mjdint)+' '+str(ms)+'\n'
                f.write(line)
                mjd = t.mjd + (ha_minutes + dha_minutes - 0.020/60.)*min2mjd # Add dha_minutes - 20 ms
                mjdint = int(mjd)   # Integer part of day
                ms = int(round((mjd-mjdint)*86400*1000))   # Time of day as ms
                line = str(rad)+' '+str(decd)+' '+str(mjdint)+' '+str(ms)+'\n'
                f.write(line)
                mjd = t.mjd + (ha_minutes + dha_minutes - 0.010/60.)*min2mjd # Add dha_minutes - 10 ms
                mjdint = int(mjd)   # Integer part of day
                ms = int(round((mjd-mjdint)*86400*1000))   # Time of day as ms
                line = str(rad)+' '+str(decd)+' '+str(mjdint)+' '+str(ms)+'\n'
                f.write(line)
                dec = src.a_dec
                adecstr = str(src.a_dec).zfill(10)
                if dec < 0:
                    adecstr = '-'+str(src.a_dec)[1:].zfill(10)
                dec = src.dec
                decstr = str(src.dec).zfill(10)
                if dec < 0:
                    decstr = '-'+str(src.dec)[1:].zfill(10)
                # Print to screen
                fmt = '{0:<4} {1:<10} {2:>12} {3:>12} {4:>12} {5:>12} {6:>12} {7:>11} {8:5.2f} {9:>12}'
                print fmt.format(src.name[2:], names[i], str(src.a_ra), adecstr, str(src.ra), decstr, \
                                 str(src.az), str(src.alt), src.mag, newt.iso[:19])
                # And write to outfile
                fmt += '\n'
                o.write(fmt.format(src.name[2:], names[i], str(src.a_ra), adecstr, str(src.ra), decstr, \
                                 str(src.az), str(src.alt), src.mag, newt.iso[:19]))
                written = True
                # print 'J written =',j
                j = i+1  # advance to next star and next time
                ha_minutes += dha_minutes
                break

        if j is nstars:
            #print 'J is equal to NSTARS, so setting to zero. Nchecked is ',nchecked
            j = 0
        if nchecked >= nstars:
            #print 'Checked all of the stars for', ha_minutes,'!'
            ha_minutes += dha_minutes
            nchecked = 0
        if (i is nstars-1) and (not written):
            #print 'Got through all the stars so far, so go back to beginning'
            j = 0
        else:
            #print 'Got our star. Next J =',j
            nchecked = 0
            written = False   # Reset the "written" flag
                
    # All done, so repeat the last line twice replacing time for end of UT day so antenna will hold on source
    f.write(line[:line.rfind(' ')]+' 86399989\n')
    f.write(line[:line.rfind(' ')]+' 86399999')
    f.close()
    o.close()

def starobs2dxeldel(filename=None,hadec=False):
    """Takes the star observations (calculated RA,Dec and measured RA,Dec) and
        calculates the differences in (Az, El).  The expected file format has the
        first three lines, below, and then a series of lines like the third:

        Num   Name        RA(J2000)    Dec(J2000)    RA(Meas)      Dec(Meas)  Time(Meas)
        ==== ========== ============ ============= ============ ============= ==========
          37 BD-18   14  0:12:10.000 -17:56:18.000  0:15:54.000 -18:35:23.000  02:38:00
        """

    if filename is None:
        filename = '/home/dgary/Dropbox/Python/Star Pointing/2012Jan09_starobs.txt'

    f = open(filename, 'r')
    lines = f.readlines()
    f.close()
    idx = filename.find('solutions')
    outfile = filename[0:idx]+'reduction'+filename[idx+9:]
    o = open(outfile,'w')
    line = lines[0]
    didx = filename.find('.txt')
    datstr = filename[didx-10:didx]
    # Loop over the lines in the file
    for line in lines[2:]:
        if line == '': break
        num = int(line[0:4])
        name = line[5:15]
        RA_J2000 = RA_Angle(line[16:28].strip(),'hms')
        Dec_J2000 = Dec_Angle(line[29:42].strip(),'dms')
        RA_Obs = RA_Angle(line[43:55].strip(),'hms')
        Dec_Obs = Dec_Angle(line[56:69].strip(),'dms')
        timstr = line[71:79]
        t = Time(datstr+' '+timstr)
        # Make a couple of fake sources, one for the nominal position, and one for the observed position
        srcnom = ephem.readdb('Nominal'+',f,'+line[16:28].strip()+','+line[29:42].strip()+','+'1.0,2000')
        srcobs = ephem.readdb('Observed'+',f,'+line[43:55].strip()+','+line[56:69].strip()+','+'1.0,2000')
        # Compute sources for EOVSA array at current time
        aa = eovsa_array()
        aa.set_jultime(t.jd)
        srcnom.compute(aa)
        srcobs.compute(aa)
        HA_Obs = eovsa_ha(srcnom,t)  # Updates srcnom for OVRO location and time Time_Obs)
        HA_Obs = eovsa_ha(srcobs,t)  # Updates srcobs for OVRO location and time Time_Obs)
        if hadec:
            ha1 = HA_Obs
            dec1 = Dec_Obs.get()
        az0 = srcnom.az
        el0 = srcnom.alt
        az1 = srcobs.az
        el1 = srcobs.alt

        # Calculate differences (in radians)
        d_xel = (az1 - az0) * cos(el1)
        d_el = (el1 - el0)

        # Write out a new file in the style of KSRBL (for now)
        timesec = int((t.mjd-int(t.mjd))*86400)
        dra = (srcobs.a_ra - srcnom.a_ra)*180/pi # Degrees
        ddec = (srcobs.a_dec - srcnom.a_dec)*180/pi # Degrees
        fmt = '{0:<8} {1:5d} {2:7.3f} {3:7.3f} {4:7.3f} {5:7.3f} {6:6.3f} {7:6.3f} {8:8.3f} {9:7.3f} {10:6.3f} {11:6.3f}'
        if hadec:
            # If the hadec switch is set, write out HA and Dec coordinates (ha1 and dec1)
            print fmt.format(name[0:8],timesec,srcnom.a_ra*180/pi, srcnom.a_dec*180/pi,
                        srcobs.a_ra*180/pi, srcobs.a_dec*180/pi, dra, ddec,
                        ha1*180/pi, dec1*180/pi, d_xel*180/pi, d_el*180/pi)
            o.write(fmt.format(name[0:8],timesec,srcnom.a_ra*180/pi, srcnom.a_dec*180/pi,
                        srcobs.a_ra*180/pi, srcobs.a_dec*180/pi, dra, ddec,
                        ha1*180/pi, dec1*180/pi, d_xel*180/pi, d_el*180/pi)+'\n')
        else:
            # Otherwise, write out AZ and EL coordinates (az1 and el1)
            print fmt.format(name[0:8],timesec,srcnom.a_ra*180/pi, srcnom.a_dec*180/pi,
                        srcobs.a_ra*180/pi, srcobs.a_dec*180/pi, dra, ddec,
                        az1*180/pi, el1*180/pi, d_xel*180/pi, d_el*180/pi)
            o.write(fmt.format(name[0:8],timesec,srcnom.a_ra*180/pi, srcnom.a_dec*180/pi,
                        srcobs.a_ra*180/pi, srcobs.a_dec*180/pi, dra, ddec,
                        az1*180/pi, el1*180/pi, d_xel*180/pi, d_el*180/pi)+'\n')
    o.close()

def do_stars(yr, mo, da, hr, mn, npts=25, mount='azel'):
    userpass = 'admin:observer@'
    srcs = readbsc()
    t = Time(dt.datetime(yr,mo,da,hr,mn),format='datetime')
    ids = selectbsc(t, srcs, [5.2,5.3])
    num = zeros(len(ids))
    for i,idx in enumerate(ids):
        num[i] = int(srcs[idx].name[2:])
    names = getbscnames(num)

    startracktable(t, names, srcs, ids, npts, mount)

    # Connect to ACC /parm directory and transfer tracktable file
    try:
        acc = FTP('acc.solar.pvt')
        acc.login('admin','observer')
        acc.cwd('parm')
        # Send tracktable file to ACC
        filename = '../Star Pointing/startracktable.radec'
        f = open(filename,'r')
        acc.storlines('STOR startracktable.radec',f)
        f.close()
    except:
        print 'Could not transfer startracktable.radec file.  ACC is down?'

def analyze_stars(yr, mo, da, radius=3):
    import os, glob, time
    from astropy.io import fits
    import subprocess
    t = Time(dt.datetime(yr,mo,da,0,0),format='datetime')
    fileloc = '/home/sched/Dropbox/PythonCode/Star Pointing/'
    apikey = 'cqjhjsprirwttecb'
    radeg = 180./pi
    datestr = t.iso[:10]
    os.chdir(fileloc+datestr)
    # Read startable
    f = open(fileloc+'startable-'+datestr+'.txt','r')
    table = f.readlines()
    f.close()
    # Remove and save header lines of table
    header = array([table[0].strip(),table[1].strip()])
    table = table[2:]
    # Read times, RA and Dec for each line of table
    times = []
    RA_J2000 = []
    Dec_J2000 = []
    for line in table:
        times.append(line.strip()[-19:])
        RA_J2000.append(RA_Angle(line[16:28].strip(),'hms'))
        Dec_J2000.append(Dec_Angle(line[29:42].strip(),'dms'))
    # Convert times to Time object
    times = Time(array(times))
    # Read image file list
    filelist = sort(glob.glob('*.fts'))
    idxlist = []   # List of indexes into table that were processed
    ftimes = []    # List of times for files that were processed
    for file in filelist:
        # Gather information from file header
        hdulist = fits.open(file)
        f_datestr = hdulist[0].header['date-obs']
        f_timestr = hdulist[0].header['time-obs']
        f_time = Time(f_datestr+' '+f_timestr)
        # Identify star from time
        try:
            idx = where(f_time > times)[0][-1]
            skip = False
        except:
            print f_time.iso,'is before first time in file...skipping'
            skip = True
        if not skip:
            # Submit this image to Astrometry.net, and wait for processing
            # to complete.
            command = ['python','/common/python/current/astronet.py','--apikey='+apikey,'--upload='+file,
                       '--ra='+str(RA_J2000[idx].radians*radeg),'--dec='+str(Dec_J2000[idx].radians*radeg),
                       '--radius='+str(radius),'--wcs=wcs_'+file[:-3]+'fits']
            # print 'Sending:',command
            p = subprocess.Popen(command,stdout=subprocess.PIPE)
            tstart = time.time()
            lines = p.stdout.readlines()
            print 'Result is',lines[-1].strip()
            print 'Took',time.time() - tstart,'seconds.'
            if lines[-1][:12] == 'Wrote to wcs': 
                idxlist.append(idx)
                ftimes.append(f_timestr)
    # Write output file
    f = open(fileloc+'starsolutions-'+datestr+'.txt','w')
    f.write('Num   Name        RA(J2000)    Dec(J2000)    RA(Meas)      Dec(Meas)  Time(Meas)\n')
    f.write('==== ========== ============ ============= ============ ============= ==========\n')
    wcslist = sort(glob.glob('*.fits'))
    for i,file in enumerate(wcslist):
        hdulist = fits.open(file)
        ra = RA_Angle(hdulist[0].header['crval1']/radeg)
        dec = Dec_Angle(hdulist[0].header['crval2']/radeg)
        outline = table[idxlist[i]][:41]+'  '+ra.get('hms')[:-1]+'  '+dec.get('dms')[:-2]+'    '+ftimes[i]
        f.write(outline+'\n')
    f.write('\n')
    f.close()
         
        
            

    
