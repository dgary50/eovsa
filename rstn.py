'''
   Module for retrieving RSTN/Penticton quiet Sun flux densities and evaluating
   expected solar flux densities expected for 2.1-m antennas, as a function of 
   frequency. '''
#
# History:
#  2014-Dec-13  DG
#    Written earlier, but updated for changes to NOAA ftp location and new
#    archive file with dates older than 45 days.  Also added lots of error-handling.
#  2015-May-29  DG
#    Converted from using datime() to using Time() based on astropy
#

import urllib2
import numpy as np
from util import Time
import sun_pos

def rd_rstnflux(t=None,f=None,recur=False):
    ''' Reads the RSTN/Penticton quiet Sun solar flux density for the date specified
        in the Time() object t.  Reads from file handle f, if supplied, or else
        attempts to retrieve the data from NOAA if the date is within 45 days of today.
        Otherwise, defaults to archive file /common/tmp/txt/radioflux.noa.  Under certain
        conditions, the routine calls itself after setting recur to True.
        
        On success, returns a 9-element list of frequencies and 
                            a 9-element list of corresponding flux densities
        
        On failure, returns None, None
        
        Although this routine could work with only the archive file, accessing the
        NOAA database is retained since it allows the routine to be used anywhere,
        not just on machines with access to the archfile location.
    '''    
    # Update these if locations change
    archfile = '/common/tmp/txt/radioflux.noa'
    #noaa_url = 'http://legacy_www.swpc.noaa.gov/ftpdir/lists/radio/45day_rad.txt'
    noaa_url = 'ftp://ftp.swpc.noaa.gov/pub/lists/radio/45day_rad.txt'
    
    if t is None:
        t = Time.now()
    datstr = t.datetime.strftime("%Y %b %d")
    #if datstr[9] == '0': datstr = datstr[:9]+datstr[10:]
    if f is None:
        today = Time.now()
        if today.mjd - t.mjd > 44:
            recur = True  # Set this to prevent unneeded recursive call
        else:
            try:
                f = urllib2.urlopen(noaa_url)
                print 'Data will be retrieved from NOAA 45-day file.'
            except:
                print 'NOAA 45-day file not reachable at',noaa_url
                recur = True  # Set this to prevent unneeded recursive call
        if recur:
            try:
                f = open(archfile,'r')
                print 'Data will be retrieved from archive file',archfile,'.'
            except:
                print 'Error: Archive file',archfile,'not found.'
                return None, None

    lines = f.readlines()
    f.close()
    frq = np.zeros(9,'int')
    flux = np.zeros(9,'float')
    for i,line in enumerate(lines):
        if line.find(datstr) != -1:
            for j in range(9):
                dat = lines[i+j+1].split()
                frq[j] = int(dat.pop(0))
                flxarr = np.array(dat,'int')
                good = np.where(flxarr != -1)[0]
                if len(good) != 0:
                    flux[j] = np.median(flxarr[good])
                else:
                    flux[j] = np.nan
            break
    if frq[0] == 0:
        if recur:
            print 'Error: Date',datstr,'not found.'
            return None, None
        else:
            print 'Warning: Date',datstr,'not found in 45-day file at',noaa_url
            print 'Data will be retrieved from archive file.'
            try:
                f = open(archfile,'r')
                frq, flux = rd_rstnflux(t,f,recur=True)
            except:
                print 'Error: File',archfile,'not found.'
                return None, None
    return frq, flux
    
def rstn2ant(frq,flux,fmhz,t=None,twometer=True):
    ''' Takes 9-element list of frequencies and corresponding flux densities from
        a call to rd_rstnflux(), and fits a 2nd-degree polynomial to the last 6
        frequency-flux density pairs.  
        
        Returns the values of the polynomial fit, adjusted for the 2.1-m antenna 
        nominal beam size (if twometer ==True), evaluated at the frequencies given 
        in the supplied frequency list fmhz.  
        
        Optional input parameter t is a Time() object with the date of the RSTN 
        fluxes, which is used to determine the solar disk radius, used in the 
        beam-size adjustment.  If omitted, today's date is used.
    '''
    if t is None:
        t = Time.now()
    # Perform 2nd-degree polynomial fit to input flux density values
    # Note that p is a polynomial object that returns the flux values
    # at the frequencies in its argument 
    idx = np.isfinite(frq[3:]) & np.isfinite(flux[3:])
    p = np.poly1d(np.polyfit(frq[3:][idx],flux[3:][idx],2))
    # Get size of this day's solar disk
    pa,b0,r = sun_pos.get_pb0r(t.mjd,arcsec=True)
    rp = (r+10.)/960.   # Radius of radio Sun in units of nominal 960" (10" added for radio limb)
    
    ''' Convert total flux to flux measured by the antenna, as follows:
    ;
    ;   Assuming we are pointed exactly at disk center, and that the Sun
    ;   is a nice, flat disk (valid only at high frequencies, but that
    ;   is where it makes the most difference, and then only for 2-m ants)
    ;   the actual flux measured is a gaussian (primary beam) truncated
    ;   by the solar disk.  The measured flux S is related to the total
    ;   flux, S_0, by
    ;                        / R_0
    ;      S = 2 S_0 / R_0^2 \  r exp[-(r/alpha)^2] dr
    ;                        / 0
    ;
    ;   where R_0 is the angular solar radius, alpha is the 1/e primary
    ;   beam width, and the integral over azimuthal angle:
    ;       / 2 pi
    ;       \  dphi = 2 pi
    ;       / 0
    ;   has already been performed.
    ;
    ;   The integral can be evaluated to yield:
    ;      S = S_0 (alpha/R_0)^2 { 1 - exp[-(R_0/alpha)^2]}
    ;
    ;   The half-power-beam-width HPBW = 2 sqrt(ln 2) alpha = 46.5'/f_GHz
    ;   for the 27m dishes, or 13.5 times this for 2m dishes.  Thus,
    ;
    ;         alpha/R_0 = 27.9'/f_GHz/ 16.0' r' = 1.745/(f_GHz r')
    ;
    ;   for the 27 m antennas, or 12.857 times this for 2.1m dishes, where
    ;   r' is the angular solar radius in units of the nominal 16.0'(960").
    '''
    const = 1.745
    if twometer: const = const*12.857
    arg = (const/(fmhz*rp/1000.))**2
    s = p(fmhz)*arg*(1.0-np.exp(-1./arg))
    return s
    
from datetime import datetime
def noaa2db(infile,outfile):
    ''' Read a NOAA text file containing RSTN and Penticton flux density
        measurements and pack the information into the given binary 
        database file.  If the date in the text file is found in the
        binary file, the information is replaced/updated--otherwise it
        is inserted or appended in date order.
        
        This is incomplete and untested--abandoned for now, but may be
        revisited.
    '''
    def line2mjd(line):
        ''' Return the modified Julian day number for a line that starts with a
            date of the form, e.g. 1999 Jan  1.  If the line does not start with
            that form, return -1
        '''
        try:
            mjd = Time(datetime.strptime(line.strip(),'%Y %b %d')).mjd
        except:
            mjd = -1
        return mjd
        
    try:
        f = open(infile,'r')
        lines = np.array(f.readlines())
        f.close()
    except:
        print 'Cannot find or open file',infile
        return
    fline = -1
    # Find the first line with a date
    for i,line in enumerate(lines):
        try:
            mjd = line2mjd(line)
            fline = i
            break
        except:
            pass
    if fline == -1:
        print 'No lines with dates found!'
        return
    idx = np.arange(len(lines)/11)*11 + fline
    mjds = []
    frq = np.zeros(9,'int')
    flux = np.zeros(9,'float')
    # Loop over lines with times
    for i,line in enumerate(lines[idx]):
        # Convert time to mjd
        mjds.append(line2mjd(line))
        # Read next 9 lines into frequency and flux arrays
        for j in range(9):
            dat = lines[idx[i]+j+1].split()
            frq[j] = int(dat.pop(0))
            flxarr = np.array(dat,'int')
            good = np.where(flxarr != -1)[0]
            if len(good) != 0:
                # Use median of good values
                flux[j] = np.median(flxarr[good])
            else:
                flux[j] = np.nan
            break
        
    mjds = np.array(mjds)
    
