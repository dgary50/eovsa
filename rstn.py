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
#  2016-May-07  DG
#    All of a sudden, closing a urlopen() file handle is hanging, so skip it.
#  2020-Jun-04  OG
#    Added new funtions to extract the RSTN data from NOAA or from the archive
#    file and write them to the SQL database
#    These are: 
#       - rstnfluxfromnoaa(dt) extracts flux data from the NOAA website for the
#         specified date
#       - rstnfluxfromtextarchive(startdt, enddt) extracts data from the old
#         noaa text file.
#       - writerstn2sql(data,t=None) writes the rstn flux data stored in data
#         to SQL at the specified time or the current time if none given.
#       - writerstnprev2sql() writes the previous days flux data to SQL
#       - rstntext2sql(startdt, enddt, logfile = None) writes flux values
#         from the old RSTN flux text archive to the SQL database in the 
#         date range startdt upt to but not including enddt. logfile is
#         is a file that will contain the dates of missing data.
#       - sql2rstn(t): extracts the rstn flux values and SQL_timestamp from
#         the SQL database
#
#    rd_rstndata() has been updated to allow reading from SQL.
#  2020-Jun-13  DG
#    Fixed some problems with rstnfluxfromnoaa() no longer being defined, and
#    rd_rstnflux now allows data to be returned from "yesterday" if current
#    date has not been written yet.

import urllib2
import numpy as np
from util import Time
import sun_pos
from time import strptime
import cal_header as ch
from stateframe import extract

def rd_rstnflux(t=None,f=None,recur=False):
    ''' Reads the RSTN/Penticton quiet Sun solar flux density for the date specified
        in the Time() object t.  Reads from file handle f, if supplied, or else
        attempts to retrieve the data from the SQL database. If this fails it
        will attempt to retrieve the data from NOAA if the date is within 45 days of
        today. Otherwise, defaults to archive file /common/tmp/txt/radioflux.noa.  
        Under certain conditions, the routine calls itself after setting recur to True.
        
        On success, returns a 9-element list of frequencies and 
                            a 9-element list of corresponding flux densities
        
        On failure, returns None, None
        
        If it finds data from another source other than SQL it will update the SQL
    '''    

    today = Time.now()
    if t is None:
        t = today
        
    datstr = t.datetime.strftime("%Y %b %d")
    if f is None:
        #try to get data from SQL
        data, sqlt = sql2rstn(t)
        if sqlt is None:
            print "RD_RSTNFLUX: SQL database error!"
        elif today.mjd - sqlt.mjd < 2:
            print "RD_RSTNFLUX: Today's flux data are not available yet--using yesterday's"
        else:
            # Data for this date is not in data base
            print "RD_RSTNFLUX: Error: Could not find ", datstr, " in SQL."
            if today.mjd - t.mjd > 44:
                recur = True  # Set this to prevent unneeded recursive call
            # Try to get data from NOAA
            data = rstnfluxfromnoaa(t)
            if data is not None:
                sqlt = Time(np.floor(data[0].mjd)+0.125,format='mjd')
                writerstn2sql(data, sqlt)
            else:
                print "RD_RSTNFLUX: Error: Could not retrieve ", datstr, " from NOAA"
                recur = True  # Set this to prevent unneeded recursive call

        if recur:
            enddt = t + TimeDelta(86400.0, format='sec')
            d = rstnfluxfromtextarchive(t, enddt)
            if d is None:
                print "RD_RSTNFLUX: Error: Could not find ", datstr, "in text archive"
                return None, None
            else:
                data=d[0]
                sqlt = Time(np.floor(data[0].mjd)+0.125,format='mjd')
                writerstn2sql(data, sqlt)
    else:
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
            print 'RD_RSTNFLUX: Error: Date',datstr,'not found.'
            return None, None
        else:
            return frq, flux
        
    frq = data[1] * 1000
    arrsize = data[2].shape
    n = arrsize[0]
    flux = np.zeros(n,'float')
    for j in range(n):
        good = np.where(data[2][j] !=-1 )[0]
        if len(good) != 0:
            flux[j] = np.median(data[2][j][good])
        else:
            flux[j] = np.nan
    
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
    
freq = np.array([0.245, 0.41, 0.61, 1.415, 2.695, 2.8, 4.995, 8.8, 15.4], dtype = np.float32)

def rstnfluxfromnoaa(t):
    ''' If given date is today, attempts to read from current NOAA file, otherwise
        reads from 45day file. Fails (returns None) if date are more than 45 days ago.
    '''
    if Time.now().mjd - t.mjd < 1:
        return rstnfluxfromcurrentnoaa()
    else:
        return rstnfluxfrom45daynoaa(t)
        
def rstnfluxfrom45daynoaa(dt):
    """This extracts the requested day's RSTN noon flux data from:
        ftp://ftp.swpc.noaa.gov/pub/lists/radio/45day_rad.txt
    A list is returned with each element as follows:
    
    0 - timestamp: Astropy Time which is the date on which the data was
        collected. This should match dt
        
    1 - freq: A float32 numpy array containing the 9 frequencies in GHz
    
    2 - data: The flux data which is a 9x7 int16 numpy array.
    
    If the ftp failed or the specified date is not found None is returned."""
    
    noaa_url = 'ftp://ftp.swpc.noaa.gov/pub/lists/radio/45day_rad.txt'
    try:
        f = urllib2.urlopen(noaa_url)
        lines = f.readlines()
    except:
        print "Could not read from ",noaa_url
        return None
    
    # The next several blocks below can be shortened by the more pythonic code:
    for i,line in enumerate(lines):     # Provides counter "i"
        line = line.strip()
        if line[:3] == 'MHZ':
            break
        
    if i+1 == len(lines):
        print "No data found in ",noaa_url
        return None
        
    lines = lines[i+1:]
    lines = np.array(lines)             # Converts to numpy array, for where()
    clean_lines = lines[np.where(lines != '\n')]  # Eliminates all empty lines
    clean_lines = [l.replace('\n','') for l in clean_lines] # Removes \n in remaining lines
    
    for i in range(0, len(clean_lines), 10):
        datestr = "%04d-%02d-%02d" % (int(clean_lines[i][0:4]), strptime(clean_lines[i][5:8], '%b').tm_mon, int(clean_lines[i][9:]))
        timestamp = Time(datestr, out_subfmt = 'date')
        if np.floor(dt.mjd) == np.floor(timestamp.mjd):
            data = np.zeros((9, 7), dtype = np.int16)
            for j in range(1, 10):
                d = clean_lines[i + j].split()
                data[j - 1] = list(map(int, d[1:8]))
                    
            print "Data successfully read from ",noaa_url," for date ",dt.iso
            return [timestamp, freq, data]
    
    print "No data found for data ",dt.iso," in ",noaa_url
    return None
    
def rstnfluxfromcurrentnoaa():
    """This extracts the first RSTN data block from:
        ftp://ftp.swpc.noaa.gov/pub/lists/radio/rad.txt
    A list is returned with each element as follows:
    
    0 - timestamp: Astropy Time which is the date on which the data was
        collected. This should match dt
        
    1 - freq: A float32 numpy array containing the 9 frequencies in GHz
    
    2 - data: The flux data which is a 9x7 int16 numpy array.
    
    If the ftp failed or the specified date is not found None is returned."""
    
    data = []
    success = True
    noaa_url = 'ftp://ftp.swpc.noaa.gov/pub/lists/radio/rad.txt'
    try:
        f = urllib2.urlopen(noaa_url)
        lines = f.readlines()
    except:
        print "Could not read from ",noaa_url
        return None
        
    for i,line in enumerate(lines):     # Provides counter "i"
        line = line.strip()
        if line[:3] == 'MHZ':
            break

    if i+1 == len(lines):
        print "No data found in ",noaa_url
        return None

    lines = lines[i+1:]
    lines = np.array(lines)             # Converts to numpy array, for where()
    clean_lines = lines[np.where(lines != '\n')]  # Eliminates all empty lines
    clean_lines = [l.replace('\n','') for l in clean_lines] # Removes \n in remaining lines
    
    datestr = "%04d-%02d-%02d" % (int(clean_lines[0][0:4]), strptime(clean_lines[0][5:8], '%b').tm_mon, int(clean_lines[0][9:]))
    timestamp = Time(datestr, out_subfmt = 'date')
    
    t = Time.now()
    if np.floor(timestamp.mjd) != np.floor(t.mjd-1.0):
        print "No data found for previous day (", Time(t.mjd-1.0, 'mjd', out_subfmt = 'date'), ")"
        return None
        
    data = np.zeros((9, 7), dtype = np.int16)
    for j in range(1, 10):
        d = clean_lines[j].split()
        data[j - 1] = list(map(int, d[1:8]))
    
    return [timestamp, freq, data]
    
def rstnfluxfromtextarchive(startdt, enddt):
    """This function extracts RSTN noon flux data from the old archive
    text file (/common/tmp/txt/radioflux.noa). It extracts data from the
    dates spanning startdt up to but not including enddt. It returns a
    list of the data in the date range. Each element in the list is
    another list as follows:
    
    The returned list contains the extracted data with each element
    as follows:
     
    0 - timestamp: Astropy Time which is the date on which the data was
        collected. This should match dt
        
    1 - freq: A float32 numpy array containing the 9 frequencies in GHz
    
    2 - data: The flux data which is a 9x7 int16 numpy array.
    
    If no data in the specified date range is found then None is returned."""
    
    archfile = '/common/tmp/txt/radioflux.noa'
    data = []
    try:
        f = open(archfile,'r')
        #print 'Data will be retrieved from archive file',archfile,'.'
        lines = f.readlines()
        f.close()
    except:
        print "Could not open",archfile
        return None
    
    for i,line in enumerate(lines):     # Provides counter "i"
        line = line.strip()
        if line[:5] == '# MHz':
            break
        
    if i+1 == len(lines):
        print "No data found in ",archfile
        return None
        
    lines = lines[i+1:]
    lines = np.array(lines)             # Converts to numpy array, for where()
    clean_lines = lines[np.where(lines != '\r\n')]  # Eliminates all empty lines
    clean_lines = [l.replace('\r\n','') for l in clean_lines] # Removes \n in remaining lines
    
    for i in range(0, len(clean_lines), 10):
        print clean_lines[i][:17]
        if clean_lines[i][:17] == ':Solar_Radio_Flux':
            datestr = "%04d-%02d-%02d" % (int(clean_lines[i][19:23]), strptime(clean_lines[i][24:27], '%b').tm_mon, int(clean_lines[i][28:]))
            timestamp = Time(datestr, out_subfmt = 'date')
            if np.floor(timestamp.mjd) >= np.floor(startdt.mjd) and np.floor(timestamp.mjd) < np.floor(enddt.mjd):
                fluxarr = np.zeros((9,7), dtype = np.int16)
                for j in range(1, 10):
                    d = clean_lines[i+j].split()
                    fluxarr[j - 1] = list(map(int, d[1:8]))
                data.append([timestamp, freq, fluxarr])
    
    if len(data)==0:
        print "No data found in specified date range."
        return None
    else:
        return data

def writerstn2sql(data,t=None):
    """This routine tries to write RSTN flux data stored in data to SQL.
    If t is supplied, then this is used for the SQL_timestamp, otherwise
    the current time is used.
    data is a list with the following elements:
    
    0 - timestamp: Astropy Time which is the date on which the data was
        collected. This should match dt
        
    1 - freq: A float32 numpy array containing the 9 frequencies in GHz
    
    2 - data: The flux data which is a 9x7 int16 numpy array."""
    
    if t is None:
        t=Time.now()
    
    if ch.rstnflux2sql(data, t):
        #If data written to database, display success message
        print t.iso + ": Data successfully written to database."
    else:
        #If data failed to write, show fail message
        print t.iso + ": Failed to write data to database."
    
def writerstnprev2sql():
    """This routine extracts and writes the previous days RSTN flux values
    from NOAA and writes them to SQL."""
    
    nt = Time.now()
    t = Time(np.floor(nt.mjd) - 0.875, format = 'mjd')
    data = rstnfluxfromcurrentnoaa()
    if data is None:
        #No data extracted, display error message
        print t.iso + ": No data found."
    else:
        #See if data already in database
        pd, sqlt = sql2rstn(t)
        if sqlt is None:
            if ch.rstnflux2sql(data, t):
                #If data written to database, display success message
                print t.iso + ": Data successfully written to database."
            else:
                #If data failed to write, show fail message
                print t.iso + ": Failed to write data to database."
        else:
            print t.iso + ": Data already in database."
            
def rstntext2sql(startdt, enddt, logfile = None):
    """This routine extracts data from the old archive text file and
    writes it to SQL. It will output a list of dates that were not
    archived. The SQL time will be the date of the data at 0300.
    The program checks to see if there is current data already for
    the date range. It will NOT overwrite a record if it is already
    present.
    
    startdt and enddt are the dates that will be written to SQL from
    startdt up to but not including enddt."""
    
    data=rstnfluxfromtextarchive(startdt, enddt)
    if data is None:
        print "No RSTN data found in range ", startdt.iso, " to ", enddt.iso
        return None
    
    print "Processing data from: ", startdt.iso, " to ", enddt.iso
     
    offset = int(np.floor(startdt.mjd))
    days = int(np.floor(enddt.mjd)) - offset
    
    processed = np.zeros((days), dtype=bool)
    recordswritten=0
    existingrecords=0
    for d in data:
        print "Processing Date: ", d[0].iso
        i = int(np.floor(d[0].mjd)) - offset
        processed[i] = True
        sqltime = Time(np.floor(d[0].mjd) + 0.125, format = 'mjd')
        xml, buf = ch.read_cal(12, sqltime)
        if buf is not None:
            sqltime_read = Time(extract(buf, xml['SQL_timestamp']),format='lv')
            if np.floor(sqltime_read.mjd) != np.floor(sqltime.mjd): buf = None
        
        if buf is None:
            if ch.rstnflux2sql(d, sqltime): 
                recordswritten += 1
                print "Record Written"
            else:
                print "Record Write Failed!"
        else:
            print "Record Exists."
            existingrecords += 1
    
    if logfile is None: logfile = "/tmp/missingrstn.txt"
    
    f = open(logfile, "w")
    for i in range(days):
        if not processed[i]: f.write(Time(float(i+offset)+0.125,format='mjd').iso+"\n")
    f.close()
    
    print "Number of days searched:    ", days
    print "Number of existing records: ", existingrecords
    print "Records Written:            ", recordswritten
    print "Missing Records:            ", days-(recordswritten+existingrecords)

def sql2rstn(t=None):
    """This function extracts the RSTN data from SQL with SQL_timestamp 
    0300 on the date supplied. If the values could be extracted then the
    data is returned in a list as follows:
    
    0 - timestamp: Astropy Time which is the date on which the data was
        collected. This should match dt
        
    1 - freq: A float32 numpy array containing the 9 frequencies in GHz
    
    2 - data: The flux data which is a 9x7 int16 numpy array.
    
    SQL_timestamp is also returned."""
    
    
    if t is None: t=Time.now()
    
    sqlt=Time(np.floor(t.mjd)+0.125,format='mjd')
    xml, buf = ch.read_cal(12, sqlt)
    
    if buf is None: return None, None
    
    sqlt_read=Time(extract(buf, xml['SQL_timestamp']), format='lv')
    #if np.floor(sqlt.mjd) != np.floor(sqlt_read.mjd): return None, None
    
    data=[]
    data.append(Time(extract(buf,xml['Timestamp']),format='lv'))
    data.append(extract(buf, xml['FGHz']))
    data.append(extract(buf, xml['Flux']))
    
    return data, sqlt_read
