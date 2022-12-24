#
# Main routine for implementing the scan header.
# 
# History:
#   2014-Dec-30  DG
#     Started this history log.  Changed handling of chanmask due to
#     change in definition to a numpy array rather than simple list.
#   2015-Mar-31  JV
#     Changed datfile from being defined inside the scan_header function
#     to a passable parameter, so that Subarray2 can write to a different
#     file location. (Default is still '/tmp/scan_header.dat' if not specified.)
#   2015-May-29  DG
#     Converted from using datime() to using Time() based on astropy.
#   2015-Jun-16  DG
#     FTP to ACC now requires a username and password
#   2016-May-20  DG
#     Added a retry to tranferring files to ACC
#   2018-Jan-04  DG
#     Used correct default SK thresholds for 1792 samples, for flagging RFI
#   2018-Jan-14  DG
#     Set SK_MODE to 1, to enable automatic RFI flagging by dppxmp
#   2018-Mar-14  DG
#     Set SK_MODE back to 0, while we work on installing notch filters.
#   2019-Feb-22  DG
#     Import chan_util from new chan_util_52, which defines things for new
#     IF filters, e.g. 52 channels of 325 MHz bandwidth.
#   2019-Feb-23  DG
#     Several other changes for 52-channel mode.
#   2022-Mar-07  DG
#     Made some adjustments related to unavailability of SQL--default
#     file location is now on the /nas4 RAID disk.
#   2022-Mar-14  DG
#     Changes to use the new Chan_Info object defined in chan_info_52 to 
#     implement a fast FLARE mode.
#
import struct,sys
from sun_pos import *
from math import pi
import numpy as np
import util
import chan_info_52 as ci
from eovsa_array import *
from eovsa_lst import *
from ftplib import FTP

def scan_header(sh_dict,datfile='/common/Tables/scanheader/scan_header.dat'):
    '''Writes the state frame header file from the scan header dictionary 
       created by the schedule. Returns file names datfile and xmlfile 
       corresponding to the output files in the /tmp directory that are 
       created by this routine and are updated at the start of each scan.  
       The format string fmt can be used with struct.unpack() to read 
       the data file /tmp/scan_header.dat, although that usage is not 
       anticipated except for testing.

       This routine does something sensible even if the supplied
       dictionary sh_dict is empty (i.e. is {}).
    '''
    dtor = pi/180.
    xmlfile = datfile[:-4]+'.xml'
    f = open(datfile,'wb')
    xml = open(xmlfile,'w')
    xml.write('<Cluster>\n')
    xml.write('<Name>Scan_Header</Name>')
    xml.write('<NumElts>54</NumElts>')

    fmt = ''
    buf = ''
    
    # Scan Header Timestamp (double) [s, in LabVIEW format]
    # To be compatible with other timestamps in the stateframe, this
    # will be in LabVIEW format, which is s since 1904/01/01 (don't ask).
    item = sh_dict.get('timestamp',0.0)
    fmt = '<d'
    buf = struct.pack('d',item)
    f.write(buf)
    xml.write('<DBL>\n')
    xml.write('<Name>Timestamp</Name>\n')
    xml.write('<Val></Val>\n')
    xml.write('</DBL>\n')
    
    # Schedule version (double) [N/A]
    # Version of the schedule/scan header stateframe.
    item = sh_dict.get('Version',0.4)
    fmt += 'd'
    buf = struct.pack('d',item)
    f.write(buf)
    xml.write('<DBL>\n')
    xml.write('<Name>Version</Name>\n')
    xml.write('<Val>'+str(item)+'</Val>\n')
    xml.write('</DBL>\n')

    # Project (ascii, length 32 string)
    # Purpose of observations (correlator testing, routine observations)
    # Project is always the first thing added, so open the file
    # Default is 'NormalObserving'
    item = sh_dict.get('project','NormalObserving')
    fmt += 'I32s'
    buf = struct.pack('I32s',32,item)
    f.write(buf)
    xml.write('<Array>\n')
    xml.write('<Name>Project</Name>\n')
    xml.write('<Dimsize>32</Dimsize>\n<U8>\n<Name></Name>\n<Val></Val>\n</U8>\n')
    xml.write('</Array>\n')

    # Operator (ascii, length 16 string)
    # Default is 'Kjell Nelin'
    item = sh_dict.get('operator','Kjell Nelin')
    fmt += 'I16s'
    buf = struct.pack('I16s',16,item)
    f.write(buf)
    xml.write('<Array>\n')
    xml.write('<Name>Operator</Name>\n')
    xml.write('<Dimsize>16</Dimsize>\n<U8>\n<Name></Name>\n<Val></Val>\n</U8>\n')
    xml.write('</Array>\n')

    # Operator comments/log (ascii, length 120 string)
    # Default is 'None'
    item = sh_dict.get('comments','None')
    fmt += 'I120s'
    buf = struct.pack('I120s',120,item)
    f.write(buf)
    xml.write('<Array>\n')
    xml.write('<Name>Comments</Name>\n')
    xml.write('<Dimsize>120</Dimsize>\n<U8>\n<Name></Name>\n<Val></Val>\n</U8>\n')
    xml.write('</Array>\n')

    # Version (ascii, 2 length 8 strings)
    # Current hardware/software versions
    # Default is 1.0.0 and 1.0.0 
    item = sh_dict.get('version',['1.0.0','1.0.0'])
    fmt += 'I8sI8s'
    buf = struct.pack('I8sI8s',8,item[0],8,item[1])
    f.write(buf)
    xml.write('<Array>\n')
    xml.write('<Name>HVersion</Name>\n')
    xml.write('<Dimsize>8</Dimsize>\n<U8>\n<Name></Name>\n<Val></Val>\n</U8>\n')
    xml.write('</Array>\n')
    xml.write('<Array>\n')
    xml.write('<Name>SVersion</Name>\n')
    xml.write('<Dimsize>8</Dimsize>\n<U8>\n<Name></Name>\n<Val></Val>\n</U8>\n')
    xml.write('</Array>\n')

    # Nants (unsigned integer)
    # Number of active antennas
    # Default is 16
    item = sh_dict.get('nants',16)
    fmt += 'I'
    buf = struct.pack('I',item)
    f.write(buf)
    xml.write('<U32>')
    xml.write('<Name>Nants</Name>')
    xml.write('<Val></Val>')
    xml.write('</U32>')

    # Active antenna list (list of unsigned integers, length 16)
    item = sh_dict.get('antlist',range(1,17))
    nants = len(item)
    fmt += 'I16I'
    buf = struct.pack('I',16)
    for i in range(16):
        buf += struct.pack('I',item[i])
    f.write(buf)
    xml.write('<Array>\n')
    xml.write('<Name>Antlist</Name>\n')
    xml.write('<Dimsize>16</Dimsize>\n<U32>\n<Name></Name>\n<Val></Val>\n</U32>\n')
    xml.write('</Array>\n')

    # Antpos (3 x 16 FP array) [ns]
    # Antenna equatorial coordinates
    # Default is nominal EOVSA array positions
    item = sh_dict.get('antpos',eovsa_array())
    fmt += '2I'
    buf = struct.pack('2I',3,16)
    for i in range(16):
        p = item.ants[i].pos
        fmt += '3d'
        buf += struct.pack('3d',p[0],p[1],p[2])
    f.write(buf)
    xml.write('<Array>\n')
    xml.write('<Name>Antpos</Name>\n')
    xml.write('<Dimsize>3</Dimsize>\n<Dimsize>16</Dimsize><DBL>\n<Name></Name>\n<Val></Val>\n</DBL>\n')
    xml.write('</Array>\n')

    # Pbfwhm (16-element FP array) [arcsec]
    # Primary beam FWHM at 1 GHz
    # Default is nominal primary beam for 13 2.1m dishes, 2 27m dishes, and 0.0 (bare-feed)
    item = sh_dict.get('pbfwhm',[1.22*30*180*3600./210./pi]*13 + [1.22*30*180*3600./2700./pi]*2 + [0.0])
    fmt += 'I'
    buf = struct.pack('I',16)
    fmt += '16f'
    for i in range(16):
        buf += struct.pack('f',item[i])
    f.write(buf)
    xml.write('<Array>\n')
    xml.write('<Name>PrimaryBeam</Name>\n')
    xml.write('<Dimsize>16</Dimsize><SGL>\n<Name></Name>\n<Val></Val>\n</SGL>\n')
    xml.write('</Array>\n')

    # Mount (16-element signed integer array)
    # Type of antenna mounts
    # -1 = rf, 0 = bare-feed, 1 = 2m-azel, 2 = 27m-eq, 3 = 2m-eq
    # Default is list of 13 2m-azel, 2 27m-eq, and 1 bare-feed
    item = sh_dict.get('mount',[1]*13 + [2]*2 + [0])
    fmt += 'I16i'
    buf = struct.pack('I',16)
    for i in range(16):
        buf += struct.pack('i',item[i])
    f.write(buf)
    xml.write('<Array>\n')
    xml.write('<Name>Mount</Name>\n')
    xml.write('<Dimsize>16</Dimsize><I32>\n<Name></Name>\n<Val></Val>\n</I32>\n')
    xml.write('</Array>\n')

    # SCAN_ID (ascii, length 12) 
    # A unique string that identifies the scan (yymmddhhmmss)
    # Default is current date/time
    dt = util.Time.now()
    item = sh_dict.get('scan_id',dt.iso[2:19].replace('-','').replace(':','').replace(' ',''))
    fmt += 'I12s'
    buf = struct.pack('I12s',12,item)
    f.write(buf)
    xml.write('<Array>\n')
    xml.write('<Name>ScanID</Name>\n')
    xml.write('<Dimsize>12</Dimsize>\n<U8>\n<Name></Name>\n<Val></Val>\n</U8>\n')
    xml.write('</Array>\n')

    # SCAN_TYPE (ascii, length 12)
    # solar, calibration, etc
    # Default is 'test'
    item = sh_dict.get('scan_type','test')
    fmt += 'I12s'
    buf = struct.pack('I12s',12,item)
    f.write(buf)
    xml.write('<Array>\n')
    xml.write('<Name>ScanType</Name>\n')
    xml.write('<Dimsize>12</Dimsize>\n<U8>\n<Name></Name>\n<Val></Val>\n</U8>\n')
    xml.write('</Array>\n')

    # SOURCE_ID (ascii, length 12)
    # Name of source (e.g. 3C84, Sun, 0321+123, etc.)
    # Default is 'None'
    item = sh_dict.get('source_id','None')
    fmt += 'I12s'
    buf = struct.pack('I12s',12,item)
    f.write(buf)
    xml.write('<Array>\n')
    xml.write('<Name>SourceID</Name>\n')
    xml.write('<Dimsize>12</Dimsize>\n<U8>\n<Name></Name>\n<Val></Val>\n</U8>\n')
    xml.write('</Array>\n')

    # Coordinate type (ascii, length 6)
    # 'RADEC ' => fixed, use RA, Dec;
    # 'FIXED ' => geostationary, use HA, Dec;
    # 'PLANET' => planetary, use EPHEM
    # 'SATELL' => satellite tracking, use SAT Ephem
    # Default is 'PLANET'
    item = sh_dict.get('track_mode','PLANET')
    fmt += 'I6s'
    buf = struct.pack('I6s',6,item)
    f.write(buf)
    xml.write('<Array>\n')
    xml.write('<Name>TrackMode</Name>\n')
    xml.write('<Dimsize>6</Dimsize>\n<U8>\n<Name></Name>\n<Val></Val>\n</U8>\n')
    xml.write('</Array>\n')

    # epoch of source coordinates (ascii, length 4)
    # '2000' => J2000
    # 'DATE' => current date/time
    item = sh_dict.get('epoch','DATE')
    fmt += 'I4s'
    buf = struct.pack('I4s',4,item)
    f.write(buf)
    xml.write('<Array>\n')
    xml.write('<Name>Epoch</Name>\n')
    xml.write('<Dimsize>4</Dimsize>\n<U8>\n<Name></Name>\n<Val></Val>\n</U8>\n')
    xml.write('</Array>\n')

    # Right Ascension (double) [radians]
    # Default is current LST
    item = sh_dict.get('ra',eovsa_lst(dt))
    fmt += 'd'
    buf = struct.pack('d',item)
    f.write(buf)
    xml.write('<DBL>\n')
    xml.write('<Name>RA</Name>\n')
    xml.write('<Val></Val>\n')
    xml.write('</DBL>\n')

    # Declination (double) [radians]
    # Default is 0 Dec
    item = sh_dict.get('dec',0.0)
    fmt += 'd'
    buf = struct.pack('d',item)
    f.write(buf)
    xml.write('<DBL>\n')
    xml.write('<Name>Dec</Name>\n')
    xml.write('<Val></Val>\n')
    xml.write('</DBL>\n')

    # Right Ascension Offset of phase center (double) [radians]
    # Default is 0
    item = sh_dict.get('dra',0.0)
    fmt += 'd'
    buf = struct.pack('d',item)
    f.write(buf)
    xml.write('<DBL>\n')
    xml.write('<Name>dRA</Name>\n')
    xml.write('<Val></Val>\n')
    xml.write('</DBL>\n')

    # Declination Offset of phase center (double) [radians]
    # Default is 0
    item = sh_dict.get('ddec',0.0)
    fmt += 'd'
    buf = struct.pack('d',item)
    f.write(buf)
    xml.write('<DBL>\n')
    xml.write('<Name>dDec</Name>\n')
    xml.write('<Val></Val>\n')
    xml.write('</DBL>\n')

    # Hour angle for geostationary sources (double) [radians]
    # Default is 0
    item = sh_dict.get('ha',0.0)
    fmt += 'd'
    buf = struct.pack('d',item)
    f.write(buf)
    xml.write('<DBL>\n')
    xml.write('<Name>HA</Name>\n')
    xml.write('<Val></Val>\n')
    xml.write('</DBL>\n')

    # Hour angle offset for geostationary sources (double) [radians]
    # Default is 0
    item = sh_dict.get('dha',0.0)
    fmt += 'd'
    buf = struct.pack('d',item)
    f.write(buf)
    xml.write('<DBL>\n')
    xml.write('<Name>dHA</Name>\n')
    xml.write('<Val></Val>\n')
    xml.write('</DBL>\n')

    aa = eovsa_array()
    aa.date = str(aa.date)[:11]+'00:00'  # Set date to 0 UT
    mjd0 = aa.date + 15019.5

    # Ephemeris for planetary source ([T, RA, Dec] x 3, double) [mjd, radians, radians]
    # Coordinates are geocentric, times are UTC, and apparent coords must be calculated from them
    # Default is the ephemeris for the Sun for the current date
    item = sh_dict.get('ephem')
    if item is None:
        # Generate default ephemeris, which is that for the Sun for current date
        item = []
        aa = eovsa_array()
        aa.date = str(aa.date)[:11]+'00:00'  # Set date to 0 UT
        mjd0 = aa.date + 15019.5
        sun = ephem.Sun()
        for i in range(3):
            sun.compute(aa)
            mjd = aa.date + 15019.5
            item.append([mjd, sun.g_ra, sun.g_dec])  # Geocentric coordinates
            aa.date += 1
    fmt += 'II' + '3d'*3
    buf = struct.pack('2I',3,3)
    for i in range(3):
        for j in range(3):
            buf += struct.pack('d',item[i][j])
    f.write(buf)
    xml.write('<Array>\n')
    xml.write('<Name>Ephem</Name>\n')
    xml.write('<Dimsize>3</Dimsize>\n<Dimsize>3</Dimsize><DBL>\n<Name></Name>\n<Val></Val>\n</DBL>\n')
    xml.write('</Array>\n')

    # Calculate solar P, B0 and R (FP x 3) [radians, radians, arcsec]
    # Item is mjd of beginning of scan
    # Default is to calculate P, B0 and R for the current date
    item = sh_dict.get('sun_info',dt.mjd)
    fmt += 'Ifff'
    p, b0, r = get_pb0r(item,arcsec=True)
    buf = struct.pack('If',3,p*dtor)
    buf += struct.pack('f',b0*dtor)
    buf += struct.pack('f',r)
    f.write(buf)
    xml.write('<Array>\n')
    xml.write('<Name>SunInfo</Name>\n')
    xml.write('<Dimsize>3</Dimsize>\n<SGL>\n<Name></Name>\n<Val></Val>\n</SGL>\n')
    xml.write('</Array>\n')

    # Satellite ephemeris data ([T, RA, Dec] x nlines, double) [mjd, radians, radians]
    # This is tricky, since we want a fixed stateframe length, so we
    # have to specify as many lines as the maximum we might use.  For
    # now, set to 20 lines.
    # Default is all zero
    item = sh_dict.get('sat_ephem',[[0.0,0.0,0.0]]*20)
    nlines = len(item)
    fmt += 'II' + '3d'*20
    buf = struct.pack('2I',3,20)
    for i in range(nlines):
        for j in range(3):
            buf += struct.pack('d',item[i][j])
    if nlines < 20:
        for i in range(nlines,20):
            for j in range(3):
                buf += struct.pack('d',0.0)            
    f.write(buf)
    xml.write('<Array>\n')
    xml.write('<Name>SatEphem</Name>\n')
    xml.write('<Dimsize>3</Dimsize>\n<Dimsize>20</Dimsize><DBL>\n<Name></Name>\n<Val></Val>\n</DBL>\n')
    xml.write('</Array>\n')
    
    # Polarization list (Miriad definition) (signed int)
    #     1: Stokes I
    #     2: Stokes Q
    #     3: Stokes U
    #     4: Stokes V
    #    -1: Circular RR
    #    -2: Circular LL
    #    -3: Circular RL
    #    -4: Circular LR
    #    -5: Linear XX
    #    -6: Linear YY
    #    -7: Linear XY
    #    -8: Linear YX
    #     0: Not used
    # Zero-filled list of 4
    # Default is RR, LL
    item = sh_dict.get('pol',[-5,-6,-7,-8])
    fmt += '4i'
    buf = ''
    for i in range(4):
       buf += struct.pack('i',item[i])
    f.write(buf)
    xml.write('<I32>\n')
    xml.write('<Name>Pol1</Name>\n')
    xml.write('<Val></Val>\n')
    xml.write('</I32>\n')
    xml.write('<I32>\n')
    xml.write('<Name>Pol2</Name>\n')
    xml.write('<Val></Val>\n')
    xml.write('</I32>\n')
    xml.write('<I32>\n')
    xml.write('<Name>Pol3</Name>\n')
    xml.write('<Val></Val>\n')
    xml.write('</I32>\n')
    xml.write('<I32>\n')
    xml.write('<Name>Pol4</Name>\n')
    xml.write('<Val></Val>\n')
    xml.write('</I32>\n')

    # Timing offset (UT1-UTC) (double) [fraction of day]
    # Default is to read IERS bulletin and use (UT1-UTC) for current date
    item = sh_dict.get('ut1-utc',util.UT1_UTC(mjd0))
    fmt += 'd'
    buf = struct.pack('d',item)
    f.write(buf)
    xml.write('<DBL>\n')
    xml.write('<Name>UT1UTC</Name>\n')
    xml.write('<Val></Val>\n')
    xml.write('</DBL>\n')

    # Maximum IDB filesize (unsigned int) [MB]
    # DPP will break IDB scan data into separate files none 
    # of which will exceed this size
    # Default is 100 (MB)
    item = sh_dict.get('max_file_size',100)
    fmt += 'I'
    buf = struct.pack('I',item)
    f.write(buf)
    xml.write('<U32>\n')
    xml.write('<Name>IDBMaxFileSize</Name>\n')
    xml.write('<Val></Val>\n')
    xml.write('</U32>\n')

    # IDB filename stem (length 14 string of form yyyymmddhhmmss)
    # Value of item passed in is Time() object
    # Creates IDB filename stem value from Time() object.  
    # Filename will be extended as needed to accommodate multiple files per scan.
    # (Note: truncated at integer second, since all times should happen at
    # 1pps boundary)
    # Default is current date/time
    item = sh_dict.get('date2IDB_stem',dt)
    fmt += 'I14s'
    datestr = item.iso.replace('-','').replace(' ','').replace(':','')
    buf = struct.pack('I14s',14,datestr[:14])
    f.write(buf)
    xml.write('<Array>\n')
    xml.write('<Name>IDBStem</Name>\n')
    xml.write('<Dimsize>14</Dimsize>\n<U8>\n<Name></Name>\n<Val></Val>\n</U8>\n')
    xml.write('</Array>\n')

    # Nominal start time of scan (length 14 string of form yyyymmddhhmmss)
    # (Note: truncated at integer second, since all times should happen at
    # 1pps boundary)
    fmt += 'I14s'
    buf = struct.pack('I14s',14,datestr[:14])
    f.write(buf)
    xml.write('<Array>\n')
    xml.write('<Name>IDBStart</Name>\n')
    xml.write('<Dimsize>14</Dimsize>\n<U8>\n<Name></Name>\n<Val></Val>\n</U8>\n')
    xml.write('</Array>\n')

    # Time conversion factor (double) [mjd] 
    # (i.e., Relates the absolute time to the accumulation number for the 
    # packets).  This is the time of the next 1 pps after the ARM signal to 
    # the correlator.
    # Default is current date/time, which is an error
    item = sh_dict.get('time_at_acc0',dt)
    fmt += 'd'
    buf = struct.pack('d',item.mjd)
    f.write(buf)
    xml.write('<DBL>\n')
    xml.write('<Name>TimeAtAcc0</Name>\n')
    xml.write('<Val></Val>\n')
    xml.write('</DBL>\n')

    # Duration of each spectral frame. Nominally 1 s. (unsigned int) [ms]
    # Default is 1000 [ms]
    item = sh_dict.get('dur_spec_frame',1000)
    fmt += 'I'
    buf = struct.pack('I',item)
    f.write(buf)
    xml.write('<U32>\n')
    xml.write('<Name>DurSpecFrame</Name>\n')
    xml.write('<Val></Val>\n')
    xml.write('</U32>\n')

    # Intval Integration time (unsigned int) [ms]
    # Default is 20 [ms]
    item = sh_dict.get('intval',20)
    fmt += 'I'
    buf = struct.pack('I',item)
    f.write(buf)
    xml.write('<U32>\n')
    xml.write('<Name>Intval</Name>\n')
    xml.write('<Val></Val>\n')
    xml.write('</U32>\n')

    # Number of integration intervals in each spectral frame. (unsigned int)
    # Default is 50 [per second].
    item = sh_dict.get('nintval',50)
    fmt += 'I'
    buf = struct.pack('I',item)
    f.write(buf)
    xml.write('<U32>\n')
    xml.write('<Name>NIntval</Name>\n')
    xml.write('<Val></Val>\n')
    xml.write('</U32>\n')

    # List of frequencies of 600 MHz band start AFTER reversal (length 50 double) [GHz]
    # Input is integer band number corresponding to base frequency for
    # each integration interval, where 1 = 1 GHz, 2 = 1.5 GHz, etc.
    # (50 unsigned ints). These are the 500 MHz IF bands.  The integer
    # values n are converted to GHz as follows:
    #     21 + n/2.0 is LO frequency [GHz].
    #     LO frequency - 19.95 is frequency at low end of 600 MHz band
    #     After reversal of bins, subtract another 0.6 GHz
    #     Complete expression is equivalent to 0.45 + n/2.
    # The Hittite tuning list is passed in and parsed with comma ','
    # Default is nominal solar sequence
    #*******This needs modification for 800 MHz clock*******
    item = sh_dict.get('fsequence',
                        '1, 2, 3, 4, 5, 6, 7, 8, 9,10,'+
                        '1, 2, 3, 4,11,12,13,14,15,16,'+
                        '1, 2, 3, 4,17,18,19,20,21,22,'+
                        '1, 2, 3, 4,23,24,25,26,27,28,'+
                        '1, 2, 3, 4,29,30,31,32,33,34')
    # print item
    fseqlist = item.rsplit(',')
    nintval = len(fseqlist)
    fmt += 'I50d'
    buf = struct.pack('I',50)
    for n in fseqlist:
        val = 0.775 + float(n)*0.325
        buf += struct.pack('d',val)
    if nintval < 50:
        # Fewer than 50 values sent in, so zero-fill to 50
        for n in range(nintval,50):
            buf += struct.pack('d',0.0)
    f.write(buf)
    xml.write('<Array>\n')
    xml.write('<Name>FSeqList</Name>\n')
    xml.write('<Dimsize>50</Dimsize><DBL>\n<Name></Name>\n<Val></Val>\n</DBL>\n')
    xml.write('</Array>\n')

    # Subband channel width [GHz]
    # This is fixed by the system to 0.4 GHz/4096.
    # Default is 0.4/4096
    item = sh_dict.get('subbw',0.4/4096)
    fmt += 'd'
    buf = struct.pack('d',item)
    f.write(buf)
    xml.write('<DBL>\n')
    xml.write('<Name>SubBW</Name>\n')
    xml.write('<Val></Val>\n')
    xml.write('</DBL>\n')

    chinfo = sh_dict.get('chinfo',ci.Chan_Info())
    # Number of 'wide' spectral channels (unsigned int)
    # Default 504, max 511
    item = sh_dict.get('nchan',chinfo.tot_scichan())
    fmt += 'I'
    buf = struct.pack('I',item)
    f.write(buf)
    xml.write('<U32>\n')
    xml.write('<Name>Nchan</Name>\n')
    xml.write('<Val></Val>\n')
    xml.write('</U32>\n')

    # Wide spectral channel lower value--nominally 504 values (nchan doubles) [GHz]
    # Default is the default set of 504 frequencies valid for solar observing 
    # (zero-filled to 511)
    item = sh_dict.get('fGHz')
    if item is None:
        item = []
        for band in range(1,53):
            item += chinfo.start_freq(band)
    nchan = len(item)
    fmt += 'I511d'
    buf = struct.pack('I',511)
    for i in item:
        buf += struct.pack('d',i)
    if nchan < 511:
        # Fewer than 511 values sent in, so zero-fill to 511
        for n in range(nchan,511):
            buf += struct.pack('d',0.0)        
    f.write(buf)
    xml.write('<Array>\n')
    xml.write('<Name>fGHz</Name>\n')
    xml.write('<Dimsize>511</Dimsize><DBL>\n<Name></Name>\n<Val></Val>\n</DBL>\n')
    xml.write('</Array>\n')

    # Wide spectral channel widths--nominally 504 values (nchan doubles) [GHz]
    # Default is the list of widths valid for solar observing (zero-filled to 511)
    item = sh_dict.get('chan_widths')
    if item is None:
        item = []
        for band in range(1,53):
            item += chinfo.sci_bw(band)
    nchan = len(item)
    fmt += 'I511d'
    buf = struct.pack('I',511)
    for i in item:
        buf += struct.pack('d',i)
    if nchan < 511:
        # Fewer than 511 values sent in, so zero-fill to 511
        for n in range(nchan,511):
            buf += struct.pack('d',0.0)        
    f.write(buf)
    xml.write('<Array>\n')
    xml.write('<Name>ChanWidths</Name>\n')
    xml.write('<Dimsize>511</Dimsize><DBL>\n<Name></Name>\n<Val></Val>\n</DBL>\n')
    xml.write('</Array>\n')
    
    # Science channel assignments (4096 x 50 array of 2-byte integers)
    #   0             = subband is not used.  
    #   (1-500)       = science channel to assign subband.
    #   512 + (1-500) = integration should not be time-averaged.
    # Multiple samples <512 within a spectral frame are averaged.
    # Multiple samples >512 are not averaged and are treated
    # as separate channels.)
    # Restrictions: 
    #   1. There are a maximum of 511 science channels.
    #   2. A given subband cannot be associated with more than one
    #      science channel. 
    #   3. All subbands associated with a given science channel must
    #      be in the same 20 ms integration.
    item = sh_dict.get('chan2wide')
    if item is None:
        # Create default channel assignments for the list of frequencies (bands) in fseqlist
        item = []
        for band in fseqlist:
            ch = chinfo.chan_asmt(int(band))
            #ch[0] = 0  # Assign first subband to purgatory (unused subband)
            item += ch
    nsubchan = len(item)
    fmt += 'II204800H'
    buf = struct.pack('2I',4096,50)
    for i in item:
        buf += struct.pack('H',i)
    if nsubchan < 204800:
        # Fewer than 4096*50 values sent in, so zero-fill to 204800
        for n in range(nsubchan,204800):
            buf += struct.pack('H',0)        
    f.write(buf)
    xml.write('<Array>\n')
    xml.write('<Name>Chan2Wide</Name>\n')
    xml.write('<Dimsize>4096</Dimsize><Dimsize>50</Dimsize><U16>\n<Name></Name>\n<Val></Val>\n</U16>\n')
    xml.write('</Array>\n')

    # Mask of RFI subbands (4096 x 50 array of 1-byte flags).   
    #    0 = Subband is presumed contaminated
    #    1 = otherwise
    # Default is all 1 (good values)
    item = sh_dict.get('chanmask',np.array([1]*204800,'byte'))
    nsubchan = len(item)
    fmt += 'II204800B'
    buf = struct.pack('2I',4096,50)
    buf += struct.pack(str(nsubchan)+'B',*item)
    if nsubchan < 204800:
        # Fewer than 4096*50 values sent in, so zero-fill to 204800
        for n in range(nsubchan,204800):
            buf += struct.pack('B',0)        
    f.write(buf)
    xml.write('<Array>\n')
    xml.write('<Name>Chanmask</Name>\n')
    xml.write('<Dimsize>4096</Dimsize><Dimsize>50</Dimsize><B8>\n<Name></Name>\n<Val></Val>\n</B8>\n')
    xml.write('</Array>\n')

    # RFI guard band (4-byte int)
    # Number of subbands to reject adjacent to kurtosis-flagged channels.
    #    0 = do not reject any additional subbands
    # Default is 0
    item = sh_dict.get('sk_guard_width',0)
    fmt += 'I'
    buf = struct.pack('I',item)
    f.write(buf)
    xml.write('<U32>\n')
    xml.write('<Name>SKGuard</Name>\n')
    xml.write('<Val></Val>\n')
    xml.write('</U32>\n')

    # Spectral kurtosis strategy (4-byte int) (Details TBD.)
    # Default 0 => Standard strategy
    item = sh_dict.get('sk_mode',0)
    fmt += 'I'
    buf = struct.pack('I',item)
    f.write(buf)
    xml.write('<U32>\n')
    xml.write('<Name>SKMode</Name>\n')
    xml.write('<Val></Val>\n')
    xml.write('</U32>\n')

    # Kurtosis lower and upper limits (two doubles)
    # Min, max acceptable SK thresholds
    # Default is 0.87308624388667777, 1.1564648485565636, which are for M = 1792 
    item = sh_dict.get('sk_lims',[0.87308624388667777,1.1564648485565636])
    fmt += '2d'
    buf = struct.pack('2d',item[0],item[1])
    f.write(buf)
    xml.write('<DBL>\n')
    xml.write('<Name>SKLimLo</Name>\n')
    xml.write('<Val></Val>\n')
    xml.write('</DBL>\n')
    xml.write('<DBL>\n')
    xml.write('<Name>SKLimHi</Name>\n')
    xml.write('<Val></Val>\n')
    xml.write('</DBL>\n')

    # Reference complex gains (511 x 16 x 2 complex array) [Jy/unit]
    # These are time-independent values that correspond to the
    # 'base attenuator state'  which includes the effects of the
    # preset front-end and back-end attenuators. (Amplitude units
    # are nominally Jy, but may be modified by attenuator and
    # time-variable gain factors)
    # Note that Struct has no complex format, so buffer is written
    # with each complex value as a pair of floats, real,imag.
    # Default 1 + 0j
    item = sh_dict.get('gains',[1.0+0j]*16352)
    fmt += 'IIII16352d'
    buf = struct.pack('4I',2,511,16,2)
    for i in item:
        buf += struct.pack('2d',np.real(i),np.imag(i))
    f.write(buf)
    xml.write('<Array>\n')
    xml.write('<Name>GainTable</Name>\n')
    xml.write('<Dimsize>2</Dimsize><Dimsize>511</Dimsize><Dimsize>16</Dimsize><Dimsize>2</Dimsize><DBL>\n<Name></Name>\n<Val></Val>\n</DBL>\n')
    xml.write('</Array>\n')

    # Time variable gain factors (511 x 16 x 2 complex array) [Jy]
    # Antenna-based complex gains corresponding to time-variable
    # factors (presumably determined by a recent phase calibration).
    # Defaults to 1 + 0j if not present.
    item = sh_dict.get('gain_update',[1.0+0j]*16352)
    fmt += 'IIII16352d'
    buf = struct.pack('4I',2,511,16,2)
    for i in item:
        buf += struct.pack('2d',np.real(i),np.imag(i))
    f.write(buf)
    xml.write('<Array>\n')
    xml.write('<Name>GainUpdate</Name>\n')
    xml.write('<Dimsize>2</Dimsize><Dimsize>511</Dimsize><Dimsize>16</Dimsize><Dimsize>2</Dimsize><DBL>\n<Name></Name>\n<Val></Val>\n</DBL>\n')
    xml.write('</Array>\n')

    # Base attenuator value (4-byte int) [dB]
    # attenuator value corresponding to the base attenuator state
    # Default 0
    item = sh_dict.get('attn0',0)
    fmt += 'I'
    buf = struct.pack('I',item)
    f.write(buf)
    xml.write('<U32>\n')
    xml.write('<Name>Attn0</Name>\n')
    xml.write('<Val></Val>\n')
    xml.write('</U32>\n')

    # Attenuator step size (4-byte int) [dB]
    # Nominal step size for time-dependent variable attenuation 
    # 2-byte integer (presumably = 3 or 5 db)
    # Default 3 [dB]
    item = sh_dict.get('attn_step',3)
    fmt += 'I'
    buf = struct.pack('I',item)
    f.write(buf)
    xml.write('<U32>\n')
    xml.write('<Name>AttnStep</Name>\n')
    xml.write('<Val></Val>\n')
    xml.write('</U32>\n')

    # Attenuator gain table (500 x 16 x 2 x N complex array) 
    # complex gains (one for each science channel, polarization and antenna).
    # The values are relative to the 'base' attenuator settings, which represents
    # both the front end and back end attenuation.  Each of the N steps
    # corresponds to 3 or 5 db increments applied to all antennas.
    #item = sh_dict.get('attn_table')

    # XML encapsulation of KatADC info (array of 8)
    xml.write('<Array>\n')
    xml.write('<Name>KatADC</Name>\n')
    xml.write('<Dimsize>8</Dimsize>\n')
    # Put array dimension into data
    fmt += 'I'
    buf = struct.pack('I',8)

    # This is an array, but the definitions only go in once
    # Here are the definitions for the 8 elements in the cluster
    xml.write('<Cluster>\n')
    xml.write('<Name/>\n')
    xml.write('<NumElts>6</NumElts>\n')
    xml.write('<U32>\n')
    xml.write('<Name>Status</Name>\n')
    xml.write('<Val></Val>\n')
    xml.write('</U32>\n')
    xml.write('<SGL>\n')
    xml.write('<Name>Temp.adc0</Name>\n')
    xml.write('<Val></Val>\n')
    xml.write('</SGL>\n')
    xml.write('<SGL>\n')
    xml.write('<Name>Temp.adc1</Name>\n')
    xml.write('<Val></Val>\n')
    xml.write('</SGL>\n')
    xml.write('<SGL>\n')
    xml.write('<Name>Temp.ambient0</Name>\n')
    xml.write('<Val></Val>\n')
    xml.write('</SGL>\n')
    xml.write('<SGL>\n')
    xml.write('<Name>Temp.ambient1</Name>\n')
    xml.write('<Val></Val>\n')
    xml.write('</SGL>\n')
    xml.write('<SGL>\n')
    xml.write('<Name>BoardClock</Name>\n')
    xml.write('<Val></Val>\n')
    xml.write('</SGL>\n')
    # End KatADC Cluster
    xml.write('</Cluster>\n')

    # Get the KatADC array of dictionaries (defaults to eight empty dicts)
    katadc = sh_dict.get('katadc',[{},{},{},{},{},{},{},{}])
    brd_clk = sh_dict.get('roach_brd_clk',[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0])
    # The actual data has to go in eight times
    for i in range(8):
        # Handle case of empty dictionary
        if katadc[i] == {}:
            fmt += 'Iffff'
            buf += struct.pack('I',0)
            for j in range(4):
                buf += struct.pack('f',0.0)
        
        else:
            # Status (int) [bit list]
            # Each bit signifies whether corresponding sensor is nominal [0] or in error [1]
            # Order is alphabetical, msb to lsb.
            # Default is zero
            status = 0
            keys = sorted(katadc[i].keys())
            for key in keys:
                if key.find('status') != -1:
                    if katadc[i][key] == 'nominal':
                        status = status<<1
                    else:
                        status = (status<<1) + 1
            item = status
            fmt += 'I'
            buf += struct.pack('I',item)

            # Just add the sorted list of sensors
            # Default is zero
            for key in keys:
                if key.find('status') == -1:
                    item = katadc[i][key]
                    fmt += 'f'
                    buf += struct.pack('f',item)

        fmt += 'f'
        buf += struct.pack('f',brd_clk[i])
        
    # End KatADC Array
    xml.write('</Array>\n')

    f.write(buf)
    f.close()
    xml.write('</Cluster>\n')
    xml.close()

    # Connect to ACC /parm directory and transfer scan_header files
    try:
        acc = FTP('acc.solar.pvt')
        acc.login('admin','observer')
        acc.cwd('parm')
        # Send XML file to ACC
        f = open(xmlfile,'r')
        acc.storlines('STOR scan_header.xml',f)
        f.close()
        # Send DAT file to ACC    
        g = open(datfile,'rb')
        acc.storbinary('STOR '+datfile[5:],g)
        acc.close()
        g.close()
        sys.stdout.write('Successfully transferred scan_header files to ACC.\n')
    except:
        sys.stdout.write('Transfer of scan_header files failed.  Retrying...\n')
        try:
            acc = FTP('acc.solar.pvt')
            acc.login('admin','observer')
            acc.cwd('parm')
            # Send XML file to ACC
            f = open(xmlfile,'r')
            acc.storlines('STOR scan_header.xml',f)
            f.close()
            # Send DAT file to ACC    
            g = open(datfile,'rb')
            acc.storbinary('STOR scan_header.dat',g)
            acc.close()
            g.close()
            sys.stdout.write('Successfully transferred scan_header files to ACC.\n')
        except:
            sys.stdout.write('Could not transfer scan_header files.  ACC is down?\n')

    return fmt, datfile, xmlfile

