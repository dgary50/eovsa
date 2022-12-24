# History:
#   2015-Feb-15  DG
#      Add an alternative make_geosattable() for handling Geosynchronous satellites
#      Instead of a source name, to look up in the aa.cat catalog, the saved 
#      ephem.EarthSatellite object is passed in.
#   2015-May-29  DG
#      Converted from using datime() to using Time() based on astropy.
#   2016-Dec-18  DG
#      Fixed problems with "wraps" in either RA or time, which cause grief to
#      some of the antennas.  Now puts 100-ms increments around such wrap points
#      to make the controller tracking algorithm unhappy only momentarily.
#   2017-Jan-05  DG
#      Tried to get make_tracktable work with GEOSATS, but gave up and
#      reinstated make_geosattable()
#   2022-Mar-20  DG
#      Today the 100-ms increments were too small for the Sun (it moves slowly
#      in RA).  Modified it to use 2-s increments when src.name == 'Sun'.
#
import aipy, ephem, numpy
import util
from numpy import pi,abs
from eovsa_array import *

def make_tracktable(srcname,aa,mjd1=None,mjd2=None,dt=1/24.):
    '''Generate a tracktable of coordinates of source with name srcname, which must be in 
       the src catalog aa.cat, with 1 h separation between lines, from the time indicated
       by mjd1 to 1 hour past the time indicated by mjd2.  If called without time arguments,
       a tracktable for the entire current day is generated.
       
       Default dt is 1 hour (dt is in units of days) but you can pass a different value when
       you call this function if you would like.
       
       This should work for Solar system and sidereal sources - not so sure yet about
       geostationary satellites.
    '''
    
    # make sure the source is in the source catalog aa.cat
    try:
        src = aa.cat[srcname]
    except:
        return 'Source ',srcname,' not found in catalog!  Tracktable not generated.'
    
    # if called without start and stop times, set start and stop times to now and 1 day from now
    if mjd1 is None:
        mjd1 = int(util.Time.now().mjd)
        mjd2 = mjd1 + 1.
    
    tbl = ''
    mjd = mjd1 - 2*dt # start tracktable a couple time intervals early for purposes of interpolation
    aa.date = mjd - dt - 15019.5
    src.compute(aa)
    raprev = src.ra
    mjdprev = mjd
    while mjd < (mjd2+4*dt): # and end it a few time intervals beyond mjd2 for same reason
        jultime = mjd + 2400000.5
        
        # precess the source's RA/dec coords to jultime
        aa.set_jultime(jultime)
        src.compute(aa)
        # Check if the abs(RA-RAprev) is greater than pi (indicates "wrap" in RA, 
        #                                          which 2.1 m dishes cannot tolerate)
        if abs(src.ra - raprev) > pi:
            # There was a wrap, so get time of transition
            mjd0 = mjd-dt # previous time
            mjd1 = mjd    # current time
            mjdz = (2*pi - raprev)*(mjd1 - mjd0)/((src.ra+2*pi)-raprev) + mjd0  # time of transition
            # Create lines around transition time
            for i in range(-2,3):
                if src.name == 'Sun':
                    mjd = mjdz + i*2/86400. # Sun moves slowly... -4, -2, 0, 2, 4 s
                else:
                    mjd = mjdz + i/864000. # -0.2, -0.1, 0, 0.1, 0.2 s
                aa.date = mjd - 15019.5
                src.compute(aa)
                msec = round((mjd % 1)*86400000.)        
                tbl = tbl+'{:7d} {:7d} {:5d} {:8d}\n'.format(int(src.ra*1800000./pi),
                                                     int(src.dec*1800000./pi),
                                                     int(mjd),int(msec))
            # Put back current time
            mjd = mjd1
            aa.date = mjd - 15019.5  # Convert mdj to ephem's strange time base
            src.compute(aa)
        raprev = src.ra

        # Do the same thing for a day change, since both Ant12 and Ant14 are too stupid
        # to do the right thing during a day change...
        if (int(mjd) - int(mjdprev)) != 0:
            # There was a day change, so get time of transition
            mjd0 = mjd-dt # previous time
            mjd1 = mjd    # current time
            mjdz = int(mjd1)  # time of transition
            # Create lines -2, -1, 0, 1, and 2 s around transition time
            for i in range(-2,3):
                mjd = float(mjdz) + i/864000.
                aa.date = mjd - 15019.5
                src.compute(aa)
                msec = round((mjd % 1)*86400000.)        
                tbl = tbl+'{:7d} {:7d} {:5d} {:8d}\n'.format(int(src.ra*1800000./pi),
                                                     int(src.dec*1800000./pi),
                                                     int(mjd),int(msec))
            # Put back current time
            mjd = mjd1
            aa.date = mjd - 15019.5  # Convert mdj to ephem's strange time base
            src.compute(aa)
        mjdprev = mjd
        
        # units for tracktable:
        # angles are in units of 1/10,000 of a degree
        # times (msec) are in units of ms
        msec = round((mjd % 1)*86400.)*1000.
        tbl = tbl+'{:7d} {:7d} {:5d} {:8d}\n'.format(int(src.ra*1800000./pi),
                                                     int(src.dec*1800000./pi),
                                                     int(mjd),int(msec))
        #print src.az,src.alt
        mjd = mjd + dt

    
    return tbl

def make_geosattable(sat,aa,mjd1=None,mjd2=None,dt=1/24.):
    '''Generate a tracktable of coordinates for geosynchronous satellite contained in
       the ephem.EarthSatellite object sat.  As for other objects, there is a dt separation 
       between lines, from 2 dt earlier than the time indicated by mjd1 to 3 dt past 
       the time indicated by mjd2.  If called without time arguments, a tracktable for the 
       entire current day is generated.
       
       Default dt is 1 hour (dt is in units of days) but you can pass a different value when
       you call this function if you would like.
    '''
    
    if sat is None:
        return 'No GEOSAT defined!  Tracktable not generated.'    

    # if called without start and stop times, set start and stop times to now and 1 day from now
    if mjd1 is None:
        mjd1 = int(util.Time.now().mjd)
        mjd2 = mjd1 + 1.
    
    tbl = ''
    mjd = mjd1 - 2*dt # start tracktable a couple time intervals early for purposes of interpolation
    aa.date = mjd - dt - 15019.5
    sat.compute(aa)
    raprev = sat.ra
    mjdprev = mjd
    while mjd < (mjd2+4*dt): # and end it a few time intervals beyond mjd2 for the same reason
        # Update aa's time
        aa.date = mjd - 15019.5  # Convert mdj to ephem's strange time base
        sat.compute(aa)

        # Check if the abs(RA-RAprev) is greater than pi (indicates "wrap" in RA, 
        #                                          which 2.1 m dishes cannot tolerate)
        if abs(sat.ra - raprev) > pi:
            # There was a wrap, so get time of transition
            mjd0 = mjd-dt # previous time
            mjd1 = mjd    # current time
            mjdz = (2*pi - raprev)*(mjd1 - mjd0)/((sat.ra+2*pi)-raprev) + mjd0  # time of transition
            # Create lines -2, -1, 0, 1, and 2 s around transition time
            for i in range(-2,3):
                mjd = mjdz + i/86400.
                aa.date = mjd - 15019.5
                sat.compute(aa)
                msec = round((mjd % 1)*86400.)*1000.        
                tbl = tbl+'{:7d} {:7d} {:5d} {:8d}\n'.format(int(sat.ra*1800000./pi),
                                                     int(sat.dec*1800000./pi),
                                                     int(mjd),int(msec))
            # Put back current time
            mjd = mjd1
            aa.date = mjd - 15019.5  # Convert mdj to ephem's strange time base
            sat.compute(aa)
        raprev = sat.ra

        # Do the same thing for a day change, since both Ant12 and Ant14 are too stupid
        # to do the right thing during a day change...
        if (int(mjd) - int(mjdprev)) != 0:
            # There was a day change, so get time of transition
            mjd0 = mjd-dt # previous time
            mjd1 = mjd    # current time
            mjdz = int(mjd1)  # time of transition
            # Create lines -2, -1, 0, 1, and 2 s around transition time
            for i in range(-2,3):
                mjd = mjdz + i/86400.
                aa.date = mjd - 15019.5
                sat.compute(aa)
                msec = round((mjd % 1)*86400.)*1000.        
                tbl = tbl+'{:7d} {:7d} {:5d} {:8d}\n'.format(int(sat.ra*1800000./pi),
                                                     int(sat.dec*1800000./pi),
                                                     int(mjd),int(msec))
            # Put back current time
            mjd = mjd1
            aa.date = mjd - 15019.5  # Convert mdj to ephem's strange time base
            sat.compute(aa)
        mjdprev = mjd

        # units for tracktable:
        # angles are in units of 1/10,000 of a degree
        # times (msec) are in units of ms
        msec = round((mjd % 1)*86400.)*1000.
        
        tbl = tbl+'{:7d} {:7d} {:5d} {:8d}\n'.format(int(sat.ra*1800000./pi),
                                                     int(sat.dec*1800000./pi),
                                                     int(mjd),int(msec))
        mjd = mjd + dt

    return tbl
