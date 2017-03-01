#
# Defines the EOVSA antenna array, and also calculates Sun rise and set
# times for the array
#
# History:
#   2015-May-29  DG
#      Converted from using datime() to using Time() based on astropy.
#   2016-Mar-17  DG
#      Added bl_cor() routine to add baseline corrections in X, Y, Z
#   2016-Mar-20  DG
#      Updated antenna coordinates
#   2016-Mar-30  DG
#      Another update of antenna coordinates
#   2016-May-22  DG
#      Update of antenna coordinates based on 3C84 observations
#   2016-May-27  DG
#      Update of antenna coordinates based on 3C273 observations
#   2016-May-28  DG
#      First Bz correction of antenna coordinates, based on 3C273 and 3C286 observations
#   2016-Nov-18  DG
#      Updated ants 1-5, 9, 12 and 14 based on correction for offset-axes on new dishes.
#      Also updated ants 9,10,11,13 based on 3C84 obs later same day.
#   2016-Nov-20  DG
#      Finally!  Got updates for Bz for ants 1-5, 9-14.  Ant 11's 57 ns is just a guess
#   2016-Nov-22  DG
#      Ant 11 looks like a bad guess, but at least now it is measureable.  I reduced
#      it by 25.1 ns, to 31.9 ns.
#   2016-Nov-24  DG
#      A (final?) tweak to Bz.
#   2016-Nov-25  DG
#      Nope, not final!  I got the wrong sign on Bz correction.  Here is another try.
#   2017-Jan-06  DG
#      Another update of Bz, based on observations on Dec. 22.  Also update of Bx and By
#      based on today's observations on 3C273.
#

import aipy, ephem, numpy
from math import cos, sin
from util import Vector,Time
from numpy import pi, mat

global lat
lat = 37.233170*numpy.pi/180       # OVSA Latitude (radians)

def eovsa_array():
    ''' Define EOVSA antenna array, which consists of tabulated E,N,U
        locations of 15 antennas plus "array center", latitude,
        longitude and elevation of array center.  Returns AIPY
        AntennaArray object.
    '''
    global lat
    mperns = 0.299792458  # Meters per nanosecond

    # Define antenna ENU locations.  Ant 16 is the test input (0,0,0).
    # Divide by mperns to convert m to ns.
    ante = numpy.array([187.86, 196.15, 175.11, 197.96, 194.11, 147.42,
                        266.83, 98.95, 20.35, 167.43, -442.00, 640.22,
                        -329.06, -631.00, -213.00, 0.0])/mperns
    antn = numpy.array([71.74, 75.14, 77.39, 50.25, 108.86, 35.91,
                        67.10, 169.34, -218.49, 280.78, -138.59, -355.82,
                        861.82, -184.00, -187.00, 0.0])/mperns
    antu = numpy.zeros(16)

    lng = -118.286953*numpy.pi/180      # OVSA Longitude (radians)
    elev = 1207.0                       # OVSA Elevation (meters)
    clat = cos(lat)
    slat = sin(lat)

    # Latitude rotation matrix to convert ENU to XYZ
    latrot = mat([[0, -slat, clat],[1, 0, 0],[0, clat, slat]])
    f = numpy.array([1.0])
    beam = aipy.phs.Beam(f)
    ants = []
    for i in range(16):
        enu = Vector([ante[i], antn[i], antu[i]])
        x, y, z = enu.rotate(latrot).get()
        # Apply (add) any baseline corrections
        xp, yp, zp = bl_cor(x,y,z,i) 
        ants.append(aipy.phs.Antenna(xp,yp,zp,beam))

    aa = aipy.phs.AntennaArray(ants=ants,location=(lat, lng, elev))
    aa.horizon='10:30:00'
    aa.compute_pressure()
    # Create some standard sources for the source catalog
    srcs = []
    srcs.append(aipy.amp.RadioSpecial('Sun'))
    srcs.append(aipy.amp.RadioSpecial('Moon'))
    cat = aipy.amp.SrcCatalog(srcs)
    cat.compute(aa)
    # Attach catalog to aa object
    aa.cat = cat
    return aa

def bl_cor(x, y, z, iant):

    mperns = 0.299792458  # Meters per nanosecond

    # Initial baseline corrections (based on Satellite obs. on 2016 Mar 20)
    dx = numpy.array([ 0.00, 0.08, 0.30, 0.67, 0.35, -0.13, -0.09, 0.94, -6.37, 6.51, 1.15,-12.50, 13.31,  0.0, 0.0, 0.0])
    dy = numpy.array([ 0.00,-0.45, 0.17,-0.39, 0.18, -0.79, -0.64, 1.47,-23.56,21.54,-2.50,-38.75, 65.38,  0.0, 0.0, 0.0])
    dz = numpy.array([ 0.00, 0.00, 0.00, 0.00, 0.00,  0.00,  0.00, 0.00,  0.00, 0.00, 0.00,  0.00,  0.00,  0.0, 0.0, 0.0])
    # Update based on Satellite obs. on 2016 Mar 29 (adds Ant14)
    dx += numpy.array([0.00, 0.61,-0.22, 0.14,-0.06, -0.19,  0.28,-0.19,  2.38,-1.84,-0.91,  4.58, -7.61,-2.44, 0.0, 0.0])
    dy += numpy.array([0.00, 0.22, 0.29, 0.44, 0.00,  0.38,  0.17,-0.29,  1.72,-2.97, 1.64,  5.10, -8.89, 3.56, 0.0, 0.0])
    # Update based on 3C84 obs. on 2016 May 22 -- these are in m, hence the divieion by mperns
    #                   A1     2     3     4     5     6      7     8      9     10    11     12    13    14
    dx += numpy.array([0.00, 0.00, 0.00, 0.00, 0.00,  0.00,  0.00, 0.00,  0.00,-4.13,-3.13,  0.0, -3.75,-10.01, 0.0, 0.0])/mperns
    dy += numpy.array([0.00, 0.00, 0.00, 0.00, 0.00,  0.00,  0.00, 0.00,  0.00,-2.01,-2.01,  0.0, -3.94, -4.88, 0.0, 0.0])/mperns
    # Update based on 3C273 obs. on 2016 May 27 -- these are in ns
    dx += numpy.array([0.00,-1.20, 0.66,-0.19, 0.86,  1.33,  0.37,-0.42, 2.41, 0.36,-3.19, 25.91, -1.44,-0.69, 0.0, 0.0])
    dy += numpy.array([0.00, 0.19,-0.13, 0.05, 0.15,  0.39,  0.68,-0.90, 2.29, 0.22,-0.25, 14.35,  0.00, 0.15, 0.0, 0.0])
    # Update based on 3C273 and 3C286 obs. on 2016 May 28 -- again in ns
    dz += numpy.array([0.00,-0.01, 0.98,-0.15,-2.19, -0.10,  0.05, 0.05,19.08,11.80, 7.80,-15.20,  8.75, 6.80, 0.0, 0.0])
    # Update based on multiple sources on 2016 Nov 17, after axis-offset correction (ants 5-8 not in service) 
    # Further update 2016 Nov 18
    #                   A1     2      3      4      5      6     7    8     9      10    11     12     13     14
    dx += numpy.array([0.00, 0.018, 0.018,-0.035,-0.054,  0.0,  0.0, 0.0, 0.624, 0.388, 0.272,-0.132, 0.582, 0.353, 0.0, 0.0])
    dy += numpy.array([0.00,-0.018, 0.015, 0.036, 0.017,  0.0,  0.0, 0.0,-0.111,-0.035,-0.025,-0.097,-0.106,-0.035, 0.0, 0.0])    
    # Update based on 3C273 and 3C286 obs. on 2016 Nov 20 -- again in ns (ants 5-8 not in service)
    dz += numpy.array([0.00,-0.171,-0.258, 0.114,-0.403,  0.0,  0.0, 0.0, 1.323, 0.363, 31.90, 1.463, 0.798, 0.637, 0.0, 0.0])
    # Update based on 2136+006 and 2253+161 obs. on 2016 Nov 23, in ns (ants 5-8 not in service).
    dz += numpy.array([0.00, 0.034, 0.065,-0.106,-0.002,  0.0,  0.0, 0.0,-0.211, 0.101,-3.702,-1.730,-0.206,-0.168, 0.0, 0.0])
    # Update based on 2136+006 and 2253+161 obs. on 2016 Nov 25, in ns (ants 5-8 not in service).
    dz += numpy.array([0.00, 0.000,-0.085, 0.000,-0.061,  0.0,  0.0, 0.0, 0.680,-0.039, 7.696,-0.186, 0.607, 0.645, 0.0, 0.0])
    # Update based on obs. on 2016 Dec 22, in ns (ants 7-8 not in service).
    dz += numpy.array([0.00,-0.090, 0.010, 0.093, 0.010, 0.190, 0.0, 0.0,-0.394,-0.210,-0.437,-0.193,-0.357,-0.444, 0.0, 0.0])
    # Update based on obs. on 2017 Jan 06, in ns (ant 7 not in service).  Ant 13 not a good fit.
    dx += numpy.array([0.00, 0.017, 0.003,-0.037,-0.027, 0.050, 0.0, 0.160,-0.027,-0.057,-0.060,-0.010,-0.130,-0.027,0.0,0.0])
    dy += numpy.array([0.00, 0.00 ,-0.007, 0.00 , 0.00 ,-0.003, 0.0,-0.057, 0.020, 0.023, 0.027, 0.010,-0.010, 0.017,0.0,0.0])

    # Corrections are subtracted from nominal positions.
    xp = x - dx[iant]
    yp = y - dy[iant]
    zp = z - dz[iant]
    return xp, yp, zp
    
def suntimes(t=None,out=None):
    '''Returns the rise and set times of the Sun for EOVSA (10.5-degree
       horizon) for the day of the given Time() object, or today, 
       if not supplied.  The values are returned as a pair of mjd 
       values giving times of sunrise and sunset, by default, or 
       standard time strings if keyword out='str'.
    '''
    aa = eovsa_array()
    if t is None:
        t = Time()
    date = t.iso.replace('-','/')[:10]+' 20:00'
    risestr = str(aa.previous_rising(ephem.Sun(),start=date)).replace('/','-')
    setstr = str(aa.next_setting(ephem.Sun(),start=date)).replace('/','-')
    if out is 'str':
        return risestr, setstr
    trange = Time([risestr,setstr])

    return trange[0].mjd, trange[1].mjd

def sun_risetime(t,limit_deg=10):
    '''Calculate the times that the Sun rises above a given elevation
       limit for the OVRO site.

       t is a Time() object representing the day in question.
       limit_deg is the elevation limit in degrees (default 10 degrees)
       returns a Time() object for the sun risetime above limit_deg
    '''
    limit_rad = limit_deg*pi/180  
    aa = eovsa_array()
    s1 = aipy.phs.RadioSpecial('Sun')
    mjd_day = int(t.mjd)
    # Loop over hour of day
    for h in range(13,18):
        jultime = 2400000.5 + mjd_day + h/24.
        aa.set_jultime(jultime)
        s1.compute(aa)
        if s1.alt > limit_rad:
            h -= 1
            print s1.alt, h
            # Loop over minute of hour
            for m in range(60):
                jultime = 2400000.5 + mjd_day + (h + m/60.)/24.
                aa.set_jultime(jultime)
                s1.compute(aa)
                if s1.alt > limit_rad:
                    m -= 1
                    print s1.alt, h, m
                    # Loop over second of minute
                    for s in range(60):
                        jultime = 2400000.5 + mjd_day + (h + (m + s/60.)/60.)/24.
                        aa.set_jultime(jultime)
                        s1.compute(aa)
                        if s1.alt > limit_rad:
                            print s1.alt, h, m, s
                            risetime = Time(mjd_day + (h + (m + s/60.)/60.)/24.,format='mjd')
                            return risetime

def sun_settime(t,limit_deg=10):
    '''Calculate the times that the Sun sets below a given elevation
       limit for the OVRO site.

       t is a Time() object representing the day in question.
       limit_deg is the elevation limit in degrees (default 10 degrees)
       returns a Time() object for the sun settime below limit_deg
    '''
    limit_rad = limit_deg*pi/180  
    aa = eovsa_array()
    s1 = aipy.phs.RadioSpecial('Sun')
    mjd_day = int(t.mjd)
    # Loop over hour of day
    for h in range(22,28):
        jultime = 2400000.5 + mjd_day + h/24.
        aa.set_jultime(jultime)
        s1.compute(aa)
        if s1.alt < limit_rad:
            h -= 1
            # Loop over minute of hour
            for m in range(60):
                jultime = 2400000.5 + mjd_day + (h + m/60.)/24.
                aa.set_jultime(jultime)
                s1.compute(aa)
                if s1.alt < limit_rad:
                    m -= 1
                    # Loop over second of minute
                    for s in range(60):
                        jultime = 2400000.5 + mjd_day + (h + (m + s/60.)/60.)/24.
                        aa.set_jultime(jultime)
                        s1.compute(aa)
                        if s1.alt < limit_rad:
                            settime = Time(mjd_day + (h + (m + s/60.)/60.)/24.,format='mjd')
                            return settime

