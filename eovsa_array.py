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

    # Initial baseline corrections (based on Satellite obs. on 2016 Mar 20)
    dx = numpy.array([ 0.00, 0.08, 0.30, 0.67, 0.35, -0.13, -0.09, 0.94, -6.37, 6.51, 1.15,-12.50, 13.31,  0.0, 0.0, 0.0])
    dy = numpy.array([ 0.00,-0.45, 0.17,-0.39, 0.18, -0.79, -0.64, 1.47,-23.56,21.54,-2.50,-38.75, 65.38,  0.0, 0.0, 0.0])
    dz = numpy.array([ 0.00, 0.00, 0.00, 0.00, 0.00,  0.00,  0.00, 0.00,  0.00, 0.00, 0.00,  0.00,  0.00,  0.0, 0.0, 0.0])
    # Update based on Satellite obs. on 2016 Mar 29 (adds Ant14)
    dx += numpy.array([0.00, 0.61,-0.22, 0.14,-0.06, -0.19,  0.28,-0.19,  2.38,-1.84,-0.91,  4.58, -7.61,-2.44, 0.0, 0.0])
    dy += numpy.array([0.00, 0.22, 0.29, 0.44, 0.00,  0.38,  0.17,-0.29,  1.72,-2.97, 1.64,  5.10, -8.89, 3.56, 0.0, 0.0])
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

