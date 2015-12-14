import aipy, ephem, numpy
import util
from math import cos, sin
from numpy import pi
from eovsa_array import *
from readvla import readvlacaldb
import urllib2

def solar_tracktable(mjd1=None,mjd2=None):
    '''Generate a solar tracktable with 1 h separation of solar coordinates
       between lines, from the time indicated by mjd1 to 1 hour past the time
       indicated by mjd2.  If called without arguments, a tracktable for the
       entire current day is generated.
    '''
    if mjd1 is None:
        d = util.datime()
        mjd1 = int(d.get())
        mjd2 = mjd1 + 1.
    aa = eovsa_array()
    aa.compute_pressure()
    srcs = aipy.phs.RadioSpecial('Sun')
    cat = aipy.phs.SrcCatalog(srcs)
    tbl = ''
    dt = 1/24.
    mjd = mjd1 - 2*dt
    while mjd < (mjd2+4*dt):
        jultime = mjd + 2400000.5
        aa.set_jultime(jultime)
        cat.compute(aa)
        msec = round((mjd % 1)*86400.)*1000.
        tbl = tbl+'{:7d} {:7d} {:5d} {:8d}\n'.format(int(cat['Sun'].ra*1800000./pi),
                                                     int(cat['Sun'].dec*1800000./pi),
                                                     int(mjd),int(msec))
        mjd = mjd + dt

    return tbl

def calibrator_tracktable(srcname,mjd1=None,mjd2=None):
    '''Generate a tracktable for the given calibrator source (quasar).
       This is just three lines giving the precessed RA, Dec coordinates.
    '''
    if mjd1 is None:
        d = util.datime()
        mjd1 = int(d.get())
        mjd2 = mjd1 + 1.

    cal = readvlacaldb()
    for c in cal:
        if c.name.find(srcname) == 0:
            break

    if srcname in c.name:
        aa = eovsa_array()
        aa.compute_pressure()
        src = aipy.phs.RadioFixedBody(c.ra[0],c.dec[0],name=c.name,epoch='2000')
        cat = aipy.phs.SrcCatalog(src)
        tbl = ''
        dt = 0.5/24.
        mjd = mjd1 - 2*dt
        while mjd < (mjd2+4*dt):
            jultime = mjd + 2400000.5
            aa.set_jultime(jultime)
            cat.compute(aa)
            msec = round((mjd % 1)*86400.)*1000.
            tbl = tbl+'{:7d} {:7d} {:5d} {:8d}\n'.format(int(cat[c.name].ra*1800000./pi),
                                                 int(cat[c.name].dec*1800000./pi),
                                                 int(mjd),int(msec))
            mjd = mjd + dt

        return tbl
    else:
        return 'Source '+srcname+' not found'

def geosat_tracktable(srcname,mjd1=None,mjd2=None):
    '''Generate a tracktable for the given geostationary satellite.
       If MDJ1 and MJD2 are not set, the tracktable covers just over
       one day (the current UT day) with entries 1/2 hour apart. 
    '''
    if mjd1 is None:
        d = util.datime()
        mjd1 = int(d.get())
        mjd2 = mjd1 + 1.

    # Retrieve TLE file for geostationary satellites from Celestrak site.
    f = urllib2.urlopen('http://www.celestrak.com/NORAD/elements/geo.txt')
    lines = f.readlines()
    f.close()

    # Names in schedule have underscore ('_') in place of spaces, so we have
    # to convert back to spaces
    name = srcname.replace('_',' ')

    # Find the source name in the GEOSAT TLE file
    for i,line in enumerate(lines):
        if line.find(name) == 0:
            break

    if name in line:
        aa = eovsa_array()
        aa.compute_pressure()
        src = ephem.readtle(lines[i], lines[i+1], lines[i+2])
        tbl = ''
        dt = 0.5/24.
        mjd = mjd1 - 2*dt
        while mjd < (mjd2+4*dt):
            jultime = mjd + 2400000.5
            aa.set_jultime(jultime)
            src.compute(aa)
            msec = round((mjd % 1)*86400.)*1000.
            tbl = tbl+'{:7d} {:7d} {:5d} {:8d}\n'.format(int(src.ra*1800000./pi),
                                                 int(src.dec*1800000./pi),
                                                 int(mjd),int(msec))
            #print src.alt, src.az
            mjd = mjd + dt

        return tbl
    else:
        return 'Source '+name+' not found'

