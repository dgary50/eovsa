#  2015-May-29  DG
#    Converted from using datime() to using Time() based on astropy.

##from datime import *
##from Astrometry import *
##from Angle import *
##from pyslalib import *
##
# This is the first version of the old code
##def eovsa_lst(datin=None):
##
##    if datin is None:
##        datin = datime()
##    # OVRO Longitude
##    longitude = Dec_Angle()
##    longitude.set(-118.286952892965,'degrees')
##    secondsperday_ = 86400.
##
##    ast = Astrometry()
##    # Get values of ut1-utc for today and two following days
##    mjd = int(datin.get())
##    mjd_ = mjd
##
##    for mjd_ in [mjd,mjd+1,mjd+2]:
##       ast.eqneqx.extend(mjd_,slalib.sla_eqeqx(mjd_)) # Radians
##
##    return ast.mjdutc2lst(datin.get(),longitude, ut1utc=-0.42)

#from math import *
#from util import *
#import ephem

# This is the second version of the old code
#def eovsa_lst(datin=None):
#    # Uses ephem to find LST.  Note that ephem does not directly
#    # provide LST, so it is calculated from Alt and Az of Sun for
#    # the OVSA site at the current time (or time given as argument,
#    # with datin a datime() object).  LST = RA + HA
#    #
#    sun = ephem.Sun()
#    hasun = eovsa_ha(sun,datin)
#    lst = RA_Angle()
#    lst.set(sun.ra + hasun)
#    return lst

#def eovsa_ha(src,datin=None):
#    # Returns the hour angle of the provided ephem src an observer at OVRO, 
#    # for the time in the datime() object datin, or if not given, for the current
#    # moment.
#    if datin is None:
#        datin = datime()
#    ovsa = ephem.Observer()
#    ovsa.date = datin.get() - 15019.5  # Converts MJD to ephem date
#    ovsa.lon = '-118.286953'
#    ovsa.lat = '37.233170'
#    ovsa.elevation = 1200
#    ovsa.pressure = 0  # Eliminates refraction due to atmosphere
#    src.compute(ovsa)
#    # Calculate HA of source based on Alt and Az
#    ha = atan2(-sin(src.az)*cos(src.alt),-cos(src.az)*sin(ovsa.lat)*cos(src.alt)
#                                     + sin(src.alt)*cos(ovsa.lat))
#    return ha

from eovsa_array import *
from math import pi
from util import Time

# New code is ridiculously simple
def eovsa_lst(tin=None):
    ''' Input is a Time() object (or None to use current time).
        Returns local sidereal time for EOVSA. NB: Now returns LST in radians,
        not as an RA_Angle()
    '''
    if tin is None:
        tin = Time.now()
    aa = eovsa_array()
    aa.set_jultime(tin.jd)
    return aa.sidereal_time()

def eovsa_ha(src,tin=None):
    ''' Input is a Time() object (or None to use current time).
        Returns the hour angle of the provided src for an observer at OVRO, 
        for the time in the Time() object tin, or if not given, for the current
        moment.
    '''
    if tin is None:
        tin = Time.now()
    ha = eovsa_lst(tin) - src.ra
    if ha > pi:
        ha = ha - 2*pi
    elif ha < -pi:
        ha = 2*pi + ha
    return ha

