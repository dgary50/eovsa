#   2015-May-29  DG
#      Converted from using datime() to using Time() based on astropy.
#   2018-Jan-16  DG
#      Added azel2radec() function.  Hmm.  In order to get the routines
#      to agree I had to wrestle with some changes. I think the earlier
#      routines worked--in any case they work now.
#   2018-Jan-19  DG
#      More changes to get dradec2dazel() to work.  I tested it in
#      detail to verify it.  One change--inputs and outputs are
#      in radians.

from eovsa_lst import *
from numpy import pi, sin, cos, arcsin, arccos, arctan2
lat = 37.233170*pi/180
lng = -118.286953*pi/180
def radec2azel(ra, dec, t):
    ''' Adapted from idlastro routine hadec2altaz '''
    # Provide RA, Dec in radians, and time as Time() object t
    #lat = 37.233170*pi/180
    #lng = -118.286953*pi/180
    ha = eovsa_lst(t) - ra
    
    sinel = sin(dec)*sin(lat) + cos(dec)*cos(lat)*cos(ha)
    el = arcsin(sinel)

    az = arctan2(-cos(dec)*sin(ha), sin(dec)*cos(lat) - cos(dec)*cos(ha)*sin(lat))

    if az < 0:
        az += 2*pi
    
    return az, el#, a
    
def azel2radec(az, el, t):
    ''' Adapted from idlastro routine altaz2hadec '''
    # Provide Az, El in radians, and time as Time() object t
    #lat = 37.233170*pi/180
    #lng = -118.286953*pi/180

    dec = arcsin(sin(el)*sin(lat) + cos(el)*cos(lat)*cos(az))
    ha = arctan2(-sin(az)*cos(el), sin(el)*cos(lat) - cos(el)*cos(az)*sin(lat))
    ra = eovsa_lst(t) - ha
    if ra < 0:
        ra += 2*p1
    return ra, dec    
    
def dradec2dazel(ra,dec,t,dra,ddec):
    # Provide RA, Dec in radians, and time as Time() object t
    # Provide dra and ddec in radians (NB: Does NOT work if in degrees)
    # Returns dXel = daz/cos(el) and d_el in radians
    ha = eovsa_lst(t) - ra

    az1, el1 = radec2azel(ra, dec, t)
    az2, el2 = radec2azel(ra+dra, dec+ddec, t)
    
    daz = az2 - az1
    d_el = el2 - el1
    
    return daz*cos(el1), d_el
    
def old_radec2azel(ra, dec, t):
    # Provide RA, Dec in radians, and time as Time() object t
    lat = 37.233170*pi/180
    lng = -118.286953*pi/180
    ha = eovsa_lst(t) - ra
    
    sinel = sin(dec)*sin(lat) + cos(dec)*cos(lat)*cos(ha)
    el = arcsin(sinel)

    cosa = (sin(dec) - sinel*sin(lat)) / (cos(el)*cos(lat))
    a = arccos(cosa)

    if ha <= 0:
        az = a
    else:
        az = 2*pi - a

    return az, el, a
    
def old_dradec2dazel(ra,dec,t,dra,ddec):
    # Provide RA, Dec in radians, and time as Time() object t
    # Provide dra and ddec in radians
    # Returns dXel = daz/cos(el) and d_el
    ha = eovsa_lst(t) - ra

    az, el, a = old_radec2azel(ra, dec, t)
    #az2, el2 = radec2azel(ra+dra, dec+ddec, t)
    
    #da = az2 - az1
    #d_el = el2 - el1
    
    d_el = (ddec*(cos(dec)*sin(lat) - sin(dec)*cos(lat)*cos(ha)) + dra*cos(dec)*cos(lat)*sin(ha)) / cos(el)

    da = (d_el*cos(el)*sin(lat) - ddec*cos(dec)) / (cos(el)*cos(lat)*sin(a)) \
         + (d_el*(sin(el)*sin(lat)-sin(dec))*sin(el)*cos(lat)) / ((cos(el)*cos(lat))**2*sin(a))

    if ha <= 0:
        daz = da
    else:
        daz = -da

    return daz*cos(el), d_el