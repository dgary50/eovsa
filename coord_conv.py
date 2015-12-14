#   2015-May-29  DG
#      Converted from using datime() to using Time() based on astropy.

from eovsa_lst import *
from numpy import pi, sin, cos, arcsin, arccos
def radec2azel(ra, dec, t):
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

def dradec2dazel(ra,dec,t,dra,ddec):
    # Provide RA, Dec in radians, and time as Time() object t
    # Provide dra and ddec in radians
    # Returns dXel = daz/cos(el) and d_el
    lat = 37.233170*pi/180
    lng = -118.286953*pi/180
    ha = eovsa_lst(t) - ra

    az, el, a = radec2azel(ra, dec, t)
    
    d_el = (ddec*(cos(dec)*sin(lat) - sin(dec)*cos(lat)*cos(ha)) + dra*cos(dec)*cos(lat)*sin(ha)) / cos(el)

    da = (d_el*cos(el)*sin(lat) - ddec*cos(dec)) / (cos(el)*cos(lat)*sin(a)) \
         + (d_el*(sin(el)*sin(lat)-sin(dec))*sin(el)*cos(lat)) / ((cos(el)*cos(lat))**2*sin(a))

    if ha <= 0:
        daz = da
    else:
        daz = -da

    return daz*cos(el), d_el
    
