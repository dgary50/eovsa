# -*- coding: utf-8 -*-
"""
Created on Wed Jun 18 18:14:41 2014

@author: jackie
"""
#   2015-May-29  DG
#     Converted from using datime() to using Time() based on astropy.
#   2016-Dec-10  DG
#     Added check for restricted range of "old" 2m antennas.
#   2017-Jan-31  DG
#     27-m HA limit has to be less that -51 for now, due to cable tension
#     problem discovered yesterday, so a check for that is added to
#     check_27m_visible().
#   2017-Jul-13  DG
#     27-m cable tension was fixed, so now we can go back to the 55-degree limit.

from eovsa_lst import eovsa_ha
from numpy import rad2deg
import util

def get_27m_HAlim(dec):
    # dec should be in degrees and must be between 90 and -90
    # HA limit is returned in degrees
    # note that the HA limit is symmetric about zero so this function
    # returns only the positive limit
    #
    # declimit1 through 4 and HAlim match the parameters on page 150 of
    # the 27-m O&M manual so that if these parameters get changed, it is
    # easy to change this function
    
    declim1 = 87.5    
    declim2 = -23.
    declim3 = -45.
    declim4 = -57.6
    HAlim = 55.
    
    if dec > declim1 or dec < declim3:
        return 0
    elif dec > declim2:
        return 55
    else:
        slope = HAlim/(declim2-declim4)
        return slope*(dec-declim4)

def check_27m_visible(ha,dec):
    # ha and dec should be in degrees
    
    # note: get_27m_HAlim returns 0 outside the min and max dec, so we
    #       only need to check ha against HAlim    
    HAlim = get_27m_HAlim(dec)
    if abs(ha)<HAlim:
        #if ha < -51:
        #    # Temporary HA limit due to problem with HA cable on 27-m, 2017-01-31 (Commented out 2017-07-13)
        #    return False
        return True
    else:
        return False

def check_2m_visible(az,alt):
    # az and alt should be in degrees
    # limits: alt between 10 and 88 degrees
    #         az between 30 and 330 degrees
    
    if az > 30 and az < 330 and alt > 10 and alt < 88:
        return True
    else:
        return False  

def check_2meq_visible(ha,dec):
    # ha and dec should be in degrees
    # limits: ha  between -58 and 58 degrees
    #         dec between -24 and 45 degrees
    
    if ha >= -58 and ha <= 58 and dec >= -24 and dec <= 45:
        return True
    else:
        return False  

def check_visible(src,aa,t=None,check27m=True,check2m=True):
    ''' Check whether a source is observable from the location of aa
        given the limits on az/el of the 2-m's and HA/dec of the 27-m's, at
        date/time dt. Returns True if visible, False otherwise.
        
        If no time is given, it checks for the current time.
        
        By default, returns True only if it is visible from both the 2-m's and the
        27-m's; use check27m=False or check2m=False to ignore either the 27-m's or
        the 2-m's, respectively.
    '''
    
    if t is None:
        t = util.Time.now()
    
    # get alt, az, ha, dec at given date/time    
    aa.set_jultime(t.jd)
    src.compute(aa)
    alt = rad2deg(src.alt)
    az  = rad2deg(src.az)
    ha  = rad2deg(eovsa_ha(src,t))
    dec = rad2deg(src.dec)
    
    if check27m and check2m:
        return check_27m_visible(ha,dec) and check_2m_visible(az,alt) and check_2meq_visible(ha,dec)
    elif check27m:
        return check_27m_visible(ha,dec)
    elif check2m:
        return check_2m_visible(az,alt) and check_2meq_visible(ha,dec)
    else:
        return True

def scan_visible(src,aa,trange,check27m=True,check2m=True):
    ''' Check whether a source is observable from the location of aa
        given the limits on az/el of the 2-m's and HA/dec of the 27-m's, at both
        times in trange. (This tells you if the source is visible over
        the duration of a scan, as long as the scan is short enough.)
        Returns True if visible at both times, False otherwise.
        
        By default, returns True only if it is visible from both the 2-m's and the
        27-m's; use check27m=False or check2m=False to ignore either the 27-m's or
        the 2-m's, respectively.
    '''
    
    visible_start = check_visible(src,aa,trange[0],check27m,check2m)
    visible_stop  = check_visible(src,aa,trange[1],check27m,check2m)
    return visible_start and visible_stop