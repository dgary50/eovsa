# Last edited by: Jackie Villadsen, 4/25/14
# Original version: Jackie Villadsen, 4/24/14

import aipy, ephem, numpy
import numpy
import util
from math import cos, sin
from numpy import pi
from eovsa_array import *

class Source:
    pass

def read_starburst_srclist():
    ''' Read in the file starburst_srclist.txt and return it as a list of
        objects with attributes name, ra, and dec
    '''
    f = open('starburst/starburst_srclist.txt')
    temp_linelist = f.readlines() # identifies ends of lines using \n
    f.close()
    
    # For some reason, when I make the source list by copying from Excel, lines
    # are separated by \r instead of \n, so readlines doesn't split them apart.
    # In order to make this flexible (in case later the source list has \n in it),
    # I am setting this up to split lines by both \r and \n.
    linelist = []
    for t in temp_linelist:
        lines = t.split('\r')
        for l in lines:
            linelist.append(l)
    
    srclist = []
    for l in linelist:
        s = Source()
        # strip() gets rid of the \n on the end
        # split breaks apart fields separated by tabs
        s.name,s.ra,s.dec = l.strip().split()
        srclist.append(s)
    
    return srclist


def import_source(srcname):
    # takes the source name (such as 'V*ADLEO') and returns an aipy RadioFixedBody object
    # with the attributes name, ra, and dec set (not magnitude)
    
    srclist = read_starburst_srclist()
    for s in srclist:
        if s.name.find(srcname) == 0:
            break
    
    if srcname in s.name:
        print s.ra, s.dec
        src = aipy.phs.RadioFixedBody(s.ra,s.dec,name=s.name,epoch='2000')
        print src
        return src
    else:
        return 'Source '+srcname+' not found in Starburst catalog. Not added to schedule!'


def star_tracktable(srcname,mjd1=None,mjd2=None):
    '''Generate a tracktable for the given Starburst source (flare star).
       This is just three lines giving the RA, Dec coordinates and times.
       Important: aipy is used to precess the coordinates to the current epoch
       from J2000 - this can be an effect of a couple tenths of a degree.
    '''
    if mjd1 is None:
        d = util.datime()
        mjd1 = int(d.get())
        mjd2 = mjd1 + 1.

    src = import_source(srcname)

    if type(src) is str:
        return 'Source '+srcname+' not found'
    else:
        aa = eovsa_array()
        aa.compute_pressure()
        cat = aipy.phs.SrcCatalog(src)
        tbl = ''
        dt = 0.5/24.
        mjd = mjd1 - 2*dt
        while mjd < (mjd2+4*dt):
            jultime = mjd + 2400000.5
            aa.set_jultime(jultime)
            cat.compute(aa)
            msec = round((mjd % 1)*86400.)*1000.
            tbl = tbl+'{0:7d} {1:7d} {2:5d} {3:8d}\n'.format(int(cat[srcname].ra*1800000./pi),
                                                 int(cat[srcname].dec*1800000./pi),
                                                 int(mjd),int(msec))
            mjd = mjd + dt
        return tbl
