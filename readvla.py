#   Routines to read VLA calibrator file
#
# 2015-May-29  DG
#    Converted from using datime() to using Time() based on astropy
#

class Calibrator():
    def setra(self, rastr):
        rastr = rastr.replace('h', ':')
        rastr = rastr.replace('m', ':')
        rastr = rastr.replace('s', ':')
        self.ra.append(rastr)
    def setdec(self, decstr):
        decstr = decstr.replace('d', ':')
        decstr = decstr.replace('\'', ':')
        decstr = decstr.replace('"', ':')
        self.dec.append(decstr)
        

def readvlacaldb():
    #filename = 'C:\Home\OVSA Expansion\Design\Calibration\\vla_calibrator_list.txt'
    filename = 'vla_calibrator_list.txt'  # Must be in current directory
    f = open(filename)
    # read lines and close file
    lines = f.readlines()
    f.close()
    cal = []
    k = -1
    nextline = False
    for line in lines[3:]:
        # print line[11:16]
        if line[11:16] == 'J2000':
            cal.append(Calibrator())
            k += 1
            cal[k].name = line[0:8]
            cal[k].ra = []
            cal[k].dec = []
            cal[k].setra(line[20:32])
            cal[k].setdec(line[37:50])
            cal[k].band = []
            cal[k].qual = []
            cal[k].flux = []
            cal[k].uvmin = []
            cal[k].uvmax = []
            if line[63:] != '':
                cal[k].altname = line[63:]
        elif line[0:1] is '=':
            nextline = True
        elif nextline:
            if len(line) == 2:
                # Empty line, so set up for a new calibrator
                nextline = False
            else:
                # A good line, should have what we need
                cal[k].band.append(line[9:10])
                cal[k].qual.append(line[16:17])
                try:
                    cal[k].flux.append(float(line[24:31]))
                except ValueError:
                    cal[k].flux.append(0.0)
                if line[41:43] != '' and line[41:43] != '  ':
                    try:
                        cal[k].uvmin.append(float(line[41:43]))
                    except ValueError:
                        pass
                if line[52:55] != '' and line[52:55] != '   ' and line[52:55] != ' vi':
                    try:
                        cal[k].uvmax.append(float(line[52:55]))
                    except ValueError:
                        pass
    return cal

from numpy import array
from math import pi
from util import *
from eovsa_array import *
import ephem

def compute_cal(cal,t):
    '''Given a cal object representing one of the VLA calibrators, 
       convert to an ephem object observed at OVSA for the time in 
       the given Time() object, t.  Returns the computed source
       and the antenna array object.
    '''
    time = t.iso.replace('-','/')[:16]
    aa = eovsa_array()
    aa.set_ephemtime(time)
    aa.compute_pressure()
    flux = '1.0'   # Just a dummy flux since no band is specified
    name = cal.name
    ra = cal.ra[0]
    dec = cal.dec[0]
    src = ephem.readdb(name+',f,'+ra+','+dec+','+flux+',2000')
    src.compute(aa)
    return src, aa

def findcal(cal, band='C', qual='P', t=None, dtheta=[20,40]):
    '''Searches the VLA calibrator database [cal list returned by
       a call to readvlacaldb()] to find sources matching the band
       and quality codes provided:
           band = 'L', 'C', 'X', 'U', 'Q' (default 'C')
           qual = 'P', 'S', 'X', 'W' (default 'P')
       If t and dtheta are provided, only sources within an annulus 
       from dtheta[0] to dtheta[1] degrees from the Sun are selected, 
       for the time specified in Time() object t. Default is 20-40 deg.
       Returns a list of ephem source objects.  If t is not provided, 
       the current time is used.  In all cases, the first source in the 
       list is the Sun.
    '''
    srclist = []
    dtor = pi/180.
    if t is None:
        time = Time.now().iso.replace('-','/')[:16]
    else:
        time = t.iso.replace('-','/')[:16]
    aa = eovsa_array()
    aa.set_ephemtime(time)
    aa.compute_pressure()
    sun = ephem.Sun()
    sun.compute(aa)
    ## Add sun to source list in first location
    srclist.append(sun)
        
    for i in range(len(cal)):
        if band in cal[i].band:
            j = (array(cal[i].band) == band) & (array(cal[i].qual) == qual)
            if j.any():
                flux = str(array(cal[i].flux)[j][0])
                name = cal[i].name
                ra = cal[i].ra[0]
                dec = cal[i].dec[0]
                src = ephem.readdb(name+',f,'+ra+','+dec+','+flux+',2000')
                src.compute(aa)
                if t is None:
                    # No selection based on distance to Sun, so append
                    srclist.append(src)
                else:
                    sep = ephem.separation(sun,src)
                    if sep > dtheta[0]*dtor and sep < dtheta[1]*dtor:
                        # This source qualifies so append it to list
                        srclist.append(src)

    return srclist, aa

