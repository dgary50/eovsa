#
# History:
#  2014-Dec-07  DG
#    The schedule was dying today due to one of the GEOSAT satellites
#    having 0.000 for its inclination.  This was tricky to find!  The
#    code has been changed to adjust such lines to 0.0001 (and update
#    the checksum in that line).  This is a very rare occurrence...
#
#  2015-Apr-21  JV
#    The schedule was crashing due to an error in readtle b/c of a formatting error
#    in the geosat file online at http://www.celestrak.com/NORAD/elements/geo.txt.
#    Added error checking to load_geosats() function so that it prints an error w/o
#    crashing the schedule.
#  2015-May-30  DG
#    Removed unneeded import of util's Vector and datime routines.
#  2015-Oct-24  DG
#    Added Venus
#  2016-Dec-14 BC
#    Changed the method for locating the calibrator source files in load_sidereal_cats()
#  2019-Feb-26  DG
#    Added GPS satellites, now that we have L band!
#  2019-Apr-06  DG
#    Apparently I had the wrong GPS elements text file.  It should
#    have been https://celestrak.com/NORAD/elements/gps-ops.txt.  Now fixed.
#  2024-Apr-25  DG
#    Apparently celestrak.com is now celestrak.org.  The calling sequence for the
#    files also changed, apparently (although the direct URL to text files also
#    works in some cases).  Here are the URLs:
#      GEO: https://celestrak.org/NORAD/elements/gp.php?GROUP=geo&FORMAT=tle
#      GPS: https://celestrak.org/NORAD/elements/gp.php?GROUP=gps-ops&FORMAT=tle
#      O3B: https://celestrak.org/NORAD/elements/gp.php?GROUP=other-comm&FORMAT=tle

import aipy, ephem, numpy
from math import cos, sin
from numpy import pi, mat
from readvla import readvlacaldb
from eovsa_array import *
import urllib2
import re
import os

global lat
lat = 37.233170*numpy.pi/180       # OVSA Latitude (radians)

class RadioGeosat(aipy.phs.RadioBody, object):
    '''A geostationary satellite.  Combines ephem versions of these objects with RadioBody.
       Modeled off of aipy.phs class RadioSpecial (which is for major Solar System bodies).
    '''
    def __init__(self,geosat_body, mfreq=.150,
            ionref=(0.,0.), srcshape=(0.,0.,0.), **kwargs):
        """`name' is used to lookup appropriate ephem celestial object."""
        aipy.phs.RadioBody.__init__(self, geosat_body.name, mfreq, ionref, srcshape)
        self.Body = geosat_body
    def __getattr__(self, nm):
        """First try to access attribute from this class, but if that fails, 
        try to get it from the underlying ephem object."""
        try: return object.__getattr__(self, nm)
        except(AttributeError): return self.Body.__getattribute__(nm)
    def __setattr__(self, nm, val):
        """First try to set attribute for this class, buf if that fails, 
        try to set it for the underlying ephem object."""
        try: object.__setattr__(self, nm, val)
        except(AttributeError): return setattr(self.Body, nm, val)
    def compute(self, observer):
        self.Body.compute(observer)
        aipy.phs.RadioBody.compute(self, observer)

def load_geosats():
    ''' Read the list of geostationary satellites from the Celestrak site and create a list
        of RadioGeosat objects containing all satellites. (List contains 399 sats as of 6/19/14.)
    '''
    # Retrieve TLE file for geostationary satellites from Celestrak site.
    try:
        f = urllib2.urlopen('https://celestrak.org/NORAD/elements/gp.php?GROUP=geo&FORMAT=tle')
    except urllib2.URLError as err:
        print 'Error reading GEO satellite web file:', err
        return []
        
    lines = f.readlines()
    f.close()
    nlines = len(lines)
    
    # use every 3 lines to create another RadioGeosat object
    satlist = []
    for i in range(0,nlines,3):
        if lines[i+2][9:16] == ' 0.0000':
            # aa.compute() hangs for a satellite with zero inclination!
            # Change to 0.0001 degrees, and do not forget to change the checksum.
            chksum = str(int(lines[i+2][-4:]) + 1)
            lines[i+2] = lines[i+2][:9]+' 0.0001'+lines[i+2][16:-4]+chksum
        try:
            geosat_body = ephem.readtle(lines[i], lines[i+1], lines[i+2])
            src = RadioGeosat(geosat_body) # convert from an ephem Body object to a RadioGeosat object
            satlist.append(src)
        except:
            print 'Error in ephem.readtle: Geosat', lines[i].strip(), 'not added to source catalog.'
    return satlist
    
def load_gpssats():
    ''' Read the list of global positioning satellites from the Celestrak site and create a list
        of RadioGeosat objects containing all satellites. (List contains 31 sats as of 2/26/2019.)
    '''
    # Retrieve TLE file for geostationary satellites from Celestrak site.
    try:
        f = urllib2.urlopen('https://celestrak.org/NORAD/elements/gp.php?GROUP=gps-ops&FORMAT=tle')
    except urllib2.URLError as err:
        print 'Error reading GPS satellite web file:', err
        return []
    lines = f.readlines()
    f.close()
    nlines = len(lines)
    
    # use every 3 lines to create another RadioGeosat object
    satlist = []
    for i in range(0,nlines,3):
        if lines[i+2][9:16] == ' 0.0000':
            # aa.compute() hangs for a satellite with zero inclination!
            # Change to 0.0001 degrees, and do not forget to change the checksum.
            chksum = str(int(lines[i+2][-4:]) + 1)
            lines[i+2] = lines[i+2][:9]+' 0.0001'+lines[i+2][16:-4]+chksum
        try:
            geosat_body = ephem.readtle(lines[i], lines[i+1], lines[i+2])
            src = RadioGeosat(geosat_body) # convert from an ephem Body object to a RadioGeosat object
            satlist.append(src)
        except:
            print 'Error in ephem.readtle: Geosat', lines[i].strip(), 'not added to source catalog.'
    return satlist

def load_o3bsats():
    ''' Read the list of ob3 satellites from the Celestrak site and create a list
        of RadioGeosat objects containing all satellites.  
    '''
    # Retrieve TLE file for o3b satellites from Celestrak site.
    try:
        f = urllib2.urlopen('https://celestrak.org/NORAD/elements/gp.php?GROUP=other-comm&FORMAT=tle')
    except urllib2.URLError as err:
        print 'Error reading ob3 satellite web file:', err
        return []

    lines = f.readlines()
    f.close()
    nlines = len(lines)
    
    # use every 3 lines to create another RadioGeosat object
    satlist = []
    for i in range(0,nlines,3):
        if lines[i+2][9:16] == ' 0.0000':
            # aa.compute() hangs for a satellite with zero inclination!
            # Change to 0.0001 degrees, and do not forget to change the checksum.
            chksum = str(int(lines[i+2][-4:]) + 1)
            lines[i+2] = lines[i+2][:9]+' 0.0001'+lines[i+2][16:-4]+chksum
        try:
            geosat_body = ephem.readtle(lines[i], lines[i+1], lines[i+2])
            src = RadioGeosat(geosat_body) # convert from an ephem Body object to a RadioGeosat object
            satlist.append(src)
        except:
            print 'Error in ephem.readtle: o3bsat', lines[i].strip(), 'not added to source catalog.'
    return satlist

def load_VLAcals():
    ''' Read the list of VLA calibrators and create a list of RadioFixedBody objects containing
        all calibrators.
        
        The list contains 1865 calibrators as of 6/19/14. However, there are four with duplicate
        names so when I add them to the aipy SrcCatalog it ends up containing only 1861 VLA
        calibrators since you can't have two dictionary entries with the same key (the first one
        just gets overwritten by the second).  The sources with duplicate names are:
            '1914+166', '0354+801', '0632+159', '1300+142'
        The code to find these duplicate names is:
            src_names = [s.src_name for s in srclist]
            from collections import Counter
            [s for s,n in Counter(src_names).items() if n>1]
        I haven't looked at the calibrator list to figure out why there are duplicates - are these
        distinct sources? If someone has an issue with this they will have to figure it out.
    '''
    cal_list = readvlacaldb()
    srclist = []
    for c in cal_list:
        src = aipy.phs.RadioFixedBody(c.ra[0],c.dec[0],name=c.name,epoch='2000')
        srclist.append(src)
    return srclist

def load_sidereal_cats():
    ''' Read all files in directory Dropbox/PythonCode/Current/SourceCat with extension .srclist
        and for each line make a RadioFixedBody object and return a list of these objects.
        
        The SourceCat files should have name, RA (h:m:s), dec (d:m:s) in the first 3 columns,
        and can have an optional 4th column with flux in Jy.
    '''
    # both \r and \n can be used in files to mark new lines and sometimes both
    # so split lines based on either and don't split twice if there is more than one in a row
    if not os.getenv('EOVSAPY'):
        cmd='cat SourceCat/*.srclist' # from the current directory
    else:
        cmd='cat '+os.path.expandvars('$EOVSAPY')+'/SourceCat/*.srclist'
    tmp = os.popen(cmd).read().strip().replace('\r','\n')
    lines = re.split('\n+',tmp)
    
    srclist = []
    for l in lines:
        props = l.strip().split()
        if len(props) == 3:
            srcname, ra, dec = props
            src = aipy.phs.RadioFixedBody(ra,dec,name=srcname,epoch='2000')
        elif len(props) == 4:
            srcname, ra, dec, fluxJy = props
            src = aipy.amp.RadioFixedBody(ra,dec,name=srcname,jys=fluxJy,epoch='J2000')
        srclist.append(src)
    return srclist

def load_cat():
    ''' Create standard cat with Sun, Moon, all VLA calibrators (N~2000),
        all geosats from Celestrak (N~400),
        and all sidereal sources whose coords are listed in a .srclist file in the
        src_cat directory (so you can add a .srclist file there and it will automatically
        be added to the catalog created when this function is called).
        
        Note: this function does not run compute yet - the catalog returned is generic to all
        Observer locations.
    '''
    srclist = load_VLAcals() + load_geosats() + load_sidereal_cats() + load_o3bsats() + load_gpssats()
    
    # append Sun and Moon
    # use aipy.amp RadioSpecial objects and SrcCatalog object - they are extensions
    # of the aipy.phs classes by the same names, with the addition that they allow
    # you to set and retrieve source fluxes (in Jy)
    srclist.append(aipy.amp.RadioSpecial('Sun'))
    srclist.append(aipy.amp.RadioSpecial('Moon'))
    srclist.append(aipy.amp.RadioSpecial('Venus'))
    cat = aipy.amp.SrcCatalog(srclist)
    return cat

def eovsa_array_with_cat():
    ''' Return an aa object created by the eovsa_array module but with a source
        catalog in the .cat attribute.
    '''
    aa = eovsa_array()
    cat = load_cat()
    cat.compute(aa)
    aa.cat = cat
    return aa
