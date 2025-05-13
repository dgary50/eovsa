#!/usr/bin/env python

#
# Routines for reading, interpreting, and displaying the stateframe
#
# Updated to use read_xml2
#
# History:
#   2014-Nov-29  DG
#      Added a mew routine azel_from_sqldict(), which does the same thing as 
#      azel_from_stateframe().
#   2014-Dec-06  DG
#      Fixed to work on dpp as well as helios...This has to be generalized.
#      Also now gets datime() from util.datime()
#   2014-Dec-18  DG
#      Rather than look for specific host names, now call socket.getfqdn(),
#      and look for 'solar.pvt'.  This allows it to work on any machine in
#      the local domain.
#   2015-Jan-30  DG
#      Add a second attempt to read weather station on failure, hopefully
#      eliminating so many "problem parsing" messages.
#   2015-Mar-28  DG
#      Argh.  Magnum Energy changed the format of the solar power station
#      web page, so I had to make changes to rd_solpwr() to read it
#      correctly.
#   2015-May-29  DG
#      Converted from using datime() to using Time() based on astropy
#   2015-Jun-16  DG
#      FTP to ACC now requires a username and password
#   2015-Jun-25  DG
#      Now that Ant13's solar power station is online, changed rd_solpwr() to 
#      read from either power station depending on supplied url.  Now returns
#      a single dictionary. 
#   2016-Jan-15  DG
#      Cleaned up azel_from_stateframe() code to use extract().
#   2016-Aug-12  DG
#      Doh!  Subtle bug in parallactic angle calculation, doing lat**dtor 
#      instead of lat*dtor!  All of the stateframe values up to now are
#      wrong!  Now fixed...
#   2016-Oct-20  DG
#      Added hadec2altaz() routine, so that I can find parallactic angle
#      more easily vs. HA, DEC.
#   2016-Nov-15  DG
#      Added PA_adjust() routine to continuously adjust the PA of the 
#      27-m feed to track the parallactic angle of a given antenna.
#   2016-Nov-16  DG
#      Corrected long-standing problem with parallactic angle--apparently
#      my hadec2altaz() routine was never completed, and never used to
#      convert HA, Dec to Az, El.  Now it should work as intended.
#   2016-Nov-20  DG
#      Had to update hadec2altaz(), because my changes made it no longer
#      work for array arguments.
#   2016-Dec-10  DG
#      Added PA_sweep() routine to sweep PA of 27-m feed from -PA to PA,
#      at a specified rate.
#   2016-Dec-12  DG
#      Tweaked PA_sweep() to wait until initial PA is acquired before
#      starting the sweep (or times out and starts sweep after 2 minutes).
#   2017-Feb-01  DG
#      Increase timeout for solar power stations from 0.2 to 0.4 s.
#   2017-Jun-28  DG
#      Add "crossed" keyword to PA_adjust() to orient the feed 90 degrees
#      from the parallactic angle
#   2017-Aug-09  DG
#      Fixed a bug in TrackFlag and dAz, in azel_from_stateframe() and 
#      azel_from_sqldict()
#   2018-Jan-10  DG
#      Added control_room_temp() function to return the ambient temperature in
#      the control room.
#   2018-Mar-01  DG
#      Added temporary get_median_wind() routine to get around current
#      weather station glitches.  Replaces average wind with median wind speed.
#      Must change weather() back to original when this has been fixed.
#   2019-Feb-12  DG
#      Discovered a bug in reading the solar power temperatures, when FET temp
#      was less than 10 C.  Fixed by reading either two or one digit.
#   2019-Apr-27  DG
#      Timeout of control_room_temp() was throwing an error and crashing.
#      This is now fixed by putting readlines() call inside the try: except:
#      clause.
#   2021-Feb-03  DG
#      Added a tracksrcflag to indicate when the antennas are supposed to be
#      tracking the source (i.e. no intentional offsets).
#   2022-Mar-07  DG
#      Oops--"temporary" change in 2018 (4 years ago!) was never reversed.  
#      I have taken it out now, since SQL is not working...
#   2024-Dec-09  DG
#      Fix bug in PA_sweep() to avoid a crash when the stateframe is not
#      successfully read.
#   2025-Apr-23  DG
#      Changes to reflect removal of equatorially mounted antennas and replacement
#      with AzEl ones for ants 9-13.  This was just a simple change in 
#      azel_from_stateframe().

import struct, sys
import socket
import urllib2
import numpy as np
from read_xml2 import xml_ptrs
import copy
from util import Time
from Tkinter import Tk
from tkFileDialog import *
import xml.etree.ElementTree as ET
import Queue
q = Queue.Queue()

#============================
def control_room_temp():
    '''Read the 'http://192.168.24.233/state.xml' page and return the
       ambient temperature, in C.  If reading data fails, returns an
       impossible number, -99 C.
    '''
    try:
        f = urllib2.urlopen('http://192.168.24.233/state.xml',timeout=0.4)
        lines = f.readlines()
        f.close()
    except:
        # Timeout error
        print Time.now().iso,'Control room temperature connection timed out'
        return -99.0
    try:
        return int((float(lines[3][13:17]) - 32)*50/9.)/10.
    except:
        return -99.0
        
#============================
def weather(attempt=0):
    '''Read the http://wx.cm.pvt/latestsampledata.xml page and
    take the title and the information and put it in dictionary form'''

    try:
        f = urllib2.urlopen('http://wx.cm.pvt/latestsampledata.xml',timeout=0.4)
    except:
        # Timeout error
        print Time.now().iso,'Weather connection timed out'
        return {}
    try:
        #tree = ET.parse(f)
        line = f.readline()
        if line.find('</oriondata>') == -1:
           # Line is often truncated, so fix it if possible
           line = line[:line.find('</o')]+'</oriondata>'
           print Time.now().iso,'Fixed Weather info'
        #tree = ET.XML(line)
    except:
        # Error reading weather info, so return blank dictionary
        print Time.now().iso,'Problem reading Weather info'
        return {}
    f.close()

    #root = tree.getroot()
    try:
        root = ET.XML(line)
    except:
        if attempt == 0:
            # Try again, then bail if it doesn't work
            return weather(attempt=1)
        # Error reading weather info, so return blank dictionary
        print Time.now().iso,'Problem parsing Weather info'
        return {}
    index = 0

    ovro_dict = {}
    for element in root.findall('meas'):
        name = element.get('name')
        text = root[index].text
        ovro_dict.update({name : text})
        index = index + 1
    
    # Convert pressure in inches Hg to mBar
    try:
        temp = ovro_dict['mtRawBaromPress']
        temp = float(temp) * 33.8637526
    except:
        return ovro_dict
    ovro_dict['mtRawBaromPress'] = str(temp)
#    return get_median_wind(ovro_dict)   # Removed due to loss of SQL
    return ovro_dict
    
def get_median_wind(wthr):
    '''  Temporary work-around for mis-behaving weather station.
         Given the weather dictionary, query the SQL database for
         the last 120-s of wind data, and calculate median rather
         than average.  I hope this does not take too long!  Returns
         the same dictionary, with median replacing average wind.
    '''
    import dbutil as db
    cursor = db.get_cursor()
    ver = db.find_table_version(cursor,  int(Time.now().lv))
    
    query = 'select top 120 Timestamp,Sche_Data_Weat_Wind from fV'+ver+'_vD1 order by Timestamp desc'
    data, msg = db.do_query(cursor,query)
    if msg == 'Success':
        try:
            medwind = np.median(data['Sche_Data_Weat_Wind'])
            wthr.update({'mt2MinRollAvgWindSpeed': medwind})
        except:
            pass
    cursor.close()
    return wthr

#============================
def rd_solpwr(url='http://data.magnumenergy.com/MW5127'):
    '''Reads the data from the solar power station at Ant 12 or 13, which is sent
       to the Magnum Energy web site and they then serve it to us at the
       address: Ant12: http://data.magnumenergy.com/MW5127.
                Ant13: http://data.magnumenergy.com/MW5241.
       Now returns a single dictionary, for whichever station is pointed to by url
    '''
    # Read and decode the information from the power station at 12
    try:
        f = urllib2.urlopen(url,timeout=0.4)
    except:
        # Timeout error
        print Time.now().iso,'Solar Power connection timed out'
        solpwr = {}
        return solpwr
    try:
        lines = f.readlines()
    except:
        print Time.now().iso,'Solar Power readlines timed out'
        lines = None
    f.close()
    solpwr = {}
    if lines is None:
        return solpwr

    for i,line in enumerate(lines):
        if i > 185:
            break
        if line.find('Data Date:<') > 0:
            t = Time(lines[i+2][23:23+19])
            solpwr.update({'Time':t.lv})
        elif line.find('State of Charge:<') > 0:
            idx = lines[i+2].find('%')
            solpwr.update({'Charge':int(lines[i+2][:idx])})
        elif line.find('volts / amps') > 0:
            args = lines[i+2].split(' ')
            solpwr.update({'Volts':float(args[0])})
            solpwr.update({'Amps':float(args[3])})
        elif line.find('Amp Hours') > 0:
            solpwr.update({'AmpHours':int(lines[i+2].split(' ')[0])})
        elif line.find('Battery Temperature') > 0:
            btemp = lines[i][lines[i].find('<td>')+4:lines[i].find('&deg;C')]
            try:
                solpwr.update({'BatteryTemp':int(btemp)})
            except:
                # When value is below 0 C, web page returns < 0 C, which cannot be set to int,
                # so just set to -1 C on error.
                solpwr.update({'BatteryTemp':-1})
        elif line.find('Transformer Temperature') >0:
            try:
                solpwr.update({'TransformerTemp':int(lines[i].split('&deg;C')[0][-2:])})
            except:
                solpwr.update({'TransformerTemp':int(lines[i].split('&deg;C')[0][-1:])})
        elif line.find('FET Temperature') >0:
            try:
                solpwr.update({'FETTemp':int(lines[i].split('&deg;C')[0][-2:])})
            except:
                solpwr.update({'FETTemp':int(lines[i].split('&deg;C')[0][-1:])})

    return solpwr

#============================
def rd_ACCfile():
    '''Reads key variables from ACC.ini file on ACC (using urllib2)
    '''
    # List of strings to search for
    s0 = '[Stateframe]'
    s1 = 'bin size = '
    s2 = 'template path = '
    n0 = '[Network]'
    n1 = 'TCP.schedule.port = '
    n2 = 'TCP.stateframe.port = '
    n3 = 'TCP.schedule.stateframe.port = '
    r0 = '[ROACH]'
    r1 = 'boffile = '
    
    userpass = 'admin:observer@'
    ACCfile = None
    if socket.getfqdn().find('solar.pvt') != -1:
        try:
            ACCfile = urllib2.urlopen('ftp://'+userpass+'acc.solar.pvt/ni-rt/startup/acc.ini',timeout=0.5)
        except:
            # Timeout error
            print Time.now().iso,'FTP connection to ACC timed out'
        # Since this is the HELIOS machine, make a disk copy of ACC.ini in the
        # current (dropbox) directory.  This will be used by other instances of
        # sf_display() on other machines that do not have access to acc.solar.pvt.
        try:
            lines = ACCfile.readlines()
            o = open('acc.ini','w')
            for line in lines:
                o.write(line+'\n')
            o.close()
            ACCfile.close()
            ACCfile = urllib2.urlopen('ftp://'+userpass+'acc.solar.pvt/ni-rt/startup/acc.ini',timeout=0.5)
            # Also read XML file for stateframe from ACC, and decode template for later use
            sf, version = xml_ptrs()
        except:
            pass
    if ACCfile is None:
        # ACC not reachable?  Try reading static files.
        print 'Cannot ftp ACC.ini.  Reading static acc.ini and stateframe.xml from current directory instead.'
        ACCfile = open('acc.ini','r')
        # Also read XML file for stateframe from static file, and decode template for later use
        sf, version = xml_ptrs('stateframe.xml')

    for line in ACCfile:
        if s0 in line:    # String s0 ([Stateframe]) found
            for line in ACCfile:
                if s1 in line:
                    binsize = int(line[len(s1):])
                elif s2 in line:
                    xmlpath = line[len(s2):]
                    break
                elif line == '':
                    break
        if n0 in line:    # String n0 ([Network]) found
            for line in ACCfile:
                if n1 in line:
                    scdport = int(line[len(n1):])
                elif n2 in line:
                    sfport = int(line[len(n2):])
                    print '\nConnecting to ACC at port:',sfport
                elif n3 in line:
                    scdsfport = int(line[len(n3):])
                    break
                elif not line:
                    break
        if r0 in line:    # String r0 ([ROACH]) found
            for line in ACCfile:
                if r1 in line:
                    boffile = line[len(r1):].strip()
                elif not line:
                    break
    ACCfile.close()
    accdict = {'host':'acc.solar.pvt','binsize':binsize,'xmlpath':xmlpath,
               'scdport':scdport,'sfport':sfport,'scdsfport':scdsfport,'sf':sf,'version':version,'boffile':boffile}
    #if socket.gethostname() != 'helios':
        # The host is not OVSA, so assume port forwarding of stateframe port
        # to localhost port 6341
        #accdict['host'] = 'localhost'
    return accdict

#============================
def rd_stateframe(s,sf_num,n_expected):
    '''Does multiple reads of opened connection s until
       n_expected bytes are read.  sf_num is sent to the
       ACC to indicate which stateframe to read (1 normally)
    '''
    totlen = 0; totdata = []; data = ''
    sf_pck = struct.pack(">i",sf_num)
    s.settimeout(0.5)
    #sys.stdout.write('+')
    #sys.stdout.flush()  # Flush stdout (/tmp/schedule.log) so we can see the output.
    try:
        s.send(sf_pck)
        #sys.stdout.write('.')
        #sys.stdout.flush()  # Flush stdout (/tmp/schedule.log) so we can see the output.
        while totlen < n_expected:
            data = s.recv(n_expected)
            totdata.append(data)
            totlen = sum([len(i) for i in totdata])
        #sys.stdout.write('-')
        #sys.stdout.flush()  # Flush stdout (/tmp/schedule.log) so we can see the output.
    except socket.timeout:
        print Time.now().iso,'Socket time-out when reading stateframe from ACC'
    return ''.join(totdata)

#============================
def get_stateframe(accini):
    '''Connects to the ACC's stateframe port and reads the
       current stateframe data.  Returns both data and a message.
       If the port cannot be opened, or cannot be read, data is None,
       and an appropriate message is returned.
    '''
    try:
        s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        try:
            s.connect((accini['host'],accini['sfport']))
            data = rd_stateframe(s,1,accini['binsize'])
            if len(data) == accini['binsize']:
                return data, 'No Error'
            else:
                return data, 'Incorrect stateframe size returned from ACC'
        except:
            return None, 'Cannot read from port '+str(accini['sfport'])
        s.close()
    except:
        return None, 'Cannot open socket to port '+str(accini['sfport'])

#============================
def get_stateframefromfile(filename,f=None,recsiz=None):
    if not f:
        f = open(filename,'rb')
        data = f.read(100)
        recsiz = struct.unpack_from('<i',data,16)[0]
        f.close()
        f = open(filename,'rb')
        data = f.read(recsiz)
        return data,'No Error',f,recsiz
    else:
        try:
            data = f.read(recsiz)
            return data,'No Error',f,recsiz
        except:
            f.close()
            return None, 'End of file reached.',None,recsiz
    
#============================
def extract(data,k):
    '''Helper function that extracts a value from data, based on stateframe
       info pair k (k[0] is fmt string, k[1] is byte offset into data)
    '''
    if len(k) == 3:
       k[2].reverse()
       val = np.array(struct.unpack_from(k[0],data,k[1]))
       val.shape = k[2]
       k[2].reverse()
    else:
       val = struct.unpack_from(k[0],data,k[1])[0]
    return val

#============================
def azel_from_stateframe(sf, data, antlist=None):
    '''Given a stateframe dictionary and a data record, calculate
       the actual and requested azimuth and elevation for each antenna, as well as the 
       difference between them, and a track flag, all as a dictionary of numpy float arrays.
    '''
    daz = []
    delv = []
    az_act = []
    el_act = []
    az_req = []
    el_req = []
    chi = []
    tracksrcflag = []
    dtor = np.pi/180.
    if antlist is None:
        # No antlist, so assume all antennas
        antlist = range(15)

    for i, ant in enumerate(antlist):
        c = sf['Antenna'][ant]['Controller']
        # True if antenna is supposed to be tracking the source (no offsets)
        tracksrcflag.append((c['RAOffset'] + c['DecOffset'] + c['ElOffset'] + c['AzOffset']) == 0)
        az1 = extract(data,c['Azimuth1'])/10000.
        az_corr = extract(data,c['AzimuthPositionCorrected'])/10000.
        el1 = extract(data,c['Elevation1'])/10000.
        el_corr = extract(data,c['ElevationPositionCorrected'])/10000.
        rm = extract(data,c['RunMode'])
        if rm == 4:
            # Track mode
            az_req.append(extract(data,c['AzimuthVirtualAxis'])/10000.)
            el_req.append(extract(data,c['ElevationVirtualAxis'])/10000.)
        else:
            # All other modes
            az_req.append(extract(data,c['AzimuthPosition'])/10000.)
            el_req.append(extract(data,c['ElevationPosition'])/10000.)

        if rm == 1 or ant in [8,9,10,11,12]:    # New telescopes work differently
            # Position mode
            daz.append(az1 - az_corr)
            az_act.append(az_req[i] + daz[i])
            delv.append(el1 - el_corr)
            el_act.append(el_req[i] + delv[i])
        else:
            # All other modes
            daz.append(az1 - az_req[i])
            az_act.append(az1)
            delv.append(el1 - el_req[i])
            el_act.append(el1)

        if ant == 13:
            # Case of equatorial mount antennas, convert HA, Dec to El, Az
            eqel, eqaz = hadec2altaz(az_act[i]*dtor,el_act[i]*dtor)
            chi.append(par_angle(eqel, eqaz))
        else:
            chi.append(par_angle(el_act[i]*dtor,az_act[i]*dtor))

    daz = np.array(az_act) - np.array(az_req)
    # Track limit is set at 1/10th of primary beam at 18 GHz
    tracklim = np.array([0.0555]*13+[0.0043]*2)       # 15-element array
    trackflag = (np.abs(daz) <= tracklim) & (np.abs(np.array(delv)) <= tracklim)
    trackflag = np.append(trackflag,False)   # Ant 16 is never tracking

    return {'dAzimuth':daz,   'ActualAzimuth':np.array(az_act),  'RequestedAzimuth':np.array(az_req),
            'dElevation':np.array(delv),'ActualElevation':np.array(el_act),'RequestedElevation':np.array(el_req),
            'ParallacticAngle':np.array(chi)/dtor, 'TrackFlag':trackflag, 'TrackSrcFlag':tracksrcflag}

#============================
def par_angle(alt, az):
    '''Calculate the parallactic angle for a sky location given by
       altitude alt [radians] and azimuth az [radians].  This is the
       "nominal" parallactic angle for an X feed exactly aligned with
       the meridian.  It is defined to be the angle of the feed on the
       sky relative to the local line of constant hour angle, +ve east
       of north.

       This is likely reversed in sign for a feed at prime focus.
    '''
    dtor = np.pi/180.
    lat = 37.233170*dtor
    chi = np.arctan2(-np.cos(lat)*np.sin(az),
                  np.sin(lat)*np.cos(alt) - np.cos(lat)*np.sin(alt)*np.cos(az))
    return chi

#============================
def hadec2altaz(ha, dec):
    ''' Given an hour angle and declination, both in radians, return
        the corresponding altitude and azimuth for OVRO.
        
        This gives the same result as radec2azel() in coord_conv.py,
        but uses HA as input, and the order of the outputs is swapped.
    '''
    lat = 37.233170*np.pi/180.
    salt = np.sin(dec)*np.sin(lat) + np.cos(dec)*np.cos(lat)*np.cos(ha)
    alt = np.arcsin(salt)
    caz = (np.sin(dec) - np.sin(alt)*np.sin(lat)) / (np.cos(alt)*np.cos(lat))
    if type(caz) is np.ndarray:
        az = np.zeros(caz.shape,float)
        for i,c in enumerate(caz):
            if c >= 1 or c <= -1:
                az[i] = np.pi
            else:
                az[i] = np.arccos(c)
        if np.sin(ha[i]) > 0: 
            az[i] = 2*np.pi - az[i]
    else:
        if caz >= 1 or caz <= -1:
            az = np.pi
        else:
            az = np.arccos(caz)
        if np.sin(ha) > 0: return alt, 2*np.pi - az
    return alt, az

#============================
def azel_from_sqldict(sqldict, antlist=None):
    '''Given a dictionary read from a dimension-15 SQL stateframe query, calculate
       the actual and requested azimuth and elevation for each antenna, as well as the 
       difference between them, and a track flag, all as a dictionary of numpy float arrays.
       
       Added track source flag, which summarizes intentional offsets
    '''
    dtor = np.pi/180.
    if antlist is None:
        # No antlist, so assume all antennas
        antlist = range(15)

    az1 = copy.deepcopy(sqldict['Ante_Cont_Azimuth1'].astype('float'))/10000.
    az_corr = copy.deepcopy(sqldict['Ante_Cont_AzimuthPositionCorre'].astype('float'))/10000.
    el1 = copy.deepcopy(sqldict['Ante_Cont_Elevation1'].astype('float'))/10000.
    el_corr = copy.deepcopy(sqldict['Ante_Cont_ElevationPositionCor'].astype('float'))/10000.
    az_req = copy.deepcopy(sqldict['Ante_Cont_AzimuthPosition'].astype('float'))/10000.
    el_req = copy.deepcopy(sqldict['Ante_Cont_ElevationPosition'].astype('float'))/10000.
    # Use alternate source of requested positions where RunMode is 4
    rm = copy.deepcopy(sqldict['Ante_Cont_RunMode'].astype('int'))
    rms = rm.shape
    rm.shape = np.prod(rms)
    good = np.where(rm == 4)[0]
    if len(good) != 0:
        az_req_alt = copy.deepcopy(sqldict['Ante_Cont_AzimuthVirtualAxis'].astype('float'))/10000.
        el_req_alt = copy.deepcopy(sqldict['Ante_Cont_ElevationVirtualAxis'].astype('float'))/10000.
        az_req.shape = el_req.shape = az_req_alt.shape = el_req_alt.shape = np.prod(rms)
        az_req[good] = copy.deepcopy(az_req_alt[good])
        el_req[good] = copy.deepcopy(el_req_alt[good])
        az_req.shape = el_req.shape = rms
        
    daz = copy.deepcopy(az1 - az_req)
    az_act = copy.deepcopy(az1)
    delv = copy.deepcopy(el1 - el_req)
    el_act = copy.deepcopy(el1)
    # Set antenna 12 to RunMode 1 for this next selection, since
    # new S. Pole telescope works differently
    rm.shape = rms
    rm[:,11] = 1
    rm.shape = np.prod(rms)
    good = np.where(rm == 1)[0]
    if len(good) != 0:
        daz.shape = delv.shape = az_act.shape = el_act.shape = az1.shape = np.prod(rms)
        az_corr.shape = az_req.shape = el1.shape = el_corr.shape = el_req.shape = np.prod(rms)
        daz[good] = copy.deepcopy(az1[good] - az_corr[good])
        az_act[good] = copy.deepcopy(az_req[good] + daz[good])
        delv[good] = copy.deepcopy(el1[good] - el_corr[good])
        el_act[good] = copy.deepcopy(el_req[good] + delv[good])
        daz.shape = delv.shape = az_req.shape = el_req.shape = az_act.shape = el_act.shape = rms
    chi = par_angle(el_act*dtor,az_act*dtor)
    # Override equatorial antennas
    for iant in [8,9,10,12,13,14]:
        # Case of equatorial mount antennas, convert HA, Dec to El, Az
        eqel, eqaz = hadec2altaz(az_act[:,iant]*dtor,el_act[:,iant]*dtor)
        chi[:,iant] = par_angle(eqel, eqaz)

    daz = az_act - az_req
    # Track limit is set at 1/10th of primary beam at 18 GHz
    tracklim = np.array([0.0555]*13+[0.0043]*2)       # 15-element array
    trackflag = np.zeros(rms,'bool')
    for i in range(rms[0]):
        trackflag[i,:] = (np.abs(daz[i,:]) <= tracklim) & (np.abs(delv[i,:]) <= tracklim)

    trackflag[:,14] = False   # Ant 15 is never tracking
    
    # Check offsets to see if the antennas are intentionally not tracking the source
    tracksrcflag = np.ones(rms,bool)
    offsource = (sqldict['Ante_Cont_RAOffset'] + sqldict['Ante_Cont_DecOffset'] + sqldict['Ante_Cont_AzOffset'] +
                 sqldict['Ante_Cont_ElOffset']).nonzero()
    tracksrcflag[offsource] = False
    
    return {'dAzimuth':daz,   'ActualAzimuth':az_act,  'RequestedAzimuth':az_req,
            'dElevation':delv,'ActualElevation':el_act,'RequestedElevation':el_req,
            'ParallacticAngle':chi/dtor, 'TrackFlag':trackflag, 'TrackSrcFlag':tracksrcflag}
            
def PA_adjust(ant=None, crossed=False, offset_angle=0):
    ''' Spawned task to check the changing parallactic angle of given
        antenna and rotate the position angle of the focus rotation
        mechanism on Ant14 to counteract it.  Checks for Abort message
        once per second, and updates PA once per minute (if needed).
        
        This routine is invoked with $PA-TRACK command in the schedule,
        and aborts with $PA-STOP command, or if Ant14 is removed from
        the current subarray.
        
        Optional keyword:
          crossed    Boolean. If True, rotates the FRM to be 90-degrees
                       from the nominal parallactic angle. Default is False
    '''
    import time
    import adc_cal2
    if ant is None:
        q.put_nowait('No antenna specified. Exiting...')
        return
    accini = rd_ACCfile()
    acc = {'host': accini['host'], 'scdport':accini['scdport']}
    sf = accini['sf']
    sub1 = sf['LODM']['Subarray1']
    chikey = sf['Schedule']['Data']['Chi']
    timekey = sf['Schedule']['Data']['Timestamp']
    pakey = sf['FEMA']['FRMServo']['PositionAngle']['Position']
    while 1:
        # Read stateframe from ACC
        data, sfmsg = get_stateframe(accini)
        if extract(data,sf['Timestamp']) != 0 and extract(data,sub1) >> 13 == 0:
            # Stateframe has a valid Timestamp, and Antenna 14 is not in the subarray, so exit
            break
        if extract(data,timekey) != 0:
            # Do this only if stateframe timestamp is valid--otherwise just skip this update
            # Get Chi for this antenna from the stateframe, converted to degrees
            chi = extract(data,chikey)[ant]*180/np.pi
            
            # If the crossed keyword is set, the orientation angle is 90-degrees from chi
            if crossed: chi += 90
            
            # Make sure it is in range...
            if chi > 90.:
                chi = 180. - chi
            elif chi < -90:
                chi = 180. + chi
            pa_to_send = -np.int(chi)+offset_angle   # Desired rotation angle is -chi
            current_pa = np.int(extract(data,pakey)+0.5)+offset_angle
            if pa_to_send != current_pa:
                # Current PA is different from new one, so rotate feed to new position.
                adc_cal2.send_cmds(['frm-set-pa '+str(pa_to_send)+' ant14'],acc)
        # Sleep for 1 minute (but checking for Abort message every second), 
        # and then repeat
        for i in range(60):
            try:
                msg = q.get_nowait()
                if msg == 'Abort':
                    # Got abort message, so exit.
                    adc_cal2.send_cmds(['frm-set-pa 0 ant14'],acc)
                    return
            except:
                pass
            time.sleep(1)
    # To get here, either Ant14 is not in the subarray, or else we got 
    # an Abort message.  In either case, reset the PA to 0 and exit.
    adc_cal2.send_cmds(['frm-set-pa 0 ant14'],acc)
        
def PA_sweep(PA=80,rate=3):
    ''' Spawned task to rotate the 27-m focus rotation mechanism to
        value given by negative of PA argument (waits up to 2 minutes to 
        reach it), and then rotate the position angle at a rate given by 
        the rate argument (units = s/deg) until PA is reached. Checks 
        for Abort message once per second.
        
        This routine is invoked with $PA-SWEEP command in the schedule,
        and aborts with $PA-STOP command, or if Ant14 is removed from
        the current subarray.
        
        PA:  Initial PA is negative of this, and sweeps until PA is reached.
               Default = 80, for full sweep from -80 to 80
        rate: Rate of rotation, in seconds/degree. Default is 3, or 
               20 degrees/minute (can acquire and complete -80 to 80 
               sweep in 10 minutes)
    '''
    import time
    import adc_cal2
    if PA > 90:
        # Make sure PA is not too big
        PA = 90
    # Initial PA is negative of argument given
    pa_to_send = -PA
    accini = rd_ACCfile()
    acc = {'host': accini['host'], 'scdport':accini['scdport']}
    sf = accini['sf']
    sub1 = sf['LODM']['Subarray1']
    timekey = sf['Schedule']['Data']['Timestamp']
    pakey = sf['FEMA']['FRMServo']['PositionAngle']['Position']
    # Send FRM to initial position, and wait up to two minutes, checking every 5 s, until there
    adc_cal2.send_cmds(['frm-set-pa '+str(pa_to_send)+' ant14'],acc)
    current_pa = -999  # Start with impossible value for current_pa
    msg = ''
    for i in range(24):
        data, sfmsg = get_stateframe(accini)
        if sfmsg == 'No Error':
            if extract(data,sf['Timestamp']) != 0 and extract(data,sub1) >> 13 == 0:
                # Stateframe has a valid Timestamp, and Antenna 14 is not in the subarray, so exit
                msg = 'Abort'
                break
            if extract(data,timekey) != 0:
                current_pa = np.round(extract(data,pakey)+0.5)
            #print 'Current and target PAs:', current_pa, pa_to_send,
            if pa_to_send != current_pa:
                try:
                    msg = q.get_nowait()
                    if msg == 'Abort':
                        # Got abort message, so exit.
                        break
                except:
                    pass
                # Current PA is different from new one, so sleep 5 minutes
                #print 'Not equal, so sleeping 5 s'
                time.sleep(5)
            else:
                # We are on the desired PA, so proceed
                #print 'Target PA reached...'
                break
    # When we get here, either we are on the desired PA, or the 2 min is up, or
    # an Abort message was received.
    while 1:
        if msg == 'Abort':
            # Handle case of abort while positioning in above loop
            break
        # Read stateframe from ACC
        data, sfmsg = get_stateframe(accini)
        if sfmsg == 'No Error':
            if extract(data,sf['Timestamp']) != 0 and extract(data,sub1) >> 13 == 0:
                # Stateframe has a valid Timestamp, and Antenna 14 is not in the subarray, so exit
                break
        #print 'Sending command','frm-set-pa '+str(pa_to_send)+' ant14'
        adc_cal2.send_cmds(['frm-set-pa '+str(pa_to_send)+' ant14'],acc)
        # Sleep for number of seconds given by rate (but checking for Abort message every second), 
        # and then repeat
        for i in range(rate):
            try:
                msg = q.get_nowait()
                if msg == 'Abort':
                    # Got abort message, so exit.
                    break
            except:
                pass
            time.sleep(1)
        # Time to increment the PA
        pa_to_send += 1
        # If PA is reached, then exit
        if pa_to_send == PA:
            break
    # To get here, either Ant14 is not in the subarray, or we got 
    # an Abort message, or the FRM has reached PA.  In any of these cases, 
    # reset the PA to 0 and exit.
    adc_cal2.send_cmds(['frm-set-pa 0 ant14'],acc)
