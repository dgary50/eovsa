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
#

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
    return ovro_dict

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
        f = urllib2.urlopen(url,timeout=0.2)
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
            solpwr.update({'TransformerTemp':int(lines[i].split('&deg;C')[0][-2])})
        elif line.find('FET Temperature') >0:
            solpwr.update({'FETTemp':int(lines[i].split('&deg;C')[0][-2])})

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
        return None, 'Cannot connect to port '+str(accini['sfport'])

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
    dtor = np.pi/180.
    if antlist is None:
        # No antlist, so assume all antennas
        antlist = range(15)

    for i, ant in enumerate(antlist):
        c = sf['Antenna'][ant]['Controller']
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

        if rm == 1 or ant == 11:    # New S. Pole telescope works differently
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

        chi.append(par_angle(el_act[i]*dtor,az_act[i]*dtor))

    # Track limit is set at 1/10th of primary beam at 18 GHz
    tracklim = np.array([0.0555]*13+[0.0043]*2)       # 15-element array
    trackflag = (np.array(daz) <= tracklim) & (np.array(delv) <= tracklim)
    trackflag = np.append(trackflag,False)   # Ant 16 is never tracking

    return {'dAzimuth':np.array(daz),   'ActualAzimuth':np.array(az_act),  'RequestedAzimuth':np.array(az_req),
            'dElevation':np.array(delv),'ActualElevation':np.array(el_act),'RequestedElevation':np.array(el_req),
            'ParallacticAngle':np.array(chi)/dtor, 'TrackFlag':trackflag}

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
        the corresponding altitude and azimuth for OVRO
    '''
    lat = 37.233170*pi/180.
    salt = sin(dec)*sin(lat) + cos(dec)*cos(lat)*cos(ha)
    alt = arcsin(salt)
    caz = (sin(dec) - sin(alt)*sin(lat)) / (cos(alt)*cos(lat))
    az = arccos(caz)
    if sin(ha) > 0: return alt, 2*pi - az
    return alt, az

#============================
def azel_from_sqldict(sqldict, antlist=None):
    '''Given a dictionary read from a dimension-15 SQL stateframe query, calculate
       the actual and requested azimuth and elevation for each antenna, as well as the 
       difference between them, and a track flag, all as a dictionary of numpy float arrays.
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
        daz.shape = delv.shape = az_act.shape = el_act.shape = rms
    chi = par_angle(el_act*dtor,az_act*dtor)

    # Track limit is set at 1/10th of primary beam at 18 GHz
    tracklim = np.array([0.0555]*13+[0.0043]*2)       # 15-element array
    trackflag = np.zeros(rms,'bool')
    for i in range(rms[0]):
        trackflag[i,:] = (daz[i,:] <= tracklim) & (delv[i,:] <= tracklim)

    return {'dAzimuth':daz,   'ActualAzimuth':az_act,  'RequestedAzimuth':az_req,
            'dElevation':delv,'ActualElevation':el_act,'RequestedElevation':el_req,
            'ParallacticAngle':chi/dtor, 'TrackFlag':trackflag}
