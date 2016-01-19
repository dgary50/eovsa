#!/usr/bin/env python
#
# History:
#   2016-Jan-16  DG
#     Started this history log.  Added ant_toggle() routine to toggle the
#     power of various devices in the field at each antenna.
#   2016-Jan-19  DG
#     Added Queue use to ant_toggle so that calling program (schedule.py)
#     can receive messages when spawned by threading module.

import requests
from requests.auth import HTTPDigestAuth
import telnetlib, time
import Queue

q = Queue.Queue()

def ant_toggle(antnum, device=None, wait=None, cycle=True):
    ''' Toggles power to one of the devices attached to the Viking Relay switch
        in each antenna controller box.  If cycle=True, then end state is ALWAYS
        for the device to be turned ON (relay turned OFF).
        
        The possible devices are 'antenna' or 'ant', 'frontend' or 'fem', 'crio'.
        There is also a relay 4 that can be switched by specifying 'other',
        although this will not affect anything unless a hardware change is made.
        These names are not case sensitive.
        
        The requested device is powered off for 15 s, then back on, unless
        wait is set to a number of seconds to wait.
        
        All messages go into Queue q, which can be read by the calling program.
    '''
    relaydef = {'ANTENNA':1, 'ANT':1, 'FEM':2, 'FRONTEND':2, 'CRIO':3, 'OTHER':4}
    devstr = {1: 'Antenna', 2: 'Frontend', 3: 'CRIO', 4: 'Relay 4'}
    if device is None:
        # Default to cycling the antenna controller
        relay = 1
    else:
        try:
            relay = relaydef[device.upper()]
        except:
            q.put('Ant'+str(antnum)+' Error interpreting device '+device)
            return
        
    # Try to log in 3 times, with 1 s delay between each
    for i in range(3):
        url = 'http://vik'+str(antnum)+'.solar.pvt/protect/'
        try:
            q.put('Attempt '+str(i+1)+' to login.')
            test = requests.get(url,auth=HTTPDigestAuth('admin','pwr4me'))
        except requests.ConnectionError as e:
            q.put('Ant'+str(antnum)+' Error could not connect to '+url+'.  Message: '+e.message.message)
            return
        if test.status_code == 200:
            q.put('Ant'+str(antnum)+' Login successful.')
            break
        else:
            q.put('Ant'+str(antnum)+' Login failed with status code: '+str(test.status_code))
        time.sleep(1)
    if test.status_code != 200:
        q.put('Ant'+str(antnum)+' Error HTML status code: '+str(test.status_code))
        return
        
    url = 'http://vik'+str(antnum)+'.solar.pvt/protect/relays.cgi?relay='+str(relay)+'&state=toggle'
    dur = 15
    if wait:
        if type(wait) != int:
            q.put('Ant'+str(antnum)+' Warning: Could not interpret wait duration',wait,'.  Must be an integer type.  Will use 15 s')
            dur = 15
        else:
            dur = wait
    try:
        r = requests.get(url,auth=HTTPDigestAuth('admin','pwr4me'))
        if cycle == True:
            if r.text == devstr[relay]+' now on':
                # Request to turn relay on (to turn device off) worked,
                # so wait for requested length of time and then toggle
                # again to turn device back on.
                time.sleep(dur)
                r = requests.get(url,auth=HTTPDigestAuth('admin','pwr4me'))
        # Close the connection
        r.close()
        # This fixes the peculiarity that the relay is off when power is on
        # and vice versa.  The "result" gives the state of the power, not the relay
        if r.text[-3:] == 'off':
            q.put('Ant'+str(antnum)+' '+r.text.replace('off','on'))
        else:
            q.put('Ant'+str(antnum)+' '+r.text.replace('off','on'))
        return
    except:
        q.put('Ant'+str(antnum)+' Error communicating with Viking Relay '+devstr)
        return

def pwr_cycle(host,loadn,user='admin',passwd='pwr4me',wait=None):
    ''' Connect to a Tripp Lite PDU (Power Distribution Unit) and power cycle one of the loads.
           host    One of the PDUs, one of 'pdunetwork.solar.pvt', 'pduanalog.solar.pvt', or 'pdudigital.solar.pvt'
           user    Username, currently set to 'admin' for all three PDUs
           passwd  Password
           loadn   Integer load number to power cycle
           wait    Duration [s] to wait before turning load back on 
                     if omitted, "Cycle" is used instead of "Off" followed by "On
        Returns True if successful, False otherwise
    '''
    if wait is None:
        #           Username, Password,    Devices, Device, Actions, Loads, Load #,         Cycle, X,    X,    X,    Logout
        term_str = ['login:', 'Password:', '> ',    '> ',   '> ',    '> ',  '> ',           '> ',  '> ', '> ', '> ', '> ']
        response = [user+'\n',passwd+'\n', '1\n',   '1\n',  '2\n',   '2\n', str(loadn)+'\n','3\n', 'X\n','X\n','X\n','X\n']
        dur = False
    else:
        #           Username, Password,    Devices, Device, Actions, Loads, Load #,         Off,   On,   X,    X,    X,    Logout
        term_str = ['login:', 'Password:', '> ',    '> ',   '> ',    '> ',  '> ',           '> ',  '> ', '> ', '> ', '> ', '> ']
        response = [user+'\n',passwd+'\n', '1\n',   '1\n',  '2\n',   '2\n', str(loadn)+'\n','2\n', '2\n','X\n','X\n','X\n','X\n']
        if type(wait) != int:
            print 'Warning: Could not interpret wait duration',wait,'.  Must be an integer type.  Will use 30 s'
            dur = 30
        else:
            dur = wait

    # Initiate connection
    tn = telnetlib.Telnet(host)
    time.sleep(1)

    # Loop over response
    for i in range(len(term_str)):
        out = tn.read_until(term_str[i],1)
        if out[-len(term_str[i]):] != term_str[i]:
            print 'Telnet connection to',host,'timed out.'
            tn.close()
            return False
            break
        if dur and i == 8:
            # Wait dur [s] before turning load back on
            time.sleep(dur)
        tn.write(response[i])
    tn.close()
    return True

def pwr_off(host,user,passwd,loadn):
    ''' Connect to a Tripp Lite PDU (Power Distribution Unit) and power down one of the loads.
           host    One of the PDUs, one of 'pdunetwork.solar.pvt', 'pduanalog.solar.pvt', or 'pdudigital.solar.pvt'
           user    Username, currently set to 'admin' for all three PDUs
           passwd  Password
           loadn   Integer load number to power down
        Returns True if successful, False otherwise
    '''
    #           Username, Password,    Devices, Device, Actions, Loads, Load #,         Off,   X,    X,    X,    Logout
    term_str = ['login:', 'Password:', '> ',    '> ',   '> ',    '> ',  '> ',           '> ',  '> ', '> ', '> ', '> ']
    response = [user+'\n',passwd+'\n', '1\n',   '1\n',  '2\n',   '2\n', str(loadn)+'\n','2\n', 'X\n','X\n','X\n','X\n']

    # Initiate connection
    tn = telnetlib.Telnet(host)
    time.sleep(1)

    # Loop over response
    for i in range(len(term_str)):
        out = tn.read_until(term_str[i],1)
        if out[-len(term_str[i]):] != term_str[i]:
            print 'Telnet connection to',host,'timed out.'
            tn.close()
            return False
            break
        tn.write(response[i])
    tn.close()
    return True

if __name__ == "__main__":
    pwr_cycle('pdudigital.solar.pvt',14)

