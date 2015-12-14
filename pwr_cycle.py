#!/usr/bin/env python

import telnetlib, time

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

