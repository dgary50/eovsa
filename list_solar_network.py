'''
   Routines for searching the solar.pvt network and providing information
   about the active IP addresses.'''
   
# History:
#   2015-Mar-03  DG
#      Added the list_unassigned_dhcp() routine, which lists the IP
#      addresses in the top half of the solar subnet.

import subprocess

def list_solar_network(mac=False):
    # Uses NSLOOKUP to find all defined machines on the 192.168.24.* subnet and
    # lists the DNS name.  Optionally pings and looks at /proc/net/arp file to
    # determine MAC address.
    for i in range(1,129):
        ipaddr = "192.168.24."+str(i)
        res = subprocess.check_output("nslookup "+ipaddr,shell=True)
        if res.find('name') != -1:
            name = res[res.find('name')+7:res.find('solar')+9]
            print ipaddr+':',name,
            if mac:
                try:
                    res = subprocess.check_output("ping -c 1 "+ipaddr,shell=True)
                    f = open('/proc/net/arp','r')
                    for line in f.readlines():
                        if line.find(ipaddr) != -1:
                            print ' MAC:',line[41:58]
                            break
                except:
                    print ' MAC: None--ping failed.'
            else:
                print '\n'

def list_nonassigned_dhcp():
    # Uses ping to find active IP addresses in the "non-assigned" range of
    # solar IP addresses 192.168.24.129-255.
    for i in range(129,255):
        ipaddr = "192.168.24."+str(i)
        try:
            res = subprocess.check_output("ping -c 1 "+ipaddr,shell=True)
            if res.find('transmitted') != -1:
                print ' IP: '+ipaddr+' found.'
            else:
                print 'Unexpected result:\n',res
        except:
            print ' IP: '+ipaddr+' --ping failed.'
