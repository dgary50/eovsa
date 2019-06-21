# History
# 2019-04-10  DG
#   Several changes to dpp_fix_packets() to that it never exits (CTRL-C to kill it),
#   and to fix a bug so that the timeout is honored correctly.  Also, the procstat()
#   routine is no longer called during a timeout.  Since the interfaces only need 
#   resetting every 15 minutes, the timeout is increased to 10 minutes.

import numpy as np
import time
from util import Time

def procstat(k1,k2):
    f = open('/proc/net/softnet_stat','r')
    lines = f.readlines()
    f.close()
    return np.array([int(lines[k1].split()[0],16),int(lines[k2].split()[0],16)])

def plot_packets(n, cpu=[22,23], ax=None):
    import matplotlib.pylab as plt
    import subprocess
    if ax is None:
        f, ax = plt.subplots(1,1)
    time.sleep(1 - (time.time() % 1))
    val0 = procstat(cpu[0],cpu[1])
    time.sleep(1 - (time.time() % 1))
    t0 = time.time()
    tiso = Time.now().iso
    ax.set_title(tiso)
    ax.set_ylim(0,200000)
    ax.set_xlim(-10,n+10)
    ax.set_xlabel('Time [s after '+tiso+']')
    ax.set_ylabel('Packets/s')
    z = np.zeros(n,dtype=float)
    ld1 = z*np.nan
    ld2 = z*np.nan
    line1, = ax.plot(z,ld1,'x',color='C0')
    line2, = ax.plot(z,ld2,'.',color='C1')
    t = np.zeros(n+1,dtype=float)
    t[0] = -1.0
    i = 0
    timeout = 5
    while t[i] < n:
        try:
            t[i+1] = time.time() - t0
        except:
            print 'All done.'
            return
        val = procstat(cpu[0],cpu[1])
        if1 = (val-val0)[0]/(t[i+1]-t[i])  # Packets/s on interface 1
        if2 = (val-val0)[1]/(t[i+1]-t[i])  # Packets/s on interface 2
        ld1[i] = if1
        ld2[i] = if2
        if timeout == 0 and (100000 < if1 < 130000) or (100000 < if2 < 130000):
            print Time.now().iso,'Packet loss detected!', if1, if2, 'Resetting interfaces'
            command = ['sudo','/usr/sbin/netplan','apply']
            proc = subprocess.Popen(command)
            timeout = 5   # Leave a 5-s window to avoid resetting too often
        else:
            timeout = max(timeout - 1, 0)
        line1.set_data(t[1:],ld1)
        line2.set_data(t[1:],ld2)
        val0 = val
        plt.pause(0.001)
        i += 1
        time.sleep(1 - (time.time() % 1))
        
def fix_packets(cpu=[22,23]):
    ''' Continuously checks /proc/net/softnet_stat once each second,
        and calculates the number of packets/s on each interface.
        
        If the number of packets on either interface drops below 130,000,
        it sends a command to reset the interfaces.
    '''
    import  subprocess
    time.sleep(1 - (time.time() % 1))
    val0 = procstat(cpu[0],cpu[1])
    time.sleep(1 - (time.time() % 1))
    t0 = time.time()
    tiso = Time.now().iso
    i = 0
    timeout = 5   # Leave a 5-s window to avoid resetting too often
    print Time.now().iso, 'Started...'
    tpre = -1
    while 1:
        t = time.time() - t0
        if timeout == 0:
            val = procstat(cpu[0],cpu[1])
            if1 = (val-val0)[0]/(t-tpre)
            if2 = (val-val0)[1]/(t-tpre)
            val0 = val
            tpre = t
            if (100000 < if1 < 130000) or (100000 < if2 < 130000):
                print Time.now().iso,'Packet loss detected!', np.round(if1), np.round(if2), 'Resetting interfaces'
                command = ['sudo','/usr/sbin/netplan','apply']
                proc = subprocess.Popen(command)
                timeout = 600   # Leave a 10-min window to avoid resetting too often
        else:
            timeout = max(timeout - 1, 0)
        if int(round(t) % 60) == 0:
            if timeout == 0:
                print Time.now().iso, np.round(if1), np.round(if2)
            else:
                print Time.now().iso, 'Waiting for end of timeout, now',timeout
            
        time.sleep(1 - (time.time() % 1))
    
