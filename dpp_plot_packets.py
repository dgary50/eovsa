# History
# 2019-04-10  DG
#   Several changes to dpp_fix_packets() to that it never exits (CTRL-C to kill it),
#   and to fix a bug so that the timeout is honored correctly.  Also, the procstat()
#   routine is no longer called during a timeout.  Since the interfaces only need 
#   resetting every 15 minutes, the timeout is increased to 10 minutes.
# 2021-08-08  DG
#   Owen wrote alternative version, fix_packets2(), which rewrites the SMP_AFFINITY.sh
#   script on the fly when packets go to zero.  I added a couple of lines to remove
#   the empty /tmp/netplan_* folders created during a network reset.
# 2021-08-10  DG
#   Removing the files with a wild-card was not working, so now filenames are found
#   with glob() before removing by explicit name.

import numpy as np
import time
import glob
from util import Time
import subprocess

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
    #import  subprocess
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
                files = glob.glob('/tmp/netplan*')
                for file in files:
                    command = ['sudo','rm','-rf',file]  # Resetting leaves an empty folder, so delete it
                    proc = subprocess.Popen(command)
        else:
            timeout = max(timeout - 1, 0)
        if int(round(t) % 60) == 0:
            if timeout == 0:
                print Time.now().iso, np.round(if1), np.round(if2)
            else:
                print Time.now().iso, 'Waiting for end of timeout, now',timeout
            
        time.sleep(1 - (time.time() % 1))

def init_fix_packets(cpu):
    print "Reinitilizing fix_packets() using CPUs "+str(cpu[0]) +' and '+ str(cpu[1])
    
    #first need to edit the SMP_AFFINITY.sh
    smpaff='/home/user/test_svn/shell_scripts/SMP_AFFINITY.sh'
    f = open(smpaff,'r')
    lines = f.readlines()
    f.close()
    
    #comment out all the CPU lines
    for i,l in enumerate(lines):
        if '#(CPU' in l:
            if l[0] != '#':
                lines[i] = '#'+lines[i]
    
    #now uncomment the appropriate cpu lines
    for c in cpu:
        substring = '#(CPU '+str(c)+')'
        for i,l in enumerate(lines):
            if substring in l:
                if l[0] == '#':
                    lines[i] = lines[i][1:]
                    break
    
#    for l in lines:
#        print l.strip()
    
    f = open(smpaff,'w')
    for l in lines:
        f.write(l)
    f.close()
    
    #Now run the script
    command = ['/home/user/test_svn/shell_scripts/SMP_AFFINITY.sh']
    proc = subprocess.Popen(command)

def update_log(msg):
    f = open("/tmp/dpp_fix_packets_log.txt",'a')
    f.write(Time.now().iso[:19]+": "+msg+'\n')
    f.close()
    
def fix_packets2():
    ''' Continuously checks /proc/net/softnet_stat once each second,
        and calculates the number of packets/s on each interface.
        
        If the number of packets on either interface drops below 130,000,
        it sends a command to reset the interfaces.
    '''
    
    #set up a list of cpu pairs from cpu 18 to 23
    cpu=[]
    for i in range(18,23,2):
        cpu.append([i,i+1])
    
    reboot_required = False
    started = False
    while 1:
        for i,c in enumerate(cpu):
            init_fix_packets(c)
        
            time.sleep(1 - (time.time() % 1))
            val0 = procstat(c[0],c[1])
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
                    val = procstat(c[0],c[1])
                    if1 = (val-val0)[0]/(t-tpre)
                    if2 = (val-val0)[1]/(t-tpre)
                    val0 = val
                    tpre = t
                    if if1 == 0.0 or if2 == 0.0:
                        if i == 0:
                            if started == False:
                                started = True
                            else:
                                reboot_required = True
                        break
                        
                    if (100000 < if1 < 130000) or (100000 < if2 < 130000):
                        print Time.now().iso,'Packet loss detected!', np.round(if1), np.round(if2), 'Resetting interfaces'
                        update_log('Packet loss detected! ' + str(np.round(if1)) + " " + str(np.round(if2)) + ' Resetting interfaces')
                        command = ['sudo','/usr/sbin/netplan','apply']
                        proc = subprocess.Popen(command)
                        timeout = 600   # Leave a 10-min window to avoid resetting too often
                        files = glob.glob('/tmp/netplan*')
                        for file in files:
                            command = ['sudo','rm','-rf',file]  # Resetting leaves an empty folder, so delete it
                            proc = subprocess.Popen(command)
                            update_log('Removed '+file)
                else:
                    timeout = max(timeout - 1, 0)
                if int(round(t) % 60) == 0:
                    if timeout == 0:
                        print Time.now().iso, np.round(if1), np.round(if2)
                    else:
                        print Time.now().iso, 'Waiting for end of timeout, now',timeout
            
                time.sleep(1 - (time.time() % 1))
        
        update_log('Resetting CPUs')
        if reboot_required:
            print "DPP REBOOT IS REQUIRED!"
            update_log("DPP REBOOT IS REQUIRED!")
