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
# 2021-08-25  OG
#   Modified update_log so that it writes to the /common/webplots/dpp_fix_packets_log.txt
#   file. This is so it can be accessed by the sf_display program. It only saves the
#   last 5 lines. This can be changed if necessary.

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
        
def init_fix_packets(cpu):
    #print "Reinitilizing fix_packets() using CPUs "+str(cpu[0]) +' and '+ str(cpu[1])
    
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
        
    f = open(smpaff,'w')
    for l in lines:
        f.write(l)
    f.close()
    
    #Now run the script
    command = ['bash','/home/user/test_svn/shell_scripts/SMP_AFFINITY.sh']
    proc = subprocess.Popen(command)

def update_log(msg):
    '''Reads each line of the file /common/webplots/dpp_fix_packets_log.txt into
    a list and appends the msg to the list. Only the last 20 lines are written
    back to the file'''
    
    logfile="/common/webplots/dpp_fix_packets_log.txt"
    f = open(logfile,'r')
    lines = f.readlines()
    f.close()
    lines.append(Time.now().iso[:19]+": "+msg+'\n')
    if len(lines) > 20:
        lines = lines[-20:]
    
    f = open(logfile,'w')
    for l in lines:
        f.write(l)
    f.close()

def update_cpu_file(cpu):
    '''writes out the cpu's currently being used to /common/webplots/dpp_cpu.txt'''
    
    cpufile="/common/webplots/dpp_cpu.txt"
    f = open(cpufile,'w')
    for c in cpu:
        f.write(str(c)+'\n')
        
    f.close()
    
#Continuously checks /proc/net/softnet_stat once each second,
#and calculates the number of packets/s on each interface.
#        
#If the number of packets on either interface drops below 130,000,
#it sends a command to reset the interfaces.
    
#set up a list of cpu pairs from cpu 18 to 23
cpu=[]
for i in range(18,23,2):
    cpu.append([i,i+1])
    
reboot_required = False
started = False
update_log('dpp_fix_packets.py has been restarted')
while 1:
    for i,c in enumerate(cpu):
        update_cpu_file(c)
        update_log('Resetting CPUs - Using CPUs '+str(c[0])+' and '+str(c[1])+'.')
        if reboot_required:
            update_log('DPP REBOOT REQUIRED!')
            
        init_fix_packets(c)
        time.sleep(1 - (time.time() % 1))
        val0 = procstat(c[0],c[1])
        time.sleep(1 - (time.time() % 1))
        t0 = time.time()
        tiso = Time.now().iso
        i = 0
        timeout = 5   # Leave a 5-s window to avoid resetting too often
        #print Time.now().iso, 'Started...'
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
                    # ==== This section can be removed once testing is complete ====
                    print Time.now().iso,'Packet loss detected!', np.round(if1), np.round(if2), 'Resetting interfaces'
                    update_log('Packet loss detected! ' + str(np.round(if1)) + " " + str(np.round(if2)) + ' Resetting interfaces')
                    # ==============================================================
                    command = ['sudo','/usr/sbin/netplan','apply']
                    proc = subprocess.Popen(command)
                    timeout = 600   # Leave a 10-min window to avoid resetting too often
                    files = glob.glob('/tmp/netplan*')
                    for file in files:
                        command = ['sudo','rm','-rf',file]  # Resetting leaves an empty folder, so delete it
                        proc = subprocess.Popen(command)
                        #update_log('Removed '+file)
            else:
                timeout = max(timeout - 1, 0)
            
            # ==== This section can be removed once testing is complete ====
            if int(round(t) % 60) == 0:
                if timeout == 0:
                    print Time.now().iso, np.round(if1), np.round(if2)
                else:
                    print Time.now().iso, 'Waiting for end of timeout, now',timeout
            #        update_log('Waiting for end of timeout, now '+str(timeout))
            # ==============================================================
            time.sleep(1 - (time.time() % 1))
