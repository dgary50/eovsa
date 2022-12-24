#
# ADC_PLOT reads the Roach ADC boards [adc_monitor()] and records the results
#          in a file in the /tmp folder named ADCplotYYYYMMDD_HHMM.npy, and can
#          also plot [adc_plot()] the results of such a file.
#
# Optionally uses the results to make a DCM Master Table for DCM attn calibration
# (if attnval is set to a specific value)
#
# History
#
#  Originally written by Owen Giersch.
#  2021-07-14 DG
#    Modified to its current form.
#  2021-07-15 DG
#    Added display of FEMATTN level to plot
#  2021-07-16  DG
#    Lots of development culminating in actually sending FEMATTN level updates
#    to keep the ADC level properly set.  It seems to work well, but needs some
#    easier way to start it and keep it running full time.
#  2021-07-24  DG
#    Removed all plotting from the ADC monitoring routine and just save the gathered information
#    as a .npy file for reading and plotting asynchronously.  The routine is renamed adc_monitor().
#    Then added a new adc_plot() routine to read that file and do the plotting.
#  2021-07-26  DG
#    Finally made this callable via the command line so that it can run as a cron job.
#  2021-07-27  DG
#    Now determines the number of 60 s loops to do from the current time.  Will exit
#    if the current time is in the gap 03:00-11:00 UT, unless a number of loops is
#    explicitly given.
#  2021-08-02  DG
#    Added code to use the previously saved levels in NormalObserving mode when a
#    new scan is started.  This should minimize the ramp-up of levels as the system
#    adjusts to the Sun over several minutes.
#  2021-08-05  DG
#    I keep having crashes due to the handling of savlevs.  One more try, avoiding
#    starting with savlevs as None.
#  2021-08-08  DG
#    Change the way adc_monitor() writes to the log file, so that it gets flushed once/min.
#    It also writes a status line every minute.
#  2021-08-09  DG
#    Change adc_plot's exit strategy--reading the file needs to fail three times
#    in a row before it exits.
#  2021-08-11  DG
#    Make adc_monitor() tolerant of failures to read from SQL
#  2021-08-16  DG
#    The adc_monitor() is working, but we are getting problems with the ACC and with
#    the schedule hanging.  Try adding a 1-s pause between FEMATTN commands.
#  2021-08-19  DG
#    This is a new version with extensive changes to eliminate most reading from SQL.  All
#    threads now read the stateframe directly from the ACC, which should be more robust
#    except for ACC failures.
#  2021-09-03  DG
#    Added option to replay ADC data from a file by specifying a start time in adc_plot().
#    The replay proceeds at a 2 s / minute rate.
#  2022-01-27  DG
#    Added FLARETEST project ID (in addition to NormalObserving) for adjusting attenuation.
#  2022-Mar-07  DG
#    Changed due to loss of SQL.  Blocks changes are preceded by
#            # ************ This block commented out due to loss of SQL **************
#
import roach as r
import struct
import numpy as np
import matplotlib.pyplot as plt
import threading
import adc_cal2 as adc2
import dbutil as db
from time import sleep,time
from util import Time
from copy import copy
import sys
import stateframe as stf

#set up threading
class AGC_Thread (threading.Thread):
    def __init__(self, threadID):
        threading.Thread.__init__(self)
        self.threadID = threadID                               #set thread ID
        self.name = "roach" + str(threadID)                    #set the thread name
        self.sd = np.empty((4,50),float)                       #set up the standard deviation array
        self.levels = np.zeros((4), int)
        self.agc = np.zeros((2), int)
        self.stop = False
        self.lock = threading.Lock()
        
    def get_FEM_level(self):
        # Read FEM level for current moment from ACC stateframe
        rn = int(self.name[5:])
        a1 = (rn-1)*2
        a2 = a1+1
        accini = self.accini
        buf, msg = stf.get_stateframe(accini)
        if msg == 'No Error':
            hpol1 = stf.extract(buf,accini['sf']['Antenna'][a1]['Frontend']['FEM']['HPol']['Regime']['Level'])
            vpol1 = stf.extract(buf,accini['sf']['Antenna'][a1]['Frontend']['FEM']['VPol']['Regime']['Level'])
            hpol2 = stf.extract(buf,accini['sf']['Antenna'][a2]['Frontend']['FEM']['HPol']['Regime']['Level'])
            vpol2 = stf.extract(buf,accini['sf']['Antenna'][a2]['Frontend']['FEM']['VPol']['Regime']['Level'])
            agc1 = stf.extract(buf,accini['sf']['Antenna'][a1]['Frontend']['FEM']['AGC']['Active'])
            agc2 = stf.extract(buf,accini['sf']['Antenna'][a2]['Frontend']['FEM']['AGC']['Active'])
            return [hpol1, vpol1, hpol2, vpol2], [agc1, agc2]
        else:
            return [-1, -1, -1, -1], [0, 0]
            
    def run(self):
        self.lock.acquire()                                    #acquire the lock
        while not self.stop:
            start = time()                                     #used for computing execution time
            f, s = np.modf(start)
            waittime = (np.floor(s/60)+1)*60 - start
            sleep(waittime)
            start=time()
            outdata, levels, agc = self.grab_all()
            self.sd = np.std(outdata,axis=2)           #calculate standard deviation from the roach
            self.levels = levels  # Current FEM levels, one for each channel of this ROACH board
            self.agc = agc  # Whether AGC is active, one for each antenna of this ROACH board
            if self.verbose:
                print self.name+" Execution Time= "+str(time()-start)   #display the execution time
        self.lock.release()
            
    def grab_all(self):
        roachname = self.name
        buf = ''                                               #adc buffer for single slot and channel
        slots = []
        levs = []
        rn = r.Roach(roachname)

        # Grab 10 slots in one channel in one second (every fifth slot) starting at the given slot.
        # This will be called five times starting with slot 0, 1, 2, 3, and 4.
        def grab10(rn,slot,chan):
            ''' Returns the appended raw byte buffer for all 10 slots, and the list of slots
                actually measured.
            '''
            f, s = np.modf(time())
            twait = 1-f+0.02*slot-0.01
            if twait > 0:
                sleep(twait)
            frac = []
            rn.fpga.write_int('swreg_snap_select',chan)
            buf = ''
            for i in range(10):
                desired_frac = (i*5 + slot)*0.02
                rn.fpga.write_int('adc_data_adc_ctrl',0)
                tcapture = time()
                rn.fpga.write_int('adc_data_adc_ctrl',7)
                buf += rn.fpga.read('adc_data_adc_bram', 2048*4, 0)
                f, s = np.modf(tcapture)
                frac.append(f)
                twait = s + desired_frac + 0.0938 - time()
                if twait > 0: sleep(twait)
            frac = np.array(frac)
            return buf, (frac/0.02).astype(int)

        agc = np.zeros(2,int)
        for chan in range(4):
            for slot in range(5):
               buf1, s = grab10(rn, slot, chan)
               buf += buf1
               slots.append(s)
            lvs, agcvals = self.get_FEM_level()
            levs.append(lvs)
            agc += np.array(agcvals)
        nbytes = 8192*4*50
        udata = np.array(struct.unpack('>'+str(nbytes)+'b',buf))
        # Reorder axes and reshape
        udata.shape = (4,50,8192)
        levs = np.array(levs)  # Shape (4,4)
        levels = np.zeros(4,int)
        # Get median over level measurements, but handle case where some measurements
        # failed on reading the stateframe
        for i in range(4):
            good, = np.where(levs[:,i] != -1)
            if len(good) == 0:
                levels[i] = -1     # All measurements were bad (ACC may be down)
            else:
                levels[i] = np.median(levs[good,i])
        outdata = np.zeros((4,50,8192),float)     #numpy array of data to return in form of [chan,slot,data]
        # The data are a bit scrambled in udata, and some bands may be missing (and others measured twice)
        # although hopefully that will be rare.  This should unscramble them.
        slots = np.array(slots)  # Shape (4,5,10)
        slots.shape = (4,50)
#        def find_missing(x):
#            bad = []
#            ptr = 0
#            for i in range(50):
#                try:
#                    if x[ptr] != i:
#                        bad.append(i)
#                    else:
#                        ptr += 1
#                except IndexError:
#                    bad.append(i)
#            return bad

        sleep(int(roachname[-1])*0.02)
        for c in range(4):
            slotarray, idx = np.unique(slots[c],return_index=True)
            if self.verbose:
                print roachname, c, '   Missing:',50-len(idx)
#            if len(idx) != 50:
#                print roachname,c,find_missing(slotarray)
            outdata[c,slotarray] = udata[c,idx]
            
        rn.fpga.stop()
        
        return outdata, levels, agc

def adc_monitor(nloop=None, verbose=False):
    ''' Performs an ADC measurement on all ROACH boards nloop times.
        
        Inputs:
            nloop    The number of times to make the measurement (once per minute)
                       If not provided, the number of loops is determined from the
                       current time such that it stops at 03:00 UT, and will exit
                       if the current time is between 03:00 - 11:00 UT.  If you
                       really want it to run during this gap, provide a value for nloop.
            verbose  If True, quite a lot of information will be printed
                       to the terminal
                       
        Output is written to /tmp/ADCplot<yyyymmdd_hhmm>.npy as a series
        of dictionaries, each with keys:
            time      time of record, as an ISO string
            fseqfile  the filename of the frequency sequence in progress (string)
            projid    the project ID in progress (string)
            levels    the FEM levels at the time, of size (7, 4)
            needs     the FEM level increment needed, of size (7, 4)
            stdev     the measured standard deviations, of size (7, 4, 50)
            
        NB: This routine also adjusts the front end gain (FEMATTN level) if in 
            NormalObserving mode, to maintain the optimal ADC levels.
    '''

    # Convert current time to string suitable for a file name <yyyymmdd_hhmm>
    t = Time.now().iso[:16].replace('-','').replace(':','').replace(' ','_')
    # Save the data as a npy file.
    f2name = '/tmp/ADCplot'+t+'.npy'
    sys.stdout = open('/tmp/adc_plot'+t+'.log','w')
    sys.stdout.write('Output saved in file '+f2name+'\n')

    tt=time()
    if not nloop:
        mjdnow = Time.now().mjd
        # Find difference between now and 03:00 UT
        dmjd = int(mjdnow) + 0.125 - mjdnow
        if dmjd <= -0.3333:
            # Time is at least 10:59:30 UT, so it is safe to start
            dmjd += 1
            mjdstop = mjdnow + dmjd  # Future time at which to stop and exit
        nloop = int(dmjd*24*60)
        if nloop < 0:
            sys.stdout.write('Current time is in the gap 03:00-11:00 UT.  Will exit.\n')
            exit()
    else:
        mjdstop = Time.now().mjd + nloop/24./60.
    #list of threads
    threads = []

    accini = stf.rd_ACCfile()
    acc = {'host': accini['host'], 'scdport':accini['scdport']}
    
    #set up the threads
    for t in range(1,8):            
        threads.append(AGC_Thread(t))
        threads[-1].verbose = verbose
        threads[-1].accini = accini

    #Start new Threads 
    for t in threads:
        t.start()

    # Give threads some time to work
    sleep(45)

    savlevs = np.zeros((7,4),int)
    newlevs = np.zeros((7,4),int)

    while Time.now().mjd < mjdstop:
        t0 = time()
        # Read once per minute on the 45 s
        waittime = 45 - t0 % 60
        if waittime < 0: waittime += 60
        if verbose: 
            print "waittime: "+str(waittime)
            sys.stdout.flush()

        sleep(waittime)
        t = Time.now().iso
        if verbose: 
            print t
            sys.stdout.flush()
        
    #    stdev[0:-1]=stdev[1:]
    #    for i,t in enumerate(threads):
    #        stdev[-1,i]=t.sd

        # Grab current stateframe
        buf, msg = stf.get_stateframe(accini)
        if msg == 'No Error':
            fseqfile = stf.extract(buf,accini['sf']['LODM']['LO1A']['FSeqFile']).strip('\0')
        else:
            fseqfile = 'Unknown'

        # Loop over threads and get sd (4,50) for each
        stdall = np.zeros((7,4,50),float)
        needs = np.zeros((7,4), float)
        agc = np.zeros((7,2), int)
        levs = np.zeros((7,4), int)
        for i in range(7):
            stdev = threads[i].sd
            levs[i] = threads[i].levels
            agc[i] = threads[i].agc
            stdall[i] = stdev
            for chan in range(4):
                # Use highest half of points to calculate FEM level needed to achieve target
                needs[i,chan] = np.log10(np.median(np.sort(stdev[chan])[24:]/34.))*5.
        agc.shape = (14,)
        # Save stdevs for this time
#        fh = open(fname,'ab')
#        fh.write(stdall)
#        fh.close()

        # ************ This block commented out due to loss of SQL **************
        ## Only do this if we are in a "Normal Observing" scan (as read from latest SQL header)
        #cursor = db.get_cursor()
        #ver = db.find_table_version(cursor, Time.now().lv, scan_header=True)
        #query = 'select top 1 * from hV'+ver+'_vD1 order by Timestamp desc'
        #result, msg = db.do_query(cursor, query)
        #cursor.close()
        #if msg == 'Success':
        #    projid = result['Project'][0].strip('\0')
        # Replaced by:
        # Read scan_header.dat file at /nas4/Tables/scanheader
        f = open('/common/Tables/scanheader/scan_header.dat')
        buf = f.read()
        f.close()
        msg = 'Success'
        if msg == 'Success':
            projid = buf[20:52].strip('\0')
    
            if projid == 'NormalObserving' or projid == 'FLARETEST':
                # Get new FEMATTN level, ensuring that it is not less than 0
                needs = np.round(needs).astype(int)
                for i in range(7):
                    if levs[i,0] == -1:
                        # This ROACH has all failed attempts to read from stateframe, 
                        # so use previous levels and skip any update
                        newlevs[i] = savlevs[i]
                    else:
                        newlevs[i] = np.clip(levs[i]+needs[i],0,None)
                if np.sum(levs[:6]) == 0:
                    # All of the levels are zero and we are on NormalObserving, so probably
                    # this is the start of a new scan. Hence, restore maximum of needs and 
                    # savlevs from previous time if available
                    newlevs = np.maximum(needs,copy(savlevs))
                savlevs = copy(newlevs)
                for i in range(13):
                    if agc[i] == 0:
                        # AGC is not active over the entire 30 s period, so adjust FEMATTN level toward target
                        nroach = i/2
                        nchan = (i % 2)*2
                        if needs[nroach,nchan] !=0 or needs[nroach,nchan+1] !=0:
                            # If either channel needs adjustment, send the command
                            femcmd = 'FEMATTN {} {} ANT{}'.format(newlevs[nroach,nchan],newlevs[nroach,nchan+1],i+1)
                            t = Time.now().iso
                            sys.stdout.write(t[11:19]+' '+femcmd+'\n')

                            adc2.send_cmds([femcmd], acc)
                            sleep(1)   # Leave a 1-s pause between sending of these commands
        
        outdict = {'time':t, 'fseqfile':fseqfile, 'projid':projid, 'levels':levs, 'needs':needs, 'stdev':stdall}
        fh2 = open(f2name,'ab')   # Open file for appending
        np.save(fh2,outdict)
        fh2.close()
        sys.stdout.write('Time now: '+t[:19]+'. Stops at: '+Time(mjdstop,format='mjd').iso[:19]+'\n')
        sys.stdout.flush()
            
    for i,t in enumerate(threads):
        t.stop=True
        
    #Wait for all threads to complete
    for t in threads:
        t.join()
        
    sys.stdout.write("Total Time= "+str(time()-tt))
    sys.stdout.close()
    exit()
    
def adc_plot(fname,tplot=None):

    import os
    fh = open(fname,'rb')

    #set up plot
    plt.ion()
    figure, ax = plt.subplots(4,7,figsize=(15, 8))
    plt.suptitle("ADC Standard Deviations", fontsize=20)
    polstr = [' X',' Y',' X',' Y']
    for i in range(7):
        for chan in range(4):
            ax[chan,i].text(25,75,'Ant '+str(i*2 + chan/2 + 1)+polstr[chan],horizontalalignment='center')
            ax[chan,i].text(25,65,'FEMATTN -',horizontalalignment='center',fontsize=9)
            ax[chan,i].plot([0,50],[32,32],'k--')
            ax[chan,i].text(0,34,'target',fontsize=9)
            ax[chan,i].set_ylim(0,100)
            if chan == 3:
                ax[chan,i].set_xlabel('Slot')
            if i == 0:
                ax[chan,i].set_ylabel('St. Dev.')
    figure.canvas.draw()
    figure.canvas.flush_events()

    # Find last record in the file
    while fh.tell() < os.fstat(fh.fileno()).st_size:
        out = np.load(fh)
        if tplot:
            t = Time(out.item()['time'])
            if (t - tplot) > 60./86400:
                break
    levs = out.item()['levels']
    needs = out.item()['needs']
    t = out.item()['time']
    projid = out.item()['projid']
    fseqfile = out.item()['fseqfile']
    stdev = out.item()['stdev']    
    strike = 0

    while(1):
        for i in range(7):
            for chan in range(4):
                ax[chan,i].plot(stdev[i,chan],'.')
                ax[chan,i].texts[1].set_text('FEMATTN {} {:+3.1f}'.format(levs[i,chan],needs[i,chan]))
                # Pop oldest line (lines[1]) from the plot (note lines[0] is the "target" line, which we want to keep)
                if len(ax[chan,i].lines) > 4: ax[chan,i].lines[1].remove()
        ax[0,6].set_title(t[11:19])
        ax[0,3].set_title(fseqfile)
        ax[0,0].set_title(projid)    

        # Save stdevs for this time
        # drawing updated values
        figure.canvas.draw()
        figure.canvas.flush_events()
        
        # If no time was given, we want a "real-time" display, otherwise
        # we are replaying past data
        if tplot:
            # Replay past data, sleeping only 1 s between plots
            sleep(1)
        else:
            # Real-time display
            # Make sure we are at the last record in the file before sleeping
            while fh.tell() < os.fstat(fh.fileno()).st_size:
                out = np.load(fh)
            # Sleep until the next record is expected (at 45 s past the minute)
            tnowsec = float(Time.now().iso[17:])
            sleep(60-(tnowsec - 47.0))  # Sleep until 47 s past the minute
        try:
            out = np.load(fh)
            strike = 0
            levs = out.item()['levels']
            needs = out.item()['needs']
            t = out.item()['time']
            projid = out.item()['projid']
            fseqfile = out.item()['fseqfile']
            stdev = out.item()['stdev']
        except IOError:
            # IOError might be a glitch, but 3 in a row means exit
            strike += 1
            if strike == 3: 
                # File is not growing, so end the loop
                break
    print 'File ended at time',t

def adc2master_table(fname=None, time=None, navg=10, attnval=None):
    ''' Reads the ADC measurements from a file in the /tmp folder and,
        given a fixed attenuation value analyzes the results and creates
        a master DCM table for calibrating the DCM attenuations.  Also
        asks the user to confirm before sending the table to SQL.
    '''
    # Read standard deviations and produce a DCM master table
    if fname is None:
        print 'Please include name of ADC data file (.npy file in /tmp)'
        return
    if time is None:
        print 'Please enter a time as Time() object.'
        return
    
    fh = open(fname,'rb')
#    buf = fh.read()
#    fh.close()
    if attnval is None:
        ans = raw_input('Enter the fixed DCM attn value used (integer dB, or q to quit): ')
        if ans.upper() == 'Q':
            return
        try:
            attnval = int(ans)
        except:
            print 'Response',ans,'did not translate to an integer dB.'
            return
        
    while (1):
        try:
            out = np.load(fh)
            if out.item()['time'] > time:
                break
        except IOError:
            print 'Reached end of file before time found.'
            fh.close()
            return
            
    stdev = [out.item()['stdev']]
    for i in range(navg-1):
        try:
            out = np.load(fh)
        except IOError:
            print 'Reached end of file before navg records read.'
            fh.close()
            return
        stdev.append(out.item()['stdev'])
    fh.close()
    data = np.array(stdev)
     
    data.shape = (navg,28,50)
    data[np.where(data < 0.00001)] = np.nan
    attn14 = np.log10(((np.nanmedian(data,0)/32.)**2))*10.
    new_table = np.zeros((52,30),int)
    for i in range(50):
        if i < 2:
            j = i
        else:
            j = i+2
        new_table[j,:28] = np.clip((np.round(attn14[:,i]+attnval)/2.).astype(int)*2,0,30)
    newtbl = []
    newtbl.append('#      Ant1  Ant2  Ant3  Ant4  Ant5  Ant6  Ant7  Ant8  Ant9 Ant10 Ant11 Ant12 Ant13 Ant14 Ant15')
    newtbl.append('#      X  Y  X  Y  X  Y  X  Y  X  Y  X  Y  X  Y  X  Y  X  Y  X  Y  X  Y  X  Y  X  Y  X  Y  X  Y')
    newtbl.append('#     ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- -----')
    for i in range(52):
        newtbl.append('{:2} : {:3}{:3}{:3}{:3}{:3}{:3}{:3}{:3}{:3}{:3}{:3}{:3}{:3}{:3}{:3}{:3}{:3}{:3}{:3}{:3}{:3}{:3}{:3}{:3}{:3}{:3}{:3}{:3}{:3}{:3}'.format(i+1,*new_table[i]))
    for line in newtbl:
        print line

    # ************ This block commented out due to loss of SQL **************
    ans = raw_input('Do you want to write this to the SQL database DCM master table? [y/n]?')
    if ans.upper() == 'Y':
        import cal_header
        print cal_header.dcm_master_table2sql(newtbl)
    # Replaced by:
    ans = raw_input('Do you want to save this to a DCM master table file? [y/n]?')
    if ans.upper() == 'Y':
        timestr = time.iso
        datstr = timestr[:19].replace('-','').replace(' ','_').replace(':','')
        filename = '/common/Tables/DCM_master/DCM_master_table_'+datstr+'.txt'
        f = open(filename, 'w')
        for line in newtbl:
            f.write(line+'\n')
        f.close()
        sleep(1)  # Make sure file has time to be closed.
        from ftplib import FTP
        f = open(filename, 'r')
        acc = FTP('acc.solar.pvt')
        acc.login('admin', 'observer')
        acc.cwd('parm')
        # Send DCM table lines to ACC
        print acc.storlines('STOR DCM_master_table.txt', f)
        f.close()
        print 'Successfully wrote DCM_master_table to ACC'


def another_pid(mypid):
    # check whether this routine is running in another instance of adc_plot.py (case-insensitive).
    # return PID if it is running, -1 if it is not
    import subprocess
    pidlist = subprocess.check_output(["pidof","python"]).split() # list of PIDs for all python processes
    for pid in pidlist:
        if mypid != pid:
            ps_out = subprocess.check_output(["ps","-lfp",pid])
            ind = ps_out.find('adc_plot.py') # position in string at which 'adc_plot.py' can be found (-1 if not found)
            if ind != -1:
                # This PID is running adc_plot.py
                return pid
    return False

if __name__ == "__main__":
    import os
    mypid = str(os.getpid())
    dup_pid = another_pid(mypid)
    if dup_pid:
        print "Another instance of adc_plot.py is already running at PID:", dup_pid
    else:
        adc_monitor()