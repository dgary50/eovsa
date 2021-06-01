import roach as r
import struct
import numpy as np
import matplotlib.pyplot as plt
import threading
from time import sleep,time
from util import Time

#set up threading
class AGC_Thread (threading.Thread):
    def __init__(self, threadID):
        threading.Thread.__init__(self)
        self.threadID = threadID                                        #set thread ID
        self.name = "roach" + str(threadID)                             #set the thread name
        self.sd = np.empty((4,50),float)                                #set up the standard deviation array
        self.stop = False
    
    def run(self):
        threadLock.acquire(0)                                           #set threadlock
        while not self.stop:
            start = time()                                              #used for computing execution time
            f, s = np.modf(start)
            waittime = (np.floor(s/60)+1)*60 - start
            sleep(waittime)
            start=time()
            self.sd = np.std(self.grab_all(),axis=2)                    #calculate standard deviation from the roach
            print self.name+" Execution Time= "+str(time()-start)           #display the execution time
            
    def grab_all(self):
        roachname = self.name
        buf = ''                                                            #adc buffer for single slot an channel
        slots = []
        rn = r.Roach(roachname)

        # Grab 10 slots in one channel in one second (every fifth slot) starting at the given slot.
        # This will be called five times starting with slot 0, 1, 2, 3, and 4.
        def grab10(rn,slot,chan):
            ''' Returns the appended raw byte buffer for all 10 slots, and the list of slots
                actually measured.
            '''
            f, s = np.modf(time())
            sleep(1-f+0.02*slot-0.01)
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

        for chan in range(4):
            for slot in range(5):
               buf1, s = grab10(rn, slot, chan)
               buf += buf1
               slots.append(s)
        nbytes = 8192*4*50
        udata = np.array(struct.unpack('>'+str(nbytes)+'b',buf))
        udata.shape = (4,50,8192)
        outdata = np.zeros((4,50,8192),float)                                 #numpy array of data to return in form of [chan,slot,data]
        # The data are a bit scrambled in udata, and some bands may be missing (and others measured twice)
        # although hopefully that will be rare.  This should unscramble them.
        slots = np.array(slots)
        slots.shape = (4,50)
        def find_missing(x):
            bad = []
            ptr = 0
            for i in range(50):
                try:
                    if x[ptr] != i:
                        bad.append(i)
                    else:
                        ptr += 1
                except IndexError:
                    bad.append(i)
            return bad

        sleep(int(roachname[-1])*0.02)
        for c in range(4):
            slotarray, idx = np.unique(slots[c],return_index=True)
            print roachname, c, '   Missing:',50-len(idx)
#            if len(idx) != 50:
#                print roachname,c,find_missing(slotarray)
            outdata[c,slotarray] = udata[c,idx]
        rn.fpga.stop()
        
        return outdata

#set up plot
plt.ion()
figure, ax = plt.subplots(4,7,figsize=(15, 8))
plt.suptitle("ADC Standard Deviations", fontsize=20)
polstr = [' X',' Y',' X',' Y']
for i in range(7):
    for chan in range(4):
        ax[chan,i].text(25,75,'Ant '+str(i*2 + chan/2 + 1)+polstr[chan],horizontalalignment='center')
        ax[chan,i].plot([0,50],[32,32],'k--')
        ax[chan,i].text(0,34,'target',fontsize=9)
        ax[chan,i].set_ylim(0,100)
        if chan == 3:
            ax[chan,i].set_xlabel('Band')
        if i == 0:
            ax[chan,i].set_ylabel('St. Dev.')
figure.canvas.draw()
figure.canvas.flush_events()

# Open file for saving output stdevs
t = Time.now().iso[:16].replace('-','').replace(':','').replace(' ','_')
fname = '/tmp/stdout'+t+'.dat'
fh = open(fname,'wb')

tt=time()

started=False

#set threading
threadLock = threading.Lock()                                           
threads = []

#set the threads
for t in range(1,8):            
    threads.append(AGC_Thread(t))

#Start new Threads 
for t in threads:
    t.start()

# Give threads some time to work
sleep(45)

for l in range(30):
    t0 = time()
    # Read once per minute on the 45 s
    waittime = 45 - t0 % 60
    if waittime < 0: waittime += 60
    print "waittime: "+str(waittime)

    sleep(waittime)
    t = Time.now().iso
    print t
    
#    stdev[0:-1]=stdev[1:]
#    for i,t in enumerate(threads):
#        stdev[-1,i]=t.sd
    
    # Loop over threads and get sd (4,50) for each
    stdall = np.zeros((7,4,50),float)
    for i in range(7):
        stdev = threads[i].sd
        stdall[i] = stdev
        for chan in range(4):
            ax[chan,i].plot(stdev[chan],'.')
            # Pop oldest line (lines[1]) from the plot (note lines[0] is the "target" line, which we want to keep)
            if len(ax[chan,i].lines) > 4: ax[chan,i].lines[1].remove()
    ax[0,6].set_title(t[11:19])
            
    # Save stdevs for this time
    fh.write(stdall)
    # drawing updated values
    figure.canvas.draw()
    
    # This will run the GUI event
    # loop until all UI events
    # currently waiting have been processed
    figure.canvas.flush_events()
        
# All done saving stdevs, so close file
fh.close()
print 'Output saved in file',fname,'\n'

for i,t in enumerate(threads):
    t.stop=True
    
#Wait for all threads to complete
for t in threads:
    t.join()
    
plt.close(figure)
print "Total Time= "+str(time()-tt)
print "*************************" 
    
ans = 'Y'
ans = raw_input('Do you want to form a DCM master table? (say yes only if using a fixed DCM attn) [y/n]?')
if ans.upper() == 'Y':
    ans = raw_input('Enter the fixed DCM attn value used (integer dB, or Q to quit): ')
    if ans.upper() == 'Q':
        exit()
    try:
        attnval = int(ans)
    except:
        print 'Response',ans,'did not translate to an integer dB.  Try again or Q to quit.'
        ans = raw_input('Enter the fixed DCM attn value used (integer dB, or Q to quit): ')
        if ans.upper() == 'Q':
            exit()
        try:
            attnval = int(ans)
        except:
            print 'Response',ans,'still did not translate to an integer dB.  Cannot continue.'
            exit()
else:
    exit()
    
# Read standard deviations and produce a DCM master table
fh = open(fname,'rb')
buf = fh.read()
fh.close()
data = array(struct.unpack('14000d',buf))
data.shape = (10,28,50)
data[where(data < 0.00001)] = nan
attn14 = np.log10(((np.nanmedian(data,0)/32.)**2))*10.
new_table = zeros((52,30),int)
for i in range(50):
    if i < 2:
        j = i
    else:
        j = i+2
    new_table[j,:28] = np.clip((np.round(attn14[:,i]+attnval)/2.).astype(int)*2,0,30)
newtbl = []
newtbl = []
newtbl.append('#      Ant1  Ant2  Ant3  Ant4  Ant5  Ant6  Ant7  Ant8  Ant9 Ant10 Ant11 Ant12 Ant13 Ant14 Ant15')
newtbl.append('#      X  Y  X  Y  X  Y  X  Y  X  Y  X  Y  X  Y  X  Y  X  Y  X  Y  X  Y  X  Y  X  Y  X  Y  X  Y')
newtbl.append('#     ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- -----')
for i in range(52):
    newtbl.append('{:2} : {:3}{:3}{:3}{:3}{:3}{:3}{:3}{:3}{:3}{:3}{:3}{:3}{:3}{:3}{:3}{:3}{:3}{:3}{:3}{:3}{:3}{:3}{:3}{:3}{:3}{:3}{:3}{:3}{:3}{:3}'.format(i+1,*new_table[i]))
for line in newtbl:
    print line

ans = 'Y'
ans = raw_input('Do you want to write this to the SQL database DCM master table? [y/n]?')
if ans.upper() == 'Y':
    import cal_header
    print cal_header.dcm_master_table2sql(newtbl)
