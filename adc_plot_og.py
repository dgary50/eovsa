import roach as r
import struct
import numpy as np
import matplotlib.pyplot as plt
import threading
from time import sleep,time
from util import Time
import matplotlib.animation as animation

#set up threading
class AGC_Thread (threading.Thread):
    def __init__(self, threadID):
        threading.Thread.__init__(self)
        self.threadID = threadID                                        #set thread ID
        self.name = "roach" + str(threadID)                             #set the thread name
        self.sd = np.empty((4,50),float)                                #set up the standard deviation array
        self.iter_n = 0
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
        buf = []                                                            #adc buffer for single slot an channel
        self.iter_n += 1
        fh = open(roachname+'-'+str(self.iter_n)+'.dat','wb')
        udata = np.zeros((4,50,8192),float)                                 #numpy array of data to return in form of [chan,slot,data]
        done = np.zeros((4,50),bool)                                        #values in done array are true if a given slot and channel have been processed    
        lasttm = time()                                                     #last time is needed to ensure processing does not hang (gets stuck on slot 29 for some reason)
        rn = r.Roach(roachname)                                             #set up the roach
        
        #loop until all slots and channels processed or stop signal sent
        while not np.all(done) and not self.stop:
            f, s = np.modf(time())                                          #get fractional and integer part of the time
            tfrac = (np.floor(f*50)+1.0)/50.0                               #compute the fractional part of the time for the next adc read
            slot = int(tfrac*50) % 50                                       #get the slot number for the next read
            sleep(tfrac-f)                                                  #wait till next read
            chan, = np.where(done[:,slot] == False)                         # get list of unprocessed channels
            if (chan.size>0):                                               #if still channels to process
                #print "Processing Slot "+str(slot)+" Chan "+str(chan[0])    #display slot and channel (used for error checking) 
                
                #read in the adc values
                rn.fpga.write_int('swreg_snap_select',chan[0])
                rn.fpga.write_int('adc_data_adc_ctrl',0)
                # Calculate slot to be read, from time at instant before trigger is written
                t1 = time()
                rn.fpga.write_int('adc_data_adc_ctrl',7)
                t2 = time()
                f, s = np.modf(t1)
                slot_read = int(f*50)
                buf = rn.fpga.read('adc_data_adc_bram', 2048*4, 0)
                udata[chan[0],slot_read] = np.array(struct.unpack('>8192b', buf))
                
                done[chan[0],slot] = True                                   #set value of done for current slot and channel to true
                lasttm = time()                                             #update lasttm
                x = np.array([slot_read, t1, t2])
                fh.write(x)
            
            if time()-lasttm>4.0:                                           #if more than 4 seconds since last read, the exit the loop
                break
           
        #go through each slot to see if any left to process
        if not self.stop:
            for slot in range(50):
                chan, = np.where(done[:,slot] == False)                         #get channels to process for current slot
                if chan.size>0:                                                 #if channels left to be processed
                    for c in chan:                                              #go through each channel and get get adc values
                        f, s = np.modf(time())                                  #get fractional and integer part of the time 
                        dt = slot*0.02                                          #get fractional part of time for current slot 
                    
                        #wait for time to read slot
                        if f<dt:
                            sleep(dt-f)
                        else:
                            sleep(dt+1.0-f)
                    
                        #process the slot as above
                        #print "Processing Slot "+str(slot)+" Chan "+str(c) 
                        rn.fpga.write_int('swreg_snap_select',c)           
                        rn.fpga.write_int('adc_data_adc_ctrl',0)                 
                        t1 = time()
                        rn.fpga.write_int('adc_data_adc_ctrl',7)                 
                        t2 = time()
                        f, s = np.modf(t1)
                        slot_read = int(f*50)
                        buf = rn.fpga.read('adc_data_adc_bram', 2048*4, 0)       
                        udata[c,slot_read] = np.array(struct.unpack('>8192b', buf))
                        done[c,slot] = True       
                        x = np.array([slot_read, t1, t2])
                        fh.write(x)
         
        rn.fpga.stop()
        fh.close()
        
        return udata

tt=time()
stdev=np.empty((120,7,4,50),float)
stdev[:,:,:,:]=np.nan

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

#set up plot
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

# Open file for saving output stdevs
t = Time.now().iso[:19].replace('-','').replace(':','').replace(' ','')
fh = open('stdout_'+t+'.dat','wb')
start=time()
waittime = np.ceil(start/60)*60 - start + 45
print "Waiting for "+str(waittime)+" seconds"
sleep(waittime)

def animate(j):
    t = Time.now().iso
    print t
    
    # Loop over threads and get sd (4,50) for each
    stdall = np.zeros((7,4,50),float)
    for i in range(7):
        stdev = threads[i].sd
        stdall[i] = stdev
        for chan in range(4):
            ax[chan,i].plot(stdev[chan],'.')
            #Remove previous n-3 line (note lines[0] is the "target" line, which we want to keep)
            if len(ax[chan,i].lines) > 4: ax[chan,i].lines[1].remove()
            
    ax[0,6].set_title(t[11:19])
    
            
    # Save stdevs for this time
    fh.write(stdall)
            
# Create animation, new frame every 60 seconds (60,000 ms)
ani = animation.FuncAnimation(figure, animate, interval=60000)
plt.show(block=True)    #block=True is to prevent program proceeding until animation window closed 

#close the file
fh.close()

for i,t in enumerate(threads):
    t.stop=True
    
#Wait for all threads to complete
for t in threads:
    t.join()
    
print "Total Time= "+str(time()-tt)
print "*************************" 
