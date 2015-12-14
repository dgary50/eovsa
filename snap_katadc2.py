#!/usr/bin/env python

import corr, os, struct, time
import matplotlib
matplotlib.use('TkAgg')

import matplotlib.pyplot as plt
import numpy as np
import sys
import datetime

'''Arguments: Roach, ADC board, Channel, RecordData
'''

# Get command line arguments
if len(sys.argv) < 2:
    # Called without arguments, so set defaults
    roach_ip = 'roach1.solar.pvt'
    board = 0
    channel = 'both'
    recdata = None
elif len(sys.argv) == 2:
    roach_ip = sys.argv[1]
    board = 0
    channel = 'both'
    recdata = None
elif len(sys.argv) == 3:
    roach_ip = sys.argv[1]
    board = sys.argv[2]
    channel = 'both'
    recdata = None
elif len(sys.argv) == 4:
    roach_ip = sys.argv[1]
    board = sys.argv[2]
    channel = sys.argv[3]   # Options are 'I', 'Q', or 'both'
    recdata = None
else:
    roach_ip = sys.argv[1]
    board = sys.argv[2]
    channel = sys.argv[3]   # Options are 'I', 'Q', or 'both'
    recdata = sys.argv[4]   # Record if equal to 'rec'
    
if recdata  == 'rec' or recdata == 'Rec':
    rec = True
else:
    rec = False
  
roach = corr.katcp_wrapper.FpgaClient(roach_ip)
time.sleep(1)
devs = roach.listdev()


def trigger():
    ''' Trigger the snapshot of data by writing 0 to _ctrl registers, then 7.
    '''
    for dev in devs:
        if dev.find('_ctrl') != -1:
            roach.write_int(dev,0)

    #time.sleep(0.001)

    for dev in devs:
        if dev.find('_ctrl') != -1:
            roach.write_int(dev,7)

    #time.sleep(0.001)


def grab_spec(roach,n=100,adc=0,chan=0,):
    ''' Grab n (default 100) snapshots of data and create summed spectrum
    '''
    spec = np.zeros(4096,'float')
    bram = 'data_adc'+str(adc)+'_ch'+str(chan)+'_bram'

    for i in range(n):
        trigger()
        data = roach.read(bram,2048*4,0)
        blah = np.array(struct.unpack('>8192b',data))
        spec += abs(np.fft.fft(blah))[0:4096]
        if rec:
            df.write(data)
    # Return last snapshot and accumulated spectrum
    return blah, spec

def plot_callback(channel, board, roach):
    ''' Called every few seconds by canvas manager
    '''
    #print channel, board, roach
    global df, rec,  roach_ip

    if channel != 'I' and channel != 'Q' and channel != 'both':
        print 'Channel input not recognized, set to both'
        channel = 'both'
    if int(board) != 1 and int(board) != 0:
        print board, 'Board input not recognized, set to 0'
        board = 0
    
    try:
        if channel == 'I' or channel == 'both':
            if rec:
                df = open('/archive/eovsa/Katadc_snapshot_data/'+roach_ip+'_'+board+'_'+datetime.datetime.now().strftime('%Y_%m_%d_%H:%M:%S')+'I.dat', 'w')
                
            data, spec = grab_spec(roach,100,board,0)
            if rec:
                df.close()
        if channel == 'Q' or channel == 'both':
            if rec:
                df = open('/archive/eovsa/Katadc_snapshot_data/'+roach_ip+'_'+board+'_'+datetime.datetime.now().strftime('%Y_%m_%d_%H:%M:%S')+'Q.dat', 'w')
            data2, spec2 = grab_spec(roach,100,board,1)
            if rec:
                df.close()
        

        subplots[0].cla()
        subplots[0].set_xticks(range(-130, 131, 20))
        
        if channel == 'I' or channel == 'both':
            histData, bins, patches = subplots[0].hist(data, bins = 256, range = (-128,128))
        
        if rec:
            yesno = 'Recording'
        else:
            yesno = ''
        
        subplots[0].set_title('Roach: '+roach_ip+', ADC Board #: '+board+', Channel: '+channel+'        '+yesno)
        subplots[0].set_ylabel('Counts')
        subplots[0].set_xlabel('ADC sample bins.')

        if channel == 'Q' or channel == 'both':
            histData, bins, patches = subplots[0].hist(data2, bins=256, range=(-128,128))
            
        plt.ylim(ymax = (max(histData) * 1.05))
        subplots[1].cla()
        
        if channel == 'I' or channel == 'both':
            subplots[1].plot(data,'b')
        if channel == 'Q' or channel == 'both':
            subplots[1].plot(data2,'g')

        subplots[2].cla()
        
        if channel == 'I' or channel == 'both':
            subplots[2].plot(spec,'b')
        if channel == 'Q' or channel == 'both':
            subplots[2].plot(spec2,'g')
            
        plt.yscale('log')
        plt.ylim(ymin=1.e4, ymax=1.e6)
        fig.canvas.draw()
        fig.canvas.manager.window.after(100, plot_callback, channel, board, roach)

        
    except KeyboardInterrupt:
        
        if rec:
            startstop = 'stop'
        else:
            startstop = 'start'
        
        prompt = raw_input('\nEnter "S" to '+startstop+' recording, "X" to change channel, "quit" to exit program: ')
        if prompt == 'S':
            rec = not rec
            if rec:
                print 'Started recording'
            else:
                print 'Stopped recording'            
        elif prompt == 'quit':
            sys.exit()
        elif prompt == 'X':
       
            #Change Channel
            channel = raw_input('\nChannel I or Q (Enter "both" for both, "X" to change adc board, "quit" to exit program): ')
            if channel == 'quit':
                sys.exit()
            elif channel == 'X':
            
                #Change Board
                board = raw_input('Board 0 or 1 (Enter "X" to change Roach IP address, "quit" to exit program): ')
                if board == 'quit':
                    sys.exit()
                elif board == 'X':
                
                    #Change Roach  
                    roach_ip = raw_input('Enter Roach IP address ("quit" to exit program): ')
                    if roach_ip == 'quit':
                        sys.exit()
                    roach = corr.katcp_wrapper.FpgaClient(roach_ip)
                    board = raw_input('Board 0 or 1:')
                    channel = raw_input('Channel I or Q (Enter "both" for both): ')
                    fig.canvas.manager.window.after(100, plot_callback, channel, board, roach)
                    
                elif int(board) == 0 or int(board) == 1:
                    channel = raw_input('Channel I or Q (Enter "both" for both): ')
                
                fig.canvas.manager.window.after(100, plot_callback, channel, board, roach) 
               
            fig.canvas.manager.window.after(100, plot_callback, channel, board, roach)
        
        else:
            print 'input not recognized'
        fig.canvas.manager.window.after(100, plot_callback, channel, board, roach)
                    
if __name__ == '__main__':

    # create the subplots
    fig = plt.figure()
    subplots = []

    n_subplots=3

    for p in range(n_subplots):
        subPlot = fig.add_subplot(n_subplots, 1, p + 1)
        subplots.append(subPlot)

    # start the process
    print 'Starting plots...'
    fig.subplots_adjust(hspace=0.8)
    fig.canvas.manager.window.after(100, plot_callback, channel, board, roach)
    plt.show()

