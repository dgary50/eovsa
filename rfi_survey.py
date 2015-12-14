#
# Routines to set up system for RFI survey, capture and display pcap
# packet capture files
#
#   2014-Dec-20  DG
#     First assembled routines from elsewhere.
#
import pcapture as p

def acc_tune(band):
    import time, socket, stateframe
    if type(band) is int:
        fsqfile = 'BAND'+str(band)+'.FSQ'
    elif type(band) is str:
        if band.lower() == 'solar.fsq' or band.lower() == 'pcal.fsq':
            fsqfile = band.lower()
    else:
        print 'Error: Unknown band',band
        return
        
    try:
        accini = stateframe.rd_ACCfile()
    except:
        print 'Error: Could not access ACC.'
        return
    cmds = ['FSEQ-OFF','FSEQ-INIT','WAIT','FSEQ-FILE '+fsqfile.lower(), 'FSEQ-ON']
    for cmd in cmds:
        print 'Command:',cmd
        s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        try:
            s.connect((accini['host'],accini['scdport']))
            s.send(cmd)
            time.sleep(0.01)
            s.close()
        except:
            print 'Error: Could not send command',cmd,' to ACC.'
    return

def set_roach_attn():
    '''Auto-set the ROACH attenuations, plus set fftshift to 31'''
    import roach as r
    r1 = r.Roach('roach1')
    r1.set_attn(update=True)
    pktfft = r1.fpga.read_int('swreg_pkt_fft')
    r1.fpga.write_int('swreg_pkt_fft',(pktfft & 0xff0000)+31)
    r2 = r.Roach('roach2')
    r2.set_attn(update=True)
    pktfft = r2.fpga.read_int('swreg_pkt_fft')
    r2.fpga.write_int('swreg_pkt_fft',(pktfft & 0xff0000)+31)
    r1.fpga.stop()
    r2.fpga.stop()

def grab_band(band):
    '''Capture 3 seconds of data on each of interfaces eth2 and eth3, saving the
       results in files in /home/user/Python'''
    outfile = '/home/user/Python/band'+str(band)+'_2.pcap'
    command = '/usr/sbin/tcpdump -i eth2 -c '+str(19200*3)+' -w '+outfile+' -s 2000'
    p.sendcmd(command)
    outfile = '/home/user/Python/band'+str(band)+'_3.pcap'
    command = 'tcpdump -i eth3 -c '+str(19200*3)+' -w '+outfile+' -s 2000'
    p.sendcmd(command)

def do_rfi_survey():
    '''Go through each band, tuning, setting attenuation, and capturing packets.'''
    import time
    for band in range(1,35):
        acc_tune(band)
        set_roach_attn()
        time.sleep(1)
        grab_band(band)

def show_rfi_sk(band):
    '''Display the result of an SK measurement for the given band
       (an integer between 1 and 34).'''
    import matplotlib.pyplot as plt
    import numpy as np
    pcapfile = '/home/user/Python/band'+str(band)+'_2.pcap'
    out1 = p.rd_spec(pcapfile)    
    out1.shape = (200,4096,8)
    good = np.where(out1[:,1000,0] != 0)[0]
    n = len(good)
    out1 = out1[good,:,:]
    sk = (1792+1)/(1792.-1) * (1792*out1[:,:,4:8]/(out1[:,:,0:4]**2) - 1)
    plt.figure()
    title = ['Ant1 X','Ant1 Y', 'Ant2 X','Ant2 Y']
    for i in range(4):
        plt.subplot(4,1,i+1)
        plt.imshow(np.logical_or(sk[:,:,i]>1.15,sk[:,:,i]<0.85),aspect='auto',interpolation='nearest')
        plt.title(title[i])
    plt.set_cmap('gray')
    pcapfile = '/home/user/Python/band'+str(band)+'_3.pcap'
    out2 = p.rd_spec(pcapfile)
    out2.shape = (200,4096,8)
    good = np.where(out2[:,1000,0] != 0)[0]
    n = len(good)
    out2 = out2[good,:,:]
    sk = (1792+1)/(1792.-1) * (1792*out2[:,:,4:8]/(out2[:,:,0:4]**2) - 1)
    plt.figure()
    title = ['Ant3 X','Ant3 Y', 'Ant4 X','Ant4 Y']
    for i in range(4):
        plt.subplot(4,1,i+1)
        plt.imshow(np.logical_or(sk[:,:,i]>1.15,sk[:,:,i]<0.85),aspect='auto',interpolation='nearest')
        plt.title(title[i])
    return out1,out2

def qlook(file):
    ''' Routine to show all four channels (2 ants, 2 poln) for
        a pcapture file with a name of pattern 'bandnn_e.pcap'
        where nn is band number (can be 1 digit) and e is 
        ethernet interface (eth) number (2 or 3)
    '''
    import matplotlib.pyplot as plt
    out = p.rd_spec(file).reshape(200,4096,8)
    f, ax = plt.subplots(4,1,sharex=True)
    ax[0].imshow(out[:,:,0])
    ax[1].imshow(out[:,:,1])
    ax[2].imshow(out[:,:,2])
    ax[3].imshow(out[:,:,3])
    band = file.split('band')[1].split('_')[0]
    if file.find('_2') != -1: 
        ax[0].set_title('Ants 1 & 2, band '+band)
    else:
        ax[0].set_title('Ants 3 & 4, band '+band)
    return out