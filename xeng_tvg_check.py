from pcapfile import savefile
import struct
import numpy as np
from corr import sim
import sys


def reorder_glitch():
    '''Returns channel order for glitch in 16-ant correlator circa 2016 Feb 13
    '''
    kern = range(8)
    kern += (np.array(kern)+128).tolist()
    kern += (np.array(kern)+2048).tolist()
    kern1 = np.rollaxis(np.array(kern).reshape(4,8),1).reshape(32).tolist()
    ad = range(0,128,8)
    ad += (np.array(ad)+1024).tolist()
    ad += (np.array(ad)+256).tolist()
    ad += (np.array(ad)+512).tolist()
    chan = []
    for val in ad:
        chan += (np.array(kern1) + val).tolist()
    chan = np.array(chan)
    return chan

def get_order(nroach, nchans=4096):
    '''
    Returns the channel order after the troublesome
    firmware block when the firmware is configured
    for nroach ROACH boards.
    '''
    v = []
    for k in range(nchans/nroach/2/2):
        for j in range(2):
            for i in range(nroach*2):
                v += [(2*i + j)*nchans/nroach/4 + k]
    return v

def get_rel_order():
    '''
    Returns a mapping of channel label to real channel,
    based on the diagnosed misconfiguration of one of the
    firmware reorder blocks.
    '''
    v2 = get_order(2) #reordering based on firmware for 2 roaches (what we have)
    v8 = get_order(8) #reordering based on firmware for 8 roaches (what we want)
    
    #print len(v2), v2[:40]
    #print len(v8), v8[:80]
    
    rel_order = [0]*4096
    
    for i in range(4096):
        rel_order[v8[i]] = v2[i] #the channel we label as v8[i] is really v2[i]

    return rel_order

def decode_x_header(p):
    h = struct.unpack('<QLHHQLHHQQQQQQQ', p[0:8*11])
    mcnt = h[-1]
    chan = h[7]
    xeng = h[2]
    acc_num = h[1]
    return mcnt, chan, xeng, acc_num

def get_expected_output(n_ants=16):
    '''
    Generate the correlation matrix based on the
    test vectors specified in firmware.
    '''
    bl_order = sim.get_bl_order(n_ants) # Get a list of antenna pairs
    n_bls = len(bl_order)
    oput = np.zeros([4096, n_bls, 4], dtype=np.complex128) # n_chans * n_baselines * n_stokes
    iput = np.zeros([n_ants, 2, 4096], dtype=np.complex128) # n_antennas * n_polarizations * n_chans
    # First generate the inputs
    # y pol first.
    counter = np.arange(4096)
    im = counter&0xf
    re = (counter>>4)&0xf
    # convert to 2's complement binary
    re[re>7] -= 16
    im[im>7] -= 16
    # vectors are the same for all antennas
    for ant in range(n_ants):
        iput[ant, 1, :] = re + 1j*im

    # xpol
    im = (counter>>8) & 0xf
    # convert to 2's complement binary
    im[im>7] -= 16
    # xpol -- odd numbered antennas (i.e., second input of each ROACH) have 1 in real part
    for ant in range(n_ants):
        if ant % 2 == 0:
            iput[ant, 0, :] = 0 + 1j*im
        elif ant % 2 == 1:
            iput[ant, 0, :] = 1 + 1j*im

    # compute correlation matrix
    for bn,bl in enumerate(bl_order):
        oput[:, bn, 0] = iput[bl[0],0] * np.conj(iput[bl[1],0]) #XX
        oput[:, bn, 1] = iput[bl[0],1] * np.conj(iput[bl[1],1]) #YY
        oput[:, bn, 2] = iput[bl[0],0] * np.conj(iput[bl[1],1]) #XY
        oput[:, bn, 3] = iput[bl[0],1] * np.conj(iput[bl[1],0]) #YX

    return oput

from optparse import OptionParser
parser = OptionParser(usage='%prog [options] pcap_file1 pcap_file2 ...')
parser.add_option("-r", "--reorder", dest="reorder", action="store_true", default=False,
                  help="Use this flag to reorder frequency channels based on a known firmware issue. This issue was fixed in firmware version eovsa_corr_2016_Feb_20_1600.bof. This option should be unnecessary for pcap files generated with this, or later, firmware.")
parser.add_option("-p", "--npackets", dest="npackets", type="int", default=100000,
                  help="Maximum number of packets to read (per pcap files). Default:100000")
(opts, args) = parser.parse_args()

if opts.reorder:
    #rel_order[x] = y -> the channel labelled x really y
    # i.e. real_channel = rel_order[channel_label]
    rel_order = get_rel_order()
else:
    rel_order = range(4096)

bl_order = sim.get_bl_order(16)

limit = opts.npackets

hdr_len = 42+88
# we only capture 1000 bytes per packet --
# header is 42+88 bytes => 870 bytes data
# => 108 complex values
# => 27 dual-pol baselines (+ some leftover)
n_bls = 27
#n_bls = 16*17/2 * 4

n_windows = 100
packet_cnts = np.zeros(n_windows)
time_slots = np.ones(n_windows) * -1
spectra = np.ones([n_windows, 4096, n_bls*4], dtype=complex)*2**31

for fn in args:
    with open(fn, 'r') as fh:
        print 'loading %s'%fn
        capfile = savefile.load_savefile(fh)
        for pn, p in enumerate(capfile.packets[::1]):
            if pn % 10000 == 0:
                print 'Read %d packets'%pn
            # check this is an X packet
            if p.packet_len == 4482:
                mcnt, chan, xeng, acc_num = decode_x_header(p.raw()[42:])
                t = mcnt >> 12 #the lower bits of mcnt are a channel ID (and should be fixed per roach)
                t = acc_num
                window = t%n_windows
                data = struct.unpack('>%dl'%(n_bls*2*4), p.raw()[hdr_len:hdr_len + (4*n_bls*8)])
                data_c = np.array(data[::2]) + 1j*np.array(data[1::2])
                data_c /= 1792.
                #print pn, acc_num, mcnt, xeng, mcnt>>12, '%4d'%(mcnt&(2**12-1)), 'chan: %4d'%chan, data_c[0]
                time_slots[window] = t
                packet_cnts[window] += 1
                last_t = t
                spectra[window, rel_order[chan], :] = data_c
            else:
                pass
            if pn == limit:
                break

print ''
for i in range(n_windows):
    if time_slots[i] != -1:
        print 'time %d had %d packets'%(time_slots[i], packet_cnts[i])


badcnt = 0
totcnt = 0
goodcnt = 0
print 'Checking outputs against TVG'
expected = get_expected_output()

# iterate over time slots
tn = 0
for t in time_slots:
    # skip if this timeslot isn't valid
    if t == -1:
        continue
    else:
        tn += 1
    print 'Checking time slot %d (time %d).'%(tn,t),
    vldpkts = 0
    # iterate over channels (= packets)
    for cn,cspec in enumerate(spectra[t%n_windows]):
        # iterate over baselines
        for bn,blspec in enumerate(cspec):
            # check if a packet was received (else the value will be that at initialisation
            if blspec != 2**31:
                totcnt += 1
                if bn == 0:
                    vldpkts += 1
                if blspec != expected[cn, bn//4, bn%4]:
                    badcnt += 1
                    print i, t, bn//4, bn%4, cn, cspec, expected[cn, bn//4, bn%4]
                else:
                    goodcnt += 1
    print 'Valid packets in this time slice: %d'%vldpkts

print '%d values checked -- %d good, %d bad'%(totcnt, goodcnt, badcnt)
