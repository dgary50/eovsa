# pcapture.py
#    Routines to enable packet capture on the 10-Gbe interfaces on the dpp,
#    and also to read the "pcap" files and return or display the data in
#    the packets.
#
# History:
#   2014-Dec-05  DG
#     The routines are now matured to usefulness, so this is the start of
#     the official history record.  TODO: Read stateframe information from 
#     the SQL server to display on the pcapture.png file indicating LO sweep 
#     status.
#   2014-Dec-13  DG
#     Updated figure plot to use better labeling.
#   2015-Jan-21  DG
#     Updated to work with 8 antennas (2 pair of ROACHes)
#   2015-Jan-24  DG
#     Added code to allow any pair of ROACHes to be plotted with capture_fig()
#   2015-Jan-30  DG
#     Eliminated inversion of plot of overflow points, now that plots are 
#     "right-side-up".
#   2015-Feb-24  DG
#     Oops.  Change to ROACH setup on Feb. 18 changed the interface the packets
#     were sent to, which broke this routine.  Now fixed.
#   2015-Mar-03  DG
#     The dppcapture.png file used to only include 4 antennas--now fixed to plot
#     all 8 currently in use.
#   2015-Mar-15  DG
#     Fix a bug in capture_fig that accessed the same ADC board twice instead of
#     switching.
#   2015-May-15  DG
#     Improve attenuation printout to allow list of backend attenuations, and
#     improve handling of negative attenuations.
#   2015-May-29  DG
#      Converted from using datime() to using Time() based on astropy
#   2015-Jun-13  DG
#      Changed output of rd_spec() to concatenate first two dimensions into one
#   2015-Jun-16  DG
#      Added new routine capture(), which captures on all interfaces and
#      returns a complete set of data for either P or X (not P^2)
#   2015-Aug-23  DG
#      Changed capture() output to be in "standard" order (nant or nbl,npol,nf,nt)
#   2015-Aug-28  DG
#      I realized that ovfl values in plot_fig() were plotted wrong--each ROACH
#      only has one ovfl value per band, i.e. the two ADC boards produce only one
#      combined ovfl value.
#   2016-Jan-18  DG
#      Added code to handle data from 16-antenna "production" correlator.
#   2016-Feb-26 BC
#      Updated routine capture() to reflect the new 16-antenna correlator design. Changed the
#      output to a dictionary {p, p^2} for ptype 'P', or {auto-corr, x-corr} for ptype 'X'
#      Slightly modified rd_spec() to fix the time dimension
#   2016-Mar-01  DG
#      Added rd_jspec() routine to read a 1-s packet dump file from Jim McT
#   2016-Mar-01 BC
#      Added get_spec() routine as a wrapper of capture() and rd_jspec()
#      Added bl_mapper() to find out baseline indices for selected antenna in the x-corr output
#      Added plot_auto_corr() to show auto correlation results
#      Added plot_x_corr() to show cross correlation results
#   2016-Mar-06  DG
#      Added bl_list() routine to return bl2ord lookup table
#      Added summary_plot() routine to create an overview plot of cross-correlations
#   2016-Mar-16  DG
#      Added get_bl_order() from the corr.sim module, so that it is not necessary to
#      import the entire, complex corr module!
#   2016-Apr-01  DG
#      Added ant_str2list()

import numpy as np
import pdb

def ant_str2list(ant_str):
    ant_list = []
    try:
        grps = ant_str.split()
        for grp in grps:
            antrange = grp[3:].split('-')
            if len(antrange) == 1:
                if antrange != '':
                    ant_list.append(int(antrange[0])-1)
            elif len(antrange) == 2:
                ant_list += range(int(antrange[0])-1,int(antrange[1]))
    except:
        print 'Error: cannot interpret ant_str',ant_str
        return None
    return np.array(ant_list)

def sudo(command, password=None, prompt="Enter password "):
    ''' Issues <command> string as sudo, and prompts for password
        if not supplied.
        
        use with command = 'tcpdump -i eth3 -c 19200 -w eth3.pcap -s 2000'
        where 19200 gives 1 s of packets from interface eth3.  Or use eth2.
        Here 'eth3.pcap' is the filename that will contain the packets.
    '''
    import pexpect

    if not password:
        import getpass
        password = getpass.getpass(prompt)

    # See if the command is consistent with tcpdump -c nnnn
    # where nnnn is a number of packets to capture.  19200 packets take 1 s
    # so set timeout to 10 s more than that the nominal time required.
    timeout = 10
    if command.find('tcpdump') != -1:
        tok = command.split('-c ')
        if len(tok) == 2:
            # Yep, there is a -c in the line
            try:
                timeout += int(tok[1].split()[0])/19200
            except:
                pass
    child = pexpect.spawn('sudo '+command)
    child.expect(['ssword', pexpect.EOF])
    child.sendline(password)
    child.expect(pexpect.EOF, timeout=timeout)
    # is this necessary?
    child.close()

def sendcmd(command):
    ''' Issues <command> string
        Use with command = 'tcpdump -i eth3 -c 19200 -w eth3.pcap -s 2000'
        where 19200 gives 1 s of packets from interface eth3.  Or use eth2.
        Here 'eth3.pcap' is the filename that will contain the packets.
        Output is list of strings output from command.  If tcpdump command
        is sent, typical output is:
        ['tcpdump: listening on eth3, link-type EN10MB (Ethernet), capture size 2000 bytes',
         '19200 packets captured',
         '19204 packets received by filter',
         '0 packets dropped by kernel',
         '']
    '''
    import subprocess

    out = subprocess.check_output(command.split(),stderr=subprocess.STDOUT).split('\n')
    return out


def list_header(filename,ptype='P',boardID=None,verbose=False):
    ''' Reads the header of each record of a file and returns a list of header lines, one
        line per spectrum.  Returns a vector of header lines.  If boardID is not None,
        only header lines from the provided boardID are returned (if any).
    '''
    import dpkt, struct

    f = open(filename,'rb')
    pcap = dpkt.pcap.Reader(f)
    hdr = '<HHHHIHHHHHHHHHHHHHH4I4Q'
    khdr = ['HeaderLength','PacketNum','FFTShift','AccumLength','GlobalAccumNum',
            'BoardID','AccumNum','DataType','PolType','Ai','Aj','ADCOverflow',
            'QuantClipNum','NSubbands','iFreq','Delay0','Delay1','Delay2','Delay3',
            'PX0','PY0','PX1','PY1','P2X0','P2Y0','P2X1','P2Y1']
    lines = []
    # Start reading records
    out = pcap.readpkts()
    npkt = len(out)
    for j in range(npkt):
        if verbose: print 'working on packet',j,'\r',
        buf = out[j][1]
        if len(buf) == 898:
            # This is a P packet
            if ptype == 'P':
                header = struct.unpack(hdr,buf[-856:88-856])
                h = dict(zip(khdr,header))
                l = h['HeaderLength']
                if l == 88:
                    if boardID != None and h['BoardID'] != boardID:
                        pass
                    else:
                        n = h['PacketNum']
                        a = h['AccumNum']
                        if n == 0:
                            lines.append('{:4d}'.format(h['AccumNum'])
                                +' {:11d}'.format(h['GlobalAccumNum'])
                                +' {:4d}'.format(h['FFTShift'])
                                +' {:2d}'.format(h['BoardID'])
                                +' {:5d}'.format(h['AccumLength'])
                                +' {:5d}'.format(h['ADCOverflow'])
                                +'  {:5d} {:5d}  {:5d} {:5d}'.format(h['PX0'],h['PY0'],h['PX1'],h['PY1'])
                                +'  {:4d} {:4d} {:4d} {:4d}'.format(h['Delay0'],h['Delay1'],h['Delay2'],h['Delay3']))
                            if a == 0:
                                # At every 0 Acc#, pop the last line, add a "header" and then add line back
                                line = lines.pop()
                                lines.append('Acc# Global-Acc# FFTs ID --M--  OvFl    P1x, P1y    P2x, P2y    ---Delays [nsec]---')
                                lines.append(line)
        elif len(buf) > 898:
            # This is an X packet
            if ptype == 'X':
                header = struct.unpack(hdr,buf[42:88+42])
                h = dict(zip(khdr,header))
                l = h['HeaderLength']
                if l == 88:
                    if boardID != None and h['BoardID'] != boardID:
                        pass
                    else:
                        n = h['PacketNum']
                        a = h['AccumNum']
                        if n == 0:
                            lines.append('{:4d}'.format(h['AccumNum'])
                                +' {:11d}'.format(h['GlobalAccumNum'])
                                +' {:4d}'.format(h['FFTShift'])
                                +' {:2d}'.format(h['BoardID'])
                                +' {:5d}'.format(h['AccumLength'])
                                +' {:5d}'.format(h['ADCOverflow'])
                                +'  {:5d} {:5d}  {:5d} {:5d}'.format(h['PX0'],h['PY0'],h['PX1'],h['PY1'])
                                +'  {:4d} {:4d} {:4d} {:4d}'.format(h['Delay0'],h['Delay1'],h['Delay2'],h['Delay3']))
                            if a == 0:
                                # At every 0 Acc#, pop the last line, add a "header" and then add line back
                                line = lines.pop()
                                lines.append('Acc# Global-Acc# FFTs ID --M--  OvFl    P1x, P1y    P2x, P2y    ---Delays [nsec]---')
                                lines.append(line)
    return lines

def rd_spec(filename,ptype='P',boardID=0,nboards=2,verbose=False):
    ''' Read all spectra in a packet capture file according to ptype and board ID.
        If ptype = 'P', read P and P^2 packets only, and returns a float
        array of size [nsec, 50, 4096, 8], where nsec is one more than
        the number of seconds in the packet capture file, and nbds is the
        number of boards currently active.  The 4 columns for each board
        are p1x, p1y, p2x, p2y, P1x, P1y, P2x, P2y (lower case p = P,
        upper case P = P^2).
        If ptype = 'X', read X packets for the channels captured on the
        interface.  This will be for all antennas, but only a subset
        of channels.  However, the array returned is of size
        [nsec,50,4096,n(n-1)], where nsec is as above, and n is the
        number of antenans in the array.  Of the 4096 channels, only
        those in the packets are non-zero.  The outputs from multiple
        interfaces can thus be added to fill in more channels.
        
        Now that we have the production correlator, I have added the
        proto=True keyword, so that if it is set to False the X data
        will be interpreted as from the production correlator.
    '''
    import dpkt, struct

    def pspectra(power):
        p = np.array(power)
        # idx is array([0,1, 16,17, 32,33, 48,49, 64,65, 80,81, 96,97, 112,113])
        idx = np.array(zip(np.arange(0,128,16),np.arange(0,128,16)+1)).flatten()
        # Distribute P, P2 into their spectra
        p1x = p[idx]/(2.**14)
        p1y = p[idx+2]/(2.**14)
        p2x = p[idx+4]/(2.**14)
        p2y = p[idx+6]/(2.**14)
        P1x = p[idx+8]/(2.**44)
        P1y = p[idx+10]/(2.**44)
        P2x = p[idx+12]/(2.**44)
        P2y = p[idx+14]/(2.**44)
        return p1x,p1y,p2x,p2y,P1x,P1y,P2x,P2y

    f = open(filename,'rb')
    pcap = dpkt.pcap.Reader(f)
    hdr = '<HHHHIHHHHHHHHHHHHHH4I4Q'
    pp2 = '8I8Q8I8Q8I8Q8I8Q8I8Q8I8Q8I8Q8I8Q'
    xfmt = '704h'
    khdr = ['HeaderLength','PacketNum','FFTShift','AccumLength','GlobalAccumNum',
            'BoardID','AccumNum','DataType','PolType','Ai','Aj','ADCOverflow',
            'QuantClipNum','NSubbands','iFreq','Delay0','Delay1','Delay2','Delay3',
            'PX0','PY0','PX1','PY1','P2X0','P2Y0','P2X1','P2Y1']
    # Start reading records
    out = pcap.readpkts()
    npkt = len(out)
    # Find number of P packets, and use it to determine how many accumulations were
    # captured in the file.
    t0 = out[0][0]  # initial X packet timestamp
    t1 = out[-1][0] # final packet timestamp
    nsec = int(t1) - int(t0) + 1
    # Find number of P packets, and use it to determine how many accumulations were
    # captured in the file.
    nPpkt = 0
    xlen = 0
    for a,b in out:
        if len(b) == 898: 
            nPpkt += 1
        elif xlen == 0:
            # Get length of X packets and value of first accumulation number
            xlen = len(b)
            header = struct.unpack(hdr,b[42:130])
            h = dict(zip(khdr,header))
            xaccum = h['AccumNum']
            pktnum = h['PacketNum']
            t0 = a  # Initial X packet timestamp
    if ptype is 'P':
        outarr = np.zeros([nsec,50,4096,8],'float')
        p1x = np.zeros(4096,'float')
        p1y = np.zeros(4096,'float')
        p2x = np.zeros(4096,'float')
        p2y = np.zeros(4096,'float')
        P1x = np.zeros(4096,'float')
        P1y = np.zeros(4096,'float')
        P2x = np.zeros(4096,'float')
        P2y = np.zeros(4096,'float')
    elif ptype is 'X':
        # In case of short capture packets, create the proper xfmt string.
        # The 130 is the 42-byte tcp packet header + 88 byte CASPER header.
        # The if statement is in case the capture buffer is not an even
        # number of 4-byte values, in which case stick on another half-word.
        nx = (xlen-130)/8
        x0 = 130  # Start location in packet buffer
        x1 = 130 + nx*8  # End location in packet buffer
        xfmt = '>'+str(2*nx)+'i'  # Read buffer as 32-bit signed integers
        # The format is all baselines, all poln products, for one channel
        outarr = np.zeros([nsec,50,4096,nx],'complex')
    else:
        print 'Invalid ptype: must be "P" or "X"'
        f.close()
        return False
    
    isec = 0
    a = -1
    for j in range(npkt):
        if verbose: print 'working on packet',j,'\r',
        t, buf = out[j]
        if ptype is 'P' and len(buf) == 898:
            # This is a P packet
            header = struct.unpack(hdr,buf[42:130])
            h = dict(zip(khdr,header))
            n = h['PacketNum']
            if  h['BoardID'] == boardID:
                if n == 0:
                    # Beginning of accumulation sample, so save the old one
                    # unless this is the first time through, indicated when a 
                    # is an impossible number (-1)
                    if a != -1:
                        outarr[isec,a,:,0] = p1x
                        outarr[isec,a,:,1] = p1y
                        outarr[isec,a,:,2] = p2x
                        outarr[isec,a,:,3] = p2y
                        outarr[isec,a,:,4] = P1x
                        outarr[isec,a,:,5] = P1y
                        outarr[isec,a,:,6] = P2x
                        outarr[isec,a,:,7] = P2y
                        p1x = np.zeros(4096,'float')
                        p1y = np.zeros(4096,'float')
                        p2x = np.zeros(4096,'float')
                        p2y = np.zeros(4096,'float')
                        P1x = np.zeros(4096,'float')
                        P1y = np.zeros(4096,'float')
                        P2x = np.zeros(4096,'float')
                        P2y = np.zeros(4096,'float')
                        if h['AccumNum'] == 0: 
                            # Beginning of a new second
                            isec += 1
                a = h['AccumNum']
                power = struct.unpack(pp2,buf[-768:])
                p1x_,p1y_,p2x_,p2y_,P1x_,P1y_,P2x_,P2y_ = pspectra(power)
                p1x[n*16:(n+1)*16] = p1x_ 
                p2x[n*16:(n+1)*16] = p2x_ 
                p1y[n*16:(n+1)*16] = p1y_ 
                p2y[n*16:(n+1)*16] = p2y_ 
                P1x[n*16:(n+1)*16] = P1x_ 
                P1y[n*16:(n+1)*16] = P1y_ 
                P2x[n*16:(n+1)*16] = P2x_ 
                P2y[n*16:(n+1)*16] = P2y_
        elif ptype is 'X' and len(buf) > 898:
            # This is an X packet from the production correlator
            header = struct.unpack(hdr,buf[42:130])
            h = dict(zip(khdr,header))
            n = h['PacketNum']
            # Override AccumNum with value based on packet timestamp, since
            # AccumNum is currently messed up.
            # This assumes packets coming out at 1-s mark were accumulated in
            # previous 20 ms (hence -2)
            aprev = a
            a = int(t*100 - 2) / 2 % 50  # This is calculated AccumNum
            if a == 0 and aprev == 49:
                isec += 1
                if isec == nsec:
                    sout = outarr.shape
                    outarr.shape = (sout[0]*sout[1],sout[2],sout[3])
                    return outarr                    
            iFreq = h['iFreq']
            xdata = np.array(struct.unpack(xfmt,buf[x0:x1]))
            outarr[isec,a,iFreq,:] = (xdata[::2] + 1j*xdata[1::2]).astype('complex64')
        else:
            # Some unknown packet?
            pass
    if verbose: print '\n'
    sout = outarr.shape
    outarr.shape = (sout[0]*sout[1],sout[2],sout[3])
    return outarr

def rd_jspec(filename):
    ''' Read all spectra in a Jim McT capture file
    
        Returns a dictionary with keys:
                                               ants pol chan slots
          'p': total power in p-packets, size [ 16,  2, 4096, 50  ]
          'p2': power-squared in p-packets,   [ 16,  2, 4096, 50  ]
          'a': auto-correlations in x-pkts,   [ 16,  4, 4096, 50  ]
                                              blines pol chan slots
          'x': cross-correlations in x-pkts, [ 120,   4, 4096, 50 ]
    '''
    import struct, os
    bl_order = get_bl_order(16)
    iauto = []
    icross = np.zeros((16,16),dtype='I')
    for i, bl in enumerate(bl_order):
        if bl[0] == bl[1]:
            iauto.append(i*4)
        else:
            icross[bl[0],bl[1]] = i*4
    iauto = np.array(iauto)
    icross = icross[icross.nonzero()]
    
    def pspectra(power):
        p = np.array(power)
        # idx is array([0,1, 16,17, 32,33, 48,49, 64,65, 80,81, 96,97, 112,113])
        idx = np.array(zip(np.arange(0,128,16),np.arange(1,128,16))).flatten()
        # Distribute P, P2 into their spectra
        p1x = p[idx]/(2.**14)
        p1y = p[idx+2]/(2.**14)
        p2x = p[idx+4]/(2.**14)
        p2y = p[idx+6]/(2.**14)
        P1x = p[idx+8]/(2.**44)
        P1y = p[idx+10]/(2.**44)
        P2x = p[idx+12]/(2.**44)
        P2y = p[idx+14]/(2.**44)
        return p1x,p1y,p2x,p2y,P1x,P1y,P2x,P2y

    hdr = '<HHHHIHHHHHHHHHHHHHH4I4Q'
    pp2 = '8I8Q8I8Q8I8Q8I8Q8I8Q8I8Q8I8Q8I8Q'
    # Read buffer as 136 channels of 4 polns, both real & imaginary 32-bit signed integers
    reclen = 136*4*2
    xfmt = '>'+str(reclen)+'i'
    khdr = ['HeaderLength','PacketNum','FFTShift','AccumLength','GlobalAccumNum',
            'BoardID','AccumNum','DataType','PolType','Ai','Aj','ADCOverflow',
            'QuantClipNum','NSubbands','iFreq','Delay0','Delay1','Delay2','Delay3',
            'PX0','PY0','PX1','PY1','P2X0','P2Y0','P2X1','P2Y1']
    # Get number of packets from file size
    npkt = os.stat(filename).st_size/4440
    outauto  = np.zeros(( 16, 4, 4096, 50),dtype='complex64')
    outcross = np.zeros((120, 4, 4096, 50),dtype='complex64')
    outp     = np.zeros(( 16, 2, 4096, 50),dtype='float'    )
    outp2    = np.zeros(( 16, 2, 4096, 50),dtype='float'    )
    # Start reading records
    with open(filename,'rb') as f:
        for i in range(npkt):
            buf = f.read(reclen*4+88)
            header = struct.unpack(hdr,buf[:88])
            h = dict(zip(khdr,header))
            if h['DataType'] == 0:
                # This is a P packet
                n = h['PacketNum']
                bid = h['BoardID']
                a = h['AccumNum']
                power = struct.unpack(pp2,buf[88:88+768])
                p1x,p1y,p2x,p2y,P1x,P1y,P2x,P2y = pspectra(power)
                outp[ bid*2,   0, n*16:(n+1)*16, a] = p1x
                outp[ bid*2,   1, n*16:(n+1)*16, a] = p1y
                outp[ bid*2+1, 0, n*16:(n+1)*16, a] = p2x
                outp[ bid*2+1, 1, n*16:(n+1)*16, a] = p2y
                outp2[bid*2,   0, n*16:(n+1)*16, a] = P1x
                outp2[bid*2,   1, n*16:(n+1)*16, a] = P1y
                outp2[bid*2+1, 0, n*16:(n+1)*16, a] = P2x
                outp2[bid*2+1, 1, n*16:(n+1)*16, a] = P2y
            elif h['DataType'] == 1:
                # This is an X packet
                a = h['AccumNum']
                iFreq = h['iFreq']
                xdata = np.array(struct.unpack(xfmt,buf[88:]))
                cxdata = (xdata[::2] + 1j*xdata[1::2]).astype('complex64')
                outauto[ :, 0, iFreq, a] = cxdata[iauto]
                outauto[ :, 1, iFreq, a] = cxdata[iauto+1]
                outauto[ :, 2, iFreq, a] = cxdata[iauto+2]
                outauto[ :, 3, iFreq, a] = cxdata[iauto+3]
                outcross[:, 0, iFreq, a] = cxdata[icross]
                outcross[:, 1, iFreq, a] = cxdata[icross+1]
                outcross[:, 2, iFreq, a] = cxdata[icross+2]
                outcross[:, 3, iFreq, a] = cxdata[icross+3]
            else:
                # Some unknown packet?
                pass
        out = {'p':outp,'p2':outp2,'a':outauto,'x':outcross}
    return out

def get_bl_order(n_ants):
    """Return the order of baseline data output by a CASPER correlator
    X engine."""
    order1, order2 = [], []
    for i in range(n_ants):
        for j in range(int(n_ants/2),-1,-1):
            k = (i-j) % n_ants
            if i >= k: order1.append((k, i))
            else: order2.append((i, k))
    order2 = [o for o in order2 if o not in order1]
    return tuple([o for o in order1 + order2])

def bl_list(nant=16):
    ''' Returns a two-dimensional array bl2ord that will translate
        a pair of antenna indexes (antenna number - 1) to the ordinal
        number of the baseline in the 'x' key.  Note bl2ord(i,j) = bl2ord(j,i),
        and bl2ord(i,i) = -1.
    '''
    bl2ord = np.ones((nant,nant),dtype='int')*(-1)
    k = 0
    for i in range(nant-1):
        for j in range(i+1,nant):
            bl2ord[i,j] = k
            bl2ord[j,i] = k
            k+=1
    return bl2ord
    
def summary_plot(out,ant_str='ant1-13',ptype='phase'):
    ''' Makes a summary amplitude or phase plot for all antennas from 0:nant
        in out dictionary.
    '''
    import matplotlib.pyplot as plt
    ant_list = ant_str2list(ant_str)
    nant = len(ant_list)
    if ptype != 'amp' and ptype != 'phase':
        print "Invalid plot type.  Must be 'amp' or 'phase'."
        return
    f, ax = plt.subplots(nant,nant)
    f.subplots_adjust(hspace=0,wspace=0)
    bl2ord = bl_list()
    for axrow in ax:
        for a in axrow:
            a.xaxis.set_visible(False)
            a.yaxis.set_visible(False)
    for i in range(nant-1):
        ai = ant_list[i]
        for j in range(i+1,nant):
            aj = ant_list[j]
            if ptype == 'phase':
                ax[i,j].imshow(np.angle(out['x'][bl2ord[ai,aj],0].reshape(64,64,10,5).sum(1).sum(2)/(64*5)))
                ax[j,i].imshow(np.angle(out['x'][bl2ord[ai,aj],1].reshape(64,64,10,5).sum(1).sum(2)/(64*5)))
            elif ptype == 'amp':
                ax[i,j].imshow(np.abs(out['x'][bl2ord[ai,aj],0].reshape(64,64,10,5).sum(1).sum(2)/(64*5)))
                ax[j,i].imshow(np.abs(out['x'][bl2ord[ai,aj],1].reshape(64,64,10,5).sum(1).sum(2)/(64*5)))
    for i in range(nant):
        ai = ant_list[i]
        ax[i,i].text(0.5,0.5,str(ai+1),ha='center',va='center',transform=ax[i,i].transAxes,fontsize=14)

def capture(filename='dump',nsec=1,ptype='P',snaplen=5000,overwrite=True):
    ''' All-in-one command to capture a given number of seconds of packets on
        both interfaces and return the result in a single array.  
        Arguments: filename   create or read a file "eth<n>_<filename>.pcap"
                                where <n> is 2 and 3 for the two interfaces,
                                and <filename> is the given string.
                   nsec       number of seconds to capture packets
                   ptype      either 'P' (default) or 'X'
                   overwrite  if True (default), new packets will be grabbed,
                                otherwise an existing file will be read.
        Returns:   out        a directionary that contains two numpy arrays whose size depends of ptype
                                ptype = 'P'   P packages
                    'p': total power with dimension (nant,npol,nt,nf) = (16,2,50*(nsec+1),4096)
                    'p2': total power square with same dimension as above
                                ptype = 'X'   X packets
                    'a': auto-correlation (nant,npol,nt,nf) = (16,4,50*(nsec+2),4096)
                    'x': cross-correlation (nbl,npol,nt,nf) = (120,4,50*(nsec+2),4096)
    '''
    import glob
    if ptype == 'P': snaplen = 1000
    if overwrite or glob.glob('eth2_'+filename+'.pcap') == []:
        command = 'tcpdump -i eth2 -c '+str(153600*nsec)+' -w eth2_'+filename+'.pcap -s '+str(snaplen)
        sendcmd(command)
    if overwrite or glob.glob('eth3_'+filename+'.pcap') == []:
        command = 'tcpdump -i eth3 -c '+str(153600*nsec)+' -w eth3_'+filename+'.pcap -s '+str(snaplen)
        sendcmd(command)
    if ptype == 'P':
        out1 = rd_spec('eth2_'+filename+'.pcap',ptype=ptype,boardID=0)
        out2 = rd_spec('eth2_'+filename+'.pcap',ptype=ptype,boardID=1)
        out5 = rd_spec('eth2_'+filename+'.pcap',ptype=ptype,boardID=4)
        out6 = rd_spec('eth2_'+filename+'.pcap',ptype=ptype,boardID=5)
        out3 = rd_spec('eth3_'+filename+'.pcap',ptype=ptype,boardID=2)
        out4 = rd_spec('eth3_'+filename+'.pcap',ptype=ptype,boardID=3)
        out7 = rd_spec('eth3_'+filename+'.pcap',ptype=ptype,boardID=6)
        out8 = rd_spec('eth3_'+filename+'.pcap',ptype=ptype,boardID=7)
        outall=[out1,out2,out3,out4,out5,out6,out7,out8]
        nt2, nf, nif = out2.shape #captured package dimensions for interface eth2
        nt3, nf, nif = out3.shape #captured package dimensions for interface eth3
        if nt2 == nt3:
            print 'packets from eth2 and eth3 have the same size '+str(nt2)+', happily proceeding...'
        if nt2 < nt3:
            print 'packets from eth2 have less time points '+str(nt2)+', resizing to that from eth3 '+str(nt3)+'...' 
            out1.resize((nt3,nf,nif),refcheck=False)
            out2.resize((nt3,nf,nif),refcheck=False)
            out5.resize((nt3,nf,nif),refcheck=False)
            out6.resize((nt3,nf,nif),refcheck=False)
        if nt2 > nt3:
            print 'packets from eth3 have less time points '+str(nt3)+', resizing to that from eth2 '+str(nt2)+'...'
            out3.resize((nt2,nf,nif),refcheck=False)
            out4.resize((nt2,nf,nif),refcheck=False)
            out7.resize((nt2,nf,nif),refcheck=False)
            out8.resize((nt2,nf,nif),refcheck=False)
        for i in range(8):
            outall[i] = np.rollaxis(outall[i],2,0)
            #out = np.rollaxis(out,2,1)
            nif, nt, nf = outall[i].shape
            outall[i].shape = (2,nif/4,2,nt,nf)
        out = np.concatenate(outall,1)
        outp=out[0,:,:,:,:]
        outp2=out[1,:,:,:,:]
        return {'p':outp,'p2':outp2}
    if ptype == 'X':
        out2 = rd_spec('eth2_'+filename+'.pcap',ptype=ptype)
        out3 = rd_spec('eth3_'+filename+'.pcap',ptype=ptype)
        nt2, nf, nch = out2.shape
        nt3, nf, nch = out3.shape
        if nt2 == nt3:
            print 'packets from eth2 and eth3 have the same size '+str(nt2)+', happily proceeding...'
        if nt2 < nt3:
            print 'packets from eth2 have less time points '+str(nt2)+', resizing to that from eth3 '+str(nt3)+'...' 
            out2.resize((nt3,nf,nch),refcheck=False)
        if nt2 > nt3:
            print 'packets from eth3 have less time points '+str(nt3)+', resizing to that from eth2 '+str(nt2)+'...'
            out3.resize((nt2,nf,nch),refcheck=False)
        out = out2+out3
        out = np.rollaxis(out,2,0)
        #out = np.rollaxis(out,2,1)
        nch, nt, nf = out.shape
        nbl=nch/4
        out.shape = (nbl,4,nt,nf)
        #find out which ones are auto-correlation and cross-correlation
        idxa=[]
        bls=get_bl_order(16)   
        for (i,bl) in enumerate(bls):
            if bl[0] == bl[1]:
                idxa.append(i) #all indices for auto-correlation
        #find out which indices are actually in the range of the captured data
        idxa=np.array(idxa)
        ix=np.where(idxa < (nbl-1))
        outa=out[idxa[ix],:,:,:]
        outx=np.delete(out,idxa[ix],0)
        return {'a':outa,'x':outx}

def get_spec(capfile=None, mode='jcap'):
    '''basic wrapper for capture and rd_jspec to read existing capture/Jim capture files and get unified output
        OUTPUT: po     total power
                ac     auto correlation
                xc     cross correlation
        example: (po, ac, xc) = get_spec(capfile=capfile, mode='jcap')
    '''
    po=[]
    ac=[]
    xc=[]
    if not isinstance(capfile,list):
        print 'capfile is required to be a list of file names'
        return
    nfile=len(capfile)
    if nfile < 1:
        print 'capfile has 0 element'
        return
    for n in range(nfile):
        if mode=='cap':
            outp=capture(filename=capfile[n], ptype='P', overwrite=False)
            outx=capture(filename=capfile[n], ptype='X', snaplen=1000, overwrite=False)
            po.append(outp['p'])
            ac.append(outx['a'])
            xc.append(outx['x'])
        if mode=='jcap':
            out=rd_jspec(capfile[n])
            po.append(out['p']) #power
            ac.append(out['a']) #auto-correlation
            xc.append(out['x']) #auto-correlation
    return po, ac, xc

def bl_mapper(nant=16, refant=1, show=False, exclude_auto=True):
    '''For selected antenna, calculate baseline index and return baseline pair name in the x-correlation output 
    '''
    nbl=nant*(nant-1)/2
    blprs=[]
    if exclude_auto:
        ncks=[nant-1-i for i in range(nant-1)]
        for ant in range(nant-1):
            for i in range(ncks[ant]):
                blprs.append([ant+1,ant+2+i])
    else:
        ncks=[nant-i for i in range(nant)]
        for ant in range(nant):
            for i in range(ncks[ant]):
                blprs.append([ant+1,ant+1+i])

    bl_idx=[]
    bl_list=[]
    for i in range(nbl):
        blpr=blprs[i]
        if refant in blpr:
            bl_idx.append(i)
            bl_list.append(str(blpr[0])+' & '+str(blpr[1]))
            if show:
                print blpr
    return {'idx':bl_idx, 'name':bl_list}

def plot_auto_corr(ac=None, nant=16, nx=2, chran='0-4095'):
    '''Reads the outputs from get_spec and plot auto correlations on selected antenna
    '''
    import matplotlib.pyplot as plt
    nfile=len(ac) 
    if isinstance(chran, basestring):
        (ch1,ch2) = (int(s) for s in chran.split('-'))
        if ch1 > ch2 or ch1 < 0 or ch2 > 4095:
            print 'start channel no must be less than end channel #'
            print 'start and end channel must be >= 0 and <= 4096'
            print 'use the default range 0-4095 instead'
            ch1 = 0
            ch2 = 4095
    else:
        print 'chran has to be a string'
        return
    for n in range(nfile):
        ac_=ac[n]
        nbl_,nx_,nch,nt=ac_.shape
        #calculate the baseline numbers for all antennas
        f,ax=plt.subplots(nx,nant)
        for j in range(nx):
            for i in range(nant): 
                ax[j,i].imshow(abs(ac_[i,j,ch1:ch2,:]),extent=[0,nt-1,ch1,ch2])
                ax[j,i].set_title('ant '+str(i+1))
                if i > 0:
                    ax[j,i].get_yaxis().set_ticks([])

def plot_x_corr(xc=None, plot_type='pha', nbl=15, nx=2, refant=1, chran='2048-4095',lineplot=False,tidx=25):
    '''Read the outputs from get_spec and plot cross-correlated phases or amplitudes on selected antenna
        REQUIRED INPUTS: xc           x-correlation output from get_spec
        OPTIONAL INPUTS: plot_type    'pha' or 'amp' for phase or amplitude display
                         nbl          number of baselines to display
                         nx           number of polarizations to display
                         refant       which antenna to plot? refant=1 means antenna 1, etc.
                         chran        channel range to display (within the 4096 raw channels)
                         lineplot     do you wish to do a lineplot (amp/phase vs channel)?
                         tidx         time index for the time plot, default to 25
    '''
    import matplotlib.pyplot as plt
    if isinstance(chran, basestring):
        (ch1,ch2) = (int(s) for s in chran.split('-'))
        if ch1 > ch2 or ch1 < 0 or ch2 > 4095:
            print 'start channel no must be less than end channel #'
            print 'start and end channel must be >= 0 and <= 4096'
            print 'use the default range 2048-4095 instead'
            ch1 = 2048
            ch2 = 4095
    else:
        print 'chran has to be a string'
        return
    if not xc:
        print 'please provide a cross-power spectrum file'
        print "by (po,ac,xc)=get_spec(capfile,'jcap')"
    nfile=len(xc) 
    bl=bl_mapper(refant=refant)
    bli=bl['idx']
    bln=bl['name']
    for n in range(nfile):
        xc_=xc[n]
        nbl_,nx_,nch,nt=xc_.shape
        #calculate the baseline numbers for all antennas
        f,ax=plt.subplots(nx,nbl)
        for j in range(nx):
            for i in range(nbl): 
                if plot_type=='pha':
                    ax[j,i].imshow(np.angle(xc_[bli[i],j,ch1:ch2,:]),vmin=-np.pi,vmax=np.pi,
                                   extent=[0,nt-1,ch1,ch2])
                    ax[j,i].set_title(bln[i])
                    if i > 0:
                        ax[j,i].get_yaxis().set_ticks([])
                if plot_type=='amp':
                    ax[j,i].imshow(abs(xc_[bli[i],j,ch1:ch2,:]),
                                   extent=[0,nt-1,ch1,ch2])
                    ax[j,i].set_title(bln[i])
                    if i > 0:
                        ax[j,i].get_yaxis().set_ticks([])

        if lineplot:
            f,ax=plt.subplots(nx,nbl)
            chan=np.arange(4096)[ch1:ch2]
            for j in range(nx):
                for i in range(nbl): 
                    if plot_type=='pha':
                        ax[j,i].plot(np.angle(xc_[bli[i],j,ch1:ch2,tidx]),chan,'.')
                        ax[j,i].set_xlim([-np.pi,np.pi])
                        ax[j,i].set_ylim([ch1,ch2])
                        ax[j,i].set_title(bln[i])
                        if i > 0:
                            ax[j,i].get_yaxis().set_ticks([])
                    if plot_type=='amp':
                        ax[j,i].plot(abs(xc_[bli[i],j,ch1:ch2,tidx]),chan,'.')
                        ax[j,i].set_ylim([ch1,ch2])
                        ax[j,i].set_title(bln[i])
                        if i > 0:
                            ax[j,i].get_yaxis().set_ticks([])

   
def capture_fig(useroach=[1,2],print_attn=False):

    import matplotlib.pyplot as plt
    from util import Time
    
    # Select which pair of ROACH boards is used for plot

    # Grab 1.02 s of data (19200 + 384 packets)*2 to make sure we have at least 1 s of good data
    npkts = (19200 + 384)*3
    if useroach == [1,2,5]:
        iface = 'eth2'
    else:
        iface = 'eth3'
    ret = sendcmd('/usr/sbin/tcpdump -i '+iface+' -c '+str(npkts)+' -w /home/user/Python/'+iface+'.pcap -s 2000')

    out1 = rd_spec(iface+'.pcap',boardID=useroach[0]-1)
    out2 = rd_spec(iface+'.pcap',boardID=useroach[1]-1)
    out3 = rd_spec(iface+'.pcap',boardID=useroach[2]-1)

    lines = np.array(list_header(iface+'.pcap',boardID=useroach[0]-1))
    idx1 = np.where(np.char.startswith(lines,'Acc'))[0][0]
    lines = np.append(lines[idx1+1:],lines[:idx1])
    if idx1 == 0:
        out1 = out1[0:50,:,:]
    else:
        out1 = out1[50-idx1:100-idx1,:,:]
    ovfl1 = []
    for line in lines:
        ovfl1.append(int(line[30:36]))
    lines = np.array(list_header(iface+'.pcap',boardID=useroach[1]-1))
    idx2 = np.where(np.char.startswith(lines,'Acc'))[0][0]
    lines = np.append(lines[idx2+1:],lines[:idx2])
    if idx2 == 0:
        out2 = out2[0:50,:,:]
    else:
        out2 = out2[50-idx2:100-idx2,:,:]
    ovfl2 = []
    for line in lines:
        ovfl2.append(int(line[30:36]))
    lines = np.array(list_header(iface+'.pcap',boardID=useroach[2]-1))
    idx3 = np.where(np.char.startswith(lines,'Acc'))[0][0]
    lines = np.append(lines[idx3+1:],lines[:idx3])
    if idx3 == 0:
        out3 = out3[0:50,:,:]
    else:
        out3 = out3[50-idx3:100-idx3,:,:]
    ovfl3 = []
    for line in lines:
        ovfl3.append(int(line[30:36]))

    # Reorganize data order to correspond to the 34 bands of solar.fsq,
    # putting the repeated bands 1,2,3,4 together for averaging
    idx = np.array([0,10,20,30,40,1,11,21,31,41,2,12,22,32,42,3,13,23,33,43,
                 4,5,6,7,8,9,14,15,16,17,18,19,24,25,26,27,28,29,34,35,36,
                 37,38,39,44,45,46,47,48,49])
    out = np.concatenate((out1[(idx+idx1)%50,:,0:4],out2[(idx+idx1)%50,:,0:4],out3[(idx+idx1)%50,:,0:4]),axis=2)
    # Average repeated channels and place into slots 16-19
    out[19,:,:] = out[15:20,:,:].sum(0)/5.
    out[18,:,:] = out[10:15,:,:].sum(0)/5.
    out[17,:,:] = out[5:10,:,:].sum(0)/5.
    out[16,:,:] = out[0:5,:,:].sum(0)/5.
    # Truncate to slots 16-49 (34 slots)
    out = out[16:,:,:]
    
    if print_attn:
        # At this point, for solar sequence, out[0:34,:,0:12] corresponds to
        # the power in each antenna/feed on these three ROACHes.  Print out the sum
        # in a nice tabular format:
        attn = []
        for i in range(12):
            attn_list = (10*np.log10(out[:,:,i].sum(1)/7.0)).astype('int')/2 * 2
            #bad = np.where(attn_list < 0)[0]
            #if len(bad) > 0:
            #    attn_list[bad] = 0
            # Append list of 34 attenuations
            attn.append(attn_list)
        attn = np.array(attn)
        attn.shape = (12,34)
        return attn
    ovfl1 = np.array(ovfl1)[idx]
    ovfl2 = np.array(ovfl2)[idx]
    ovfl3 = np.array(ovfl3)[idx]
    ovfl1[19] = ovfl1[15:20].sum()/5.
    ovfl1[18] = ovfl1[10:15].sum()/5.
    ovfl1[17] = ovfl1[5:10].sum()/5.
    ovfl1[16] = ovfl1[0:5].sum()/5.
    ovfl2[19] = ovfl2[15:20].sum()/5.
    ovfl2[18] = ovfl2[10:15].sum()/5.
    ovfl2[17] = ovfl2[5:10].sum()/5.
    ovfl2[16] = ovfl2[0:5].sum()/5.
    ovfl3[19] = ovfl3[15:20].sum()/5.
    ovfl3[18] = ovfl3[10:15].sum()/5.
    ovfl3[17] = ovfl3[5:10].sum()/5.
    ovfl3[16] = ovfl3[0:5].sum()/5.
    ovfl = np.concatenate((ovfl1,ovfl2,ovfl3))
    ovfl.shape = (50,3)
    return Time.now(),out,ovfl

def plot_fig(t,out1,out2,ovfl1,ovfl2):
    out = np.concatenate((out1[:,:,0:8],out2[:,:,0:8],out1[:,:,8:12],out2[:,:,8:12]),axis=2)
    ovfl = np.concatenate((ovfl1[:,0:2],ovfl2[:,0:2],ovfl1[:,2],ovfl2[:,2]),axis=1)
    f, ax = plt.subplots(6, 4, sharex='col', sharey='row')
    f.set_size_inches(10,12,forward=True)
    f.suptitle('Packet Capture at '+t.iso[:19]+' UT',fontsize=18)
    #a1 = str((useroach[0]-1)*2+1)
    #a2 = str((useroach[0]-1)*2+2)
    #a3 = str((useroach[1]-1)*2+1)
    #a4 = str((useroach[1]-1)*2+2)
    titstr = np.array([['Ant 1 X','Ant 1 Y','Ant 2 X','Ant 2 Y'],['Ant 3 X','Ant 3 Y','Ant 4 X','Ant 4 Y'],
              ['Ant 5 X','Ant 5 Y','Ant 6 X','Ant 6 Y'],['Ant 7 X','Ant 7 Y','Ant 8 X','Ant 8 Y'],
              ['Ant 9 X','Ant 9 Y','Ant 10 X','Ant 10 Y'],['Ant 11 X','Ant 11 Y','Ant 12 X','Ant 12 Y']])
    for i in range(4):
        ax[3,i].set_xlabel('Band Number')
    for j in range(6):
        ax[j,0].set_ylabel('Channel')
    for j in range(6):
        for i in range(4):
            ax[j,i].plot(ovfl[16:,j]/10,'o',color='white')
    ax[0,0].imshow(np.clip(np.transpose(out1[:,:,0]),0,np.median(out1[:,:,0])*5),aspect='auto',interpolation='nearest')
    ax[0,0].set_title(titstr[0,0])
    ax[0,1].plot(ovfl3[16:]/10,'o',color='white')
    ax[0,1].imshow(np.clip(np.transpose(out2[:,:,0]),0,np.median(out2[:,:,0])*5),aspect='auto',interpolation='nearest')
    ax[0,1].set_title(titstr[0,1])
    ax[0,2].plot(ovfl1[16:]/10,'o',color='white')
    ax[0,2].imshow(np.clip(np.transpose(out1[:,:,1]),0,np.median(out1[:,:,1])*5),aspect='auto',interpolation='nearest')
    ax[0,2].set_title(titstr[0,2])
    ax[0,3].plot(ovfl3[16:]/10,'o',color='white')
    ax[0,3].imshow(np.clip(np.transpose(out2[:,:,1]),0,np.median(out2[:,:,1])*5),aspect='auto',interpolation='nearest')
    ax[0,3].set_title(titstr[0,3])
    ax[1,0].plot(ovfl1[16:]/10,'o',color='white')
    ax[1,0].imshow(np.clip(np.transpose(out1[:,:,2]),0,np.median(out1[:,:,2])*5),aspect='auto',interpolation='nearest')
    ax[1,0].set_title(titstr[1,0])
    ax[1,1].plot(ovfl3[16:]/10,'o',color='white')
    ax[1,1].imshow(np.clip(np.transpose(out2[:,:,2]),0,np.median(out2[:,:,2])*5),aspect='auto',interpolation='nearest')
    ax[1,1].set_title(titstr[1,1])
    ax[1,2].plot(ovfl1[16:]/10,'o',color='white')
    ax[1,2].imshow(np.clip(np.transpose(out1[:,:,3]),0,np.median(out1[:,:,3])*5),aspect='auto',interpolation='nearest')
    ax[1,2].set_title(titstr[1,2])
    ax[1,3].plot(ovfl3[16:]/10,'o',color='white')
    ax[1,3].imshow(np.clip(np.transpose(out2[:,:,3]),0,np.median(out2[:,:,3])*5),aspect='auto',interpolation='nearest')
    ax[1,3].set_title(titstr[1,3])
    ax[2,0].plot(ovfl2[16:]/10,'o',color='white')
    ax[2,0].imshow(np.clip(np.transpose(out1[:,:,4]),0,np.median(out1[:,:,4])*5),aspect='auto',interpolation='nearest')
    ax[2,0].set_title(titstr[2,0])
    ax[2,1].plot(ovfl4[16:]/10,'o',color='white')
    ax[2,1].imshow(np.clip(np.transpose(out2[:,:,4]),0,np.median(out2[:,:,4])*5),aspect='auto',interpolation='nearest')
    ax[2,1].set_title(titstr[2,1])
    ax[2,2].plot(ovfl2[16:]/10,'o',color='white')
    ax[2,2].imshow(np.clip(np.transpose(out1[:,:,5]),0,np.median(out1[:,:,5])*5),aspect='auto',interpolation='nearest')
    ax[2,2].set_title(titstr[2,2])
    ax[2,3].plot(ovfl4[16:]/10,'o',color='white')
    ax[2,3].imshow(np.clip(np.transpose(out2[:,:,5]),0,np.median(out2[:,:,5])*5),aspect='auto',interpolation='nearest')
    ax[2,3].set_title(titstr[2,3])
    ax[3,0].plot(ovfl2[16:]/10,'o',color='white')
    ax[3,0].imshow(np.clip(np.transpose(out1[:,:,6]),0,np.median(out1[:,:,6])*5),aspect='auto',interpolation='nearest')
    ax[3,0].set_title(titstr[3,0])
    ax[3,1].plot(ovfl4[16:]/10,'o',color='white')
    ax[3,1].imshow(np.clip(np.transpose(out2[:,:,6]),0,np.median(out2[:,:,6])*5),aspect='auto',interpolation='nearest')
    ax[3,1].set_title(titstr[3,1])
    ax[3,2].plot(ovfl2[16:]/10,'o',color='white')
    ax[3,2].imshow(np.clip(np.transpose(out1[:,:,7]),0,np.median(out1[:,:,7])*5),aspect='auto',interpolation='nearest')
    ax[3,2].set_title(titstr[3,2])
    ax[3,3].plot(ovfl4[16:]/10,'o',color='white')
    ax[3,3].imshow(np.clip(np.transpose(out2[:,:,7]),0,np.median(out2[:,:,7])*5),aspect='auto',interpolation='nearest')
    ax[3,3].set_title(titstr[3,3])
#    for i in range(8):
#        if i < 4:
#            ax[i/2,i % 2].plot(ovfl[16:]/10,'o',color='white')
#        else:
#            ax[i/2,i % 2].plot(ovfl2[16:]/10,'o',color='white')
#        ax[i/2,i % 2].imshow(np.clip(np.transpose(out[:,:,i]),0,np.median(out[:,:,i])*5),aspect='auto',interpolation='nearest')
#        ax[i/2,i % 2].set_title(titstr[i])
    plt.savefig('/common/tmp/dppcapture.png',bbox_inches='tight')
    return f

if __name__ == "__main__":
    # For non-interactive use, use a backend that does not require a display
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import sys
    if len(sys.argv) == 3:
        if sys.argv[1] == 'print_attn':
            attn1 = capture_fig(useroach=[1,2,5],print_attn=True)
            attn2 = capture_fig(useroach=[3,4,6],print_attn=True)
            attn = np.zeros((30,34),dtype='int')
            val = sys.argv[2].split()
            if len(val) == 16:
                # An array of 16 attenuation values was provided, so add them to each line of output
                vals = np.zeros(len(val),'int')
                for i,v in enumerate(val):
                    vals[i] = int(v)
                for j in range(34):
                    attn[0:8,j] = attn1[0:8,j] + vals[0:8]
                    attn[8:16,j] = attn2[0:8,j] + vals[8:16]
                    attn[16:20,j] = attn1[8:12,j] + vals[16:20]
                    attn[20:24,j] = attn2[8:12,j] + vals[20:24]
                    attn[24:,j] = 10
            else:
                # Assume the second argument is only a single number
                attn[0:8] = attn1[0:8] + int(sys.argv[2])
                attn[8:16] = attn2[0:8] + int(sys.argv[2])
                attn[16:20] = attn1[8:12] + int(sys.argv[2])
                attn[20:24] = attn2[8:12] + int(sys.argv[2])
                attn[24:] = 10
            attn.shape = (30*34)
            bad = np.where(attn < 0)[0]
            if len(bad) > 0:
                attn[bad] = 0
            attn.shape = (30,34)

            print '       Ant1  Ant2  Ant3  Ant4  Ant5  Ant6  Ant7  Ant8  Ant9 Ant10 Ant11 Ant12 Ant13 Ant14 Ant15'
            print '       X  Y  X  Y  X  Y  X  Y  X  Y  X  Y  X  Y  X  Y  X  Y  X  Y  X  Y  X  Y  X  Y  X  Y  X  Y'
            print '      ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- -----'
            for i in range(34):
                print '{:2} :  {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2}'.format(i+1,*attn[:,i])
        else:
            print 'Command line argument',sys.argv[1],'not understood.  Exiting.'
    else:
        t1,out1,ovfl1 = capture_fig(useroach=[1,2,5])
        t2,out2,ovfl2 = capture_fig(useroach=[3,4,6])
        plt.close(plot_fig(t1,out1,out2,ovfl1,ovfl2))
