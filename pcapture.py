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
#

import numpy as np

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

def rd_spec(filename,ptype='P',boardID=0,nboards=2,verbose=False,proto=True):
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

    def xspectra(xdata):
        # Case of prototype X spectra
        x = np.array(xdata)
        # idx is array([0,44,88,132,176,220,264,308,352,396,440,484,528,572,616,660])
        idx = np.arange(0,704,44)
        # Distribute 16 channels of complex correlation data into their spectra
        #idx += 4   # Uncomment this to read xy and yx data
        a12x = x[idx]+x[idx+1]*1j
        a12y = x[idx+2]+x[idx+3]*1j
        a13x = x[idx+8]+x[idx+9]*1j
        a13y = x[idx+10]+x[idx+11]*1j
        a14x = x[idx+16]+x[idx+17]*1j
        a14y = x[idx+18]+x[idx+19]*1j
        a23x = x[idx+24]+x[idx+25]*1j
        a23y = x[idx+26]+x[idx+27]*1j
        a24x = x[idx+32]+x[idx+33]*1j
        a24y = x[idx+34]+x[idx+35]*1j
        #idx -= 4   # Uncomment this to read xy and yx data
        a34x = x[idx+40]+x[idx+41]*1j
        a34y = x[idx+42]+x[idx+43]*1j
        return a12x,a12y,a13x,a13y,a14x,a14y,a23x,a23y,a24x,a24y,a34x,a34y
    
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
    t1 = a #Final packet timestamp
            
    naccum = nPpkt / 256 / nboards
    nsec = int(t1 - t0) + 1
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
    elif ptype is 'X' and proto:
        outarr = np.zeros([nsec,50,4096,12],'complex')
        a12x = np.zeros(4096,'complex')
        a12y = np.zeros(4096,'complex')
        a13x = np.zeros(4096,'complex')
        a13y = np.zeros(4096,'complex')
        a14x = np.zeros(4096,'complex')
        a14y = np.zeros(4096,'complex')
        a23x = np.zeros(4096,'complex')
        a23y = np.zeros(4096,'complex')
        a24x = np.zeros(4096,'complex')
        a24y = np.zeros(4096,'complex')
        a34x = np.zeros(4096,'complex')
        a34y = np.zeros(4096,'complex')
    elif ptype is 'X':
        # In case of short capture packets, create the proper xfmt string.
        # The 130 is the 42-byte tcp packet header + 88 byte CASPER header.
        # The if statement is in case the capture buffer is not an even
        # number of 4-byte values, in which case stick on another half-word.
        nx = (xlen-130)/8
        x0 = 130  # Start location in packet buffer
        x1 = 130 + nx*8  # End location in packet buffer
        xfmt = str(2*nx)+'i'  # Read buffer as 32-bit signed integers
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
            if proto:
                # This is an X packet from the prototype
                header = struct.unpack(hdr,buf[42:130])
                h = dict(zip(khdr,header))
                n = h['PacketNum']
                # Case of prototype packets
                if  h['BoardID'] == boardID:
                    if n == 0:
                        # Beginning of accumulation sample, so save the old one
                        # unless this is the first time through, indicated when a 
                        # is an impossible number (-1)
                        if a != -1:
                            outarr[isec,a,:,0] = a12x
                            outarr[isec,a,:,1] = a12y
                            outarr[isec,a,:,2] = a13x
                            outarr[isec,a,:,3] = a13y
                            outarr[isec,a,:,4] = a14x
                            outarr[isec,a,:,5] = a14y
                            outarr[isec,a,:,6] = a23x
                            outarr[isec,a,:,7] = a23y
                            outarr[isec,a,:,8] = a24x
                            outarr[isec,a,:,9] = a24y
                            outarr[isec,a,:,10] = a34x
                            outarr[isec,a,:,11] = a34y
                            a12x = np.zeros(4096,'complex')
                            a12y = np.zeros(4096,'complex')
                            a13x = np.zeros(4096,'complex')
                            a13y = np.zeros(4096,'complex')
                            a14x = np.zeros(4096,'complex')
                            a14y = np.zeros(4096,'complex')
                            a23x = np.zeros(4096,'complex')
                            a23y = np.zeros(4096,'complex')
                            a24x = np.zeros(4096,'complex')
                            a24y = np.zeros(4096,'complex')
                            a34x = np.zeros(4096,'complex')
                            a34y = np.zeros(4096,'complex')
                            if h['AccumNum'] == 0: 
                                # Beginning of a new second
                                isec += 1
                                if isec == nsec:
                                    sout = outarr.shape
                                    outarr.shape = (sout[0]*sout[1],sout[2],sout[3])
                                    return outarr                    
                    a = h['AccumNum']
                    xdata = struct.unpack(xfmt,buf[130:])
                    a12x_,a12y_,a13x_,a13y_,a14x_,a14y_,a23x_,a23y_,a24x_,a24y_,a34x_,a34y_ = xspectra(xdata)
                    bid = h['BoardID'] % 2
                    ns = n*16 + bid*2048
                    ne = (n+1)*16 + bid*2048
                    a12x[ns:ne] = a12x_ 
                    a12y[ns:ne] = a12y_ 
                    a13x[ns:ne] = a13x_ 
                    a13y[ns:ne] = a13y_ 
                    a14x[ns:ne] = a14x_ 
                    a14y[ns:ne] = a14y_ 
                    a23x[ns:ne] = a23x_ 
                    a23y[ns:ne] = a23y_ 
                    a24x[ns:ne] = a24x_ 
                    a24y[ns:ne] = a24y_ 
                    a34x[ns:ne] = a34x_ 
                    a34y[ns:ne] = a34y_ 
            else:
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
  
def capture(filename='dump',nsec=1,ptype='P',overwrite=True):
    ''' All-in-one command to capture a given number of seconds of packets on
        both interfaces and return the result in a single array.  
        Arguments: filename   create or read a file "eth<n>_<filename>.pcap"
                                where <n> is 2 and 3 for the two interfaces,
                                and <filename> is the given string.
                   nsec       number of seconds to capture packets
                   ptype      either 'P' (default) or 'X'
                   overwrite  if True (default), new packets will be grabbed,
                                otherwise an existing file will be read.
        Returns:   out        a numpy array whose size depends of ptype
                                ptype = 'P' => (nant,npol,nf,nt) = (8,2,4096,50*(nsec+1))
                                ptype = 'X' => (nbl,npol,nf,nt) = (12,2,4096,50*(nsec+1))
    '''
    import glob
    if overwrite or glob.glob('eth2_'+filename+'.pcap') == []:
        command = 'tcpdump -i eth2 -c '+str(38400*nsec)+' -w eth2_'+filename+'.pcap -s 2000'
        sendcmd(command)
    if overwrite or glob.glob('eth3_'+filename+'.pcap') == []:
        command = 'tcpdump -i eth3 -c '+str(38400*nsec)+' -w eth3_'+filename+'.pcap -s 2000'
        sendcmd(command)
    if ptype == 'P':
        out1 = rd_spec('eth2_'+filename+'.pcap',ptype=ptype,boardID=0)
        out2 = rd_spec('eth2_'+filename+'.pcap',ptype=ptype,boardID=1)
        out3 = rd_spec('eth3_'+filename+'.pcap',ptype=ptype,boardID=2)
        out4 = rd_spec('eth3_'+filename+'.pcap',ptype=ptype,boardID=3)
        out1[:,:,4:] = out2[:,:,:4]
        out3[:,:,4:] = out4[:,:,:4]
        del out2,out4
        out = np.concatenate((out1,out3),2)
    else:
        out1 = rd_spec('eth2_'+filename+'.pcap',ptype=ptype,boardID=0)
        out2 = rd_spec('eth2_'+filename+'.pcap',ptype=ptype,boardID=1)
        out3 = rd_spec('eth3_'+filename+'.pcap',ptype=ptype,boardID=2)
        out4 = rd_spec('eth3_'+filename+'.pcap',ptype=ptype,boardID=3)
        out1[:,2048:,:] = out2[:,2048:,:]
        out3[:,2048:,:] = out4[:,2048:,:]
        del out2,out4
        out = np.concatenate((out1,out3),2)
    out = np.rollaxis(out,2,0)
    out = np.rollaxis(out,2,1)
    nif, nf, nt = out.shape
    out.shape = (nif/2,2,nf,nt)
    return out
   
def capture_fig(useroach=[1,2],print_attn=False):

    import matplotlib.pyplot as plt
    from util import Time
    
    # Select which pair of ROACH boards is used for plot

    # Grab 1.02 s of data (19200 + 384 packets)*2 to make sure we have at least 1 s of good data
    if useroach == [1,2]:
        iface = 'eth2'
    else:
        iface = 'eth3'
    ret = sendcmd('/usr/sbin/tcpdump -i '+iface+' -c 39168 -w /home/user/Python/'+iface+'.pcap -s 2000')
    out = rd_spec(iface+'.pcap',boardID=useroach[0]-1)
    #out.shape = (100,4096,8)
    out2 = rd_spec(iface+'.pcap',boardID=useroach[1]-1)
    #out2.shape = (100,4096,8)
    lines = np.array(list_header(iface+'.pcap',boardID=useroach[0]-1))
    idx1 = np.where(np.char.startswith(lines,'Acc'))[0][0]
    lines = np.append(lines[idx1+1:],lines[:idx1])
    if idx1 == 0:
        out = out[0:50,:,:]
    else:
        out = out[50-idx1:100-idx1,:,:]
    ovfl = []
    for line in lines:
        ovfl.append(int(line[30:36]))
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

    # Reorganize data order to correspond to the 34 bands of solar.fsq,
    # putting the repeated bands 1,2,3,4 together for averaging
    idx = np.array([0,10,20,30,40,1,11,21,31,41,2,12,22,32,42,3,13,23,33,43,
                 4,5,6,7,8,9,14,15,16,17,18,19,24,25,26,27,28,29,34,35,36,
                 37,38,39,44,45,46,47,48,49])
    out = out[(idx+idx1)%50,:,:]
    out2 = out2[(idx+idx2)%50,:,:]
    # Put board n+1 P data into last four (P-squared) slots of board n, for convenience 
    out[:,:,4:8] = out2[:,:,0:4]
    # Average repeated channels and place into slots 16-19
    out[19,:,:] = out[15:20,:,:].sum(0)/5.
    out[18,:,:] = out[10:15,:,:].sum(0)/5.
    out[17,:,:] = out[5:10,:,:].sum(0)/5.
    out[16,:,:] = out[0:5,:,:].sum(0)/5.
    # Truncate to slots 16-49 (34 slots)
    out = out[16:,:,:]
    
    if print_attn:
        # At this point, for solar sequence, out[0:34,:,0:8] corresponds to
        # the power in each antenna/feed on these two ROACHes.  Print out the sum
        # in a nice tabular format:
        attn = []
        for i in range(8):
            attn_list = (10*np.log10(out[:,:,i].sum(1)/7.0)).astype('int')/2 * 2
            #bad = np.where(attn_list < 0)[0]
            #if len(bad) > 0:
            #    attn_list[bad] = 0
            # Append list of 34 attenuations
            attn.append(attn_list)
        attn = np.array(attn)
        attn.shape = (8,34)
        return attn
    ovfl = np.array(ovfl)[idx]
    ovfl2 = np.array(ovfl2)[idx]
    ovfl[19] = ovfl[15:20].sum()/5.
    ovfl[18] = ovfl[10:15].sum()/5.
    ovfl[17] = ovfl[5:10].sum()/5.
    ovfl[16] = ovfl[0:5].sum()/5.
    ovfl2[19] = ovfl2[15:20].sum()/5.
    ovfl2[18] = ovfl2[10:15].sum()/5.
    ovfl2[17] = ovfl2[5:10].sum()/5.
    ovfl2[16] = ovfl2[0:5].sum()/5.
    return Time.now(),out,ovfl,ovfl2

def plot_fig(t,out1,out2,ovfl1,ovfl2,ovfl3,ovfl4):
    f, ax = plt.subplots(4, 4, sharex='col', sharey='row')
    f.set_size_inches(10,12,forward=True)
    f.suptitle('Packet Capture at '+t.iso[:19]+' UT',fontsize=18)
    #a1 = str((useroach[0]-1)*2+1)
    #a2 = str((useroach[0]-1)*2+2)
    #a3 = str((useroach[1]-1)*2+1)
    #a4 = str((useroach[1]-1)*2+2)
    titstr = np.array([['Ant 1 X','Ant 5 X','Ant 1 Y','Ant 5 Y'],['Ant 2 X','Ant 6 X','Ant 2 Y','Ant 6 Y'],
              ['Ant 3 X','Ant 7 X','Ant 3 Y','Ant 7 Y'],['Ant 4 X','Ant 8 X','Ant 4 Y','Ant 8 Y']])
    ax[3,0].set_xlabel('Band Number')
    ax[3,1].set_xlabel('Band Number')
    ax[3,2].set_xlabel('Band Number')
    ax[3,3].set_xlabel('Band Number')
    ax[0,0].set_ylabel('Channel')
    ax[1,0].set_ylabel('Channel')
    ax[2,0].set_ylabel('Channel')
    ax[3,0].set_ylabel('Channel')
    ax[0,0].plot(ovfl1[16:]/10,'o',color='white')
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
            attn1 = capture_fig(useroach=[1,2],print_attn=True)
            attn2 = capture_fig(useroach=[3,4],print_attn=True)
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
                    attn[16:,j] = 10
            else:
                # Assume the second argument is only a single number
                attn[0:8,:] = attn1[0:8,:] + int(sys.argv[2])
                attn[8:16,:] = attn2[0:8,:] + int(sys.argv[2])
                attn[16:,:] = 10
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
        t1,out1,ovfl1,ovfl2 = capture_fig(useroach=[1,2])
        t2,out2,ovfl3,ovfl4 = capture_fig(useroach=[3,4])
        plt.close(plot_fig(t1,out1,out2,ovfl1,ovfl2,ovfl3,ovfl4))
    
