import struct, os, time
from pylab import *

def show_image(filename,chan=0, bid=0, norm=False, sort=False, logrange=[None,None]):
    f = open(filename,'rb')
    fsize = os.stat(filename).st_size
    hdr = '<HHHHIHHHHHHHHHHHHHH4I4Q'
    pp2 = '8I8Q8I8Q8I8Q8I8Q8I8Q8I8Q8I8Q8I8Q'
    junk = '320H'
    reclen = struct.calcsize(hdr+pp2+junk)
    nrec = fsize/reclen
    nspec = nrec/256
    print 'This file contains',nspec,'spectra.'
    khdr = ['HeaderLength','PacketNum','FFTShift','AccumLength','GlobalAccumNum',
            'BoardID','AccumNum','DataType','PolType','Ai','Aj','ADCOverflow',
            'QuantClipNum','NSubbands','iFreq','Delay0','Delay1','Delay2','Delay3',
            'PX0','PY0','PX1','PY1','P2X0','P2Y0','P2X1','P2Y1']

    img = np.zeros([4096,50],'float')
    #img2 = np.zeros([4096,nspec],'float')
    accum = np.zeros(50,'int')
    specx = np.zeros(50,'int')
    specy = np.zeros(50,'int')
    # Start reading spectra
    for j in range(nspec):
        n = 0
        # Read records until n=0 (beginning of new accumulation)
        while n != 255:
            # Read an entire record
            try:
                rec = f.read(reclen)
            except:
                print 'End of file'
                break
            # Unpack the header portion and convert to dictionary
            if len(rec) < 88:
                break
            header = struct.unpack(hdr,rec[0:88])
            h = dict(zip(khdr,header))
            dtype = h['DataType']
            bdid = h['BoardID']
            if dtype == 0 and bdid == bid:
                n = h['PacketNum']
                a = h['AccumNum']
                # Unpack power/power-squared
                power = struct.unpack(pp2,rec[88:856])
                # Extract desired channel from this record
                for i in range(8):
                    img[n*16 + i*2,a] = power[i*16 + chan*2]
                    img[n*16 + i*2 + 1,a] = power[i*16 + chan*2 + 1]
        accum[a] = a
        specx[a] = h['PX0']
        specy[a] = h['PY0']

    if chan < 4:
        img = img/(2.**14)
    else:
        img = img/(2.**44)

    if norm:
        bspec = median(img,1)
        for j in range(nspec):
            img[:,a] = img[:,a]/bspec
    for i in range(4):
        # Change 10, 20, 30, ... to 0, and 11, 21, 31, ... to 1, etc. for 0-3
        accum[where(accum % 10 == i)] = i
    idx = accum.argsort()
    #figure()
    #plot(accum,specx,'b.')
    #plot(accum,specy,'r.')
    #yscale('log')
    #title('Blue = X [H], Red = Y [V]')
    figure()
    if logrange[0] is None:
        # Scale to data
        logrange = [log10(img[:4090,:]).min(),log10(img[:4090,:]).max()]
    if sort:
        imp = imshow(log10(img[:,idx]),aspect='auto',vmin=logrange[0],vmax=logrange[1],interpolation='nearest')
    else:
        imp = imshow(log10(img),aspect='auto',vmin=logrange[0],vmax=logrange[1],interpolation='nearest')
    return imp, img

def pspectra(power):

    p = array(power)
    # idx is array([0,1, 16,17, 32,33, 48,49, 64,65, 80,81, 96,97, 112,113])
    idx = array(zip(arange(0,128,16),arange(0,128,16)+1)).flatten()
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

    x = array(xdata)
    # idx is array([0,44,88,132,176,220,264,308,352,396,440,484,528,572,616,660])
    idx = arange(0,704,44)
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

def list_header(filename):
    ''' Reads the header of each record of a file and returns a list of header lines, one
        line per spectrum.  Returns a vector of header lines.
    '''
    f = open(filename,'rb')
    fsize = os.stat(filename).st_size
    hdr = '<HHHHIHHHHHHHHHHHHHH4I4Q'
    pp2 = '8I8Q8I8Q8I8Q8I8Q8I8Q8I8Q8I8Q8I8Q'
    junk = '320H'
    reclen = struct.calcsize(hdr+pp2+junk)
    nrec = fsize/reclen
    nspec = nrec/256
    x = '704H'
    khdr = ['HeaderLength','PacketNum','FFTShift','AccumLength','GlobalAccumNum',
            'BoardID','AccumNum','DataType','PolType','Ai','Aj','ADCOverflow',
            'QuantClipNum','NSubbands','iFreq','Delay0','Delay1','Delay2','Delay3',
            'PX0','PY0','PX1','PY1','P2X0','P2Y0','P2X1','P2Y1']
    lines = []
    badheader = False
    eof = False
    # Start reading records
    for j in range(nspec):
        n = 0
        # Read records until n=0 (beginning of new accumulation)
        while n != 500:
            # Read an entire record
            try:
                rec = f.read(reclen)
            except:
                print 'End of file'
                eof = True
                break
            # Unpack the header portion and convert to dictionary
            try:
                header = struct.unpack(hdr,rec[0:88])
                h = dict(zip(khdr,header))
                a = h['AccumNum']
                n = h['PacketNum']
                #print a, n, h['DataType'],
            except:
                print 'Error reading header.'
                badheader = True
                print 'Expected 88 bytes, but got:',len(rec)
                break
            if n == 0:
                lines.append('{:4d}'.format(h['AccumNum'])
                    +' {:11d}'.format(h['GlobalAccumNum'])
                    +' {:4d}'.format(h['FFTShift'])
                    +' {:2d}'.format(h['BoardID'])
                    +' {:5d}'.format(h['AccumLength'])
                    +' {:5d}'.format(h['ADCOverflow'])
                    +'  {:4d} {:4d}  {:4d} {:4d}'.format(h['PX0'],h['PY0'],h['PX1'],h['PY1'])
                    +'  {:4d} {:4d} {:4d} {:4d}'.format(h['Delay0'],h['Delay1'],h['Delay2'],h['Delay3']))
            if n == 255:
                # Signal that this is the end of PP2 data
                n = 500
        if a == 0:
            # At every 0 Acc#, pop the last line, add a "header" and then add line back
            line = lines.pop()
            lines.append('Acc# Global-Acc# FFTs --M--  OvFl   P1x, P1y   P2x, P2y  ---Delays [nsec]---')
            lines.append(line)
        if badheader:
           break
        if eof:
           break
    return lines

def print_packet_str(filename):
    ''' Prints to the screen a summary of the packet order in the file.
    '''
    f = open(filename,'rb')
    hdr = '<HHHHIHHHHHHHHHHHHHH4I4Q'
    pp2 = '8I8Q8I8Q8I8Q8I8Q8I8Q8I8Q8I8Q8I8Q'
    x = '704h'
    junk = '320H'
    reclen = struct.calcsize(hdr+pp2+junk)
    khdr = ['HeaderLength','PacketNum','FFTShift','AccumLength','GlobalAccumNum',
            'BoardID','AccumNum','DataType','PolType','Ai','Aj','ADCOverflow',
            'QuantClipNum','NSubbands','iFreq','Delay0','Delay1','Delay2','Delay3',
            'PX0','PY0','PX1','PY1','P2X0','P2Y0','P2X1','P2Y1']
    n = 0
    cur_dtype = -1
    cur_bid = -1
    cur_accn = -1
    while 1:
        try:
            rec = f.read(reclen)
            n += 1
        except:
            print 'End of file reached'
            break
        header = struct.unpack(hdr,rec[0:88])
        h = dict(zip(khdr,header))
        this_bid = h['BoardID']
        this_dtype = h['DataType']
        this_accn = h['AccumNum']
        if this_dtype != cur_dtype or this_accn != cur_accn:
            nend = n-1
            if cur_dtype != -1:
                dtype = 'P'
                if cur_dtype:
                    dtype = 'X'
                board = 'R1'
                if cur_bid:
                    board = 'R2'
                print '{:5d}:{:5d} -- '.format(nstart,nend),board,dtype,cur_accn
            nstart = n
            cur_bid = this_bid
            cur_dtype = this_dtype
            cur_accn = this_accn

def show_capture(filename,source):
    f = open(filename,'rb')
    hdr = '<HHHHIHHHHHHHHHHHHHH4I4Q'
    pp2 = '8I8Q8I8Q8I8Q8I8Q8I8Q8I8Q8I8Q8I8Q'
    xfmt = '704h'
    junk = '320H'
    reclen = struct.calcsize(hdr+pp2+junk)
    khdr = ['HeaderLength','PacketNum','FFTShift','AccumLength','GlobalAccumNum',
            'BoardID','AccumNum','DataType','PolType','Ai','Aj','ADCOverflow',
            'QuantClipNum','NSubbands','iFreq','Delay0','Delay1','Delay2','Delay3',
            'PX0','PY0','PX1','PY1','P2X0','P2Y0','P2X1','P2Y1']

    p1ytot = np.zeros(4096)
    if source == 'X':
        # Read only cross-correlation (x) data, skipping any P packets
        # Start reading records
        while 1:
            # Clear spectra to zeros
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

            n = 0
            # Read records until n=0 (beginning of new accumulation)
            while n != 500:
                # Read an entire record
                try:
                    rec = f.read(reclen)
                except:
                    print 'End of file'
                    break
                # Unpack the header portion and convert to dictionary
                header = struct.unpack(hdr,rec[0:88])
                h = dict(zip(khdr,header))
                n = h['PacketNum']
                M = h['AccumLength']
                dtype = h['DataType']
                if dtype:
                    bdid = h['BoardID']
                    if bdid:
                        n = n+128
                    # This is X-engine data
                    # Unpack X data
                    xdata = struct.unpack(xfmt,rec[88:])
                    print n,xdata[0:4]
                    a12x_,a12y_,a13x_,a13y_,a14x_,a14y_,a23x_,a23y_,a24x_,a24y_,a34x_,a34y_ = xspectra(xdata)
                    a12x[n*16:(n+1)*16] = a12x_ 
                    a12y[n*16:(n+1)*16] = a12y_ 
                    a13x[n*16:(n+1)*16] = a13x_ 
                    a13y[n*16:(n+1)*16] = a13y_ 
                    a14x[n*16:(n+1)*16] = a14x_ 
                    a14y[n*16:(n+1)*16] = a14y_ 
                    a23x[n*16:(n+1)*16] = a23x_ 
                    a23y[n*16:(n+1)*16] = a23y_ 
                    a24x[n*16:(n+1)*16] = a24x_ 
                    a24y[n*16:(n+1)*16] = a24y_ 
                    a34x[n*16:(n+1)*16] = a34x_ 
                    a34y[n*16:(n+1)*16] = a34y_ 
                    if n == 255:
                        # Signal that this is the end of 1 accumulation of data
                        n = 500

            # Plot X-engine data
            response = raw_input("Plot What? [a12,a13,a14,a23,a24,a34,p12,etc. else exit]:")
            if response == '':
                response = rsav
            else:
                rsav = response
            if response == 'a12':
                clf()
                plot(abs(a12x),'r')
                plot(abs(a12y),'b')
                ylabel('Amplitude')
                title('Ant12 Amp: X-Chan (Red); Y-Chan (Blu)')
            elif response == 'p12':
                clf()
                plot(angle(a12x),'r.')
                plot(angle(a12y),'b.')
                ylabel('Phase')
                title('Ant12 Phase X-Chan (Red); Y-Chan (Blu)')
            elif response == 'a13':
                clf()
                plot(abs(a13x),'r')
                plot(abs(a13y),'b')
                ylabel('Amplitude')
                title('Ant13 Amp: X-Chan (Red); Y-Chan (Blu)')
            elif response == 'p13':
                clf()
                plot(angle(a13x),'r.')
                plot(angle(a13y),'b.')
                ylabel('Phase')
                title('Ant13 Phase: X-Chan (Red); Y-Chan (Blu)')
            elif response == 'a14':
                clf()
                plot(abs(a14x),'r')
                plot(abs(a14y),'b')
                ylabel('Amplitude')
                title('Ant14 Amp: X-Chan (Red); Y-Chan (Blu)')
            elif response == 'p14':
                clf()
                plot(angle(a14x),'r.')
                plot(angle(a14y),'b.')
                ylabel('Phase')
                title('Ant14 Phase: X-Chan (Red); Y-Chan (Blu)')
            elif response == 'a23':
                clf()
                plot(abs(a23x),'r')
                plot(abs(a23y),'b')
                ylabel('Amplitude')
                title('Ant23 Amp: X-Chan (Red); Y-Chan (Blu)')
            elif response == 'p23':
                clf()
                plot(angle(a23x),'r.')
                plot(angle(a23y),'b.')
                ylabel('Phase')
                title('Ant23 Phase: X-Chan (Red); Y-Chan (Blu)')
            elif response == 'a24':
                clf()
                plot(abs(a24x),'r')
                plot(abs(a24y),'b')
                ylabel('Amplitude')
                title('Ant24 Amp: X-Chan (Red); Y-Chan (Blu)')
            elif response == 'p24':
                clf()
                plot(angle(a24x),'r.')
                plot(angle(a24y),'b.')
                ylabel('Phase')
                title('Ant24 Phase: X-Chan (Red); Y-Chan (Blu)')
            elif response == 'a34':
                clf()
                plot(abs(a34x),'r')
                plot(abs(a34y),'b')
                ylabel('Amplitude')
                title('Ant34 Amp: X-Chan (Red); Y-Chan (Blu)')
            elif response == 'p34':
                clf()
                plot(angle(a34x),'r.')
                plot(angle(a34y),'b.')
                ylabel('Phase')
                title('Ant34 Phase: X-Chan (Red); Y-Chan (Blu)')
            else:
                break
            if response[0] == 'p':
                yscale('linear')
            else:
                yscale('log')
            xlabel('Accum '+str(h['AccumNum'])+' Channel Number')

    else:
        # Read only P data, skipping any X packets
        while 1:
            p1x = np.zeros(4096,'float')
            p1y = np.zeros(4096,'float')
            p2x = np.zeros(4096,'float')
            p2y = np.zeros(4096,'float')
            p3x = np.zeros(4096,'float')
            p3y = np.zeros(4096,'float')
            p4x = np.zeros(4096,'float')
            p4y = np.zeros(4096,'float')
            P1x = np.zeros(4096,'float')
            P1y = np.zeros(4096,'float')
            P2x = np.zeros(4096,'float')
            P2y = np.zeros(4096,'float')
            P3x = np.zeros(4096,'float')
            P3y = np.zeros(4096,'float')
            P4x = np.zeros(4096,'float')
            P4y = np.zeros(4096,'float')

            n = 0
            # Read records until n=0 (beginning of new accumulation)
            while n != 500:
                # Read an entire record
                try:
                    rec = f.read(reclen)
                except:
                    print 'End of file'
                    break
                # Unpack the header portion and convert to dictionary
                header = struct.unpack(hdr,rec[0:88])
                h = dict(zip(khdr,header))
                n = h['PacketNum']
                M = h['AccumLength']
                dtype = h['DataType']
                if dtype:
                    pass
                else:
                    # This is P-P2 data
                    # Unpack power/power-squared
                    power = struct.unpack(pp2,rec[88:856])
                    p1x_,p1y_,p2x_,p2y_,P1x_,P1y_,P2x_,P2y_ = pspectra(power)
                    bdid = h['BoardID']
                    if bdid:
                        # This is from ROACH2
                        p3x[n*16:(n+1)*16] = p1x_ 
                        p4x[n*16:(n+1)*16] = p2x_ 
                        p3y[n*16:(n+1)*16] = p1y_ 
                        p4y[n*16:(n+1)*16] = p2y_ 
                        P3x[n*16:(n+1)*16] = P1x_ 
                        P4y[n*16:(n+1)*16] = P1y_ 
                        P3x[n*16:(n+1)*16] = P2x_ 
                        P4y[n*16:(n+1)*16] = P2y_
                    else:
                        # This is from ROACH1
                        p1x[n*16:(n+1)*16] = p1x_ 
                        p2x[n*16:(n+1)*16] = p2x_ 
                        p1y[n*16:(n+1)*16] = p1y_ 
                        p2y[n*16:(n+1)*16] = p2y_ 
                        P1x[n*16:(n+1)*16] = P1x_ 
                        P1y[n*16:(n+1)*16] = P1y_ 
                        P2x[n*16:(n+1)*16] = P2x_ 
                        P2y[n*16:(n+1)*16] = P2y_
                    if bdid and n == 255:
                        # Signal that this is the end of 1 accumulation of data
                        n = 500

            response = raw_input("Plot What? [px,py,p1,p2,p3,p4 for power, sk*, P* for power^2, else exit]:")
            if response == '':
                response = rsav
            else:
                rsav = response
            if response == 'px':
                clf()
                plot(p1x,'r')
                plot(p2x,color='orange')
                plot(p3x,'g')
                plot(p4x,'b')
                ylabel('Power')
                title('X-Chan Power: Ant1 (Red); Ant2 (Org); Ant3 (Grn); Ant4 (Blu)')
            elif response == 'py':
                clf()
                plot(p1y,'r')
                plot(p2y,color='orange')
                plot(p3y,'g')
                plot(p4y,'b')
                ylabel('Power')
                title('Y-Chan Power: Ant1 (Red); Ant2 (Org); Ant3 (Grn); Ant4 (Blu)')
            elif response == 'p1':
                clf()
                plot(p1x,'r')
                plot(p1y,'g')
                ylabel('Power')
                title('Antenna 1 X-Chan Power (Red); Y-Chan Power (Grn)')
            elif response == 'p2':
                clf()
                plot(p2x,'r')
                plot(p2y,'g')
                ylabel('Power')
                title('Antenna 2 X-Chan Power (Red); Y-Chan Power (Grn)')
            elif response == 'p3':
                clf()
                plot(p3x,'r')
                plot(p3y,'g')
                ylabel('Power')
                title('Antenna 3 X-Chan Power (Red); Y-Chan Power (Grn)')
            elif response == 'p4':
                clf()
                plot(p4x,'r')
                plot(p4y,'g')
                ylabel('Power')
                title('Antenna 4 X-Chan Power (Red); Y-Chan Power (Grn)')
            elif response == 'Px':
                clf()
                plot(P1x,'r')
                plot(P2x,color='orange')
                plot(P3x,'g')
                plot(P4x,'b')
                ylabel('Power Squared')
                title('X-Chan Power Squared: Ant1 (Red); Ant2 (Org); Ant3 (Grn); Ant4 (Blu)')
            elif response == 'Py':
                clf()
                plot(P1y,'r')
                plot(P2y,color='orange')
                plot(P3y,'g')
                plot(P4y,'b')
                ylabel('Power Squared')
                title('Y-Chan Power Squared: Ant1 (Red); Ant2 (Org); Ant3 (Grn); Ant4 (Blu)')
            elif response == 'P1':
                clf()
                plot(P1x,'r')
                plot(P1y,'g')
                ylabel('Power Squared')
                title('Antenna 1 X-Chan Power Squared (Red); Y-Chan Power Squared (Grn)')
            elif response == 'P2':
                clf()
                plot(P2x,'r')
                plot(P2y,'g')
                ylabel('Power Squared')
                title('Antenna 2 X-Chan Power Squared (Red); Y-Chan Power Squared (Grn)')
            elif response == 'P3':
                clf()
                plot(P3x,'r')
                plot(P3y,'g')
                ylabel('Power Squared')
                title('Antenna 3 X-Chan Power Squared (Red); Y-Chan Power Squared (Grn)')
            elif response == 'P4':
                clf()
                plot(P4x,'r')
                plot(P4y,'g')
                ylabel('Power')
                title('Antenna 4 X-Chan Power (Red); Y-Chan Power (Grn)')
            elif response == 'sk1':
                clf()
                plot(((M+1.)/(M-1.))*(M*P1x/p1x**2 - 1),'r')
                plot(((M+1.)/(M-1.))*(M*P1y/p1y**2 - 1),'g')
                ylabel('SK Ant 1')
                title('Antenna 1 SK; X-Channel (Red); Y-Channel (Grn)')
            elif response == 'sk2':
                clf()
                plot(((M+1.)/(M-1.))*(M*P2x/p2x**2 - 1),'r')
                plot(((M+1.)/(M-1.))*(M*P2y/p2y**2 - 1),'g')
                ylabel('SK Ant 2')
                title('Antenna 2 SK; X-Channel (Red); Y-Channel (Grn)')
            elif response == 'sk3':
                clf()
                plot(((M+1.)/(M-1.))*(M*P3x/p3x**2 - 1),'r')
                plot(((M+1.)/(M-1.))*(M*P3y/p3y**2 - 1),'g')
                ylabel('SK Ant 3')
                title('Antenna 3 SK; X-Channel (Red); Y-Channel (Grn)')
            elif response == 'sk4':
                clf()
                plot(((M+1.)/(M-1.))*(M*P4x/p4x**2 - 1),'r')
                plot(((M+1.)/(M-1.))*(M*P4y/p4y**2 - 1),'g')
                ylabel('SK Ant 4')
                title('Antenna 4 SK; X-Channel (Red); Y-Channel (Grn)')
            elif response == 'skx':
                clf()
                plot(((M+1.)/(M-1.))*(M*P1x/p1x**2 - 1),'r')
                plot(((M+1.)/(M-1.))*(M*P2x/p2x**2 - 1),color='orange')
                plot(((M+1.)/(M-1.))*(M*P3x/p3x**2 - 1),'g')
                plot(((M+1.)/(M-1.))*(M*P4x/p4x**2 - 1),'b')
                ylabel('X-Channel SK')
                title('X-Channel SK: Ant1 (Red); Ant2 (Org); Ant3 (Grn); Ant4 (Blu)')
            elif response == 'sky':
                clf()
                plot(((M+1.)/(M-1.))*(M*P1y/p1y**2 - 1),'r')
                plot(((M+1.)/(M-1.))*(M*P2y/p2y**2 - 1),color='orange')
                plot(((M+1.)/(M-1.))*(M*P3y/p3y**2 - 1),'g')
                plot(((M+1.)/(M-1.))*(M*P4y/p4y**2 - 1),'b')
                ylabel('Y-Channel SK')
                title('Y-Channel SK: Ant1 (Red); Ant2 (Org); Ant3 (Grn); Ant4 (Blu)')
            else:
                break
            yscale('log')
            xlabel('Accum '+str(h['AccumNum'])+' Channel Number')
            
