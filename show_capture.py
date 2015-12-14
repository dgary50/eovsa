import struct, os, time
from pylab import *

def show_image(filename,chan=0,norm=False, sort=False, logrange=[None,None]):
    f = open(filename,'rb')
    fsize = os.stat(filename).st_size
    hdr = 'HHHHIHHHHHHHHHHHHHH4I4Q'
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

    img = np.zeros([4096,nspec],'float')
    img2 = np.zeros([4096,nspec],'float')
    accum = np.zeros(nspec,'int')
    specx = np.zeros(nspec,'int')
    specy = np.zeros(nspec,'int')
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
            header = struct.unpack(hdr,rec[0:88])
            h = dict(zip(khdr,header))
            n = h['PacketNum']
            a = h['AccumNum']
            # Unpack power/power-squared
            power = struct.unpack(pp2,rec[88:856])

            # Extract desired channel from this record
            for i in range(8):
                img[n*16 + i*2,j] = power[i*16 + chan*2]
                img[n*16 + i*2 + 1,j] = power[i*16 + chan*2 + 1]
        accum[j] = a
        specx[j] = h['PX0']
        specy[j] = h['PY0']

    if chan < 4:
        img = img/(2.**14)
    else:
        img = img/(2.**44)

    if norm:
        bspec = median(img,1)
        for j in range(nspec):
            img[:,j] = img[:,j]/bspec
    for i in range(4):
        # Change 10, 20, 30, ... to 0, and 11, 21, 31, ... to 1, etc. for 0-3
        accum[where(accum % 10 == i)] = i
    idx = accum.argsort()
    figure()
    plot(accum,specx,'b.')
    plot(accum,specy,'r.')
    yscale('log')
    title('Blue = X [H], Red = Y [V]')
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

    p1x = np.zeros(16,'float')
    p1y = np.zeros(16,'float')
    p2x = np.zeros(16,'float')
    p2y = np.zeros(16,'float')
    P1x = np.zeros(16,'float')
    P1y = np.zeros(16,'float')
    P2x = np.zeros(16,'float')
    P2y = np.zeros(16,'float')
    # Distribute P, P2 into their spectra
    for i in range(8):
        p1x[i*2] = power[i*16]
        p1x[i*2 + 1] = power[i*16 + 1]
        p1y[i*2] = power[i*16 + 2]
        p1y[i*2 + 1] = power[i*16 + 3]
        p2x[i*2] = power[i*16 + 4]
        p2x[i*2 + 1] = power[i*16 + 5]
        p2y[i*2] = power[i*16 + 6]
        p2y[i*2 + 1] = power[i*16 + 7]
        P1x[i*2] = power[i*16 + 8]
        P1x[i*2 + 1] = power[i*16 + 9]
        P1y[i*2] = power[i*16 + 10]
        P1y[i*2 + 1] = power[i*16 + 11] 
        P2x[i*2] = power[i*16 + 12]
        P2x[i*2 + 1] = power[i*16 + 13]
        P2y[i*2] = power[i*16 + 14]
        P2y[i*2 + 1] = power[i*16 + 15]

    # We now have one accumulation, so make a plot (after scaling)
    p1x = p1x/(2.**14)
    p1y = p1y/(2.**14)
    p2x = p2x/(2.**14)
    p2y = p2y/(2.**14)
    P1x = P1x/(2.**44)
    P1y = P1y/(2.**44)
    P2x = P2x/(2.**44)
    P2y = P2y/(2.**44)

    return p1x,p1y,p2x,p2y,P1x,P1y,P2x,P2y

def show_capture(filename,x=False):
    f = open(filename,'rb')
    hdr = 'HHHHIHHHHHHHHHHHHHH4I4Q'
    pp2 = '8I8Q8I8Q8I8Q8I8Q8I8Q8I8Q8I8Q8I8Q'
    x = '704H'
    junk = '320H'
    reclen = struct.calcsize(hdr+pp2+junk)
    khdr = ['HeaderLength','PacketNum','FFTShift','AccumLength','GlobalAccumNum',
            'BoardID','AccumNum','DataType','PolType','Ai','Aj','ADCOverflow',
            'QuantClipNum','NSubbands','iFreq','Delay0','Delay1','Delay2','Delay3',
            'PX0','PY0','PX1','PY1','P2X0','P2Y0','P2X1','P2Y1']

    p1ytot = np.zeros(4096)
    # Start reading records
    while 1:
        # Clear spectra to zeros
        p1x = np.zeros(4096,'float')
        p1y = np.zeros(4096,'float')
        p2x = np.zeros(4096,'float')
        p2y = np.zeros(4096,'float')
        P1x = np.zeros(4096,'float')
        P1y = np.zeros(4096,'float')
        P2x = np.zeros(4096,'float')
        P2y = np.zeros(4096,'float')

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
                # This is X-engine data
                # Unpack X data
                xdata = struct.unpack(x,rec[88:1496])
                print h
            else:
                # This is P-P2 data
                # Unpack power/power-squared
                power = struct.unpack(pp2,rec[88:856])
                p1x_,p1y_,p2x_,p2y_,P1x_,P1y_,P2x_,P2y_ = pspectra(power)
                p1x[n*16:(n+1)*16] = p1x_ 
                p2x[n*16:(n+1)*16] = p2x_ 
                p1y[n*16:(n+1)*16] = p1y_ 
                p2y[n*16:(n+1)*16] = p2y_ 
                P1x[n*16:(n+1)*16] = P1x_ 
                P1y[n*16:(n+1)*16] = P1y_ 
                P2x[n*16:(n+1)*16] = P2x_ 
                P2y[n*16:(n+1)*16] = P2y_ 
                if n == 255:
                    # Signal that this is the end of PP2 data
                    n = 500

        if dtype:
            # Plot X-engine data
            pass
        else:
            response = raw_input("Plot What? [px,py,p1,p2 for power, sk*, P* for power^2, else exit]:")
            if response == '':
                response = rsav
            else:
                rsav = response
            if response == 'px':
                clf()
                plot(p1x,'r')
                plot(p2x,'b')
#                plot(p1x,'r.')
#                plot(p2x,'b.')
                ylabel('Power')
                title('Ant '+str(h['Ai'])+' X-Chan Power (Red); Ant '+str(h['Aj'])+' X-Chan Power (Blu)')
            elif response == 'py':
                clf()
                plot(p1y,'g')
                plot(p2y,color='orange')
#                plot(p1y,'g.')
#                plot(p2y,'.',color='orange')
                ylabel('Power')
                title('Ant '+str(h['Ai'])+' Y-Chan Power (Grn); Ant '+str(h['Aj'])+' Y-Chan Power (Org)')
            elif response == 'p1':
                clf()
                plot(p1x,'r')
                plot(p1y,'g')
#                plot(p1x,'r.')
#                plot(p1y,'g.')
                ylabel('Power')
                title('Ant '+str(h['Ai'])+' X-Chan Power (Red); Ant '+str(h['Ai'])+' Y-Chan Power (Grn)')
#                if h['AccumNum'] > 14 and h['AccumNum'] < 20:
#                    p1ytot += p1x
#                if h['AccumNum'] == 20:
#                    fout = open(filename[22:],'wb')
#                    fout.write(p1ytot)
#                    fout.close()
#                    p1ytot = np.zeros(4096)
#                    time.sleep(1)
#                    fin = open(filename[22:],'rb')
#                    buf = fin.read()
#                    fin.close()
#                    new = struct.unpack('4096d',buf)
#                    figure()
#                    plot(new,'b')
                    
            elif response == 'p2':
                clf()
                plot(p2x,'b')
                plot(p2y,color='orange')
#                plot(p2x,'b.')
#                plot(p2y,'.',color='orange')
                ylabel('Power')
                title('Ant '+str(h['Aj'])+' X-Chan Power (Blu); Ant '+str(h['Aj'])+' Y-Chan Power (Org)')
            elif response == 'Px':
                clf()
                plot(P1x,'r')
                plot(P2x,'b')
#                plot(P1x,'r.')
#                plot(P2x,'b.')
                ylabel('Power^2')
                title('Ant '+str(h['Ai'])+' X-Chan Power^2 (Red); Ant '+str(h['Aj'])+' X-Chan Power^2 (Blu)')
            elif response == 'Py':
                clf()
                plot(P1y,'g')
                plot(P2y,color='orange')
#                plot(P1y,'g.')
#                plot(P2y,'.',color='orange')
                ylabel('Power^2')
                title('Ant '+str(h['Ai'])+' Y-Chan Power^2 (Grn); Ant '+str(h['Aj'])+' Y-Chan Power^2 (Org)')
            elif response == 'P1':
                clf()
                plot(P1x,'r')
                plot(P1y,'b')
#                plot(P1x,'r.')
#                plot(P1y,'b.')
                ylabel('Power^2')
                title('Ant '+str(h['Ai'])+' X-Chan Power^2 (Red); Ant '+str(h['Ai'])+' Y-Chan Power^2 (Blu')
            elif response == 'P2':
                clf()
                plot(P2x,'g')
                plot(P2y,color='orange')
#                plot(P2x,'g.')
#                plot(P2y,'.',color='orange')
                ylabel('Power^2')
                title('Ant '+str(h['Aj'])+' X-Chan Power^2 (Grn); Ant '+str(h['Aj'])+' Y-Chan Power^2 (Org)')
            elif response == 'sk1':
                clf()
                plot(((M+1.)/(M-1.))*(M*P1x/p1x**2 - 1))
                plot(((M+1.)/(M-1.))*(M*P1y/p1y**2 - 1))
                ylabel('SK Ant 1')
                title('Ant '+str(h['Ai'])+' X-Channel SK; Ant '+str(h['Ai'])+' Y-Channel SK')
            elif response == 'sk2':
                clf()
                plot(((M+1.)/(M-1.))*(M*P2x/p2x**2 - 1))
                plot(((M+1.)/(M-1.))*(M*P2y/p2y**2 - 1))
                ylabel('SK Ant 1')
                title('Ant '+str(h['Aj'])+' X-Channel SK; Ant '+str(h['Aj'])+' Y-Channel SK')
            elif response == 'skx':
                clf()
                plot(((M+1.)/(M-1.))*(M*P1x/p1x**2 - 1))
                plot(((M+1.)/(M-1.))*(M*P2x/p2x**2 - 1))
                ylabel('SK Ant 1')
                title('Ant '+str(h['Ai'])+' X-Channel SK; Ant '+str(h['Aj'])+' X-Channel SK')
            elif response == 'sky':
                clf()
                plot(((M+1.)/(M-1.))*(M*P1y/p1y**2 - 1))
                plot(((M+1.)/(M-1.))*(M*P2y/p2y**2 - 1))
                ylabel('SK Ant 1')
                title('Ant '+str(h['Ai'])+' Y-Channel SK; Ant '+str(h['Aj'])+' Y-Channel SK')
            else:
                break

            if response[0:2] == 'sk':
                # These are relaxed limits to flag fewer points, to account for union of flags
                # across all antennas.
                plot([0,4096],[0.84,0.84],'r')
                plot([0,4096],[1.2,1.2],'r')
                ylim([0.1,10])
            yscale('log')
            xlabel('Channel Number')
            print 'Accum#',h['AccumNum'],', FFTs',h['FFTShift'],', M',h['AccumLength'],', ADCOvr',h['ADCOverflow'],', P1x,y',h['PX0'],h['PY0'],', Delays',h['Delay0'],h['Delay1'],h['Delay2'],h['Delay3']
