def get_xy_corr(npzlist=None, doplot=True):
    ''' Analyze a pair of parallel and cross polarization calibration scans and
        return the X vs. Y delay phase corrections on all antennas 1-14.
        
        Required keyword:
           npzlist   a list of 2 NPZ filenames, the first being the 
                       parallel-feed scan, and the second being the
                       crossed-feed scan.
        Optional keyword:
           doplot    True => plot the final result, False => no plot
    '''
    if npzlist is None:
        print 'Must provide a list of 2 NPZ files.'
        return None, None
    import read_idb as ri
    import numpy as np
    from util import lobe
        
    if doplot: import matplotlib.pylab as plt

    out0 = ri.read_npz([npzlist[0]])  # Parallel scan
    out1 = ri.read_npz([npzlist[1]])  # Perpendicular scan
    ph0 = np.angle(np.sum(out0['x'][ri.bl2ord[:13,13]],3))
    ph1 = np.angle(np.sum(out1['x'][ri.bl2ord[:13,13]],3))
    ph0[:,2:] = ph1[:,2:]  # Insert crossed-feed phases from ph1 into ph0

    fghz = out0['fghz']
    nf = len(fghz)
    dph = np.zeros((14,nf),np.float)
    # Form differential delay phase from channels, and average them
    # dph14 = XY - XX and YY - YX + pi
    dph14 = np.concatenate((lobe(ph0[:,2] - ph0[:,0] + np.pi/2),lobe(ph0[:,1] - ph0[:,3] - np.pi/2)))  # 26 values for Ant 14
    dph[13] = np.angle(np.sum(np.exp(1j*dph14),0))  # Very clever average does not suffer from wrapping issues
    # dphi = XX - YX and XY - YY + pi 
    dphi = np.array((lobe(ph0[:,0] - ph0[:,3] - np.pi/2),lobe(ph0[:,2] - ph0[:,1] + np.pi/2)))  # 2 values for Ant 14
    dph[:13] = np.angle(np.sum(np.exp(1j*dphi),0))
    
    if doplot:
        f, ax = plt.subplots(4,4)
        ax.shape = (16,)
        for i in range(13): 
            ax[i].plot(fghz,dphi[0,i],'.')
            ax[i].plot(fghz,dphi[1,i],'.')
            ax[i].plot(fghz,dph[i],'k.')
        for i in range(26):
            ax[13].plot(fghz,dph14[i],'.')
        ax[13].plot(fghz,dph[13],'k.')
        for i in range(14): ax[i].set_ylim(-4,4)
        f.suptitle('Multicolor: Measurements, Black: Final Results')
    np.savez('/common/tmp/Feed_rotation/'+npzlist[0].split('/')[-1][:14]+'_delay_phase.npz',fghz=fghz,dph=dph)
    return dph

def apply_xy_corr(out,dph):
    ''' Does not actually change the data, only calculates and displays it
    '''
    import copy
    import matplotlib.pylab as plt
    ph0 = np.angle(np.sum(out['x'][ri.bl2ord[:13,13]],3))
    ph1 = copy.deepcopy(ph0)
    for i in range(13):
        ph1[i,1] += dph[i] - dph[13]
        ph1[i,2] += -dph[13] + np.pi/2
        ph1[i,3] += dph[i] + np.pi/2
    f, ax = plt.subplots(4,13)
    for i in range(13):
        for j in range(4):
            ax[j,i].cla()
            ax[j,i].plot(fghz,lobe(ph0[i,j]),'.',color='lightgreen')
            ax[j,i].plot(fghz,lobe(ph1[i,j]),'.',color='black')
            ax[j,i].set_ylim(-4,4)

def apply_unrot(filename):

    import read_idb as ri
    import dbutil as db
    import copy
    from util import lobe, Time
    import matplotlib.pylab as plt
    import numpy as np
    blah = np.load('/common/tmp/Feed_rotation/20170702121949_delay_phase.npz')
    dph = blah['dph']
    fghz = blah['fghz']
    out = ri.read_npz([filename])
    nbl, npol, nfrq, nt = out['x'].shape
    # Correct data for phase
    #n = [0,0,0,1,1,0,1,0,1,1,0,0,0]
    for i in range(13):
        a1 = lobe(dph[i] - dph[13])
        a2 = -dph[13] + np.pi/2
        a3 = dph[i] - np.pi/2
        for j in range(nt):
            out['x'][ri.bl2ord[i,13],1,:,j] *= np.exp(1j*a1)
            out['x'][ri.bl2ord[i,13],2,:,j] *= np.exp(1j*a2) 
            out['x'][ri.bl2ord[i,13],3,:,j] *= np.exp(1j*a3)
        
    trange = Time(out['time'][[0,-1]],format='jd')
    times, chi = db.get_chi(trange)
    nskip = len(times)/nt
    chi = np.transpose(chi[::nskip+1])
    chi[[8,9,10,12]] = 0.0
    outp = copy.deepcopy(out)
    for i in range(nt):
        for k in range(13):
            outp['x'][ri.bl2ord[k,13],0] = out['x'][ri.bl2ord[k,13],0]*np.cos(chi[k,i]) + out['x'][ri.bl2ord[k,13],3]*np.sin(chi[k,i])
            outp['x'][ri.bl2ord[k,13],2] = out['x'][ri.bl2ord[k,13],2]*np.cos(chi[k,i]) + out['x'][ri.bl2ord[k,13],1]*np.sin(chi[k,i])
            outp['x'][ri.bl2ord[k,13],3] = out['x'][ri.bl2ord[k,13],3]*np.cos(chi[k,i]) - out['x'][ri.bl2ord[k,13],0]*np.sin(chi[k,i])
            outp['x'][ri.bl2ord[k,13],1] = out['x'][ri.bl2ord[k,13],1]*np.cos(chi[k,i]) - out['x'][ri.bl2ord[k,13],2]*np.sin(chi[k,i])
    amp0 = np.abs(np.sum(out['x'][ri.bl2ord[:,13]],3))
    amp2 = np.abs(np.sum(outp['x'][ri.bl2ord[:,13]],3))
    f, ax = plt.subplots(4,13)
    for i in range(13):
        for j in range(4):
            ax[j,i].cla()
            ax[j,i].plot(fghz, amp0[i,j],'.',color='lightgreen')
            ax[j,i].plot(fghz, amp2[i,j],'k.')
    ph0 = np.angle(np.sum(out['x'][ri.bl2ord[:,13]],3))
    ph2 = np.angle(np.sum(outp['x'][ri.bl2ord[:,13]],3))
    f, ax = plt.subplots(4,13)
    for i in range(13):
        for j in range(4):
            ax[j,i].cla()
            ax[j,i].plot(fghz, ph0[i,j],'.',color='lightgreen')
            ax[j,i].plot(fghz, ph2[i,j],'k.')
            