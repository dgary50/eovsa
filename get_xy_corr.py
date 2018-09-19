#
# Routines for handling XY delay phase measurements.
# 
# History:
#   2018-Jun-13  DG
#     Started this history log, on the occasion of adding routine xydelay_anal()
#
import numpy as np

def get_xy_corr(npzlist=None, doplot=True, npzlist2=None):
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
    from util import lobe, Time
        
    if doplot: import matplotlib.pylab as plt

    out0 = ri.read_npz([npzlist[0]])  # Parallel scan
    out1 = ri.read_npz([npzlist[1]])  # Perpendicular scan
    if npzlist2 is None:
        pass
    else:
        # Interpret second list as a set of additional files to be concatenated to the first
        for file in npzlist2:
            outx = ri.read_npz([file])
            out0['x'] = np.concatenate((out0['x'],outx['x']),3)
            out1['x'] = np.concatenate((out1['x'],outx['x']),3)
    ph0 = np.angle(np.sum(out0['x'][ri.bl2ord[:13,13]],3))
    ph1 = np.angle(np.sum(out1['x'][ri.bl2ord[:13,13]],3))
    ph0[:,2:] = ph1[:,2:]  # Insert crossed-feed phases from ph1 into ph0

    fghz = out0['fghz']
    nf = len(fghz)
    dph = np.zeros((14,nf),np.float)
    # Determine xi_rot
    xi2 = ph0[:,2] - ph0[:,0] + ph0[:,3] - ph0[:,1]  # This is 2 * xi, measured separately on each of 13 antennas
    xi_rot = np.unwrap(np.angle(np.sum(np.exp(1j*xi2),0)))/2.   # Very clever average does not suffer from wrapping issues
    # Form differential delay phase from channels, and average them
    # dph14 = XY - XX and YY - YX + pi
    #dph14 = np.concatenate((lobe(ph0[:,2] - ph0[:,0] + np.pi/2),lobe(ph0[:,1] - ph0[:,3] - np.pi/2)))  # 26 values for Ant 14
    #dph[13] = np.angle(np.sum(np.exp(1j*dph14),0))  # Very clever average does not suffer from wrapping issues
    # dphi = XX - YX and XY - YY + pi 
    #dphi = np.array((lobe(ph0[:,0] - ph0[:,3] - np.pi/2),lobe(ph0[:,2] - ph0[:,1] + np.pi/2)))  # 2 values for Ant 14
    #dph[:13] = np.angle(np.sum(np.exp(1j*dphi),0))
    # dph14 = XY - XX - xi_rot and YY - YX + xi_rot
    dph14 = np.concatenate((lobe(ph0[:,2] - ph0[:,0] - xi_rot),lobe(ph0[:,1] - ph0[:,3] + xi_rot)))  # 26 values for Ant 14
    dph[13] = np.angle(np.sum(np.exp(1j*dph14),0))  # Very clever average does not suffer from wrapping issues
    # dphi = XX - YX + xi_rot and XY - YY - xi_rot 
    dphi = np.array((lobe(ph0[:,0] - ph0[:,3] + xi_rot),lobe(ph0[:,2] - ph0[:,1] - xi_rot)))  # 2 values for Ant 14
    dph[:13] = np.angle(np.sum(np.exp(1j*dphi),0))
    
    if doplot:
        f, ax = plt.subplots(4, 4, num='XY_Phase')
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
    np.savez('/common/tmp/Feed_rotation/'+npzlist[0].split('/')[-1][:14]+'_delay_phase.npz',fghz=fghz,dph=dph,xi_rot=xi_rot)
    xy_phase = {'timestamp':Time(out0['time'][0],format='jd').lv,'fghz':fghz,'xyphase':dph,'xi_rot':xi_rot}
    return xy_phase

def xydelay_anal(npzfiles):
    ''' Analyze a "standard" X vs. Y delay calibration, consisting of four observations
        on a strong calibrator near 0 HA, in the order:
           90-degree  Low-frequency  receiver,
           90-degree  High-frequency receiver,
            0-degree  High-frequency receiver,
            0-degree  Low-frequency  receiver
    '''
    import matplotlib.pylab as plt
    from util import common_val_idx
    npzfiles = np.array(npzfiles)
    dph_lo = get_xy_corr(npzfiles[[3,0]], doplot=False)
    dph_hi = get_xy_corr(npzfiles[[2,1]])
    ax = plt.figure('XY_Phase').get_axes()
    for i in range(14): 
        ax[i].plot(dph_lo['fghz'],dph_lo['xyphase'][i],'r.')
        ax[i].set_xlim(0,20)
    fghz = np.union1d(dph_lo['fghz'],dph_hi['fghz'])
    nf, = fghz.shape
    flo_uniq = np.setdiff1d(dph_lo['fghz'],dph_hi['fghz'])  # List of frequencies in LO not in HI
    idx_lo_not_hi, idx2 = common_val_idx(fghz, flo_uniq)    # List of indexes of unique LO frequencies
    # Make empty arrays with enough frequencies
    xyphase = np.zeros((14,nf),dtype=float)
    xi_rot = np.zeros((nf),dtype=float)
    idx_hi, idx2 = common_val_idx(fghz,dph_hi['fghz'])  # List of indexes of HI receiver frequencies
    xyphase[:14,idx_hi] = dph_hi['xyphase']       # Insert all high-receiver xyphases
    xyphase[:14,idx_lo_not_hi] = dph_lo['xyphase'][:14,idx_lo_not_hi]  # For unique low-receiver frequencies, insert LO xyphases
    
    xi_rot[idx_hi] = dph_hi['xi_rot']   # Insert all high-receiver xi_rot
    xi_rot[idx_lo_not_hi] = dph_lo['xi_rot'][idx_lo_not_hi]   # For unique low-receiver frequencies, insert LO xi_rot
    dph_hi.update({'xi_rot':xi_rot, 'xyphase':xyphase, 'fghz':fghz})
    print 'Referring to the output of this routine as "xyphase,"'
    print 'run cal_header.xy_phasecal2sql(xyphase) to write the SQL record.' 
    return dph_hi

def apply_xy_corr(out,dph,dphnew=None):
    ''' Does not actually change the data, only calculates and displays it
    '''
    import copy
    import matplotlib.pylab as plt
    from util import lobe
    fghz = out['fghz']
    ph0 = np.angle(np.sum(out['x'][ri.bl2ord[:13,13]],3))
    ph1 = copy.deepcopy(ph0)
    for i in range(13):
        ph1[i,1] += dph[i] - dph[13]
        ph1[i,2] += -dph[13] + np.pi/2
        ph1[i,3] += dph[i] - np.pi/2
    if not dphnew is None:
        ph2 = copy.deepcopy(ph0)
        dph, xi_rot = dphnew
        for i in range(13):
            ph2[i,1] += dph[i] - dph[13]
            ph2[i,2] += -dph[13] - xi_rot[i]
            ph2[i,3] += dph[i] - xi_rot[i] + np.pi
    f, ax = plt.subplots(4,13)
    f.suptitle('Old Correction Applied')
    for i in range(13):
        for j in range(4):
            ax[j,i].cla()
            ax[j,i].plot(fghz,lobe(ph0[i,j]),'.',color='lightgreen')
            ax[j,i].plot(fghz,lobe(ph1[i,j]),'.',color='black')
            ax[j,i].set_ylim(-4,4)
    f, ax = plt.subplots(4,13)
    f.suptitle('New Correction Applied')
    for i in range(13):
        for j in range(4):
            ax[j,i].cla()
            ax[j,i].plot(fghz,lobe(ph0[i,j]),'.',color='lightgreen')
            ax[j,i].plot(fghz,lobe(ph2[i,j]),'.',color='black')
            ax[j,i].set_ylim(-4,4)

def apply_unrot(filename):

    import read_idb as ri
    import dbutil as db
    import copy
    from util import lobe, Time
    import matplotlib.pylab as plt
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
    amp0 = np.abs(np.sum(out['x'][ri.bl2ord[:13,13]],3))
    amp2 = np.abs(np.sum(outp['x'][ri.bl2ord[:13,13]],3))
    f, ax = plt.subplots(4,13)
    for i in range(13):
        for j in range(4):
            ax[j,i].cla()
            ax[j,i].plot(fghz, amp0[i,j],'.',color='lightgreen')
            ax[j,i].plot(fghz, amp2[i,j],'k.')
    ph0 = np.angle(np.sum(out['x'][ri.bl2ord[:13,13]],3))
    ph2 = np.angle(np.sum(outp['x'][ri.bl2ord[:13,13]],3))
    f, ax = plt.subplots(4,13)
    for i in range(13):
        for j in range(4):
            ax[j,i].cla()
            ax[j,i].plot(fghz, ph0[i,j],'.',color='lightgreen')
            ax[j,i].plot(fghz, ph2[i,j],'k.')
            
def apply_unrot_new(filename):

    import read_idb as ri
    import dbutil as db
    import copy
    from util import lobe, Time
    import matplotlib.pylab as plt
    import numpy as np
    blah = np.load('/common/tmp/Feed_rotation/20171223001448_delay_phase.npz')
    dph = blah['dph']
    fghz = blah['fghz']
    xi_rot = blah['xi_rot']
    out = ri.read_npz([filename])
    nbl, npol, nfrq, nt = out['x'].shape
    # Correct data for phase
    #n = [0,0,0,1,1,0,1,0,1,1,0,0,0]
    for i in range(13):
        a1 = lobe(dph[i] - dph[13])
        a2 = -dph[13] - xi_rot
        a3 = dph[i] - xi_rot + np.pi
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
            outp['x'][ri.bl2ord[k,13],0,:,i] = out['x'][ri.bl2ord[k,13],0,:,i]*np.cos(chi[k,i]) + out['x'][ri.bl2ord[k,13],3,:,i]*np.sin(chi[k,i])
            outp['x'][ri.bl2ord[k,13],2,:,i] = out['x'][ri.bl2ord[k,13],2,:,i]*np.cos(chi[k,i]) + out['x'][ri.bl2ord[k,13],1,:,i]*np.sin(chi[k,i])
            outp['x'][ri.bl2ord[k,13],3,:,i] = out['x'][ri.bl2ord[k,13],3,:,i]*np.cos(chi[k,i]) - out['x'][ri.bl2ord[k,13],0,:,i]*np.sin(chi[k,i])
            outp['x'][ri.bl2ord[k,13],1,:,i] = out['x'][ri.bl2ord[k,13],1,:,i]*np.cos(chi[k,i]) - out['x'][ri.bl2ord[k,13],2,:,i]*np.sin(chi[k,i])
    amp0 = np.abs(np.sum(out['x'][ri.bl2ord[:13,13]],3))
    amp2 = np.abs(np.sum(outp['x'][ri.bl2ord[:13,13]],3))
    f, ax = plt.subplots(4,13)
    for i in range(13):
        for j in range(4):
            ax[j,i].cla()
            ax[j,i].plot(fghz, amp0[i,j],'.',color='lightgreen')
            ax[j,i].plot(fghz, amp2[i,j],'k.')
            ax[j,i].set_ylim(0,10)
    ph0 = np.angle(np.sum(out['x'][ri.bl2ord[:13,13]],3))
    ph2 = np.angle(np.sum(outp['x'][ri.bl2ord[:13,13]],3))
    f, ax = plt.subplots(4,13)
    for i in range(13):
        for j in range(4):
            ax[j,i].cla()
            ax[j,i].plot(fghz, ph0[i,j],'.',color='lightgreen')
            ax[j,i].plot(fghz, ph2[i,j],'k.')
            
