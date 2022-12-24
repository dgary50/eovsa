#
# Routines for handling XY delay phase measurements.
# 
# History:
#   2018-Jun-13  DG
#     Started this history log, on the occasion of adding routine xydelay_anal()
#   2020-Jun-21  DG
#     Updated xydelay_anal() to adjust the low-frequency receiver xyphase and xi_rot
#     to make them consistent with the high-frequency receiver by correcting for the
#     pi ambiguity.
#   2020-Aug-09  DG
#     More adjustments to xydelay_anal() for making low- and high-frequency receiver 
#     corrections agree.  To fix a glitch in the low-frequency receiver observations
#     that was sometimes occurring, I changed xydelay_anal() to read the npzfiles and
#     correct the offending scan, and then changed get_xy_corr() accordingly to take
#     the already read-in data as input.
#
import numpy as np
from util import lobe, Time, ant_str2list
import read_idb as ri

def get_xy_corr(out, ant_str='ant1-13',doplot=True):
    ''' Analyze a pair of parallel and cross polarization calibration scans and
        return the X vs. Y delay phase corrections on all antennas 1-14.
        
        Required keyword:
           out   a 2-element array of dicts representing two scans, the first 
                       being the parallel-feed scan, and the second being the
                       crossed-feed scan.
        Optional keyword:
           doplot    True => plot the final result, False => no plot
    '''
    if doplot: import matplotlib.pylab as plt

    tstr = Time(out[0]['time'][0],format='jd').iso[:19].replace('-','').replace(' ','').replace(':','')
    ph0 = np.angle(np.sum(out[0]['x'][ri.bl2ord[:13,13]],3))
    ph1 = np.angle(np.sum(out[1]['x'][ri.bl2ord[:13,13]],3))
    ph0[:,2:] = ph1[:,2:]  # Insert crossed-feed phases from ph1 into ph0
    #import pdb; pdb.set_trace()

    fghz = out[0]['fghz']
    nf = len(fghz)
    dph = np.zeros((14,nf),np.float)
    dphi = np.zeros((2,13,nf),np.float)
    # Determine xi_rot
    xi2 = ph0[:,2] - ph0[:,0] + ph0[:,3] - ph0[:,1]  # This is 2 * xi, measured separately on each of 13 antennas
    antlist = ant_str2list(ant_str)
    xi_rot = np.unwrap(np.angle(np.sum(np.exp(1j*xi2[antlist]),0)))/2.   # Very clever average does not suffer from wrapping issues
    # Form differential delay phase from channels, and average them
    # dph14 = XY - XX and YY - YX + pi
    #dph14 = np.concatenate((lobe(ph0[:,2] - ph0[:,0] + np.pi/2),lobe(ph0[:,1] - ph0[:,3] - np.pi/2)))  # 26 values for Ant 14
    #dph[13] = np.angle(np.sum(np.exp(1j*dph14),0))  # Very clever average does not suffer from wrapping issues
    # dphi = XX - YX and XY - YY + pi 
    #dphi = np.array((lobe(ph0[:,0] - ph0[:,3] - np.pi/2),lobe(ph0[:,2] - ph0[:,1] + np.pi/2)))  # 2 values for Ant 14
    #dph[:13] = np.angle(np.sum(np.exp(1j*dphi),0))
    # dph14 = XY - XX - xi_rot and YY - YX + xi_rot
    dph14 = np.concatenate((lobe(ph0[antlist,2] - ph0[antlist,0] - xi_rot),lobe(ph0[antlist,1] - ph0[antlist,3] + xi_rot)))  # 26 values for Ant 14
    dph[13] = np.angle(np.sum(np.exp(1j*dph14),0))  # Very clever average does not suffer from wrapping issues
    # dphi = XX - YX + xi_rot and XY - YY - xi_rot 
    dphi[:,antlist] = np.array((lobe(ph0[antlist,0] - ph0[antlist,3] + xi_rot),lobe(ph0[antlist,2] - ph0[antlist,1] - xi_rot)))  # 2 values for Ant 14
    dph[antlist] = np.angle(np.sum(np.exp(1j*dphi[:,antlist]),0))
    
    if doplot:
        f, ax = plt.subplots(4, 4, num='XY_Phase')
        ax.shape = (16,)
        for i in range(13): 
            ax[i].plot(fghz,dphi[0,i],'.')
            ax[i].plot(fghz,dphi[1,i],'.')
            ax[i].plot(fghz,dph[i],'k.')
        for i in range(len(dph14[:,0])):
            ax[13].plot(fghz,dph14[i],'.')
        ax[13].plot(fghz,dph[13],'k.')
        for i in range(14): ax[i].set_ylim(-4,4)
        f.suptitle('Multicolor: Measurements, Black: Final Results')
    np.savez('/common/tmp/Feed_rotation/'+tstr+'_delay_phase.npz',fghz=fghz,dph=dph,xi_rot=xi_rot)
    xy_phase = {'timestamp':Time(out[0]['time'][0],format='jd').lv,'fghz':fghz,'xyphase':dph,'xi_rot':xi_rot, 'dphi':dphi, 'dph14':dph14}
    return xy_phase

def xydelay_anal(npzfiles, ant_str='ant1-13', fix_tau_lo=None):
    ''' Analyze a "standard" X vs. Y delay calibration, consisting of four observations
        on a strong calibrator near 0 HA, in the order:
           90-degree  Low-frequency  receiver,
           90-degree  High-frequency receiver,
            0-degree  High-frequency receiver,
            0-degree  Low-frequency  receiver
            
        It has happened that the low-frequency receiver delays were not set for one
        of the observation sets.  This can be fixed by reading the delay-center table
        for the relevant date and applying a phase correction according to the delay
        difference.  Setting fix_tau_lo to "first" means correct first scan, and to
        "last" means correct last scan.
    '''
    import matplotlib.pylab as plt
    from util import common_val_idx

    npzfiles = np.array(npzfiles)
    out = []
    for file in npzfiles:
        out.append(ri.read_npz([file]))
    out = np.array(out)
    if fix_tau_lo != None:
        # Correct for low-frequency delay error, if requested
        import cal_header as ch
        from stateframe import extract
        if fix_tau_lo == 'first':
            icorr=0
        elif fix_tau_lo == 'last':
            icorr=3
        else:
            print 'Invalid value for fix_tau_lo.  Must be "first" or "last"'
            return
        xml, buf = ch.read_cal(4,t=Time(out[icorr]['time'][0],format='jd'))
        dlatbl = extract(buf,xml['Delaycen_ns'])
        dtau_x, dtau_y = dlatbl[14] - dlatbl[13]
        dp_x = out[icorr]['fghz']*2*np.pi*dtau_x
        dp_y = out[icorr]['fghz']*2*np.pi*dtau_y
        nt, = out[icorr]['time'].shape
        for i in range(nt):
            for iant in range(13):
                out[icorr]['x'][ri.bl2ord[iant,13],0,:,i] *= np.exp(1j*dp_x)
                out[icorr]['x'][ri.bl2ord[iant,13],1,:,i] *= np.exp(1j*dp_y)
                out[icorr]['x'][ri.bl2ord[iant,13],2,:,i] *= np.exp(1j*dp_y)
                out[icorr]['x'][ri.bl2ord[iant,13],3,:,i] *= np.exp(1j*dp_x)
    dph_lo = get_xy_corr(out[[3,0]], ant_str=ant_str, doplot=False)
    dph_hi = get_xy_corr(out[[2,1]], ant_str=ant_str)
    fghz = np.union1d(dph_lo['fghz'],dph_hi['fghz'])
    # Check for LO and HI being off by pi due to pi-ambiguity in xi_rot
    lo_com, hi_com = common_val_idx(dph_lo['fghz'],dph_hi['fghz'])
    # Average xi_rot angle difference over common frequencies
    a = lobe(dph_hi['xi_rot'][hi_com]-dph_lo['xi_rot'][lo_com]) # angle difference
    xi_rot_diff = np.angle(np.sum(np.exp(1j*a)))  # Average angle difference
    if np.abs(xi_rot_diff) > np.pi/2:
        # Looks like shifting by pi will get us closer, so shift both xyphase and xi_rot
        # This does not actually change any phase relationships, it only makes the plots
        # and HI/LO data comparison more consistent.
        dph_lo['xyphase'] += np.pi
        dph_lo['xi_rot'] += np.pi
        dph_lo['dphi'] += np.pi
        dph_lo['dph14'] += np.pi
    ax = plt.figure('XY_Phase').get_axes()
    for i in range(13):
        ax[i].plot(dph_lo['fghz'], lobe(dph_lo['dphi'][0,i]), '.',color='C0')    
        ax[i].plot(dph_lo['fghz'], lobe(dph_lo['dphi'][1,i]), '.',color='C1')    
        ax[i].plot(dph_lo['fghz'],lobe(dph_lo['xyphase'][i]),'r.')
        ax[i].set_xlim(0,20)
    for i in range(len(dph_lo['dph14'][:,0])):
        ax[13].plot(dph_lo['fghz'],lobe(dph_lo['dph14'][i]),'.')
    ax[13].plot(dph_lo['fghz'],lobe(dph_lo['xyphase'][13]),'r.')
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
    xi_rot[idx_lo_not_hi] = lobe(dph_lo['xi_rot'][idx_lo_not_hi])   # For unique low-receiver frequencies, insert LO xi_rot
    ax[14].plot(fghz,xi_rot)
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
            
def sat_xy_corr(out0, out1, band=0, ant_str='ant1-13', doplot=True):
    ''' Analyze a pair of parallel and cross polarization calibration packet captures
        on a Geosat in K-band (bands 33, 34, 35, 36, 37) and
        return the X vs. Y delay phase corrections on all antennas 1-14.
        
        Required keyword:
           prtlist   a list of 2 PRT filenames, the first being the 
                       parallel-feed scan, and the second being the
                       crossed-feed scan.
           band      the band index (0-4) corresponding to the above 5 bands
        Optional keyword:
           doplot    True => plot the final result, False => no plot
    '''
    import pcapture2 as p
    from util import bl2ord, ant_str2list
        
    if doplot: import matplotlib.pylab as plt

    antlist = ant_str2list(ant_str)
    nant = len(antlist)
#    out0 = p.rd_jspec(prtlist[0])  # Parallel scan
#    out1 = p.rd_jspec(prtlist[1])  # Perpendicular scan
    # Integrate over (10) repeated records for the desired band
    ph0 = np.angle(np.sum(out0['x'][:,:,:,10*band:10*(band+1)],3))
    ph1 = np.angle(np.sum(out1['x'][:,:,:,10*band:10*(band+1)],3))
    # Determine secular change in phase at the two times, relative to ant 1
    for i in antlist:
        if i == antlist[0]:
            dp = np.zeros_like(ph0[0,0])
        else:
            dp = lobe(ph1[bl2ord[antlist[0],i],0] - ph0[bl2ord[antlist[0],i],0])
        ph0[bl2ord[i,13],2:] = lobe(ph1[bl2ord[i,13],2:]-dp)     # Insert crossed-feed phases from ph1 into ph0, corrected for secular change

    ph0 = ph0[bl2ord[antlist,13]]          # Now restrict to only baselines with ant 14
    fstart = (band+32)*0.325 + 1.1 - 0.025
    fghz = np.linspace(fstart,fstart+0.400,4096)
    nf = len(fghz)
    dph = np.zeros((nant+1,nf),np.float)
    # Determine xi_rot
    xi2 = ph0[:,2] - ph0[:,0] + ph0[:,3] - ph0[:,1]  # This is 2 * xi, measured separately on each of 13 antennas
    xi_rot = lobe(np.unwrap(np.angle(np.sum(np.exp(1j*xi2),0)))/2.)   # Very clever average does not suffer from wrapping issues
    #xi_rot = np.zeros_like(xi_rot) + np.pi/2.    # *********** Zero out xi_rot for now ****************
    # Form differential delay phase from channels, and average them
    # dph14 = XY - XX - xi_rot and YY - YX + xi_rot
    dph14 = np.concatenate((lobe(ph0[:,2] - ph0[:,0] - xi_rot),lobe(ph0[:,1] - ph0[:,3] + xi_rot)))  # 26 values for Ant 14
    dph[nant] = np.angle(np.sum(np.exp(1j*dph14),0))  # Very clever average does not suffer from wrapping issues
    # dphi = XX - YX + xi_rot and XY - YY - xi_rot 
    dphi = np.array((lobe(ph0[:,0] - ph0[:,3] + xi_rot),lobe(ph0[:,2] - ph0[:,1] - xi_rot)))  # 2 values for Ant 14
    dph[:nant] = np.angle(np.sum(np.exp(1j*dphi),0))
    
    if doplot:
        figlabel = 'XY_Phase_'+str(band+33)
        if figlabel in plt.get_figlabels():
            f = plt.figure(figlabel)
            ax = f.get_axes()
        else:
            f, ax = plt.subplots(4, 4, num=figlabel)
            ax.shape = (16,)
        for i in range(nant): 
            ax[antlist[i]].plot(fghz,dphi[0,i],',')
            ax[antlist[i]].plot(fghz,dphi[1,i],',')
            ax[antlist[i]].plot(fghz,dph[i],'k,')
            ax[antlist[i]].set_title('Ant '+str(antlist[i]+1),fontsize=9)
        for i in range(2*nant):
            ax[13].plot(fghz,dph14[i],',')
            ax[13].set_title('Ant 14',fontsize=9)
        ax[13].plot(fghz,dph[nant],'k,')
        for i in range(14): ax[i].set_ylim(-4,4)
        f.suptitle('Multicolor: Measurements, Black: Final Results')
        ax[14].plot(fghz,xi_rot)
        for i in range(nant):
            ax[15].plot(fghz,xi2[i],',')
#    np.savez('/common/tmp/Feed_rotation/' + npzlist[0].split('/')[-1][:14] + '_delay_phase.npz', fghz=fghz, dph=dph, xi_rot=xi_rot)
    time = Time.now().lv
    xy_phase = {'antlist':antlist, 'timestamp':time, 'fghz':fghz, 'xyphase':dph, 'xi_rot':xi_rot, 'dphi':dphi, 'dph14':dph14}
    return xy_phase

def sat_unrot(data,xyphase,band=0):
    from util import bl2ord
    xi_rot = xyphase['xi_rot']
    dph = xyphase['xyphase']
    nf = 4096
    fidx2 = np.arange(nf)
    for i in range(13):
        for j in range(i + 1, 14):
            k = bl2ord[i, j]
            if j == 13:                  # xi_rot was applied for all antennas, but this
                xi = xi_rot[fidx2]       # is wrong.  Now it is only done for ant14.
            else:
                xi = 0.0                 # xi_rot for other antennas is just zero.
            a1 = lobe(dph[i, fidx2] - dph[j, fidx2])
            a2 = -dph[j, fidx2] - xi
            a3 = dph[i, fidx2] - xi + np.pi
            data['x'][k, 1, fidx2, 10*band:10*(band+1)] *= np.repeat(np.exp(1j * a1), 10).reshape(nf, 10)
            data['x'][k, 2, fidx2, 10*band:10*(band+1)] *= np.repeat(np.exp(1j * a2), 10).reshape(nf, 10)
            data['x'][k, 3, fidx2, 10*band:10*(band+1)] *= np.repeat(np.exp(1j * a3), 10).reshape(nf, 10)
    return data