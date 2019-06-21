#
# Pipeline Calibration Routines
#
# These routines are used to apply calibration to the output of the
# udb_util.py's readXdata() routine.
#
# History
#  2017-08-06  DG
#    Initial start of file
#  2017-08-07  DG
#    Further additions and debugging
#  2017-08-14  DG
#    Implemented flagging for non-tracking antennas, in unrot().
#  2017-09-02  DG
#    Change to allow concatenation of a list of files in udb_corr()
#  2017-09-09  DG
#    Added apply_fem_level() routine, similar to apply_gain_state()
#    except it takes account of non-uniform attenuation vs. frequency
#    as measured by GAINCALTEST.  Also added get_calfac() and apply_calfac()
#    routines, which can optionally be applied to the data, although
#    this is likely not going to be done to IDB/UDB data, but rather
#    will be applied as a CASA bandpass calibration table.
#  2017-09-11 DG
#    Updated udb_corr() to optionally apply the new scheme of calibration
#    using the GAINCALTEST results if "new" is True (default).
#  2017-10-13  DG
#    Changed apply_fem_level() to read attn values from SQL rather than
#    calculating it on the fly.
#  2017-12-03  DG
#    Updated the calibrate=True part of udb_corr() to give an error and
#    do nothing if no record for TP calibration has been entered in SQL.
#  2017-12-05  DG
#    The above did not work, because the SQL timestamp was not available
#    from cal_header's read_cal() routine.  I fixed that, and changed get_calfac()
#    to add the SQL_timestamp as the key sqltime.
#  2018-01-03  DG
#    Updated udb_corr() to work with new xi_rot.
#  2018-03-22  NK
#    Added an optional keyword attncal in udb_corr() that, when set to False, 
#    will skip the attenuation correction.
#  2018-06-13  DG
#    Added code to get_sql_info() and unrot() to determine which 27-m
#    receiver is in use, and correct baselines with Ant 14 accordingly.
#    The X vs. Y delay is in general very different for the two receivers.
#  2019-04-28  DG
#    I discovered a major screw up in applying xi_rot when unrotating
#    phases.  The xi_rot term should only apply to ant14, never to any
#    other antennas.  So now xi_rot is set to zero for all except ant14.
#  2019-05-22  DG
#    Updated apply_fem_level to subtract receiver noise before scaling,
#    if a skycal or gaincal is available.  The call to skycal_anal()
#    occurs in udb_corr()
#  2019-06-02  DG
#    Important change to apply_cal, which now does offsun subtraction and
#    zeroes-out the power-squared since it is no longer valid after calibration
#  2019-06-20  DG
#    Another change to apply_cal, which now does auto-correlation offsun
#    subtraction.  Note that this only works for XX and YY, and XY and YX
#    are left with no subtraction.  Also added subtraction of auto-correlation
#    receiver noise in apply_fem_level().
#

import dbutil as db
import numpy as np
from util import Time, nearest_val_idx, common_val_idx, lobe, bl2ord
import stateframe
import cal_header as ch



def get_sql_info(trange):
    ''' Get all antenna information from the SQL database for a given
        timerange, including TrackFlag and Parallactic Angle
        
        Also determines if the RFSwitch state (i.e. which 27-m receiver 
        is being used).
    '''
    cursor = db.get_cursor()
    sqldict = db.get_dbrecs(cursor, dimension=15, timestamp=trange)
    azeldict = stateframe.azel_from_sqldict(sqldict)
    time = Time(sqldict['Timestamp'][:, 0].astype(int), format='lv')
    azeldict.update({'Time': time})
    sqldict = db.get_dbrecs(cursor, dimension=1, timestamp=trange)
    azeldict.update({'RFSwitch':sqldict['FEMA_Powe_RFSwitchStatus']})
    azeldict.update({'LF_Rcvr':sqldict['FEMA_Rece_LoFreqEnabled']})
    if np.median(azeldict['RFSwitch']) == 0.0 and np.median(azeldict['LF_Rcvr']) == 1.0:
        azeldict.update({'Receiver':'Low'})
    elif np.median(azeldict['RFSwitch']) == 1.0 and np.median(azeldict['LF_Rcvr']) == 0.0:
        azeldict.update({'Receiver':'High'})
    else:
        azeldict.update({'Receiver':'Unknown'})
    cursor.close()
    return azeldict


def apply_fem_level(data, gctime=None, skycal=None):
    ''' Applys the FEM level corrections to the given data dictionary.
        
        Inputs:
          data     A dictionary such as that returned by readXdata().
          gctime   A Time() object whose date specifies which GAINCALTEST
                     measurements to use.  If omitted, the date of the data
                     is used.
          skycal   Optional array of receiver noise from SKYCAL or GAINCAL
                     calibration.  Only the receiver noise is applied (subtracted)

        Output:
          cdata    A dictionary with the level-corrected data.  The keys
                     p, x, p2, and a are all updated.
    '''
    import attncal as ac
    from gaincal2 import get_fem_level
    import copy

    # Get timerange from data
    trange = Time([data['time'][0], data['time'][-1]], format='jd')
    if gctime is None:
        gctime = trange[0]
    # Get time cadence
    dt = np.int(np.round(np.median(data['time'][1:] - data['time'][:-1]) * 86400))
    if dt == 1: dt = None
    # Get the FEM levels of the requested timerange
    src_lev = get_fem_level(trange, dt)  # solar gain state for timerange of file
    nf = len(data['fghz'])
    nt = len(src_lev['times'])
    attn = ac.read_attncal(gctime)[0]  # Reads attn from SQL database (returns a list, but use first, generally only, one)
    # attn = ac.get_attncal(gctime)[0]   # Analyzes GAINCALTEST (returns a list, but use first, generally only, one)
    antgain = np.zeros((15, 2, nf, nt), np.float32)  # Antenna-based gains [dB] vs. frequency
    # Find common frequencies of attn with data
    idx1, idx2 = common_val_idx(data['fghz'], attn['fghz'], precision=4)
    # Currently, GAINCALTEST measures 8 levels of attenuation (16 dB).  I assumed this would be enough,
    # but the flare of 2017-09-10 actually went to 10 levels (20 dB), so we have no choice but to extend
    # to higher levels using only the nominal, 2 dB steps above the 8th level.  This part of the code
    # extends to the maximum 16 levels.
    a = np.zeros((16, 13, 2, nf), float)  # Extend attenuation to 14 levels
    a[1:9, :, :, idx1] = attn['attn'][:, :13, :, idx2]  # Use GAINCALTEST results in levels 1-9 (bottom level is 0dB)
    for i in range(8, 15):
        # Extend to levels 9-15 by adding 2 dB to each previous level
        a[i + 1] = a[i] + 2.
    a[15] = 62.  # Level 15 means 62 dB have been inserted.
    #print 'Attn list (dB) for ant 1, pol xx, lowest frequency:',a[:,0,0,0]
    for i in range(13):
        for k, j in enumerate(idx1):
            antgain[i, 0, j] = a[src_lev['hlev'][i], i, 0, j]
            antgain[i, 1, j] = a[src_lev['vlev'][i], i, 0, j]
    cdata = copy.deepcopy(data)
    nblant = 136
    blgain = np.zeros((nf, nblant, 4, nt), float)  # Baseline-based gains vs. frequency

    for i in range(14):
        for j in range(i, 14):
            k = bl2ord[i, j]
            blgain[:, k, 0] = 10 ** ((antgain[i, 0] + antgain[j, 0]) / 20.)
            blgain[:, k, 1] = 10 ** ((antgain[i, 1] + antgain[j, 1]) / 20.)
            blgain[:, k, 2] = 10 ** ((antgain[i, 0] + antgain[j, 1]) / 20.)
            blgain[:, k, 3] = 10 ** ((antgain[i, 1] + antgain[j, 0]) / 20.)
    # Reorder antgain axes to put frequencies in first slot, to match data
    antgain = np.swapaxes(np.swapaxes(antgain, 1, 2), 0, 1)
    antgainf = 10 ** (antgain / 10.)

    idx = nearest_val_idx(data['time'], src_lev['times'].jd)
    nt = len(idx)  # New number of times
    # If a skycal dictionary exists, subtract auto-correlation receiver noise before scaling (clip to 0)
    if skycal:
        sna, snp, snf = skycal['rcvr_bgd_auto'].shape
        bgd = skycal['rcvr_bgd_auto'].repeat(nt).reshape((sna,snp,snf,nt))
        # Reorder axes
        bgd = np.swapaxes(bgd,0,2)
        bslice = bgd[:,:,:,idx]
        for i in range(13):
            cdata['x'][:, bl2ord[i,i], 0] = np.clip(cdata['x'][:, bl2ord[i,i], 0] - bslice[:,0,i],0,None)
            cdata['x'][:, bl2ord[i,i], 1] = np.clip(cdata['x'][:, bl2ord[i,i], 1] - bslice[:,1,i],0,None)
    # Correct the auto- and cross-correlation data
    cdata['x'] *= blgain[:, :, :, idx]
    # Reshape px and py arrays
    cdata['px'].shape = (nf, 16, 3, nt)
    cdata['py'].shape = (nf, 16, 3, nt)
    # If a skycal dictionary exists, subtract total power receiver noise before scaling (clip to 0)
    # NB: This will break SK!
    if skycal:
        sna, snp, snf = skycal['rcvr_bgd'].shape
        bgd = skycal['rcvr_bgd'].repeat(nt).reshape((sna,snp,snf,nt))
        # Reorder axes
        bgd = np.swapaxes(bgd,0,2)
        bslice = bgd[:,:,:,idx]
        bgnd = np.rollaxis(bslice,3)
        cdata['px'][:, :13, 0] = np.clip(cdata['px'][:, :13, 0] - bslice[:,0],0,None)
        cdata['py'][:, :13, 0] = np.clip(cdata['py'][:, :13, 0] - bslice[:,1],0,None)
    # Correct the power
    cdata['px'][:, :15, 0] *= antgainf[:, :, 0, idx]
    cdata['py'][:, :15, 0] *= antgainf[:, :, 1, idx]
    # Correct the power-squared
    cdata['px'][:, :15, 1] *= antgainf[:, :, 0, idx] ** 2
    cdata['py'][:, :15, 1] *= antgainf[:, :, 1, idx] ** 2
    # Reshape px and py arrays back to original
    cdata['px'].shape = (nf * 16 * 3, nt)
    cdata['py'].shape = (nf * 16 * 3, nt)
    return cdata


def apply_attn_corr(data, tref=None):
    ''' NB: This routine has been superseded by apply_fem_level()
    
        Applies nominal attenuator state corrections to the given data dictionary,
        corrected to the gain-state at time given by Time() object tref.
        
        Inputs:
          data     A dictionary returned by udb_util.py's readXdata().
          tref     A Time() object with the reference time, or if None,
                     the gain state of the nearest earlier REFCAL is 
                     used.
          skycal   Optional array of receiver noise from SKYCAL or GAINCAL
                     calibration.  Only the receiver noise is applied (subtracted)

        Output:
          cdata    A dictionary with the gain-corrected data.  The keys
                     px, py, and x, are updated.
                     
        NB: This is the same routine as in gaincal2.py, but modified
        to handle the different ordering/format of data from udb_util.py's
        readXdata() routine.
    '''
    from gaincal2 import get_gain_state
    import copy
    if tref is None:
        # No reference time specified, so get nearest earlier REFCAL
        trange = Time(data['time'][[0, -1]], format='jd')
        xml, buf = ch.read_cal(8, t=trange[0])
        tref = Time(stateframe.extract(buf, xml['Timestamp']), format='lv')
    # Get the gain state at the reference time (actually median over 1 minute)
    trefrange = Time([tref.iso, Time(tref.lv + 60, format='lv').iso])
    ref_gs = get_gain_state(trefrange)  # refcal gain state for 60 s
    # Get median of refcal gain state (which should be constant anyway)
    ref_gs['h1'] = np.median(ref_gs['h1'], 1)
    ref_gs['h2'] = np.median(ref_gs['h2'], 1)
    ref_gs['v1'] = np.median(ref_gs['v1'], 1)
    ref_gs['v2'] = np.median(ref_gs['v2'], 1)

    # Get timerange from data
    trange = Time([data['time'][0], data['time'][-1]], format='jd')
    # Get time cadence
    dt = np.int(np.round(np.median(data['time'][1:] - data['time'][:-1]) * 86400))
    if dt == 1: dt = None
    # Get the gain state of the requested timerange
    src_gs = get_gain_state(trange, dt)  # solar gain state for timerange of file
    nt = len(src_gs['times'])
    antgain = np.zeros((15, 2, 34, nt), np.float32)  # Antenna-based gains vs. band
    for i in range(15):
        for j in range(34):
            antgain[i, 0, j] = src_gs['h1'][i] + src_gs['h2'][i] - ref_gs['h1'][i] - ref_gs['h2'][i] + \
                               src_gs['dcmattn'][i, 0, j] - ref_gs['dcmattn'][i, 0, j]
            antgain[i, 1, j] = src_gs['v1'][i] + src_gs['v2'][i] - ref_gs['v1'][i] - ref_gs['v2'][i] + \
                               src_gs['dcmattn'][i, 1, j] - ref_gs['dcmattn'][i, 1, j]

    cdata = copy.deepcopy(data)
    # Create giant array of baseline-based gains, translated to baselines and frequencies
    fghz = data['fghz']
    nf = len(fghz)
    blist = (fghz * 2 - 1).astype(int) - 1  # Band list corresponding to frequencies in data
    nblant = 136
    blgain = np.zeros((nf, nblant, 4, nt), float)  # Baseline-based gains vs. frequency
    for i in range(14):
        for j in range(i, 14):
            k = bl2ord[i, j]
            blgain[:, k, 0] = 10 ** ((antgain[i, 0, blist] + antgain[j, 0, blist]) / 20.)
            blgain[:, k, 1] = 10 ** ((antgain[i, 1, blist] + antgain[j, 1, blist]) / 20.)
            blgain[:, k, 2] = 10 ** ((antgain[i, 0, blist] + antgain[j, 1, blist]) / 20.)
            blgain[:, k, 3] = 10 ** ((antgain[i, 1, blist] + antgain[j, 0, blist]) / 20.)
    # Reorder antgain axes to put frequencies in first slot, to match data
    antgain = np.swapaxes(np.swapaxes(antgain, 1, 2), 0, 1)
    antgainf = 10 ** (antgain[blist] / 10.)

    idx = nearest_val_idx(data['time'], src_gs['times'].jd)
    nt = len(idx)  # New number of times
    # Correct the auto- and cross-correlation data
    cdata['x'] *= blgain[:, :, :, idx]
    # Reshape px and py arrays
    cdata['px'].shape = (nf, 16, 3, nt)
    cdata['py'].shape = (nf, 16, 3, nt)
    # Correct the power
    cdata['px'][:, :15, 0] *= antgainf[:, :, 0, idx]
    cdata['py'][:, :15, 0] *= antgainf[:, :, 1, idx]
    # Correct the power-squared
    cdata['px'][:, :15, 1] *= antgainf[:, :, 0, idx] ** 2
    cdata['py'][:, :15, 1] *= antgainf[:, :, 1, idx] ** 2
    # Reshape px and py arrays back to original
    cdata['px'].shape = (nf * 16 * 3, nt)
    cdata['py'].shape = (nf * 16 * 3, nt)
    return cdata


def get_calfac(t=None):
    ''' Read total power and auto-correlation calibration factors from the SQL
        database, for the time specified by Time() object t, or if None, at the
        next earlier calibration time to the current time.
    '''
    tpcal_type = 10  # Calibration type specified in cal_header.py
    if t is None:
        t = Time.now()
    xml, buf = ch.read_cal(tpcal_type, t=t)
    fghz = stateframe.extract(buf, xml['FGHz'])
    nf = len(fghz)
    tpcalfac = np.zeros((13, 2, nf), np.float)
    tpoffsun = np.zeros((13, 2, nf), np.float)
    accalfac = np.zeros((13, 2, nf), np.float)
    acoffsun = np.zeros((13, 2, nf), np.float)
    nant = len(xml['Antenna'])
    for i in range(nant):
        iant = stateframe.extract(buf, xml['Antenna'][i]['Antnum']) - 1
        tpcalfac[iant] = stateframe.extract(buf, xml['Antenna'][i]['TPCalfac'])
        accalfac[iant] = stateframe.extract(buf, xml['Antenna'][i]['ACCalfac'])
        tpoffsun[iant] = stateframe.extract(buf, xml['Antenna'][i]['TPOffsun'])
        acoffsun[iant] = stateframe.extract(buf, xml['Antenna'][i]['ACOffsun'])
    try:
        sqltime = stateframe.extract(buf, xml['SQL_timestamp'])
    except:
        sqltime = None
    return {'fghz': fghz, 'timestamp': stateframe.extract(buf, xml['Timestamp']), 'sqltime': sqltime,
            'tpcalfac': tpcalfac, 'accalfac': accalfac, 'tpoffsun': tpoffsun, 'acoffsun': acoffsun}


def apply_calfac(data, calfac):
    ''' Applies calibration factors in calfac dictionary returned by get_calfac(),
        to the data and returns the calibrated data in the same form.  No calibration
        can be applied for antennas/frequencies in data that are not included in
        calfac, so for those the data are returned unchanged, except the correlated
        data for missing frequencies are flagged.
        
        Inputs:
            data    The data to be calibrated
            calfac  The calfac dictionary returned by a call to get_calfac().
            
        Output:
            cdata   The calibrated data
            
        N.B.:  Spectral Kurtosis is no longer valid after applying calibration
    '''
    import copy
    fghz = data['fghz']
    nfin = len(fghz)
    # Find common frequencies
    idx1, idx2 = common_val_idx(fghz, calfac['fghz'], precision=4)
    nf = len(idx1)
    missing = np.setdiff1d(np.arange(len(fghz)), idx1)
    nant = 16
    nblant = nant * (nant - 1) / 2 + nant
    blfac =  np.ones((nf, nblant, 4), float)  # Factors for missing antennas/frequencies are set to unity
    bloff = np.zeros((nf, nblant, 4), float)  # Offsun values (will be zero except for auto-correlations)
    fac = calfac['accalfac'][:, :, idx2]  # Extract only common frequencies
    acoffsun = calfac['acoffsun'][:, :, idx2]  # Extract only common frequencies
    for i in range(13):
        for j in range(i, 13):
            blfac[idx1, bl2ord[i, j], 0] = np.sqrt(fac[i, 0] * fac[j, 0])
            blfac[idx1, bl2ord[i, j], 1] = np.sqrt(fac[i, 1] * fac[j, 1])
            blfac[idx1, bl2ord[i, j], 2] = np.sqrt(fac[i, 0] * fac[j, 1])
            blfac[idx1, bl2ord[i, j], 3] = np.sqrt(fac[i, 1] * fac[j, 0])
            if i == j:
                bloff[idx1, bl2ord[i, j], 0] = acoffsun[j, 0]
                bloff[idx1, bl2ord[i, j], 1] = acoffsun[j, 1]
    antfac = np.ones((nf, nant, 2), float)  # Factors for missing antennas/frequencies are set to unity
    offsun = np.zeros((nf, nant, 2), float) # Offsun for missing antennas/frequencies are set to zero
    for i in range(13):
        for j in range(2):
            antfac[idx1, i, j] = calfac['tpcalfac'][i, j, idx2]
            offsun[idx1, i, j] = calfac['tpoffsun'][i, j, idx2]
    cdata = copy.deepcopy(data)
    nt = len(data['time'])
    # Calibrate the auto- and cross-correlation data
    # Note that bloff is non-zero only for the real part of auto-correlations
    for i in range(nt):
        cdata['x'][idx1, :, :, i] = (cdata['x'][idx1, :, :, i]-bloff)*blfac
    # Reshape px and py arrays
    cdata['px'].shape = (nfin, 16, 3, nt)
    cdata['py'].shape = (nfin, 16, 3, nt)
    for i in range(nt):
        # Correct the power
        cdata['px'][idx1, :, 0, i] = (cdata['px'][idx1, :, 0, i] - offsun[:, :, 0])*antfac[:, :, 0]
        cdata['py'][idx1, :, 0, i] = (cdata['py'][idx1, :, 0, i] - offsun[:, :, 1])*antfac[:, :, 1]
        # Zero-out the power-squared, since SK is no longer valid after applying calibration
        cdata['px'][idx1, :, 1, i] = 0.0
        cdata['py'][idx1, :, 1, i] = 0.0
    # Reshape px and py arrays back to original
    cdata['px'].shape = (nfin * 16 * 3, nt)
    cdata['py'].shape = (nfin * 16 * 3, nt)
    # Set flags for any missing frequencies (hopefully this also works when "missing" is np.array([]))
    cdata['x'][missing] = np.ma.masked
    return cdata


def unrot(data, azeldict=None):
    ''' Apply the correction to differential feed rotation to data, and return
        the corrected data.  This also applies flags to data whose antennas are
        not tracking.

        Inputs:
          data     A dictionary returned by udb_util.py's readXdata().
          azeldict The dictionary returned from get_sql_info(), or if None, the appropriate
                     get_sql_info() call is done internally.

        Output:
          cdata    A dictionary with the phase-corrected data.  Only the key
                     x is updated.
    '''
    import copy
    trange = Time(data['time'][[0, -1]], format='jd')

    if azeldict is None:
        azeldict = get_sql_info(trange)
    chi = azeldict['ParallacticAngle'] * np.pi / 180.  # (nt, nant)
    # Correct parallactic angle for equatorial mounts, relative to Ant14
    chi[:, [8, 9, 10, 12, 13]] = 0  # Currently 0, but can be measured and updated

    # Which antennas are tracking
    track = azeldict['TrackFlag']  # True if tracking

    # Ensure that nearest valid parallactic angle is used for times in the data
    good = np.where(azeldict['ActualAzimuth'] != 0)
    tidx = []  # List of arrays of indexes for each antenna
    for i in range(14):
        gd = good[0][np.where(good[1] == i)]
        tidx.append(nearest_val_idx(data['time'], azeldict['Time'][gd].jd))

    # Read X-Y Delay phase from SQL database and get common frequencies
    xml, buf = ch.read_cal(11, t=trange[0])
    fghz = stateframe.extract(buf, xml['FGHz'])
    good, = np.where(fghz != 0.)
    fghz = fghz[good]
    dph = stateframe.extract(buf, xml['XYphase'])
    dph = dph[:, good]
    xi_rot = stateframe.extract(buf, xml['Xi_Rot'])
    xi_rot = xi_rot[good]
    fidx1, fidx2 = common_val_idx(data['fghz'], fghz, precision=4)
    missing = np.setdiff1d(np.arange(len(data['fghz'])), fidx1)

    nf, nbl, npol, nt = data['x'].shape
    nf = len(fidx1)
    # Correct data for X-Y delay phase
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
            data['x'][fidx1, k, 1] *= np.repeat(np.exp(1j * a1), nt).reshape(nf, nt)
            data['x'][fidx1, k, 2] *= np.repeat(np.exp(1j * a2), nt).reshape(nf, nt)
            data['x'][fidx1, k, 3] *= np.repeat(np.exp(1j * a3), nt).reshape(nf, nt)

    # Correct data for differential feed rotation
    cdata = copy.deepcopy(data)
    for n in range(nt):
        for i in range(13):
            for j in range(i + 1, 14):
                k = bl2ord[i, j]
                ti = tidx[i][n]
                tj = tidx[j][n]
                if track[ti, i] and track[tj, j]:
                    dchi = chi[ti, i] - chi[tj, j]
                    cchi = np.cos(dchi)
                    schi = np.sin(dchi)
                    cdata['x'][:, k, 0, n] = data['x'][:, k, 0, n] * cchi + data['x'][:, k, 3, n] * schi
                    cdata['x'][:, k, 2, n] = data['x'][:, k, 2, n] * cchi + data['x'][:, k, 1, n] * schi
                    cdata['x'][:, k, 3, n] = data['x'][:, k, 3, n] * cchi - data['x'][:, k, 0, n] * schi
                    cdata['x'][:, k, 1, n] = data['x'][:, k, 1, n] * cchi - data['x'][:, k, 2, n] * schi
                else:
                    cdata['x'][:, k, :, n] = np.ma.masked

    # Set flags for any missing frequencies (hopefully this also works when "missing" is np.array([]))
    cdata['x'][missing] = np.ma.masked
    return cdata


def udb_corr(filelist, outpath='./', calibrate=False, new=True, gctime=None, attncal=True):
    ''' Complete routine to read in an existing idb or udb file and output
        a new file of the same name in the local directory, with all corrections
        applied.
        
        Inputs:
          filelist  List of files to read.  The output of all files in the list
                        are concatenated into a single output file based on the first
                        file in the list.  Use an external loop over files if this is
                        not wanted.
          calibrate A boolean flag to indicate whether calibration factors should be
                        applied to the data.  The default is False, in anticipation that
                        such calibration will be applied as a CASA bandpass table, but
                        if set to True, the SOLPNTCAL analysis is used.
          new       If True (default), the "new" scheme of attenuation corrections
                        based on GAINCALTEST results is applied.  Otherwise, only the
                        nominal attenuation corrections are applied.
          gctime    A Time() object whose date is used to find the GAINCALTEST data.
                        If None (default), then the date of the data is used.  Note that
                        gctime is only used if parameter new is True.
          attncal   If False, the attenuation correction is skipped - expected to be
                        applied manually in post-processing (e.g. 2017-09-10 X8 flare)          
    '''
    import sys
    import os
    import udb_util as uu
    import time
    from pathlib2 import Path
    from copy import deepcopy
    if type(filelist) is str or type(filelist) is np.string_:
        # Convert input filename to list if not already a list
        filelist = [filelist]

    for idx, file in enumerate(filelist):
        if file[-1] == '/':
            filelist[idx] = file[:-1]

    filecount = 0
    for filename in filelist:
        t1 = time.time()
        out = uu.readXdata(filename)
        # temp
        nf, nt = np.array(out['px'].shape)/np.array([3*16,1])
        out['px'].shape = (nf, 16, 3, nt)
        f = open('pctemp.dat','ab')
        f.write(deepcopy(out['px'][:,0,0,25]))
        out['px'].shape = (nf*16*3,nt)
        f.close()
        ## end temp
        print 'Reading file took', time.time() - t1, 's'
        sys.stdout.flush()
        trange = Time(out['time'][[0, -1]], format='jd')
        t1 = time.time()
        azeldict = get_sql_info(trange)
        print 'Reading SQL info took', time.time() - t1, 's'
        sys.stdout.flush()
        ## Correct data for attenuation changes
        if attncal:
            from calibration import skycal_anal
            t1 = time.time()
            skycal = skycal_anal(t=trange[0],do_plot=False)
            if new:
                # Subtract receiver noise, then correct for front end attenuation
                cout = apply_fem_level(out, gctime, skycal=skycal)
            else:
                cout = apply_attn_corr(out)
            print 'Applying attn correction took', time.time() - t1, 's'
            sys.stdout.flush()
            ## temp
            nf, nt = np.array(out['px'].shape)/np.array([3*16,1])
            cout['px'].shape = (nf, 16, 3, nt)
            f = open('pctemp.dat','ab')
            f.write(deepcopy(cout['px'][:,0,0,25]))
            cout['px'].shape = (nf*16*3,nt)
            f.close()
            ## end temp
            t1 = time.time()
        else:
            cout = out
        # Correct data for differential feed rotation
        coutu = unrot(cout, azeldict)
        print 'Applying feed rotation correction took', time.time() - t1, 's'
        sys.stdout.flush()
        # Optionally apply calibration to convert to solar flux units
        if calibrate:
            t1 = time.time()
            if trange[0].datetime.hour < 7:
                # Data time is earlier than 7 UT (i.e. on previous local day) so
                # use previous date at 20 UT.
                mjd = int(trange[0].mjd) - 1 + 20. / 24
            else:
                # Use current date at 20 UT
                mjd = int(trange[0].mjd) + 20. / 24
            calfac = get_calfac(Time(mjd, format='mjd'))
            if Time(calfac['sqltime'], format='lv').mjd == mjd:
                coutu = apply_calfac(coutu, calfac)
                ## temp
                nf, nt = np.array(out['px'].shape)/np.array([3*16,1])
                coutu['px'].shape = (nf, 16, 3, nt)
                f = open('pctemp.dat','ab')
                f.write(deepcopy(coutu['px'][:,0,0,25]))
                coutu['px'].shape = (nf*16*3,nt)
                f.close()
                ## end temp
                print 'Applying calibration took', time.time() - t1, 's'
            else:
                print 'Error: no TP calibration for this date.  Skipping calibration.'
        sys.stdout.flush()
        filecount += 1
        if filecount == 1:
            x = coutu
        else:
            x = uu.concatXdata(x, coutu)
    ufilename = outpath + filelist[0].split('/')[-1]
    while Path(ufilename).exists():
        # Handle case of existing file, by appending _n, where n increments (up to 9)
        if ufilename[-2] == '_':
            ufilename = ufilename[:-1] + str(int(ufilename[-1]) + 1)
        else:
            ufilename += '_1'
    ufile_out = uu.udbfile_write(x, filelist[0], ufilename)
    return ufilename

