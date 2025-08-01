# gaincal2
# Modification History
#  2017-05-13  DG
#    First wrote routine get_gain_state()
#  2017-05-15  DG
#    Added get_fseqfile() and fseqfile2bandlist() to help with conversion
#    of slots to bands.  Also ensure that all arrays returned by get_gain_state() 
#    are in canonical order [nant, npol, nf/nband, nt].
#  2017-05-22  DG
#    Some changes to apply_gain_corr() to make it more general.
#  2017-06-12  DG
#    Important change to apply_gain_corr() so that it does not drop unmeasured
#    points.  Instead, it uses the gain state of the nearest time in case of a missing one
#  2017-07-13  DG
#    Fixed a bug in apply_gain_corr() -- defined trange, when tref not given.
#  2017-07-15  DG
#    Rather extensive changes to allow for correct gain correction for averaged
#    data.  When apply_gain_corr() is called with data with a cadence dt longer than
#    1 s (e.g. UDB data), it will detect it and apply the approporiate average
#    correction.
#  2017-08-06  DG
#    Changed get_gain_state() to return data for inclusive timerange, i.e. it
#    includes the begin and end times.
#  2017-09-09  DG
#    Added get_fem_level() and apply_fem_level() routines, similar to get_gain_state()
#    and apply_gain_state(), except these take account of non-uniform attenuation
#    vs. frequency, based on GAINCALTEST measurements. 
#  2018-01-26  DG
#    First try to read from SQL in apply_fem_level(), and then go to data only if
#    that fails.
#  2019-02-25  DG
#    Fixed the bugs where the number of bands was hardcoded as 34. Now the number of
#    bands is determined from the arrays returned from the SQL database.
#  2019-07-22  DG
#    Initial writing of eq_anal(), which will ultimately be used to set the EQ
#    coefficients of the correlator.
#  2019-10-24  DG
#    Strangely, when a SQL attncal was not found, and it was then calculated from the
#    data, the code did not then write it to SQL.  It seems obvious that once it is
#    calculated a new SQL record should be written, so that is now added.
#  2020-08-24  DG
#    Change apply_gain_corr() to use the first PHASECAL of the requested date to
#    determine reference gain state.  This will be wrong if $FEM_INIT is issued
#    subsequent to that PHASECAL...
#  2021-01-21  DG
#    Found a bug in apply_fem_level(), where polarization index for Vpol was set to 0.
#    Also changed get_fem_level() to return times from the middle of dt rather than
#    the start of dt.
#  2021-01-23  DG
#    Important change to get_fem_level() for non-None dt case (mainly UDB data that
#    have a 60-s integration).  Now determines the proportion of time for each fem
#    attenuation state within integration time dt and returns a dictionary for each 
#    antenna and integrated time sample.  This required corresponding changes to
#    apply_fem_level().
#  2021-10-30  DG
#    A bug occurs in get_gain_state() during times when the 27-m is not working, so 
#    that we have no reference calibrations.  The apply_gain_corr() routine uses
#    an early time (13:30 UT) as reference time, but the SQL query fails if there
#    are no records at that time.  Now I have introduced a fall-back option that
#    uses the keyword relax=True in the call to get_gain_state(), which will just
#    take the first nt records after 13:30 UT.  This is not guaranteed to be a
#    reference gain state, but it has some chance to be so.  This is a very rare
#    occurrence.
#  2025-01-11  DG
#    Changes needed because of changes to SQL tables (16 ants and DCMs).  This
#    version simply truncates the output to the usual 15 channels, but when the
#    new dishes are online this will have to a changed again.
#  2025-05-22 DG
#    Further changes to handle the real 16 antennas.  Note that there is
#    parallel code in gaincal2.py and pipeline_cal.py, but the arrays they work
#    with are entirely different.  Here the arrays are, e.g. out['x'] which are
#    only cross-correlations (120 baselines), while in pipeline_cal.py the arrays
#    are all correlations (136 channels including auto-correlations).  This caused
#    some confusion initially, but hopefully I have it straight now.
#
import dbutil as db
import read_idb as ri
import cal_header as ch
from util import Time, nearest_val_idx
import stateframe as stf
import numpy as np

def get_fseqfile(t=None):
    if t is None:
        # 10 s ago...
        tlv = Time.now().lv - 10
    else: 
        tlv = int(t.lv)
    cursor = db.get_cursor()
    ver = db.find_table_version(cursor,tlv)
    # Get front end attenuator states
    query = 'select Timestamp,LODM_LO1A_FSeqFile from fV'+ver+'_vD1 where Timestamp between '+str(tlv)+' and '+str(tlv+1)+' order by Timestamp'
    data, msg = db.do_query(cursor, query)
    if msg == 'Success':
        fseqfile = data['LODM_LO1A_FSeqFile'][0].replace('\x00','')
        if fseqfile == 'none':
            fseqfile = None
    else:
        print 'Error: ',msg
        fseqfile = None
    cursor.close()
    return fseqfile
    
def fseqfile2bandlist(fseqfile=None):
    ''' Reads named fseqfile from ACC and returns the list of bands for
        the 50 slots of each 1-s sequence.
        
        Input:
           fseqfile    string filename (must exist in ACC:/parm folder.
        
        Returns:
           bandlist    numpy 50-element integer array of band numbers, 1-nbands
    '''
    import urllib2
    if fseqfile is None:
        print 'Must specify a frequency sequence.'
        return None
    userpass = 'admin:observer@'
    fseq_handle = urllib2.urlopen('ftp://'+userpass+'acc.solar.pvt/parm/'+fseqfile,timeout=0.5)
    lines = fseq_handle.readlines()
    fseq_handle.close()
    for line in lines:
        if line.find('LIST:SEQUENCE') != -1:
            line = line[14:]
            bands = np.array(map(int,line.split(',')))
        elif line.find('LIST:DWELL') != -1:
            line = line[11:]
            dwellist = np.array(map(float,line.split(',')))
    bandlist = []
    for band in bands:
        bandlist += [band]*int(np.round(dwellist[band-1]/0.02))
    return np.array(bandlist)
    
def get_fem_level(trange, dt=None):
    ''' Get FEM attenuation levels for a given timerange.  Returns a dictionary
        with keys as follows:

        times:     A Time object containing the array of times, size (nt)
        hlev:      The FEM attenuation level for HPol, size (nt, 15) 
        vlev:      The FEM attenuation level for VPol, size (nt, 15)
        dcmattn:   The base DCM attenuations for nbands (34 or 52) bands x 15 antennas x 2 Poln, size (nbands,30)
                      The order is Ant1 H, Ant1 V, Ant2 H, Ant2 V, etc.
        dcmoff:    If DPPoffset-on is 0, this is None (meaning there are no changes to the
                      above base attenuations).  
                   If DPPoffset-on is 1, then dcmoff is a table of offsets to the 
                      base attenuation, size (nt, 50).  The offset applies to all 
                      antennas/polarizations.
                      
        Optional keywords:
           dt      Seconds between entries to read from SQL stateframe database. 
                     If omitted, 1 s is assumed.
        
    '''
    def proportion(seq):
        '''Return dict of proportion of each value in "seq".'''
        hist = {}
        n = len(seq)
        for i in seq:
            hist[i] = hist.get(i, 0) + 1./n
        return hist

    if dt is None:
        tstart,tend = [str(i) for i in trange.lv]
    else:
        # Expand time by 1/2 of dt before and after
        tstart = str(np.round(trange[0].lv - dt/2))
        tend = str(np.round(trange[1].lv + dt/2))
    cursor = db.get_cursor()
    ver = db.find_table_version(cursor,trange[0].lv)
    if int(ver) > 66:
        tdim = 16   # Antenna info is in dimension 16 table starting with version 67
    else:
        tdim = 15
    if trange[0] < Time('2025-05-22'):
        nant = 15
    else:
        nant = 16
    # Get front end attenuator states
    query = 'select Timestamp,Ante_Fron_FEM_Clockms,' \
            +'Ante_Fron_FEM_HPol_Regi_Level,Ante_Fron_FEM_VPol_Regi_Level from fV' \
            +ver+'_vD'+str(tdim)+' where Timestamp >= '+tstart+' and Timestamp <= '+tend+' order by Timestamp'
    data, msg = db.do_query(cursor, query)
    if msg == 'Success':
        if dt:
            # If we want other than full cadence, get new array shapes and times
            n = len(data['Timestamp'])  # Original number of times
            new_n = (n/tdim/dt)*tdim*dt     # Truncated number of times equally divisible by dt
            new_shape = (n/tdim/dt,tdim) # New shape of truncated arrays
            times = Time(data['Timestamp'][:new_n].astype('int')[int(tdim*dt/2)::tdim*dt],format='lv')
        else:
            times = Time(data['Timestamp'].astype('int')[::tdim],format='lv')
        hlev = data['Ante_Fron_FEM_HPol_Regi_Level'].astype(int)
        vlev = data['Ante_Fron_FEM_VPol_Regi_Level'].astype(int)
        ms = data['Ante_Fron_FEM_Clockms']
        nt = len(hlev)/tdim
        hlev.shape = (nt,tdim)
        vlev.shape = (nt,tdim)
        ms.shape = (nt,tdim)
        # Find any entries for which Clockms is zero, which indicates where no
        # gain-state measurement is available.
        for i in range(nant):
            bad, = np.where(ms[:,i] == 0)
            if bad.size != 0 and bad.size != nt:
                # Find nearest adjacent good value
                good, = np.where(ms[:,i] != 0)
                idx = nearest_val_idx(bad,good)
                hlev[bad,i] = hlev[good[idx],i]
                vlev[bad,i] = vlev[good[idx],i]
        if dt:
            # If we want other than full cadence, find proportion of each level
            # during the dt time interval
            hlevout = []
            vlevout = []
            # Determine proportion of each state within each dt time integration
            for i in range(new_n/dt/tdim):
                for j in range(tdim):
                    hlevout.append(proportion(hlev[i*dt:(i+1)*dt,j]))
                    vlevout.append(proportion(vlev[i*dt:(i+1)*dt,j]))
            # Note, for the dt case hlev and vlev are an array of dicts with keys being
            # the level and values being the proportion of that level within the integration
            hlev = np.array(hlevout).reshape(new_shape)
            vlev = np.array(vlevout).reshape(new_shape)
        # Put results in canonical order [nant, nt]
        hlev = hlev.T
        vlev = vlev.T
    else:
        print 'Error reading FEM levels:',msg
        return {}
    # Get back end attenuator states
    xml, buf = ch.read_cal(2, t=trange[0])
    dcmattn = stf.extract(buf,xml['Attenuation'])
    nbands = dcmattn.shape[0]
    nants = dcmattn.shape[1]//2
    dcmattn.shape = (nbands, nants, 2)
    # Put into canonical order [nant, npol, nband]
    dcmattn = np.moveaxis(dcmattn,0,2)
    # See if DPP offset is enabled
    query = 'select Timestamp,DPPoffsetattn_on from fV' \
            +ver+'_vD1 where Timestamp >= '+tstart+' and Timestamp <= '+tend+'order by Timestamp'
    data, msg = db.do_query(cursor, query)
    if msg == 'Success':
        dppon = data['DPPoffsetattn_on']
        if np.where(dppon > 0)[0].size == 0:
            dcm_off = None
        else:
            query = 'select Timestamp,DCMoffset_attn from fV' \
                    +ver+'_vD50 where Timestamp >= '+tstart+' and Timestamp <= '+tend+' order by Timestamp'
            data, msg = db.do_query(cursor, query)
            if msg == 'Success':
                otimes = Time(data['Timestamp'].astype('int')[::50],format='lv')
                dcmoff = data['DCMoffset_attn']
                dcmoff.shape = (nt, 50)
                # We now have a time-history of offsets, at least some of which are non-zero.
                # Offsets by slot number do us no good, so we need to translate to band number.
                # Get fseqfile name at mean of timerange, from stateframe SQL database
                fseqfile = get_fseqfile(Time(int(np.mean(trange.lv)),format='lv')) 
                if fseqfile is None:
                    print 'Error: No active fseq file.'
                    dcm_off = None
                else:
                    # Get fseqfile from ACC and return bandlist
                    bandlist = fseqfile2bandlist(fseqfile)
                    # Use bandlist to covert nt x 50 array to nt x nbands array of DCM attn offsets
                    # Note that this assumes DCM offset is the same for any multiply-sampled bands
                    # in the sequence.
                    nbands = len(bandlist)
                    dcm_off = np.zeros((nt,nbands),float)
                    dcm_off[:,bandlist - 1] = dcmoff
                    # Put into canonical order [nband, nt]
                    dcm_off = dcm_off.T
                    if dt:
                        # If we want other than full cadence, find mean over dt measurements
                        new_nt = len(times)
                        dcm_off = dcm_off[:,:new_nt*dt]
                        dcm_off.shape = (nbands,dt,new_nt)
                        dcm_off = np.mean(dcm_off,1)
            else:
                print 'Error reading DCM attenuations:',msg
                dcm_off = None
    else:
        print 'Error reading DPPon state:',msg
        dcm_off = None
    cursor.close()
    # tdim is the number of dimensions in SQL tables, while nant is the number of antennas. These are both
    # 15 prior to V67 SQL tables, but tdim=16 (while nant remains 15) in V67, at least until we get 16 ants...
    # The statement below returns only the first nant attenuations.  Now we have 16 ants, but nant is set
    # accordingly so it will work in either case.
    return {'times':times,'hlev':hlev[:nant],'vlev':vlev[:nant],'dcmattn':dcmattn,'dcmoff':dcm_off}
    
def get_gain_state(trange, dt=None, relax=False):
    ''' Get all gain-state information for a given timerange.  Returns a dictionary
        with keys as follows:
        
        times:     A Time object containing the array of times, size (nt)
        h1:        The first HPol attenuator value for 15 antennas, size (nt, 15) 
        v1:        The first VPol attenuator value for 15 antennas, size (nt, 15) 
        h2:        The second HPol attenuator value for 15 antennas, size (nt, 15) 
        v2:        The second VPol attenuator value for 15 antennas, size (nt, 15)
        dcmattn:   The base DCM attenuations for nbands x 15 antennas x 2 Poln, size (34 or 52,30)
                      The order is Ant1 H, Ant1 V, Ant2 H, Ant2 V, etc.
        dcmoff:    If DPPoffset-on is 0, this is None (meaning there are no changes to the
                      above base attenuations).  
                   If DPPoffset-on is 1, then dcmoff is a table of offsets to the 
                      base attenuation, size (nt, 50).  The offset applies to all 
                      antennas/polarizations.
                      
        Optional keywords:
           dt      Seconds between entries to read from SQL stateframe database. 
                     If omitted, 1 s is assumed.
           relax   Used for gain of reference time, in case there are no SQL data for the
                     requested time.  In that case it finds the data for the nearest later time.
    '''
    if dt is None:
        tstart,tend = [str(i) for i in trange.lv]
    else:
        # Expand time by 1/2 of dt before and after
        tstart = str(np.round(trange[0].lv - dt/2))
        tend = str(np.round(trange[1].lv + dt/2))
    cursor = db.get_cursor()
    ver = db.find_table_version(cursor,trange[0].lv)
    if int(ver) > 66:
        tdim = 16   # Antenna info is in dimension 16 table starting with version 67
    else:
        tdim = 15
    if trange[0] < Time('2025-05-22'):
        nant = 15
    else:
        nant = 16

    # Get front end attenuator states
    # Attempt to solve the problem if there are no data 
    if relax:
        # Special case of reference gain, where we want the first nt records after tstart, in case there
        # are no data at time tstart
        nt = int(float(tend) - float(tstart) - 1)*tdim
        query = 'select top '+str(nt)+' Timestamp,Ante_Fron_FEM_HPol_Atte_First,Ante_Fron_FEM_HPol_Atte_Second,' \
            +'Ante_Fron_FEM_VPol_Atte_First,Ante_Fron_FEM_VPol_Atte_Second,Ante_Fron_FEM_Clockms from fV' \
            +ver+'_VD'+str(tdim)+' where Timestamp >= '+tstart+' order by Timestamp'
    else:
        query = 'select Timestamp,Ante_Fron_FEM_HPol_Atte_First,Ante_Fron_FEM_HPol_Atte_Second,' \
            +'Ante_Fron_FEM_VPol_Atte_First,Ante_Fron_FEM_VPol_Atte_Second,Ante_Fron_FEM_Clockms from fV' \
            +ver+'_VD'+str(tdim)+' where Timestamp >= '+tstart+' and Timestamp < '+tend+' order by Timestamp'
    #if dt:
    #    # If dt (seconds between measurements) is set, add appropriate SQL statement to query
    #    query += ' and (cast(Timestamp as bigint) % '+str(dt)+') = 0 '
    data, msg = db.do_query(cursor, query)
    if msg == 'Success':
        if dt:
            # If we want other than full cadence, get new array shapes and times
            n = len(data['Timestamp'])  # Original number of times
            new_n = (n/tdim/dt)*tdim*dt     # Truncated number of times equally divisible by dt
            new_shape = (n/tdim/dt,dt,tdim) # New shape of truncated arrays
            times = Time(data['Timestamp'][:new_n].astype('int')[::tdim*dt],format='lv')
        else:
            times = Time(data['Timestamp'].astype('int')[::tdim],format='lv')
        # Change tstart and tend to correspond to actual times from SQL
        tstart, tend = [str(i) for i in times[[0,-1]].lv]
        h1 = data['Ante_Fron_FEM_HPol_Atte_First']
        h2 = data['Ante_Fron_FEM_HPol_Atte_Second']
        v1 = data['Ante_Fron_FEM_VPol_Atte_First']
        v2 = data['Ante_Fron_FEM_VPol_Atte_Second']
        ms = data['Ante_Fron_FEM_Clockms']
        nt = len(h1)/tdim
        h1.shape = (nt,tdim)
        h2.shape = (nt,tdim)
        v1.shape = (nt,tdim)
        v2.shape = (nt,tdim)
        ms.shape = (nt,tdim)
        # Find any entries for which Clockms is zero, which indicates where no
        # gain-state measurement is available.
        for i in range(tdim):
            bad, = np.where(ms[:,i] == 0)
            if bad.size != 0 and bad.size != nt:
                # Find nearest adjacent good value
                good, = np.where(ms[:,i] != 0)
                idx = nearest_val_idx(bad,good)
                h1[bad,i] = h1[good[idx],i]
                h2[bad,i] = h2[good[idx],i]
                v1[bad,i] = v1[good[idx],i]
                v2[bad,i] = v2[good[idx],i]
        if dt:
            # If we want other than full cadence, find mean over dt measurements
            h1 = np.mean(h1[:new_n/tdim].reshape(new_shape),1)
            h2 = np.mean(h2[:new_n/tdim].reshape(new_shape),1)
            v1 = np.mean(v1[:new_n/tdim].reshape(new_shape),1)
            v2 = np.mean(v2[:new_n/tdim].reshape(new_shape),1)
        # Put results in canonical order [nant, nt]
        h1 = h1.T
        h2 = h2.T
        v1 = v1.T
        v2 = v2.T
    else:
        print 'Error reading FEM attenuations:',msg
        return {}
    # Get back end attenuator states
    xml, buf = ch.read_cal(2, t=trange[0])
    dcmattn = stf.extract(buf,xml['Attenuation'])
    nbands = dcmattn.shape[0]
    dcmattn.shape = (nbands, dcmattn.shape[1]//2, 2)
    # Put into canonical order [nant, npol, nband]
    dcmattn = np.moveaxis(dcmattn,0,2)
    # See if DPP offset is enabled
    query = 'select Timestamp,DPPoffsetattn_on from fV' \
            +ver+'_vD1 where Timestamp >= '+tstart+' and Timestamp <= '+tend+'order by Timestamp'
    data, msg = db.do_query(cursor, query)
    if msg == 'Success':
        dppon = data['DPPoffsetattn_on']
        if np.where(dppon > 0)[0].size == 0:
            dcm_off = None
        else:
            query = 'select Timestamp,DCMoffset_attn from fV' \
                    +ver+'_vD50 where Timestamp >= '+tstart+' and Timestamp <= '+tend
            #if dt:
            #    # If dt (seconds between measurements) is set, add appropriate SQL statement to query
            #    query += ' and (cast(Timestamp as bigint) % '+str(dt)+') = 0 '
            query += ' order by Timestamp'
            data, msg = db.do_query(cursor, query)
            if msg == 'Success':
                otimes = Time(data['Timestamp'].astype('int')[::50],format='lv')
                dcmoff = data['DCMoffset_attn']
                dcmoff.shape = (nt, 50)
                # We now have a time-history of offsets, at least some of which are non-zero.
                # Offsets by slot number do us no good, so we need to translate to band number.
                # Get fseqfile name at mean of timerange, from stateframe SQL database
                fseqfile = get_fseqfile(Time(int(np.mean(trange.lv)),format='lv')) 
                if fseqfile is None:
                    print 'Error: No active fseq file.'
                    dcm_off = None
                else:
                    # Get fseqfile from ACC and return bandlist
                    bandlist = fseqfile2bandlist(fseqfile)
                    nbands = len(bandlist)
                    # Use bandlist to covert nt x 50 array to nt x nbands array of DCM attn offsets
                    # Note that this assumes DCM offset is the same for any multiply-sampled bands
                    # in the sequence.
                    dcm_off = np.zeros((nt,nbands),float)
                    dcm_off[:,bandlist - 1] = dcmoff
                    # Put into canonical order [nband, nt]
                    dcm_off = dcm_off.T
                    if dt:
                        # If we want other than full cadence, find mean over dt measurements
                        new_nt = len(times)
                        dcm_off = dcm_off[:,:new_nt*dt]
                        dcm_off.shape = (nbands,dt,new_nt)
                        dcm_off = np.mean(dcm_off,1)
            else:
                print 'Error reading DCM attenuations:',msg
                dcm_off = None
    else:
        print 'Error reading DPPon state:',msg
        dcm_off = None
    cursor.close()
    # tdim is the number of dimensions in SQL tables, while nant is the number of antennas. These are both
    # 15 prior to V67 SQL tables, but tdim=16 (while nant remains 15) in V67, at least until we get 16 ants...
    # The statement below returns only the first nant attenuations.  Now we have 16 ants, but nant is set
    # accordingly so it will work in either case.
    return {'times':times,'h1':h1[:nant],'v1':v1[:nant],'h2':h2[:nant],'v2':v2[:nant],'dcmattn':dcmattn,'dcmoff':dcm_off}

def apply_fem_level(data, skycal={}, gctime=None):
    ''' Applys the FEM level corrections to the given data dictionary.
        
        Inputs:
          data     A dictionary such as that returned by read_idb().
          skycal   A dictionary returned by skycal_anal() in calibration.py.  This is
                     used to subtract a small "receiver background" before scaling for
                     fem level, and then adding it back.
          gctime   A Time() object whose date specifies which GAINCALTEST
                     measurements to use.  If omitted, the date of the data
                     is used.

        Output:
          cdata    A dictionary with the level-corrected data.  The keys
                     p, x, p2, and a are all updated.
    '''
    from util import common_val_idx, nearest_val_idx
    import attncal as ac
    import copy

    # Get timerange from data
    trange = Time([data['time'][0],data['time'][-1]],format='jd')
    if trange[0] < Time('2025-05-22'):
        nsolant = 13
        nant = 15
    else:
        nsolant = 15
        nant = 16
    if gctime is None:
        gctime = trange[0]
    # Get time cadence
    dt = np.int(np.round(np.median(data['time'][1:] - data['time'][:-1]) * 86400))
    if dt == 1: dt = None
    cdata = copy.deepcopy(data)
    # Get the FEM levels of the requested timerange
    src_lev = get_fem_level(trange,dt)   # solar gain state for timerange of file
    if src_lev == {}:
        print 'APPLY_FEM_LEVEL: No GAINCALTEST scans for this date, so no FEM level correction applied.'
        return cdata

    nf = len(data['fghz'])
    nt = len(src_lev['times'])
    # First attempt to read from the SQL database.  If that fails, read from the IDB file itself
    try:
        attn = ac.read_attncal(gctime)[0]   # Attn from SQL
        if (gctime.mjd - attn['time'].mjd) > 1:
            # SQL entry is too old, so analyze the GAINCALTEST
            attn = ac.get_attncal(gctime)[0]   # Attn measured by GAINCALTEST (returns a list, but use first, generally only, one)
            if attn != {}:
                ch.fem_attn_val2sql([attn])   # Go ahead and write it to SQL
    except:
        attn = ac.get_attncal(gctime)[0]   # Attn measured by GAINCALTEST (returns a list, but use first, generally only, one)
    antgain = np.zeros((nant,2,nf,nt),np.float32)   # Antenna-based gains [dB] vs. frequency
    # Find common frequencies of attn with data
    idx1, idx2 = common_val_idx(data['fghz'],attn['fghz'],precision=4)
    # Currently, GAINCALTEST measures 8 levels of attenuation (16 dB).  I assumed this would be enough,
    # but the flare of 2017-09-10 actually went to 10 levels (20 dB), so we have no choice but to extend
    # to higher levels using only the nominal, 2 dB steps above the 8th level.  This part of the code
    # extends to the maximum 16 levels.
    a = np.zeros((16, nsolant, 2, nf), float)  # Extend attenuation to 14 levels
    a[1:9, :, :, idx1] = attn['attn'][:, :nsolant, :, idx2]  # Use GAINCALTEST results in levels 1-9 (bottom level is 0dB)
    for i in range(8, 15):
        # Extend to levels 9-15 by adding 2 dB to each previous level
        a[i + 1] = a[i] + 2.
    #a[15] = 62.  # Level 15 means 62 dB have been inserted.
    a[15] = 0.  # Do not try to correct for full attenuation--it just makes huge noise...
    if dt:
        # For this case, src_lev is an array of dictionaries where keys are levels and
        # values are the proportion of that level for the given integration
        for i in range(nsolant):
            for k,j in enumerate(idx1):
                for m in range(nt):
                    for lev, prop in src_lev['hlev'][i,m].items():
                        antgain[i,0,j,m] += prop*a[lev,i,0,idx2[k]]
                    for lev, prop in src_lev['vlev'][i,m].items():
                        antgain[i,1,j,m] += prop*a[lev,i,1,idx2[k]]
    else:
        # For this case, src_lev is just an array of levels
        for i in range(nsolant):
            for k,j in enumerate(idx1):
                antgain[i,0,j] = a[src_lev['hlev'][i],i,0,idx2[k]]
                antgain[i,1,j] = a[src_lev['vlev'][i],i,1,idx2[k]]
    blgain = np.zeros((120,4,nf,nt),float)     # Baseline-based gains vs. frequency
    # Loop over cross-correlations (no auto-correlations)
    for i in range(nant-1):
         for j in range(i+1,nant):
             blgain[ri.bl2ord[i,j],0] = 10**((antgain[i,0] + antgain[j,0])/20.)
             blgain[ri.bl2ord[i,j],1] = 10**((antgain[i,1] + antgain[j,1])/20.)
             blgain[ri.bl2ord[i,j],2] = 10**((antgain[i,0] + antgain[j,1])/20.)
             blgain[ri.bl2ord[i,j],3] = 10**((antgain[i,1] + antgain[j,0])/20.)
    antgainf = 10**(antgain/10.)

    #idx1, idx2 = common_val_idx(data['time'],src_gs['times'].jd)
    idx = nearest_val_idx(data['time'],src_lev['times'].jd)
    # Apply corrections (some times may be eliminated from the data)
    # Correct the cross-correlation data
    cdata['x'] *= blgain[:,:,:,idx]
    # If a skycal dictionary exists, subtract receiver noise before scaling
    # NB: This will break SK!
    if skycal != {}:
        sna, snp, snf = skycal['rcvr_bgd'].shape
        bgd = skycal['rcvr_bgd'].repeat(nt).reshape((sna,snp,snf,nt))
        bgd_auto = skycal['rcvr_bgd_auto'].repeat(nt).reshape((sna,snp,snf,nt))
        cdata['p'][:nsolant] -= bgd[:,:,:,idx]
        cdata['a'][:nsolant,:2] -= bgd_auto[:,:,:,idx]
    # Correct the power,
    cdata['p'][:nant] *= antgainf[:,:,:,idx]
    # Correct the autocorrelation
    cdata['a'][:nant,:2] *= antgainf[:,:,:,idx]
    # If a skycal dictionary exists, add back the receiver noise
    #if skycal:
    #    cdata['p'][:nsolant] += bgd[:,:,:,idx]
    #    cdata['a'][:nsolant,:2] += bgd_auto[:,:,:,idx]
    cross_fac = np.sqrt(antgainf[:,0]*antgainf[:,1])
    cdata['a'][:nant,2] *= cross_fac[:,:,idx]
    cdata['a'][:nant,3] *= cross_fac[:,:,idx]
    # Correct the power-squared -- this should preserve SK
    cdata['p2'][:nant] *= antgainf[:,:,:,idx]**2
    # Remove any uncorrected times before returning
    #cdata['time'] = cdata['time'][idx1]
    #cdata['p'] = cdata['p'][:,:,:,idx1]
    #cdata['a'] = cdata['a'][:,:,:,idx1]
    #cdata['p2'] = cdata['p2'][:,:,:,idx1]
    #cdata['ha'] = cdata['ha'][idx1]
    #cdata['m'] = cdata['m'][:,:,:,idx1]
    return cdata

def apply_gain_corr(data, tref=None):
    ''' Applys the gain_state() corrections to the given data dictionary,
        corrected to the gain-state at time given by Time() object tref.
        
        Inputs:
          data     A dictionary such as that returned by read_idb().
          tref     A Time() object with the reference time, or if None,
                     the gain state of the nearest earlier REFCAL is 
                     used.
        Output:
          cdata    A dictionary with the gain-corrected data.  The keys
                     p, x, p2, and a are all updated.
    '''
    from util import fname2mjd, common_val_idx, nearest_val_idx
    import copy
    if tref is None:
        # No reference time specified, so get nearest PHASECAL scan (should guarantee femauto-off state)
        mjd = int(Time(data['time'][0],format='jd').mjd)
        trange = Time([mjd+10./24.,mjd+1.], format='mjd')
        import pcal_anal as pa
        try:
            # Get filename of first PHASECAL on this day, and use it to get reference time
            scanfile = pa.findfile(trange)['scanlist'][0][0]
            tref = Time(fname2mjd(scanfile)+60./86400,format='mjd')  # Add one minute to ensure scan is active
            # Get the gain state at the reference time (actually median over 1 minute)
            trefrange = Time([tref.iso,Time(tref.lv+61,format='lv').iso])
            ref_gs =  get_gain_state(trefrange)  # refcal gain state for 60 s
        except:
            # No phasecal in timerange, so just use an early time as reference
            tref = Time(trange[0].iso[:10]+' 13:30')
            # Get the gain state at the reference time (actually median over 1 minute)
            trefrange = Time([tref.iso,Time(tref.lv+61,format='lv').iso])
            ref_gs =  get_gain_state(trefrange, relax=True)  # refcal gain state for 60 s after ref time
        if trange[0].mjd - tref.mjd > 2:
            # Reference calibration is too old, so just use an early time as reference
            tref = Time(trange[0].iso[:10]+' 13:30')
            # Get the gain state at the reference time (actually median over 1 minute)
            trefrange = Time([tref.iso,Time(tref.lv+61,format='lv').iso])
            ref_gs =  get_gain_state(trefrange, relax=True)  # refcal gain state for 60 s after ref time
    if tref < Time('2025-05-22'):
        nant = 15
    else:
        nant = 16
    # Get median of refcal gain state (which should be constant anyway)
    ref_gs['h1'] = np.median(ref_gs['h1'],1)
    ref_gs['h2'] = np.median(ref_gs['h2'],1)
    ref_gs['v1'] = np.median(ref_gs['v1'],1)
    ref_gs['v2'] = np.median(ref_gs['v2'],1)

    # Get timerange from data
    trange = Time([data['time'][0],data['time'][-1]],format='jd')
    # Get time cadence
    dt = np.int(np.round(np.median(data['time'][1:] - data['time'][:-1]) * 86400))
    if dt == 1: dt = None
    # Get the gain state of the requested timerange
    src_gs = get_gain_state(trange,dt)   # solar gain state for timerange of file
    nt = len(src_gs['times'])
    nbands = src_gs['dcmattn'].shape[2]
    if nbands != ref_gs['dcmattn'].shape[2]:
        # Reference gain state is incompatible with this one, so set to src gain state (no correction will be applied)
        ref_gs = src_gs
        print 'GAINCAL2 Warning: Data and reference gain states are not compatible. No correction applied!'
    antgain = np.zeros((nant,2,nbands,nt),np.float32)   # Antenna-based gains vs. band
    for i in range(nant):
        for j in range(nbands):
            antgain[i,0,j] = src_gs['h1'][i] + src_gs['h2'][i] - ref_gs['h1'][i] - ref_gs['h2'][i] + src_gs['dcmattn'][i,0,j] - ref_gs['dcmattn'][i,0,j]
            antgain[i,1,j] = src_gs['v1'][i] + src_gs['v2'][i] - ref_gs['v1'][i] - ref_gs['v2'][i] + src_gs['dcmattn'][i,1,j] - ref_gs['dcmattn'][i,1,j]

    cdata = copy.deepcopy(data)
    # Frequency list is provided, so produce baseline-based gain table as well
    # Create giant array of gains, translated to baselines and frequencies
    fghz = data['fghz']
    nf = len(fghz)
    blist = (fghz*2 - 1).astype(int) - 1
    blgain = np.zeros((120,4,nf,nt),float)     # Baseline-based gains vs. frequency
    for i in range(nant-1):
         for j in range(i+1,nant):
             blgain[ri.bl2ord[i,j],0] = 10**((antgain[i,0,blist] + antgain[j,0,blist])/20.)
             blgain[ri.bl2ord[i,j],1] = 10**((antgain[i,1,blist] + antgain[j,1,blist])/20.)
             blgain[ri.bl2ord[i,j],2] = 10**((antgain[i,0,blist] + antgain[j,1,blist])/20.)
             blgain[ri.bl2ord[i,j],3] = 10**((antgain[i,1,blist] + antgain[j,0,blist])/20.)
    antgainf = 10**(antgain[:,:,blist]/10.)

    #idx1, idx2 = common_val_idx(data['time'],src_gs['times'].jd)
    idx = nearest_val_idx(data['time'],src_gs['times'].jd)
    # Apply corrections (some times may be eliminated from the data)
    # Correct the cross-correlation data
    cdata['x'] *= blgain[:,:,:,idx]
    # Correct the power
    cdata['p'][:nant] *= antgainf[:,:,:,idx]
    # Correct the autocorrelation
    cdata['a'][:nant,:2] *= antgainf[:,:,:,idx]
    cross_fac = np.sqrt(antgainf[:,0]*antgainf[:,1])  # Gains from one ant, two poln
    cdata['a'][:nant,2] *= cross_fac[:,:,idx]
    cdata['a'][:nant,3] *= cross_fac[:,:,idx]
    # Correct the power-squared -- this should preserve SK
    cdata['p2'][:nant] *= antgainf[:,:,:,idx]**2
    # Remove any uncorrected times before returning
    #cdata['time'] = cdata['time'][idx1]
    #cdata['p'] = cdata['p'][:,:,:,idx1]
    #cdata['a'] = cdata['a'][:,:,:,idx1]
    #cdata['p2'] = cdata['p2'][:,:,:,idx1]
    #cdata['ha'] = cdata['ha'][idx1]
    #cdata['m'] = cdata['m'][:,:,:,idx1]
    return cdata
    
def get_gain_corr(trange, tref=None, fghz=None):
    ''' Calls get_gain_state() for a timerange and a reference time,
        and returns the gain difference table to apply to data in the
        given timerange.  If no reference time is provided, the gain
        state is referred to the nearest earlier REFCAL.
        
        Returns a dictionary containing:
          antgain    Array of size (nant, npol, nbands, nt)
          times      A Time() object corresponding to the times in 
                       antgain
    '''
    if tref is None:
        # No reference time specified, so get nearest earlier REFCAL
        xml, buf = ch.read_cal(8,t=trange[0])
        tref = Time(stf.extract(buf,xml['Timestamp']),format='lv')
    # Get the gain state at the reference time (actually median over 1 minute)
    trefrange = Time([tref.iso,Time(tref.lv+61,format='lv').iso])
    ref_gs =  get_gain_state(trefrange)  # refcal gain state for 60 s
    # Get median of refcal gain state (which should be constant anyway)
    ref_gs['h1'] = np.median(ref_gs['h1'],1)
    ref_gs['h2'] = np.median(ref_gs['h2'],1)
    ref_gs['v1'] = np.median(ref_gs['v1'],1)
    ref_gs['v2'] = np.median(ref_gs['v2'],1)
    if trange[0] > Time('2025-05-22'):
        nant = 16
    else:
        nant = 15

    # Get the gain state of the requested timerange
    src_gs = get_gain_state(trange)   # solar gain state for timerange of file
    nt = len(src_gs['times'])
    nbands = src_gs['dcmattn'].shape[2]
    antgain = np.zeros((nant,2,nbands,nt),np.float32)   # Antenna-based gains vs. band
    for i in range(nant):
        for j in range(nbands):
            antgain[i,0,j] = src_gs['h1'][i] + src_gs['h2'][i] - ref_gs['h1'][i] - ref_gs['h2'][i] + src_gs['dcmattn'][i,0,j] - ref_gs['dcmattn'][i,0,j]
            antgain[i,1,j] = src_gs['v1'][i] + src_gs['v2'][i] - ref_gs['v1'][i] - ref_gs['v2'][i] + src_gs['dcmattn'][i,1,j] - ref_gs['dcmattn'][i,1,j]

    return {'antgain': antgain, 'times': src_gs['times']}

def eq_anal(out):
    ''' Analyze 3 consecutive GAINCALTEST observations taken with EQ settings 8, 4 and 2, to
        determine the best settings for each antenna/polarization/frequency.
    '''
    import chan_util_52 as cu
    blah = cu.freq2bdname(out['fghz'])
    if out['time'][0] < Time('2025-05-22').jd:
        nsolant = 13
    else:
        nsolant = 15
    nfpb = []
    for i in range(1,53):
        nfpb.append(len(where(blah == i)[0]))
    nfb = (array(nfpb)[blah-1]/8.)
    nf = len(out['fghz'])
    ratio = np.zeros((3,nsolant,2,nf))
    idx = np.array([30,40])
    for k in range(nsolant):
        for j in range(2):
            ratio[0,k,j,:] = np.mean(out['p'][k,j,:,    idx]/np.abs(out['a'][k,j,:,    idx]),0)*8**2*nfb
            ratio[1,k,j,:] = np.mean(out['p'][k,j,:,115+idx]/np.abs(out['a'][k,j,:,115+idx]),0)*4**2*nfb
            ratio[2,k,j,:] = np.mean(out['p'][k,j,:,230+idx]/np.abs(out['a'][k,j,:,230+idx]),0)*2**2*nfb
    foo = abs(ratio - 400).argmin(axis=0)
    # This routine appears to be unfinished, and is certainly not used...