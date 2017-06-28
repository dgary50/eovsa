# gaincal2
# Modification History
#  2017-05-13  DG
#    First wrote routine get_gain_state()
#  2017-05-15  DG
#    Added get_fseqfile() and fseqfile2bandlist() to help with conversion
#    of slots to bands.  Also ensure that all arrays returned by get_gain_state() 
#    are in canonical order [nant, npol, nf/nband, nt].
#  2017-05-21  DG
#    Some changes to apply_gain_corr() to make it more general.
#  2017-06-12  DG
#    Important change to apply_gain_corr() so that it does not drop unmeasured
#    points.  Instead, it uses the gain state of the nearest time in case of a missing one
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
           bandlist    numpy 50-element integer array of band numbers, 1-34
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
    
def get_gain_state(trange):
    ''' Get all gain-state information for a given timerange.  Returns a dictionary
        with keys as follows:
        
        times:     A Time object containing the array of times, size (nt)
        h1:        The first HPol attenuator value for 15 antennas, size (nt, 15) 
        v1:        The first VPol attenuator value for 15 antennas, size (nt, 15) 
        h2:        The second HPol attenuator value for 15 antennas, size (nt, 15) 
        v2:        The second VPol attenuator value for 15 antennas, size (nt, 15)
        dcmattn:   The base DCM attenuations for 34 bands x 15 antennas x 2 Poln, size (34,30)
                      The order is Ant1 H, Ant1 V, Ant2 H, Ant2 V, etc.
        dcmoff:    If DPPoffset-on is 0, this is None (meaning there are no changes to the
                      above base attenuations).  
                   If DPPoffset-on is 1, then dcmoff is a table of offsets to the 
                      base attenuation, size (nt, 50).  The offset applies to all 
                      antennas/polarizations.
    '''
    tstart,tend = [str(i) for i in trange.lv]
    cursor = db.get_cursor()
    ver = db.find_table_version(cursor,trange[0].lv)
    # Get front end attenuator states
    query = 'select Timestamp,Ante_Fron_FEM_HPol_Atte_First,Ante_Fron_FEM_HPol_Atte_Second,' \
            +'Ante_Fron_FEM_VPol_Atte_First,Ante_Fron_FEM_VPol_Atte_Second,Ante_Fron_FEM_Clockms from fV' \
            +ver+'_vD15 where Timestamp > '+tstart+' and Timestamp < '+tend+'order by Timestamp'
    data, msg = db.do_query(cursor, query)
    if msg == 'Success':
        times = Time(data['Timestamp'].astype('int'),format='lv')[::15]
        h1 = data['Ante_Fron_FEM_HPol_Atte_First']
        h2 = data['Ante_Fron_FEM_HPol_Atte_Second']
        v1 = data['Ante_Fron_FEM_VPol_Atte_First']
        v2 = data['Ante_Fron_FEM_VPol_Atte_Second']
        ms = data['Ante_Fron_FEM_Clockms']
        nt = len(h1)/15
        h1.shape = (nt,15)
        h2.shape = (nt,15)
        v1.shape = (nt,15)
        v2.shape = (nt,15)
        ms.shape = (nt,15)
        # Find any entries for which Clockms is zero, which indicates where no
        # gain-state measurement is available.
        for i in range(15):
            bad, = np.where(ms[:,i] == 0)
            if bad.size != 0 and bad.size != nt:
                # Find nearest adjacent good value
                good, = np.where(ms[:,i] != 0)
                idx = nearest_val_idx(bad,good)
                h1[bad,i] = h1[good[idx],i]
                h2[bad,i] = h2[good[idx],i]
                v1[bad,i] = v1[good[idx],i]
                v2[bad,i] = v2[good[idx],i]
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
    dcmattn.shape = (34, 15, 2)
    # Put into canonical order [nant, npol, nband]
    dcmattn = np.moveaxis(dcmattn,0,2)
    # See if DPP offset is enabled
    query = 'select Timestamp,DPPoffsetattn_on from fV' \
            +ver+'_vD1 where Timestamp > '+tstart+' and Timestamp < '+tend+'order by Timestamp'
    data, msg = db.do_query(cursor, query)
    if msg == 'Success':
        dppon = data['DPPoffsetattn_on']
        if np.where(dppon > 0)[0].size == 0:
            dcm_off = None
        else:
            query = 'select Timestamp,DCMoffset_attn from fV' \
                    +ver+'_vD50 where Timestamp > '+tstart+' and Timestamp < '+tend+'order by Timestamp'
            data, msg = db.do_query(cursor, query)
            if msg == 'Success':
                otimes = Time(data['Timestamp'].astype('int'),format='lv')[::15]
                dcmoff = data['DCMoffset_attn']
                dcmoff.shape = (nt, 50)
                # We now have a time-history of offsets, at least some of which are non-zero.
                # Offsets by slot number do us no good, so we need to translate to band number.
                # Get fseqfile name at mean of timerange, from stateframe SQL database
                fseqfile = get_fseqfile(Time(int(np.mean(trange.lv)),format='lv')) 
                if fseqfile is None:
                    return {}
                else:
                    # Get fseqfile from ACC and return bandlist
                    bandlist = fseqfile2bandlist(fseqfile)
                    # Use bandlist to covert nt x 50 array to nt x 34 band array of DCM attn offsets
                    # Note that this assumes DCM offset is the same for any multiply-sampled bands
                    # in the sequence.
                    dcm_off = np.zeros((nt,34),float)
                    dcm_off[:,bandlist - 1] = dcmoff
                    # Put into canonical order [nband, nt]
                    dcm_off = dcm_off.T
            else:
                print 'Error reading DCM attenuations:',msg
                return {}
    else:
        print 'Error reading DPPon state:',msg
        return {}
    cursor.close()
    return {'times':times,'h1':h1,'v1':v1,'h2':h2,'v2':v2,'dcmattn':dcmattn,'dcmoff':dcm_off}

def apply_gain_corr(data, tref=None):
    ''' Applys the gain_state() corrections to the given data dictionary,
        corrected to the gain-state at time given by Time() object tref.
        
        Inputs:
          data     A dictionary such as that returned by read_idb().
          tref     A Time() object with the reference time, or if None,
                     the gain state of the nearest earlier SOLPNTCAL is 
                     used.
        Output:
          cdata    A dictionary with the gain-corrected data.  The keys
                     p, x, p2, and a are all updated.
    '''
    from util import common_val_idx, nearest_val_idx
    import copy
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

    # Get timerange from data
    trange = Time([data['time'][0],data['time'][-1]],format='jd')
    # Get the gain state of the requested timerange
    src_gs = get_gain_state(trange)   # solar gain state for timerange of file
    nt = len(src_gs['times'])
    antgain = np.zeros((15,2,34,nt),np.float32)   # Antenna-based gains vs. band
    for i in range(15):
        for j in range(34):
            antgain[i,0,j] = src_gs['h1'][i] + src_gs['h2'][i] - ref_gs['h1'][i] - ref_gs['h2'][i] + src_gs['dcmattn'][i,0,j] - ref_gs['dcmattn'][i,0,j]
            antgain[i,1,j] = src_gs['v1'][i] + src_gs['v2'][i] - ref_gs['v1'][i] - ref_gs['v2'][i] + src_gs['dcmattn'][i,1,j] - ref_gs['dcmattn'][i,1,j]

    cdata = copy.deepcopy(data)
    # Frequency list is provided, so produce baseline-based gain table as well
    # Create giant array of gains, translated to baselines and frequencies
    fghz = data['fghz']
    nf = len(fghz)
    blist = (fghz*2 - 1).astype(int) - 1
    blgain = np.zeros((120,4,nf,nt),float)     # Baseline-based gains vs. frequency
    for i in range(14):
         for j in range(i+1,15):
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
    cdata['p'][:15] *= antgainf[:,:,:,idx]
    # Correct the autocorrelation
    cdata['a'][:15,:2] *= antgainf[:,:,:,idx]
    cross_fac = np.sqrt(antgainf[:,0]*antgainf[:,1])
    cdata['a'][:15,2] *= cross_fac[:,:,idx]
    cdata['a'][:15,3] *= cross_fac[:,:,idx]
    # Correct the power-squared -- this should preserve SK
    cdata['p2'][:15] *= antgainf[:,:,:,idx]**2
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
          antgain    Array of size (15, 2, 34, nt) = (nant, npol, nbands, nt)
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

    # Get the gain state of the requested timerange
    src_gs = get_gain_state(trange)   # solar gain state for timerange of file
    nt = len(src_gs['times'])
    antgain = np.zeros((15,2,34,nt),np.float32)   # Antenna-based gains vs. band
    for i in range(15):
        for j in range(34):
            antgain[i,0,j] = src_gs['h1'][i] + src_gs['h2'][i] - ref_gs['h1'][i] - ref_gs['h2'][i] + src_gs['dcmattn'][i,0,j] - ref_gs['dcmattn'][i,0,j]
            antgain[i,1,j] = src_gs['v1'][i] + src_gs['v2'][i] - ref_gs['v1'][i] - ref_gs['v2'][i] + src_gs['dcmattn'][i,1,j] - ref_gs['dcmattn'][i,1,j]

    return {'antgain': antgain, 'times': src_gs['times']}
