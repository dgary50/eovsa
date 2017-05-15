# gaincal2
# Modification History
#  2017-05-13  DG
#    First wrote routine get_gain_state()
#  2017-05-15  DG
#    Added get_fseqfile() and fseqfile2bandlist() to help with conversion
#    of slots to bands.  Also ensure that all arrays returned by get_gain_state() 
#    are in canonical order [nant, npol, nf/nband, nt].
#
import dbutil as db
import cal_header as ch
from util import Time
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
    query = 'select Timestamp,Ante_Fron_FEM_HPol_Atte_First,Ante_Fron_FEM_HPol_Atte_Second,Ante_Fron_FEM_VPol_Atte_First,Ante_Fron_FEM_VPol_Atte_Second from fV' \
            +ver+'_vD15 where Timestamp > '+tstart+' and Timestamp < '+tend+'order by Timestamp'
    data, msg = db.do_query(cursor, query)
    if msg == 'Success':
        times = Time(data['Timestamp'].astype('int'),format='lv')[::15]
        h1 = data['Ante_Fron_FEM_HPol_Atte_First']
        h2 = data['Ante_Fron_FEM_HPol_Atte_Second']
        v1 = data['Ante_Fron_FEM_VPol_Atte_First']
        v2 = data['Ante_Fron_FEM_VPol_Atte_Second']
        nt = len(h1)/15
        h1.shape = (nt,15)
        h2.shape = (nt,15)
        v1.shape = (nt,15)
        v2.shape = (nt,15)
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
