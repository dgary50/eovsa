#
# Routines to set up system for ADC calibration/DCM attenuation setting
#
#   2016-Feb-20  DG
#     First written.
#   2016-Mar-19  DG
#     Converted to adc_cal2.py (for now), and added set_fem_attn() routine.
#   2016-May-05  DG
#     Added a new routine DCM_master_attn_cal() to do the quickest possible
#     calibration--takes only 5 minutes.  It returns DCM_lines, ready to be
#     inserted into the SQL table by
#        import cal_header 
#        cal_header.dcm_master_table2sql(DCM_lines)
#   2016-May-21  DG
#     Quite a few changes in an attempt to get set_dcm_attn() to work.
#   2016-Aug-01  DG
#     Important change to DCM_master_attn_cal() (completely rewritten) to
#     use a new scheme involving capture of packets on the dpp.  This
#     routine can only be run on the dpp...
#   2016-Aug-02  DG
#     Added a new routine, DCM_attn_anal(), which analyzes an observation
#     made using DCMATTNTEST.ctl.
#   2016-Aug-06  DG
#     Added gain_state() routine, which returns the FEM and DCM attenuations
#     for a given timerange, as complex arrays in dB units.
#   2016-Aug-07  DG
#     Added output of timestamp array (in Julian date) to gain_state(). Also
#     added routines and changed DCM_master_attn_cal() to ensure that the
#     desired tuning sequence and DCM-switching is active.
#

import time
import numpy as np
import roach as r
import stateframe as stf

def acc_tune(band,acc):
    if type(band) is int:
        fsqfile = 'BAND'+str(band)+'.FSQ'
    elif type(band) is str:
        if band.lower() == 'solar.fsq' or band.lower() == 'pcal.fsq':
            fsqfile = band.lower()
    else:
        print 'Error: Unknown band',band
        return
    cmds = ['FSEQ-OFF','FSEQ-INIT','WAIT','FSEQ-FILE '+fsqfile.lower(), 'FSEQ-ON']
    send_cmds(cmds,acc)

def send_cmds(cmds,acc):
    ''' Sends a series of commands to ACC.  The sequence of commands
        is not checked for validity!
        
        cmds   a list of strings, each of which must be a valid command
    '''
    import socket, stateframe

    for cmd in cmds:
        #print 'Command:',cmd
        s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        try:
            s.connect((acc['host'],acc['scdport']))
            s.send(cmd)
            time.sleep(0.01)
            s.close()
        except:
            print 'Error: Could not send command',cmd,' to ACC.'
    return

def ant_str2list(ant_str):
    ant_list = []
    try:
        grps = ant_str.split()
        for grp in grps:
            antrange = grp[3:].split('-')
            if len(antrange) == 1:
                if antrange != '':
                    ant_list.append(int(antrange[0])-1)
            elif len(antrange) == 2:
                ant_list += range(int(antrange[0])-1,int(antrange[1]))
    except:
        print 'Error: cannot interpret ant_str',ant_str
        return None
    return np.array(ant_list)
                    
def set_fem_attn(level=3,ant_str='ant1-15'):
    ''' Read current power and attenuation values in FEMs, and set the attenuation
        such that the power level will be as close to "level" as possible.
    '''
    accini = stf.rd_ACCfile()
    acc = {'host': accini['host'], 'scdport':accini['scdport']}
    
    ant_list = ant_str2list(ant_str)
    if ant_list is None:
        return 'Bad antenna list'
    accini = stf.rd_ACCfile()
    hatn1 = np.zeros((10,15),dtype='int')
    vatn1 = np.zeros((10,15),dtype='int')
    hatn2 = np.zeros((10,15),dtype='int')
    vatn2 = np.zeros((10,15),dtype='int')
    hpwr = np.zeros((10,15),dtype='float')
    vpwr = np.zeros((10,15),dtype='float')
    # Read 10 instances of attenuation and power, and take the median to avoid
    # glitches
    for i in range(10):
        # Read the frontend attenuations and powers for each antenna
        data, msg = stf.get_stateframe(accini)
        for iant in range(15):
            fem = accini['sf']['Antenna'][iant]['Frontend']['FEM']
            hatn1[i,iant] = stf.extract(data,fem['HPol']['Attenuation']['First'])
            vatn1[i,iant] = stf.extract(data,fem['VPol']['Attenuation']['First'])
            hatn2[i,iant] = stf.extract(data,fem['HPol']['Attenuation']['Second'])
            vatn2[i,iant] = stf.extract(data,fem['VPol']['Attenuation']['Second'])
            hpwr[i,iant] = stf.extract(data,fem['HPol']['Power'])
            vpwr[i,iant] = stf.extract(data,fem['VPol']['Power'])
        time.sleep(1)
    hatn1 = np.median(hatn1,0).astype('int')
    vatn1 = np.median(vatn1,0).astype('int')
    hatn2 = np.median(hatn2,0).astype('int')
    vatn2 = np.median(vatn2,0).astype('int')
    hpwr = np.median(hpwr,0)
    vpwr = np.median(vpwr,0)
    hatn2 = np.clip(hatn2 - (level - hpwr).astype('int'),0,31)
    vatn2 = np.clip(vatn2 - (level - vpwr).astype('int'),0,31)
    # Send attenuation commands to each antenna in ant_list
    for iant in ant_list:
        hatn = str(hatn1[iant])+' '+str(hatn2[iant])+' ant'+str(iant+1)
        vatn = str(vatn1[iant])+' '+str(vatn2[iant])+' ant'+str(iant+1)
        send_cmds(['HATTN '+hatn,'VATTN '+vatn],acc)
    return 'Success'

def chk_lo1a(accini, band, iteration=1):
    data, msg = stf.get_stateframe(accini)
    errstr = stf.extract(data,accini['sf']['LODM']['LO1A']['ERR']).split('"')[1]
    if iteration == 1 and errstr != 'No error':
        # Looks like a reboot of LO1A is needed!
        print '10-s delay while attempting to reboot LO1A'
        send_cmds(['LO1A-REBOOT'],acc)
        time.sleep(10)
        acc = {'host': accini['host'], 'scdport':accini['scdport']}
        acc_tune(band+1,acc)
        time.sleep(5)
        errstr = chk_lo1a(accini,band,iteration=2)
        if errstr != 'No error':
            errstr = 'Reboot attempt failed.'
    return errstr
    
    
def set_dcm_attn(roach_list,fem_level=5,nd_state='on',adc_nom=25,ant_list='ant1-13',do_plot=False):
    ''' Set the FEM attenuation to the given value, switch ND state to the given state
        and cycle through the IF bands, and find the optimum DCM attenuation needed to 
        attain the given ADC level.
        
        Returns the corresponding DCM table
    '''
    import copy
    accini = stf.rd_ACCfile()
    acc = {'host': accini['host'], 'scdport':accini['scdport']}
    
    # Switch noise diode to requested state
    if nd_state == 'on':
        send_cmds(['ND-ON '+ant_list],acc)
    else:
        send_cmds(['ND-OFF '+ant_list],acc)
    # Set FEM power level to requested level
    if set_fem_attn(fem_level,ant_list) == 'Failure':
        return
    time.sleep(5)
    dcm_table = np.zeros((34,32), dtype='int')
    # Set DCM state to standard values
    send_cmds(['DCMAUTO-OFF '+ant_list],acc)
    time.sleep(1)
    send_cmds(['DCMATTN 12 12 '+ant_list],acc)
    time.sleep(1)
    # Cycle through bands to get ADC levels
    for band in range(34):
        # Start with nominal DCM attenuation
        #send_cmds(['DCMATTN 12 12 '+ant_list],acc)
        #time.sleep(1)
        print 'Band:',band+1
        acc_tune(band+1,acc)
        time.sleep(5)
        errstr = chk_lo1a(accini,band)
        if errstr != 'No error':
            print errstr
            return None
        # Go through ROACH list
        for i,ro in enumerate(roach_list):
            # Get ADC levels at nominal setting
            dcm_base = np.array([12,12,12,12],dtype='int')
            r.adc_levels([ro])
            # Calculate new attenuation to achieve nominal level
            ch_atn = np.clip(((20*np.log10(ro.adc_levels/adc_nom)+dcm_base + 1)/2).astype('int')*2,0,30)
            # Set new attenuation levels (twice, for good measure, and check result
            # send_cmds(['DCMATTN '+str(ch_atn[0])+' '+str(ch_atn[1])+' ant'+str(ro.ants[0])],acc)
            # time.sleep(1)
            # send_cmds(['DCMATTN '+str(ch_atn[0])+' '+str(ch_atn[1])+' ant'+str(ro.ants[0])],acc)
            # time.sleep(1)
            # send_cmds(['DCMATTN '+str(ch_atn[2])+' '+str(ch_atn[3])+' ant'+str(ro.ants[1])],acc)
            # time.sleep(1)
            # send_cmds(['DCMATTN '+str(ch_atn[2])+' '+str(ch_atn[3])+' ant'+str(ro.ants[1])],acc)
            # time.sleep(1)
            # r.adc_levels([ro])
            # ch_atn2 = np.clip(((20*np.log10(ro.adc_levels/adc_nom)+ch_atn + 1)/2).astype('int')*2,0,30)
            # print '  ',ro.roach_ip[:6],'Attn:',ch_atn,'Check:',ch_atn2
            print '  ',ro.roach_ip[:6],'Attn:',ch_atn
            dcm_table[band,np.array(((ro.ants[0]-1)*2,(ro.ants[0]-1)*2+1,(ro.ants[1]-1)*2,(ro.ants[1]-1)*2+1))] = copy.copy(ch_atn)
    return dcm_table
    
# Insert 62 dB into FEMs, cycle through bands, get ADC levels (optionally plot results)
# Set FEM power between 2 and 3 dBm with ND-ON, cycle through bands, get ADC levels (optionally plot results)
# Use ADC level tests to "guess" best DCM attenuation settings

def adc_cal(roach_list,ant_list='ant1-15',do_plot=False):
    ''' Perform a sequence of FEM settings, using ADC levels to 
        deduce optimum DCM attenuation settings for all 34 bands.
        This can also reveal problems in FEM or DCM hardware.
        TAKES ABOUT 17 MINUTES TO COMPLETE
        
        roach_list  a set of ROACH objects created with roach.py
        ant_list    a list of antennas in the form of a string,
                      e.g. "ant1-5 ant7" on which to adjust FEMs
                      Default is all antennas, and an empty string
                      means all antennas in current subarray.
        do_plot     if True, makes a summary plot of results
        
        Returns numpy arrays :
                adc_nosig[34, nroach, 4] (no-signal ADC levels)
                adc_ndoff[34, nroach, 4] (ADC levels for ND-OFF)
                adc_ndon [34, nroach, 4] (ADC levels for ND-ON)
    '''
    accini = stf.rd_ACCfile()
    acc = {'host': accini['host'], 'scdport':accini['scdport']}

    n = len(roach_list)
    adc_nosig = np.zeros((34,n,4),dtype='float')
    adc_ndoff = np.zeros((34,n,4),dtype='float')
    adc_ndon = np.zeros((34,n,4),dtype='float')
    # Set DCM state to standard values
    send_cmds(['DCMAUTO-OFF '+ant_list,'DCMATTN 12 12 '+ant_list],acc)
    # Set FEM attenuations to maximum
    send_cmds(['FEMATTN 15 '+ant_list],acc)
    # Cycle through bands to get "zero-input" ADC levels
    for band in range(34):
        acc_tune(band+1,acc)
        time.sleep(3)
        # Go through roaches and find optimum ADC levels
        for i,ro in enumerate(roach_list):
            r.adc_levels([ro])
            
            adc_nosig[band,i] = ro.adc_levels
    # Set FEM attenuations to nominal
    send_cmds(['FEMATTN 0 '+ant_list],acc)
    # Cycle through bands to get "nd-on" ADC levels
    send_cmds(['ND-ON '+ant_list],acc)
    for band in range(34):
        acc_tune(band+1,acc)
        time.sleep(3)
        r.adc_levels(roach_list)
        for i,ro in enumerate(roach_list):
            adc_ndon[band,i] = ro.adc_levels
    # Cycle through bands to get "nd-off" ADC levels
    send_cmds(['ND-OFF '+ant_list],acc)
    for band in range(34):
        acc_tune(band+1,acc)
        time.sleep(3)
        r.adc_levels(roach_list)
        for i,ro in enumerate(roach_list):
            adc_ndoff[band,i] = ro.adc_levels
    if do_plot: 
        plot_adc_cal(roach_list, adc_nosig, adc_ndoff, adc_ndon)
    return adc_nosig, adc_ndoff, adc_ndon

def fseq_is_running(fseqfile,accini=None):
    ''' Check current stateframe to see if the given sequence file is
        running.  Returns True if so, False otherwise
    '''
    import stateframe as stf
    if accini is None:
        accini = stf.rd_ACCfile()
    # Make sure this sequence is actually running, or start it if not
    buf, msg = stf.get_stateframe(accini)
    if msg != 'No Error':
        print 'Error reading stateframe.'
        return None
    fseq = stf.extract(buf,accini['sf']['LODM']['LO1A']['FSeqFile'])
    fseq = fseq.strip('\x00')  # strip nulls from name
    result = fseq == fseqfile and stf.extract(buf,accini['sf']['LODM']['LO1A']['SweepStatus']) == 8
    return result

def bandlist2dcmtable(bandlist):
    '''Use list of bands representing a frequency sequence, to set 
       dcmtable.txt from the DCM_master_table, and send to ACC
       The sequence numbers are 1-based band numbers
    '''
    import stateframe as stf
    import cal_header as ch
    from ftplib import FTP
    # Convert from 1-based bandlist to zero-based band numbers
    bands = bandlist-1
    # Read master table from SQL server
    dcm, buf = ch.read_cal(2)
    dcm_m_attn = stf.extract(buf,dcm['Attenuation'])
    dcm_attn = dcm_m_attn[bands]
    lines = []
    g = open('DCM_table.txt','w')
    for line in dcm_attn:
        l = ' '.join(map(str,line))
        lines.append(l)
        g.write(l+'\n')
    g.close()
    ch.dcm_table2sql(lines)
    # Connect to ACC /parm directory and transfer scan_header files
    try:
        g = open('DCM_table.txt','r')
        acc = FTP('acc.solar.pvt')
        acc.login('admin','observer')
        acc.cwd('parm')
        # Send DCM table lines to ACC
        print acc.storlines('STOR dcm.txt',g)
        g.close()
        print 'Successfully wrote dcm.txt to ACC'
    except:
        print 'Cannot FTP dcm.txt to ACC'
    
def DCM_master_attn_cal(fseqfile=None,dcmattn=None,update=False):
    ''' New version of this command, which uses the power values in
        the 10gbe packet headers instead of the very slow measurement
        of the ADC levels themselves.  This version only takes about 8 s!
        
        If update is True, it writes the results to the SQL database.
        
        Returns the DCM_master_table in the form of lines of text
        strings, with labels (handy for viewing).
    '''
    import pcapture2 as p
    import dbutil as db
    import cal_header as ch
    import stateframe as stf
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
            bandlist = np.array(map(int,line.split(',')))

    # Make sure this sequence is actually running, or start it if not
    accini = stf.rd_ACCfile()
    if not fseq_is_running(fseqfile,accini):
        # Sequence is not running, so send ACC commands to start it
        send_cmds(['FSEQ-OFF'],accini)
        send_cmds(['FSEQ-INIT'],accini)
        send_cmds(['FSEQ-FILE '+fseqfile],accini)
        send_cmds(['FSEQ-ON'],accini)
        bandlist2dcmtable(bandlist)
        time.sleep(3)
        if not fseq_is_running(fseqfile,accini):
            print 'Frequency sequence could not be started.'
            return None
        else:
            print 'Successfully started frequency sequence.'
        send_cmds(['dcmtable dcm.txt'],accini)
        send_cmds(['dcmauto-on'],accini)

    pwr = np.zeros((50,8,4),'int')
    # Capture on eth2 interface
    command = 'tcpdump -i eth2 -c 155000 -w /home/user/Python/dcm2.pcap -s 1000'
    p.sendcmd(command)
    # Capture on eth3 interface
    command = 'tcpdump -i eth3 -c 155000 -w /home/user/Python/dcm3.pcap -s 1000'
    p.sendcmd(command)
    headers = p.list_header('/home/user/Python/dcm2.pcap')
    for line in headers:
        try:
            j, id, p1,p2,p3,p4 = np.array(map(int,line.split()))[[0,3,6,7,8,9]]
            pwr[j,id] = (p1, p2, p3, p4)
        except:
            # This is to skip the non-data header lines in the list
            pass
    headers = p.list_header('/home/user/Python/dcm3.pcap')
    for line in headers:
        try:
            j, id, p1,p2,p3,p4 = np.array(map(int,line.split()))[[0,3,6,7,8,9]]
            pwr[j,id] = (p1, p2, p3, p4)
        except:
            # This is to skip the non-data header lines in the list
            pass
    # Reshape to (slot, nant, npol)
    pwr.shape = (50,16,2)
#    # Read current frequency sequence from database
#    cursor = db.get_cursor()
#    query = 'select top 50 FSeqList from hV37_vD50 order by Timestamp desc'
#    fseq, msg = db.do_query(cursor, query)
#    if msg == 'Success':
#        fseqlist = fseq['FSeqList'][::-1]  # Reverse the order
#        bandlist = ((np.array(fseqlist)-0.44)*2).astype(int)
#    cursor.close()
    if dcmattn is None:
        # Read current DCM_master_table from database
        xml, buf = ch.read_cal(2)
        orig_table = stf.extract(buf,xml['Attenuation'])
    else:
        # DCM attenuation is set to a constant value so create a table of such values.
        orig_table = np.zeros((34,30)) + dcmattn
        orig_table[:,26:] = 0
    # Order pwr values according to bandlist, taking median of any repeated values
    new_pwr = np.zeros((34,16,2))
    for i in range(34):
        idx, = np.where(bandlist-1 == i)
        if len(idx) > 0:
            new_pwr[i] = np.median(pwr[idx],0)
    new_pwr.shape = (34,32)
    # Now determine the change in attenuation needed to achieve a target
    # value of 1600.  Eliminate last two entries, corresponding to Ant16
    attn = np.log10(new_pwr[:,:-2]/1600.)*10.
    new_table = (np.clip(orig_table + attn,0,30)/2).astype(int)*2
    DCMlines = []
    DCMlines.append('#      Ant1  Ant2  Ant3  Ant4  Ant5  Ant6  Ant7  Ant8  Ant9 Ant10 Ant11 Ant12 Ant13 Ant14 Ant15')
    DCMlines.append('#      X  Y  X  Y  X  Y  X  Y  X  Y  X  Y  X  Y  X  Y  X  Y  X  Y  X  Y  X  Y  X  Y  X  Y  X  Y')
    DCMlines.append('#     ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- -----')
    for band in range(1,35):
        DCMlines.append('{:2} :  {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2}'.format(band,*new_table[band-1]))
    if update:
        msg = ch.dcm_master_table2sql(DCMlines)
        if msg:
            print 'Success'
        else:
            print 'Error writing table to SQL database!'
    return DCMlines
    
#def DCM_master_attn_cal(roach_list,ant_list='ant1-15'):
    #''' Perform a sequence of FEM settings, using ADC levels to 
        #deduce optimum DCM attenuation settings for all 34 bands
        #and return the master table DCM attenuation lines.
        #TAKES ABOUT 5 MINUTES TO COMPLETE

        #roach_list  a set of ROACH objects created with roach.py
        #ant_list    a list of antennas in the form of a string,
                      #e.g. "ant1-5 ant7" on which to adjust FEMs
                      #Default is all antennas, and an empty string
                      #means all antennas in current subarray.
        
        #Returns DCM_lines list of strings, one for each band
    #'''
    #accini = stf.rd_ACCfile()
    #acc = {'host': accini['host'], 'scdport':accini['scdport']}

    #n = len(roach_list)
    #adc_ndon = np.zeros((34,n,4),dtype='float')
    ## Set DCM state to standard values
    #send_cmds(['DCMAUTO-OFF '+ant_list,'DCMATTN 12 12 '+ant_list],acc)
    ## Set FEM attenuations to nominal
    #send_cmds(['FEMATTN 0 '+ant_list],acc)
    ## Cycle through bands to get "nd-on" ADC levels
    #send_cmds(['ND-ON '+ant_list],acc)
    #for band in range(34):
        #acc_tune(band+1,acc)
        #time.sleep(3)
        #r.adc_levels(roach_list)
        #for i,ro in enumerate(roach_list):
            #adc_ndon[band,i] = ro.adc_levels
    #send_cmds(['ND-OFF '+ant_list],acc)
    #DCM_lines = make_DCM_table(roach_list,adc_ndon,dcm_base=12,adc_nom=30)
    #return DCM_lines
    
def DCM_attn_anal(filename):
    ''' Analyze a DCMATTNTEST observation to determine the 2- and 4-bit
        attenuation values.  Input is a Miriad file.  Returns two arrays, 
           at2 and at4 of size (nant,npol) = (13,2)
        representing the attenuation, in dB, of the 2- and 4-bit, resp.
    '''
    import read_idb as ri
    import dbutil as db
    import cal_header as ch
    import stateframe as stf
    import copy
    from util import Time
    import matplotlib.pylab as plt
    
    out = ri.read_idb([filename])
    ts = int(Time(out['time'][0],format='jd').lv + 0.5)
    te = int(Time(out['time'][-1],format='jd').lv + 0.5)
    query = 'select Timestamp,DCM_Offset_Attn from fV65_vD15 where Timestamp between '+str(ts)+' and '+str(te)+' order by Timestamp'
    cursor = db.get_cursor()
    data, msg = db.do_query(cursor,query)
    cursor.close()
    dcm_offset = data['DCM_Offset_Attn'].reshape(len(data['DCM_Offset_Attn'])/15,15)
    dcm_offset = dcm_offset[:,0]   # All antennas are the same
    t = Time(out['time'][0],format='jd')
    xml, buf = ch.read_cal(2,t)
    table = stf.extract(buf,xml['Attenuation'])
    bandlist = ((out['fghz']-0.5)*2).astype(int)
    tbl = table[bandlist-1]
    tbl.shape = (len(bandlist),15,2)
    tbl = np.swapaxes(np.swapaxes(tbl,0,-1),0,1)
    tbl2 = np.broadcast_to(tbl,(out['time'].shape[0],15,2,134))
    tbl = copy.copy(np.rollaxis(tbl2,0,4))  # Shape (nant,npol,nf,nt)
    pwr = out['p'][:15]  # Shape (nant,npol,nf,nt)
    # Add value of dcm_offset to table
    for i,offset in enumerate(dcm_offset):
        tbl[:,:,:,i] += offset
    # Clip to valid attenuations
    tbl = np.clip(tbl,0,30)
    # Isolate good times in various attn states
    goodm2, = np.where(dcm_offset == -2)
    goodm2 = goodm2[2:-3]
    good2, = np.where(dcm_offset == 2)
    good2 = good2[2:-3]
    good0, = np.where(dcm_offset[goodm2[-1]:good2[0]] == 0)
    good0 += goodm2[-1]
    good0 = good0[2:-3]
    good4, = np.where(dcm_offset == 4)
    good4 = good4[2:-3]
    good6, = np.where(dcm_offset == 6)
    good6 = good6[2:-3]
    goodbg = good6 + 30  # Assumes FEMATTN 15 follows good6 30 s later
    # Perform median over good times and create pwrmed with medians
    # The 5 indexes correspond to dcm_offsets -2, 0, 2, 4 and 6
    nant,npol,nf,nt = pwr.shape
    pwrmed = np.zeros((nant,npol,nf,5))
    # Do not forget to subtract the background
    bg = np.median(pwr[:,:,:,goodbg],3)
    pwrmed[:,:,:,0] = np.median(pwr[:,:,:,goodm2],3) - bg
    pwrmed[:,:,:,1] = np.median(pwr[:,:,:,good0],3) - bg
    pwrmed[:,:,:,2] = np.median(pwr[:,:,:,good2],3) - bg
    pwrmed[:,:,:,3] = np.median(pwr[:,:,:,good4],3) - bg
    pwrmed[:,:,:,4] = np.median(pwr[:,:,:,good6],3) - bg
    good = np.array([goodm2[0],good0[0],good2[0],good4[0],good6[0]])
    tbl = tbl[:,:,:,good]
    at2 = np.zeros((13,2),float)
    at4 = np.zeros((13,2),float)
    at8 = np.zeros((13,2),float)
    f1, ax1 = plt.subplots(2,13)
    f2, ax2 = plt.subplots(2,13)
    f3, ax3 = plt.subplots(2,13)
    for ant in range(13):
        for pol in range(2):
            pts = []
            for i in range(4):
                for v in [0,4,8,12,16,20,24,28]:
                    idx, = np.where(tbl[ant,pol,:,i] == v)
                    if len(idx) != 0:
                        good, = np.where((tbl[ant,pol,idx,i] + 2) == tbl[ant,pol,idx,i+1])
                        if len(good) != 0:
                            pts.append(pwrmed[ant,pol,idx[good],i]/pwrmed[ant,pol,idx[good],i+1])
            pts = np.concatenate(pts)
            ax1[pol,ant].plot(pts,'.')
            ax1[pol,ant].set_ylim(0,2)
            at2[ant,pol] = np.log10(np.median(pts))*10.
            pts = []
            for i in range(3):
                for v in [0,2,8,10,16,18,24,26]:
                    idx, = np.where(tbl[ant,pol,:,i] == v)
                    if len(idx) != 0:
                        good, = np.where((tbl[ant,pol,idx,i] + 4) == tbl[ant,pol,idx,i+2])
                        if len(good) != 0:
                            pts.append(pwrmed[ant,pol,idx[good],i]/pwrmed[ant,pol,idx[good],i+2])
            pts = np.concatenate(pts)
            ax2[pol,ant].plot(pts,'.')
            ax2[pol,ant].set_ylim(0,3)
            at4[ant,pol] = np.log10(np.median(pts))*10.
            pts = []
            i = 0
            for v in [0,2,4,6,16,18,20,22]:
                idx, = np.where(tbl[ant,pol,:,i] == v)
                if len(idx) != 0:
                    good, = np.where((tbl[ant,pol,idx,i] + 8) == tbl[ant,pol,idx,i+4])
                    if len(good) != 0:
                        pts.append(pwrmed[ant,pol,idx[good],i]/pwrmed[ant,pol,idx[good],i+4])
            try:
                pts = np.concatenate(pts)
            except:
                # Looks like there were no points for this antenna/polarization, so set to nominal attn
                pts = [6.30957,6.30957,6.30957]
            ax3[pol,ant].plot(pts,'.')
            ax3[pol,ant].set_ylim(5,8)
            at8[ant,pol] = np.log10(np.median(pts))*10.
    plt.show()
    # Generate output table, a complex array of size (nant,npol,nbits)
    attn = np.zeros((16,2,4),np.complex)
    # Set to nominal values, then overwrite with measured ones
    for i in range(16):
        for j in range(2):
            attn[i,j] = [2.0+0j, 4.0+0j, 8.0+0j, 16.0+0j]
    attn[:13,:,0] = at2 + 0j
    attn[:13,:,1] = at4 + 0j
    attn[:13,:,2] = at8 + 0j
    return attn
    
def plot_adc_cal(roach_list,adc_nosig,adc_ndoff,adc_ndon):
    import matplotlib.pylab as plt
    n = len(roach_list)
    chans = ['X','Y']
    f, ax = plt.subplots(n,4)
    f.set_size_inches(10,2.5*n, forward=True)
    for i in range(n):
        rstr = 'Roach'+roach_list[i].roach_ip[5:6]
        for j in range(4):
            ant = roach_list[i].ants[j / 2]
            chan = chans[j % 2]
            astr = ' Ant '+str(ant)+chan+':'
            ax[i,j].plot(adc_nosig[:,i,j],'.')
            ax[i,j].plot(adc_ndoff[:,i,j],'.')
            ax[i,j].plot(adc_ndon[:,i,j],'.')
            ax[i,j].set_ylim(0, 60)
            ax[i,j].text(5,50,rstr+astr,fontsize=10)
    plt.show()

def make_DCM_table(roach_list,adc_ndon,dcm_base=12,adc_nom=30):

    DCMlines = []
    DCMlines.append('#      Ant1  Ant2  Ant3  Ant4  Ant5  Ant6  Ant7  Ant8  Ant9 Ant10 Ant11 Ant12 Ant13 Ant14 Ant15')
    DCMlines.append('#      X  Y  X  Y  X  Y  X  Y  X  Y  X  Y  X  Y  X  Y  X  Y  X  Y  X  Y  X  Y  X  Y  X  Y  X  Y')
    DCMlines.append('#     ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- -----')
    for band in range(1,35):
        out = np.zeros(32,dtype='int') + 12  # Default to 12 dB if not present
        for i in range(len(roach_list)):
            # Calculate DCM attenuation for the 4 channels on this roach at this band
            # The target standard deviation is adc_nom (default is 30), and the base
            # attenuation at which the observations were made is dcm_base (default is 12 dB).
            # This uses the ratio of standard deviations to determine the factor in dB
            # needed to change it to the target standard deviation.  The division by two,
            # conversion to integer, and multiplication by 2 is because the attenuation steps
            # are in units of 2 dB.  The result is clipped to be between 0 and 30 dB.
            ch_atn = np.clip(((10*np.log(adc_ndon[band-1,i,:]/adc_nom)+dcm_base + 1)/2).astype('int')*2,0,30)
            # Determine the two antennas on this roach (-1 converts to 0-based index)
            ant1,ant2 = roach_list[i].ants - 1
            # Use indexes to assign the 4 channels to the right place in the array
            out[np.array((ant1*2,ant1*2+1,ant2*2,ant2*2+1))] = ch_atn
        DCMlines.append('{:2} :  {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2}'.format(band,*out[:30]))
    return DCMlines
        
def adc_check(roach_list,dcmlines,ant_list='ant1-15',do_plot=False):
    ''' Perform a sequence of FEM settings, using ADC levels to 
        deduce optimum DCM attenuation settings for all 34 bands.
        This can also reveal problems in FEM or DCM hardware.
        TAKES ABOUT 17 MINUTES TO COMPLETE
        
        roach_list  a set of ROACH objects created with roach.py
        ant_list    a list of antennas in the form of a string,
                      e.g. "ant1-5 ant7" on which to adjust FEMs
                      Default is all antennas, and an empty string
                      means all antennas in current subarray.
        do_plot     if True, makes a summary plot of results
        
        Returns numpy arrays :
                adc_nosig[34, nroach, 4] (no-signal ADC levels)
                adc_ndoff[34, nroach, 4] (ADC levels for ND-OFF)
                adc_ndon [34, nroach, 4] (ADC levels for ND-ON)
    '''
    accini = stf.rd_ACCfile()
    acc = {'host': accini['host'], 'scdport':accini['scdport']}

    n = len(roach_list)
    adc_nosig = np.zeros((34,n,4),dtype='float')
    adc_ndoff = np.zeros((34,n,4),dtype='float')
    adc_ndon = np.zeros((34,n,4),dtype='float')
    # Set DCM state to standard values
    send_cmds(['DCMAUTO-OFF '+ant_list,'DCMATTN 12 12 '+ant_list],acc)
    # Set FEM attenuations to maximum
    send_cmds(['FEMATTN 15 '+ant_list],acc)
    # Cycle through bands to get "zero-input" ADC levels
    for band in range(34):
        acc_tune(band+1,acc)
        line = dcmlines[band+3]
        for ant in range(1,16):
            send_cmds(['DCMATTN'+line[ant*6-1:(ant+1)*6-1]+' ant'+str(ant)],acc)
            time.sleep(1)
        r.adc_levels(roach_list)
        for i,ro in enumerate(roach_list):
            adc_nosig[band,i] = ro.adc_levels
    # Set FEM attenuations to nominal
    send_cmds(['FEMATTN 0 '+ant_list],acc)
    # Cycle through bands to get "nd-on" ADC levels
    send_cmds(['ND-ON '+ant_list],acc)
    for band in range(34):
        acc_tune(band+1,acc)
        line = dcmlines[band+3]
        for ant in range(1,16):
            send_cmds(['DCMATTN'+line[ant*6-1:(ant+1)*6-1]+' ant'+str(ant)],acc)
            time.sleep(1)
        r.adc_levels(roach_list)
        for i,ro in enumerate(roach_list):
            adc_ndon[band,i] = ro.adc_levels
    # Cycle through bands to get "nd-off" ADC levels
    send_cmds(['ND-OFF '+ant_list],acc)
    for band in range(34):
        acc_tune(band+1,acc)
        line = dcmlines[band+3]
        for ant in range(1,16):
            send_cmds(['DCMATTN'+line[ant*6-1:(ant+1)*6-1]+' ant'+str(ant)],acc)
            time.sleep(1)
        r.adc_levels(roach_list)
        for i,ro in enumerate(roach_list):
            adc_ndoff[band,i] = ro.adc_levels
    if do_plot: 
        plot_adc_cal(roach_list, adc_nosig, adc_ndoff, adc_ndon)
    return adc_nosig, adc_ndoff, adc_ndon

def gain_state(trange=None):
    ''' Read and assemble the gain state for the given timerange from 
        the SQL database, or for the last 10 minutes if trange is None.
        
        Returns the complex attenuation of the FEM for the timerange
        as an array of size (nant, npol, ntimes) [not band dependent],
        and the complex attenuation of the DCM for the same timerange
        as an array of size (nant, npol, nbands, ntimes).  Also returns
        the time as a Time() object array.
    '''
    from util import Time
    import dbutil as db
    from fem_attn_calib import fem_attn_update
    import cal_header as ch
    
    if trange is None:
        t = Time.now()
        t2 = Time(t.jd - 600./86400.,format='jd')
        trange = Time([t2.iso,t.iso])
    ts = trange[0].lv   # Start timestamp
    te = trange[1].lv   # End timestamp
    cursor = db.get_cursor()
    # First get FEM attenuation for timerange
    D15dict=db.get_dbrecs(cursor,dimension=15,timestamp=trange)
    DCMoffdict=db.get_dbrecs(cursor,dimension=50,timestamp=trange)
    DCMoff_v_slot = DCMoffdict['DCMoffset_attn']
#    DCMoff_0 = D15dict['DCM_Offset_Attn'][:,0]  # All ants are the same
    fem_attn={}
    fem_attn['timestamp']=D15dict['Timestamp'][:,0]
    nt=len(fem_attn['timestamp'])
    junk=np.zeros([nt,1],dtype='int') #add the non-existing antenna 16
    fem_attn['h1']=np.append(D15dict['Ante_Fron_FEM_HPol_Atte_First'],junk,axis=1) #FEM hpol first attn value
    fem_attn['h2']=np.append(D15dict['Ante_Fron_FEM_HPol_Atte_Second'],junk,axis=1) #FEM hpol second attn value
    fem_attn['v1']=np.append(D15dict['Ante_Fron_FEM_VPol_Atte_First'],junk,axis=1) #FEM vpol first attn value
    fem_attn['v2']=np.append(D15dict['Ante_Fron_FEM_VPol_Atte_Second'],junk,axis=1) #FEM vpol second attn value
    fem_attn['ants']=np.append(D15dict['I15'][0,:],[15])
    # Add corrections from SQL database for start time of timerange
    fem_attn_corr = fem_attn_update(fem_attn,trange[0])
    # Next get DCM attenuation for timerange
    # Getting next earlier scan header
    ver = db.find_table_version(cursor,ts, True)
    query = 'select top 50 Timestamp,FSeqList from hV'+ver+'_vD50 where Timestamp <= '+str(ts)+' order by Timestamp desc'
    fseq, msg = db.do_query(cursor, query)
    if msg == 'Success':
        fseqlist = fseq['FSeqList'][::-1]  # Reverse the order
        bandlist = ((np.array(fseqlist)-0.44)*2).astype(int)
    cursor.close()
    # Read current DCM_table from database
    xml, buf = ch.read_cal(3,trange[0])
    orig_table = stf.extract(buf,xml['Attenuation']).astype('int')
    orig_table.shape = (50,15,2)
    xml, buf = ch.read_cal(6,trange[0])
    dcm_attn_bitv=np.nan_to_num(stf.extract(buf, xml['DCM_Attn_Real'])) + np.nan_to_num(stf.extract(buf, xml['DCM_Attn_Imag'])) * 1j
#    # Add one more bit (all zeros) to take care of unit bit
#    dcm_attn_bitv = np.concatenate((np.zeros((16,2,1),'int'),dcm_attn_bitv),axis=2)
    # We now have:
    #   orig_table     the original DCM at start of scan, size (nslot, nant=15, npol)
    #   DCMoff_0       the offset applied to all antennas and slots (ntimes)
    #   DCMoff_v_slot  the offest applied to all antennas but varies by slot (ntimes, nslot)
    #   dcm_attn_bitv  the measured (non-nominal) attenuations for each bit value (nant=16, npol, nbit) -- complex
    # Now I need to convert slot to band, add appropriately, and organize as (nant=16, npol, nband, ntimes)
    # Add one more antenna (all zeros) to orig_table
    orig_table = np.concatenate((orig_table,np.zeros((50,1,2),'int')),axis=1)
    ntimes, nslot = DCMoff_v_slot.shape
    dcm_attn = np.zeros((16,2,34,ntimes),np.int)
    for i in range(ntimes):
        for j in range(50):
            idx = bandlist[j]-1
            # This adds attenuation for repeated bands--hopefully the same value for each repeat
            dcm_attn[:,:,idx,i] += orig_table[j,:,:]+DCMoff_v_slot[i,j]
    # Normalize repeated bands by finding number of repeats and dividing.
    for i in range(1,35):
        n = len(np.where(bandlist == i)[0])
        if n > 1:
            dcm_attn[:,:,i-1,:] /= n
    # Make sure attenuation is in range
    dcm_attn = np.clip(dcm_attn,0,30)
    # Finally, correct for non-nominal (measured) bit values
    # Start with 0 attenuation as reference
    dcm_attn_corr = dcm_attn*(0+0j)
    att = np.zeros((16,2,34,ntimes,5),np.complex)
    # Calculate resulting attenuation based on bit attn values (2,4,8,16)
    for i in range(4):
        # Need dcm_attn_bitv[...,i] to be same shape as dcm_attn
        bigger_bitv = np.broadcast_to(dcm_attn_bitv[...,i],(ntimes,34,16,2))
        bigger_bitv = np.swapaxes(np.swapaxes(np.swapaxes(bigger_bitv,0,3),1,2),0,1)
        att[...,i] = (np.bitwise_and(dcm_attn,2**(i+1))>>(i+1))*bigger_bitv
        dcm_attn_corr = dcm_attn_corr + att[...,i]

    # Move ntimes column to next to last position, and then sum over last column (the two attenuators)
    fem_attn_corr = np.sum(np.rollaxis(fem_attn_corr,0,3),3)
    # Output is FEM shape (nant, npol, ntimes) = (16, 2, ntimes)
    #           DCM shape (nant, npol, nband, ntimes) = (16, 2, 34, ntimes)
    # Arrays are complex, in dB units
    tjd = Time(fem_attn['timestamp'].astype('int'),format='lv').jd
    return fem_attn_corr, dcm_attn_corr, tjd
    
