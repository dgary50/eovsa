import pcapture2 as p
import dbutil as db
import cal_header as ch
import stateframe as stf
import numpy as np

def DCM_master_attn_cal(update=False):
    ''' New version of this command, which uses the power values in
        the 10gbe packet headers instead of the very slow measurement
        of the ADC levels themselves.  This version only takes about 8 s!
        
        If update is True, it writes the results to the SQL database.
        
        Returns the DCM_master_table in the form of lines of text
        strings, with labels (handy for viewing).
    '''
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
    # Read current frequency sequence from database
    cursor = db.get_cursor()
    query = 'select top 50 FSeqList from hV37_vD50 order by Timestamp desc'
    fseq, msg = db.do_query(cursor, query)
    if msg == 'Success':
        fseqlist = fseq['FSeqList'][::-1]  # Reverse the order
        bandlist = ((np.array(fseqlist)-0.44)*2).astype(int)
    cursor.close()
    # Read current DCM_master_table from database
    xml, buf = ch.read_cal(2)
    orig_table = stf.extract(buf,xml['Attenuation'])
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
    
if __name__ == "__main__":

    import sys
    print len(sys.argv)
    if len(sys.argv) == 2:
        if sys.argv[1] == 'update':
            lines = DCM_master_attn_cal(True)
            for line in lines:
                print line
    else:
            lines = DCM_master_attn_cal()
            for line in lines:
                print line
