#    ATTNCAL
# History:
#  2017-09-07  DG
#    First written
#  2017-09-09  DG
#    Changed the sign of the output, so that it really is attenuation,
#    for consistency with other gaincal2.py routines.
#  2017-10-13  DG
#    Added read_attncal() routine, to read attenuation values from SQL.
#    Also added a __main__ entry point for calling from the command line.
#  2017-10-23  DG
#    Added rcvr key to get_attncal() output, representing the power
#    when 62 dB is inserted in the front end (i.e. is the receiver
#    contribution to the noise).
#
from util import Time
import numpy as np
import dump_tsys as dt
import read_idb as ri

def get_attncal(trange, do_plot=False):
    ''' Finds GAINCALTEST scans from FDB files corresponding to the days
        present in trange Time() object (can be multiple days), calculates
        the attenuation differences for the various FEMATTN states 1-8 
        relative to FEMATTN state 0, and optionally plots the results for
        states 1 and 2 (the most commonly used).  To analyze only a single
        day, trange Time() object can have the same time repeated, or can
        be a single time.
        
        Returns a list of dictionaries, each pertaining to one of the days
        in trange, with keys defined as follows:
           'time':      The start time of the GAINCALTEST scan, as a Time() object
           'fghz':      The list of frequencies [GHz] at which attenuations are measured
           'attn':      The array of attenuations [dB] of size (nattn, nant, npol, nf), 
                           where nattn = 8, nant = 13, npol = 2, and nf is variable
           'rcvr':      The array of receiver noise level (raw units) of size 
                           (nant, npol, nf), where nant = 13, npol = 2, and nf is variable
           'rcvr_auto': Same as rcvr, but for auto-correlation (hence it is complex)
                           
        N.B.: Ignores days with other than one GAINCALTEST measurement, e.g. 0 or 2,
              the first is obvious, while the second is because there is no way to
              tell which of the 2 are good.
    '''
    import os
    if type(trange.mjd) == np.float:
        # Interpret single time as both start and end time
        mjd1 = int(trange.mjd)
        mjd2 = mjd1
    else: 
        mjd1, mjd2 = trange.mjd.astype(int)
    if do_plot:
        import matplotlib.pylab as plt
        f, ax = plt.subplots(4,13)
        f.set_size_inches((14,5))
        ax[0,0].set_ylabel('Atn1X [dB]')
        ax[1,0].set_ylabel('Atn1Y [dB]')
        ax[2,0].set_ylabel('Atn2X [dB]')
        ax[3,0].set_ylabel('Atn2Y [dB]')
        for i in range(13):
            ax[0,i].set_title('Ant '+str(i+1))
            ax[3,i].set_xlabel('Freq [GHz]')
            for j in range(2):
                ax[j,i].set_ylim(1,3)
                ax[j+2,i].set_ylim(3,5)
    outdict = []
    for mjd in range(mjd1,mjd2+1):
        fdb = dt.rd_fdb(Time(mjd,format='mjd'))
        gcidx, = np.where(fdb['PROJECTID'] == 'GAINCALTEST')
        if len(gcidx) == 1:
            print fdb['FILE'][gcidx]

            datadir=os.getenv('EOVSADB')
            if not datadir:
                # go to default directory on pipeline
                datadir='/data1/eovsa/fits/IDB/'+fdb['FILE'][gcidx][0][3:11]+'/'

            file = datadir+fdb['FILE'][gcidx][0]
            out = ri.read_idb([file])
            vx = np.mean(out['p'][:13,:,:,6:12],3)
            va = np.mean(out['a'][:13,:2,:,6:12],3)
            val0 = np.median(out['p'][:13,:,:,16:22],3) - vx
            val1 = np.median(out['p'][:13,:,:,26:32],3) - vx
            val2 = np.median(out['p'][:13,:,:,36:42],3) - vx
            val3 = np.median(out['p'][:13,:,:,46:52],3) - vx
            val4 = np.median(out['p'][:13,:,:,56:62],3) - vx
            val5 = np.median(out['p'][:13,:,:,66:72],3) - vx
            val6 = np.median(out['p'][:13,:,:,76:82],3) - vx
            val7 = np.median(out['p'][:13,:,:,86:92],3) - vx
            val8 = np.median(out['p'][:13,:,:,96:102],3) - vx
            attn1 = np.log10(val0/val1)*10.
            attn2 = np.log10(val0/val2)*10.
            attn3 = np.log10(val0/val3)*10.
            attn4 = np.log10(val0/val4)*10.
            attn5 = np.log10(val0/val5)*10.
            attn6 = np.log10(val0/val6)*10.
            attn7 = np.log10(val0/val7)*10.
            attn8 = np.log10(val0/val8)*10.
            if do_plot:
                for i in range(13):
                    for j in range(2):
                        ax[j,i].plot(out['fghz'],attn1[i,j],'.')
                        ax[j+2,i].plot(out['fghz'],attn2[i,j],'.')
            outdict.append({'time': Time(out['time'][0],format='jd'),'fghz': out['fghz'], 'rcvr':vx, 'rcvr_auto':va,
                            'attn': np.array([attn1, attn2, attn3, attn4, attn5, attn6, attn7, attn8])})
    return outdict
    
def read_attncal(trange=None):
    ''' Given a timerange as a Time() object, read FEM attenuation
        records for each date from the SQL database, and return then
        as a list of attn dictionaries.  To get values for only a single
        day, the trange Time() object can have the same time repeated, or can
        be a single time.
        
        Returns a list of dictionaries, each pertaining to one of the days
        in trange, with keys defined as follows:
           'time':      The start time of the GAINCALTEST scan, as a Time() object
           'fghz':      The list of frequencies [GHz] at which attenuations are measured
           'attn':      The array of attenuations [dB] of size (nattn, nant, npol, nf), 
                           where nattn = 8, nant = 13, npol = 2, and nf is variable

    '''
    import cal_header as ch
    import stateframe as stf
    if trange is None:
        trange = Time.now()
    if type(trange.mjd) == np.float:
        # Interpret single time as both start and end time
        mjd1 = int(trange.mjd)
        mjd2 = mjd1
    else: 
        mjd1, mjd2 = trange.mjd.astype(int)
    attn = []
    for mjd in range(mjd1,mjd2+1):
        # Read next earlier SQL entry from end of given UT day (mjd+0.999) 
        xml, buf = ch.read_cal(7,t=Time(mjd+0.999,format='mjd'))
        t = Time(stf.extract(buf,xml['Timestamp']),format='lv')
        fghz = stf.extract(buf,xml['FGHz'])
        nf = len(np.where(fghz != 0.0)[0])
        fghz = fghz[:nf]
        attnvals = stf.extract(buf,xml['FEM_Attn_Real'])[:,:,:,:nf]
        attn.append({'time':t,'fghz':fghz,'attn':attnvals})
    return attn
    
if __name__ == "__main__":
    ''' Entry for cron job.  Command line can be a single time string, or two time strings
        representing a time range.  GAINCALTEST scans for each day in the time range are
        analyzed and the results are written to the SQL database.
    '''
    import sys
    import cal_header as ch
    t = Time.now()
    if len(sys.argv) == 2:
        try:
            t = Time(sys.argv[1])
        except:
            print 'Cannot interpret',sys.argv[1],'as a valid date/time string.'
            exit()
    if len(sys.argv) == 3:
        try:
            t = Time([sys.argv[1],sys.argv[2]])
        except:
            print 'Cannot interpret',sys.argv[1],'and/or',sys.argv[2],'as a valid date/time string.'
            exit()    
    print t.iso
    attn_list = get_attncal(t)
    for attn in attn_list:
        ch.fem_attn_val2sql([attn])
        
