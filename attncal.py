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
#  2020-05-10  DG
#    Updated get_attncal() to use util.get_idbdir() to find IDB root path.
#  2020-05-11  DG
#    Further update to make this work on the DPP.
#  2021-07-13  DG
#    I noticed that my guessing the indexes of attenuation transitions was
#    sometimes failing, so now it is "done right" by using the SQL record
#    of attenuation settings to determine the different attenuation states.
#  2021-07-27  DG
#    The new code was still failing due to bad SQL records (erroneous 0 levels)
#    so I now eliminate those before searching for transitions.  I also
#    commented out the attna calculation, which is not used anywhere and
#    cannot possibly be useful (due to uncorrected saturation).
#
from util import Time
import numpy as np
import dump_tsys as dt
import read_idb as ri

def get_attncal(trange, do_plot=False, dataonly=False):
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
        
        The dataonly parameter tells the routine to skip calculating the attenuation
        and only return the IDB data from the (first) gaincal.
    '''
    from util import get_idbdir, fname2mjd, nearest_val_idx
    import socket
    import dbutil
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

            datadir = get_idbdir(Time(mjd,format='mjd'))
            # Add date path if on pipeline
            # if datadir.find('eovsa') != -1: datadir += fdb['FILE'][gcidx][0][3:11]+'/'

            host = socket.gethostname()
            if host == 'pipeline': datadir += fdb['FILE'][gcidx][0][3:11]+'/'

            file = datadir + fdb['FILE'][gcidx][0]
            out = ri.read_idb([file])
            if dataonly:
                return out
            # Get time from filename and read 120 records of attn state from SQL database
            filemjd = fname2mjd(fdb['FILE'][gcidx][0])
            cursor = dbutil.get_cursor()
            d15 = dbutil.get_dbrecs(cursor, dimension=15, timestamp=Time(filemjd,format='mjd'), nrecs=120)
            cursor.close()
            # Find time indexes of the 62 dB attn state
            # Uses only ant 1 assuming all are the same
            dtot = (d15['Ante_Fron_FEM_HPol_Atte_Second'] + d15['Ante_Fron_FEM_HPol_Atte_First'])[:,0]
            good, = np.where(dtot != 0)
            #import pdb; pdb.set_trace()
            # Indexes into SQL records where a transition occurred.
            transitions, = np.where(dtot[good] - np.roll(dtot[good],1) != 0)
            # Eliminate any zero-index transition (if it exists)
            if transitions[0] == 0:
                transitions = transitions[1:]
            # These now have to be translated into indexes into the data, using the times
            idx = nearest_val_idx(d15['Timestamp'][good,0][transitions],Time(out['time'],format='jd').lv)
            vx = np.nanmedian(out['p'][:13,:,:,np.arange(idx[0]+1,idx[1]-1)],3)
            va = np.mean(out['a'][:13,:2,:,np.arange(idx[0]+1,idx[1]-1)],3)
            vals = []
            attn = []
            for i in range(1,10):
                vals.append(np.nanmedian(out['p'][:13,:,:,np.arange(idx[i]+1,idx[i+1]-1)],3) - vx)
                attn.append(np.log10(vals[0]/vals[-1])*10.)
            #vals = []
            #attna = []
            #for i in range(1,10):
            #    vals.append(np.median(out['a'][:13,:2,:,np.arange(idx[i],idx[i+1])],3) - va)
            #    attna.append(np.log10(vals[0]/vals[-1])*10.)
            
            if do_plot:
                for i in range(13):
                    for j in range(2):
                        ax[j,i].plot(out['fghz'],attn[1][i,j],'.',markersize=3)
                        #ax[j,i].plot(out['fghz'],attna[1][i,j],'.',markersize=1)
                        ax[j+2,i].plot(out['fghz'],attn[2][i,j],'.',markersize=3)
                        #ax[j+2,i].plot(out['fghz'],attna[2][i,j],'.',markersize=1)
            outdict.append({'time': Time(out['time'][0],format='jd'),'fghz': out['fghz'], 
                            'rcvr_auto':va, # 'attna': np.array(attna[1:]), 
                            'rcvr':vx, 'attn': np.array(attn[1:])})
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
        
