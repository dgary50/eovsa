#    ATTNCAL
# History:
#  2017-09-07  DG
#    First written
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
                           
        N.B.: Ignores days with other than one GAINCALTEST measurement, e.g. 0 or 2,
              the first is obvious, while the second is because there is no way to
              tell which of the 2 are good.
    '''
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
                ax[j,i].set_ylim(-3,-1)
                ax[j+2,i].set_ylim(-5,-3)
    outdict = []
    for mjd in range(mjd1,mjd2+1):
        fdb = dt.rd_fdb(Time(mjd,format='mjd'))
        gcidx, = np.where(fdb['PROJECTID'] == 'GAINCALTEST')
        if len(gcidx) == 1:
            print fdb['FILE'][gcidx]
            file =  '/data1/eovsa/fits/IDB/'+fdb['FILE'][gcidx][0][3:11]+'/'+fdb['FILE'][gcidx][0]
            out = ri.read_idb([file])
            vx = np.mean(out['p'][:13,:,:,6:12],3)
            val0 = np.median(out['p'][:13,:,:,16:22],3) - vx
            val1 = np.median(out['p'][:13,:,:,26:32],3) - vx
            val2 = np.median(out['p'][:13,:,:,36:42],3) - vx
            val3 = np.median(out['p'][:13,:,:,46:52],3) - vx
            val4 = np.median(out['p'][:13,:,:,56:62],3) - vx
            val5 = np.median(out['p'][:13,:,:,66:72],3) - vx
            val6 = np.median(out['p'][:13,:,:,76:82],3) - vx
            val7 = np.median(out['p'][:13,:,:,86:92],3) - vx
            val8 = np.median(out['p'][:13,:,:,96:102],3) - vx
            attn1 = np.log10(val1/val0)*10.
            attn2 = np.log10(val2/val0)*10.
            attn3 = np.log10(val3/val0)*10.
            attn4 = np.log10(val4/val0)*10.
            attn5 = np.log10(val5/val0)*10.
            attn6 = np.log10(val6/val0)*10.
            attn7 = np.log10(val7/val0)*10.
            attn8 = np.log10(val8/val0)*10.
            if do_plot:
                for i in range(13):
                    for j in range(2):
                        ax[j,i].plot(out['fghz'],attn1[i,j],'.')
                        ax[j+2,i].plot(out['fghz'],attn2[i,j],'.')
            outdict.append({'time': Time(out['time'][0],format='jd'),'fghz': out['fghz'], 
                            'attn': np.array([attn1, attn2, attn3, attn4, attn5, attn6, attn7, attn8])})
    return outdict
