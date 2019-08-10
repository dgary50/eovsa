import pipeline_cal as pc
import gaincal2 as gc
import attncal as ac
import numpy as np
from util import Time, nearest_val_idx, ant_str2list
import matplotlib.pylab as plt
from matplotlib.dates import DateFormatter

def autocorrect(out,ant_str='ant1-13',brange=[0,300]):
    nt = len(out['time'])
    nf = len(out['fghz'])
    pfac1 = (out['p'][:,:,:,:-1] - out['p'][:,:,:,1:])/out['p'][:,:,:,:-1]
    trange = Time(out['time'][[0,-1]],format='jd')
    src_lev = gc.get_fem_level(trange)   # Read FEM levels from SQL
    # Match times with data
    tidx = nearest_val_idx(out['time'],src_lev['times'].jd)
    # Find attenuation changes
    for ant in range(13):
        for pol in range(2):
            if pol == 0:
                lev = src_lev['hlev'][ant,tidx]
            else:
                lev = src_lev['vlev'][ant,tidx]
            jidx, = np.where(abs(lev[:-1] - lev[1:]) == 1)
            for freq in range(nf):
                idx, = np.where(np.logical_and(abs(pfac1[ant,pol,freq]) > 0.05,abs(pfac1[ant,pol,freq]) < 0.95))
                for i in range(len(idx-1)):
                    if idx[i] in jidx or idx[i] in jidx-1:
                        out['p'][ant,pol,freq,idx[i]+1:] /= (1-pfac1[ant,pol,freq,idx[i]])
    # Time of total power calibration is 20 UT on the date given
    tptime = Time(np.floor(trange[0].mjd) + 20./24.,format='mjd')
    calfac = pc.get_calfac(tptime)
    tpcalfac = calfac['tpcalfac']
    tpoffsun = calfac['tpoffsun']
    hlev = src_lev['hlev'][:13,0]
    vlev = src_lev['vlev'][:13,0]
    attn_dict = ac.read_attncal(trange[0])[0]   # Read GCAL attn from SQL
    attn = np.zeros((13,2,nf))
    for i in range(13):
        attn[i,0] = attn_dict['attn'][hlev[i],0,0]
        attn[i,1] = attn_dict['attn'][vlev[i],0,1]
        print 'Ant',i+1,attn[i,0,20],attn[i,1,20]
    attnfac = 10**(attn/10.)
    for i in range(13): print attnfac[i,0,20],attnfac[i,1,20]
    for i in range(nt):
        out['p'][:13,:,:,i] = (out['p'][:13,:,:,i]*attnfac - tpoffsun)*tpcalfac 
    antlist = ant_str2list(ant_str)
    bg = np.zeros_like(out['p'])
    # Subtract background for each antenna/polarization
    for ant in antlist:
        for pol in range(2):
            bg[ant,pol] = np.median(out['p'][ant,pol,:,brange[0]:brange[1]],1).repeat(nt).reshape(nf,nt)
            #out['p'][ant,pol] -= bg
    # Form median over antennas/pols
    med = np.mean(np.median((out['p']-bg)[antlist],0),0)
    # Do background subtraction once more for good measure
    bgd = np.median(med[:,brange[0]:brange[1]],1).repeat(nt).reshape(nf,nt)
    med -= bgd
    pdata = np.log10(med)
    f, ax = plt.subplots(1,1)
    vmax = np.median(np.nanmax(pdata,1))
    im = ax.pcolormesh(Time(out['time'],format='jd').plot_date,out['fghz'],pdata,vmin=1,vmax=vmax)
    ax.axvspan(Time(out['time'][brange[0]],format='jd').plot_date,Time(out['time'][brange[1]],format='jd').plot_date,color='w',alpha=0.3)
    plt.colorbar(im,ax=ax,label='Log Flux Density [sfu]')
    ax.xaxis_date()
    ax.xaxis.set_major_formatter(DateFormatter("%H:%M"))
    ax.set_ylim(out['fghz'][0], out['fghz'][-1])
    ax.set_xlabel('Time [UT]')
    ax.set_ylabel('Frequency [GHz]')
    return {'caldata':out, 'med_sub':med, 'bgd':bg}
    
# def tp_bgnd(tpdata):
    # ''' Create time-variable background from ROACH inlet temperature
        # This may be a crude correction, but it has been seen to work.
        
        # Inputs:
          # tpdata   dictionary returned by read_idb()  NB: tpdata is not changed.
          
        # Returns:
          # bgnd     The background fluctuation array of size (nf,nt) to be 
                     # subtracted from any antenna's total power (or mean of
                     # antenna total powers)
    # '''
    # import dbutil as db
    # import read_idb as ri
    # from util import Time
    # outfghz = tpdata['fghz']
    # try:
        # outtime = tpdata['time']
        # trange = Time(outtime[[0,-1]],format='jd')
    # except:
        # outtime = tpdata['ut_mjd']
        # trange = Time(outtime[[0,-1]],format='mjd')
    
    # tstr = trange.lv.astype(int).astype(str)
    # nt = len(outtime)
    # nf = len(outfghz)
    # outpd = Time(outtime,format='jd').plot_date
    # cursor = db.get_cursor()
    # version = db.find_table_version(cursor,int(tstr[0]))
    # query = 'select * from fV'+version+'_vD8 where (Timestamp between '+tstr[0]+' and '+tstr[1]+')'
    # data, msg = db.do_query(cursor, query)
    # pd = Time(data['Timestamp'][::8].astype(int),format='lv').plot_date
    # inlet = data['Sche_Data_Roac_TempInlet'].reshape(len(pd),8)   # Inlet temperature variation
    # sinlet = np.sum(inlet.astype(float),1)
    # sint = np.interp(outpd,pd,sinlet)
    # sint = np.roll(sint,-90)  # Shift phase of variation by 90 s earlier
    # sint -= np.mean(sint)     # Remove offset, to provide zero-mean fluctuation
    # bgnd = np.zeros((nf, nt), float)
    # for i in range(nf):
        # if 13.5 < outfghz[i] < 14.0:
            # bgnd[i] = sint*2
        # elif 14.0 < outfghz[i] < 15.0:
            # bgnd[i] = sint*5        
        # elif outfghz[i] > 15.0:
            # bgnd[i] = sint*3
    # return bgnd

def smooth(x,window_len=11,window='hanning'):
    """smooth the data using a window with requested size.

    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.

    input:
        x: the input signal
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal

    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)

    see also:

    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter

    TODO: the window parameter could be the window itself if an array instead of a string
    NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.
    """

    if x.ndim != 1:
        raise ValueError, "smooth only accepts 1 dimension arrays."

    if x.size < window_len:
        raise ValueError, "Input vector needs to be bigger than window size."


    if window_len<3:
        return x


    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"


    s = np.r_[x[window_len-1:0:-1], x, x[-2:-window_len-1:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w = np.ones(window_len,'d')
    else:
        w = eval('np.'+window+'(window_len)')

    y = np.convolve(w/w.sum(), s, mode='valid')
    return y


def tp_bgnd(tpdata):
    ''' Create time-variable background from ROACH inlet temperature
        This version is far superior to the earlier, crude version, but
        beware that it works best for a long timerange of data, especially
        when there is a flare in the data.
        
        Inputs:
          tpdata   dictionary returned by read_idb()  NB: tpdata is not changed.
          
        Returns:
          bgnd     The background fluctuation array of size (nf,nt) to be 
                     subtracted from any antenna's total power (or mean of
                     antenna total powers)
    '''
    import dbutil as db
    from util import Time, nearest_val_idx
    outfghz = tpdata['fghz']
    try:
        outtime = tpdata['time']
        trange = Time(outtime[[0,-1]],format='jd')
    except:
        outtime = tpdata['ut_mjd']
        trange = Time(outtime[[0,-1]],format='mjd')
    
    tstr = trange.lv.astype(int).astype(str)
    nt = len(outtime)
    if nt < 2000:
        print 'TP_BGND: Error, timebase too small.  Must have at least 2000 time samples.'
        return None
    nf = len(outfghz)
    outpd = Time(outtime,format='jd').plot_date
    cursor = db.get_cursor()
    version = db.find_table_version(cursor,int(tstr[0]))
    query = 'select * from fV'+version+'_vD8 where (Timestamp between '+tstr[0]+' and '+tstr[1]+')'
    data, msg = db.do_query(cursor, query)
    pd = Time(data['Timestamp'][::8].astype(int),format='lv').plot_date
    inlet = data['Sche_Data_Roac_TempInlet'].reshape(len(pd),8)   # Inlet temperature variation
    sinlet = np.sum(inlet.astype(float),1)
    # Eliminate 0 values in sinlet by replacing with nearest good value
    bad, = np.where(sinlet == 0)
    good, = np.where(sinlet != 0)
    idx = nearest_val_idx(bad,good)  # Find locations of nearest good values to bad ones
    sinlet[bad] = sinlet[good[idx]]  # Overwrite bad values with good ones
    # Interpolate sinlet values to the times in the data
    sint = np.interp(outpd,pd,sinlet)
    sint = np.roll(sint,-90)  # Shift phase of variation by 90 s earlier
    sint -= np.mean(sint)     # Remove offset, to provide zero-mean fluctuation
    sdev = np.std(sint)
    sint_ok = np.abs(sint) < sdev
    bgnd = np.zeros((nf, nt), float)
    for i in range(nf):
        # Subtract smooth trend from data
        sig = tpdata['p'][i] - smooth(tpdata['p'][i],2000,'blackman')[1000:-999]
        # Eliminate the worst outliers and repeat
        stdev = np.std(sig)
        good, = np.where(np.abs(sig) < stdev)
        sig = tpdata['p'][i,good] - smooth(tpdata['p'][i,good],2000,'blackman')[1000:-999]
        sint_i = sint[good]
        stdev = np.std(sig)
        # Final check for data quality
        good, = np.where(np.logical_and(sig < 2*stdev, sint_ok[good]))
        p = np.polyfit(sint_i[good],sig[good],1)
        # Apply correction for this frequency
        bgnd[i] = sint*p[0] + p[1]
    return bgnd