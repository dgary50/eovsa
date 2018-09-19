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
    calfac = pc.get_calfac(trange[0])
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
    # Subtract background for each antenna/polarization
    for ant in antlist:
        for pol in range(2):
            bg = np.median(out['p'][ant,pol,:,brange[0]:brange[1]],1).repeat(nt).reshape(nf,nt)
            out['p'][ant,pol] -= bg
    # Form median over antennas/pols
    med = np.mean(np.median(out['p'][antlist],0),0)
    # Do background subtraction once more for good measure
    bg = np.median(med[:,0:300],1).repeat(nt).reshape(nf,nt)
    med -= bg
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
    return out, med
