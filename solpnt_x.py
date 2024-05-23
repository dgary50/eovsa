'''
   Module for accessing, analyzing and plotting SOLPNTCAL X-pattern data
   valid since 2023 Oct 20. In addition to using a different approach to
   observations (an X sweep pattern instead of a + pattern), this code
   does not separately analyze auto-correlation data.  For compatability,
   it still writes out auto-correlation numbers, but these are just copied
   from the total power.'''
#
# History:
#   2023-Oct-24  DG
#     First working version.
#   2023-Nov-03  DG
#     I think I finally got all the signs right...
#   2023-Nov-09  DG
#     Changes needed to deal with special case of Ant 12.  It is an
#     AzEl dish, but its offsets have to be in RA-Dec.  The only change
#     is to calculate a third set of offsets and adjust the code for
#     making the mask.  All else is unchanged (I hope).
#   2023-Nov-13  DG
#     Rather large change to add back autocorrelation calibration, which
#     unfortunately nearly doubles the running time.  I thought it was not
#     really needed, but it is critical for the correct running of udb_corr().
#   2023-Dec-02  DG
#     Fixed a problem where the apply_fem_level() was not applied (with the side
#     effect that no attncal record was ever written)!
#   2024-Jan-15  DG
#     Cosmetic change to plot offset points with transparency (alpha=0.2).
#   2024-Apr-29  DG
#     Add attncal_from_solpnt_to_sql() routine to analyze and write to SQL
#     the "new" GAINCAL information, which is now appended to a SOLPNTCAL scan.
#   2024-May-22  DG
#     Fixed a number of display glitches (some really stupid mistakes...).  Also
#     fixed problem in attncal_from_solpnt_to_sql() that resulted in a new
#     attncal record being written even when one existed already. Wow!  I just
#     discovered that putting an Az offset into a 2-m azel antenna is really
#     a cross-El offset, so I have been calculating it all wrong.  I have
#     to update my trajectories (and the code here) to account for this.
#
import dbutil
import stateframe
import cal_header as ch
from disk_conv import disk_conv
import urllib2
from util import Time, common_val_idx, ant_str2list
from eovsa_array import eovsa_array
from eovsa_tracktable import make_trajtables
import read_idb
from coord_conv import radec2azel
import gaincal2 as gc
import numpy as np
from time import time
import warnings

# This block reads the SQL metadata for the time of interest
nsec = 420  # Number of seconds required for completion of a SOLPNTCAL
nant = 13
antidx = np.arange(nant)                    # List of all 13 antennas
oldidx = np.array([8, 9, 10, 12])           # Indexes of the equatorial antennas
newidx = np.array([0, 1, 2, 3, 4, 5, 6, 7]) # Indexes of the azel antennas
ant12idx = np.array([11])                   # Ant 12 is a special case

def solpnt_meta(tsolpnt):
    ''' Reads the SOLPNT metadata for 7-minutes past the given time, which is the
        start time of the SOLPNTCAL scan.
        
        Input:
           tsolpnt    Time object representing the start time of the SOLPNTCAL scan
        
        Returns:
           meta       A dictionary encapsulating the metadata.
    '''
    t1 = time()
    cursor = dbutil.get_cursor()
    solpntdict = dbutil.get_dbrecs(cursor,dimension=15,timestamp=tsolpnt,nrecs=nsec)
    antlist = np.array([0,1,2,3,4,5,6,7,8,9,10,11,12])
    bgmask = solpntdict['Ante_Fron_FEM_VPol_Atte_First']==31
    ra = (solpntdict['Ante_Cont_RAVirtualAxis'][:,antlist]*np.pi/10000./180.).astype('float')
    dec = (solpntdict['Ante_Cont_DecVirtualAxis'][:,antlist]*np.pi/10000./180.).astype('float')
    rao = (solpntdict['Ante_Cont_RAOffset'][:,antlist]).astype('float')
    deco = (solpntdict['Ante_Cont_DecOffset'][:,antlist]).astype('float')
    azo = (solpntdict['Ante_Cont_AzOffset'][:,antlist]).astype('float')
    elo = (solpntdict['Ante_Cont_ElOffset'][:,antlist]).astype('float')
    times = solpntdict['Timestamp'][:,0].astype('int64').astype('float')
    outdict = stateframe.azel_from_sqldict(solpntdict)
    trk = outdict['TrackFlag'][:,antlist]
    meta = {'Timestamp':times[0],'tstamps':times,'antlist':antlist,'ra':ra,'dec':dec,'rao':rao,
           'deco':deco, 'azo':azo, 'elo':elo, 'trk':trk, 'bgmask':bgmask}
    print 'solpnt_meta took',time()-t1,'s'
    return meta

def solpnt_trajectory(tsolpnt):
    ''' Creates TRAJECTORY files based on the SOLPNTCAL start date/time and reads them.
    
        Innput:
           tsolpnt    Time object representing the start time of the SOLPNTCAL scan
        
        Returns:
           trajdict   A dictionary encapsulating the trajectory offsets.
    '''
    t1 = time()
    aa = eovsa_array()
    res = make_trajtables('Sun',aa,'caltraj',tsolpnt.mjd)
    trjradec = open('/tmp/caltraj.radec','r')
    trjlines = trjradec.readlines()
    trjradec.close()
    trjrao = []
    trjdeco = []
    for line in trjlines:
        rao, deco, junk = line.split(' ')
        trjrao.append(int(rao))
        trjdeco.append(int(deco))
    trjazel = open('/tmp/caltraj.azel','r')
    trjlines = trjazel.readlines()
    trjazel.close()
    trjazo = []
    trjelo = []
    for line in trjlines:
        azo, elo, junk = line.split(' ')
        trjazo.append(int(azo))
        trjelo.append(int(elo))
    # Special case of antenna 12
    trjradec = open('/tmp/caltraj12.radec','r')
    trjlines = trjradec.readlines()
    trjradec.close()
    trjrao12 = []
    trjdeco12 = []
    for line in trjlines:
        rao, deco, junk = line.split(' ')
        trjrao12.append(int(rao))
        trjdeco12.append(int(deco))
    print 'solpnt_trajectory took',time()-t1,'s'
    return {'time':tsolpnt,'trjrao':trjrao,'trjdeco':trjdeco,'trjrao12':trjrao12,'trjdeco12':trjdeco12,'trjazo':trjazo,'trjelo':trjelo}

def solpnt_addmask(meta,trajdict):
    ''' Calculates the mask which, where True, indicates good SOLPNTCAL data
        (i.e. antennas are tracking the correct location).
        
        Inputs:
           meta       The SQL metadata dictionary for the scan (returned from solpnt_meta())
           trajdict   The TRAJECTORY information of offsets (returned from solpnt_trajextory())
           
        Returns:
           mask       The array of size (420, 13, 28) indicating the good SOLPNTCAL data
    '''
    t1 = time()
    npts = len(trajdict['trjrao'])
    mask = np.zeros((nsec,nant,npts),dtype='bool')
    for i in range(nsec):
        # Loop over each RAO,DECO position and calculate "good data" mask for each
        # desired position
        for j in range(npts):
            # True if RAO and DECO matches this line (within 0.01 degrees)
            m1 = np.logical_and(abs(meta['rao'][i,antidx]  - trajdict['trjrao'][j])<100.,
                                abs(meta['deco'][i,antidx] - trajdict['trjdeco'][j])<100.)
            m2 = np.logical_and(abs(meta['azo'][i,antidx]  - trajdict['trjazo'][j])<100.,
                                abs(meta['elo'][i,antidx]  - trajdict['trjelo'][j])<100.)
            m3 = np.logical_and(abs(meta['rao'][i,antidx]  - trajdict['trjrao12'][j])<100.,
                                abs(meta['deco'][i,antidx] - trajdict['trjdeco12'][j])<100.)
            # True if above is true AND antenna is tracking
            mask[i,oldidx,j] = np.logical_and(m1[oldidx],meta['trk'][i,oldidx])
            mask[i,newidx,j] = np.logical_and(m2[newidx],meta['trk'][i,newidx])
            mask[i,ant12idx,j] = np.logical_and(m3[ant12idx],meta['trk'][i,ant12idx])
    ra0 = np.median(meta['ra'])
    dec0 = np.median(meta['dec'])
    az0, el0 = radec2azel(ra0, dec0, Time(meta['Timestamp'],format='lv'))
    meta.update({'mask':mask, 'ra0':ra0, 'dec0':dec0, 'az0':az0, 'el0':el0})
    print 'solpnt_addmask took',time()-t1,'s'

def attncal_from_solpnt_to_sql(out, do_plot=False):
    ''' Examines data from a "new" SOLPNTCAL scan to make sure it includes
        GAINCAL data, and if so, calculates the attenuation differences for 
        the various FEMATTN states 1-8 relative to FEMATTN state 0, and optionally 
        plots the results for states 1 and 2 (the most commonly used).
        
        After a successful analysis, writes the result as an ATTNCAL SQL record.
    '''
    from util import get_idbdir, fname2mjd, nearest_val_idx
    import attncal as ac
    import socket
    import dbutil

    if len(out['time']) < 500:
        # This scan is too short--probably taken before the era of GAINCALs as part of SOLPNT,
        # or other data problem.  Return an empty dictionary.
        return None

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

    # Get time from data (at 400th index, which should be the start of the GAINCAL data, 
    # and read 120 records of attn state from SQL database
    t = Time(out['time'][395],format='jd')
    cursor = dbutil.get_cursor()
    d15 = dbutil.get_dbrecs(cursor, dimension=15, timestamp=t, nrecs=120)
    cursor.close()
    # Find time indexes of the 62 dB attn state
    # Uses only ant 1 assuming all are the same
    dtot = (d15['Ante_Fron_FEM_HPol_Atte_Second'] + d15['Ante_Fron_FEM_HPol_Atte_First'])[:,0]
    # Use system clock day number to identify bad SQL entries and eliminate them
    good, = np.where(d15['Ante_Cont_SystemClockMJDay'][:,0] != 0)
    #import pdb; pdb.set_trace()
    # Indexes into SQL records where a transition occurred.
    transitions, = np.where(dtot[good] - np.roll(dtot[good],1) != 0)
    # Eliminate any zero-index transition (if it exists)
    if transitions[0] == 0:
        transitions = transitions[1:]
    # These now have to be translated into indexes into the data, using the times
    idx = nearest_val_idx(d15['Timestamp'][good,0][transitions],Time(out['time'],format='jd').lv)
    #import pdb; pdb.set_trace()
    vx = np.nanmedian(out['p'][:13,:,:,np.arange(idx[0]+1,idx[1]-1)],3)
    va = np.mean(out['a'][:13,:2,:,np.arange(idx[0]+1,idx[1]-1)],3)
    vals = []
    attn = []
    for i in range(1,10):
        try:
            vals.append(np.nanmedian(out['p'][:13,:,:,np.arange(idx[i]+1,idx[i+1]-1)],3) - vx)
            attn.append(np.log10(vals[0]/vals[-1])*10.)
        except IndexError:
            if i > 6:
                print 'Warning: GAINCAL for dB step',i,'missing. Nominal 2dB attenuation added.'
                attn.append(attn[-1]+2.)  # The gaincal was short.  Just add 2 dB to previous.
            else:
                print 'Error: GAINCAL scan is present but mangled.'
                raise
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
    outdict = {'time': Time(out['time'][0],format='jd'),'fghz': out['fghz'], 
                        'rcvr_auto':va, # 'attna': np.array(attna[1:]), 
                        'rcvr':vx, 'attn': np.array(attn[1:])}

    attn = ac.read_attncal(t)   # Attempt to read attn from SQL
    if attn[0]['time'].iso == outdict['time'].iso:
        print 'Attenuation record already exists in SQL--not updated.'
    else:
        print 'Writing attenuations to SQL for',outdict['time'].iso
        ch.fem_attn_val2sql([outdict])

    return outdict


def solpnt_read(meta):
    ''' Reads the data from the disk and removes the background.

        Inputs:
           meta       The SQL metadata dictionary for the scan (returned from solpnt_meta())
           
        Returns:
           otp        A dictionary containing background-subtracted data
           Also adds fghz to meta.
    '''
    warnings.filterwarnings("ignore", category=RuntimeWarning)
    #import pdb; pdb.set_trace()
    t1 = time()
    trange = Time([meta['Timestamp'],meta['Timestamp']+nsec],format='lv')
    # Determine skycal from the solpntcal file itself
    filename = read_idb.get_trange_files(trange)
    out = read_idb.readXdata(filename[0],desat=True)
    nt = len(out['time'])
    # This is the number of accumulations for each frequency (should be the same for all
    # antennas, polarizations, and times, but this does a median over all just to be
    # sure there are no glitches.
    mval = np.nanmedian(np.nanmedian(np.nanmedian(out['m'],3),0),0)
    # Conversion factor to scale power to autocorrelation
    # NB: The factor 102 was determined empirically and will likely change if the 
    # equalizer coefficient (current 2) is changed. (2023-Nov-23)
    p2afac = np.max(mval)/(102*mval)
    tsecsql = meta['tstamps']
    # The miriad times are slightly off from exact times, so round to the msec
    tsecdata = (Time(out['time']-2400000.5,format='mjd').lv + 0.001).astype('int64').astype('float')
    sqlidx, dataidx = common_val_idx(tsecsql, tsecdata)
    bg = np.zeros_like(out['p'][:nant,:,:,0])
    bgauto = np.zeros_like(out['p'][:nant,:,:,0])
    for ant in range(nant):
        bg[ant] = np.nanmedian(out['p'][ant,:,:,dataidx[np.where(meta['bgmask'][sqlidx,ant])]],0)
        bgauto[ant,0] = bg[ant,0]*p2afac
        bgauto[ant,1] = bg[ant,1]*p2afac
    skycal = {'timestamp':meta['Timestamp'], 'fghz':out['fghz'], 'rcvr_bgd':bg, 'rcvr_bgd_auto':bgauto}
    # Check if skycal record exists in SQL database, and write it if not
    xml, buf = ch.read_cal(13, trange[0])
    sqlmjd = Time(stateframe.extract(buf, xml['Timestamp']), format='lv').mjd  # Date of existing skycal record
    if int(sqlmjd) != int(trange[0].mjd):
        # Day numbers of existing and new records do not agree, so write a new record
        ch.skycal2sql(skycal)
        print 'Skycal record written to SQL.'
    else:
        print 'Skycal record already exists in SQL database--not updated.'
    # Analyze current gain state (this also writes attncal record to SQL if needed, 
    # so return value is not used)
    attn = attncal_from_solpnt_to_sql(out)
    if attn is None:
        print 'Warning: SOLPNTCAL has no appended GAINCAL'
    cout = gc.apply_fem_level(out, skycal=skycal)
    # Subtract background receiver level from data (already done in apply_fem_level...)
    #sna, snp, snf = bg.shape
    #bgd  =  bg.repeat(nt).reshape((sna,snp,snf,nt))
    #out['p'][:13] -= bgd
    #sna, snp, snf = bgauto.shape
    #bgdauto  =  bgauto.repeat(nt).reshape((sna,snp,snf,nt))
    #out['a'][:13,:2] -= bgdauto
    meta.update({'fghz':out['fghz']})
    print 'solpnt_read took',time()-t1,'s'
    return {'source':out['source'], 'fghz':out['fghz'], 'ut_mjd':out['time']-2400000.5, 
            'tsys':cout['p'], 'auto':np.real(cout['a'][:,:2]), 'dataidx':dataidx, 
            'sqlidx':sqlidx, 'skycal':skycal, 'p2afac':p2afac}


def solpnt_applymask(otp, meta):
    ''' Isolate the good data and obtain the median for each pointing
    
        Inputs:
           otp         The dictionary containing background-subtracted data (returned from solpnt_read())
           meta        The SQL metadata dictionary for the scan (returned from solpnt_addmask())
           
        Returns nothing, but it updates the otp dictionary to add tsyspts key
    '''
    t1 = time()
    nallants, npol, nf, nt = otp['tsys'].shape   # Note, nallants is 16 here, but nant remains 13
    npol = 2
    dataidx = otp['dataidx']
    sqlidx = otp['sqlidx']
    tsys0 = otp['tsys'][:nant,0,:,dataidx]               # Size nt, 13, 451
    tsys1 = otp['tsys'][:nant,1,:,dataidx]
    m = meta['mask'][sqlidx]   # Size nt, 13, 28
    npt = len(m[0,0,:])
    tsyspts = np.zeros([nf,npt,nant,npol],'float')    # The actual fluxes for each frequency, offset, and antenna
    # Step through pointings, applying mask
    for ipnt in range(npt):
        for iant in range(nant):
            good = np.where(m[:,iant,ipnt])[0]
            if len(good) == 0:
                # This pointing location had no good tracking, so set values to Nan
                tsyspts[:,ipnt,iant,:] = np.nan
            else:
                # Set values to median of good-tracking period
                tsyspts[:,ipnt,iant,0] = np.median(tsys0[good,iant,:],0)
                tsyspts[:,ipnt,iant,1] = np.median(tsys1[good,iant,:],0)
    otp.update({'tsyspts':tsyspts})
    print 'solpnt_applymask took',time()-t1,'s'

def gausfit(X,data, bounds=None):
    ''' Given a set of locations X and amplitude data at each location
        fit the data with a gaussian and return the fit parameters and
        x,y arrays of smooth, fitted data '''
    from scipy.optimize import curve_fit
    def func(x, a, b, c, d):
        return a*np.exp(-((x-b)/c)**2) + d
    if bounds is None:
        bounds=(-np.inf,np.inf)

    # This dies horribly if there are nans in the data, so remove them first
    X = X[~np.isnan(data)]
    data = data[~np.isnan(data)]
    if len(X) < 4:
        # Do not attempt fit if too many nans--just return all zeros
        return [0.0,0.0,0.0,0.0],np.zeros(100),np.zeros(100)
    xmax = X.max()
    xmin = X.min()
    try:
        popt, pcov = curve_fit(func, X, data, bounds=bounds)
    except RuntimeError:
        # Parameters not found--just return all zeros
        return [0.0,0.0,0.0,0.0],np.zeros(100),np.zeros(100)
    x = np.linspace(xmin,xmax,100)
    y = func(x, *popt)
    #plt.plot(x,peval(x, plsq[0]), X, data, 'o')
    return popt,x,y

def solpnt_gaussfit(otp):
    ''' Perform Gaussian fits to the SOLPNT data for each antenna, frequency, and 
        polarization, and return the parameters of the fit (and the flux measurement
        at each offset).
        
        Inputs:
           otp         The dictionary containing background-subtracted data (returned from 
                         solpnt_read() after solpnt_applymask() has been run).
        
        Outputs:
           params      The dictionary of fit parameters and flux measurements at each offset.
                         Note that px and py are the list of 4 Gaussian fit parameters for each
                         frequency and antenna, for the narrower beam direction, and are far
                         more accurate.  pcrossx and pcrossy are the same for the wider beam
                         direction, for completeness, but should not be used for quantitative
                         measurement.  Likewise, xo and yo are the flux measurements in the 
                         narrower beam direction, and xcrosso, ycrosso are the corresponding 
                         measurements in the wider beam direction.  The fits parameter order is
                         px[0]   Amplitude (arbitrary units) of the Gaussian.  For a well-pointed
                                   dish this is the Flux density of the Sun.
                         px[1]   Offset in 10000ths of a degree of the center of the Gaussian 
                                   along the sweep direction.
                         px[2]   Width in 10000ths of a degree of the Gaussian along the sweep
                                   direction.
                         px[3]   Constant (off-Sun) background on which the Gaussian sits.
    '''
    warnings.filterwarnings("ignore", category=RuntimeWarning)
    t1 = time()
    fghz, a, aout = disk_conv(otp['fghz'])
    aout *= 10000./(2*np.sqrt(np.log(2)))

    xpts = np.array([-100000, -50000, -20000, -10000, -5000, -2000, -1000, 0, 1000, 2000, 5000, 10000, 20000, 50000])
    tsyspts = otp['tsyspts']
    nf, npnt, nant, npol = tsyspts.shape
    px = np.zeros([4,nf,nant],'float')  # Array of gaussian fits
    pcrossx = np.zeros([4,nf,nant],'float')  # Array of gaussian fits
    crossxo = tsyspts[:,:14,:,0]       # First 14 pointings are Cross-X Offsets
    xo = tsyspts[:,14:,:,0]            # Second 14 pointings are X Offsets
    for iant in range(nant):
        xdata = xpts
        for ifreq in range(nf):
            ydata = xo[ifreq,:,iant]
            allzero = len(ydata.nonzero()[0]) == 0  # True if all data are zero
            ymax = np.nanmax(ydata)
            ymin = np.nanmin(ydata)
            yrange = ymax - ymin
            low_bounds =  [0.1*yrange, -30000., aout[ifreq]*0.5, ymin - 0.1*yrange]
            high_bounds = [1.0*yrange,  30000., aout[ifreq]*3.0, ymin + 0.1*yrange]    
            if allzero:
                high_bounds[0] = 1  # Ensure that high bounds are greater than low bounds
                high_bounds[3] = 1
            bounds = (low_bounds,high_bounds)
            #if iant == 1 and ifreq == 30:
            #    import pdb; pdb.set_trace()
            px[:,ifreq,iant], x, y = gausfit(xdata, ydata, bounds=bounds)

            ydata = crossxo[ifreq,:,iant]
            allzero = len(ydata.nonzero()[0]) == 0  # True if all data are zero
            ymax = np.nanmax(ydata)
            ymin = np.nanmin(ydata)
            yrange = ymax - ymin
            # width bounds are 0.9 to 1.1 nominally
            low_bounds =  [0.1*yrange, -30000., aout[ifreq]*0.5, ymin - 0.1*yrange]
            high_bounds = [1.0*yrange,  30000., aout[ifreq]*3.0, ymin + 0.1*yrange]    
            if allzero:
                high_bounds[0] = 1  # Ensure that high bounds are greater than low bounds
                high_bounds[3] = 1
            bounds = (low_bounds,high_bounds)            
            pcrossx[:,ifreq,iant], x, y = gausfit(xdata, ydata, bounds=bounds)
    # Repeat for other feed
    py = np.zeros([4,nf,nant],'float')  # Array of gaussian fits
    pcrossy = np.zeros([4,nf,nant],'float')  # Array of gaussian fits
    yo = tsyspts[:,:14,:,1]       # First 14 pointings are Y Offsets
    crossyo = tsyspts[:,14:,:,1]  # Second 14 pointings are Cross-Y Offsets
    for iant in range(nant):
        xdata = xpts
        for ifreq in range(nf):
            ydata = crossyo[ifreq,:,iant]
            allzero = len(ydata.nonzero()[0]) == 0  # True if all data are zero
            ymax = np.nanmax(ydata)
            ymin = np.nanmin(ydata)
            yrange = ymax - ymin
            low_bounds =  [0.1*yrange, -30000., aout[ifreq]*0.5, ymin - 0.1*yrange]
            high_bounds = [1.0*yrange,  30000., aout[ifreq]*3.0, ymin + 0.1*yrange]    
            if allzero:
                high_bounds[0] = 1  # Ensure that high bounds are greater than low bounds
                high_bounds[3] = 1
            bounds = (low_bounds,high_bounds)
            #if iant == 1 and ifreq == 30:
            #    import pdb; pdb.set_trace()
            pcrossy[:,ifreq,iant], x, y = gausfit(xdata, ydata, bounds=bounds)

            ydata = yo[ifreq,:,iant]
            allzero = len(ydata.nonzero()[0]) == 0  # True if all data are zero
            ymax = np.nanmax(ydata)
            ymin = np.nanmin(ydata)
            yrange = ymax - ymin
            # width bounds are 0.9 to 1.1 nominally
            low_bounds =  [0.1*yrange, -30000., aout[ifreq]*0.5, ymin - 0.1*yrange]
            high_bounds = [1.0*yrange,  30000., aout[ifreq]*3.0, ymin + 0.1*yrange]    
            if allzero:
                high_bounds[0] = 1  # Ensure that high bounds are greater than low bounds
                high_bounds[3] = 1
            bounds = (low_bounds,high_bounds)            
            py[:,ifreq,iant], x, y = gausfit(xdata, ydata, bounds=bounds)
    params = {'px':px, 'pcrossx':pcrossx, 'py':py, 'pcrossy':pcrossy, 
              'xo':xo, 'crossxo':crossxo, 'yo':yo, 'crossyo':crossyo}            
    print 'solpnt_gaussfit took',time()-t1,'s'
    return params

def solpnt_xanal(tsolpnt):
    ''' Complete analysis of a SOLPNTCAL scan, combining all of the routines in the proper
        order.
        
        Inputs:
           tsolpnt      Time object representing the start time of the SOLPNTCAL scan
           
        Outputs:
           solpnt       A dictionary containing all of the products of the SOLPNTCAL analysis.
    '''
    from copy import deepcopy
    meta = solpnt_meta(tsolpnt)
    trajdict = solpnt_trajectory(tsolpnt)
    solpnt_addmask(meta, trajdict)
    otp = solpnt_read(meta)
    solpnt_applymask(otp, meta)
    params = solpnt_gaussfit(otp)
    mjd = Time(meta['Timestamp'].astype(int),format='lv').mjd
    if mjd < Time('2024-05-22 21:10').mjd:
        # Prior to this date, there was an error in the trajectories for
        # azel antennas 1-8 (not 12), where azoff was really cross-eloff.
        # This is a rough correction for the fit parameters for those antennas.
        params['px'][1:2,:,:8] /= np.cos(meta['el0'])
        params['py'][1:2,:,:8] /= np.cos(meta['el0'])
        params['pcrossx'][1:2,:,:8] /= np.cos(meta['el0'])
        params['pcrossy'][1:2,:,:8] /= np.cos(meta['el0'])
    autoparams = deepcopy(params)
    for ant in antidx:
        # Create autoparams by scaling flux parameters from params
        autoparams['px'][0,:,ant] *= otp['p2afac']
        autoparams['px'][3,:,ant] *= otp['p2afac']
        autoparams['pcrossx'][0,:,ant] *= otp['p2afac']
        autoparams['pcrossx'][3,:,ant] *= otp['p2afac']
        autoparams['py'][0,:,ant] *= otp['p2afac']
        autoparams['py'][3,:,ant] *= otp['p2afac']
        autoparams['pcrossy'][0,:,ant] *= otp['p2afac']
        autoparams['pcrossy'][3,:,ant] *= otp['p2afac']
    solpnt = {'meta':meta, 'trajdict':trajdict, 'otp':otp, 
              'params':params, 'autoparams':autoparams}
    return solpnt

def solpnt_offsets(inparams, meta=None, savefig=True):
    ''' Calculates offsets in RA, Dec, Az, and El from the X pattern offsets, and
        optionally plots them.
        
        Inputs:
           inparams  Parameters of the Gaussian fits (returned from solpnt_gaussfit()).
           meta      The SQL metadata dictionary for the scan (returned from solpnt_read()).
                       If meta is None, params is interpreted as a solpnt dictionary as
                       output from solpnt_xanal()
           savefig  If true, saves the offsets as a solar map of mean over frequency, for 
                       all antennas, to the web folder /common/webplots/PTG/.
        
        Returns:
           offsets   A dictionary of offsets for each solar antenna.  NB: Old (equatorial) 
                       antennas use raoff and decoff, while new (azel) antennas use azoff
                       and eloff.
    '''
    t1 = time()
    import matplotlib.pylab as plt
    if meta is None:
        params = inparams['params']
        meta = inparams['meta']
    else:
        params = inparams
    px = params['px']
    py = params['py']
    fghz = meta['fghz']
    mjd = Time(solout['meta']['Timestamp'].astype(int),format='lv').mjd
    f, ax = plt.subplots(1,1)
    xoff = np.zeros(nant)
    yoff = np.zeros(nant)
    dx = np.zeros(nant)
    dy = np.zeros(nant)
    for i in antidx:
        bestx = np.sort(px[1,:,i].flatten())[50:400]
        xoff[i] = np.median(bestx) # Median x-feed offset over best 350 frequencies
        dx[i] = np.std(bestx)  # Std error in mean of x-feed offset over best 350 frequencies
        besty = np.sort(py[1,:,i].flatten())[50:400]
        yoff[i] = np.median(besty) # Median y-feed offset over best 350 frequencies
        dy[i] = np.std(besty)  # Std error in mean of y-feed offset over best 350 frequencies
    cp4 = np.cos(np.pi/4)
    cross_eloff = (xoff - yoff)*cp4
    cross_decoff = -cross_eloff
    eloff = -(xoff + yoff)*cp4   # Positive sweep offsets are negative elevation offsets and vice versa
    decoff = eloff
    azoff = cross_eloff/np.cos(meta['el0'])
    raoff = cross_decoff/np.cos(meta['dec0'])
    for i in antidx:
        ax.plot((px[1,200:300,i]-py[1,200:300,i])*cp4,
            -(px[1,200:300,i]+py[1,200:300,i])*cp4,'.',alpha=0.2)
    for i in antidx:
        ax.plot(cross_eloff[i],eloff[i],'k*')
        ax.plot([cross_eloff[i]-dx[i]*cp4,cross_eloff[i]+dx[i]*cp4],
                [eloff[i]-dx[i]*cp4,eloff[i]+dx[i]*cp4],'k')
        ax.plot([cross_eloff[i]-dy[i]*cp4,cross_eloff[i]+dy[i]*cp4],
                [eloff[i]+dy[i]*cp4,eloff[i]-dy[i]*cp4],'k')
        ax.text(cross_eloff[i],eloff[i]+300,str(i+1),fontsize=14,ha='center')
    tstr = Time(meta['Timestamp'],format='lv').iso[:19]
    ax.set_title('Pointing Offsets for '+tstr)
    ax.set_xlabel('Az or RA Offset')
    ax.set_ylabel('El or Dec Offset')
    th = np.linspace(0, 2*np.pi, 360)
    ax.fill(2500*np.cos(th),2500*np.sin(th),color='lightyellow')
    ax.plot(2500*np.cos(th),2500*np.sin(th),color='r')
    ax.axis('square')
    ax.set_xlim(-5000,5000)
    ax.set_ylim(-5000,5000)
    if savefig:
        tstr = tstr.replace('-','').replace(':','').replace(' ','')[:14]
        plt.savefig('/common/webplots/PTG/PTG'+tstr+'.png',bbox_inches='tight')

    f, ax = plt.subplots(4,4,figsize=[15,8.5])
    ax = ax.flatten()
    for ant in antidx:
        ax[ant].set_ylim(-5000,5000)
        if ant in oldidx:
            ax[ant].set_title('Ant '+str(ant+1)+' [blue=RA, red=Dec]')
            xdecoff = (px[1,:,ant]-py[1,:,ant])*cp4
            raoff_ = -xdecoff/np.cos(meta['dec0']) # NB: RA has opposite sign
            decoff_ = -(px[1,:,ant]+py[1,:,ant])*cp4
            ax[ant].plot(fghz, raoff_,'.',color='skyblue')  
            ax[ant].plot(fghz, decoff_,'.',color='salmon')
            ax[ant].plot([fghz[0],fghz[-1]],[raoff[ant],raoff[ant]],'b')
            ax[ant].plot([fghz[0],fghz[-1]],[decoff[ant],decoff[ant]],'r')
        else:
            ax[ant].set_title('Ant '+str(ant+1)+' [blue=Az, red=El]')
            xeloff = (px[1,:,ant]-py[1,:,ant])*cp4
            azoff_ = xeloff/np.cos(meta['el0'])
            eloff_ = -(px[1,:,ant]+py[1,:,ant])*cp4
            ax[ant].plot(fghz, azoff_,'.',color='skyblue')
            ax[ant].plot(fghz, eloff_,'.',color='salmon')
            ax[ant].plot([fghz[0],fghz[-1]],[azoff[ant],azoff[ant]],'b')
            ax[ant].plot([fghz[0],fghz[-1]],[eloff[ant],eloff[ant]],'r')
        if ant % 4 == 0:
            ax[ant].set_ylabel('Offset [10000th-deg]')
    ax[12].set_xlabel('Frequency [GHz]')
    ax[13].set_xlabel('Frequency [GHz]')
    # I am hazy about why the raoff needs to be opposite sign, but experience shows it's true.
    offsets = {'raoff':-raoff, 'decoff':decoff, 'azoff':azoff, 'eloff':eloff, 'Timestamp':meta['Timestamp']}
    print 'solpnt_offsets took',time()-t1,'s'
    return offsets
    
def offsets2ants(offsets,ant_str=None):
    ''' Given the list of offsets output by solpnt_offsets() for 13 antennas, convert to 
        pointing coefficients, update the existing coefficients, and send to the relevant 
        antennas.  The antennas to update are specified with ant_str (defaults to no 
        antennas, for safety). 

        Inputs:
           offsets      Dictionary of offsets returned by solpnt_offsets().
           ant_str      The standard antenna list string (e.g. 'ant1-3 ant5-8') giving the
                          antenna list to be updated.

        Outputs:
           None, but commands are sent to each antenna indicated by ant_str, to update 
           their pointing coefficients 1 and 7.        
    '''
    t1 = time()
    from calibration import send_cmds
    if ant_str is None:
        print 'No antenna list specified, so there is nothing to do!'
        return

    timestamp = offsets['Timestamp']
    antlist = ant_str2list(ant_str)
    if antlist is None:
        return
    cursor = dbutil.get_cursor()
    # Read 10 stateframe records at timestamp of SOLPNTCAL and find first good ones
    D15data = dbutil.get_dbrecs(cursor,dimension=15,timestamp=timestamp,nrecs=10)
    p1_cur = []
    p7_cur = []
    good_antlist = []
    for ant in antlist:
        # Check that current pointing parameters are good (signified by non-zero clock)
        good, = np.where(D15data['Ante_Cont_SystemClockms'][:,ant] != 0)
        if good.tolist() != []:
            p1_cur.append(D15data['Ante_Cont_PointingCoefficient1'][good[0],ant])
            p7_cur.append(D15data['Ante_Cont_PointingCoefficient7'][good[0],ant])
            good_antlist.append(ant)
        else:
            p1_cur.append(np.nan)
            p7_cur.append(np.nan)
            print 'Ant',str(ant+1),'clock is zero--pointing coefficients cannot be updated.'
    if good_antlist != []:
        accini = stateframe.rd_ACCfile()
        acc = {'host': accini['host'], 'scdport':accini['scdport']}
        for i in good_antlist:
            j, = np.where(i == antlist)[0]
            if i in oldidx:
                p1_inc = -int(round(offsets['raoff'][i]))   # These apparently need to be negative
                p7_inc = int(round(offsets['decoff'][i]))
            else:
                p1_inc = -int(round(offsets['azoff'][i]))   # These apparently need to be negative
                p7_inc = int(round(offsets['eloff'][i]))
            p1_new = p1_cur[j] + p1_inc
            p7_new = p7_cur[j] + p7_inc
            print 'Updating P1 for Ant',i+1,'P1_old =',p1_cur[j],'P1_inc =',p1_inc,'P1_new =',p1_new
            cmd1 = 'pointingcoefficient1 '+str(p1_new)+' ant'+str(i+1)
            print 'Updating P7 for Ant',i+1,'P7_old =',p7_cur[j],'P7_inc =',p7_inc,'P7_new =',p7_new
            cmd7 = 'pointingcoefficient7 '+str(p7_new)+' ant'+str(i+1)
            print 'Commands to be sent:'
            print cmd1
            print cmd7,'\n'
            send_cmds([cmd1],acc)
            send_cmds([cmd7],acc)
    print 'offsets2ants',time()-t1,'s'

def solpnt_bsize(inparams, meta=None):
    ''' Makes a plot of the beam size from the SOLPNTCAL measurements
        
        Inputs:
           params    Parameters of the Gaussian fits (returned from solpnt_gaussfit()).
           meta      The SQL metadata dictionary for the scan (returned from solpnt_read()).
                       If meta is None, params is interpreted as a solpnt dictionary as
                       output from solpnt_xanal()
        
        Returns:
           None, but a plot is made of the beam size vs. frequency for each antenna and feed.
    '''
    import matplotlib.pylab as plt
    if meta is None:
        params = inparams['params']
        meta = inparams['meta']
    else:
        params = inparams
    fghz, a, aout = disk_conv(meta['fghz'])
    aout *= 10000./(2*np.sqrt(np.log(2)))
    px = params['px']
    py = params['py']
    f, ax = plt.subplots(4,4,figsize=[15,8.5])
    ax = ax.flatten()
    for ant in antidx:
        ax[ant].set_ylim(2000,100000)
        ax[ant].set_yscale('log')
        ax[ant].set_title('Ant '+str(ant+1)+' [blue=X-feed, red=Y-feed]')
        xwid = px[2,:,ant]
        ywid = py[2,:,ant]
        ax[ant].plot(fghz, xwid,'.',color='skyblue')  
        ax[ant].plot(fghz, ywid,'.',color='salmon')
        ax[ant].plot(fghz, aout, color='k')
        if ant % 4 == 0:
            ax[ant].set_ylabel('Beam Width [10000th-deg]')
    ax[12].set_xlabel('Frequency [GHz]')
    ax[13].set_xlabel('Frequency [GHz]')

def params2calfac(s, params, fghz):
    ''' Helper routine for solpnt_calfac, to repeat all steps for both total power
        and autocorrelation.
    '''
    fg, a, aout = disk_conv(fghz)
    aout *= 10000./(2*np.sqrt(np.log(2)))
    a *= 10000.

    nfrq = len(fghz)
    nant = len(antidx)
    calfac = np.zeros((2,nfrq,nant),'float')
    offsun = np.zeros((2,nfrq,nant),'float')
    for ant in antidx:
        # Do flux calculation for X feed
        s1 = params['px'][0,:,ant]      # X feed peak flux in arb. units
        s2 = params['pcrossx'][0,:,ant] # X feed cross sweep peak flux in arb. units
        o1 = params['px'][1,:,ant]      # X feed offset (1/10000th of a degree)
        o2 = params['pcrossx'][1,:,ant] # X feed cross sweep offset (1/10000th of a degree)
        # Use nominal beamwidth for narrow direction, and 2*nominal width for wide direction
        s01 = s1*np.exp((o2/(aout*2))**2)  # Corrects X feed flux for cross sweep off-pointing
        s02 = s2*np.exp((o1/aout)**2)      # Corrects X feed cross sweep flux for off-pointing
        s0 = (s01 + s02)/2.             # S01 should equal S02, take the mean
        s0[np.where(s0 == 0)] = np.nan
        calfac[0,:,ant] = s/s0          # sfu/unit
        offsun[0,:,ant] = (params['px'][3,:,ant] + params['pcrossx'][3,:,ant])/2  # arb. units

        # Repeat for Y feed
        s1 = params['py'][0,:,ant]      # Y feed peak flux in arb. units
        s2 = params['pcrossy'][0,:,ant] # Y feed cross sweep peak flux in arb. units
        o1 = params['py'][1,:,ant]      # Y feed offset (1/10000th of a degree)
        o2 = params['pcrossy'][1,:,ant] # Y feed cross sweep offset (1/10000th of a degree)
        # Use nominal beamwidth for narrow direction, and 2*nominal width for wide direction
        s01 = s1*np.exp((o2/(aout*2))**2)  # Corrects Y feed flux for cross sweep off-pointing
        s02 = s2*np.exp((o1/aout)**2)      # Corrects Y feed cross sweep flux for off-pointing
        s0 = (s01 + s02)/2.             # S01 should equal S02, take the mean
        s0[np.where(s0 == 0)] = np.nan
        calfac[1,:,ant] = s/s0          # sfu/unit
        offsun[1,:,ant] = (params['py'][3,:,ant] + params['pcrossy'][3,:,ant])/2  # arb. units
    return calfac, offsun

def solpnt_calfac(inparams, do_plot=False, prompt=True):
    ''' Reads the RSTN/Penticton flux, fits to the observed frequencies, and applies
        them to the antenna solar response in input dictionaries x and y to return
        the calibration factors with which to MULTIPLY solar data to convert to solar
        flux units.  The offsun (background) spectrum for each polarization and antenna
        is also returned.  
        
        if do_plot is True, also make a nice plot of the factors applied to the data,
        for total power only (not autocorrelation)
        
        TODO: These need to be scaled for gain state
    '''
    import matplotlib.pylab as plt
    import rstn
    params = inparams['params']
    autoparams = inparams['autoparams']
    meta = inparams['meta']

    t = Time(meta['Timestamp'],format='lv')
    fghz = meta['fghz']
    frq, flux = rstn.rd_rstnflux(t)
    if frq is None:
        print 'Cannot read RSTN flux--cannot continue.'
        return None
    fghz = meta['fghz']
    fmhz = fghz*1000.
    nfrq = len(fghz)
    nant = len(antidx)
    s = rstn.rstn2ant(frq, flux, fmhz, t)

    if do_plot:
        # Set up summary plot
        f, ax = plt.subplots(4, 4)
        ax = ax.flatten()
        f.set_size_inches(15,8.5)
        f.suptitle('Calibration for SOLPNT scan at '+t.iso[:19]+' UT',fontsize=18)
        for ant in antidx:
            ax[ant].set_title('Ant '+str(ant+1)+' Solar Spectrum')
            if ant % 4 == 0:
                ax[ant].set_ylabel('Solar Flux [sfu]')
        ax[12].set_xlabel('Frequency [GHz]')
        ax[13].set_xlabel('Frequency [GHz]')

    tpcalfac, tpoffsun = params2calfac(s, params, fghz)
    accalfac, acoffsun = params2calfac(s, autoparams, fghz)   # I could replace this by scaling from params
    if do_plot:
        for ant in antidx:
            s1 = params['px'][0,:,ant]
            ax[ant].plot(fghz,s1*tpcalfac[0,:,ant],'.',color='skyblue',label='X-Feed')
            s1 = params['py'][0,:,ant]
            ax[ant].plot(fghz,s1*tpcalfac[1,:,ant],'.',color='salmon',label='Y-Feed')
            #ax[ant].plot(fghz,s2*calfac[1,:,ant],'.',color='khaki',label='Y-cross-Feed')
            ax[ant].set_ylim(0,600)
            ax[ant].set_xlim(0,19)
            ax[ant].legend(loc='upper left',fontsize='small')

    # The auto-correlation values are just copies of the total power values
    tpcal_dict = {'fghz':fghz,'timestamp':meta['Timestamp'],
                  'tpcalfac':tpcalfac,'accalfac':accalfac,'tpoffsun':tpoffsun,'acoffsun':acoffsun}
    tsql = Time(int(t.mjd) + 20/24.,format='mjd')
    if prompt:
        ok = raw_input('Okay to write result to SQL database? [Y/N]: ')
    else:
        ok = 'Y'
    if ok.upper() == 'Y':
        ch.tpcal2sql(tpcal_dict,t=tsql)
        print 'Result was written to the SQL database'
    else:
        print 'Result was NOT written to the SQL database'

    return tpcal_dict

def solpnt_check_offsets(solpnt,offsets,ant):
    ''' Plots actual measurements on the sky (circles with size proportional to signal strength), 
        derived offsets from the fit parameters (orange dots), and the offsets calculated by 
        solpnt_offsets() (black + symbols) to verify consistency.
    '''
    import matplotlib.pylab as plt
    f, ax = plt.subplots(3,3,figsize=(7,7))
    ax = ax.flatten()
    if ant in oldidx:
        xoff = solpnt['trajdict']['trjrao']
        yoff = solpnt['trajdict']['trjdeco']
        y0 = solpnt['meta']['dec0']
    else:
        xoff = solpnt['trajdict']['trjazo']
        yoff = solpnt['trajdict']['trjelo']
        y0 = solpnt['meta']['el0']
        xsign = 1
    antfac = abs(np.float64(xoff[0])/yoff[0])
    px = solpnt['params']['px'][1,200:300,ant]
    py = solpnt['params']['py'][1,200:300,ant]
    f.suptitle(Time(solpnt['meta']['Timestamp'],format='lv').iso+' Ant '+str(ant+1))
    for k,j in enumerate([200,210,220,230,240,250,260,270,280]):
        adx = solpnt['params']['xo'][j,:,ant]-np.nanmin(solpnt['params']['xo'][j,:,ant])
        ady = solpnt['params']['yo'][j,:,ant]-np.nanmin(solpnt['params']['yo'][j,:,ant])
        adxmax = np.nanmax(adx)
        adymax = np.nanmax(ady)
        for i in range(14):
            ax[k].plot(xoff[i], yoff[i],'.',markersize=adx[i]*30/adxmax,color='C0',alpha=0.5)
        for i in range(14,28):
            ax[k].plot(xoff[i], yoff[i],'.',markersize=ady[i-14]*30/adymax,color='C0',alpha=0.5)
        #ax[k].axis('square')
        ax[k].set_xlim(-5000*antfac,5000*antfac)
        ax[k].set_ylim(-5000,5000)
        if k in [6, 7, 8]:
            ax[k].set_xlabel('-Az or RA')
        if k in [0,3,6]:
            ax[k].set_ylabel('El or Dec')
        if ant in oldidx:
            ax[k].plot(-(px-py)*np.cos(np.pi/4)/np.cos(y0),-(px+py)*np.cos(np.pi/4),'.',color='C1')
            ax[k].plot(-offsets['raoff'][ant],offsets['decoff'][ant],'k+')
        else:
            ax[k].plot((px-py)*np.cos(np.pi/4)/np.cos(y0),-(px+py)*np.cos(np.pi/4),'.',color='C1')
            ax[k].plot(offsets['azoff'][ant],offsets['eloff'][ant],'k+')

def sp_check_qual(solpnt):
    ''' Does a sanity check for quality of SOLPNTCAL, based on comparison
        of measured beamsize with that expected from a 2.1-m dish.  For 
        some reason, the new feeds return a much narrower beam than 
        expected (around 0.8 of nominal).
    '''
    fghz = solpnt['meta']['fghz']
    fout,a,aout = disk_conv(fghz[:400])
    aout = aout*10000/(2*np.sqrt(np.log(2)))
    nparms,nf,nant = solpnt['params']['px'].shape
    qual = np.zeros((2,nant),'bool')
    for i in range(nant):
        # Beamsize must be within 10 percent of nominal value
        val = np.median(solpnt['params']['px'][2,:400,i]/aout)
        qual[0,i] = val < 1.0 and val > 0.7
        val = np.median(solpnt['params']['py'][2,:400,i]/aout)
        qual[1,i] = val < 1.0 and val > 0.7
    return qual

# Problem--the autocorr parameters are time-consuming to calculate and can be different than 
# the total power parameters, leading to an incorrect calibration.  It turns out that I can 
# predict the autocorr values from total power by

# solout['autoparams']['px'][0,:,ant] = solout['params']['px'][0,:,ant]/104./(outc['m'][ant,0,:,0]/745472.)

# where outc['m'] are the m values for the data.  

# In other words, these two commands give essentially the same plot:
#plot(solout['autoparams']['px'][0,:,2]+solout['autoparams']['px'][3,:,2],'.')
#plot((solout['params']['px'][0,:,2]+solout['params']['px'][3,:,2])/104./(outc['m'][2,0,:,0]/745472.),'.')
