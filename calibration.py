'''
   Module for accessing, analyzing and plotting calibration data'''
#
# History:
#   2014-Dec-11  DG
#     First written (to analyze SOLPNTCAL scans)
#   2014-Dec-12  DG
#     Lots of tweaks to make and beautify the plots.  Lots more to do!
#   2014-Dec-14  DG
#     Added delay after dump_tsys(), to allow for files to finish writing.
#     Also added error return for when RSTN data are not available.
#   2014-Dec-15  DG
#     Changed to work with UDB files as well as IDB files.
#   2014-Dec-16  DG
#     Added ability to run as cron job (__main__ call at the end).
#     Also added a routine sp_check_qual() to do quality/sanity check,
#     to verify that the calibration is okay.
#   2014-Dec-17  DG
#     Converted s0 zeros to np.nan so division no longer sets a warning.
#   2014-Dec-26  DG
#     Added write of today's calibration file to tomorrow, so that
#     tomorrow's initial scans will be calibrated.
#   2015-Jan-25  DG
#     Fixed all of the routines to work with any number of antennas (?)
#     At least it should work for 4 or 8...
#   2015-May-29  DG
#      Converted from using datime() to using Time() based on astropy.
#   2015-Jun-27  DG
#    Changed the code to standardize on names, content, and order of indices
#    of outputs for rd_miriad_tsys.  Names ut_mjd, fghz, and tsys will be used, 
#    with units hopefully made obvious.  Order of indices will be 
#          (nant/nbl, npol, nfreq, ntimes).
#    This required some changes to several of the routines here.
#   2016-May-05  DG
#    Add test in solpntanal() for 16-ant correlator data (based on date),
#    and read with rd_miriad_tsys_16() routine if so.
#   2016-May-07  DG
#    Fix some plotting bugs for when number of antennas is odd.
#   2016-Jun-15  DG
#    Added new versions sp_apply_cal2() and sp_bg_subtract2() to work
#    with new data format.
#   2016-Aug-07  DG
#    Change sp_offsets() to give RA and Dec offsets for old ants.  Also
#    wrote offsets2ants() routine (to replace one that I wrote earlier,
#    but somehow got lost...)
#   2016-Nov-20  DG
#    There was a slight bug in offsets2ants(), now fixed.
#   2016-Nov-21  DG
#    Had the wrong sign of the correction in offsets2ants()!   Yikes!  It
#    was only x that was the wrong sign. (I worry that this is only for the
#    equatorial dishes...)
#   2016-Nov-27  DG
#    New routine to analyze interferometric (cal) pointing, called calpntanal().
#   2016-Dec-09  DG
#    I had determined RA offsets in units of "cross-Dec," but in fact I want
#    them in units of RA, so do not multiply by cos(Dec).  Also, I now write
#    out HA and Dec in degrees instead of radians.
#

if __name__ == "__main__":
    ''' For cron-job use, make sure we use a backend that does not require
        X-server    
    '''
    import matplotlib
    matplotlib.use('Agg')

import numpy as np
import solpnt
from util import Time
import struct, time, glob, sys, socket
from disk_conv import *
import dump_tsys

def calpntanal(t,ant_str='ant1-13',do_plot=True):
    ''' Does a complete analysis of CALPNTCAL, reading information from the SQL
        database, finding the corresponding Miriad IDB data, and doing the 
        gaussian fit to the beam, to return the beam and offset parameters.
    '''
    import matplotlib.pyplot as plt
    import read_idb
    bl2ord = read_idb.p.bl_list()
    # Read pointing data (timerange t must be accurate)
    out = read_idb.read_idb(t, navg=30)
    # Determine wanted baselines with ant 14 from ant_str
    idx = read_idb.p.ant_str2list(ant_str)
    # Do appropriate sums over frequency and baseline
    pntdata = np.sum(np.abs(np.sum(out['x'][bl2ord[idx,13],0,:,:48],1)),0)
    # Measurements are 90 s long, hence 3 consecutive 30 s points, so do final
    # sum over these
    pntdata.shape = (16,3)
    stdev = np.std(pntdata,1)
    pntdata = np.sum(pntdata,1)
    # Multiply by cos(dec) to convert to "Cross-dec" sky coordinates -- NO, leave in RA units
    #rao = np.array([-1.00, -0.20, -0.10, -0.05, 0.00, 0.05, 0.10, 0.20])*np.cos(out['dec'])
    rao = np.array([-1.00, -0.20, -0.10, -0.05, 0.00, 0.05, 0.10, 0.20])
    radat = pntdata[:8]
    deco = np.array([-0.20, -0.10, -0.05, 0.00, 0.05, 0.10, 0.20, 1.00])
    decdat = pntdata[8:]
    plsqr, xr, yr = solpnt.gausfit(rao, radat)
    plsqd, xd, yd = solpnt.gausfit(deco, decdat)
    midtime = Time((t[1].lv + t[0].lv)/2.,format='lv')
    if (do_plot):
        f, ax = plt.subplots(1,2)
        f.set_size_inches(7,3.5,forward=True)
        ax[0].errorbar(rao,radat,yerr=stdev[:8],fmt='.')
        ax[0].plot(xr,yr)
        ax[0].axvline(x=0,color='k')
        ax[0].axvline(x=plsqr[1],linestyle='--')
        ax[1].errorbar(deco,decdat,yerr=stdev[8:],fmt='.')
        ax[1].plot(xd,yd)
        ax[1].axvline(x=0,color='k')
        ax[1].axvline(x=plsqd[1],linestyle='--')
        for j in range(2):
            ax[j].set_xlim(-0.3, 0.3)
            ax[j].grid()
        ax[0].text(0.6,0.9,'RAO :'+str(plsqr[1])[:5],transform=ax[0].transAxes)
        ax[0].text(0.6,0.85,'FWHM:'+str(plsqr[2])[:5],transform=ax[0].transAxes)
        ax[0].set_xlabel('RA Offset [deg]')
        ax[1].text(0.6,0.9,'DECO:'+str(plsqd[1])[:5],transform=ax[1].transAxes)
        ax[1].text(0.6,0.85,'FWHM:'+str(plsqd[2])[:5],transform=ax[1].transAxes)
        ax[1].set_xlabel('Dec Offset [deg]')
        f.suptitle('Pointing on '+out['source']+' at '+midtime.iso)
    return {'source':out['source'],'ha':out['ha'][24]*180./np.pi,'dec':out['dec']*180/np.pi,
               'rao':plsqr[1],'deco':plsqd[1],'time':midtime}
    
    
def solpntanal(t,udb=False):
    ''' Does a complete analysis of SOLPNTCAL, reading information from the SQL
        database, finding and dumping the corresponding Miriad IDB data, and 
        doing gaussian fit to beam to return the beam parameters.  The outputs
        are identical dictionaries for x and y polarizations.
    '''
    
    def fname2mjd(filename):
        fstem = filename.split('/')[-1]
        fstr = fstem[3:7]+'-'+fstem[7:9]+'-'+fstem[9:11]+' '+fstem[11:13]+':'+fstem[13:15]+':'+fstem[15:17]
        return Time(fstr).mjd

    pnt = solpnt.get_solpnt(t)
    proc = solpnt.process_solpnt(pnt)
    trange = Time([pnt['Timestamp'],pnt['Timestamp']+300.],format='lv')
    if trange[0].mjd < 57450:
        otp = dump_tsys.rd_miriad_tsys(trange,udb=udb)
    else:
        otp = dump_tsys.rd_miriad_tsys_16(trange,udb=udb)
#    if udb:
#        fstr = t1.iso
#        folder = '/data1/UDBTXT/'+fstr[:4]
#        files = glob.glob(folder+'/UDB'+fstr.replace('-','').split()[0]+'*_xtsys.txt')
#        files.sort()
#        for filename in files:
#            mjd = fname2mjd(filename)
#            if mjd >= t1.mjd:
#                xfile = filename
#                yfile = filename.replace('xtsys','ytsys')
#                sfile = filename.replace('xtsys','sfreq')
#                break
#    else:
#        xfiles, yfiles = solpnt.dmp_tsys(t)
#        # Give it some time for the files to be written and closed.
#        time.sleep(5)
#        xfile = xfiles[0]
#        yfile = yfiles[0]
#        sfile = '/common/tmp/txt/sfreq.txt'
#    otp = solpnt.rd_tsys(xfile,sfile)
    fghz = otp['fghz']
    xra,xdec,xrao,xdeco = solpnt.process_tsys(otp,proc,pol=0)
    x = {'ut_mjd':otp['ut_mjd'],'fghz':fghz,'ra0':proc['ra0'],'dec0':proc['dec0'],'raparms':xra,'decparms':xdec,'rao':xrao,'deco':xdeco}
#    otp = solpnt.rd_tsys(yfile,sfile)
    yra,ydec,yrao,ydeco = solpnt.process_tsys(otp,proc,pol=1)
    y = {'ut_mjd':otp['ut_mjd'],'fghz':fghz,'ra0':proc['ra0'],'dec0':proc['dec0'],'raparms':yra,'decparms':ydec,'rao':yrao,'deco':ydeco}
    qual = sp_check_qual(x,y)
    return x,y, qual
        
def sp_get_calfac(x,y, do_plot=True):
    ''' Reads the RSTN/Penticton flux, fits to the observed frequencies, and applies
        them to the antenna solar response in input dictionaries x and y to return
        the calibration factors with which to MULTIPLY solar data to convert to solar
        flux units.  The offsun (background) spectrum for each polarization and antenna
        is also returned.  
        
        if do_plot is True, also make a nice plot of the factors applied to the data
        
        TODO: These need to be scaled for gain state
    '''
    import rstn
    import matplotlib.pyplot as plt

    t = Time(x['ut_mjd'][0],format='mjd')
    frq, flux = rstn.rd_rstnflux(t)
    if frq is None:
        print 'Cannot continue.'
        return None
    nfrq, npnt, nant = x['rao'].shape
    fmhz = x['fghz']*1000.
    fghz = x['fghz']
    s = rstn.rstn2ant(frq, flux, fmhz, t)
    calfac = np.zeros((2,nfrq,nant),'float')
    offsun = np.zeros((2,nfrq,nant),'float')
    
    if do_plot:
        # Set up summary plot
        nrow = 2
        ncol = (nant+1)/2
        f, ax = plt.subplots(nrow, ncol, sharex='col', sharey='row')
        f.set_size_inches(2*nant,7,forward=True)
        f.suptitle('Calibration for SOLPNT scan at '+t.iso[:19]+' UT',fontsize=18)
        for ant in range(nant):
            ax[ant % nrow, ant/nrow].set_title('Ant '+str(ant+1)+' Solar Spectrum')
            if ant % ncol == 0:
                ax[ant/ncol,0].set_ylabel('Solar Flux [sfu]')
            if ant >= nant/nrow:
                ax[1,ant - nant/nrow].set_xlabel('Frequency [GHz]')

    for ant in range(nant):
        # Do flux calculation for X feed
        a1 = x['raparms'][2,:,ant]   # 1/e half-width of RA  beam  (1/10000th of a degree)
        a2 = x['decparms'][2,:,ant]  # 1/e half-width of Dec beam  (1/10000th of a degree)
        s1 = x['raparms'][0,:,ant]   # RA  peak flux in arb. units
        s2 = x['decparms'][0,:,ant]  # Dec peak flux in arb. units
        o1 = x['raparms'][1,:,ant]   # RA  offset (1/10000th of a degree)
        o2 = x['decparms'][1,:,ant]  # Dec offset (1/10000th of a degree)
        s01 = s1*np.exp((o2/a1)**2)  # Corrects RA flux for off-pointing in Dec
        s02 = s2*np.exp((o1/a2)**2)  # Corrects Dec flux for off-pointing in RA
        s0 = (s01 + s02)/2.          # S01 should equal S02, take the mean
        s0[np.where(s0 == 0)] = np.nan
        calfac[0,:,ant] = s/s0       # sfu/unit
        offsun[0,:,ant] = (x['raparms'][3,:,ant] + x['decparms'][3,:,ant])/2  # arb. units
        if do_plot:
            ax[ant % 2, ant/2].plot(fghz,s1*calfac[0,:,ant],'.',label='X-RA')
            ax[ant % 2, ant/2].plot(fghz,s2*calfac[0,:,ant],'.',label='X-Dec')
        # Repeat for Y feed
        a1 = y['raparms'][2,:,ant]
        a2 = y['decparms'][2,:,ant]
        s1 = y['raparms'][0,:,ant]
        s2 = y['decparms'][0,:,ant]
        o1 = y['raparms'][1,:,ant]
        o2 = y['decparms'][1,:,ant]
        s01 = s1*np.exp((o2/a1)**2)  # Corrects RA flux for off-pointing in Dec
        s02 = s2*np.exp((o1/a2)**2)  # Corrects Dec flux for off-pointing in RA
        s0 = (s01 + s02)/2.
        s0[np.where(s0 == 0)] = np.nan
        calfac[1,:,ant] = s/s0       # sfu/unit
        offsun[1,:,ant] = (y['raparms'][3,:,ant] + y['decparms'][3,:,ant])/2  # arb. units
        if do_plot:
            ax[ant % nrow, ant/nrow].plot(fghz,s1*calfac[1,:,ant],'.',label='Y-RA')
            ax[ant % nrow, ant/nrow].plot(fghz,s2*calfac[1,:,ant],'.',label='Y-Dec')
            ax[ant % nrow, ant/nrow].set_ylim(0,600)
            ax[ant % nrow, ant/nrow].set_xlim(0,19)
            ax[ant % nrow, ant/nrow].legend(loc='upper left',fontsize='small')
    return calfac,offsun

def sp_check_qual(x, y):
    ''' Does a sanity check for quality of SOLPNTCAL, based on comparison
        of measured beamsize with that expected from a 2.1-m dish
    '''
    fghz = x['fghz']
    fout,a,aout = disk_conv(fghz)
    aout = aout*10000/(2*np.sqrt(np.log(2)))
    nparms,nf,nant = x['raparms'].shape
    qual = np.zeros((4,nant),'bool')
    for i in range(nant):
        # Beamsize must be within 10 percent of nominal value
        val = np.median(x['raparms'][2,:,i]/aout)
        qual[0,i] = val < 1.1 and val > 0.9
        val = np.median(x['decparms'][2,:,i]/aout)
        qual[1,i] = val < 1.1 and val > 0.9
        val = np.median(y['raparms'][2,:,i]/aout)
        qual[2,i] = val < 1.1 and val > 0.9
        val = np.median(y['decparms'][2,:,i]/aout)
        qual[3,i] = val < 1.1 and val > 0.9
    return qual

def sp_write_calfac(x, y, calfac, offsun):
    ''' Given input dictionaries x, y, and the calibration factors calfac and offsun,
        generated by sp_get_calfac(), write standard binary files with the name
            TPCALyyyymmdd_a_bbb_c.dat, where a = number of feeds (2), 
                                           bbb = number of frequencies,  
                                             c = number of antennas.
        These values are needed in order to know how to read the file...
        
        The file contains a list of bbb frequencies, and the arrays calfac and offsun,
        both of dimensions [a,bbb,c].
    '''
    fghz = x['fghz']
    nf = len(fghz)
    t1 = Time(x['ut_mjd'][0],format='mjd')
    datstr = t1.iso[:10].replace('-','')
    dims = calfac.shape
    siz = np.prod(dims)
    buf = struct.pack(str(nf)+'f',*fghz)
    buf += struct.pack(str(siz)+'f',*calfac.reshape(siz))
    buf += struct.pack(str(siz)+'f',*offsun.reshape(siz))
    filename = '/data1/TPCAL/TPCAL'+datstr+'_'+str(dims[0])+'_'+str(dims[1])+'_'+str(dims[2])+'.dat'
    f = open(filename,'wb')
    f.write(buf)
    f.close()
    return
    
def sp_read_calfac(t):
    ''' Read the contents of a SOLPNT calibration file and return
        fghz, calfac, offsun arrays
    '''
    files = glob.glob('/data1/TPCAL/*')
    for file in files:
        if file.find(t.iso.replace('-','')[:8]) != -1:
            npol,nf,nant = file.replace('.dat','').split('_')[1:]
            nsiz = int(npol)*int(nf)*int(nant)
            nfi = int(nf)
            f = open(file,'rb')
            data = f.read()
            f.close()
            fghz = np.array(struct.unpack_from(nf+'f',data,0))
            calfac = np.array(struct.unpack_from(str(nsiz)+'f',data,nfi*4)).reshape(int(npol),int(nf),int(nant))
            offsun = np.array(struct.unpack_from(str(nsiz)+'f',data,(nfi+nsiz)*4)).reshape(int(npol),int(nf),int(nant))
            return fghz, calfac, offsun
    print 'Calibration file not found for date',t.iso
    return None, None, None

def sp_apply_cal2(out,calfac,offsun):
    ''' Given "standard" output of tp_display.rd_tsys_multi(),
        and corresponding calibration data from sp_read_calfac(),
        apply the calibration and return the calibrated output.
    '''
    fghz = out['fghz']
    good = np.where(np.logical_and(fghz != 0,fghz > 2.0))
    nant, npol, nf, nt = out['p'].shape
    # Interchange antenna and frequency axes in offsun and calfac
    osun = np.rollaxis(offsun,2,0)
    cfac = np.rollaxis(calfac,2,0)
    nant,npol,nf = cfac.shape
    for i in range(nt):
        out['p'][:nant,:,good,i] = (out['p'][:nant,:,good,i]-osun[:,:,good])*cfac[:,:,good]
#        out['tsys'][:,1,good,i] = (out['tsys'][:,1,good,i]-offsun[1,good,:])*calfac[1,good,:]
    return out
        
def sp_apply_cal(out,fghz,calfac,offsun):
    ''' Given "standard" output of tp_display.rd_tsys_multi(),
        and corresponding calibration data from sp_read_calfac(),
        apply the calibration and return the calibrated output.
    '''
    good = np.where(np.logical_and(fghz != 0,fghz > 2.0))
    nant, npol, nf, nt = out['tsys'].shape
    # Interchange antenna and frequency axes in offsun and calfac
    osun = np.rollaxis(offsun,2,0)
    cfac = np.rollaxis(calfac,2,0)
    for i in range(nt):
        out['tsys'][:,:,good,i] = (out['tsys'][:,:,good,i]-osun[:,:,good])*cfac[:,:,good]
#        out['tsys'][:,1,good,i] = (out['tsys'][:,1,good,i]-offsun[1,good,:])*calfac[1,good,:]
    return out
    
def sp_bg_subtract2(out,idx):
    ''' Given "standard" output of tp_display.rd_tsys_multi(),
        and a range of time indexes idx into the arrays, use the
        median of data at times indicated by idx as background.
        Subtract that background from both xtsys and ytsys
    '''
    nant, npol, nf, nt = out['p'].shape
    # Create the backgrounds (shape nf,nant)
    if len(idx) == 1:
        # Only one index in idx, so no median necessary
        bg = out['p'][:,:,:,idx]
    else:
        # Multiple values in idx, so do median
        bg = np.median(out['p'][:,:,:,idx],3)
    # Perform the background subtraction
    for i in range(nt):
        out['p'][:,:,:,i] -= bg
    return out
    
def sp_bg_subtract(out,idx):
    ''' Given "standard" output of tp_display.rd_tsys_multi(),
        and a range of time indexes idx into the arrays, use the
        median of data at times indicated by idx as background.
        Subtract that background from both xtsys and ytsys
    '''
    nant, npol, nf, nt = out['tsys'].shape
    # Create the backgrounds (shape nf,nant)
    if len(idx) == 1:
        # Only one index in idx, so no median necessary
        bg = out['tsys'][:,:,:,idx]
    else:
        # Multiple values in idx, so do median
        bg = np.median(out['tsys'][:,:,:,idx],3)
    # Perform the background subtraction
    for i in range(nt):
        out['tsys'][:,:,:,i] -= bg
    return out
    
def sp_bsize(x,y):
    ''' Make nice plots of the beamsize.
    '''
    import matplotlib.pyplot as plt

    nfrq, npnt, nant = x['rao'].shape
    nrow = 2
    ncol = (nant+1)/2
    f, ax = plt.subplots(nrow, ncol, sharex='col', sharey='row')
    f.set_size_inches(2*nant,7,forward=True)
    t1 = Time(x['ut_mjd'][0],format='mjd')
    f.suptitle('X-feed Beam Widths for SOLPNT scan at '+t1.iso[:19]+' UT',fontsize=18)
    fout,a,aout = disk_conv()
    fgood = np.where(x['fghz'] > 2.48)[0]
    for ant in range(nant):
        ax[ant % nrow, ant/nrow].set_title('Ant '+str(ant+1)+' [blue=RA, red=Dec]')
        if ant % ncol == 0:
            ax[ant/ncol,0].set_ylabel('Beam FWHM [deg]')
        if ant >= ncol:
            ax[1,ant - ncol].set_xlabel('Frequency [GHz]')
        ax[ant % nrow, ant/nrow].set_ylim(0.5,5)
        ax[ant % nrow, ant/nrow].set_xlim(1,20)
        ax[ant % nrow, ant/nrow].set_yscale('log')
        ax[ant % nrow, ant/nrow].set_xscale('log')
        ax[ant % nrow, ant/nrow].plot(x['fghz'][fgood],x['raparms'][2,fgood,ant]*2*np.sqrt(np.log(2))/10000,'b.')
        ax[ant % nrow, ant/nrow].plot(fout,aout,color='black')
        ax[ant % nrow, ant/nrow].plot(x['fghz'][fgood],x['decparms'][2,fgood,ant]*2*np.sqrt(np.log(2))/10000,'r.')
    plt.draw()
    f, ax = plt.subplots(nrow, ncol, sharex='col', sharey='row')
    f.set_size_inches(2*nant,7,forward=True)
    t1 = Time(x['ut_mjd'][0],format='mjd')
    f.suptitle('Y-feed Beam Widths for SOLPNT scan at '+t1.iso[:19]+' UT',fontsize=18)
    for ant in range(nant):
        ax[ant % nrow, ant/nrow].set_title('Ant '+str(ant+1)+' [blue=RA, red=Dec]')
        if ant % ncol == 0:
            ax[ant/ncol,0].set_ylabel('Beam FWHM [deg]')
        if ant >= ncol:
            ax[1,ant - ncol].set_xlabel('Frequency [GHz]')
        ax[ant % nrow, ant/nrow].set_ylim(0.5,5)
        ax[ant % nrow, ant/nrow].set_xlim(1,20)
        ax[ant % nrow, ant/nrow].set_yscale('log')
        ax[ant % nrow, ant/nrow].set_xscale('log')
        ax[ant % nrow, ant/nrow].plot(y['fghz'][fgood],y['raparms'][2,fgood,ant]*2*np.sqrt(np.log(2))/10000,'b.')
        ax[ant % nrow, ant/nrow].plot(fout,aout,color='black')
        ax[ant % nrow, ant/nrow].plot(y['fghz'][fgood],y['decparms'][2,fgood,ant]*2*np.sqrt(np.log(2))/10000,'r.')
    plt.draw()
    
def sp_offsets(x,y,save_plot=False):
    import matplotlib.pyplot as plt

    oldant = [8,9,10,12]
    nfrq, npnt, nant = x['rao'].shape
    nrow = (nant+1)/2
    ncol = 2
    f, ax = plt.subplots(nrow,ncol, sharex='col', sharey='row')
    f.set_size_inches(16,2*((nant+1)/2),forward=True)
    t1 = Time(x['ut_mjd'][0],format='mjd')
    f.suptitle('X-feed Offsets for SOLPNT scan at '+t1.iso[:19]+' UT',fontsize=18)
    user2rad = np.pi/10000./180.
    cosdec = np.cos(x['dec0'])
    fgood = np.where(x['fghz'] > 2.48)[0]
    xelx = []
    elx = []
    for ant in range(nant):
        ax[ant/ncol,ant % ncol].set_ylim(-0.5,0.5)
        ax[ant/ncol,ant % ncol].set_title('Ant '+str(ant+1)+' [blue=RA, red=Dec]')
        if ant >= nant-2:
            ax[nant/ncol-1,ant % ncol].set_xlabel('Frequency [GHz]')
        if ant % ncol == 0:
            ax[ant/ncol,ant % ncol].set_ylabel('Offset [deg]')
        ramed = np.median(x['raparms'][1,fgood,ant])
        ax[ant/ncol,ant % ncol].plot(x['fghz'][fgood],x['raparms'][1,fgood,ant]/10000,'b.')
        ax[ant/ncol,ant % ncol].plot([0,18],np.array([1,1])*ramed/10000,'b')
        decmed = np.median(x['decparms'][1,fgood,ant])
        ax[ant/ncol,ant % ncol].plot(x['fghz'][fgood],x['decparms'][1,fgood,ant]/10000,'r.')
        ax[ant/ncol,ant % ncol].plot([0,18],np.array([1,1])*decmed/10000,'r')
        if ant in oldant:
            ax[ant/ncol,ant % ncol].text(0.05,0.80,'RA = {:6.3f}'.format(ramed/10000.),color='red',transform=ax[ant/2,ant % 2].transAxes)
            ax[ant/ncol,ant % ncol].text(0.05,0.65,'Dec = {:6.3f}'.format(decmed/10000.),color='red',transform=ax[ant/2,ant % 2].transAxes)
        else:
            ax[ant/ncol,ant % ncol].text(0.05,0.80,'RA = {:6.3f}'.format(ramed/10000.),color='gray',transform=ax[ant/2,ant % 2].transAxes)
            ax[ant/ncol,ant % ncol].text(0.05,0.65,'Dec = {:6.3f}'.format(decmed/10000.),color='gray',transform=ax[ant/2,ant % 2].transAxes)
        xel, el = solpnt.dradec2dazel(x['ra0'],x['dec0'],t1,ramed*user2rad/cosdec,decmed*user2rad)
        if ant in oldant:
            ax[ant/ncol,ant % ncol].text(0.65,0.80,'XEL = {:6.3f}'.format(xel*180./np.pi),color='gray',transform=ax[ant/2,ant % 2].transAxes)
            ax[ant/ncol,ant % ncol].text(0.65,0.65,'EL = {:6.3f}'.format(el*180./np.pi),color='gray',transform=ax[ant/2,ant % 2].transAxes)
        else:
            ax[ant/ncol,ant % ncol].text(0.65,0.80,'XEL = {:6.3f}'.format(xel*180./np.pi),color='red',transform=ax[ant/2,ant % 2].transAxes)
            ax[ant/ncol,ant % ncol].text(0.65,0.65,'EL = {:6.3f}'.format(el*180./np.pi),color='red',transform=ax[ant/2,ant % 2].transAxes)
        if ant in oldant:
            xelx.append(ramed/10000)
            elx.append(decmed/10000)
        else:
            xelx.append(xel*180./np.pi)
            elx.append(el*180./np.pi)
    plt.draw()
    f, ax = plt.subplots(nrow, ncol, sharex='col', sharey='row')
    f.set_size_inches(16,2*((nant+2)/2),forward=True)
    f.suptitle('Y-feed Offsets for SOLPNT scan at '+t1.iso[:19]+' UT',fontsize=18)
    user2rad = np.pi/10000./180.
    xely = []
    ely = []
    for ant in range(nant):
        ax[ant/ncol,ant % ncol].set_ylim(-0.5,0.5)
        ax[ant/ncol,ant % ncol].set_title('Ant '+str(ant+1)+' [blue=RA, red=Dec]')
        if ant >= nant-2:
            ax[nant/ncol-1,ant % ncol].set_xlabel('Frequency [GHz]')
        if ant % ncol == 0:
            ax[ant/ncol,ant % ncol].set_ylabel('Offset [deg]')
        ramed = np.median(y['raparms'][1,fgood,ant])
        ax[ant/ncol,ant % ncol].plot(y['fghz'][fgood],y['raparms'][1,fgood,ant]/10000,'b.')
        ax[ant/ncol,ant % ncol].plot([0,18],np.array([1,1])*ramed/10000,'b')
        decmed = np.median(y['decparms'][1,fgood,ant])
        ax[ant/ncol,ant % ncol].plot(y['fghz'][fgood],y['decparms'][1,fgood,ant]/10000,'r.')
        ax[ant/ncol,ant % ncol].plot([0,18],np.array([1,1])*decmed/10000,'r')
        if ant in oldant:
            ax[ant/ncol,ant % ncol].text(0.05,0.80,'RA = {:6.3f}'.format(ramed/10000.),color='red',transform=ax[ant/2,ant % 2].transAxes)
            ax[ant/ncol,ant % ncol].text(0.05,0.65,'Dec = {:6.3f}'.format(decmed/10000.),color='red',transform=ax[ant/2,ant % 2].transAxes)
        else:
            ax[ant/ncol,ant % ncol].text(0.05,0.80,'RA = {:6.3f}'.format(ramed/10000.),color='gray',transform=ax[ant/2,ant % 2].transAxes)
            ax[ant/ncol,ant % ncol].text(0.05,0.65,'Dec = {:6.3f}'.format(decmed/10000.),color='gray',transform=ax[ant/2,ant % 2].transAxes)
        xel, el = solpnt.dradec2dazel(y['ra0'],y['dec0'],t1,ramed*user2rad/cosdec,decmed*user2rad)
        if ant in oldant:
            ax[ant/ncol,ant % ncol].text(0.65,0.80,'XEL = {:6.3f}'.format(xel*180./np.pi),color='gray',transform=ax[ant/2,ant % 2].transAxes)
            ax[ant/ncol,ant % ncol].text(0.65,0.65,'EL = {:6.3f}'.format(el*180./np.pi),color='gray',transform=ax[ant/2,ant % 2].transAxes)
        else:
            ax[ant/ncol,ant % ncol].text(0.65,0.80,'XEL = {:6.3f}'.format(xel*180./np.pi),color='red',transform=ax[ant/2,ant % 2].transAxes)
            ax[ant/ncol,ant % ncol].text(0.65,0.65,'EL = {:6.3f}'.format(el*180./np.pi),color='red',transform=ax[ant/2,ant % 2].transAxes)
        if ant in oldant:
            xely.append(ramed/10000)
            ely.append(decmed/10000)
        else:
            xely.append(xel*180./np.pi)
            ely.append(el*180./np.pi)
    plt.draw()
    xout = (np.array(xelx) + np.array(xely))/2.
    yout = (np.array(elx) + np.array(ely))/2.
    dxout = abs((np.array(xelx) - np.array(xely)))/2.
    dyout = abs((np.array(elx) - np.array(ely)))/2.
    plt.figure()
    plt.errorbar(xout,yout,xerr=dxout,yerr=dyout,fmt=' ')
    plt.axis('equal')
    plt.plot(0.25*np.cos(np.linspace(0,2*np.pi)),0.25*np.sin(np.linspace(0,2*np.pi)))
    plt.xlim(-0.5,0.5)
    plt.ylim(-0.5,0.5)
    plt.xlabel('XDec or XEl Offset [deg]')
    plt.xlabel('Dec or El Offset [deg]')
    for i in range(nant):
        plt.text(xout[i],yout[i],str(i+1))
    plt.title('EOVSA Pointing for '+t1.iso)
    if save_plot:
        tstr = t1.iso.replace('-','').replace(':','').replace(' ','')[:14]
        plt.savefig('/common/webplots/PTG/PTG'+tstr+'.png',bbox_inches='tight')
    return xout,yout,dxout,dyout

def offsets2ants(t,xoff,yoff,ant_str=None):
    ''' Given a start time (Time object) and a list of offsets output by sp_offsets()
        for 13 antennas, convert to pointing coefficients (multiply by 10000),
        add to coefficients listed in stateframe, and send to the relevant 
        antennas.  The antennas to update are specified with ant_str 
        (defaults to no antennas, for safety).        
    '''
    if ant_str is None:
        print 'No antenna list specified, so there is nothing to do!'
        return
    try:
        timestamp = int(Time(t,format='mjd').lv)
    except:
        print 'Error interpreting time as Time() object'
        return
    import pcapture2 as p
    import dbutil as db
    import adc_cal2 as ac
    import stateframe as stf
    accini = stf.rd_ACCfile()
    acc = {'host': accini['host'], 'scdport':accini['scdport']}
    antlist = p.ant_str2list(ant_str)
    if antlist is None:
        return
    cursor = db.get_cursor()
    # Read current stateframe data (as of 10 s ago)
    D15data = db.get_dbrecs(cursor,dimension=15,timestamp=timestamp,nrecs=1)
    p1_cur, = D15data['Ante_Cont_PointingCoefficient1']
    p7_cur, = D15data['Ante_Cont_PointingCoefficient7']
    
    for i in antlist:
        p1_inc = int(xoff[i]*10000)
        p7_inc = int(yoff[i]*10000)
        p1_new = p1_cur[i] - p1_inc
        p7_new = p7_cur[i] + p7_inc
        print 'Updating P1 for Ant',i+1,'P1_old =',p1_cur[i],'P1_inc =',p1_inc,'P1_new =',p1_new
        cmd1 = 'pointingcoefficient1 '+str(p1_new)+' ant'+str(i+1)
        print 'Updating P7 for Ant',i+1,'P7_old =',p7_cur[i],'P7_inc =',p7_inc,'P7_new =',p7_new
        cmd7 = 'pointingcoefficient7 '+str(p7_new)+' ant'+str(i+1)
        print 'Commands to be sent:'
        print cmd1
        print cmd7
        ac.send_cmds([cmd1],acc)
        ac.send_cmds([cmd7],acc)

if __name__ == "__main__":
    ''' Run automatically via cron job, or at command line.
        Usage: python /common/python/current/calibration.py "2014-12-15 18:30"
        where the time string is optional.  If omitted, the current time is used.
        
        The logic is to check whether there is a new SOLPNTCAL available, and
        analyze it if so.  Compares times from solpnt.find_solpnt() with current
        time and analyzes any that are between 5 and 10 minutes old.
    '''
    arglist = str(sys.argv)
    t = Time.now()
    if len(sys.argv) == 2:
        try:
            t = Time(sys.argv[1])
        except:
            print 'Cannot interpret',sys.argv[1],'as a valid date/time string.'
            exit()

    timestamp = t.lv  # Current timestamp
    times, tstamp = solpnt.find_solpnt(t)
    # Find first SOLPNTCAL occurring after timestamp (time given by Time() object)
    if times == []:
        # No SOLPNTCAL scans (yet)
        print t.iso[:19]+': No SOLPNTCAL scans for today'
        exit()
    elif type(times[0]) is np.ndarray:
        # Annoyingly necessary when only one time in tstamps
        times = times[0]
    igt5 = np.where((timestamp - times) > 300)[0]
    if len(igt5) == 0:
        # SOLPNTCAL scan in progress
        print t.iso[:19]+': SOLPNTCAL scan still in progress'
        exit()
    if (timestamp - times[igt5[-1]]) > 600:
        # Latest SOLPNTCAL scan is too old, so nothing to do
        print t.iso[:19]+': Last SOLPNTCAL scan too old. Age:',int(timestamp - times[igt5[-1]])/60,'minutes'
        exit()
    # Looks like this is the right "age" and ready to be analyzed
    # Wait 30 s to ensure scan is done.
    time.sleep(30)
    t = Time(times[igt5[-1]],format='lv')
    if socket.gethostname() == 'pipeline':
        x, y, qual = solpntanal(t,udb=True)
    elif socket.gethostname() == 'dpp':
        x, y, qual = solpntanal(t,udb=False)
        xout,yout,dxout,dyout = sp_offsets(x,y,save_plot=True)
    else:
        print 'CALIBRATION Error: This routine only runs on dpp or pipeline.'
        exit()
    status = np.zeros(qual.shape,'S4')
    status[np.where(qual)] = 'Good'
    status[np.where(qual==False)] = '*Bad'
    nparms,nf,nant = x['raparms'].shape
    print t.iso[:19],': Quality of TP Calibration'
    print '    Ant      X-Feed        Y-Feed'
    print '            RA    DEC     RA    DEC'
    for i in range(nant):
        print '    ',i+1,' ',status[:,i]
    percent_good = len(np.where(qual)[0])*100/(4*nant)
    if percent_good > 50:
        calfac, offsun = sp_get_calfac(x,y)
        # If another file for today's date already exists, this will overwrite it
        sp_write_calfac(x,y,calfac,offsun)
        print 'Calibration file successfully written'
        if t.iso[:10] == Time.now().iso[:10]:
           # The calibration file is for today, so rewrite as tomorrow's file also
           x['ut_mjd'][0] = x['ut_mjd'][0]+1.   # Add one day to first timestamp
           sp_write_calfac(x,y,calfac,offsun)
           print "Also wrote tomorrow's file"
    else:
        print 'Calibration file not written--too many bad values.'
    exit()  
