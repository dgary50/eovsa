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
    otp = dump_tsys.rd_miriad_tsys(trange,udb=udb)
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
        f, ax = plt.subplots(2, nant/2, sharex='col', sharey='row')
        f.set_size_inches(2*nant,7,forward=True)
        f.suptitle('Calibration for SOLPNT scan at '+t.iso[:19]+' UT',fontsize=18)
        for ant in range(nant):
            ax[ant % 2, ant/2].set_title('Ant '+str(ant+1)+' Solar Spectrum')
            if ant % 4 == 0:
                ax[ant/4,0].set_ylabel('Solar Flux [sfu]')
            if ant >= nant/2:
                ax[1,ant - nant/2].set_xlabel('Frequency [GHz]')

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
            ax[ant % 2, ant/2].plot(fghz,s1*calfac[1,:,ant],'.',label='Y-RA')
            ax[ant % 2, ant/2].plot(fghz,s2*calfac[1,:,ant],'.',label='Y-Dec')
            ax[ant % 2, ant/2].set_ylim(0,600)
            ax[ant % 2, ant/2].set_xlim(0,19)
            ax[ant % 2, ant/2].legend(loc='upper left',fontsize='small')
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
    filename = '/common/tmp/TPCAL'+datstr+'_'+str(dims[0])+'_'+str(dims[1])+'_'+str(dims[2])+'.dat'
    #filename = '/data1/TPCAL/TPCAL'+datstr+'_'+str(dims[0])+'_'+str(dims[1])+'_'+str(dims[2])+'.dat'
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
    ''' Make nice plots of the beamsize, currently only for the first eight
        antennas.  This routine will have to be modified to work with more.
    '''
    import matplotlib.pyplot as plt

    nfrq, npnt, nant = x['rao'].shape
    f, ax = plt.subplots(2, nant/2, sharex='col', sharey='row')
    f.set_size_inches(2*nant,7,forward=True)
    t1 = Time(x['ut_mjd'][0],format='mjd')
    f.suptitle('X-feed Beam Widths for SOLPNT scan at '+t1.iso[:19]+' UT',fontsize=18)
    fout,a,aout = disk_conv()
    fgood = np.where(x['fghz'] > 2.48)[0]
    for ant in range(nant):
        ax[ant % 2, ant/2].set_title('Ant '+str(ant+1)+' [blue=RA, red=Dec]')
        if ant % 4 == 0:
            ax[ant/4,0].set_ylabel('Beam FWHM [deg]')
        if ant >= nant/2:
            ax[1,ant - nant/2].set_xlabel('Frequency [GHz]')
        ax[ant % 2, ant/2].set_ylim(0.5,5)
        ax[ant % 2, ant/2].set_xlim(1,20)
        ax[ant % 2, ant/2].set_yscale('log')
        ax[ant % 2, ant/2].set_xscale('log')
        ax[ant % 2, ant/2].plot(x['fghz'][fgood],x['raparms'][2,fgood,ant]*2*np.sqrt(np.log(2))/10000,'b.')
        ax[ant % 2, ant/2].plot(fout,aout,color='black')
        ax[ant % 2, ant/2].plot(x['fghz'][fgood],x['decparms'][2,fgood,ant]*2*np.sqrt(np.log(2))/10000,'r.')
    plt.draw()
    f, ax = plt.subplots(2, nant/2, sharex='col', sharey='row')
    f.set_size_inches(2*nant,7,forward=True)
    t1 = Time(x['ut_mjd'][0],format='mjd')
    f.suptitle('Y-feed Beam Widths for SOLPNT scan at '+t1.iso[:19]+' UT',fontsize=18)
    for ant in range(nant):
        ax[ant % 2, ant/2].set_title('Ant '+str(ant+1)+' [blue=RA, red=Dec]')
        if ant % 4 == 0:
            ax[ant/4,0].set_ylabel('Beam FWHM [deg]')
        if ant >= nant/2:
            ax[1,ant - nant/2].set_xlabel('Frequency [GHz]')
        ax[ant % 2, ant/2].set_ylim(0.5,5)
        ax[ant % 2, ant/2].set_xlim(1,20)
        ax[ant % 2, ant/2].set_yscale('log')
        ax[ant % 2, ant/2].set_xscale('log')
        ax[ant % 2, ant/2].plot(y['fghz'][fgood],y['raparms'][2,fgood,ant]*2*np.sqrt(np.log(2))/10000,'b.')
        ax[ant % 2, ant/2].plot(fout,aout,color='black')
        ax[ant % 2, ant/2].plot(y['fghz'][fgood],y['decparms'][2,fgood,ant]*2*np.sqrt(np.log(2))/10000,'r.')
    plt.draw()
    
def sp_offsets(x,y):
    import matplotlib.pyplot as plt

    nfrq, npnt, nant = x['rao'].shape
    f, ax = plt.subplots(nant/2, 2, sharex='col', sharey='row')
    f.set_size_inches(16,2*((nant+1)/2),forward=True)
    t1 = Time(x['ut_mjd'][0],format='mjd')
    f.suptitle('X-feed Offsets for SOLPNT scan at '+t1.iso[:19]+' UT',fontsize=18)
    user2rad = np.pi/10000./180.
    cosdec = np.cos(x['dec0'])
    fgood = np.where(x['fghz'] > 2.48)[0]
    for ant in range(nant):
        ax[ant/2,ant % 2].set_ylim(-0.5,0.5)
        ax[ant/2,ant % 2].set_title('Ant '+str(ant+1)+' [blue=RA, red=Dec]')
        if ant >= nant-2:
            ax[nant/2-1,ant % 2].set_xlabel('Frequency [GHz]')
        if ant % 2 == 0:
            ax[ant/2,ant % 2].set_ylabel('Offset [deg]')
        ramed = np.median(x['raparms'][1,fgood,ant])
        ax[ant/2,ant % 2].plot(x['fghz'][fgood],x['raparms'][1,fgood,ant]/10000,'b.')
        ax[ant/2,ant % 2].plot([0,18],np.array([1,1])*ramed/10000,'b')
        ax[ant/2,ant % 2].text(0.05,0.80,'RA = {:6.3f}'.format(ramed/10000.),transform=ax[ant/2,ant % 2].transAxes)
        decmed = np.median(x['decparms'][1,fgood,ant])
        ax[ant/2,ant % 2].plot(x['fghz'][fgood],x['decparms'][1,fgood,ant]/10000,'r.')
        ax[ant/2,ant % 2].plot([0,18],np.array([1,1])*decmed/10000,'r')
        ax[ant/2,ant % 2].text(0.05,0.65,'Dec = {:6.3f}'.format(decmed/10000.),transform=ax[ant/2,ant % 2].transAxes)
        xel, el = solpnt.dradec2dazel(x['ra0'],x['dec0'],t1,ramed*user2rad/cosdec,decmed*user2rad)
        ax[ant/2,ant % 2].text(0.65,0.80,'XEL = {:6.3f}'.format(xel*180./np.pi),transform=ax[ant/2,ant % 2].transAxes)
        ax[ant/2,ant % 2].text(0.65,0.65,'EL = {:6.3f}'.format(el*180./np.pi),transform=ax[ant/2,ant % 2].transAxes)
    plt.draw()
    f, ax = plt.subplots(nant/2, 2, sharex='col', sharey='row')
    f.set_size_inches(16,2*((nant+2)/2),forward=True)
    f.suptitle('Y-feed Offsets for SOLPNT scan at '+t1.iso[:19]+' UT',fontsize=18)
    user2rad = np.pi/10000./180.
    for ant in range(nant):
        ax[ant/2,ant % 2].set_ylim(-0.5,0.5)
        ax[ant/2,ant % 2].set_title('Ant '+str(ant+1)+' [blue=RA, red=Dec]')
        if ant >= nant-2:
            ax[nant/2-1,ant % 2].set_xlabel('Frequency [GHz]')
        if ant % 2 == 0:
            ax[ant/2,ant % 2].set_ylabel('Offset [deg]')
        ramed = np.median(y['raparms'][1,fgood,ant])
        ax[ant/2,ant % 2].plot(y['fghz'][fgood],y['raparms'][1,fgood,ant]/10000,'b.')
        ax[ant/2,ant % 2].plot([0,18],np.array([1,1])*ramed/10000,'b')
        ax[ant/2,ant % 2].text(0.05,0.80,'RA = {:6.3f}'.format(ramed/10000.),transform=ax[ant/2,ant % 2].transAxes)
        decmed = np.median(y['decparms'][1,fgood,ant])
        ax[ant/2,ant % 2].plot(y['fghz'][fgood],y['decparms'][1,fgood,ant]/10000,'r.')
        ax[ant/2,ant % 2].plot([0,18],np.array([1,1])*decmed/10000,'r')
        ax[ant/2,ant % 2].text(0.05,0.65,'Dec = {:6.3f}'.format(decmed/10000.),transform=ax[ant/2,ant % 2].transAxes)
        xel, el = solpnt.dradec2dazel(y['ra0'],y['dec0'],t1,ramed*user2rad/cosdec,decmed*user2rad)
        ax[ant/2,ant % 2].text(0.65,0.80,'XEL = {:6.3f}'.format(xel*180./np.pi),transform=ax[ant/2,ant % 2].transAxes)
        ax[ant/2,ant % 2].text(0.65,0.65,'EL = {:6.3f}'.format(el*180./np.pi),transform=ax[ant/2,ant % 2].transAxes)
    plt.draw()
    
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
