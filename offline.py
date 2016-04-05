'''
   Module for accessing, analyzing and plotting data offlikne'''
#
# History:
#   2014-Dec-19  DG
#     First written.  Creating calibration requires working online
#     (needs access to SQL database/stateframe), so really there is
#     not much to be done unless saved calibration data are available.
#   2014-Dec-20  DG
#     Cleaned up some labeling of tsys_show_dynspec().  Also added
#     idx argument for restricting the timerange to be shown.
#   2015-Apr-10  DG
#     Fixed a problem introduced earlier where I converted out['ut'] to
#     integer, which killed the date plotting.  I covert it back to float
#     now in tsys_show_dynspec()
#   2015-Apr-18  DG
#     In tsys_show_dynspec(), subtract day1 instead of using % 1 to allow
#     times to span a day.  This is still not entirely satisfactory,
#     and a considerable rewrite of date/time code is necessary
#   2015-May-29  DG
#     Converted from using datime() to using Time() based on astropy.
#   2015-May-31  DG
#     Fixed rd_tsys_multi() to read files correctly when trange spans a
#     day change.  (Still does not work for more than two consecutive days.)
#   2015-Jun-27  DG
#     Changed the code to standardize on names, content, and order of indices
#     of outputs for rd_miriad_tsys.  Names ut_mjd, fghz, and tsys will be used, 
#     with units hopefully made obvious.  Order of indices will be 
#          (nant/nbl, npol, nfreq, ntimes).
#     This required some changes to several of the routines here.
#

import numpy as np
import solpnt
import struct, time, glob, sys, socket
from disk_conv import *


def read_calfac(t):
    ''' Read the contents of a SOLPNT calibration file and return
        fghz, calfac, offsun arrays
    '''
    #files = glob.glob('/data1/TPCAL/*')
    files = glob.glob('/common/tmp/TPCAL*')
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

def apply_cal(out,fghz,calfac,offsun):
    ''' Given "standard" output of tp_display.rd_tsys_multi(),
        and corresponding calibration data from sp_read_calfac(),
        apply the calibration and return the calibrated output.
    '''
    nant, npol, nf, nt = out['tsys'].shape
    # Find the common frequencies (to within 1 MHz) in the calibration and in the data
    idx1, idx2 = solpnt.common_val_idx((fghz*1000).astype('int'),(out['fghz']*1000).astype('int'))
    # Interchange antenna and frequency axes in offsun and calfac
    osun = np.rollaxis(offsun,2,0)
    cfac = np.rollaxis(calfac,2,0)
    for i in range(nt):
        out['tsys'][:,:,idx2,i] = (out['tsys'][:,:,idx2,i]-osun[:,:,idx1])*cfac[:,:,idx1]
    # Restrict output to only calibrated frequencies
    out['tsys'] = out['tsys'][:,:,idx2,:]
    out['fghz'] = out['fghz'][idx2]
    return out
    
def bg_subtract(out,idx):
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

#
# tp_display routines
#
import matplotlib.pyplot as plt

def qlook(file):
    ''' Routine to show all four channels (2 ants, 2 poln) for
        a pcapture file with a name of pattern 'bandnn_e.pcap'
        where nn is band number (can be 1 digit) and e is 
        ethernet interface (eth) number (2 or 3)
    '''
    out = p.rd_spec(file).reshape(200,4096,8)
    f, ax = plt.subplots(4,1)
    ax[0].imshow(out[:,:,0])
    ax[1].imshow(out[:,:,1])
    ax[2].imshow(out[:,:,2])
    ax[3].imshow(out[:,:,3])
    band = file.split('band')[1].split('_')[0]
    if file.find('_2') != -1: 
        ax[0].set_title('Ants 1 & 2, band '+band)
    else:
        ax[0].set_title('Ants 3 & 4, band '+band)

  
def rd_tsys_multi(trange):
 
    ''' Given 2-element Time() object specifying a timerange, dump and
        return the tsys data for all IDB files in that timerange.
    '''
    import glob,sys, os
    if sys.platform[:3] == 'win':
        path = 'c:\\users\\gary\\dropbox\\PythonCode\\current\\tmp\\txt'
    else:
        path = '/common/tmp/txt'
        
    def fname2mjd(filename):
        from util import Time
        fstem = filename.split(os.sep)[-1]
        fstr = fstem[2:6]+'-'+fstem[6:8]+'-'+fstem[8:10]+' '+fstem[10:12]+':'+fstem[12:14]+':'+fstem[14:16]
        print fstr
        return Time(fstr).mjd

    datstr = trange[0].iso.replace('-','')[:8]+'*.txt'
    xfiles = glob.glob(path+os.sep+'xt'+datstr)
    if int(trange[0].mjd) != int(trange[1].mjd):
        datstr2 = trange[1].iso.replace('-','')[:8]+'*.txt'
        xfiles += glob.glob(path+os.sep+'xt'+datstr2)
    xfiles.sort()
    xtsyslist = []  # List of xtsys arrays
    utlist = []     # List of times
    ytsyslist = []  # List of ytsys arrays
    for file in xfiles:
        mjd = fname2mjd(file)
        if mjd >= trange[0].mjd and mjd < trange[1].mjd:
            otp = solpnt.rd_tsys(file,path+os.sep+'sfreq.txt')
            xtsyslist.append(otp['tsys'])
            utlist.append(otp['ut_mjd'])
            otp = solpnt.rd_tsys(file.replace('/xt','/yt'),path+os.sep+'sfreq.txt')
            ytsyslist.append(otp['tsys'])
    # We have read all of the files, now concatenate into two dicts with single arrays
    xtsys = np.concatenate(xtsyslist,2)
    ytsys = np.concatenate(ytsyslist,2)
    ut = np.concatenate(utlist)
    tsys = np.rollaxis(np.array((xtsys,ytsys)),1,0)
    out = {'fghz':otp['fghz'],'tsys':tsys,'ut_mjd':ut}
    return out
    
def tsys_show_dynspec(out,idx=None,ampscl=None,domedian=True,frq='linear'):
    ''' Given "standard" output of rd_tsys_multi(), possibly
        calibrated using calibration.sp_apply_cal() and/or 
        background-subtracted using calibration.sp_bg_subtract(),
        make a nice image plot of the dynamic spectrum.  The
        plot can contain multiple panels if domedian is False,
        or plot a single spectrum representing the median of
        multiple antennas.  Only linear frequency scale is supported
        at this time.
    '''
    from matplotlib.image import NonUniformImage
    from matplotlib.dates import AutoDateLocator, DateFormatter
    from matplotlib import colors
    from matplotlib.pyplot import locator_params
    from util import Time

    nant, npol, nf, nt = out['tsys'].shape
    if idx is None:
        idx_ = np.arange(nt)
    else:
        idx_ = idx
    ut = Time(out['ut_mjd'][idx_],format='mjd')
    utd = (ut.mjd - int(ut[0].mjd) + 1).astype('float')
    good = np.where(out['fghz'] > 2.0)[0]
    if frq == 'linear':
        fghz = out['fghz'][good]
    else:
        fghz = np.log10(out['fghz'][good])

    locator = AutoDateLocator(interval_multiples=True) #maxticks=7,minticks=3,
    tsys = out['tsys'][:,:,good,idx_]
    if domedian:
        medtsys = np.nanmedian(np.nanmedian(tsys,0),0)
        fig = plt.figure()
        fig.suptitle('EOVSA Total Power Data for '+ut[0].iso[:10],fontsize=14)
        # Plot X-feed
        ax = fig.add_subplot(211)
        ax.xaxis_date()
        ax.set_ylabel('Frequency [GHz]')
        ax.set_title('Median Total Power')
        extent=[utd[0],utd[-1],fghz[0],fghz[-1]]
#        extent=[ut[0],ut[-1],fghz[0],fghz[-1]]
        im = NonUniformImage(ax,extent=extent)
        if ampscl != None:
            im.set_norm(colors.Normalize(vmin=ampscl[0],vmax=ampscl[1]))
        im.set_data(utd,fghz,medtsys)
        #fig.colorbar(im,ax)
        ax.images.append(im)
        ax.set_xlim(extent[0],extent[1])
        ax.set_ylim(extent[2],extent[3])
        ax.xaxis.set_major_locator(locator)
        #ax.xaxis.set_minor_locator(MinuteLocator(interval=10))
        # Set up date formatting
        fmt = DateFormatter('%H:%M:%S')
        ax.xaxis.set_major_formatter(fmt)
        labels = (10**ax.get_yticks()).astype('str')
        for i in range(len(labels)):
            labels[i] = labels[i][:4]
        ax.set_yticklabels(labels)
        # Repeat for Y-feed
        ax = fig.add_subplot(212)
        ax.xaxis_date()
        ax.set_xlabel('Start time '+ut[0].iso[:19]+' UT')
        ax.set_title('Median of Y-poln')
        ax.set_ylabel('Frequency [GHz]')
        im = NonUniformImage(ax,extent=extent)
        if ampscl != None:
            im.set_norm(colors.Normalize(vmin=ampscl[0],vmax=ampscl[1]))
        im.set_data(utd,fghz,medytsys)
        #fig.colorbar(im,ax)
        ax.images.append(im)
        ax.set_xlim(extent[0],extent[1])
        ax.set_ylim(extent[2],extent[3])
        ax.xaxis.set_major_locator(locator)
        # Set up date formatting
        ax.xaxis.set_major_formatter(fmt)
        ax.set_yticklabels(labels)
    plt.draw()
    plt.show()
