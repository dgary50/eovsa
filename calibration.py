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
#   2017-Jan-29  DG
#    Updated calpntanal() to take a single time and find the relevant scan
#    to analyze.  Also added calpnt_multi(), which calls calpntanal() for all
#    scans in a timerange.
#   2017-Feb-05  DG
#    Made changes to calpnt_multi() and calpntanal() to improve plots.
#   2017-Feb-09  DG
#    Still struggling with getting the right signs on updating offsets on antennas.
#    I have verified that the old antennas 9, 10, 11 and 12 need to have the
#    opposite sign for P1 correction (RA to HA involves inverting the sign).  This
#    sign reversal is implemented in offsets2ants()
#   2017-Apr-11  DG
#    Some changes to calpnt_multi() to survive a bad scan, and write results
#    to a file as it goes along.  This should allow restarting where we left
#    off in case of a crash or looking at results before the entire set of
#    observations is done.
#   2017-Apr-28  DG
#    Began work on calpntanal() to perform on either CALPNTCAL or CALPNT2M data
#    with automatic adjustment based on scans in fdb file.
#   2017-May-17  DG
#    Added 'auto' keyword to solpntanal(), for optionally analyzing solar
#    pointing based on auto-correlation instead of total power.
#   2017-Jul-30  DG
#    Added 'auto' keyword to sp_apply_cal2() and sp_bg_subtract2(), to apply
#    the calibration to auto-correlation if desired.
#   2017-Oct-27  DG
#    Added tsql keyword to solpnt2sql(), for overriding date of SQL record.
#   2017-Dec-03  DG
#    Now tsql is set by default to 20 UT on the date of observation, which is
#    the new standard.  A second analysis at any time of day will overwrite the
#    existing one.
#   2018-Jan-04  DG
#    Added best_solpnt2sql(), to aid in getting past TP calibrations into SQL.
#    Also added a find parameter to solpntanal() call.
#   2018-Jan-07  DG
#    Small bug fixes to reduce crashes for unusual situations.
#   2018-Jan-11  DG
#    imported ant_str2list() from util, instead of trying to use pcapture2 version
#   2018-Jan-15  DG
#    Changed calpntanal() to require feed-rotation corrected data (it was not
#    working well otherwise).  As part of this, an alternate directory can be
#    provided in either calpntanal() or calpnt_multi().
#   2018-Jun-14  DG
#    Added skycal_anal() routine, which together with some code in solpntanal() will
#    be used to allow fitting of lower frequencies, for which the current SOLPNTCAL
#    uses too small a range of distances for fitting the large primary beam.  More
#    to come.
#   2019-May-17  DG
#    Changes to the case of running from the command line, since the routine runs
#    fine from pipeline now.
#   2019-May-21  DG
#    Finally completed the skycal_anal() routine, which is now called in
#    solpntanal() to improve fitting by adding pseudo-10-degree pointing offset.
#    Also added code to give a "poor-man's" skycal when those are not available,
#    by using the GAINCAL full attenuation setting.
#   2019-Jun-23  DG
#    Change to calpntanal() to do feed rotation correction, and change to calpnt_multi() 
#    to append '2m' to output filename if calpnt2m = True and outfile is not specified, 
#    to avoid mixing 2-m and 27-m results in the same file.  Also fixed a bug in
#    calpnt_multi() that occurred in ax array when nrows is 1.
#  2020-May-10  DG
#    Updated skycal_anal() to use util.get_idbdir() to find IDB root path.
#  2020-May-11  DG
#    Further update to make this work on the DPP.
#  2020-May-29 SY
#    update calpntanal() to use util.get_idbdir() to find IDB root path.
#  2020-Dec-12 DG
#    Added sp_explore() routine for investigating the data and fits for
#    a given antenna and frequency.
#  2021-Jan-11 DG
#    Important change to default to correcting for correlator saturation
#    (only affects the total power auto-correlation calibration).  This
#    activates the correction in the read_idb.py function readXdata().
#  2021-Jul-25  DG
#    I discovered a bug in skycal_anal() where the times of SQL levels and data 
#    power levels could be mismatched.  Now finds the common value indexes first.
#    Argh.  Use of median() instead of nanmedian() got me again.
#  2022-Jun-28  DG
#    Found that offsets2ants() sometimes read zeros from SQL.  It now bails
#    with a message to try a different timestamp.
#

if __name__ == "__main__":
    ''' For cron-job use, make sure we use a backend that does not require
        X-server    
    '''
    import matplotlib
    matplotlib.use('Agg')

import numpy as np
import solpnt
from util import Time, nearest_val_idx, lobe, ant_str2list, get_idbdir, common_val_idx
import struct, time, glob, sys, socket, os
from disk_conv import *
import dump_tsys

def calpnt_multi(trange, fdir=None, ant_str='ant1-13', calpnt2m=False, do_plot=True, outfile=None):
    ''' Runs calpntanal() for all scans within the given time range.
    
        trange   Time object with start and end time over which scans
                   should be identified and analyzed.
                   
        Returns a list of dictionaries containing 
           Source name, Time, HA, Dec, RA offset and Dec offset
    '''
    fdb = dump_tsys.rd_fdb(trange[0])
    #calpnt2m = False   # Flag to indicate whether we are doing CALPNT2M analysis
    if calpnt2m:
        # Search for CALPNT2M
        scanidx, = np.where(fdb['PROJECTID'] == 'CALPNT2M')
        if scanidx.size == 0:
            print 'No CALPNT2M project IDs found for date'+t.iso[:10]
            return {}
    else:
        # Search for CALPNTCAL
        scanidx, = np.where(fdb['PROJECTID'] == 'CALPNTCAL')
        if scanidx.size == 0:
            print 'No CALPNTCAL project IDs found for date'+t.iso[:10]
            return {}
    scans,sidx = np.unique(fdb['SCANID'][scanidx],return_index=True)
    tlist = Time(fdb['ST_TS'][scanidx[sidx]].astype(float).astype(int),format='lv')
    telist = Time(fdb['EN_TS'][scanidx[sidx]].astype(float).astype(int),format='lv')
    idx = ant_str2list(ant_str)
    out = []
    if outfile is None:
        if calpnt2m:
            outfile = '/common/tmp/Pointing/'+tlist[0].iso[:10].replace('-','')+'_calpnt2m.txt'
        else:
            outfile = '/common/tmp/Pointing/'+tlist[0].iso[:10].replace('-','')+'_calpnt27m.txt'
    k = -1
    # First find out how many pointings there will be
    npt = 0
    for i,t in enumerate(tlist):
        if t.jd >= trange[0].jd and telist[i].jd <= trange[1].jd:
            npt += 1

    for i,t in enumerate(tlist):
        if t.jd >= trange[0].jd and telist[i].jd <= trange[1].jd:
            k += 1  # Plot-axis row pointer
            if do_plot:
                if k % 6 == 0:
                    nrows = min(npt-k,6)  # Number of rows in this plot--up to 6
                    if calpnt2m:
                        f, ax = plt.subplots(nrows,2*idx.size)
                        f.set_size_inches(2.5*idx.size,1.5*nrows,forward=True)
                        plt.subplots_adjust(left=0.03, right=0.97, top=0.89, bottom=0.3, wspace=0.1, hspace=0.3)
                    else:
                        f, ax = plt.subplots(nrows,2)
                        f.set_size_inches(2.5,1.5*nrows,forward=True)
                    k = 0  # Reset row pointer to top of new plot
                else:
                    # Kill previous row's xaxis tick labels so that title is visible
                    if calpnt2m:
                        for j in range(idx.size*2):
                            ax[k-1,j].xaxis.set_ticklabels([])
                            ax[k-1,j].set_xlabel('')
                    else:
                        for j in range(2):
                            ax[k-1,j].xaxis.set_ticklabels([])
                            ax[k-1,j].set_xlabel('')
                try:
                    if nrows == 1:
                        out.append(calpntanal(t, fdir=fdir, ant_str=ant_str, calpnt2m=calpnt2m, do_plot=do_plot, ax=ax))
                    else:
                        out.append(calpntanal(t, fdir=fdir, ant_str=ant_str, calpnt2m=calpnt2m, do_plot=do_plot, ax=ax[k]))
                except:
                    out.append({})
            else:
                try:
                    out.append(calpntanal(t, fdir=fdir, ant_str=ant_str, calpnt2m=calpnt2m))
                except:
                    out.append({})
            src = out[-1]
            if src == {}:
                print 'Scan starting at',t.iso,'failed.'
            else:
                if calpnt2m:
                    if k == 0:
                        # Print heading to screen and to file
                        heading = ' Source     Date     Time       HA     Dec   '
                        for iant in idx:
                            heading += '     Ant {:2d}     '.format(iant+1)
                        print heading
                        # Write heading to file only if the file has not been opened/created yet.
                        if not os.path.isfile(outfile):
                            f = open(outfile,'w')
                            f.write(heading+'\n')
                            f.close()
                    rao_deco = ''
                    for j in range(idx.size):
                        rao_deco += ' {:7.3f} {:7.3f}'.format(src['rao'][j],src['deco'][j])
                    line = '{:s} {:s}  {:6.2f}  {:6.2f}'.format(src['source'],
                                                                src['time'].iso[:19],
                                                                src['ha'],src['dec'])
                    line += rao_deco
                else:
                    if k == 0:
                        # Print heading to screen and to file
                        heading = ' Source     Date     Time       HA     Dec        Ant 14'
                        print heading
                        # Write heading to file only if the file has not been opened/created yet.
                        if not os.path.isfile(outfile):
                            f = open(outfile,'w')
                            f.write(heading+'\n')
                            f.close()
                    line = '{:s} {:s}  {:6.2f}  {:6.2f} {:7.3f} {:7.3f}'.format(src['source'],
                                                                                src['time'].iso[:19],
                                                                                src['ha'],src['dec'],
                                                                                src['rao'],src['deco'])
                print line
                f = open(outfile,'a')   # This should create the file if it does not exist, or append otherwise.
                f.write(line+'\n')
                f.close()
    return out
    
def calpntanal(t, fdir=None, ant_str='ant1-13', calpnt2m=False, do_plot=True, ax=None):
    ''' Does a complete analysis of CALPNTCAL, reading information from the SQL
        database, finding the corresponding Miriad IDB data, and doing the 
        gaussian fit to the beam, to return the beam and offset parameters.
        
        t   Time object with a time near the start of the desired scan (finds the scan
              that starts closest time to the given time)

        ax  If specified as a two-element array of axes, the plots will be placed
              in an already existing window, allowing reuse of the same window.
              Otherwise, a new figure is created (if do_plot is True).
        Returns a dictionary containing 
           Source name, Time, HA, Dec, RA offset and Dec offset
           
        NB: This routine has been changed to use feed-rotation-corrected data, by
            commenting out the original code with ##.  Two replacement lines were
            added, below comments marked with ##--.
    '''
    import matplotlib.pyplot as plt
    from matplotlib.transforms import Bbox
    import read_idb
    import dbutil as db
    bl2ord = read_idb.bl2ord
    tdate = t.iso.replace('-','')[:8]
    if fdir is None:
        fdir = get_idbdir(t)+tdate+'/'
    else:
        if fdir[-1] != '/': 
            fdir += '/'
    fdb = dump_tsys.rd_fdb(t)
    # Set offset coordinates appropriate to the type of PROJECTID found
    if calpnt2m:
        # Search for CALPNT2M
        scanidx, = np.where(fdb['PROJECTID'] == 'CALPNT2M')
        if scanidx.size == 0:
            print 'No CALPNT2M project IDs found for date'+t.iso[:10]
            return {}
        else:
            # Found CALPNT2M, so set offset coordinates to match calpnt2m.trj
            rao = np.array([-5.00, -2.0, -1.0, -0.5, 0.00, 0.5, 1.0, 2.0])
            deco = np.array([-2.0, -1.0, -0.5, 0.00, 0.5, 1.0, 2.0, 5.00])
            pltfac = 10.        
    else:
        # Search for CALPNTCAL
        scanidx, = np.where(fdb['PROJECTID'] == 'CALPNTCAL')
        if scanidx.size == 0:
            print 'No CALPNTCAL project IDs found for date'+t.iso[:10]
            return {}
        else:
            # Found CALPNTCAL, so set offset coordinates to match calpnt.trj
            rao = np.array([-1.00, -0.20, -0.10, -0.05, 0.00, 0.05, 0.10, 0.20])
            deco = np.array([-0.20, -0.10, -0.05, 0.00, 0.05, 0.10, 0.20, 1.00])
            pltfac = 1.
    scans,sidx = np.unique(fdb['SCANID'][scanidx],return_index=True)
    tlist = Time(fdb['ST_TS'][scanidx[sidx]].astype(float).astype(int),format='lv')
    idx, = nearest_val_idx([t.jd],tlist.jd)
    filelist = [fdir+f for f in fdb['FILE'][np.where(fdb['SCANID'] == scans[idx])]]
    #import pdb; pdb.set_trace()
    # Read pointing data (timerange t must be accurate)
    outo = read_idb.read_idb(filelist, navg=30)
    # Perform feed rotation correction
    out = read_idb.unrot(outo)
    # Determine wanted baselines with ant 14 from ant_str
    idx = ant_str2list(ant_str)
    ##idx1 = idx[idx>7]  # Ants > 8
    ##idx2 = idx[idx<8]  # Ants <= 8
    # Determine parallactic angle for azel antennas (they are all the same, so find median).  
    #    If 0 < abs(chi) < 30, use channel XX
    #    If 30 < abs(chi) < 60, use sum of channel XX and XY
    #    If 60 < abs(chi) < 90, use channel XY
    ##midtime = Time((out['time'][0] + out['time'][-1])/2.,format='jd')
    ##times, chi = db.get_chi(Time([midtime.lv+1,midtime.lv + 10],format='lv'))
    ##abschi = abs(lobe(np.median(chi[0,0:8])))
    if pltfac == 1.:
        # Case of 27m antenna pointing
        # Do appropriate sums over frequency, polarization and baseline
        ##--This line replaces all of the lines marked with ##
        pntdata = np.sum(np.abs(np.sum(np.sum(out['x'][bl2ord[idx,13],:2,:,:48],1),1)),0)
        ##if abschi >= 0 and abschi < np.pi/6:
        ##    pntdata = np.sum(np.abs(np.sum(out['x'][bl2ord[idx,13],0,:,:48],1)),0)  # Use only XX
        ##elif abschi >= np.pi/6 and abschi < np.pi/3:
        ##    pntdata1 = np.sum(np.abs(np.sum(out['x'][bl2ord[idx1,13],0,:,:48],1)),0)  # Use XX only for ants > 8
        ##    pntdata2 = np.sum(np.abs(np.sum(out['x'][bl2ord[idx2,13],:,:,:48],2)),0)  # Use sum of XX and XY for ants <= 8
        ##    pntdata = pntdata1 + np.sum(pntdata2[np.array([0,2])],0)
        ##else:
        ##    pntdata1 = np.sum(np.abs(np.sum(out['x'][bl2ord[idx1,13],0,:,:48],1)),0)  # Use XX only for ants > 8
        ##    pntdata2 = np.sum(np.abs(np.sum(out['x'][bl2ord[idx2,13],2,:,:48],1)),0)  # Use sum of XY for ants <= 8
        ##    pntdata = pntdata1 + pntdata2
        # Measurements are 90 s long, hence 3 consecutive 30 s points, so do final
        # sum over these
        pntdata.shape = (16,3)
        stdev = np.std(pntdata,1)
        pntdata = np.sum(pntdata,1)
        radat = pntdata[:8]
        decdat = pntdata[8:]
        plsqr, xr, yr = solpnt.gausfit(rao, radat)
        plsqd, xd, yd = solpnt.gausfit(deco, decdat)
        midtime = Time((out['time'][0] + out['time'][-1])/2.,format='jd')
        if (do_plot):
            if ax is None:
                f, ax = plt.subplots(1,2)
                f.set_size_inches(2.5,1.5,forward=True)
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
            ax[0].text(0.05,0.9,'RAO :'+str(plsqr[1])[:5],transform=ax[0].transAxes)
            ax[0].text(0.55,0.9,'FWHM:'+str(plsqr[2])[:5],transform=ax[0].transAxes)
            ax[0].set_xlabel('RA Offset [deg]')
            ax[1].text(0.05,0.9,'DECO:'+str(plsqd[1])[:5],transform=ax[1].transAxes)
            ax[1].text(0.55,0.9,'FWHM:'+str(plsqd[2])[:5],transform=ax[1].transAxes)
            ax[1].set_xlabel('Dec Offset [deg]')
            if ax is None:
                f.suptitle('Pointing on '+out['source']+' at '+midtime.iso)
            else:
                ax[0].set_title(out['source']+' at')
                ax[1].set_title(midtime.iso[:16])
            plt.pause(0.5)
        return {'source':out['source'],'ha':out['ha'][24]*180./np.pi,'dec':out['dec']*180/np.pi,
                   'rao':plsqr[1],'deco':plsqd[1],'time':midtime,'antidx':13}
    else:
        # Case of 2m antenna pointing
        # Do appropriate sum over frequency and polarization but not baseline
        ##--This line replaces all of the lines marked with ##
        pntdata = np.abs(np.sum(np.sum(out['x'][bl2ord[idx,13],:2,:,:48],1),1)) # Combines XX and YY
        ##if abschi >= 0 and abschi < np.pi/6:
        ##    pntdata = np.abs(np.sum(out['x'][bl2ord[idx,13],0,:,:48],1))  # Use only XX
        ##elif abschi >= np.pi/6 and abschi < np.pi/3:
        ##    pntdata1 = np.abs(np.sum(out['x'][bl2ord[idx1,13],0,:,:48],1))  # Use XX only for ants > 8
        ##    pntdata2 = np.abs(np.sum(out['x'][bl2ord[idx2,13],:,:,:48],2))  # Use sum of XX and XY for ants <= 8
        ##    pntdata = np.concatenate((pntdata1,np.sum(pntdata2[:,np.array([0,2])],1)),0)
        ##else:
        ##    pntdata1 = np.abs(np.sum(out['x'][bl2ord[idx1,13],0,:,:48],1))  # Use XX only for ants > 8
        ##    pntdata2 = np.abs(np.sum(out['x'][bl2ord[idx2,13],2,:,:48],1))  # Use sum of XY for ants <= 8
        ##    pntdata = np.concatenate((pntdata1,pntdata2),0)
        # Measurements are 90 s long, hence 3 consecutive 30 s points, so do final
        # sum over these
        pntdata.shape = (idx.size,16,3)
        stdev = np.std(pntdata,2)
        pntdata = np.sum(pntdata,2)
        rao_fit = np.zeros(idx.size,float)
        deco_fit = np.zeros(idx.size,float)
        if (do_plot):
            if ax is None:
                f, ax = plt.subplots(1,2*idx.size)
                f.set_size_inches(2.5*idx.size,1.5,forward=True)
                plt.subplots_adjust(left=0.03, right=0.97, top=0.89, bottom=0.3, wspace=0.1, hspace=0.1)
        for k in range(idx.size):
            radat = pntdata[k,:8]
            decdat = pntdata[k,8:]
            plsqr, xr, yr = solpnt.gausfit(rao, radat)
            plsqd, xd, yd = solpnt.gausfit(deco, decdat)
            midtime = Time((out['time'][0] + out['time'][-1])/2.,format='jd')
            if (do_plot):
                ax[k*2+0].errorbar(rao,radat,yerr=stdev[k,:8],fmt='.')
                ax[k*2+0].plot(xr,yr)
                ax[k*2+0].axvline(x=0,color='k')
                ax[k*2+0].axvline(x=plsqr[1],linestyle='--')
                ax[k*2+1].errorbar(deco,decdat,yerr=stdev[k,8:],fmt='.')
                ax[k*2+1].plot(xd,yd)
                ax[k*2+1].axvline(x=0,color='k')
                ax[k*2+1].axvline(x=plsqd[1],linestyle='--')
                for j in range(2):
                    ax[k*2+j].set_xlim(-3.0, 3.0)
                    ax[k*2+j].grid()
                ax[k*2+0].text(0.05,0.7,'RAO :'+str(plsqr[1])[:5],transform=ax[k*2+0].transAxes,fontsize=9)
                ax[k*2+0].text(0.05,0.5,'FWHM:'+str(plsqr[2])[:5],transform=ax[k*2+0].transAxes,fontsize=9)
                ax[k*2+0].set_xlabel('RAO [deg]',fontsize=9)
                ax[k*2+1].text(-0.2,0.85,'Ant :'+str(idx[k]+1),transform=ax[k*2+1].transAxes,fontsize=9)
                ax[k*2+1].text(0.05,0.7,'DECO:'+str(plsqd[1])[:5],transform=ax[k*2+1].transAxes,fontsize=9)
                ax[k*2+1].text(0.05,0.5,'FWHM:'+str(plsqd[2])[:5],transform=ax[k*2+1].transAxes,fontsize=9)
                ax[k*2+1].set_yticklabels([])
                ax[k*2+1].set_xlabel('DECO [deg]',fontsize=9)
                if k == idx.size-1:
                    ax[k+0].set_title(out['source']+' at',fontsize=9)
                    ax[k+1].set_title(midtime.iso[:16],fontsize=9)
                # Adjust the pairs of subplots to group RA/Dec for each antenna together 
                ax[k*2+0].set_position(Bbox(ax[k*2+0].get_position().get_points() + np.array([[0.0025,0],[0.0025,0]])))
                ax[k*2+1].set_position(Bbox(ax[k*2+1].get_position().get_points() + np.array([[-0.0025,0],[-0.0025,0]])))
                # Set to the same amplitude scale
                ymin = min(ax[k*2+0].get_ylim()[0],ax[k*2+1].get_ylim()[0])
                ymax = max(ax[k*2+0].get_ylim()[1],ax[k*2+1].get_ylim()[1])
                ax[k*2+0].set_ylim(ymin,ymax)
                ax[k*2+1].set_ylim(ymin,ymax)
                plt.pause(0.5)
            rao_fit[k] = plsqr[1]
            deco_fit[k] = plsqd[1]
        return {'source':out['source'],'ha':out['ha'][24]*180./np.pi,'dec':out['dec']*180/np.pi,
                   'rao':rao_fit,'deco':deco_fit,'time':midtime,'antidx':idx}
    
    
def solpntanal(t, udb=False, auto=False, find=True, desat=False):
    ''' Does a complete analysis of SOLPNTCAL, reading information from the SQL
        database, finding and dumping the corresponding Miriad IDB data, and 
        doing gaussian fit to beam to return the beam parameters.  The outputs
        are identical dictionaries for x and y polarizations.
        
        Added find parameter, which if False skips the find_solpnt() call, so
        the Time() object t must be an actual SOLPNTCAL scan time.
    '''
    from copy import deepcopy
    
    def fname2mjd(filename):
        fstem = filename.split('/')[-1]
        fstr = fstem[3:7]+'-'+fstem[7:9]+'-'+fstem[9:11]+' '+fstem[11:13]+':'+fstem[13:15]+':'+fstem[15:17]
        return Time(fstr).mjd

    pnt = solpnt.get_solpnt(t, find=find)
    proc = solpnt.process_solpnt(pnt)
    trange = Time([pnt['Timestamp'],pnt['Timestamp']+300.],format='lv')
    skycal = None
    if trange[0].mjd < 57450:
        if auto: print "Warning: 'auto' keyword does nothing for older data.  Keyword ignored."
        otp = dump_tsys.rd_miriad_tsys(trange,udb=udb)
    elif trange[0].mjd < 58623:
        # The SKYCALTEST scans may have been taken, but they were not setting the DCM table correctly
        # until 2019-05-20, so earlier dates just do this almost as good approach.
        poor_skycal = skycal_anal(trange[0], do_plot=False, desat=desat)
        otp = dump_tsys.rd_miriad_tsys_16(trange, udb=udb, auto=auto, tref=trange[0], skycal=poor_skycal, desat=desat)
    else:
        # New code as of 2019-May-21, tries to read from SKYCALTEST and use it
        # for a far off-Sun background measurement.
        skycal = skycal_anal(trange[0], do_plot=False, last=True, desat=desat)
        otp = dump_tsys.rd_miriad_tsys_16(trange, udb=udb, auto=auto, tref=trange[0], skycal=skycal, desat=desat)
    if skycal:
        newskycal = deepcopy(skycal)
        if auto:
            # Modify the skycal to put the auto information into the standard location
            newskycal['offsun'] = np.real(skycal['offsun_auto'])
            newskycal['rcvr_bgd'] = np.real(skycal['rcvr_bgd_auto'])
    else:
        newskycal = None
    fghz = otp['fghz']
    xra,xdec,xrao,xdeco = solpnt.process_tsys(otp,proc,pol=0,skycal=newskycal)
    x = {'ut_mjd':otp['ut_mjd'],'fghz':fghz,'ra0':proc['ra0'],'dec0':proc['dec0'],'raparms':xra,'decparms':xdec,
         'rao':xrao,'deco':xdeco,'ra':proc['rao'],'dec':proc['deco'],'antlist':proc['antlist']}
#    otp = solpnt.rd_tsys(yfile,sfile)
    yra,ydec,yrao,ydeco = solpnt.process_tsys(otp,proc,pol=1,skycal=newskycal)
    y = {'ut_mjd':otp['ut_mjd'],'fghz':fghz,'ra0':proc['ra0'],'dec0':proc['dec0'],'raparms':yra,'decparms':ydec,
         'rao':yrao,'deco':ydeco,'ra':proc['rao'],'dec':proc['deco'],'antlist':proc['antlist']}
    qual = sp_check_qual(x,y)
    return x,y, qual
        
def skycal_anal(t=None, do_plot=False, last=False, desat=False, file=''):
    ''' Reads the daily SKYCALTEST data and returns the far "off Sun" background for use with
        total power calibration.  The SKYCALTEST observation can double as a GAINCAL if needed,
        but this has not yet been implemented.
        Returns the receiver background power and the backend power (with 62 dB front-end attenuation
        applied), for both total power and auto-correlation data.
        
        Inputs:
          t        An optional Time() object with the date for which to return the SKYCAL data.  Uses
                     today's date if omitted (or None)
          do_plot  A flag requesting a summary plot of the results (default = True) 
                     NB: Only total power is plotted, not auto-correlation
          last     A flag specifying that the last SKYCALTEST scan should be used (normally only one
                     per day is obtained, so this allows a second one to be processed).  Default is
                     False, meaning only the first SKYCALTEST scan is used.  Note that if there is
                     only one SKYCALTEST scan, first and last are the same and this flag has no effect.
          
        Returns:
          Dictionary with keys (or None):
          
          'offsun'         Total power level on blank sky, with femattn level 0 (includes receiver noise)
          'rcvr_bgd'       Total power level on blank sky, with femattn level 15 (receiver noise only)
          'offsun_auto'    Auto correlation level on blank sky, with femattn level 0
          'rcvr_bgd_auto'  Auto correlation level on blank sky, with femattn level 15
          
    '''
    from gaincal2 import get_fem_level
    import read_idb as ri
    import socket


    #    lev = get_fem_level(trange,300)
    #    levs = np.zeros((13,2), dtype=int)
    #    levs[:,0] = lev['hlev'][:13,0]
    #    levs[:,1] = lev['vlev'][:13,0]
    if t is None:
        t = Time.now()
    if t.mjd < 58623:
        print 'SKYCAL_ANAL: No valid SKYCALTEST scans for dates prior to 2019-05-20'
        print 'SKYCAL_ANAL: Using GAINCALTEST 62 dB setting instead.'
        # Do a "poor-man's" skycal background by using the gaincal 62 dB setting.
        # This should be almost as good for receiver noise subtraction as a true SKYCAL
        from attncal import get_attncal
        outdict = get_attncal(t)
        return {'rcvr_bgd': outdict[0]['rcvr'], 'rcvr_bgd_auto': outdict[0]['rcvr_auto']}
    fdb = dump_tsys.rd_fdb(t)
    gcidxes, = np.where(fdb['PROJECTID'] == 'SKYCALTEST')
    host = socket.gethostname()
    if len(gcidxes) != 0:
        if last:
            gcidx = gcidxes[-1]
        else:
            gcidx = gcidxes[0]

        datadir = get_idbdir(t)
        # Add date path if on pipeline
        # if datadir.find('eovsa') != -1:
        #     datadir += fdb['FILE'][gcidx][3:11]+'/'
        if host == 'pipeline': datadir += fdb['FILE'][gcidx][3:11]+'/'

        if file == '':
            file = datadir+fdb['FILE'][gcidx]
        
        out = ri.read_idb([file], desat=desat)
        skytrange = Time(out['time'][[0,-1]],format='jd')
        nant, npol, nf, nt = out['p'].shape
        if nant > 13: nant = 13
        offsun = np.zeros((nant,npol,nf), dtype=float)
        rcvr = np.zeros((nant,npol,nf), dtype=float)
        offsun_auto = np.zeros((nant,npol,nf), dtype=complex)
        rcvr_auto = np.zeros((nant,npol,nf), dtype=complex)
        # Get the SKYCALTEST fem level vs. time and use it to define
        # the range of times for femattn 15 (62 dB) and femattn 0 (0 dB)
        hvlev = ['hlev','vlev']
        lev = get_fem_level(skytrange)
        # Get indexes to common times
        idx1, idx2 = common_val_idx(lev['times'].jd,out['time'])
        for i in range(nant):
            for j in range(npol):
                idx0, = np.where(lev[hvlev[j]][i,idx1] == 0)   # 0 dB
                idx15, = np.where(lev[hvlev[j]][i,idx1] == 15) # 62 dB
                # NB: time and frequency indexes are reordered for some reason.
                # The next two lines are a median over times, despite what it seems.
                offsun[i,j] = np.nanmedian(out['p'][i,j,:,idx2[idx0]],0)
                rcvr[i,j] = np.nanmedian(out['p'][i,j,:,idx2[idx15]],0)
                offsun_auto[i,j] = np.nanmedian(out['a'][i,j,:,idx2[idx0]],0)
                rcvr_auto[i,j] = np.nanmedian(out['a'][i,j,:,idx2[idx15]],0)
        if do_plot:
            nrow = 2
            ncol = (nant+1)/2
            f, ax = plt.subplots(nrow, ncol, sharex='col', sharey='row')
            f.set_size_inches(2*nant,7,forward=True)
            f.suptitle('HPol Receiver and Sky Noise for SKYCAL scan at '+skytrange[0].iso[:19]+' UT',fontsize=18)
            for ant in range(nant):
                ax[ant % nrow, ant/nrow].set_title('Ant '+str(ant+1)+' Noise')
                if ant % ncol == 0:
                    ax[ant/ncol,0].set_ylabel('Arbitrary Units')
                if ant >= nant/nrow:
                    ax[1,ant - nant/nrow].set_xlabel('Frequency [GHz]')
                ax[ant % nrow, ant/nrow].plot(offsun[ant,0])
                ax[ant % nrow, ant/nrow].plot(rcvr[ant,0])
            f, ax = plt.subplots(nrow, ncol, sharex='col', sharey='row')
            f.set_size_inches(2*nant,7,forward=True)
            f.suptitle('VPol Receiver and Sky Noise for SKYCAL scan at '+skytrange[0].iso[:19]+' UT',fontsize=18)
            for ant in range(nant):
                ax[ant % nrow, ant/nrow].set_title('Ant '+str(ant+1)+' Noise')
                if ant % ncol == 0:
                    ax[ant/ncol,0].set_ylabel('Arbitrary Units')
                if ant >= nant/nrow:
                    ax[1,ant - nant/nrow].set_xlabel('Frequency [GHz]')
                ax[ant % nrow, ant/nrow].plot(offsun[ant,1])
                ax[ant % nrow, ant/nrow].plot(rcvr[ant,1])
    else:
        print 'SKYCAL_ANAL: No SKYCALTEST scan for this date'
        print 'SKYCAL_ANAL: Using GAINCALTEST 62 dB setting instead.'
        # Do a "poor-man's" skycal background by using the gaincal 62 dB setting.
        # This should be almost as good for receiver noise subtraction as a true SKYCAL
        from attncal import get_attncal
        outdict = get_attncal(t)
        return {'rcvr_bgd': outdict[0]['rcvr'], 'rcvr_bgd_auto': outdict[0]['rcvr_auto']}
    return {'offsun': offsun, 'rcvr_bgd': rcvr, 'offsun_auto': offsun_auto, 'rcvr_bgd_auto': rcvr_auto}
    
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
    files = glob.glob('/dppdata1/TPCAL/*')
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

def sp_apply_cal2(out, calfac, offsun, auto=False):
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
    if auto:
        for i in range(nt):
            out['a'][:nant,:2,good,i] = (out['a'][:nant,:2,good,i]-osun[:,:,good])*cfac[:,:,good]
    else:
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
    
def sp_bg_subtract2(out, idx, auto=False):
    ''' Given "standard" output of tp_display.rd_tsys_multi(),
        and a range of time indexes idx into the arrays, use the
        median of data at times indicated by idx as background.
        Subtract that background from both xtsys and ytsys
    '''
    nant, npol, nf, nt = out['p'].shape
    # Create the backgrounds (shape nf,nant)
    if len(idx) == 1:
        # Only one index in idx, so no median necessary
        if auto:
            bg = out['a'][:,:,:,idx]
        else:
            bg = out['p'][:,:,:,idx]
    else:
        # Multiple values in idx, so do median
        if auto:
            bg = np.median(out['a'][:,:,:,idx],3)
        else:
            bg = np.median(out['p'][:,:,:,idx],3)
    # Perform the background subtraction
    if auto:
        for i in range(nt):
            out['a'][:,:,:,i] -= bg
    else:
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
    #fgood = np.where(x['fghz'] > 2.48)[0]
    fgood = np.where(x['fghz'] > 1.10)[0]
    print len(fgood)
    for ant in range(nant):
        ax[ant % nrow, ant/nrow].set_title('Ant '+str(ant+1)+' [blue=RA, red=Dec]')
        if ant % ncol == 0:
            ax[ant/ncol,0].set_ylabel('Beam FWHM [deg]')
        if ant >= ncol:
            ax[1,ant - ncol].set_xlabel('Frequency [GHz]')
        ax[ant % nrow, ant/nrow].set_ylim(0.5,20)
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
        ax[ant % nrow, ant/nrow].set_ylim(0.5,20)
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
    #fgood = np.where(x['fghz'] > 2.48)[0]
    fgood = np.where(x['fghz'] > 1.10)[0]
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
            # clip errors to 0.5 degree
            if  abs(ramed) > 5000: ramed = np.sign(ramed)*5000
            if abs(decmed) > 5000: decmed = np.sign(ramed)*5000
            xelx.append(ramed/10000)
            elx.append(decmed/10000)
        else:
            # clip errors to 0.5 degree
            halfdeg = 0.5*np.pi/180.
            if abs(xel) > halfdeg: xel = np.sign(xel)*halfdeg
            if abs(el) > halfdeg: el = np.sign(el)*halfdeg
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
            # clip errors to 0.5 degree
            if  abs(ramed) > 5000: ramed = np.sign(ramed)*5000
            if abs(decmed) > 5000: decmed = np.sign(ramed)*5000
            xely.append(ramed/10000)
            ely.append(decmed/10000)
        else:
            # clip errors to 0.5 degree
            halfdeg = 0.5*np.pi/180.
            if abs(xel) > halfdeg: xel = np.sign(xel)*halfdeg
            if abs(el) > halfdeg: el = np.sign(el)*halfdeg
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

def sp_explore(x, y, aidx, fidx, ax=None):
    ''' Makes a 2x2 plot for one antenna and frequency, showing the data as dots,
        the fit as an orange line, and the nominal beam as green for comparison.
    
        To reuse the same plot, call first as 
            ax = sp_explore(x,y,ant,frq),
        then as
            ax = sp_explore(x,y,ant,frq,ax=ax)
        It is MUCH faster and more efficient that way!
    '''
    if ax is None:
        f, ax = plt.subplots(2,2)
        import disk_conv as dc
        fghz, a, aout = dc.disk_conv(x['fghz'])
        ax[0,0].aout = aout
    anom = ax[0,0].aout[fidx]/(2*np.sqrt(np.log(2)))*10000.
    lbl = [[' X-RA',' X-DEC'],[' Y-RA',' Y-DEC']]
    axis = ['ra','dec']
    ax[0,0].get_figure().suptitle('Ant '+str(aidx+1)+' '+str(x['fghz'][fidx])[:4]+' GHz')
    for i,src in enumerate([x,y]):
        for j,a in enumerate(axis):
            ax[i,j].cla()
            ax[i,j].plot(np.concatenate(([-100000],src[a],[100000]))/10000.,src[a+'o'][fidx,:,aidx],'.',label=lbl[i][j])
            xarr = np.linspace(-100000.,100000,1000)
            s1,x0,a1,o1 = src[a+'parms'][:,fidx,aidx]   # pk flux, offset, 1/e half-width
            yarr = s1*np.exp(-((xarr-x0)/a1)**2) + o1
            ynom = s1*np.exp(-(xarr/anom)**2) + o1
            ax[i,j].plot(xarr/10000.,yarr,label='Fit')
            ax[i,j].plot(xarr/10000.,ynom,label='Nominal')
            ax[i,j].legend(loc='lower right', fontsize=9)
    return ax

def sp_gaussfit(x, y, aidx, fidx, ax=None):
    from scipy.optimize import curve_fit
    def func(x, a, b, c, d):
        return a*np.exp(-((x-b)/c)**2) + d
    if ax is None:
        f, ax = plt.subplots(2,2)
        import disk_conv as dc
        fghz, a, aout = dc.disk_conv(x['fghz'])
        ax[0,0].aout = aout
    anom = ax[0,0].aout[fidx]/(2*np.sqrt(np.log(2)))*10000.
    lbl = [[' X-RA',' X-DEC'],[' Y-RA',' Y-DEC']]
    axis = ['ra','dec']
    ax[0,0].get_figure().suptitle('Ant '+str(aidx+1)+' '+str(x['fghz'][fidx])[:4]+' GHz')
    for i,src in enumerate([x,y]):
        for j,a in enumerate(axis):
            ax[i,j].cla()
            xdata = np.concatenate(([-100000],src[a],[100000]))
            ydata = src[a+'o'][fidx,:,aidx]
            ax[i,j].plot(xdata/10000.,ydata,'.',label=lbl[i][j])
            low_bounds =  [0.1*np.max(ydata), -30000., anom*0.9, np.min(ydata)-0.1*np.max(ydata)]
            high_bounds = [1.0*np.max(ydata),  30000., anom*1.1, np.min(ydata)+0.1*np.max(ydata)]
            popt, pcov = curve_fit(func, xdata, ydata, bounds=(low_bounds,high_bounds))
            xarr = np.linspace(-100000.,100000,100)
            s1,x0,a1,o1 = src[a+'parms'][:,fidx,aidx]   # pk flux, offset, 1/e half-width
            yarr = s1*np.exp(-((xarr-x0)/a1)**2) + o1
            ynom = s1*np.exp(-(xarr/anom)**2) + o1
            yfit = func(xarr, *popt)
            ax[i,j].plot(xarr/10000.,yarr,label='Fit')
            ax[i,j].plot(xarr/10000.,yfit,label='New Fit')
            ax[i,j].legend(loc='lower right', fontsize=9)
            print 'Fit       :  {:10.3f}, {:10.3f}, {:10.3f}, {:10.3f}'.format(s1, x0, a1, o1)
            print 'New    Fit:  {:10.3f}, {:10.3f}, {:10.3f}, {:10.3f}'.format(*popt)
            print 'Min Bounds: {:10.3f}, {:10.3f}, {:10.3f}, {:10.3f}'.format(*low_bounds)
            print 'Max Bounds: {:10.3f}, {:10.3f}, {:10.3f}, {:10.3f}'.format(*high_bounds)
    return ax

    

def send_cmds(cmds,acc):
    ''' Sends a series of commands to ACC.  The sequence of commands
        is not checked for validity!
        
        cmds   a list of strings, each of which must be a valid command
    '''
    import socket

    for cmd in cmds:
        #print 'Command:',cmd
        s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        try:
            s.connect((acc['host'],acc['scdport']))
            s.send(cmd)
            time.sleep(0.01)
            s.close()
        except:
            print 'Error: Could not send command',cmd,' to ACC.'
    return
    
def offsets2ants(t,xoff,yoff,ant_str=None):
    ''' Given a start time (Time object) and a list of offsets output by sp_offsets()
        for 13 antennas, convert to pointing coefficients (multiply by 10000),
        add to coefficients listed in stateframe, and send to the relevant 
        antennas.  The antennas to update are specified with ant_str 
        (defaults to no antennas, for safety).        
    '''

    oldant = [8,9,10,12]
    if ant_str is None:
        print 'No antenna list specified, so there is nothing to do!'
        return

    try:
        timestamp = int(Time(t,format='mjd').lv)
    except:
        print 'Error interpreting time as Time() object'
        return
    from util import ant_str2list
    import dbutil as db
    antlist = ant_str2list(ant_str)
    if antlist is None:
        return
    cursor = db.get_cursor()
    # Read stateframe data at timestamp (as of 10 s ago)
    D15data = db.get_dbrecs(cursor,dimension=15,timestamp=timestamp,nrecs=1)
    p1_cur, = D15data['Ante_Cont_PointingCoefficient1']
    p7_cur, = D15data['Ante_Cont_PointingCoefficient7']
    err_rtn = False
    for i in antlist:
        if p1_cur[i] == 0 and p7_cur[i] == 0:
            # This is probably a glitch where the SQL record for this ant is all zeros
            print 'Ant'+str(i+1)+' SQL record returned zeros.  Adjust the time and try again.'
            err_rtn = True
    if err_rtn:
        print 'No updates sent to antennas'
        return
    
    import stateframe as stf
    accini = stf.rd_ACCfile()
    acc = {'host': accini['host'], 'scdport':accini['scdport']}
    for i in antlist:
        if i in oldant:
            # Change sign of RA offset to be HA, for old antennas (9, 10, 11 or 13)
            p1_inc = int(-xoff[i]*10000)
        else:
            p1_inc = int(xoff[i]*10000)
        p7_inc = int(yoff[i]*10000)
        p1_new = p1_cur[i] + p1_inc
        p7_new = p7_cur[i] + p7_inc
        print 'Updating P1 for Ant',i+1,'P1_old =',p1_cur[i],'P1_inc =',p1_inc,'P1_new =',p1_new
        cmd1 = 'pointingcoefficient1 '+str(p1_new)+' ant'+str(i+1)
        print 'Updating P7 for Ant',i+1,'P7_old =',p7_cur[i],'P7_inc =',p7_inc,'P7_new =',p7_new
        cmd7 = 'pointingcoefficient7 '+str(p7_new)+' ant'+str(i+1)
        print 'Commands to be sent:'
        print cmd1
        print cmd7
        send_cmds([cmd1],acc)
        send_cmds([cmd7],acc)
        
def solpnt2sql(t,tsql=None,prompt=True, desat=False):
    ''' Handles the complete processing of a solar calibration, from analysis through
        writing the calibration data to the SQL database.  The SQL timestamp at which
        to write the record can be given as Time() object in tsql, or writes to
        "standard" location of 20 UT on date of observation, if omitted (preferred!).
    '''
    import cal_header as ch
    x1, y1, qual1 = solpntanal(t)
    calfac1, offsun1 = sp_get_calfac(x1, y1)
    x2, y2, qual2 = solpntanal(t,auto=True, desat=desat)
    calfac2, offsun2 = sp_get_calfac(x2, y2)
    tpcal_dict = {'fghz':x1['fghz'],'timestamp':Time(x1['ut_mjd'][0],format='mjd').lv,
                  'tpcalfac':calfac1,'accalfac':calfac2,'tpoffsun':offsun1,'acoffsun':offsun2}
    if tsql is None:
        # Default to writing SQL record at 20 UT on date of observation 
        tsql = Time(int(x1['ut_mjd'][0]) + 20/24.,format='mjd')
    if prompt:
        ok = raw_input('Okay to write result to SQL database? [Y/N]: ')
    else:
        ok = 'Y'
    if ok.upper() == 'Y':
        ch.tpcal2sql(tpcal_dict,t=tsql)
        print 'Result was written to the SQL database'
    else:
        print 'Result was NOT written to the SQL database'

def check_qual(x, qual):
    status = np.zeros(qual.shape,'S4')
    status[np.where(qual)] = 'Good'
    status[np.where(qual==False)] = '*Bad'
    nparms,nf,nant = x['raparms'].shape
    t = Time(x['ut_mjd'][0],format='mjd')
    print t.iso[:19],': Quality of TP Calibration'
    print '    Ant      X-Feed        Y-Feed'
    print '            RA    DEC     RA    DEC'
    for i in range(nant):
        print '    ',i+1,' ',status[:,i]
    percent_good = len(np.where(qual)[0])*100/(4*nant)
    return percent_good

def best_solpnt2sql(t):
    ''' For post-date analysis, to examine the (typically) two SOLPNTCAL scans
        on a date given by Time() object t and send the best one to the SQL database.
    '''
    import cal_header as ch
    ts, timestamp = solpnt.find_solpnt(t)
    pct_list = []   # List of percent-good quality (best will be used)
    tp_list = []    # List of results of solpntanal() for total power
    ac_list = []    # List of results of solpntanal() for autocorrelation
    t_list = []     # List of time objects
    for tstamp in ts:
        tobj = Time(tstamp,format='lv')
        x1, y1, qual1 = solpntanal(tobj, find=False)
        pct = check_qual(x1, qual1)
        if pct > 50:
           # More than 50% good, so it is worth continuing
           tp_list.append([x1, y1, qual1])
           x2, y2, qual2 = solpntanal(tobj, auto=True, find=False)
           ac_list.append([x2, y2, qual2])
           pct_list.append(pct)
           t_list.append(tobj)
    if pct_list != []:
        best_pct = max(pct_list)
        best_idx = pct_list.index(best_pct)
        x1, y1, qual1 = tp_list[best_idx]
        calfac1, offsun1 = sp_get_calfac(x1, y1, do_plot=False)
        x2, y2, qual2 = ac_list[best_idx]
        calfac2, offsun2 = sp_get_calfac(x2, y2, do_plot=False)
        tpcal_dict = {'fghz':x1['fghz'],'timestamp':Time(x1['ut_mjd'][0],format='mjd').lv,
                      'tpcalfac':calfac1,'accalfac':calfac2,'tpoffsun':offsun1,'acoffsun':offsun2}
        # Default to writing SQL record at 20 UT on date of observation 
        tsql = Time(int(x1['ut_mjd'][0]) + 20/24.,format='mjd')
        ch.tpcal2sql(tpcal_dict,t=tsql)
        print 'Result for',t_list[best_idx].iso[:16],'UT, with',str(best_pct)+'% quality, was written to the SQL database'
    else:
        print 'No good scans for',t.iso[:10]+'. Update manually with a result from another date.'
   
if __name__ == "__main__":
    ''' Run automatically via cron job, or at command line.
        Usage: python /common/python/current/calibration.py "2014-12-15 18:30"
        where the time string is optional.  If omitted, the current time is used.
        
        The logic is to check whether there is a new SOLPNTCAL available, and
        analyze it if so.  Compares times from solpnt.find_solpnt() with current
        time and analyzes any that are between 5 and 10 minutes old.
        
        Normally this is run from the DPP, in which case the results are now
        written to the SQL database.  The old method of writing to a file is
        now to be discouraged, and probably will not work on pipeline anyway.
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
    if len(times) == 0:
        # No SOLPNTCAL scans (yet)
        print t.iso[:19]+': No SOLPNTCAL scans for today'
        exit()
    #elif type(times[0]) is np.ndarray:
    #    # Annoyingly necessary when only one time in tstamps
    #    times = times[0]
    igt5 = np.where((timestamp - times) > 300)[0]
    if (Time.now().mjd - t.mjd) < 1:
        # This has a chance of being a "real-time" calibration, so make some further checks
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
 #   if socket.gethostname() == 'pipeline':
 #       x, y, qual = solpntanal(t,udb=True)
 #   elif socket.gethostname() == 'dpp':
    x, y, qual = solpntanal(t,udb=False)  # NB: the desat keyword is irrelevant because udb=False
    xout,yout,dxout,dyout = sp_offsets(x,y,save_plot=True)
 #   else:
 #       print 'CALIBRATION Error: This routine only runs on dpp or pipeline.'
 #       exit()
    percent_good = check_qual(x, qual)
    if percent_good > 50:
#        if socket.gethostname() == 'dpp':
#            # If this is the DPP, write the results to the SQL database automatically.
        solpnt2sql(t,prompt=False, desat=True)
#        else:
#            # This is the old way of writing results to a disk file--should be removed...
#            calfac, offsun = sp_get_calfac(x,y)
#            # If another file for today's date already exists, this will overwrite it
#            sp_write_calfac(x,y,calfac,offsun)
#            print 'Calibration file successfully written'
#            if t.iso[:10] == Time.now().iso[:10]:
#               # The calibration file is for today, so rewrite as tomorrow's file also
#               x['ut_mjd'][0] = x['ut_mjd'][0]+1.   # Add one day to first timestamp
#               sp_write_calfac(x,y,calfac,offsun)
#               print "Also wrote tomorrow's file"
    else:
        print 'Calibration file not written--too many bad values.'
    exit()  
