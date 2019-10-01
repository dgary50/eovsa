#
# EOVSA FITS-handling Routines
#
# These routines are used to read and plot FITS files created for EOVSA.
#
# History
#  2019-08-05  DG
#    Initial start of file (just eovsa_fits2plot() routine)
#  2019-08-10  DG
#    Added return of dictionary of data from combined files.  Changed
#    eovsa_fits2plot() name to eovsa_combinefits().
#  2019-08-12  DG
#    Remove nans before checking vmax in scaling the imshow() colors

from astropy.io import fits
import glob
from util import Time
import numpy as np
import matplotlib.pylab as plt
from matplotlib.dates import DateFormatter

def eovsa_combinefits(files, freqgaps=True, outpath=None, ac_corr=True, doplot=True):
    ''' Reads provided list of FITS files and combines them into a single,
        all-day dynamic spectrum.  Returns a dictionary with the spectrum, 
        times and frequencies.  Optionally writes the combined FITS files 
        to a single output file.  Optionally makes a nice spectrogram plot.
        
        If the spectrum type is Cross Power (header type 2), the dictionary
        has keys:
           'x'    Cross-Power spectrum of size nf, nt
           'time' UT time as Julian Date, of size nt
           'fghz' Frequency in GHz, of size nf
           'src'  Source name from header OBJ_ID
         
        If the spectrum type is anothe other type, the dictionary has keys:
           'p'    Total-Power spectrum of size nf, nt
           'time' UT time as Julian Date, of size nt
           'fghz' Frequency in GHz, of size nf
           'src'  Source name from header OBJ_ID

        Optional keyword:
           freqgaps  Boolean. If True (default), missing frequencies will show as
                       a gap in the plot (does not affect returned spectrum).
           ac_corr   Boolean. IF True (default) and the spectrum type is Total Power
                       (header type 1) the variable background due to the air 
                       conditioning is (partially) corrected.  Has no effect if
                       header type is not 1.
           outpath   File path of output FITS file.  Subdirectories with year/month/day
                       are created if they do not exist, and output filename has 
                       standard format as EOVSA_TPall_yymmddhhmm.fts if total power,
                       or EOVSA_Xall_yymmddhhmm.fts.  If None (default), no output 
                       FITS file is generated.
           doplot    Boolean.  If True (default), a nice spectrogram plot is
                       displayed on the screen.
           
    '''
    for file in files:
        spec = fits.getdata(file,ext=0)
        freq = fits.getdata(file,ext=1)
        ut = fits.getdata(file,ext=2)
        fghz = freq['sfreq']
        time = Time(ut['mjd']+ut['time']/86400000.,format='mjd')
        jd = time.jd
        pd = time.plot_date
        if file != files[0]:
            specs = np.concatenate((specs,spec),1)
            jds = np.concatenate((jds,jd))
            pds = np.concatenate((pds,pd))
        else:
            # Things to set for the first file
            specs = spec
            jds = time.jd
            pds = time.plot_date
            date = time[0].iso[:10]
            header = fits.open(file)[0].header
            if header['TYPE'] == 0:
                typstr = 'Undefined'
            elif header['TYPE'] == 1:
                typstr = 'Total Power'
            elif header['TYPE'] == 2:
                typstr = 'Cross Power'
            src = header['OBJ_ID']
            
    # Create output dictionary
    if typstr == 'Total Power':
        out = {'p':specs,'time':jds,'fghz':fghz,'source':src}
        if ac_corr:
            # Correct for varying air conditioning temperature
            from autocorrect_tp import tp_bgnd
            bgnd = tp_bgnd(out)
            if bgnd is None:
                pass
            else:
                out['p'] -= bgnd
                specs = out['p']
    else:
        out = {'x':specs,'time':jds,'fghz':fghz,'source':src}

    if outpath:
        import os
        from xspfits2 import tp_writefits
        if typstr == 'Total Power':
            # Write an all-day FITS file
            tp_writefits(out, out['p'].astype(np.float32), filestem='TPall_',outpath=outpath)
        else:
            # Write an all-day FITS file
            tp_writefits(out, out['x'].astype(np.float32), filestem='Xall_',outpath=outpath)
            
    if doplot:
        f, ax = plt.subplots(1,1,figsize=(14,5))
        # Set any time gaps to 0
        tdif = pds[1:] - pds[:-1]
        bad, = np.where(tdif > 120./86400)  # Time gaps > 2 minutes
        specs[:,bad] = 0
        # Overwrite bottom of frequency gap with 0
        if freqgaps is True:
            fdif = fghz[1:] - fghz[:-1]
            bad, = np.where(fdif > 0.1)
            specs[bad,:] = 0
        else:
            if time[0].mjd < 58536:
                # If the date is earlier than 2019-02-22, eliminate frequency gaps
                # (for display only) by smoothing the frequency list
                nf = len(fghz)
                p = np.polyfit(np.arange(nf),fghz,6)
                fghz = np.polyval(p,np.arange(nf))
        X = np.sort(specs.flatten())   # Sorted, flattened array
        X = X[np.where(~np.isnan(X))]  # Removes any nan at end of the sorted array
        vmax = X[int(len(X)*0.95)]  # Clip at 5% of points
        
        im = ax.pcolormesh(pds,fghz,specs,vmax=vmax,vmin=0)
        #    plt.colorbar(im,ax=ax,label='Amplitude [arb. units]')
        ax.xaxis_date()
        ax.xaxis.set_major_formatter(DateFormatter("%H:%M"))
        ax.set_ylim(fghz[0], fghz[-1])
        ax.set_xlabel('Time [UT]')
        ax.set_ylabel('Frequency [GHz]')
        ax.set_title('EOVSA '+typstr+' Dynamic Spectrum for '+date)
        f.autofmt_xdate(bottom=0.15)
    return out