#
# tp_display routines
#
#   2014-Dec-18  DG
#     First written.
#   2014-Dec-20  DG
#     Fixed problems with display of nonuniform image data in 
#     tsys_show_dynspec(), especially the time axis.
#   2015-Apr-10  DG
#     Fixed a problem introduced earlier where I converted out['ut'] to
#     integer, which killed the date plotting.  I covert it back to float
#     now in tsys_show_dynspec()
#   2015-Apr-18  DG
#     In tsys_show_dynspec(), subtract day1 instead of using % 1 to allow
#     times to span a day.  This is still not entirely satisfactory,
#     and a considerable rewrite of date/time code is necessary
#   2015-May-29  DG
#      Converted from using datime() to using Time() based on astropy
#   2015-Jun-02  DG
#      Greatly simplified to just call the offline.py routines, which anyway
#      were the same.
#

import dump_tsys as dtsys     
import time
import offline

  
def rd_tsys_multi(trange):
 
    ''' Given a timerange (2-element Time() array), dump and
        return the tsys data for all IDB files in that timerange.
    '''
    # Just call dump_tsys to create the xt*.txt and yt*.txt files
    dtsys.dump_tsys(trange)
    time.sleep(3)
    # Now call the offline routine
    return offline.rd_tsys_multi(trange)
    
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
    offline.tsys_show_dynspec(out,idx,ampscl,domedian,frq)
