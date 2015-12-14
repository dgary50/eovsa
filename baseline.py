#  Baseline.py (routines for measuring and analyzing baselines)
#
# History:
#   2015-07-26  DG
#     First written
#  
import eovsa_array as ea
import eovsa_lst as el
import aipy
import numpy as np

def bl_phz(bn,be,fghz,times):
    ''' Initial test routine to play with baseline corrections.  This 
        routine takes a baseline error in N direction, bn, and baseline
        error in E direction, be, as well as a list of frequencies and
        times corresponding to an actual observation, and returns the
        phase error associated with the baseline error.  The phase errors
        are relative to the phases at the first time, for easy comparison
        with "calibrated" phases in actual data.
        Usage:
            dphz = bl_phz(bn,be,fghz,times)
        '''
    nf = len(fghz)
    nt = len(times)
    aa = ea.eovsa_array()
    src = aipy.amp.RadioSpecial('Sun')
    dphz = np.zeros((nf,nt),'float')
    for i,t in enumerate(times):
        aa.set_jultime(t.jd)
        src.compute(aa)
        ha = el.eovsa_ha(src,t)
        dec = src.dec
        bx = -bn * np.sin(aa.lat)
        by = be
        bz = bn*np.cos(aa.lat)
        dphz[:,i] = 2*np.pi*(bx*np.cos(dec)*np.cos(ha) - by*np.cos(dec)*np.sin(ha) + bz*np.sin(dec))*fghz/0.3

    # Normal phases to first time
    blah2 = copy.copy(dphz)
    for i in range(nt):
        blah2[:,i] = blah2[:,i]-dphz[:,0]
    dphz = blah2
    # Ensure the phases are between -pi and pi
    dphz = dphz % (2*np.pi)
    dphz[where(dphz > np.pi)] -= 2*np.pi
    return dphz
