'''
   Module for accessing, analyzing and plotting calibration data'''
#
# History:
#   2014-Dec-16  DG
#     Added optional fghz input, to calculate the curve at specified
#     frequencies.
#
import numpy as np
from scipy.optimize import leastsq
from solpnt import gausfit, plt

def disk_conv(fghz=None,doplot=False):
    ''' Calculate nominal beam size for 2.1 m EOVSA antennas
        and the measured solar disk size when convolved with
        a nominal solar disk of size of 0.54 degrees.
        Returns:
           fghz  Frequency range in GHz
           a     Nominal FWHM of 2.1 m primary beam
           aout  FWHM of measured solar disk
    '''        
    if fghz is None:
        fghz = np.arange(101)*17./100. + 1.
    a = 1.22*(180./np.pi)*30./(fghz*210)
    aout = np.zeros(len(fghz),'float')
    alpha = a/(2*np.sqrt(np.log(2.)))
    for i in range(len(fghz)):
        x = np.arange(-10,10.01,0.01)
        disk = np.ones(55,'float')
        g = np.exp(-(x**2)/(alpha[i]**2))
        gdisk = np.convolve(g,disk,mode='same')
        [A,x0,w,b],x0,y0 = gausfit(x,gdisk)
        aout[i] = np.abs(w)*2*np.sqrt(np.log(2.))
    if doplot:
        plt.figure()
        plt.plot(fghz,a)
        plt.plot(fghz,aout)
        plt.xlabel('Frequency [GHz]')
        plt.ylabel('Disk vs. Beam Size [degrees]')
        plt.title('Unconvolved Beam (blue); Disk-Convolved (green)')
        plt.yscale('log')
        plt.xscale('log')
        plt.ylim(0.5,10)
        plt.xlim(1,18)
    return fghz,a,aout