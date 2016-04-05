#
# Routine to read FEM calibration measurement files and fit
# log(voltage) vs. power with a quartic, in order to generate
# coefficients to be entered into the crio.ini files.
#
#   2016-Mar-19  DG
#     First written.
#

import numpy as np
import matplotlib.pylab as plt
import glob

def fem_cal(ant):
    ''' Reads most recent HPOL and VPOL fem calibration files from
        current directory for the given antenna, plots the data
        and a quartic fit, and prints the coefficients of the fit.
    '''
    hfiles = glob.glob('Antenna '+str(ant)+' H*.txt')
    if len(hfiles) > 1: hfiles = np.sort(hfiles)
    f = open(hfiles[-1],'r')
    lines = f.readlines()
    f.close()
    hp = []
    hv = []
    for i,line in enumerate(lines):
        if i > 0:
            vals = np.array(line.strip().split()).astype('float')
            hp.append(vals[0])
            hv.append(vals[4])
    hp = np.array(hp)
    hv = np.array(hv)
    pfit = np.polyfit(np.log(hv),hp,4)
    lhv = np.linspace(np.log(hv).min(),np.log(hv).max(),100)
    plt.plot(np.log(hv),hp,'o')
    plt.plot(lhv,np.polyval(pfit,lhv))
    for i,p in enumerate(pfit[::-1]):
        print 'HPOL.c{:d} = {:10.7f}'.format(i,p)
    vfiles = glob.glob('Antenna '+str(ant)+' V*.txt')
    if len(vfiles) > 1: vfiles = np.sort(vfiles)
    f = open(vfiles[-1],'r')
    lines = f.readlines()
    f.close()
    vp = []
    vv = []
    for i,line in enumerate(lines):
        if i > 0:
            vals = np.array(line.strip().split()).astype('float')
            vp.append(vals[0])
            vv.append(vals[-1])
    vp = np.array(vp)
    vv = np.array(vv)
    pfit = np.polyfit(np.log(vv),vp,4)
    lvv = np.linspace(np.log(vv).min(),np.log(vv).max(),100)
    plt.plot(np.log(vv),vp,'*')
    plt.plot(lvv,np.polyval(pfit,lvv))
    plt.ylabel('Power [dBm]')
    plt.xlabel('Log Voltage [V]')
    plt.title('Antenna '+str(ant))
    for i,p in enumerate(pfit[::-1]):
        print 'VPOL.c{:d} = {:10.7f}'.format(i,p)
    plt.show()
