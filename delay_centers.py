#
'''
   Routine to capture packets while on satellite Galaxy 3C and analyze
   them to check delay centers.  Creates a plot of phase slopes, and
   a table showing original delay centers (from acc's delay_centers.txt
   file), and updated delay centers'''
#
# History:
#   2015-Jul-23  DG
#     Started this history log.  Fixed some errors in print-out of table,
#     and improved formatting of the table.
#   2015-Aug-23  DG
#     Changed to reflect new order of dimensions returned from 
#     pcapture.capture() routine.
#

import pcapture as p
import numpy as np
import matplotlib.pylab as plt
import urllib2
from util import Time

def delay_centers(overwrite=True):
    ''' Grabs raw packets while array is on the Galaxy 3C satellite,
        and analyzes a specific range of channels in the 11-11.5 GHz
        band.  First it measures the phase slope and converts to delay
        on each baseline, and then calculates the optimum delay 
        needed to correct for the phase slope.  

        TODO: It then optionally reads the delay_centers.txt file from 
        the ACC and updates it.
    '''
    # First see if we can read delay_centers.txt from the ACC (return if failure)
    userpass = 'admin:observer@'
    try:
        # Read delay center file from ACC
        dlafile = urllib2.urlopen('ftp://'+userpass+'acc.solar.pvt/parm/delay_centers.txt',timeout=1)
        lines = dlafile.readlines()
        delay_centers = []
        dceny = []
        for line in lines:
            if line[0] != '#':
                # Skip comment lines and take second number as delay center [ns]
                delay_centers.append(float(line.strip().split()[1]))
                dceny.append(float(line.strip().split()[2]))
    except:
        print Time().iso,'ACC connection for delay centers timed out.'
        return [None]*5
        
    # Capture 1 second of raw packets
    out = p.capture('galaxy',1,'X',overwrite=overwrite)
    # ch1,ch2 = 3600,3750
    ch1,ch2 = 2050,3800
    # Calculate the frequencies corresponding to selected channels ch1-ch2
    f = np.linspace(ch1,ch2-1,ch2-ch1)*0.4/4096 + 12.15
    # Determine which times in the capture contain valid (non-zero) data
    good, = np.where(np.abs(out[0,0,100]) != 0)
    outs = out[:,:,ch1:ch2,good]
    good2, = np.where(np.abs(out[6,0,100]) != 0)
    outs[6:] = out[6:,:,ch1:ch2,good2]
    nbl,npol,nf,ngood = outs.shape
    # NB: These settings hard-code the antenna assignments--probably needs
    # to be aware of actual antenna assignments
    bl = ['1-2','1-3','1-4','2-3','2-4','3-4','5-6','5-7','5-8','6-7','6-8','7-8']
    pol = ['X','Y']
    taus = []
    c = np.zeros((nbl,npol))
    for i in range(nbl):
        # Determine the two antennas for this baseline (will break if ants > 9 -- i.e. not a single digit)
        ant1 = int(bl[i][0])
        ant2 = int(bl[i][2])
        # For each baseline, determine phase slope by fitting a straight line
        # to the mean phase over the relevant channels.  Array pm are the
        # two parameters of the fit [slope,offset] in radians
        for j in range(npol):
            phi = np.unwrap(np.angle(np.mean(outs[i,j],1)))  # Mean over time axis
            pm = np.polyfit(f,phi,1)

            # plot the slopes and fits
            plt.plot(f,phi,'.')     # data points
            y = np.polyval(pm,f)
            plt.plot(f,y)           # fit
            chi_sq = np.sum((y - phi) ** 2)
       
            # Convert slopes to delay in clock steps and in nsec
            tau = int(pm[0]*0.8/(2*np.pi)+0.5)   # delay in clock steps
            tau_ns = pm[0]/(2*np.pi)             # delay in nsec
            plt.text(f[nf/2],y[nf/2],str(tau))
            print bl[i]+pol[j],'=',tau,'(in ns:',tau_ns,'), chi_sq =',chi_sq
            taus.append(tau_ns)
            if j == 0:
                c[i,j] = delay_centers[ant2-1] - delay_centers[ant1-1] - tau_ns
            else:
                c[i,j] = dceny[ant2-1] - dceny[ant1-1] - tau_ns
    taus = np.array(taus)
    
    # For now, analyze these in two sets of four
    M = [[-1,1,0,0],[-1,0,1,0],[-1,0,0,1],[0,-1,1,0],[0,-1,0,1],[0,0,-1,1]]
    x1,_,_,_ = np.linalg.lstsq(M,c[:6,0])   # Uses only X channel
    x2,_,_,_ = np.linalg.lstsq(M,c[6:,0])   # Uses only X channel
    y1,_,_,_ = np.linalg.lstsq(M,c[:6,1])  # Uses only Y channel
    y2,_,_,_ = np.linalg.lstsq(M,c[6:,1])   # Uses only Y channel
    
    # Calculate new delays for ants 1-4 based on Ant 1 as reference
    newdlax1 = x1 - x1[0] + delay_centers[0]
    newdlay1 = y1 - y1[0] + dceny[0]
    # Calculate new delays for ants 5-8 based on Ant 5 as reference
    newdlax2 = x2 - x2[0] + delay_centers[4]
    newdlay2 = y2 - y2[0] + dceny[4]
    
    for line in lines:
        if line[0] == '#':
            print line
        else:
            ant = int(line[:6])
            if ant < 5:
                print '{:4d} {:12.3f} {:12.3f} {:12.3f} {:12.3f}'.format(ant,delay_centers[ant-1],dceny[ant-1],newdlax1[ant-1],newdlay1[ant-1])
            elif ant < 9:
                print '{:4d} {:12.3f} {:12.3f} {:12.3f} {:12.3f}'.format(ant,delay_centers[ant-1],dceny[ant-1],newdlax2[ant-5],newdlay2[ant-5])
            else:
                print line[:-2]
    
    return taus,c,M,x1,x2
