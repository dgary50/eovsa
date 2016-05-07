#
'''
   Routine to read 1-s capture files (taken on a geostationary satellite)
   and determine the equivalent delays from the measured phase slope on
   multiple baselines '''
#
# History:
#   2016-Mar-08  DG
#     First written, based on earlier version for prototype correlator
#   2016-Mar-16  DG
#     Added chans keyword, to allow specifying a restricted range of 
#     channels over which to determine the delay.
#   2016-Mar-24  DG
#     Added xy2rl() routine.
#   2016-Mar-26  DG
#     Complete rewrite of delay_centers() to make it more robust
#   2016-Mar-29  DG
#     Fixed a bug with antenna indexes in case of skipped antennas.
#   2016-Apr-08  DG
#     Added delay_centers2sql() routine, to write the delay center values
#     to the SQL table abin
#   2016-Apr-09  DG
#     Moved delay_centers2sql() to cal_header.py module, to join other,
#     similar routines.
#

import pcapture2 as p
import numpy as np
import matplotlib.pylab as plt
import urllib2
from util import Time
import get_sat_info as gs

def ant_str2list(ant_str):
    ant_list = []
    try:
        grps = ant_str.split()
        for grp in grps:
            antrange = grp[3:].split('-')
            if len(antrange) == 1:
                if antrange != '':
                    ant_list.append(int(antrange[0])-1)
            elif len(antrange) == 2:
                ant_list += range(int(antrange[0])-1,int(antrange[1]))
    except:
        print 'Error: cannot interpret ant_str',ant_str
        return None
    return np.array(ant_list)
                    
def delay_centers(filename,satname,ants='ant1-13',doplot=False):
    ''' Reads specified 1-s capture file (specified as a filename string) 
        and analyzes a specific range of channels in the 12-12.5 GHz band 
        (band 23).  First measures the phase slope and converts to delay 
        on each baseline, and then calculates the optimum delay needed to 
        correct for the phase slope.

        The keyword nants can be used to limit the solution to a smaller
        number of baselines.
        
        If doplot is True, plots an overview plot of the corrected phases,
        including the standard deviation of the fit as text in each box.
        
        The keyword chans allows adjusting which subchannels the delay is
        determined from.  It is sometimes good to restrict the channel
        range to avoid problems with multiple satellites at one HA.

        TODO: It then optionally reads the delay_centers.txt file from 
        the ACC and updates it.
    '''
    
    def best_fit(freq, data, idx):
        size = len(idx)*9
        loc_idx = np.rollaxis(np.array([idx-4,idx-3,idx-2,idx-1,idx,idx+1,idx+2,idx+3,idx+4]),1)
        local_idx = np.sort(loc_idx.reshape(size))
        p1 = np.unwrap(np.angle(data[local_idx]))
        p1.shape = (len(idx),9)
        p0 = np.zeros(p1.shape,'float')
        p2 = np.zeros(p1.shape,'float')
        slope = np.zeros(len(idx),'float')
        for i in range(len(idx)):
            slope[i] = np.polyfit(freq[loc_idx[i]],p1[i],1)[0]
            p0[i] = p1[i] - 2*np.pi*i
            p2[i] = p1[i] + 2*np.pi*i
        mslope = np.median(slope)
        p0 = p0.reshape(size)
        p1 = p1.reshape(size)
        p2 = p2.reshape(size)
        slp = [np.polyfit(freq[local_idx],p0,1)[0]]
        slp.append(np.polyfit(freq[local_idx],p1,1)[0])
        slp.append(np.polyfit(freq[local_idx],p2,1)[0])
        slp = np.array(slp)
        return slp[np.argmin(np.abs(mslope-slp))]
        
    # First see if we can read delay_centers.txt from the ACC (return if failure)
    userpass = 'admin:observer@'
    try:
        # Read delay center file from ACC
        dlafile = urllib2.urlopen('ftp://'+userpass+'acc.solar.pvt/parm/delay_centers.txt',timeout=1)
        lines = dlafile.readlines()
        dcenx = []
        dceny = []
        for line in lines:
            if line[0] != '#':
                # Skip comment lines and takes next 2 numbers as delay 
                # centers [ns] in X & Y
                dcenx.append(float(line.strip().split()[1]))
                dceny.append(float(line.strip().split()[2]))
    except:
        print Time().iso,'ACC connection for delay centers timed out.'
        return [None]*3
        
    print 'Successfully read delay_centers.txt from ACC.'
    # Get satellite frequencies and polarizations (this is to use channel centers for fitting,
    # but is not yet completed)
    sat, = gs.get_sat_info(names=[satname])
    freq = np.linspace(0,4095,4096)*0.4/4096 + 12.15
    freqmhz = (freq*10000. + 0.5).astype('int')/10.
    pol = sat['pollist']
    ridx, = np.where(pol == 'R')
    rfrq = sat['freqlist'][ridx]
    ridx = []
    for f in rfrq:
        try:
            ridx.append(np.where(f == freqmhz)[0][0])
        except:
            pass
    ridx = np.array(ridx)
    ridx = ridx[np.where(ridx > 2047)]  # Indexes to use for fitting phases
    ridx = ridx[:-2]  # Last point or two seems to be unreliable
    lidx, = np.where(pol == 'L')
    lfrq = sat['freqlist'][lidx]
    lidx = []
    for f in lfrq:
        try:
            lidx.append(np.where(f == freqmhz)[0][0])
        except:
            pass
    lidx = np.array(lidx)
    lidx = lidx[np.where(lidx > 2047)]  # Indexes to use for fitting phases
    lidx = lidx[:-1]  # Last point seems to be unreliable

    ant_list = ant_str2list(ants)
    nants = len(ant_list)
    nbl = nants*(nants-1)/2
    # Create M arrays for nants antennas, and c arrays for nbl baselines
    Mx = np.zeros((nbl,nants))
    My = np.zeros((nbl,nants))
    cx = np.zeros(nbl,dtype='float')
    cy = np.zeros(nbl,dtype='float')
    # Read 1-s capture file of raw packets
    out = p.rd_jspec(filename)
    print 'Successfully read file',filename
    # Get conversion from antenna pair (bl) to "data source" index
    bl2ord = p.bl_list()
    # Declare storage for delays in ns
    tau_xx = np.zeros((nants,nants),dtype='float')
    tau_yy = np.zeros((nants,nants),dtype='float')
    tau_xy = np.zeros((nants,nants),dtype='float')
    tau_yx = np.zeros((nants,nants),dtype='float')
    tau_ns = np.zeros((nants,nants),dtype='float')
    sigma = np.zeros((nants,nants),dtype='float')
    iblx = 0
    ibly = 0
    if doplot:
        f, ax = plt.subplots(nants,nants)
        f.subplots_adjust(hspace=0,wspace=0)
        for axrow in ax:
            for a in axrow:
                a.xaxis.set_visible(False)
                a.yaxis.set_visible(False)

    # Declare storage for auto-correlation (X vs. Y) delays in ns
    tau_a = np.zeros(nants,dtype='float')
    for i in range(nants):
        ai = ant_list[i]
        # Determine phase difference between X and Y feed using auto-correlation on each ant
        phi_ar = np.unwrap(np.angle(out['a'][ai,2,ridx].sum(1)))
        phi_al = np.unwrap(np.angle(out['a'][ai,2,lidx].sum(1)))
        tau_a[i] = (np.polyfit(freq[ridx],phi_ar,1)[0] + np.polyfit(freq[lidx],phi_al,1)[0])/(4*np.pi)
    for i in range(0,nants-1):
        ai = ant_list[i]
        for j in range(i+1,nants):
            aj = ant_list[j]
            # Form all possible phase measures, and check for consistency
            data = out['x'][bl2ord[ai,aj],0,].sum(1)
            tau_xx[i,j] = (best_fit(freq,data,ridx)+best_fit(freq,data,lidx))/(4*np.pi)
            data = out['x'][bl2ord[ai,aj],1,].sum(1)
            tau_yy[i,j] = (best_fit(freq,data,ridx)+best_fit(freq,data,lidx))/(4*np.pi)
            data = out['x'][bl2ord[ai,aj],2,].sum(1)
            tau_xy[i,j] = (best_fit(freq,data,ridx)+best_fit(freq,data,lidx))/(4*np.pi)
            data = out['x'][bl2ord[ai,aj],3,].sum(1)
            tau_yx[i,j] = (best_fit(freq,data,ridx)+best_fit(freq,data,lidx))/(4*np.pi)
            # Nominally, tau_xx and tau_yy are the same, and are what we want.  However,
            # the following two checks (residual slopes) should be small, so if either is 
            # greater than 1 then set both tau_xx and tau_yy to the same value, whichever
            # has the smaller slope residual.
            chk_xx = np.abs(tau_xx[i,j] - (tau_xy[i,j] - tau_a[j]))
            chk_yy = np.abs(tau_yy[i,j] - (tau_yx[i,j] + tau_a[j]))
            if chk_xx > chk_yy and chk_xx > 1:
                tau_xx[i,j] = tau_yy[i,j]
            if chk_yy > chk_xx and chk_yy > 1:
                tau_yy[i,j] = tau_xx[i,j]
            if doplot:
                phi0 = tau_xx[i,j]*2*np.pi*freq
                phix = np.angle(out['x'][bl2ord[ai,aj],0,].sum(1)) - phi0
                phix -= phix[2048]
                ax[i,j].plot(freq, (phix + np.pi) % (2*np.pi) - np.pi , '.')
                phi0 = tau_yy[i,j]*2*np.pi*freq
                phiy = np.angle(out['x'][bl2ord[ai,aj],1,].sum(1)) - phi0
                phiy -= phiy[2048]
                ax[j,i].plot(freq, (phiy + np.pi) % (2*np.pi) - np.pi, '.')
            tau_ns[i,j] = tau_xx[i,j]             # delay in nsec
            tau_ns[j,i] = tau_yy[i,j]             # delay in nsec
            Mx[iblx,np.array((i,j))] = -1,1
            cx[iblx] = dcenx[aj] - dcenx[ai] - tau_ns[i,j]
            iblx += 1
            My[ibly,np.array((i,j))] = -1,1
            cy[ibly] = dceny[aj] - dceny[ai] - tau_ns[j,i]
            ibly += 1

    if doplot:
        # Set axis scales and label the plots
        for i in range(nants):
            ai = ant_list[i]
            ax[i,i].text(0.5,0.5,str(ai+1),ha='center',va='center',transform=ax[i,i].transAxes,fontsize=14)
            for j in range(nants):
                ax[i,j].set_xlim(freq[0],freq[-1])
                ax[i,j].set_ylim(-np.pi,np.pi)

    print 'Successfully calculated delays from data.'
    # Obtain solution, only for those baselines with sigma < 1.0
    xdla,xr,_,_ = np.linalg.lstsq(Mx[:iblx],cx[:iblx])   # For X channel
    ydla,yr,_,_ = np.linalg.lstsq(My[:ibly],cy[:ibly])   # For Y channel
    print 'Successfully analyzed delays for mutual consistency.'
    # Calculate new delays for ants based on Ant 1 as reference
    newdlax = xdla - xdla[0] + dcenx[0]
    newdlay = ydla - ydla[0] + dceny[0] - tau_a

    f = open('/tmp/delay_centers_tmp.txt','w')
    for line in lines:
        if line[0] == '#':
            print line,
            f.write(line)
        else:
            ant = int(line[:6])
            idx, = np.where((ant-1) == ant_list)
            if len(idx) == 0:
                # No update for this antenna
                fmt = '{:4d}*{:12.3f} {:12.3f} {:12.3f} {:12.3f}'
                print fmt.format(ant,dcenx[ant-1],dceny[ant-1],dcenx[ant-1],dceny[ant-1])
                fmt = '{:4d} {:12.3f} {:12.3f} {:12.3f} {:12.3f}\n'
                f.write(fmt.format(ant,dcenx[ant-1],dceny[ant-1],dcenx[ant-1],dceny[ant-1]))
            else:
                idx = idx[0]
                fmt = '{:4d} {:12.3f} {:12.3f} {:12.3f} {:12.3f}'
                print fmt.format(ant,dcenx[ant-1],dceny[ant-1],newdlax[idx],newdlay[idx])
                fmt = '{:4d} {:12.3f} {:12.3f} {:12.3f} {:12.3f}\n'
                f.write(fmt.format(ant,dcenx[ant-1],dceny[ant-1],newdlax[idx],newdlay[idx]))
    f.close()
    fmt = '{:6.2f} '*nants
    for i in range(nants): print fmt.format(*tau_ns[i])
    return tau_ns, xr, yr

def dla_solve(nant=13):
    ''' This routine requests the user to select a file of delay measurements based
        on a set of geosynchronous satellite observations, at taken at different HA
        
        The measuremeents are in a text file, with each observation represented by
        a line with HA and Dec of the satellite, and an nant x nant matrix of delay
        errors tau_ij, where tau_ii = 0 is not used, and if i < j tau_ij = X poln,
        while if i > j tau_ij = Y poln.
    '''
    # Find and read measurement file
    # For now, hard-code the filename!
    filename = 'geosat_delays.txt'
    # Read the measurement file
    f = open(filename,'r')
    lines = f.readlines()
    nsat = len(lines)/(nant+1)
    isat = -1
    ha = np.zeros(nsat,dtype='float')
    dec = np.zeros(nsat,dtype='float')
    dla = np.zeros((nsat,nant,nant),dtype='float')
    for line in lines:
        if line.strip()[0] == '#':
            isat += 1
            j = -1
        else:
            if j == -1:
                # This is a "header" line with HA and Dec
                ha[isat], dec[isat] = np.array(line.strip().split()).astype('float')
                j += 1
            else:
                # This is the jth data line for sat isat
                dla[isat,j] = np.array(line.strip().split()).astype('float')
                j += 1
                
    # Speed of light in m/ns -- These units ensure that delay errors in ns give baseline errors in m
    c_light = 299792458./1.e9
    # Calculate the terms in delay equation
    # dtau_ij = (1/c)*[(dXj-dXi)*cos(dec)*cos(ha) - (dYj-dYi)*cos(dec)*sin(ha) + (dZj-dZi)*sin(dec) + dtau_cj - dtau_ci]
    #         = [(dXj - dXi)*cx + (dYj - dYi)*cy + (dZj - dZi)*cz + dtau_cj - dtau_ci]
    cx = cos(ha*pi/180.)*cos(dec*pi/180.)/c_light
    cy = -sin(ha*pi/180.)*cos(dec*pi/180.)/c_light
    #cz = sin(dec*pi/180.)/c_light
    #vals = np.concatenate((cx, cy, cz, np.array([1]*5)))
    vals = np.concatenate((cx, cy, np.array([1]*5)))
    vals.shape = (3,1,5)
    vals = rollaxis(rollaxis(vals,1),2)
    # Declare storage for matrix representing equations to be solved
    nbl = nant*(nant-1)/2
    #M = np.zeros((nsat,nbl,(nant-1)*4),dtype='float')
    M = np.zeros((nsat,nbl,(nant-1)*3),dtype='float')
    # Declare storage for column matrix of measurements
    cxx = np.zeros(nsat*nbl,dtype='float')
    cyy = np.zeros(nsat*nbl,dtype='float')
    k = 0    # Baseline counter
    # Loop over first antenna
    for i in range(0,nant-1):
        # Loop over second antenna
        for j in range(i+1,nant):
            # Set column vector of delays to the measurements
            cxx[k::nbl] = dla[:,i,j]
            cyy[k::nbl] = dla[:,j,i]
            if i != 0:
                # If antenna i is not the reference antenna...
                # Enter coefficient values with a minus sign
                M[:,k::nbl,(i-1)::(nant-1)] = -vals
            # For antenna j, enter coefficient values with plus sign                
            M[:,k::nbl,(j-1)::(nant-1)] = vals
            k += 1
    #M.shape = (nsat*nbl,(nant-1)*4)
    M.shape = (nsat*nbl,(nant-1)*3)
    sol_xx, xr, _, _ = np.linalg.lstsq(M,cxx)
    sol_yy, yr, _, _ = np.linalg.lstsq(M,cyy)
    sol_xy, xyr, _, _ = np.linalg.lstsq(np.concatenate((M,M)),np.concatenate((cxx,cyy)))
    
def dla_solve_enu(nant=13):
    ''' This routine requests the user to select a file of delay measurements based
        on a set of geosynchronous satellite observations, at taken at different HA
        
        The measuremeents are in a text file, with each observation represented by
        a line with HA and Dec of the satellite, and an nant x nant matrix of delay
        errors tau_ij, where tau_ii = 0 is not used, and if i < j tau_ij = X poln,
        while if i > j tau_ij = Y poln.
        
        This version solves for East-North-Up, with baselines in ns
    '''
    lat = 37.233170*numpy.pi/180       # OVSA Latitude (radians)

    # Find and read measurement file
    # For now, hard-code the filename!
    filename = '/common/tmp/geosat_delays.txt'
    # Read the measurement file
    f = open(filename,'r')
    lines = f.readlines()
    nsat = len(lines)/(nant+1)
    isat = -1
    ha = np.zeros(nsat,dtype='float')
    dec = np.zeros(nsat,dtype='float')
    dla = np.zeros((nsat,nant,nant),dtype='float')
    for line in lines:
        if line.strip()[0] == '#':
            isat += 1
            j = -1
        else:
            if j == -1:
                # This is a "header" line with HA and Dec
                ha[isat], dec[isat] = np.array(line.strip().split()).astype('float')*np.pi/180.
                j += 1
            else:
                # This is the jth data line for sat isat
                dla[isat,j] = np.array(line.strip().split()).astype('float')
                j += 1
                
    # Calculate the terms in delay equation
    # dtau_ij = [(dNj-dNi)*(cos(lat)*sin(dec) - sin(lat)*cos(dec)*cos(ha)) - (dEj-dEi)*cos(dec)*cos(ha) 
    #          + (dUj-dUi)*(cos(lat)*cos(dec)*cos(ha) + sin(lat)*sin(dec))]
    #         = [(dNj - dNi)*cx + (dEj - dEi)*cy + (dUj - dUi)*cz
    cx = cos(lat)*sin(dec) - sin(lat)*cos(ha)*cos(dec)
    cy = -sin(ha)*cos(dec)
    cz = cos(lat)*cos(dec)*cos(ha) + sin(lat)*sin(dec)
    #vals = np.concatenate((cx, cy, cz, np.array([1]*5)))
    vals = np.concatenate((cx, cy, cz))
    vals.shape = (3,1,5)
    vals = rollaxis(rollaxis(vals,1),2)
    # Declare storage for matrix representing equations to be solved
    nbl = nant*(nant-1)/2
    M = np.zeros((nsat,nbl,(nant-1)*3),dtype='float')
    # Declare storage for column matrix of measurements
    cxx = np.zeros(nsat*nbl,dtype='float')
    cyy = np.zeros(nsat*nbl,dtype='float')
    k = 0    # Baseline counter
    # Loop over first antenna
    for i in range(0,nant-1):
        # Loop over second antenna
        for j in range(i+1,nant):
            # Set column vector of delays to the measurements
            cxx[k::nbl] = dla[:,i,j]
            cyy[k::nbl] = dla[:,j,i]
            if i != 0:
                # If antenna i is not the reference antenna...
                # Enter coefficient values with a minus sign
                M[:,k::nbl,(i-1)::(nant-1)] = -vals
            # For antenna j, enter coefficient values with plus sign                
            M[:,k::nbl,(j-1)::(nant-1)] = vals
            k += 1
    M.shape = (nsat*nbl,(nant-1)*3)
    sol_xx, xr, _, _ = np.linalg.lstsq(M,cxx)
    sol_yy, yr, _, _ = np.linalg.lstsq(M,cyy)
    sol_xy, xyr, _, _ = np.linalg.lstsq(np.concatenate((M,M)),np.concatenate((cxx,cyy)))
    sol_xx.shape = (3,12)
    sol_yy.shape = (3,12)
    sol_xy.shape = (3,12)

def dla_solve_enutau(nant=13):
    ''' This routine requests the user to select a file of delay measurements based
        on a set of geosynchronous satellite observations, at taken at different HA
        
        The measuremeents are in a text file, with each observation represented by
        a line with HA and Dec of the satellite, and an nant x nant matrix of delay
        errors tau_ij, where tau_ii = 0 is not used, and if i < j tau_ij = X poln,
        while if i > j tau_ij = Y poln.
        
        This version solves for East-North-Up AND tau, with baselines in ns
    '''
    lat = 37.233170*numpy.pi/180       # OVSA Latitude (radians)

    # Find and read measurement file
    # For now, hard-code the filename!
    filename = '/common/tmp/geosat_delays.txt'
    # Read the measurement file
    f = open(filename,'r')
    lines = f.readlines()
    nsat = len(lines)/(nant+1)
    isat = -1
    ha = np.zeros(nsat,dtype='float')
    dec = np.zeros(nsat,dtype='float')
    dla = np.zeros((nsat,nant,nant),dtype='float')
    for line in lines:
        if line.strip()[0] == '#':
            isat += 1
            j = -1
        else:
            if j == -1:
                # This is a "header" line with HA and Dec
                ha[isat], dec[isat] = np.array(line.strip().split()).astype('float')*np.pi/180.
                j += 1
            else:
                # This is the jth data line for sat isat
                dla[isat,j] = np.array(line.strip().split()).astype('float')
                j += 1
                
    # Calculate the terms in delay equation
    # dtau_ij = [(dNj-dNi)*(cos(lat)*sin(dec) - sin(lat)*cos(dec)*cos(ha)) - (dEj-dEi)*cos(dec)*cos(ha) +
    #          + (dUj-dUi)*(cos(lat)*cos(dec)*cos(ha) + sin(lat)*sin(dec))] + dtau_cj - dtau_ci]
    #         = [(dNj - dNi)*cx + (dEj - dEi)*cy + (dUj - dUi)*cz
    cx = cos(lat)*sin(dec) - sin(lat)*cos(ha)*cos(dec)
    cy = -sin(ha)*cos(dec)
    cz = cos(lat)*cos(dec)*cos(ha) + sin(lat)*sin(dec)
    vals = np.concatenate((cx, cy, cz, np.array([1]*5)))
    #vals = np.concatenate((cx, cy, cz))
    vals.shape = (4,1,5)
    vals = rollaxis(rollaxis(vals,1),2)
    # Declare storage for matrix representing equations to be solved
    nbl = nant*(nant-1)/2
    M = np.zeros((nsat,nbl,(nant-1)*4),dtype='float')
    # Declare storage for column matrix of measurements
    cxx = np.zeros(nsat*nbl,dtype='float')
    cyy = np.zeros(nsat*nbl,dtype='float')
    k = 0    # Baseline counter
    # Loop over first antenna
    for i in range(0,nant-1):
        # Loop over second antenna
        for j in range(i+1,nant):
            # Set column vector of delays to the measurements
            cxx[k::nbl] = dla[:,i,j]
            cyy[k::nbl] = dla[:,j,i]
            if i != 0:
                # If antenna i is not the reference antenna...
                # Enter coefficient values with a minus sign
                M[:,k::nbl,(i-1)::(nant-1)] = -vals
            # For antenna j, enter coefficient values with plus sign                
            M[:,k::nbl,(j-1)::(nant-1)] = vals
            k += 1
    M.shape = (nsat*nbl,(nant-1)*4)
    sol_xx, xr, _, _ = np.linalg.lstsq(M,cxx)
    sol_yy, yr, _, _ = np.linalg.lstsq(M,cyy)
    sol_xy, xyr, _, _ = np.linalg.lstsq(np.concatenate((M,M)),np.concatenate((cxx,cyy)))
    sol_xx.shape = (4,12)
    sol_yy.shape = (4,12)
    sol_xy.shape = (4,12)

def dla_solve_xy(nant=13):
    ''' This routine requests the user to select a file of delay measurements based
        on a set of geosynchronous satellite observations, taken at different HA
        
        The measuremeents are in a text file, with each observation represented by
        a line with HA and Dec of the satellite, and an nant x nant matrix of delay
        errors tau_ij, where tau_ii = 0 is not used, and if i < j tau_ij = X poln,
        while if i > j tau_ij = Y poln.
        
        This version solves only for baseline errors in X and Y, in ns
    '''
    # Find and read measurement file
    # For now, hard-code the filename!
    filename = '/common/tmp/geosat_delays4.txt'
    # Read the measurement file
    f = open(filename,'r')
    lines = f.readlines()
    nsat = len(lines)/(nant+1)
    isat = -1
    ha = np.zeros(nsat,dtype='float')
    dec = np.zeros(nsat,dtype='float')
    dla = np.zeros((nsat,nant,nant),dtype='float')
    for line in lines:
        if line.strip()[0] == '#':
            isat += 1
            j = -1
        else:
            if j == -1:
                # This is a "header" line with HA and Dec
                ha[isat], dec[isat] = np.array(line.strip().split()).astype('float')
                j += 1
            else:
                # This is the jth data line for sat isat
                dla[isat,j] = np.array(line.strip().split()).astype('float')
                j += 1
                
    # Calculate the terms in delay equation (only X and Y terms)
    # dtau_ij = (1/c)*[(dXj-dXi)*cos(dec)*cos(ha) - (dYj-dYi)*cos(dec)*cos(ha)]
    #         = [(dXj - dXi)*cx - (dYj - dYi)*xy]
    cx = np.cos(ha*np.pi/180.)*np.cos(dec*np.pi/180.)
    cy = -np.sin(ha*np.pi/180.)*np.cos(dec*np.pi/180.)
    #cz = sin(dec*pi/180.)/c_light
    #vals = np.concatenate((cx, cy, cz, np.array([1]*5)))
    vals = np.concatenate((cx, cy))
    vals.shape = (2,1,5)
    vals = np.rollaxis(np.rollaxis(vals,1),2)
    # Declare storage for matrix representing equations to be solved
    nbl = nant*(nant-1)/2
    #M = np.zeros((nsat,nbl,(nant-1)*4),dtype='float')
    M = np.zeros((nsat,nbl,(nant-1)*2),dtype='float')
    # Declare storage for column matrix of measurements
    cxx = np.zeros(nsat*nbl,dtype='float')
    cyy = np.zeros(nsat*nbl,dtype='float')
    k = 0    # Baseline counter
    # Loop over first antenna
    for i in range(0,nant-1):
        # Loop over second antenna
        for j in range(i+1,nant):
            # Set column vector of delays to the measurements
            cxx[k::nbl] = dla[:,i,j]
            cyy[k::nbl] = dla[:,j,i]
            if i != 0:
                # If antenna i is not the reference antenna...
                # Enter coefficient values with a minus sign
                M[:,k::nbl,(i-1)::(nant-1)] = -vals
            # For antenna j, enter coefficient values with plus sign                
            M[:,k::nbl,(j-1)::(nant-1)] = vals
            k += 1
    #M.shape = (nsat*nbl,(nant-1)*4)
    M.shape = (nsat*nbl,(nant-1)*2)
    sol_xx, xr, _, _ = np.linalg.lstsq(M,cxx)
    sol_yy, yr, _, _ = np.linalg.lstsq(M,cyy)
    sol_xy, xyr, _, _ = np.linalg.lstsq(np.concatenate((M,M)),np.concatenate((cxx,cyy)))
    sol_xx.shape = (2,12)
    sol_yy.shape = (2,12)
    sol_xy.shape = (2,12)
    return sol_xx, sol_yy, sol_xy

def xy2rl(filename, satname):
    out = p.rd_jspec(filename)
    sat, = gs.get_sat_info([satname])
    freq = np.linspace(0,4095,4096)*0.4/4096 + 12.15
    frq = sat['freqlist']
    pol = sat['pollist']
    freqmhz = (freq*10000. + 0.5).astype('int')/10.
    ridx, = np.where(pol == 'R')
    lidx, = np.where(pol == 'L')
    rfrq = frq[ridx]
    ridx = []
    for f in rfrq:
        try:
            ridx.append(np.where(f == freqmhz)[0][0])
        except:
            pass
    ridx = np.array(ridx)

    nant = 13
    ipol = np.zeros((nant,4096),dtype='float')
    vpol = np.zeros((nant,4096),dtype='float')
    rpol = np.zeros((nant,4096),dtype='float')
    lpol = np.zeros((nant,4096),dtype='float')
    for k in range(nant):
        xx,yy,xy,yx = out['a'][k,:,:,30]

        pfitr = np.polyfit(freq[ridx],np.unwrap(np.angle(xy[ridx])),1)
        phi = np.polyval(pfitr,freq)
        print 'Ant',k+1,'xy slope:',pfitr[0],'Delay:',pfitr[0]/(2*np.pi),'ns'
        xyp = xy*(np.cos(phi+np.pi/2)-1j*np.sin(phi+np.pi/2))
        yxp = yx*(np.cos(phi+np.pi/2)+1j*np.sin(phi+np.pi/2))

        ipol[k] = np.abs(xx+yy)
        vpol[k] = np.imag(yxp - xyp)
        rpol[k] = np.real(xx+yy) + np.imag(yxp - xyp)
        lpol[k] = np.real(xx+yy) - np.imag(yxp - xyp)

    # Normalize to antenna 6 IPOL (kind of random)
    fac = ipol/ipol[5,:]
    f, ax = plt.subplots(4,1)
    for k in range(nant):
        ax[0].plot(freq,ipol[k]/fac[k])
        ax[1].plot(freq,vpol[k]/fac[k])
        ax[2].plot(freq,rpol[k]/fac[k])
        ax[3].plot(freq,lpol[k]/fac[k])

    ax[0].text(0.05,0.8,'Stokes I',transform=ax[0].transAxes,fontsize=12)
    ax[1].text(0.05,0.8,'Stokes V',transform=ax[1].transAxes,fontsize=12)
    ax[2].text(0.05,0.8,'RCP',transform=ax[2].transAxes,fontsize=12)
    ax[3].text(0.05,0.8,'LCP',transform=ax[3].transAxes,fontsize=12)

    for i in range(4):
        yrng = ax[i].yaxis.get_data_interval()
        ax[i].set_xlim(ax[i].xaxis.get_data_interval())
        for j in range(len(frq)):
            if pol[j] == 'R':
                ax[i].plot(frq[j]*np.ones(2)/1000.,yrng,color='red',linewidth=2)
            else:
                ax[i].plot(frq[j]*np.ones(2)/1000.,yrng,color='green',linewidth=2)
    ax[0].title(sat['name']+' Communication Satellite',fontsize=18)
    ax[3].set_xlabel('Frequency [GHz]')
    if sat['name'] == 'Ciel-2':
        ax[3].annotate('Beacon',(12.209,14000),xytext=(12.17,30000),arrowprops=dict(width=2,headwidth=6,frac=0.2,                facecolor='black', shrink=0.05))
    plt.show()
    return ipol,vpol,rpol,lpol
