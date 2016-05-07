'''
   Module for accessing, analyzing and plotting SOLPNTCAL data'''
#
# History:
#   2014-Nov-29  DG
#      Extensive changes to utilize SQL database.  Also made gausfit()
#      routine more general.
#   2014-Dec-02  DG
#      Completed new routines find_solpnt(), get_solpnt(), and renamed
#      the older routines to find_solpnt_log(), get_solpnt_log().  Also
#      modified process_solpnt() to work with the new routines.
#   2014-Dec-10  DG
#      Made some small adjustments to get analyze_all() working with the
#      new scheme (now it definitely does not work with the old!)
#   2014-Dec-11  DG
#      Added dmp_tsys() routine.
#   2014-Dec-13  DG
#      From beginning, convert RA offsets to "cross-Dec" by multiplying
#      by cos(dec).  Fix annoying bug when only one time is returned by 
#      find_solpnt().  Add dec0 and ra0 values to output of process_solpnt()
#   2014-Dec-15  DG
#      Changes to make the tsys part work with UDB files
#   2014-Dec-21  DG
#      Encountered weird bug that just showed up--sometimes the subarray1
#      value shows up as 0 at start of a scan--possibly something that Gelu's
#      recent change caused?  Anyway, I now check the first 100 values
#      looking for a non-zero value.
#   2015-Jan-22  DG
#      Changes required to work with any number of antennas
#   2015-Jan-25  DG
#      Fix slight glitches to allow to work with more antennas.
#   2015-Jan-27  DG
#      Provide error return in gausfit() in case of too many nans
#   2015-Mar-29  DG
#      Add v42/v45 stateframe table update.
#   2015-Mar-31  DG
#      Add v45/v46 stateframe table update.
#   2015-Apr-02  DG
#      Finally succeeded in making this version independent.
#   2015-Apr-05  DG
#      Major improvement--added common_val_idx() to find the common times
#      in proc and otp arrays.  This is critical to getting good results.
#      The change required adding a tstamps key in the values returned
#      by get_solpnt() and process_solpnt(), to be used in process_tsys().
#      Also removed the now defunct routines to read solpnt from log file.
#   2015-May-24  DG
#      Moved import of dbutil inside find_solpnt and get_solpnt, since other
#      functions are useful offline, where dbutil is unnecessary.
#   2015-May-29  DG
#      Converted from using datime() to using Time() based on astropy.  Also
#      changed tstamps in solpnt and otp dictionaries to float.
#   2015-May-31  DG
#      Fixed rd_tsys to account for day number in tsys *.txt files
#   2015-Jun-16  DG
#      FTP to ACC now requires a username and password
#   2015-Jun-27  DG
#      Changed the code to standardize on names, content, and order of indices
#      of outputs for otp.  Names ut_mjd, fghz, and tsys will be used, with units
#      hopefully made obvious.  Order of indices will be 
#          (nant/nbl, npol, nfreq, ntimes). 
#   2015-May-05  DG
#      Small change to process_tsys() to use only small antennas (max number is 13)
#

import stateframe, stateframedef, struct, os, urllib2, sys
from matplotlib.backends.backend_pdf import PdfPages
from scipy.optimize import leastsq

import numpy as np
import datetime
import matplotlib.pyplot as plt
import util

if sys.platform[:5] == 'linux':
    from coord_conv import dradec2dazel

def find_solpnt(t=None):
    ''' Makes an SQL query to find the SOLPNTCAL scans for the date given in
        the Time() object t.  A list of timestamps is returned, along with 
        the timestamp of the object provided.
    '''
    import dbutil
    if t is None: 
        # Get today's date
        t = util.Time.now()
    timestamp = int(t.lv)
    stimestamp = timestamp - (timestamp % 86400)  # Start of day
    etimestamp = stimestamp + 86399               # End of day
    # Open handle to SQL database
    cursor = dbutil.get_cursor()
    # Try to find a scan header with project SOLPNTCAL (only works after 2014 Nov. 30)
    verstr = dbutil.find_table_version(cursor,timestamp,True)
    if verstr is None:
        print 'No scan_header table found for given time.'
        return [], timestamp
    # First retrieve the Project from all scan headers for the day
    cursor.execute('select timestamp,Project from hV'+verstr+'_vD1 where timestamp between '+str(stimestamp)+' and '+str(etimestamp)+' order by timestamp')
    data = np.transpose(np.array(cursor.fetchall()))
    names = stateframedef.numpy.array(cursor.description)[:,0]
    cursor.close()
    if len(data) == 0:
        # No SOLPNTCAL found, so return empty list (and timestamp)
        return [], timestamp
    else:
        projdict = dict(zip(names,data))
        projdict['timestamp'] = projdict['timestamp'].astype('float')  # Convert timestamps from string to float
    good = np.where(projdict['Project'] == 'SOLPNTCAL')[0]
    if len(good) != 0:
        if len(good) == 1:
            return [projdict['timestamp'][good]], timestamp
        else:
            return projdict['timestamp'][good], timestamp
    else:
        return [], timestamp

def get_solpnt(t=None):
    ''' Get the SOLPNT data from the SQL database, occurring after 
        time given in the Time() object t.  If omitted, the first 
        SOLPNT scan for the current day is used (if it exists). '''
    import dbutil
    tstamps, timestamp = find_solpnt(t)
    # Find first SOLPNTCAL occurring after timestamp (time given by Time() object)
    if tstamps != []:
        print 'SOLPNTCAL scans were found at ',        
        for tstamp in tstamps:
            if type(tstamp) is np.ndarray:
                # Annoyingly necessary when only one time in tstamps
                tstamp = tstamp[0]
            t1 = util.Time(tstamp,format='lv')
            print t1.iso,';',
        print ' '
        good = np.where(tstamps >= timestamp)[0]
        # This is the timestamp of the first SOLPNTCAL scan after given time
        if good.shape[0] == 0: 
            stimestamp = tstamps[0]
        else:
            stimestamp = tstamps[good][0]
    else:
        print 'Warning: No SOLPNTCAL scan found, so interpreting given time as SOLPNTCAL time.' 
        stimestamp = timestamp

    # Grab 300 records after the start time
    cursor = dbutil.get_cursor()
    # Now version independent!
    verstr = dbutil.find_table_version(cursor,timestamp)
    if verstr is None:
        print 'No stateframe table found for the given time.'
        return {}
    solpntdict = dbutil.get_dbrecs(cursor,version=int(verstr),dimension=15,timestamp=stimestamp,nrecs=300)
    # Need dimension-1 data to get antennas in subarray -- Note: sometimes the antenna list
    # is zero (an unlikely value!) around the time of the start of a scan, so keep searching
    # first 100 records until non-zero:
    for i in range(100):
        blah = dbutil.get_dbrecs(cursor,version=int(verstr),dimension=1,timestamp=stimestamp+i,nrecs=1)
        if blah['LODM_Subarray1'][0] != 0:
            break
    cursor.close()
    sub1 = blah['LODM_Subarray1'][0]
    subarray1 = []
    antlist = []
    for i in range(16): 
        subarray1.append(sub1 & (1<<i) > 0)
        if subarray1[-1]:
            antlist.append(i)

#    print 'Antlist:',antlist
#    rao = np.zeros([15,300],dtype='int')
#    deco = np.zeros([15,300],dtype='int')
#    trk = np.zeros([15,300],dtype='bool')
#    hpol = np.zeros([15,300],dtype='float')
#    vpol = np.zeros([15,300],dtype='float')
#    ra = np.zeros(15,dtype='float')
#    dec = np.zeros(15,dtype='float')
    ra = (solpntdict['Ante_Cont_RAVirtualAxis'][:,antlist]*np.pi/10000./180.).astype('float')
    dec = (solpntdict['Ante_Cont_DecVirtualAxis'][:,antlist]*np.pi/10000./180.).astype('float')
    hpol = (solpntdict['Ante_Fron_FEM_HPol_Voltage'][:,antlist]).astype('float')
    vpol = (solpntdict['Ante_Fron_FEM_VPol_Voltage'][:,antlist]).astype('float')
    rao = (solpntdict['Ante_Cont_RAOffset'][:,antlist]).astype('float')
    deco = (solpntdict['Ante_Cont_DecOffset'][:,antlist]).astype('float')
    times = solpntdict['Timestamp'][:,0].astype('int64').astype('float')
    # Convert pointing info to track information
    outdict = stateframe.azel_from_sqldict(solpntdict)
    trk = np.logical_and(outdict['dAzimuth'][:,antlist]<0.0020,outdict['dElevation'][:,antlist]<0.0020)
    
    return {'Timestamp':stimestamp,'tstamps':times,'antlist':antlist,'trjfile':'SOLPNT.TRJ','ra':ra,'dec':dec,
             'rao':rao,'deco':deco,'trk':trk,'hpol':hpol,'vpol':vpol}
             
             
def process_solpnt(soldata,trj=None,antlist=None):
    # Read trajectory file from the ACC (unavoidably assumes it has not changed
    # since the observations were taken).
    if trj:
        try:
            trjfile = open(trj,'r')
        except:
            print 'Could not open TRJFILE',trj
            return {}
    else:
        userpass = 'admin:observer@'
        trjfile = urllib2.urlopen('ftp://'+userpass+'acc.solar.pvt/parm/'+soldata['trjfile'])
    trjlines = trjfile.readlines()
    trjfile.close()
    trjrao = []
    trjdeco = []
    for line in trjlines:
        rao, deco, junk = line.split(' ')
        trjrao.append(int(rao))
        trjdeco.append(int(deco))
    if antlist is None:
        antlist = soldata['antlist']
    # Mask array
    mask = np.zeros([300,len(antlist),len(trjlines)],dtype='bool')
    # Loop over each line of solpnt
    for i in range(300):
        # Loop over each RAO,DECO position and calculate "good data" mask for each 
        # desired position
        for j in range(len(trjlines)):
            # True if RAO and DECO matches this line
            m = np.logical_and(soldata['rao'][i,antlist] == trjrao[j],soldata['deco'][i,antlist] == trjdeco[j])
            # True if above AND is tracking
            mask[i,:,j] = np.logical_and(m,soldata['trk'][i,antlist])
    # Loop over each RAO,DECO position, and average all good HPOL and VPOL voltage data
    hpol = np.zeros([len(antlist),len(trjlines)],dtype='float')
    vpol = np.zeros([len(antlist),len(trjlines)],dtype='float')
    for j in range(len(trjlines)):
        hpol[:,j] = (soldata['hpol'][:,antlist]*mask[:,:,j]).sum(0)/mask[:,:,j].sum(0)
        vpol[:,j] = (soldata['vpol'][:,antlist]*mask[:,:,j]).sum(0)/mask[:,:,j].sum(0)
    # Convert the trajectory RA offsets to "cross-dec" (use of median assumes most antennas
    # are tracking the correct declination)
    ra0 = np.median(soldata['ra'])
    dec0 = np.median(soldata['dec'])
    return {'Timestamp':soldata['Timestamp'],'tstamps':soldata['tstamps'],'antlist':antlist,
            'ra':soldata['ra'][:,antlist],'dec':soldata['dec'][:,antlist],'ra0':ra0,'dec0':dec0,
            'rao':  np.array(trjrao[:13])*np.cos(dec0), 'deco':np.array(trjdeco[13:]), 
            'hrao': hpol[:,:13], 'vrao': vpol[:,:13],
            'hdeco':hpol[:,13:], 'vdeco':vpol[:,13:], 'mask':mask}

def gausfit(X,data):
    ''' Given a set of locations X and amplitude data at each location
        fit the data with a gaussian and return the fit parameters and
        x,y arrays of smooth, fitted data '''
    def peval(x,p):
        return p[0]*np.exp(-(x - p[1])**2/p[2]**2) + p[3]
    def residuals(p,y,x):
        A, x0, w, b = p
        return y - A*np.exp(-(x - x0)**2/w**2) - b

    # This dies horribly if there are nans in the data, so remove them first
    X = X[~np.isnan(data)]
    data = data[~np.isnan(data)]
    if len(X) < 4:
        # Do not attempt fit if too many nans--just return all zeros
        return [0.0,0.0,0.0,0.0],np.zeros(100),np.zeros(100)
    xmax = X.max()
    xmin = X.min()
    wid = (xmax-xmin)*0.1  # initial width set at a 10th of the range
    p0 = [np.max(data) - np.min(data), X[np.argmax(data)], wid, np.min(data)]

    if len(data) < len(p0):
        return np.zeros(4),[0,1],[0,1]

    plsq = leastsq(residuals, p0, args=(data, X))

    xstep = (xmax-xmin)*0.01  # step is 100th of range
    x = np.arange(xmin,xmax+xstep,xstep)
    y = peval(x,plsq[0])
    #plt.plot(x,peval(x, plsq[0]), X, data, 'o')
    return plsq[0],x,y

def fitall(proc,plot=False,pp=None):
    ''' Fits SOLPNT data for all antennas contained in 
        processed data proc, returned from process_solpnt()
        Prints results to terminal and optionally plots
        the gaussian fits.  It also returns the fitted offsets. 
        If pp is provided and is not None, the plots are
        written to a multipage PDF file.
    '''
    t = util.Time(proc['Timestamp'],format='lv')
    nant = len(proc['antlist'])

    if plot:
        # row and column sharing
        f, ax = plt.subplots(7, 2, sharex='col', sharey='row')
        f.set_size_inches(10,14,forward=True)
        f.suptitle('Fits for SOLPNT scan at '+t.iso+' UT',fontsize=18)
        ax[0,0].set_title('RA Offset [blue=HPol, red=VPol]')
        ax[0,1].set_title('Dec Offset [blue=HPol, red=VPol]')

    # Temporary statement to make this work with no more than 7 antennas.
    if nant > 7: nant = 7
    raoh =  np.zeros(nant,dtype='float')
    decoh = np.zeros(nant,dtype='float')
    raov =  np.zeros(nant,dtype='float')
    decov = np.zeros(nant,dtype='float')
    xeloh = np.zeros(nant,dtype='float')
    xelov = np.zeros(nant,dtype='float')
    eloh =  np.zeros(nant,dtype='float')
    elov =  np.zeros(nant,dtype='float')
    print 'SOLPNT solution for',t.iso
    print 'Ant       XELO (deg)             ELO (deg)'
    print '      HPol   VPol   Avg      HPol   VPol   Avg'
    print '---- ------ ------ ------   ------ ------ ------'
    for i in range(nant):
        ant = proc['antlist'][i]
        [A, hrao, w, b], x, y = gausfit(proc['rao'],proc['hrao'][i,:])
        if plot: ax[i,0].plot(x,y,proc['rao'],proc['hrao'][i,:],'o',label='Hpol')
        [A, vrao, w, b], x, y = gausfit(proc['rao'],proc['vrao'][i,:])
        if plot: 
            ax[i,0].plot(x,y,proc['rao'],proc['vrao'][i,:],'o',label='VPol')
            ax[i,0].text(0.05,0.8,'Ant '+str(ant+1),transform=ax[i,0].transAxes,fontsize=14)
            ax[i,0].text(0.05,0.65,'H = {:6.3f}'.format(hrao/10000.),transform=ax[i,0].transAxes)
            ax[i,0].text(0.05,0.5,'V = {:6.3f}'.format(vrao/10000.),transform=ax[i,0].transAxes)
        [A, hdeco, w, b], x, y = gausfit(proc['deco'],proc['hdeco'][i,:])
        if plot: ax[i,1].plot(x,y,proc['deco'],proc['hdeco'][i,:],'o',label='HPol')
        [A, vdeco, w, b], x, y = gausfit(proc['deco'],proc['vdeco'][i,:])
        if plot: 
            ax[i,1].plot(x,y,proc['deco'],proc['vdeco'][i,:],'o',label='VPol')
            ax[i,1].text(0.05,0.8,'Ant '+str(ant+1),transform=ax[i,1].transAxes,fontsize=14)
            ax[i,1].text(0.05,0.65,'H = {:6.3f}'.format(hdeco/10000.),transform=ax[i,1].transAxes)
            ax[i,1].text(0.05,0.5,'V = {:6.3f}'.format(vdeco/10000.),transform=ax[i,1].transAxes)
        if sys.platform[:5] == 'linux':
            # Convert rao, deco to azo, elo (returned in radians)
            cosdec = np.cos(proc['dec0'])
            hxelo, helo = dradec2dazel(proc['ra0'],proc['dec0'],
                                      t,hrao*np.pi/10000./180./cosdec,hdeco*np.pi/10000./180.)
            vxelo, velo = dradec2dazel(proc['ra0'],proc['dec0'],
                                      t,vrao*np.pi/10000./180./cosdec,vdeco*np.pi/10000./180.)
        else:
            hxelo, helo = hrao*np.pi/10000./180., hdeco*np.pi/10000./180.
            vxelo, velo = vrao*np.pi/10000./180., vdeco*np.pi/10000./180. 
        # Print table of XEL and EL offsets as degrees
        print ' {:2d} {:6.3f} {:6.3f} {:6.3f}   {:6.3f} {:6.3f} {:6.3f}'.format(
               proc['antlist'][i]+1, hxelo*180./np.pi, vxelo*180./np.pi,(hxelo + vxelo)*180./np.pi/2.,
                                     helo*180./np.pi,velo*180./np.pi,(helo+velo)*180./np.pi/2.)
        if plot: 
            ax[i,0].text(0.65,0.65,'Hxel = {:6.3f}'.format(hxelo*180./np.pi),transform=ax[i,0].transAxes)
            ax[i,0].text(0.65,0.5,'Vxel = {:6.3f}'.format(vxelo*180./np.pi),transform=ax[i,0].transAxes)
            ax[i,1].text(0.65,0.65,'Hel = {:6.3f}'.format(helo*180./np.pi),transform=ax[i,1].transAxes)
            ax[i,1].text(0.65,0.5,'Vel = {:6.3f}'.format(velo*180./np.pi),transform=ax[i,1].transAxes)
        raoh[i] = hrao/10000.  # degrees
        raov[i] = vrao/10000.  # degrees
        decoh[i] = hdeco/10000.  # degrees
        decov[i] = vdeco/10000.  # degrees
        xeloh[i] = hxelo*180./np.pi  # degrees
        xelov[i] = vxelo*180./np.pi  # degrees
        eloh[i] = helo*180./np.pi  # degrees
        elov[i] = velo*180./np.pi  # degrees

    if plot:
        ax[nant-1,0].set_xlabel('RA Offset [0.0001 deg]') 
        ax[nant-1,1].set_xlabel('Dec Offset [0.0001 deg]') 
        if pp:
            pp.savefig(f)
            plt.close()
        else:
            plt.show()
    return {'Timestamp':proc['Timestamp'],'antlist':proc['antlist'],
            'HPol':{'rao':raoh, 'deco':decoh, 'elo':eloh, 'xelo':xeloh},
            'VPol':{'rao':raov, 'deco':decov, 'elo':elov, 'xelo':xelov}}

def analyze_all(t,plot=False,pdffile=None):
    times, timestamp = find_solpnt(t)
    out = []
    pp = None
    antlist = range(4)
    if plot:
        if pdffile:
            pp = PdfPages(pdffile)
    if len(times) > 1:
        for tstamp in times:
            t = util.Time(tstamp,format='lv')
            soldata = get_solpnt(t)
            proc = process_solpnt(soldata,antlist=antlist)
            out.append(fitall(proc,plot,pp))#=False))
    else:
        t = util.Time(times[0],format='lv')
        soldata = get_solpnt(t)
        proc = process_solpnt(soldata,antlist=antlist)
        out.append(fitall(proc,plot,pp))#=False))
    if plot:
        axtype = [' RA, DEC',' XEL, EL']
        f, ax = plt.subplots(len(antlist), 2, sharex=True, sharey=True)
        f.set_size_inches(4,2*len(antlist),forward=True)
        t = util.Time(proc['Timestamp'],format='lv')
        f.suptitle('SOLPNT results for '+t.iso[:10],fontsize=18)
        for i in range(len(ax)):
            for j in range(len(ax[i])):
                ax[i][j].grid()
                ax[i][j].set_aspect('equal',adjustable='box-forced')
                ax[i][j].spines['right'].set_position('zero')
                ax[i][j].spines['top'].set_position('zero')
                ax[i][j].spines['left'].set_color('none')
                ax[i][j].spines['bottom'].set_color('none')
                ax[i][j].text(0.05,0.8,'Ant '+str(i+1)+axtype[j],transform=ax[i][j].transAxes,fontsize=10)
        for ant in antlist:
            for i in range(len(times)):
                if ant in out[i]['antlist']:
                    iant = out[i]['antlist'].index(ant)
                    x = [out[i]['HPol']['rao'][iant],out[i]['VPol']['rao'][iant]]
                    y = [out[i]['HPol']['deco'][iant],out[i]['VPol']['deco'][iant]]
                    ax[iant,0].plot(x,y,'.',x,y,'--')
                    x = [out[i]['HPol']['xelo'][iant],out[i]['VPol']['xelo'][iant]]
                    y = [out[i]['HPol']['elo'][iant],out[i]['VPol']['elo'][iant]]
                    ax[iant,1].plot(x,y,'.',x,y,'--')
        if pp:
            pp.savefig(f)
            plt.close()
            pp.close()
    return out

def dmp_tsys(t):
    ''' Given a start time as a Time() object t, find and dump the
        corresponding tsys values using Miriad (note that this only works
        on dpp or pipeline).  Returns the filenames of the files created.
    '''
    import socket, time
    import dump_tsys as dtsys

    hostname = socket.gethostname()
    if hostname != 'dpp' and hostname != 'pipeline':
        print 'Error: Cannot run dump_tsys() on host:',hostname
        return
    tnow = util.Time.now()
    t1 = t.lv - 10.   # 10 s before given time
    t2 = t.lv + 300.  # 300 s after given time
    trange = util.Time([t1,t2],format='lv')
    dtsys.dump_tsys(trange)
    # Check that the dump worked (wait 2 s for data)
    time.sleep(2)
    filestem ='t'+t.iso.replace('-','').replace(':','').replace(' ','')[:12]+'*.txt'
    xfiles = dtsys.glob.glob('/common/tmp/txt/x'+filestem)
    yfiles = dtsys.glob.glob('/common/tmp/txt/y'+filestem)
    return xfiles,yfiles
    
def rd_tsys(tsysfile,sfreqfile):
    ''' Read the Miriad log file output for a SOLPNT scan 
        To generate the Miriad log file, enter the commands:
           varplt vis=<infile> xaxis=time yaxis=xtsys log=<tsys_file>
           varplt vis=<infile> xaxis=freq yaxis=wfreq log=<sfreq_file>
        where <infile> is the IDB or UDB file, and <tsys_file> is the
        Miriad log file name (.txt).  Note that <sfreq_file> is a Miriad
        log file containing the starting frequencies, which is not
        yet implemented, but should be added to this routine.
        
        Major change to standardize output to ut_mdj, fghz, and order
        of indices of tsys as (nant, nf, ntimes)
    '''
    with open(sfreqfile,'r') as s:
        sheader = []
        for i in range(3):
            sheader.append(s.readline().strip())
        try:
            nfrq = int(sheader[-1].split('freq(')[1].split(')')[0])
        except:
            # This must be the format used by Jim McTiernan, so read another line
            sheader.append(s.readline().strip())
            try:
                nfrq = int(sheader[-1].split('freq(')[1].split(')')[0])
            except:
                print 'Could not successfully read sfreq file.'
                #return {}
        sfreq = np.zeros(nfrq,dtype='float')
        if (nfrq/4)*4 != nfrq:
            rlim = nfrq/4 + 1
        else:
            rlim = nfrq/4
        for i in range(rlim):
            line = s.readline()
            vals = line.split()[-4:]
            sfreq[i*4:min((i+1)*4,nfrq)] = vals
        nf = len(np.where(sfreq != 0)[0])
            
    with open(tsysfile,'r') as f:
        header = []
        for i in range(4):
            header.append(f.readline().strip())
        # Set Time() object to value read from header line 2 (header[1])
        dt = datetime.datetime.strptime(header[1].split(' ')[-1],'%y%b%d:%H:%M:%S.0')
        t = util.Time(dt,format='datetime')
        # Set base time in LabVIEW Timestamp format
        basetime = t.lv
        print 'Base Time is:',t.iso,'Timestamp:',basetime
        # Get number of columns from header line 4 (header[3])
        sline = header[3].split(',')
        nant = int(sline[0].split('(')[1])
        nfrq2 = int(sline[1].split(')')[0])
        if nfrq != nfrq2:
            print tsysfile,'frequencies are not consistent with those in',sfreqfile
            return {}
        # Get number of times in file (just number of lines - 4 header lines, 
        # divided by the number of frequencies, times nant/4
        f.seek(0)
        nlinespertime = (nant-1)/4 + 1
        nt = (sum(1 for line in f)-4)/nfrq2/nlinespertime
        print 'There are',nt,'times,',nfrq2,'frequencies, and',nant,'antennas in this file'
        # Get back to data start by rewind and skipping 4 header lines
        f.seek(0)
        for i in range(4): f.readline()
        out = np.zeros([nt,nfrq,nant],dtype='float')
        tstamp = np.zeros(nt,dtype='float')
        j = -1
        while 1:
            try:
                line = f.readline()
                vals = line.strip().split()
            except:
                if line.strip() != '':
                    print 'Error interpreting line:'
                    print 'Line is:',line.strip()+';'
                break
            # This assumes Miriad dumps no more than 4 ants/line
            nvals = len(vals)
            if nvals == min(nant,4)+2:
                hms = vals[1].split(':')
                secs = int(vals[0])*86400 + int(hms[0])*3600 + int(hms[1])*60 + int(hms[2])
                j += 1
                tstamp[j] = basetime+secs
                vals = vals[2:]
                nvals -= 2
                i = 0
                ntot = 0  # Number inserted into "out" so far
            try:
                out[j,i,ntot:ntot+nvals] = np.array(vals,dtype='float')
                ntot += nvals
                if ntot >= nant:
                    ntot = 0
                    i += 1
            except:
                if line.strip() != '':
                    print 'Error interpreting line at time,freq:',j,i
                    print 'Line is:',line.strip()+';'
                break
    otp = {'ut_mjd':util.Time(tstamp[:j+1],format='lv').mjd,'fghz':sfreq[:nf],'tsys':np.swapaxes(out[:j+1,:nf,:],0,2),'header':header}
    return otp

def process_tsys(otp, proc, pol=None):
    ''' OTP contains tsys, of size [nant, npol, nf, nt], ut_mjd of size [nt], fghz of size [nf]
             and a 4-line ascii header.  nant is 4 for prototype
        This only operates on one polarization at a time.  If npol = 2, use pol=0 and pol=1 in
            two subsequent calls to do both polarizations.
        PROC contains mask of size [nant', nt', npnt], and an array of timestamps,
             to be compared with OTP['ut'] to synchronize them.  nant' is the number of
             ants in the subarray.  This assumes the first four ants are the prototype ants.
    '''
    tsec1 = proc['tstamps']
    # The otp times, if from miriad, are slightly off from exact times, so round to the ms
    tsec2 = (util.Time(otp['ut_mjd'],format='mjd').lv + 0.001).astype('int64').astype('float')
    idx1, idx2 = common_val_idx(tsec1, tsec2)
    if pol is None:
        nant, nf, nt = otp['tsys'].shape
    else:
        nant, npol, nf, nt = otp['tsys'].shape
    nant = min(nant,13) # Make sure only small dishes are used
    # nf = len(np.where(otp['sfreq'] != 0)[0])  # Override nf with number of non-zero frequencies
    m = np.transpose(proc['mask'],(1,0,2))[0:nant,idx1,:]  # mask now same number of ants and times as tsys
    print m.shape
    npt = len(m[0,0,:])
    if pol is None:
        tsys = otp['tsys']
    else:
        tsys = otp['tsys'][:,pol,:,:]
    tsys = tsys[:,:,idx2]
    hpol = np.zeros([nf,npt,nant],'float')
    pra = np.zeros([4,nf,nant],'float')  # Array of gaussian fits
    pdec = np.zeros([4,nf,nant],'float')  # Array of gaussian fits
    # Step through pointings, applying mask
    for ipnt in range(npt):
        for iant in range(nant):
            good = np.where(m[iant,:,ipnt])[0]
            if len(good) == 0:
                # This pointing location had no good tracking, so set values to Nan
                hpol[:,ipnt,iant] = np.nan
            else:
                # Set values to median of good-tracking period
                hpol[:,ipnt,iant] = np.median(tsys[iant,:,good],0)
    # At this stage, hpol is of size [nf, npnt, nant]
    rao = hpol[:,:13,:]   # First 13 pointings are RAO
    deco = hpol[:,13:,:]  # Second 13 pointings are DECO
    for iant in range(nant):
        for ifreq in range(nf):
            pra[:,ifreq,iant], x, y = gausfit(proc['rao'],rao[ifreq,:,iant])
            pdec[:,ifreq,iant], x, y = gausfit(proc['deco'],deco[ifreq,:,iant])
    return pra, pdec, rao, deco
        
def common_val_idx(array1,array2):
    ''' Find the common values in two sorted arrays, and return the array
        of indexes of those common values in the two arrays.
    '''
    common = np.intersect1d(array1,array2,True)
    idx1 = np.searchsorted(array1,common)
    idx2 = np.searchsorted(array2,common)
    return idx1, idx2
    
    # To calculate total power calibration factors for an antenna, just do
    # import solpnt, rstn
    # Get otp and proc dicts by appropriate means...
    # pra, pdec = solpnt.process_tsys(otp, proc)
    # t = util.Time(otp['ut'][0],format='lv')
    # ant0 = pra[0,:,0]
    # fmhz = otp['sfreq']*1000.
    # frq, flux = rstn.rd_rstnflux(t)
    # s = rstn.rstn2ant(frq, flux, fmhz, t)
    # calfac0 = s/ant0
    #
    # To apply to other solar data taken around the same time
    # sun = solpnt.rd_tsys(tsysfile, sfreqfile)
    # calsun0 = sun['tsys'][:,:,0]
    # for i in range(nt):
    #    calsun0[i,:] = (sun['tsys'][i,:,0]-pra[3,:,0])*calfac0
    
    # The tsys values for the UDB files are in fits files available on
    # the web site.  To read these in IDL, use mrdfits.  The relevant data
    # are in the third extension table.  To read them in Python (assumes
    # you are using Astropy), do 
    #    from astropy.io import fits
    #    hdulist = fits.open(fitsfile)
    #    hdulist.info()  # Lists extension tables--we want XTSYS
    #    xtsys = hdulist[2].data[0][1]  # array of shape (nt,nf,nant)
