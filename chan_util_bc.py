#
# Routines for setting up scan header arrays involving channels.
# 
# History:
#   2014-Dec-30  DG
#     Started this history log.  Added the function get_chanmask(),
#     to implement RFI flagging--uses new rfi*.txt files on ACC in /parm
#     directory
#   2015-Feb-28  DG
#     Fix some bugs with reading rfi*.txt files on ACC
#   2015-May-29  DG
#      Converted from using datime() to using Time() based on astropy.
#   2015-Jun-16  DG
#      FTP to ACC now requires a username and password
#   2015-Jun-18  DG
#      Rather major change to eliminate "bad" channels in overlap region
#   2016-May-20  DG
#      Adjusted the "good IF bandwith" gifb to 370 MHz, and adjusted the
#      flagging for overlapped IF in get_chanmask() to 2064 in order to better
#      optimize our 5 science channels in each science band.
#   2016-May-21  DG
#      Looks like dppxmp does not like going past the end of the nominal
#      500 MHz band, so set gifb back to 350 MHz and 2148 in get_chanmask().
#      This results in only four good frequencies per band at higher freqs.
#

import pdb
import matplotlib.pyplot as plt
import numpy as np
global nschan, ifbw, gifbw, nschanx, nsavg
# global parameters
# nschan: number of spectral channels
# ifbw: IF bandwidth in MHz
# gifbw: usable IF bandwidth in MHz, now excluding 50 MHz at the IF edge on one side (the other side is folded over for now) 
# nschanx: channels to skip at the beginning of each IF
# nsavg: number of spectral channels to average for each science channel
nschan = 4096
ifbw = 400.
gifbw = 350.
nschanx = 0
#ifbw = 600.
#gifbw = 500.
#nschanx = 341
nsavg = [64,97,131,162,200,227,262]+[284]*27
#nsavg = [64,97,131,162,40,227,262]+[568]*27 #temporarily increase the number of science channels on band 5
#nsavg = [64,97,131,162,200,227,262]+[568]*15+[20]+[568]*21 #temporarily increase the number of science channels on band 23 

def update_nsavg():
    #update each value in nsavg in order to cover all the good if bandwidth
    new_nsavg=[]
    for i in range(34):
        n=nsavg[i]
        df=ifbw/nschan
        nscichan=int(gifbw/(n*df))
        new_nscichan=[]
        for j in range(3):
            new_nscichan.append(nscichan-j)
        new_ns=[]
        new_fgaps=[]
        for nn in new_nscichan:
            new_n=int(gifbw/nn/df)
            new_ns.append(new_n)
            new_scibw=new_n*df
            new_fgap=gifbw-new_scibw*nn
            new_fgaps.append(new_fgap)
        new_ns=np.array(new_ns)
        new_fgaps=np.array(new_fgaps)
        ind=np.argmin(new_fgaps)
        new_nsavg.append(new_ns[ind])
    return new_nsavg

nsavg=update_nsavg()

def tot_scichan():
    nscichan=[]
    for n in nsavg:
        df=ifbw/nschan
        nscichan.append(int(gifbw/(n*df)))
    print nscichan
    return sum(nscichan)

def chan_asmt(bnd):
    # Input band (bnd) ranges from 1-34, but the code below uses band,
    # which ranges from 0-33.
    band = bnd - 1
    if bnd < 1 or bnd > 34:
        # Error in band number provided
        return -1

    global nschan, ifbw, gifbw, nschanx, nsavg
    df = ifbw/nschan

    chasmt = [0]*nschanx

    # Create list of science channels for each band up to this one
    nscichan = []
    for n in nsavg[0:band+1]:
        nscichan.append(int(gifbw/(n*df)))
    # Sum number of science channels for all bands up to but
    # not including this one.
    nsci = 0
    if band is 0:
        pass
    else:
        for n in nscichan[0:band]:
            nsci += n
    for i in range(nscichan[band]):
        chasmt += [i+nsci+1]*nsavg[band]
    nrest = 4096 - len(chasmt)
    chasmt += [0]*nrest
    return chasmt

def start_freq(band):
    # Input band (bnd) ranges from 1-34
    if band < 1 or band > 34:
        # Error in band number provided
        return -1

    global nschan, ifbw, gifbw, nschanx, nsavg
    # Frequency assigned to the first spectral channel, in GHz
    fsx = 0.95+(600.-ifbw)/1000.
    df = ifbw/nschan

    sf = []

    # Create list of frequencies for each science channel in this band
    nscichan = int(gifbw/(nsavg[band-1]*df))
    for n in range(nscichan):
        sf.append(fsx + (band-1)*0.5 + (nschanx + nsavg[band-1]*n)*df/1000.)

    return sf

def sci_bw(band):
    # Input band (bnd) ranges from 1-34
    if band < 1 or band > 34:
        # Error in band number provided
        return -1

    global nschan, ifbw, gifbw, nschanx, nsavg
    df = ifbw/nschan

    nscichan = int(gifbw/(nsavg[band-1]*df))
    scibw = [nsavg[band-1]*df/1000.]*nscichan

    return scibw

def plt_chan():
    bands=np.arange(34)+1
    plt.figure(figsize=(10,6),dpi=100)
    plt.axis((1,18,0,1))
    plt.xlabel('RF Frequency (GHz)')
    plt.ylabel('Normalization')
    plt.xlim(1,18)
    for band in bands:
        chasmt=chan_asmt(band)
        sf=np.array(start_freq(band))
        bw=np.array(sci_bw(band))
        ef=sf+bw
        for i in range(len(sf)):
            plt.axvspan(sf[i],ef[i],ymin=0,ymax=1,alpha=0.5,color='b')
    plt.show()

def get_chanmask(fsequence,t=None):
    ''' Given a frequency tuning sequence (specifies which bands are used in
        each of the 50 time slots of the 1-s cycle), read the appropriate
        RFI-survey file and generate the channel mask for flagging.  Optionally
        provide a Time() object with the desired date.  Not sure why this
        would be wanted, but just in case we need to go back to older data...
        
        Returns a 204800-byte array (50*4096) with 1's where channels are kept,
        and 0's were channels are to be flagged.
    '''
    import urllib2, copy
    from util import Time
    userpass = 'admin:observer@'
    if t is None:
        # Get current date
        t = Time.now()
    now = t.iso[:10].replace('-','')
    f = urllib2.urlopen('ftp://'+userpass+'acc.solar.pvt/parm',timeout=0.5)
    files = f.readlines()
    f.close()
    goodfile = ''
    flist = []
    # Find files that start with 'rfi'
    for file in files:
        fname = file.strip().split()[-1]
        if fname.find('rfi') != -1:
            flist.append(fname)
    flist.sort()
    for file in flist:
        # This line starts with rfi, so interpret rest as a date
        datstr = file[3:-4]
        try:
            idat = int(datstr)
            if int(now) >= int(datstr):
                # Current date is same or later than file date, so this is potentially the file we want
                goodfile = file
        except:
            # Filename began with 'rfi', but date could not be converted to integer, so skip this file
            pass
    bandmask = np.ones((34,4096),'byte')  # start with no flagged channels
    if goodfile != '':
        # Got an appropriate file (last one prior to current date)
        # Format of the file is BAND: List, where BAND is integer 1-34, followed by a colon ':', 
        # and List is a comma-separated list of integers or ranges, e.g. 385, 861, 2945-2946, etc.
        # Both bands and channel lists can be in any order.
        RFIfile = urllib2.urlopen('ftp://'+userpass+'acc.solar.pvt/parm/'+goodfile,timeout=0.5)
        lines = RFIfile.readlines()
        for line in lines:
            band, rest = line.split(':')
            chanlist = rest.split(',')
            for chaninfo in chanlist:
                if chaninfo.find('-') == -1:
                    # This is a single channel
                    bandmask[int(band)-1,int(chaninfo)] = 0
                else:
                    # This is a range of channels
                    chans = chaninfo.split('-')
                    bandmask[int(band)-1,np.arange(int(chans[0]),int(chans[1])+1)] = 0
    # We have bandmask, the channel mask as a function of band.  Now use the supplied fsequence
    # to generate chanmask, the mask for the 1-s cycle.
    chanmask = np.ones((50,4096),'byte')
    bands = np.array(fsequence.split(',')).astype('int')-1
    for i,band in enumerate(bands):
        # Transfers the 4096 values for band into the ith slot of chanmask
        chanmask[i,:] = copy.copy(bandmask[band,:])
        chanmask[i,:2148] = 0  # Temporary--flag all overlapped channels for 800 MHz clock
    chanmask.shape = (204800)
    return chanmask
            
        
        
