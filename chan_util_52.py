#
# Routines for setting up scan header arrays involving channels,
# for the case of 52 bands of 325 MHz each.  This is adapted from
# chan_util_bc.py, which is the routine good for the 34-band case.
# 
# History:
#   2019-Feb-22  DG
#     First adapted from chan_util_bc.py
#   2019-Jun-26  DG
#      Fix bug in freq2bdnames() when unknown frequency was given.

import numpy as np
global nschan, ifbw, gifbw, nschanx, nsavg
# global parameters
# nschan: number of spectral channels
# ifbw: IF bandwidth in MHz
# gifbw: usable IF bandwidth in MHz, now excluding 50 MHz at the IF edge on one side (the other side is folded over for now) 
# nschanx: channels to skip at the beginning of each IF
# nsavg: number of spectral channels to average for each science channel

ifbw = 400.    # 200 MHz design
gifbw = 325.
nschanx = 256
nschan = 4096
nsavg = [110,151,175,208,256,277,302,332,369]+[416]*43

def tot_scichan():
    nscichan=[]
    for n in nsavg:
        df=ifbw/nschan
        nscichan.append(int(gifbw/(n*df)))
    print nscichan
    return sum(nscichan)

def chan_asmt(bnd):
    # Input band (bnd) ranges from 1-52, but the code below uses band,
    # which ranges from 0-51.
    band = bnd - 1
    if bnd < 1 or bnd > 52:
        # Error in band number provided
        return -1

    #global nschan, ifbw, gifbw, nschanx, nsavg
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
    # Input band (bnd) ranges from 1-52
    if band < 1 or band > 52:
        # Error in band number provided
        return -1

    global nschan, ifbw, gifbw, nschanx, nsavg
    # Frequency assigned to the first spectral channel, in GHz
    fsx = 1.075
    df = ifbw/nschan

    sf = []

    # Create list of frequencies for each science channel in this band
    nscichan = int(gifbw/(nsavg[band-1]*df))
    for n in range(nscichan):
        sf.append(fsx + (band-1)*0.325 + (nschanx + nsavg[band-1]*n)*df/1000.)

    return sf

def sci_bw(band):
    # Input band (bnd) ranges from 1-52
    if band < 1 or band > 52:
        # Error in band number provided
        return -1

    #global nschan, ifbw, gifbw, nschanx, nsavg
    df = ifbw/nschan

    nscichan = int(gifbw/(nsavg[band-1]*df))
    scibw = [nsavg[band-1]*df/1000.]*nscichan

    return scibw

def plt_chan():
    import matplotlib.pyplot as plt
    bands=np.arange(52)+1
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
#    import urllib2, copy
#    from util import Time
#    userpass = 'admin:observer@'
#    if t is None:
#        # Get current date
#        t = Time.now()
#    now = t.iso[:10].replace('-','')
#    f = urllib2.urlopen('ftp://'+userpass+'acc.solar.pvt/parm',timeout=0.5)
#    files = f.readlines()
#    f.close()
#    goodfile = ''
#    flist = []
#    # Find files that start with 'rfi'
#    for file in files:
#        fname = file.strip().split()[-1]
#        if fname.find('rfi') != -1:
#            flist.append(fname)
#    flist.sort()
#    for file in flist:
#        # This line starts with rfi, so interpret rest as a date
#        datstr = file[3:-4]
#        try:
#            idat = int(datstr)
#            if int(now) >= int(datstr):
#                # Current date is same or later than file date, so this is potentially the file we want
#                goodfile = file
#        except:
#            # Filename began with 'rfi', but date could not be converted to integer, so skip this file
#            pass
#    bandmask = np.ones((34,4096),'byte')  # start with no flagged channels
#    if goodfile != '':
#        # Got an appropriate file (last one prior to current date)
#        # Format of the file is BAND: List, where BAND is integer 1-34, followed by a colon ':', 
#        # and List is a comma-separated list of integers or ranges, e.g. 385, 861, 2945-2946, etc.
#        # Both bands and channel lists can be in any order.
#        RFIfile = urllib2.urlopen('ftp://'+userpass+'acc.solar.pvt/parm/'+goodfile,timeout=0.5)
#        lines = RFIfile.readlines()
#        for line in lines:
#            band, rest = line.split(':')
#            chanlist = rest.split(',')
#            for chaninfo in chanlist:
#                if chaninfo.find('-') == -1:
#                    # This is a single channel
#                    bandmask[int(band)-1,int(chaninfo)] = 0
#                else:
#                    # This is a range of channels
#                    chans = chaninfo.split('-')
#                    bandmask[int(band)-1,np.arange(int(chans[0]),int(chans[1])+1)] = 0
#    # We have bandmask, the channel mask as a function of band.  Now use the supplied fsequence
#    # to generate chanmask, the mask for the 1-s cycle.
#    chanmask = np.ones((50,4096),'byte')
#    bands = np.array(fsequence.split(',')).astype('int')-1
#    for i,band in enumerate(bands):
#        # Transfers the 4096 values for band into the ith slot of chanmask
#        chanmask[i,:] = copy.copy(bandmask[band,:])
#        if ifbw == 400.:
#            #chanmask[i,:512] = 0  # Temporary--flag only lower 50 MHz (test for new IF filters)
#            chanmask[i,:2148] = 0  # Temporary--flag all overlapped channels for 800 MHz clock
#    chanmask.shape = (204800)
#    # NB: This line overrides everything done above!  I.e. there are no masked channels.
#    if ifbw == 600.:
    chanmask = np.ones(204800,'byte')
    return chanmask

def freq2bdname(fghz):
    '''figure out the band name (1-52) from a given frequency in GHz (such as uv['sfreq'] values from Miriad 
       UV data, which are frequency values at the CENTER of the frequency channel). Return -1 if not found 
    '''
    bfreqs=[] #start frequency of a band
    efreqs=[] #end frequency of a band
    for b in range(52):
        bfreqs.append(start_freq(b+1)[0])
        efreqs.append(start_freq(b+1)[-1]+sci_bw(b+1)[-1])
    if isinstance(fghz,list) or isinstance(fghz,np.ndarray):
        bds = []
        for fg in fghz:
            bd, = np.where((np.array(bfreqs) < fg) & (np.array(efreqs) > fg))
            if len(bd) == 1:
                bds.append(bd[0]+1)
            else:
                print '{0:f} GHz is not found in any band'.format(fg)
                bds.append(-1)
        return np.array(bds)
    else:
        bd, = np.where((np.array(bfreqs) < fghz) & (np.array(efreqs) > fghz))
        if len(bd) == 1:
            return bd[0]+1
        else:
            print '{0:f} GHz is not found in any band'.format(fghz)
            return -1
