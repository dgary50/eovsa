#
# Routines for setting up scan header arrays involving channels,
# for the case of 52 bands of 325 MHz each.  This is adapted from
# chan_util_bc.py, which is the routine good for the 34-band case.
# 
# History:
#   2022-Mar-13  DG
#      Converted to an object from chan_util_52.py, so I can keep all of
#      its parameters together without having to recalculate on each
#      call to a function.  This version allows specification of a
#      "dwell" fseq file, which allows recording at high cadence (for flares).
#

import numpy as np

class Chan_Info():
    def __init__(self):
        from copy import copy
        self.ifbw = 400.    # Actual bandwidth 200 MHz design
        self.gifbw = 325.   # Usable bandwidth in MHz
        self.nschanx = 256  # Number of subchannels to drop at low end
        self.nschan = 4096  # Total number of subchannels
        self.fsx = 1.075    # Frequency in GHz assigned to first subchannel
        #  There are 3328 subchannels in 325 MHz => n_sci = 3328/fast_nsavg
        #  n_sci = (416, 208, 104, 52) for fast_nsavg = (8, 16, 32, 64)
        self.fast_nsavg = 16   # Number of subschannels to average in fast mode
        # Default number of subchannels to average in each of 52 IF bands
        self.std_nsavg = [110,151,175,208,256,277,302,332,369]+[416]*43
        self.nsavg = copy(self.std_nsavg)

    def fseq2nsavg(self,fseqfile=None):
        ''' A "DWELL" frequency sequence file is used to override the 
            default nsavg.  The named frequency sequence file must exist 
            in the ACC parm folder, and must start with the string 'dwell'
            (or 'DWELL').
            
            The DWELL sequence should have a list of "normal" bands and
            one band that has a much longer "dwell" time. A dwell band will
            have 416 science channels (each averaging over 8 subchannels).
            There can be at most 96 additional channels (max 512).
            
            Calling this without an fseqfile will reset to the default nsavg
        '''
        from copy import copy
        self.nsavg = copy(self.std_nsavg)  # Start with standard list
        if fseqfile:
            if fseqfile.upper().find('DWELL') != 0:
                # Not a DWELL file, so exit
                return
            # Attempt to ftp named file from ACC
            import urllib2
            userpass = 'admin:observer@'
            try:
                f = urllib2.urlopen('ftp://'+userpass+'acc.solar.pvt/parm/'+fseqfile,timeout=1)
                lines = f.readlines()
                f.close()
            except:
                print 'CHAN_UTIL_52: Could not read frequency sequence from ACC'
                print 'Default science channels unchanged.'
                return
            durs = np.array(lines[0][11:].strip().split(',')).astype(float)
            bands = np.array(lines[1][14:].strip().split(',')).astype(int)
            dwellband = np.argmax(durs)+1
            bands[np.where(bands == dwellband)] *= -1
            for i in range(52):
                if (i+1) == dwellband:
                    # Minimum averaging will yield 416 channels in this band
                    #self.nsavg[i] = 8
                    self.nsavg[i] = self.fast_nsavg  
                else:
                    if (i+1) in bands:
                        # Will use the standard science channels in this band
                        pass
                    else:
                        # This implies no science channels in this band
                        self.nsavg[i] = 10000
            if self.tot_scichan() > 512:
                print 'CHAN_UTIL_52: fseq file',fseqfile,'results in too many science channels.'
                print 'Reverting to standard list.'
                self.nsavg = copy(self.std_nsavg)  # Revert to standard list

    def tot_scichan(self, verbose=False):
        ''' Calculate and return the total number of science channels
        '''
        nscichan=[]
        for n in self.nsavg:
            df=self.ifbw/self.nschan
            nscichan.append(int(self.gifbw/(n*df)))
        if verbose: print nscichan
        return sum(nscichan)

    def chan_asmt(self, bnd):
        ''' Calculate and return the science channel assignments for each
            of the 4096 subchannels in the specified band.
        '''
        # Input band (bnd) ranges from 1-52, but the code below uses band,
        # which ranges from 0-51.
        band = bnd - 1
        if bnd < 1 or bnd > 52:
            # Error in band number provided
            return -1

        df = self.ifbw/self.nschan

        chasmt = [0]*self.nschanx

        # Create list of science channels for each band up to this one
        nscichan = []
        for n in self.nsavg[0:band+1]:
            nscichan.append(int(self.gifbw/(n*df)))
        # Sum number of science channels for all bands up to but
        # not including this one.
        nsci = 0
        if band is 0:
            pass
        else:
            for n in nscichan[0:band]:
                nsci += n
        for i in range(nscichan[band]):
            if self.nsavg[band] == self.fast_nsavg:
                # This is a "dwell" band, so mark its channels by adding 512
                chasmt += [i+nsci+1+512]*self.nsavg[band]
            else:
                chasmt += [i+nsci+1]*self.nsavg[band]
        nrest = 4096 - len(chasmt)
        chasmt += [0]*nrest
        return chasmt

    def start_freq(self, band):
        ''' Calculate and return the starting frequency for the given band.
        '''
        # Input band (bnd) ranges from 1-52
        if band < 1 or band > 52:
            # Error in band number provided
            return -1

        # Frequency assigned to the first spectral channel, in GHz
        fsx = 1.075
        df = self.ifbw/self.nschan

        sf = []

        # Create list of frequencies for each science channel in this band
        nscichan = int(self.gifbw/(self.nsavg[band-1]*df))
        for n in range(nscichan):
            sf.append(self.fsx + (band-1)*0.325 
                      + (self.nschanx + self.nsavg[band-1]*n)*df/1000.)
        return sf

    def sci_bw(self, band):
        ''' Calculate and return the science channel bandwidth for the given band.
        '''
        # Input band (bnd) ranges from 1-52
        if band < 1 or band > 52:
            # Error in band number provided
            return -1

        df = self.ifbw/self.nschan

        nscichan = int(self.gifbw/(self.nsavg[band-1]*df))
        scibw = [self.nsavg[band-1]*df/1000.]*nscichan
        return scibw

    def plt_chan(self):
        ''' This is just a diagnostic to check that the channel assignments
            are as expected.
        '''
        import matplotlib.pyplot as plt
        bands = np.arange(52) + 1
        plt.figure(figsize=(10,6),dpi=100)
        plt.axis((1,18,0,1))
        plt.xlabel('RF Frequency (GHz)')
        plt.ylabel('Normalization')
        plt.xlim(1,18)
        for band in bands:
            chasmt = self.chan_asmt(band)
            sf = np.array(self.start_freq(band))
            bw = np.array(self.sci_bw(band))
            ef = sf + bw  # end frequency
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
    ci = Chan_Info()   # Default channelization is fine for this
    bfreqs=[] #start frequency of a band
    efreqs=[] #end frequency of a band
    for b in range(52):
        bfreqs.append(ci.start_freq(b+1)[0])
        efreqs.append(ci.start_freq(b+1)[-1] + ci.sci_bw(b+1)[-1])
    if isinstance(fghz,list) or isinstance(fghz,np.ndarray):
        bds = []
        for fg in fghz:
            bd, = np.where((np.array(bfreqs) < fg) & (np.array(efreqs) > fg))
            if len(bd) == 1:
                bds.append(bd[0]+1)
            else:
                if fg > 0.0:  # Suppress warning if the frequency is zero
                    print '{0:f} GHz is not found in any band'.format(fg)
                bds.append(-1)
        return np.array(bds)
    else:
        bd, = np.where((np.array(bfreqs) < fghz) & (np.array(efreqs) > fghz))
        if len(bd) == 1:
            return bd[0]+1
        else:
            if fg > 0.0:  # Suppress warning if the frequency is zero
                print '{0:f} GHz is not found in any band'.format(fghz)
            return -1
