#
#  RoachCal routines to analyze standard set of observations done after
#  rebooting the ROACHes.
# 
#  The roachcal.scd schedule captures three sets of 1-s capture files,
#  named 
#     PRT<yyyymmddhhmmss>adc.dat   Blank sky observation with standard
#                                    observing sequence, solar.fsq,
#                                    ND-OFF, DCM attn set to 16 dB.
#                                    This is for generating the standard
#                                    DCM attenuation levels for each
#                                    antenna and IF channel.
#     PRT<yyyymmddhhmmss>ndon.dat  Blank sky observation with fixed
#                                    frequency (band23.fsq) and ND-ON.
#                                    This is for measuring relative X
#                                    vs. Y delay using autocorrelations.
#     PRT<yyyymmddhhmmss>ciel.dat  Observation on CIEL-2 Geosat, with
#                                    fixed frequency (band23.fsq).  This
#                                    is for measuring delays on each
#                                    antenna relative to Ant 1.

import pcapture2 as p
import urllib2
import numpy as np

def DCM_cal(filename=None,fseqfile='solar.fsq',dcmattn=16,update=False):

    if filename is None:
        return 'Must specify ADC packet capture filename, e.g. "/dppdata1/PRT/PRT<yyyymmddhhmmss>adc.dat"'

    adc = p.rd_jspec(filename)
    pwr = rollaxis(adc['phdr'],2)[:,:,:2]
    fseqfile = 'solar.fsq'
    userpass = 'admin:observer@'
    fseq_handle = urllib2.urlopen('ftp://'+userpass+'acc.solar.pvt/parm/'+fseqfile,timeout=0.5)
    lines = fseq_handle.readlines()
    fseq_handle.close()
    for line in lines:
        if line.find('LIST:SEQUENCE') != -1:
            line = line[14:]
            bandlist = np.array(map(int,line.split(',')))
    new_pwr = np.zeros((34,16,2))
    for i in range(34):
        idx, = np.where(bandlist-1 == i)
        if len(idx) > 0:
            new_pwr[i] = np.median(pwr[idx],0)
    new_pwr.shape = (34,32)
    
    orig_table = np.zeros((34,30)) + dcmattn
    orig_table[:,26:] = 0

    attn = np.log10(new_pwr[:,:-2]/1600.)*10.
    new_table = (np.clip(orig_table + attn,0,30)/2).astype(int)*2
    DCMlines = []
    DCMlines.append('#      Ant1  Ant2  Ant3  Ant4  Ant5  Ant6  Ant7  Ant8  Ant9 Ant10 Ant11 Ant12 Ant13 Ant14 Ant15')
    DCMlines.append('#      X  Y  X  Y  X  Y  X  Y  X  Y  X  Y  X  Y  X  Y  X  Y  X  Y  X  Y  X  Y  X  Y  X  Y  X  Y')
    DCMlines.append('#     ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- -----')
    for band in range(1,35):
        DCMlines.append('{:2} :  {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2}'.format(band,*new_table[band-1]))
    return DCMlines
