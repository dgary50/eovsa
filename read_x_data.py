import aipy
import numpy as np

def readXdata(filename):
    ibl = np.array(
       [[0,0,1,2,0,0,0,0],
       [0,0,3,4,0,0,0,0],
       [0,0,0,5,0,0,0,0],
       [0,0,0,0,0,0,0,0],
       [0,0,0,0,0,6,7,8],
       [0,0,0,0,0,0,9,10],
       [0,0,3,4,0,0,0,11]])
    # Open uv file for reading
    uv = aipy.miriad.UV(filename)
    # Read one record to get number of good frequencies
    preamble, data = uv.read()
    uv.rewind()
    nf = len(data.nonzero()[0])
    freq = uv['sfreq'][data.nonzero()[0]]
    out = np.zeros((12,nf,600,2),dtype=np.complex64)
    l = -1
    tprev = 0
    for preamble, data in uv.all():
        uvw, t, (i,j) = preamble
        if uv['pol'] == -2:
            k = 1
        else:
            k = 0
        if t != tprev:
            # New time 
            l += 1
            tprev = t
            if l == 600:
                break
        if len(data.nonzero()[0]) == nf:
            out[ibl[i,j],:,l,k] = data[data.nonzero()]
    return out, freq
    
def get_X_data(files):
    #This takes a list of IDB files (as given by get_IDBfiles()) and concatenates
    #  the data they contain into a single IDBdata variable.
    IDBdatalist = []
    for file in files:
        #This will skip any files that give us errors.
        #  The names of the bad or unreadable files will
        #  be printed.
        try:
            IDBdata, freq = readXdata(file)
            IDBdatalist.append(IDBdata)
        except:
            print 'The problematic file is: ' + file
    IDBdata = concatenate(IDBdatalist,2)
    return IDBdata

import copy
import spectrogram_fit as sp
from util import Time

trange = Time(['2015-06-22 17:35:00','2015-06-22 18:35:00'])
s = sp.Spectrogram(trange)
s.fidx = [0,213]
tsys, std = s.get_median_data()

 
flaredata = ['/data1/IDB/IDB20150622173512',
 '/data1/IDB/IDB20150622174513',
 '/data1/IDB/IDB20150622175513',
 '/data1/IDB/IDB20150622180513',
 '/data1/IDB/IDB20150622181513',
 '/data1/IDB/IDB20150622182513']
 
IDBflaredata = get_X_data(flaredata)

outflare,fghzflare = readXdata('/data1/IDB/IDB20150622173512')
#this is just to dig out the frequencies, 'outflare' is not even used

pcal = angle(IDBflaredata[:,:,1050,:])
acal = abs(IDBflaredata[:,:,1050,:])
calout = copy.copy(IDBflaredata)
# Calibrate for time 1050, just before the initial peak of the flare.
for i in range(3600):
    calout[:,:,i,:] = calout[:,:,i,:]*(cos(pcal)-1j*sin(pcal))/acal

# Normalize to the total power spectrum at the same time, with reference
# to the shortest baseline (preserves the relative amplitudes on various
# baselines and polarizations.
norm = abs(calout[:,:,1050,:])
for i in range(12):
    for j in range(2):
        norm[i,:,j] = tsys[:,1050]*abs(calout[i,:,1050,j]) / abs(calout[0,:,1050,0])

for i in range(3600):
    calout[:,:,i,:] = calout[:,:,i,:]*norm

# Multi-panel Plot
# This looks good
f, ax = subplots(4,5)

sbl = ['1-2','1-3','1-4','2-3','2-4','3-4','5-7','5-8','7-8']
for i,ibl in enumerate([0,1,2,3,4,5,7,8,11]):
    if (i > 4):
        ax[2,i % 5].imshow(  calout2[ibl,50:,:,0])
        ax[2,i % 5].text(100,10,sbl[i]+' Amp',color='white')
        ax[3,i % 5].imshow(angle(calout2[ibl,50:,:,0]))
        ax[3,i % 5].text(100,10,sbl[i]+' Phase')
    else:
        ax[0,i % 5].imshow(  abs(calout2[ibl,50:,:,0]))
        ax[0,i % 5].text(100,10,sbl[i]+' Amp',color='white')
        ax[1,i % 5].imshow(angle(calout2[ibl,50:,:,0]))
        ax[1,i % 5].text(100,10,sbl[i]+' Phase')

ax[2,4].imshow(tsys[50:,:])
ax[2,4].text(100,10,'Total Power',color='white')
subplots_adjust(left=0.02, bottom=0.03, right=0.99, top=0.97, wspace=0.20, hspace=0.20)


# Saturation Plot
# This looks not so good with ants 7-8
# I re-did it with 1-2, it looks better
figure()
plot(tsys[100,:],abs(calout[0,100,:,0]),'.',label=str(fghzflare[100])[:5]+' GHz')
plot(tsys[150,:],abs(calout[0,150,:,0]),'.',label=str(fghzflare[150])[:5]+' GHz')
plot(tsys[200,:],abs(calout[0,200,:,0]),'.',label=str(fghzflare[200])[:5]+' GHz')
xlabel('Total Power [sfu]')
ylabel('Correlated Power (Ants 1-2) [sfu]')
legend(loc='lower right')
title('EOVSA Flare 2015-06-22')

for i in range(3000):
    for j in range(213):
        for k in range(12):
            for z in range(2):
                calout2[k, j, i, z] = abs(calout[k, j, i, z]) - abs(calout1[k, j, i, z])
                
# "uncalibrated" 
#f, ax = subplots(4,5)
#
#sbl = ['1-2','1-3','1-4','2-3','2-4','3-4','5-7','5-8','7-8']
#for i,ibl in enumerate([0,1,2,3,4,5,7,8,11]):
#    if (i > 4):
#        ax[2,i % 5].imshow(  abs(IDBflaredata[ibl,50:,:,0]))
#        ax[2,i % 5].text(100,10,sbl[i]+' Amp',color='white')
#        ax[3,i % 5].imshow(angle(IDBflaredata[ibl,50:,:,0]))
#        ax[3,i % 5].text(100,10,sbl[i]+' Phase')
#    else:
#        ax[0,i % 5].imshow(  abs(IDBflaredata[ibl,50:,:,0]))
#        ax[0,i % 5].text(100,10,sbl[i]+' Amp',color='white')
#        ax[1,i % 5].imshow(angle(IDBflaredata[ibl,50:,:,0]))
#        ax[1,i % 5].text(100,10,sbl[i]+' Phase')
#
#ax[2,4].imshow(tsys[50:,:])
#ax[2,4].text(100,10,'Total Power',color='white')
#subplots_adjust(left=0.02, bottom=0.03, right=0.99, top=0.97, wspace=0.20, hspace=0.20)
