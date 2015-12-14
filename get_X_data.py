
#June 20, 2015
#RG

#Updated June 28, 2015
#RG
#Now updated, this file has functions that can pull the vital information from 
#  IDBfiles, whether get_X_data(data) is given the filenames or the timerange.
#  There is now only one plotting function. It has the capability of plotting
#  angle or abs for all 12 baselines on 12 subplots, or plotting u vs v 
#  for all 12 baselines(+ and -) on one figure. 
#
#  2015-07-03  DG
#    Added check for source name, and stop reading when source name changes.
#    Also now returns times as Time() object array.
#  2015-07-04  DG
#    Throws out any "short" files that have too few frequencies in get_X_data().
#    Also fixed bug where out was not truncated in readXdata() to actual number of times.
#  2015-07-26  DG
#    Fixed an off-by-one error in times returned by readXdata().
# 

import aipy
import os
from util import Time
import glob
import numpy as np
import matplotlib.pyplot as plt
import spectrogram_fit as sp
import copy

def readXdata(filename):
    #This routine reads the vital information from a single IDBfile.
    #  It is used in get_X_data(data).
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
    if 'source' in uv.vartable:
        src = uv['source']
    uv.rewind()
    nf = len(data.nonzero()[0])
    freq = uv['sfreq'][data.nonzero()[0]]
    out = np.zeros((12,nf,600,2),dtype=np.complex64)
    uvwfile = []
    timearray = []
    l = -1
    tprev = 0
    for preamble, data in uv.all():
        uvw, t, (i,j) = preamble
        if j == 13: 
            j = 3 # Ant 14 in Ant 4 slot for now
        if j == 8:
            j = 1 # Ant 9 in Ant 2 slot for now
        if i == 8:
            i = 1 # Again with Ant 9
        if uv['pol'] == -2: 
            k = 1
            uvwfile.append(uvw)
        else:
            k = 0
        if len(data.nonzero()[0]) == nf:
            if t != tprev:
                # New time 
                l += 1
                if l == 600:
                    break
                tprev = t
                timearray.append(t)
            out[ibl[i,j],:,l,k] = data[data.nonzero()]
    return out[:,:,:l+1,:], np.array(uvwfile), freq, np.array(timearray), src

def show_IDBdatainfo(IDBdata, datatype=None, minchannel=0, maxchannel=None, polarization=0):
    #This function will show 12 subplot. It can show either phase or absolute power,
    #  depending on the inputs.
    #
    #To make a uvw plot, put uvw in for IDBdata, and None for datatype.
    #  The minchannel and maxchannel do not matter for a UVW plot.
    p = polarization
    if maxchannel == None:
        m =(minchannel +1)
    else: 
        m =maxchannel
    if datatype == np.angle:
        type = 'Phases'
    else:
        if datatype == abs:
            type = 'Abs'
    baselines = [' 1 and 2 ', ' 1 and 3 ', ' 1 and 4 ', ' 2 and 3 ', ' 2 and 4 ', ' 3 and 4 ', ' 5 and 6 ', ' 5 and 7 ', ' 5 and 8 ', ' 6 and 7 ', ' 6 and 8 ', ' 7 and 8 ']
    polarizations = [' x', ' y']
    if datatype == None:   
       uvw = IDBdata
       for i in range(uvw.shape[0]):
            plt.plot(-uvw[i][:, 0], -uvw[i][:, 1], '.')
            plt.suptitle('uv plot', fontsize=18, fontweight='bold', color='green')
    else:
       f, ax = plt.subplots(2, 6, figsize=(21,10))
       for i in range(12): 
           ax[i/6, i % 6].plot(np.transpose(datatype(IDBdata[i, p, minchannel:m, :])), '.')
           ax[i/6, i % 6].set_title('Antennas ' + baselines[i] + polarizations[p], color='blue')
           if datatype == np.angle:
               ax[i/6, i % 6].set_ylim([-np.pi, np.pi])
           else: pass
       plt.subplots_adjust(left=0.03, bottom=0.05, right=0.97)
       plt.suptitle(type + ' for the 12 Baselines', fontsize=18, fontweight='bold', color='red')
       

def get_IDBfiles(showthelast=10):
    #This will return the most recent IDB files saved to the 
    #  directory /data1/IDB/. They will be returned in a list
    #  with the format '/data1/IDB/IDByyyymmddhhmmss'
    #  We can adjust the number of files shown with 'showthelast' optional variable.
    #  It will not show the file being written currently
    s = showthelast
    IDBfiles = os.listdir('/data1/IDB/')
    IDBfiles.sort()
    IDBfiles = IDBfiles[-(s+2):-2]
    for i in range(len(IDBfiles)):
        IDBfiles[i] = '/data1/IDB/'+IDBfiles[i] 
    return IDBfiles
    
def get_X_data(data):
    #This takes a list of IDB files (as given by get_IDBfiles() or get_trange_files(trange) ) and concatenates
    #  the abs and angle data they contain into a single IDBdata variable. It also returns the uvw coordinates, the
    #  frequencies, and the time array for the files.
    if type(data) == Time:
        data = get_trange_files(data)
    else:
        pass
    IDBdatalist = []
    uvw = []
    times = []
    src = ''
    for files in data:
        #This will skip any files which give us errors.
        #  The names of the bad or unreadable files will
        #  be printed.
        try:
            IDBdata, uvwdata, freq, timearray, src0 = readXdata(files)
            if src == '':
                # This is the first file so set source name and frequency list
                src = copy.copy(src0)
                fghz = freq
            if src != src0:
                print 'Source name:',src0,'is different from initial source name:',src
                print 'Will stop reading files.'
                break
            IDBdatalist.append(IDBdata)
            uvw.append(uvwdata)
            times.append(timearray)
        except:
            print 'The problematic file is: ' + files
    # Sometimes an IDB file will be truncated early and have a number of frequencies
    # that is incompatible with the others, so make a list of the number of frequencies 
    # in each IDBdata and keep only those with the most common number.
    nfs = []
    # Find the maximum number of frequencies
    for i in range(len(IDBdatalist)):
        nbl,nf,nt,npol = IDBdatalist[i].shape
        nfs.append(nf)
    nfs = np.array(nfs)
    nfgood = np.median(nfs)
    badidx, = np.where(nfs != nfgood)
    if len(badidx) != 0 and len(nfs) < 3:
        print 'Files have different numbers of frequencies, and result is ambiguous.'
        print 'Will take the file with the higher number of frequencies'
        nfgood = nfs.max()
    # Make a list, nfbad, of entries with fewer than the max number of frequencies
    nfbad = []
    for i in range(len(nfs)):
        if nfs[i] != nfgood: nfbad.append(i)
    # This loop will run only if length of nfbad > 0
    for i in range(len(nfbad)):
        # Eliminate some entries in the list with too few frequencies
        IDBdatalist.pop(nfbad[i])
        times.pop(nfbad[i])
        uvw.pop(nfbad[i])
        print 'File',data[nfbad[i]],'eliminated due to too few frequencies.'
    IDBdata = np.concatenate(IDBdatalist,2)
    times = Time(np.array(np.concatenate(times)).astype('float'),format='jd')
    uvw = np.concatenate(uvw,0)
    IDBdata = np.swapaxes(np.swapaxes(IDBdata, 3, 1), 2, 3)
    #IDBdata is (baselines, polarizations, frequency, time)
    return IDBdata, uvw, fghz, times

def get_trange_files(trange):
    #Given a timerange, this routine will take all relevant IDBfiles from
    #  that time range, put them in a list, and return that list.
    #  This function is used in get_X_data(data).
    fstr = trange[0].iso
    folder = '/data1/IDB'
    files = glob.glob(folder+'/IDB'+fstr.replace('-','').split()[0]+'*')
    files.sort()
    mjd1, mjd2 = trange.mjd.astype('int')
    if mjd2 != mjd1:
        if (mjd2 - 1) != mjd1:
            usage('Second date must differ from first by at most 1 day')
        else:
            fstr2 = trange[1].iso
            files2 = glob.glob(folder+'/IDB'+fstr2.replace('-','').split()[0]+'*')
            files2.sort()
            files += files2

    def fname2mjd(filename):
        fstem = filename.split('/')[-1]
        fstr = fstem[3:7]+'-'+fstem[7:9]+'-'+fstem[9:11]+' '+fstem[11:13]+':'+fstem[13:15]+':'+fstem[15:17]
        t = Time(fstr)
        return t.mjd

    filelist = []
    for filename in files:
        mjd = fname2mjd(filename)
        if mjd >= trange[0].mjd and mjd < trange[1].mjd:
            filelist.append(filename)
    return filelist
    
def show_selfcalibrated(index, trange, plot='multipanel', antennas=0):
    #for index, choose a time when the flare starting, with a large slope. 
    #the options for plot are 'multipanel' , 'saturation' , and 'uncalibrated'
    #antennas controls which pair the saturation plot shows. 
  
    data = get_trange_files(trange)
    IDBdata, uvw, freq, times = get_X_data(data)
    
    s = sp.Spectrogram(trange)
    s.fidx = [0,IDBdata.shape[2]]
    tsys, std = s.get_median_data()
    
    pcal = np.angle(IDBdata[:,:, :, index])
    acal = abs(IDBdata[:,:, :, index])
    calout = copy.copy(IDBdata)
    # Calibrate for time 'index', just before the initial peak of the flare.
    for i in range(IDBdata.shape[3]):
        calout[:,:,:,i] = calout[:,:,:,i]*(np.cos(pcal)-1j*np.sin(pcal))/acal
    # Normalize to the total power spectrum at the same time, with reference
    # to the shortest baseline (preserves the relative amplitudes on various
    # baselines and polarizations.
    norm = abs(calout[:,:, :, index])
    for i in range(12):
        for j in range(2):
            norm[i,j,:] = tsys[:,index]*abs(calout[i,j,:,index]) / abs(calout[0,0,:,index])
    for i in range(IDBdata.shape[3]):
        calout[:,:,:,i] = calout[:,:,:,i]*norm       
    if plot == 'multipanel':
        # Multi-panel Plot
        f, ax = plt.subplots(4,5)
        sbl = ['1-2','1-3','1-4','2-3','2-4','3-4','5-7','5-8','7-8']
        for i,ibl in enumerate([0,1,2,3,4,5,7,8,11]):
            if (i > 4):
                ax[2,i % 5].imshow(abs(calout[ibl,0,50:,:]))
                ax[2,i % 5].text(100,10,sbl[i]+' Amp',color='white')
                ax[3,i % 5].imshow(np.angle(calout[ibl,0,50:,:]))
                ax[3,i % 5].text(100,10,sbl[i]+' Phase')
            else:
                ax[0,i % 5].imshow(abs(calout[ibl,0,50:,:]))
                ax[0,i % 5].text(100,10,sbl[i]+' Amp',color='white')
                ax[1,i % 5].imshow(np.angle(calout[ibl,0,50:,:]))
                ax[1,i % 5].text(100,10,sbl[i]+' Phase')
            ax[2,4].imshow(tsys[50:,:])
            ax[2,4].text(100,10,'Total Power',color='white')
            plt.subplots_adjust(left=0.02, bottom=0.03, right=0.99, top=0.97, wspace=0.20, hspace=0.20)
    else:
        if plot == 'saturation':
            # Saturation Plot
            ants_ = 'Ants ' + str(antennas+1) + '-' + str(antennas+2)
            plt.figure()
            plt.plot(tsys[IDBdata.shape[2]-110,:],abs(calout[antennas, 0, IDBdata.shape[2]-110,:]),'.',label=str(freq[IDBdata.shape[2]-110])[:5]+' GHz')
            plt.plot(tsys[IDBdata.shape[2]-60,:],abs(calout[antennas, 0, IDBdata.shape[2]-60,:]),'.',label=str(freq[IDBdata.shape[2]-60])[:5]+' GHz')
            plt.plot(tsys[IDBdata.shape[2]-10,:],abs(calout[antennas, 0, IDBdata.shape[2]-10,:]),'.',label=str(freq[IDBdata.shape[2]-10])[:5]+' GHz')
            plt.xlabel('Total Power [sfu]')
            plt.ylabel('Correlated Power (' + ants_ + ') [sfu]')
            plt.legend(loc='lower right')    
        else:
            if plot == 'uncalibrated':                          
                #"uncalibrated" 
                f, ax = plt.subplots(4,5)
                sbl = ['1-2','1-3','1-4','2-3','2-4','3-4','5-7','5-8','7-8']
                for i,ibl in enumerate([0,1,2,3,4,5,7,8,11]):
                    if (i > 4):
                        ax[2,i % 5].imshow(abs(IDBdata[ibl,0,50:,:]))
                        ax[2,i % 5].text(100,10,sbl[i]+' Amp',color='white')
                        ax[3,i % 5].imshow(np.angle(IDBdata[ibl,0,50:,:]))
                        ax[3,i % 5].text(100,10,sbl[i]+' Phase')
                    else:
                        ax[0,i % 5].imshow(abs(IDBdata[ibl,0,50:,:]))
                        ax[0,i % 5].text(100,10,sbl[i]+' Amp',color='white')
                        ax[1,i % 5].imshow(np.angle(IDBdata[ibl,0,50:,:]))
                        ax[1,i % 5].text(100,10,sbl[i]+' Phase')
                ax[2,4].imshow(tsys[50:,:])
                ax[2,4].text(100,10,'Total Power',color='white')
                plt.subplots_adjust(left=0.02, bottom=0.03, right=0.99, top=0.97, wspace=0.20, hspace=0.20)
            else:
                print 'please choose valid plot type'

    
