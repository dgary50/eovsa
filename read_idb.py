''' Main routine to read the auto- and cross-correlation data from IDB files.
'''
#
#  2016-03-30  DG
#    Began writing the code, adapted from get_X_data2.py
#  2016-04-01  DG
#    Added summary_plot()

import aipy
import os
from util import Time
import glob
import numpy as np
import matplotlib.pyplot as plt
import spectrogram_fit as sp
import pcapture2 as p
import copy

def readXdata(filename,filter=True):
    #This routine reads the data from a single IDBfile.
    
    # the IDB array sorts basline pairs into correct order for an output
    # array that just has 18 baselines (2 sets of 4 antennas correlated separately
    # now just need to convert i,j into slot in 136-slot array: 15*16/2 corrs +16 auto
    # this array corresponds to the first index (auto corr) for each antenna
    # 16*(iant-1)-(iant-1)(iant-2)/2
    ibl = np.array( [ 0,16,31,45,58,70,81,91,100,108,115,121,126,130,133,135 ])

    # Open uv file for reading
    uv = aipy.miriad.UV(filename)
    good_idx = np.arange(len(uv['sfreq']))
    if filter:
        good_idx = []
        # Read a bunch of records to get number of good frequencies, i.e. those with at least 
        # some non-zero data.  Read 20 records for baseline 1-2, XX pol
        uv.select('antennae',0,2,include=True)
        uv.select('polarization',-5,-5,include=True)
        for i in range(20):
            preamble, data = uv.read()
            idx, = data.nonzero()
            if len(idx) > len(good_idx):
                good_idx = copy.copy(idx)
        uv.select('clear',0,0)
        uv.rewind()

    if 'source' in uv.vartable:
        src = uv['source']
    nf = len(good_idx)
    freq = uv['sfreq'][good_idx]
    npol = uv['npol']
    nants = uv['nants']
    nbl = nants*(nants-1)/2
    bl2ord = p.bl_list(nants)
    outa = np.zeros((nants,npol,nf,600),dtype=np.complex64)  # Auto-correlations
    outx = np.zeros((nbl,npol,nf,600),dtype=np.complex64)  # Cross-correlations
    outp = np.zeros((nants,2,nf,600),dtype=np.float)
    outp2 = np.zeros((nants,2,nf,600),dtype=np.float)
    outm = np.zeros((nants,2,nf,600),dtype=np.int)
    uvwarray = []
    timearray = []
    l = -1
    tprev = 0
    # Use antennalist if available
    if 'antlist' in uv.vartable:
        ants = uv['antlist']
        antlist = map(int, ants.split())
    else:
        antlist = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]

    for preamble, data in uv.all():
        uvw, t, (i0,j0) = preamble
        i = antlist.index(i0+1)
        j = antlist.index(j0+1)
        if i > j:
            # Reverse order of indices
            j = antlist.index(i0+1)
            i = antlist.index(j0+1)
        # Assumes uv['pol'] is one of -5, -6, -7, -8
        k = -5 - uv['pol']
        if k == 3:
            uvwarray.append(uvw)
        if filter:
            if len(data.nonzero()[0]) == nf:
                if t != tprev:
                    # New time 
                    l += 1
                    if l == 600:
                        break
                    tprev = t
                    timearray.append(t)
                    xdata = uv['xsampler'].reshape(500,nants,3)
                    ydata = uv['ysampler'].reshape(500,nants,3)
                    outp[:,0,:,l] = np.swapaxes(xdata[good_idx,:,0],0,1)
                    outp[:,1,:,l] = np.swapaxes(ydata[good_idx,:,0],0,1)
                    outp2[:,0,:,l] = np.swapaxes(xdata[good_idx,:,1],0,1)
                    outp2[:,1,:,l] = np.swapaxes(ydata[good_idx,:,1],0,1)
                    outm[:,0,:,l] = np.swapaxes(xdata[good_idx,:,2],0,1)
                    outm[:,1,:,l] = np.swapaxes(ydata[good_idx,:,2],0,1)
    
                if i0 == j0:
                    # This is an auto-correlation
                    outa[i0,k,:,l] = data[data.nonzero()]
                else:
                    outx[bl2ord[i,j],k,:,l] = data[data.nonzero()]
        else:
            if t != tprev:
                # New time 
                l += 1
                if l == 600:
                    break
                tprev = t
                timearray.append(t)
                xdata = uv['xsampler'].reshape(500,nants,3)
                ydata = uv['ysampler'].reshape(500,nants,3)
                outp[:,0,:,l] = np.swapaxes(xdata[:,:,0],0,1)
                outp[:,1,:,l] = np.swapaxes(ydata[:,:,0],0,1)
                outp2[:,0,:,l] = np.swapaxes(xdata[:,:,1],0,1)
                outp2[:,1,:,l] = np.swapaxes(ydata[:,:,1],0,1)
                outm[:,0,:,l] = np.swapaxes(xdata[:,:,2],0,1)
                outm[:,1,:,l] = np.swapaxes(ydata[:,:,2],0,1)
            if i0 == j0:
                # This is an auto-correlation
                outa[i0,k,:,l] = data
            else:
                outx[bl2ord[i,j],k,:,l] = data

    out = {'a':outa, 'x':outx, 'uvw':np.array(uvwarray), 'fghz':freq, 'time':np.array(timearray),'source':src,'p':outp,'p2':outp2,'m':outm}
    return out

def readXdatmp(filename):
    # This temporary routine reads the data from a single IDBfile where the
    # polarization was written incorrectly. (-1,-2,0,0) instead of (-5,-6,-7,-8)
    
    # the IDB array sorts basline pairs into correct order for an output
    # array that just has 18 baselines (2 sets of 4 antennas correlated separately
    # now just need to convert i,j into slot in 136-slot array: 15*16/2 corrs +16 auto
    # this array corresponds to the first index (auto corr) for each antenna
    # 16*(iant-1)-(iant-1)(iant-2)/2
    ibl = np.array( [ 0,16,31,45,58,70,81,91,100,108,115,121,126,130,133,135 ])

    # Open uv file for reading
    uv = aipy.miriad.UV(filename)
    good_idx = []
    # Read a bunch of records to get number of good frequencies, i.e. those with at least 
    # some non-zero data.  Read 20 records for baseline 1-2, XX pol
    uv.select('antennae',0,2,include=True)
    uv.select('polarization',-1,-1,include=True)
    for i in range(20):
        preamble, data = uv.read()
        idx, = data.nonzero()
        if len(idx) > len(good_idx):
            good_idx = copy.copy(idx)
    if 'source' in uv.vartable:
        src = uv['source']
    uv.select('clear',0,0)
    uv.rewind()
    nf = len(good_idx)
    freq = uv['sfreq'][good_idx]
    npol = uv['npol']
    nants = uv['nants']
    nbl = nants*(nants-1)/2
    bl2ord = p.bl_list(nants)
    outa = np.zeros((nants,npol,nf,600),dtype=np.complex64)  # Auto-correlations
    outx = np.zeros((nbl,npol,nf,600),dtype=np.complex64)  # Cross-correlations
    outp = np.zeros((nants,2,nf,600),dtype=np.float)
    outp2 = np.zeros((nants,2,nf,600),dtype=np.float)
    outm = np.zeros((nants,2,nf,600),dtype=np.int)
    uvwarray = []
    timearray = []
    l = -1
    tprev = 0
    # Use antennalist if available
    if 'antlist' in uv.vartable:
        ants = uv['antlist']
        antlist = map(int, ants.split())
    else:
        antlist = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]
    kp = 0
    for preamble, data in uv.all():
        uvw, t, (i0,j0) = preamble
        i = antlist.index(i0+1)
        j = antlist.index(j0+1)
        if i > j:
            # Reverse order of indices
            j = antlist.index(i0+1)
            i = antlist.index(j0+1)
        # Assumes uv['pol'] is -1, -2, 0, 0
        if i == 0 and j == 0:
            if uv['pol'] == -1:
                k = 0
            elif uv['pol'] == -2:
                k = 1
            elif uv['pol'] == 0:
                if kp == 1: 
                    k = 2
                    kp = 0
                else:
                    k = 3
                    kp = 1
        if k == 3:
            uvwarray.append(uvw)
        if len(data.nonzero()[0]) == nf:
            if t != tprev:
                # New time 
                l += 1
                if l == 600:
                    break
                tprev = t
                timearray.append(t)
                xdata = uv['xsampler'].reshape(500,nants,3)
                ydata = uv['ysampler'].reshape(500,nants,3)
                outp[:,0,:,l] = np.swapaxes(xdata[good_idx,:,0],0,1)
                outp[:,1,:,l] = np.swapaxes(ydata[good_idx,:,0],0,1)
                outp2[:,0,:,l] = np.swapaxes(xdata[good_idx,:,1],0,1)
                outp2[:,1,:,l] = np.swapaxes(ydata[good_idx,:,1],0,1)
                outm[:,0,:,l] = np.swapaxes(xdata[good_idx,:,2],0,1)
                outm[:,1,:,l] = np.swapaxes(ydata[good_idx,:,2],0,1)

            if i0 == j0:
                # This is an auto-correlation
                outa[i0,k,:,l] = data[data.nonzero()]
            else:
                outx[bl2ord[i,j],k,:,l] = data[data.nonzero()]
    out = {'a':outa, 'x':outx, 'uvw':np.array(uvwarray), 'fghz':freq, 'time':np.array(timearray),'source':src,'p':outp,'p2':outp2,'m':outm}
    return out
    
def summary_plot(out,ant_str='ant1-13',ptype='phase'):
    ''' Makes a summary amplitude or phase plot for all baselines from ants in ant_str
        in out dictionary.
    '''
    import matplotlib.pyplot as plt
    
    ant_list = p.ant_str2list(ant_str)
    nant = len(ant_list)
    if ptype != 'amp' and ptype != 'phase':
        print "Invalid plot type.  Must be 'amp' or 'phase'."
        return
    f, ax = plt.subplots(nant,nant)
    f.subplots_adjust(hspace=0,wspace=0)
    bl2ord = p.bl_list()
    for axrow in ax:
        for a in axrow:
            a.xaxis.set_visible(False)
            a.yaxis.set_visible(False)
    for i in range(nant-1):
        ai = ant_list[i]
        for j in range(i+1,nant):
            aj = ant_list[j]
            if ptype == 'phase':
                ax[i,j].imshow(np.angle(out['x'][bl2ord[ai,aj],0]))
                ax[j,i].imshow(np.angle(out['x'][bl2ord[ai,aj],1]))
            elif ptype == 'amp':
                ax[i,j].imshow(np.abs(out['x'][bl2ord[ai,aj],0]))
                ax[j,i].imshow(np.abs(out['x'][bl2ord[ai,aj],1]))
    for i in range(nant):
        ai = ant_list[i]
        ax[i,i].text(0.5,0.5,str(ai+1),ha='center',va='center',transform=ax[i,i].transAxes,fontsize=14)

def get_IDBfiles(showthelast=10):
    #This will return the most recent IDB files saved to the 
    #  directory /data1/IDB/. They will be returned in a list
    #  with the format '/data1/IDB/IDByyyymmddhhmmss'
    #  We can adjust the number of files shown with 'showthelast' optional variable.
    #  It will not show the file being written currently
    s = showthelast
    IDBfiles = os.listdir('/dppdata1/IDB/')
    IDBfiles.sort()
    IDBfiles = IDBfiles[-(s+2):-2]
    for i in range(len(IDBfiles)):
        IDBfiles[i] = '/dppdata1/IDB/'+IDBfiles[i] 
    return IDBfiles
    
def read_idb(trange,filter=True):
    #  This finds the IDB files within a given time range and concatenates the times into a single dictionary
    files = get_trange_files(trange)
    datalist = []
    src = ''
    for file in files:
        #This will skip any files which give us errors.
        #  The names of the bad or unreadable files will
        #  be printed.
        try:
            out = readXdata(file,filter)
            if src == '':
                # This is the first file so set source name and frequency list
                src = out['source']
            if src != out['source']:
                print 'Source name:',out['source'],'is different from initial source name:',src
                print 'Will stop reading files.'
                break
            datalist.append(out)
        except:
            print 'The problematic file is:',file
    # Sometimes an IDB file will be truncated early and have a number of frequencies
    # that is incompatible with the others, so make a list of the number of frequencies 
    # in each IDBdata and keep only those with the most common number.
    nflist = []
    # Find the maximum number of frequencies
    for out in datalist:
        nflist.append(len(out['fghz']))
    nflist = np.array(nflist)
    nfgood = np.median(nflist)
    badidx, = np.where(nflist != nfgood)
    if len(badidx) != 0 and len(nflist) < 3:
        print 'Files have different numbers of frequencies, and result is ambiguous.'
        print 'Will take the file with the higher number of frequencies'
        nfgood = nflist.max()
    # Make a list, nfbad, of entries with fewer than the max number of frequencies
    nfbad, = np.where(nfgood != nflist)
    for i in range(len(nflist)):
        if nflist[i] != nfgood: nfbad.append(i)
    # This loop will run only if length of nfbad > 0
    mybad = np.array(nfbad)
    for i in range(len(nfbad)):
        # Eliminate the entries in the list with too few frequencies
        datalist.pop(mybad[i])
        mybad -= 1
        print 'File',files[nfbad[i]],'eliminated due to too few frequencies.'
    # Have to concatenate outa, outx, uvw, and time arrays
    outa = []
    outx = []
    outp = []
    outp2 = []
    outm = []
    uvw = []
    time = []
    for out in datalist:
        outa.append(out['a'])
        outx.append(out['x'])
        outp.append(out['p'])
        outp2.append(out['p2'])
        outm.append(out['m'])
        uvw.append(out['uvw'])
        time.append(out['time'])
    out['a'] = np.concatenate(outa,3)
    out['x'] = np.concatenate(outx,3)
    out['p'] = np.concatenate(outp,3)
    out['p2'] = np.concatenate(outp2,3)
    out['m'] = np.concatenate(outm,3)
    out['uvw'] = np.concatenate(uvw)
    out['time'] = np.concatenate(time)
    return out

def get_trange_files(trange):
    #Given a timerange, this routine will take all relevant IDBfiles from
    #  that time range, put them in a list, and return that list.
    #  This function is used in get_X_data(data).
    fstr = trange[0].iso
    folder = '/dppdata1/IDB'
    try:
        os.listdir(folder)
    except:
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

    
