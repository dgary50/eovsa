# CASA 6
from casatools import table
from casatools import ms as mstool
ms = mstool()
tb = table()
ptb = table()
spwtb = table()
# CASA 5
#from taskinit import tb
#from taskinit import tb as ptb
#from taskinit import tb as spwtb
#from taskinit import ms
def rd_ms(vis):
    fillnan = None
    spwlist = []
    ms.open(vis)
    tb.open(vis)
    mdata = ms.metadata()
    nspw = mdata.nspw()
    nbl = mdata.nbaselines() + mdata.nantennas()
    nscans = mdata.nscans()
    spw_nfrq = []     # List of number of frequencies in each spw
    for i in range(nspw):
        spw_nfrq.append(mdata.nchan(i))
    spw_nfrq = np.array(spw_nfrq)
    nf = np.sum(spw_nfrq)
    smry = mdata.summary()
    scan_ntimes = []  # List of number of times in each scan
    for iscan in range(nscans):
        scan_ntimes.append(len(smry['observationID=0']['arrayID=0']['scan='+str(iscan)]['fieldID=0'].keys()) - 6)
    scan_ntimes = np.array(scan_ntimes)
    nt = np.sum(scan_ntimes)
    ms.close()
    times = tb.getcol('TIME')
    if times[nbl] - times[0] != 0:
        # This is frequency/scan sort order
        order = 'f'
    elif times[nbl*nspw] - times[0] !=0:
        # This is time sort order
        order = 't'
    ptb.open(vis+'POLARIZATION')
    npol = ptb.getcol('NUM_CORR',0,1)[0]
    ptb.close()
    spwtb.open(vis+'SPECTRAL_WINDOW')
    spec = np.zeros((npol, nf, nbl, nt), np.complex)
    freq = np.zeros(nf, float)
    time = np.zeros(nt, float)
    if order == 't':
        for j in range(nt):
            fptr = 0
            # Loop over spw
            for i in range(nspw):
                cfrq = spwtb.getcol('CHAN_FREQ',i,1)[:,0]     # Get channel frequencies for this spw (annoyingly comes out as shape (nf, 1)
                if j == 0:
                    # Only need this the first time through
                    spwlist += [i]*len(cfrq)
                if i == 0:
                    time[j] = tb.getcol('TIME',nbl*(i+nspw*j),1)  # Get the time
                spec_ = tb.getcol('DATA',nbl*(i+nspw*j),nbl)  # Get complex data for this spw
                flag = tb.getcol('FLAG',nbl*(i+nspw*j),nbl)   # Get flags for this spw
                nfrq = len(cfrq)
                # Apply flags
                if type(fillnan) in [int, float]:
                    spec_[flag] = float(fillnan)
                else:
                    spec_[flag] = 0.0
                # Insert data for this spw into larger array
                spec[:, fptr:fptr+nfrq, :, j] = spec_
                freq[fptr:fptr+nfrq] = cfrq
                fptr += nfrq
            self.progressBar.setValue(100*j/nt)
    else:
        iptr = 0
        for j in range(nscans):
            # Loop over scans
            for i in range(nspw):
                #Loop over spectral windows
                s = scan_ntimes[j]
                f = spw_nfrq[i]
                s1 = np.sum(scan_ntimes[:j])     # Start time index
                s2 = np.sum(scan_ntimes[:j+1])   # End time index
                f1 = np.sum(spw_nfrq[:i])         # Start freq index
                f2 = np.sum(spw_nfrq[:i+1])       # End freq index
                spec_ = tb.getcol('DATA',iptr,nbl*s)
                flag  = tb.getcol('FLAG',iptr,nbl*s)
                if j == 0:
                    cfrq = spwtb.getcol('CHAN_FREQ',i,1)[:,0]
                    freq[f1:f2] = cfrq
                    spwlist += [i]*len(cfrq)
                time[s1:s2] = tb.getcol('TIME', iptr, nbl*s).reshape(s, nbl)[:,0]  # Get the times
                iptr += nbl*s
                # Apply flags
                if type(fillnan) in [int, float]:
                    spec_[flag] = float(fillnan)
                else:
                    spec_[flag] = 0.0
                # Insert data for this spw into larger array
                spec[:, f1:f2, :, s1:s2] = spec_.reshape(npol,f,nbl,s)
    tb.close()
    spwtb.close()
    return {'data':spec, 'freq':freq, 'time':time, 'name':ms.name, 'spwlist':np.array(spwlist)}

                