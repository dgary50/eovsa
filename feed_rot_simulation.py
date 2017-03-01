import numpy as np
import matplotlib.pylab as plt
from pcapture2 import bl_list
bl2ord = bl_list()

def rot_sim(indict):
    ''' Simulates effect of non-ideal feed behavior on input data for a single baseline.
        If 'unrot' is true, it takes the data as output data and applies the inverse
        of the corresponding Mueller matrix to it.
    
        Input is a dictionary (indict) with the following keys (and their defaults if 
        omitted)
          'data' : 4 x ntimes array of data, corresponding to a set of XX, XY, YX, YY, 
                    in that order, for each time.  No default (error return if omitted)
          'chi1' : Parallactic (or rotation) angle of first antenna (default 0) (scalar or ntimes array)
          'chi2' : Parallactic (or rotation) angle of second antenna (default 0) (scalar or ntimes array)
          'a1'   : Relative amplitude of Y wrt X for first antenna (default 1)
          'a2'   : Relative amplitude of Y wrt X for second antenna (default 1)
          'd1'   : Relative cross-talk between X and Y for first antenna (default 0)
          'd2'   : Relative cross-talk between X and Y for second antenna (default 0)
          'unrot': Whether to rotate or unrotate the input data (default False)
          'verbose': Print some diagnostic messages
          'doplot': Create some plots of the before and after amplitudes and phases
          
        Result is a plot of input and output.  Returns the rotated or unrotated data in 
        the same form as the input data.
    '''
    # Input is a dictionary, contained needed keys.  Any missing
    # keys are filled in with defaults:
    data = indict.get('data')    # No default
    if data is None:
        print 'Must supply "data" key in input dictionary'
        return
    chi1 = indict.get('chi1',0)
    chi2 = indict.get('chi2',0)
    a1 = indict.get('a1',1)
    a2 = indict.get('a2',1)
    d1 = indict.get('d1',0)
    d2 = indict.get('d2',0)
    titl = indict.get('title','')
    verbose = indict.get('verbose',False)
    unrot = indict.get('unrot',False)
    doplot = indict.get('doplot',False)
    titl = titl+' a1:'+str(a1)+' a2:'+str(a2)+' d1:'+str(d1)+' d2:'+str(d2)
    if unrot:
        titl += ' Unrot'
    # Do some sanity checks.  Any or all inputs can have only 1 time (assumed constant at other times)
    # but any that do have times must have the same number of times
    try:
        dn, dnt = data.shape
        if dn != 4:
            print 'Number of data elements for each time must be 4 [XX, XY, YX, YY].'
            return
    except:
        if len(data) == 4:
            dnt = 1
        else:
            print 'Number of data elements for each time must be 4 [XX, XY, YX, YY].'
            return
    try:
        c1nt, = chi1.shape
    except:
        c1nt = 1
    try:
        c2nt, = chi2.shape
    except:
        c2nt = 1
    if dnt > 1 or c1nt > 1 or c2nt > 1:
        # Multiple times are requested, so ensure input is compatible
        nt = np.max(np.array([dnt,c1nt,c2nt]))
        if dnt == 1:
            # Expand data to 4 x nt times
            data = np.rollaxis(np.tile(data,(nt,1)),1)
        if c1nt == 1:
            # Expand chi1 to nt times
            chi1 = chi1*np.ones(nt)
        if c2nt == 1:
            # Expand chi2 to nt times
            chi2 = chi2*np.ones(nt)
        # Now see if they all have the same length
        if nt != len(data[0]) or nt != len(chi1) or nt != len(chi2):
            print 'Number of times in data, chi1 and chi2 are not compatible.'
            return
    else:
        nt = 1
        
    if verbose:
        print 'ntimes=',nt
        print 'shapes of data, chi1, chi2:',data.shape,chi1.shape,chi2.shape
        
    # At this point, we should have uniformity of times
    # Rotation matrix for first antenna
    R1 = np.array([[np.cos(chi1), np.sin(chi1)],[-np.sin(chi1), np.cos(chi1)]])
    if verbose:
        print 'Rotation matrix R1 shape:',R1.shape,'for first time is'
        print R1[:,:,0]
    # Amplitude matrix for first antenna
    A = np.array([[1,d1],[-d1,a1]])
    # Resultant Jones matrix for first antenna
    JA = []
    for i in range(nt):
        JA.append(np.dot(A, R1[:,:,i]))
    if verbose:
        print 'First element of Jones matrix A:\n',JA[0]
    # Rotation matrix for second antenna
    R2 = np.array([[np.cos(chi2), np.sin(chi2)],[-np.sin(chi2), np.cos(chi2)]])
    if verbose:
        print 'Rotation matrix R2 shape:',R2.shape,'for first time is'
        print R2[:,:,0]
    # Amplitude matrix for second antenna
    B = np.array([[1,d2],[-d2,a2]])
    # Resultant Jones matrix for second antenna
    JB = []
    for i in range(nt):
        JB.append(np.dot(B, R2[:,:,i]))
    if verbose:
        print 'First element of Jones matrix B:\n',JB[0]
    # Resultant Mueller matrix
    M = []
    for i in range(nt):
        M.append(np.kron(JA[i],np.conj(JB[i])))
        if verbose and i == 0:
            print 'Mueller matrix at first time is:\n',M[0]
        if unrot:
            M[i] = np.linalg.inv(M[i])
    # Apply matrix to data
    out = np.zeros_like(data)
    if nt == 1:
        out = np.dot(M[0],data)
    else:
        for i in range(nt):
            out[:,i] = np.dot(M[i],data[:,i])
    # Now plot a comparison of the input data to output data:
    if doplot:
        f, ax = plt.subplots(4,2)
        f.set_size_inches(5, 6.5)
        f.suptitle(titl+' Amp')
        pol = ['XX','XY','YX','YY']
        dirxn = [' (IN)',' (OUT)']
        for i in range(4):
            ax[i,0].plot(abs(data[i]),'.')
            ax[i,0].set_ylabel('Rel. Amp')
            ax[i,1].plot(abs(out[i]),'.')
            for j in range(2):
                if i == 3:
                    ax[i,j].set_xlabel('Time index')
                ax[i,j].set_ylim(-0.05,2)
                ax[i,j].grid()
                ax[i,j].text(0.05,0.8,pol[i]+dirxn[j],transform=ax[i,j].transAxes)
        f, ax = plt.subplots(4,2)
        f.set_size_inches(5, 6.5)
        f.suptitle(titl+' Phase')
        pol = ['XX','XY','YX','YY']
        dirxn = [' (IN)',' (OUT)']
        for i in range(4):
            ax[i,0].plot(np.angle(data[i]),'.')
            ax[i,0].set_ylabel('Phase')
            ax[i,1].plot(np.angle(out[i]),'.')
            for j in range(2):
                if i == 3:
                    ax[i,j].set_xlabel('Time index')
                ax[i,j].set_ylim(-4,4)
                ax[i,j].grid()
                ax[i,j].text(0.05,0.8,pol[i]+dirxn[j],transform=ax[i,j].transAxes)
    return out