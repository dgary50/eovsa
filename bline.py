import read_idb as ri
from util import Time
import eovsa_array as ea
import eovsa_lst as el
import eovsa_cat as ec
import numpy as np
import aipy
import copy

def amphz(out):
    bl2ord = ri.p.bl_list()
    uv = copy.copy(out['uvw'][:,:,:2])
    time = Time(out['time'],format='jd')
    # Get some of the other information we need for this source
    srclist = ec.load_VLAcals()
    cat = aipy.amp.SrcCatalog(srclist)
    aa = ea.eovsa_array()
    cat.compute(aa)
    aa.cat = cat
    src = aa.cat[out['source']]
    src.name = out['source']
    ha = np.zeros(len(time),'double')
    for i in range(1,len(time)):
        ha[i] = el.eovsa_ha(src,time[i])
    # Select only baselines with Ant 14, and eliminate XY and YX polarizations, and the first time
    data = out['x'][bl2ord[:13,13],:2,:,1:]
    phz = np.angle(data)
    a0 = abs(np.sum(data,3))    # Time averaged amplitude
    damp = np.std(abs(data),3) 
    p0 = np.angle(np.sum(data,3))  # Time averaged angle
    nbl,npol,nf,nt = data.shape
    pdif = np.zeros(data.shape,'float')
    for i in range(nt):
        pdif[:,:,:,i] = phz[:,:,:,i] - p0
    pdif[np.where(pdif > np.pi)] -= 2*np.pi
    pdif[np.where(pdif < -np.pi)] += 2*np.pi
    dphase = np.std(pdif,3)    
    outdict = {'amp':a0,'damp':damp,'phase':p0,'dphase':dphase,'time':time,'ha':np.mean(ha),'dec':src.dec,  
               'fghz':out['fghz'],'source':src,'uv':uv}
    return outdict

def bline_fit(ps):
    ''' Given a list of outdicts, ps, generate phase differences between adjacent measurements
        and fit to baseline errors.  All outdicts in the list have to have the same shape data.
        
        Expression for differences in phase of adjacent measurements:
            dphi = 2*pi*fghz*{dbx*[cos(dec_i)*cos(ha_i) - cos(dec_i-1)*cos(ha_i-1)] 
                            - dby*[cos(dec_i)*sin(ha_i) - cos(dec_i-1)*sin(ha_i-1)]) 
                            + dbz*(sin(dec_i)-sin(dec_i-1))}
        which is of the form:
            A = M deltaB
        with solution
            deltaB = inv(M) A
        where A_i are the measured phase differences, B is a vector [dbx, dby, dbz], and the matrix
            M = [[cx_1   cy_1   cz_1]], ..., [cx_N    cy_N    cz_N]] is a 3 x N matrix
        where
            cx_i =  2*pi*fghz*[cos(dec_i)*cos(ha_i) - cos(dec_i-1)*cos(ha_i-1)] 
            cy_i = -2*pi*fghz*[cos(dec_i)*sin(ha_i) - cos(dec_i-1)*sin(ha_i-1)] 
            cz_i =  2*pi*fghz*[sin(dec_i) - sin(dec_i-1)]
            
        This expression needs to be evaluated for a given dbx, dby, dbz, for all frequencies and
        both polarizations, for each of the 13 baselines.
        
        The output is the vector of baseline errors deltaB
    '''
    nm = len(ps)-1  # one less than number of measurements
    nbl, npol, nf = ps[0]['phase'].shape
    ninst = nm*npol*nf
    fghz = ps[0]['fghz']
    deltaB = np.zeros((nbl,3),'float')
    cs = np.zeros((nm,3))
    # Loop over baselines (start over for each baseline)
    for j in range(nbl):
        inst = -1  # start over with counter--each baseline is independent
        M = np.zeros((ninst,3),'float')
        A = np.zeros(ninst,'float')
        dphi = np.zeros((npol,nf,nm),'float')
        for i in range(nm):
            d1 = ps[i+1]['dec']
            d0 = ps[i]['dec']
            h1 = ps[i+1]['ha']
            h0 = ps[i]['ha']
            dphi = ps[i+1]['phase'][j] - ps[i]['phase'][j]
            cx = 2*np.pi*(np.cos(d1)*np.cos(h1)-np.cos(d0)*np.cos(h0))
            cy = -2*np.pi*(np.cos(d1)*np.sin(h1)-np.cos(d0)*np.sin(h0))
            cz = 2*np.pi*(np.sin(d1)-np.sin(d0))
            cs[i] = [cx,cy,cz]   # Save these calculations for later
            # Loop over frequencies
            for k in range(nf):
                inst += 1
                M[inst,0] = fghz[k]*cx
                M[inst,1] = fghz[k]*cy
                M[inst,2] = fghz[k]*cz
                A[inst] = dphi[0,k]
                inst += 1
                M[inst,0] = fghz[k]*cx
                M[inst,1] = fghz[k]*cy
                M[inst,2] = fghz[k]*cz
                A[inst] = dphi[1,k]
        deltaB[j], var, _, _ = np.linalg.lstsq(M,A)
    # Convert the results to phase differences and list the O-C differences
    fmt = ' {:5.2f}'*13
    for i in range(nm):
        print 'O-C differences for '+ps[i+1]['source']+'-'+ps[i]['source']
        print 'Baselines 1-14  2-14  3-14  4-14  5-14  6-14  7-14  8-14  9-14 10-14 11-14 12-14 13-14'
        for k in range(nf):
            dphi_O = ps[i+1]['phase'][:13,0,k] - ps[i]['phase'][:13,0,k]
            dphi_C = deltaB[:13,0]*cs[i,0]+deltaB[:13,1]*cs[i,1]+deltaB[:13,2]*cs[i,2]
            OmC = (dphi_O - dphi_C)
            OmC[np.where(OmC > np.pi)] -= 2*np.pi
            OmC[np.where(OmC < -np.pi)] += 2*np.pi
            print str(fghz[k])[:5]+'GHz'+fmt.format(*OmC)
    return deltaB, var
    
def brute_fit(ps):
    ''' Given a list of outdicts, ps, generate phase differences between adjacent measurements
        and fit to baseline errors.  All outdicts in the list have to have the same shape data.
        
        Expression for differences in phase of adjacent measurements:
            dphi = 2*pi*fghz*{dbx*[cos(dec_i)*cos(ha_i) - cos(dec_i-1)*cos(ha_i-1)] 
                            - dby*[cos(dec_i)*sin(ha_i) - cos(dec_i-1)*sin(ha_i-1)]) 
                            + dbz*(sin(dec_i)-sin(dec_i-1))}
        This version does a brute-force trial-and-error solution.
        
        
        
        The output is the vector of baseline errors deltaB
    '''
    nm = len(ps)-1  # one less than number of measurements
    nbl, npol, nf = ps[0]['phase'].shape
    ninst = nm*npol*nf
    fghz = ps[0]['fghz']
    cx = np.zeros((nm,nf))
    cy = np.zeros((nm,nf))
    cz = np.zeros((nm,nf))
    # Loop over baselines (start over for each baseline)
    for j in range(10,11):
#        inst = -1  # start over with counter--each baseline is independent
#        M = np.zeros((ninst,3),'float')
#        A = np.zeros(ninst,'float')
        dphi_O = np.zeros((nm,npol,nf),'float')
        for i in range(nm):
            d1 = ps[i+1]['dec']
            d0 = ps[i]['dec']
            h1 = ps[i+1]['ha']
            h0 = ps[i]['ha']
            dphi_O[i] = ps[i+1]['phase'][j] - ps[i]['phase'][j]  # nm, npol, nf
            cx[i] = 2*np.pi*(np.cos(d1)*np.cos(h1)-np.cos(d0)*np.cos(h0))*fghz
            cy[i] = -2*np.pi*(np.cos(d1)*np.sin(h1)-np.cos(d0)*np.sin(h0))*fghz
            cz[i] = 2*np.pi*(np.sin(d1)-np.sin(d0))*fghz
        # At this point, we have the measurements dphi, and coefficients we need to calculate corrections
        # Now embark on a search in a 3-d baseline error space to find the optimum baselines
        mn = -5
        mx = 5
        step = 0.1
        bl_space = np.zeros(((mx-mn)/step,(mx-mn)/step,(mx-mn)/step),'float')
        for bx in np.arange(mn,mx,step):
            for by in np.arange(mn,mx,step):
                for bz in np.arange(mn,mx,step):
                    # These would all be zero for correct bx, by, bz
                    phz = (cx*bx + cy*by + cz*bz - dphi_O[:,0,:]) % 2*np.pi  # nm, nf (only 1 poln!)
                    # Ensure that they all lie between -pi and pi
                    phz[np.where(phz > np.pi)] -= 2*np.pi
                    phz[np.where(phz < -np.pi)] += 2*np.pi
                    # A simple average (should be weighted by errors somehow)
                    pavg = np.mean(abs(phz[:,0]))  # only 1 freq
                    bl_space[(bx-mn)/step,(by-mn)/step,(bz-mn)/step] = pavg
        return bl_space
        
import get_bldata as gb
def bz(out1,out2):
    phi1 = angle(sum(out1['x'][bl2ord[:13,13]],3))
    phi2 = angle(sum(out2['x'][bl2ord[:13,13]],3))
    outdict1 = gb.src2dict(out1)
    outdict2 = gb.src2dict(out2)
    fghz = out1['fghz'][array([0,1,3,4,5,6])]
    phi1 = phi1[:,:,array([0,1,3,4,5,6])]
    phi2 = phi2[:,:,array([0,1,3,4,5,6])]
    c = 2*pi*(sin(outdict1['dec']) - sin(outdict2['dec']))
    nbl,npol,nf = phi1.shape
    f, ax = subplots(npol,nbl)
    res = zeros(npol,'float')
    for i in range(nbl):
        for j in range(npol):
            dphi = unwrap(phi1[i,j] - phi2[i,j])
            ax[j,i].plot(fghz,dphi,'o')
            ax[j,i].set_ylim(ax[j,i].yaxis.get_data_interval()+array((-1,1)))
            res[j], junk =  polyfit(fghz,dphi,1)
            print 'Pol',j+1,':',res[j],'=> Bz =',res[j]/c,
            print ' '
        print 'Baseline',i+1,'Avg Bz:',mean(res)/c
    show()