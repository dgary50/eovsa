import read_idb as ri
from util import Time
import numpy as np
from adc_cal2 import ant_str2list

def make_bl_table(trange,antstr='ant1-13'):
    files = ri.get_trange_files(trange)
    bl2ord = ri.p.bl_list()
    antidx = ant_str2list(antstr)
    fxx = open(trange[0].iso[:10].replace('-','')+'_baseline_xx.txt','w')
    fyy = open(trange[0].iso[:10].replace('-','')+'_baseline_yy.txt','w')
    fxy = open(trange[0].iso[:10].replace('-','')+'_baseline_xy.txt','w')
    fyx = open(trange[0].iso[:10].replace('-','')+'_baseline_yx.txt','w')
    for i,file in enumerate(files):
        print 'Reading',file,'('+str(i+1)+' of '+str(len(files))+')'
        out = ri.read_idb([file],navg=60)
        data = out['x'][bl2ord[antidx,13],:,:,1:]
        phase = np.angle(np.mean(data,3))*180./np.pi
        amp = np.abs(np.mean(data,3))
        x = np.mean(np.real(data),3)
        y = np.mean(np.imag(data),3)
        dx = np.std(np.real(data),3)
        dy = np.std(np.imag(data),3)
        dphase = (np.sqrt(x**2 * dy**2 + y**2 * dx**2)/amp**2)*180./np.pi
        xx = phase[:,0,0]
        yy = phase[:,1,0]
        xy = phase[:,2,0]
        yx = phase[:,3,0]
        dxx = dphase[:,0,0]
        dyy = dphase[:,1,0]
        dxy = dphase[:,2,0]
        dyx = dphase[:,3,0]
        outdict = src2dict(out)
        h = np.mean(outdict['ha'])
        fxx.write(outdict['source'].name+'  {:6.2f} {:6.2f}  '.format(h*180./np.pi,outdict['dec']*180/np.pi)+('{:6.1f} '*10).format(*xx)+'\n')
        fyy.write(outdict['source'].name+'  {:6.2f} {:6.2f}  '.format(h*180./np.pi,outdict['dec']*180/np.pi)+('{:6.1f} '*10).format(*yy)+'\n')
        fxy.write(outdict['source'].name+'  {:6.2f} {:6.2f}  '.format(h*180./np.pi,outdict['dec']*180/np.pi)+('{:6.1f} '*10).format(*xy)+'\n')
        fyx.write(outdict['source'].name+'  {:6.2f} {:6.2f}  '.format(h*180./np.pi,outdict['dec']*180/np.pi)+('{:6.1f} '*10).format(*yx)+'\n')
    fxx.close()
    fyy.close()
    fxy.close()
    fyx.close()
    
def src2dict(out):
    import aipy
    import eovsa_array
    import eovsa_cat
    import eovsa_lst
    import copy

    bl2ord = ri.p.bl_list()
    uv = copy.copy(out['uvw'][:,:,:2])
    time = Time(out['time'],format='jd')
    # Get some of the other information we need for this source
    srclist = eovsa_cat.load_VLAcals()
    cat = aipy.amp.SrcCatalog(srclist)
    aa = eovsa_array.eovsa_array()
    cat.compute(aa)
    aa.cat = cat
    src = aa.cat[out['source']]
    src.name = out['source']
    ha = np.zeros(len(time),'double')
    for i in range(len(time)):
        ha[i] = eovsa_lst.eovsa_ha(src,time[i])
    phase = np.angle(out['x'][bl2ord[:,13]])
    outdict = {'phase':phase,'time':time,'ha':ha,'dec':src.dec,'fghz':out['fghz'],'source':src,'uv':uv}
    return outdict
    
def bz(fname1,fname2):
    ''' This routine assumes that two sources were observed, and the results
        read in and pickled, using out['x'] = out['x'][bl2ord[:,13]] to
        reduce the data to only baselines with ant14.  The filenames of 
        the two pickled files is provided.
    '''
    import cPickle as pickle
    f = open(fname1,'rb')
    out1 = pickle.load(f)
    f.close()
    f = open(fname2,'rb')
    out2 = pickle.load(f)
    nant, npol, nf, nt = out1['x'].shape
    f.close()
    x1, = np.where(out1['ha'] - np.roll(out1['ha'],1) > 0.005)
    x2, = np.where(out2['ha'] - np.roll(out2['ha'],1) > 0.005)
    x1 = [0]+x1.tolist()+[len(out1['ha'])]
    x2 = [0]+x2.tolist()+[len(out2['ha'])]
    p1 = np.zeros([13,npol,nf,len(x1)-1],float)
    h1 = np.zeros(len(x1)-1,float)
    for i in range(len(x1)-1):
        p1[:,:,:,i] = np.angle(np.sum(out1['x'][:13,:,:,x1[i]+1:x1[i+1]-1],3))
        h1[i] = np.mean(out1['ha'][x1[i]+1:x1[i+1]-1])
    p2 = np.zeros([13,npol,nf,len(x2)-1],float)
    h2 = np.zeros(len(x2)-1,float)
    for i in range(len(x2)-1):
        p2[:,:,:,i] = np.angle(np.sum(out2['x'][:13,:,:,x2[i]+1:x2[i+1]-1],3))
        h2[i] = np.mean(out2['ha'][x2[i]+1:x2[i+1]-1])
