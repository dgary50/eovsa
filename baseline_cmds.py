from matplotlib import pylab as plt
import glob
import read_idb as ri
from util import lobe
files = glob.glob('/data1/eovsa/fits/UDB/2018/UDB20180826*')
files.sort()
files = files[:20]
out1o = ri.read_idb(files,navg=20)
out2o = ri.read_idb(files[1:],navg=20)
ph1 = np.angle(out1o['x'][ri.bl2ord[13,:13]])
ph2 = np.angle(out2o['x'][ri.bl2ord[13,:13]])
ph_1 = np.zeros((13,2,12,13),np.float)
ph_2 = np.zeros((13,2,12,10),np.float)
for i in range(13):
    for j in range(2):
        for k in range(12):
            if i == 0:
                ph_1[i,0,k] = lobe(ph1[i,j,k]-ph1[i,j,k,0])*out1o['fghz'][0]/out1o['fghz'][k]
                ph_2[i,1,k] = lobe(ph2[i,j,k]-ph2[i,j,k,0])*out2o['fghz'][0]/out2o['fghz'][k]
            else:
                ph_1[i,0,k] = -lobe(ph1[i,j,k]-ph1[0,j,k]-ph1[i,j,k,0]+ph1[0,j,k,0])*out1o['fghz'][0]/out1o['fghz'][k]
                ph_2[i,1,k] = -lobe(ph2[i,j,k]-ph2[0,j,k]-ph2[i,j,k,0]+ph2[0,j,k,0])*out2o['fghz'][0]/out2o['fghz'][k]
# Mean over frequencies and polarizations, scaled by cos(dec)
ph_1t = np.mean(np.mean(ph_1,1),1)/cos(out1o['dec'])
ph_2t = np.mean(np.mean(ph_2,1),1)/cos(out2o['dec'])
f, ax = plt.subplots(4,4)
ax.shape = (16,)
for i in range(13):
    ax[i].cla()
    ax[i].plot(out1o['ha'],ph_1t[i],'.',color='C0')
    ax[i].plot(out2o['ha'],ph_2t[i],'.',color='C1')
    ax[i].set_ylim(-1.4,1)
ha = np.concatenate((out1o['ha'],out2o['ha']))
ph_t = np.concatenate((ph_1t,ph_2t),1)
arg = np.argsort(ha)
has = ha[arg]
ph_ts = ph_t[:,arg]

def basefit(ha, phases, f):
    ''' Given a set of hour angles and phase offsets at each hour angle,
        corresponding to observing frequency f, fit the data for dBx and dBy
        baseline errors and return the fit parameters and x,y arrays of 
        smooth, fitted curves. NB: Phases are expected to already be divided
        by cos(dec). '''
        
    from scipy.optimize import leastsq
    def peval(x,p):
        return p[0]*cos(x) - p[1]*sin(x) + p[2]
    def residuals(p,y,x):
        dBx, dBy, poff = p
        return y - (dBx*cos(x) - dBy*sin(x)) - poff


    cpns = 0.2997924 # speed of light in m/ns
    p0 = [0.0,0.0,0.0]

    # Scale phases errors to delay in ns
    ph_in = phases/(2*np.pi*f/cpns)
    plsq = leastsq(residuals, p0, args=(ph_in, ha))

    ha_out = np.linspace(-1.05,1.05,100)
    y = peval(ha_out, plsq[0])
    # Scale output delay errors back to phase
    ph_out = y*(2*np.pi*f/cpns)
    return plsq[0],ha_out,ph_out
    
cpns = 0.2997924 # speed of light in m/ns
for i in range(13):
    params, ha_out, ph_out = basefit(has,ph_ts[i],out1o['fghz'][0])
    ax[i].plot(ha_out, ph_out)
    ax[i].text(-1,0.7,'dBx, dBy [ns]')
    ax[i].text(-1,0.5,'{:6.4f}, {:6.4f}'.format(params[0],params[1]))
    ax[i].text(-1,-0.7,'dBx, dBy [mm]')
    ax[i].text(-1,-0.9,'{:6.4f}, {:6.4f}'.format(params[0]*1000./cpns,params[1]*1000./cpns))

