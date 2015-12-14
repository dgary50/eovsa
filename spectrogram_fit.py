'''
   Module for plotting EOVSA data as a spectrogram'''
#
# History:
#   2015-May-09  DG
#     First written.
#   2015-Jun-04  DG
#     Major enhancements to add Spectrogram() class and many methods for working
#     with it.
#   2015-Jun-07  DG
#     Additional enhancements, especially to do total-power fits with the Staehli
#     function.  Also added a version number.
#   2015-Jun-08  RG
#     Additional enhancement, now finds a fit for the four parameters of the Staehli 
#     function. 
#   2015-Jun-22  DG
#     Changes to reflect new selection of science bands, which resulted in a lot
#     of frequencies with all zero data.  Just eliminate those from the arrays.
#     Also updated the spectrum-fitting routines to return physical parameters
#        f_pk, S_pk, low-f slope, high-f slope
#   2015-Jun-27  DG
#     Changed the code to standardize on names, content, and order of indices
#     of outputs for rd_miriad_tsys.  Names ut_mjd, fghz, and tsys will be used, 
#     with units hopefully made obvious.  Order of indices will be 
#          (nant/nbl, npol, nfreq, ntimes).
#     This required some changes to several of the routines here.
#   2015-Jun-30  DG
#     Added suptitle() to explore plot.
#   2015-Jul-03  DG
#     Changes to plot_spectrogram() to handle case of plotting xdata (cross-correlation 
#     amplitude)
#
__version__ = '0.2'

def log_sample(fghz, ut, tsys):
    ''' Resamples a spectrogram from an irregular sampling in linear frequency space
        to a regular grid in log frequency space.
    '''
    from scipy import interpolate
    import numpy as np
    nf = len(fghz)
    fghzl = np.logspace(np.log10(fghz[0]),np.log10(fghz[-1]),nf)
    x, y = np.meshgrid(ut,fghzl)
    fint = interpolate.interp2d(ut,fghz,tsys,kind='cubic')
    out = fint(ut,fghzl)
    return fghzl,out
    
def lin_sample(fghz, ut, tsys):
    ''' Resamples a spectrogram from an irregular sampling in linear frequency space
        to a regular grid in linear frequency space.
    '''
    from scipy import interpolate
    import numpy as np
    nf = len(fghz)
    fghzl = np.linspace(fghz[0],fghz[-1],nf)
    x, y = np.meshgrid(ut,fghzl)
    fint = interpolate.interp2d(ut,fghz,tsys,kind='cubic')
    out = fint(ut,fghzl)
    return fghzl,out

def plot_spectrogram(fghz, ut, tsys, ax=None, cbar=True, logsample=False, **kwargs):
    ''' Creates standard spectrogram plot for EOVSA data, using axes supplied
        or creates a new single axis plot if None.  The intensities are log-scaled,
        the xaxis is interpreted as time (ut can be timestamps or plot_date format)
        
        kwargs:
            dmin    Clip data to this minimum value [sfu] (default = 10 sfu)
            dmax    Clip data to this maximum value [sfu] (default = tsys.max())
            xlabel  String to use as xlabel 
                      (default is 'Time [UT on YYYY-MM-DD]' if ax is supplied--none otherwise)
            ylabel  String to use as ylabel 
                      (default is 'Frequency [GHz]' if ax is supplied--none otherwise)
            title   String to use as plot title
                      (default is no title)
            logsample  Boolean. If True, resample the data on a logarithmic 
                       frequency space, if False resample on a linear space,
                       if None, do no resampling in frequency
            xdata   Boolean.  If True, label graph for cross-correlation.  If False, 
                       or omitted, lable graph for Total Power
    '''
    import matplotlib.pylab as plt
    import matplotlib.dates
    import numpy as np
    utd = ut.plot_date
    datstr = ut[0].iso[:10]
        
    if ax is None:
        # No axes supplied, so create one (and assume labels are wanted)
        f, ax = plt.subplots(1,1)
        ax.set_xlabel('Time [UT on '+datstr+']')
        ax.set_ylabel('Frequency [GHz]')
        ax.set_title('EOVSA Total Power for '+datstr)
        if 'xdata' in kwargs.keys():
            if kwargs['xdata'] is True:
                ax.set_title('EOVSA Summed Cross-Correlation Amplitude for '+datstr)

    ax.xaxis.set_tick_params(width=1.5,size=10,which='both')
    ax.yaxis.set_tick_params(width=1.5,size=10,which='both')
    if logsample:
        # Sample data only a uniform logarithmic frequency space
        fghzl, tsysl = log_sample(fghz, utd, tsys)
        ax.set_yscale('log')
        minorFormatter = plt.LogFormatter(base=10, labelOnlyBase=False)
        ax.yaxis.set_minor_formatter(minorFormatter)
    elif logsample is None:
        # No sampling in frequency space
        fghzl, tsysl = fghz, tsys
    else:
        # logsample is not True or None (presumably False), so resample 
        # data on a uniform linear frequency space
        fghzl, tsysl = lin_sample(fghz, utd, tsys)

    dmin = 1.
    if 'dmin' in kwargs.keys():
        if kwargs['dmin'] is not None:
            dmin = kwargs['dmin']
    dmax = tsys.max()
    if 'dmax' in kwargs.keys():
        if kwargs['dmax'] is not None:
            dmax = kwargs['dmax']

    # Take logarithm if TP, but not for Cross-correlation amplitude
    data = np.log10(np.clip(tsysl,dmin,dmax))
    if 'xdata' in kwargs.keys():
        if kwargs['xdata'] is True:
            data = np.clip(tsysl,dmin,dmax)

    im = ax.imshow(data,origin='lower',extent=[utd[0],utd[-1],fghzl[0],fghzl[-1]],
                     aspect='auto',interpolation='nearest')

    if cbar: 
        cbar_label = 'Log Flux Density [sfu]'
        if 'xdata' in kwargs.keys():
            if kwargs['xdata'] is True:
                cbar_label = 'Amplitude [arb. units]'
        plt.colorbar(im,ax=ax,label=cbar_label)
    ax.xaxis_date()
    ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter("%H:%M:%S"))
    
    # Set labels if requested
    if 'xlabel' in kwargs.keys():
        if kwargs['xlabel'] == 'auto':
            ax.set_xlabel('Time [UT on '+datstr+']')
        else:
            ax.set_xlabel(kwargs['xlabel'])
    if 'ylabel' in kwargs.keys():
        if kwargs['ylabel'] == 'auto':
            ax.set_ylabel('Frequency [GHz]')
        else:
            ax.set_ylabel(kwargs['ylabel'])
    if 'title' in kwargs.keys():
        ax.set_title(kwargs['title'])
    return ax

import numpy as np
import dump_tsys
from util import Time, common_val_idx
import offline

class Spectrogram():
    
    def __init__(self,trange):
        ''' Create the object for the specified timerange specified by the 
            2-element Time() trange.  The timerange is used to create a list 
            of Miriad database files to read, and the data are read.
        '''
        # Read data
        out = dump_tsys.rd_miriad_tsys(trange)
        nant, npol, nf, nt = out['tsys'].shape
        self.xdata = out['tsys'][:,0,:,:]
        self.ydata = out['tsys'][:,1,:,:]
        self.fghz = out['fghz']
        self.time = Time(out['ut_mjd'],format='mjd')
        self.tidx = [0,len(self.time)]
        # Read calibration
        fghz, self.calfac, self.offsun = offline.read_calfac(trange[0])
        # Make sure frequencies in data and calibration agree
        calidx, dataidx = common_val_idx((fghz*1000).astype('int'),(self.fghz*1000).astype('int'))
        self.bidx = [0,100]
        self.fghz = self.fghz[dataidx]
        self.xdata = self.xdata[:,dataidx,:]
        self.ydata = self.ydata[:,dataidx,:]
        if self.calfac is not None:
            # Select frequencies and swap axes to put into standard form
            self.calfac = np.rollaxis(self.calfac[:,calidx,:],2)
            self.offsun = np.rollaxis(self.offsun[:,calidx,:],2)
            fghz = fghz[calidx]
        # Set frequency index range (fidx) to default to show only frequencies > 2.5 GHz
        lowf, = np.where(self.fghz > 2.5)
        self.fidx = [lowf[0],len(self.fghz)]
        self.drange = [None,None]
        self.antlist = range(nant)
        self.cbar = True
        self.showants = range(nant)
        self.domedian = True
        self.docal = True
        self.dolog = False
        self.dosub = True
        self.ax = None
        self.version = __version__

    def show(self):
        ''' Create a spectrogram plot of the data
        '''
        tsys, stdtsys = self.get_data()
        if self.domedian:
            # This results in only a single plot
            self.ax = plot_spectrogram(self.fghz[self.fidx[0]:self.fidx[1]], self.time[self.tidx[0]:self.tidx[1]], tsys, ax=self.ax, cbar=self.cbar, logsample=self.dolog, dmin=self.drange[0], dmax=self.drange[1])
        else:
            print 'Cannot (yet) plot data for each anteanna separately.  Please set <self>.domedian = True first'
        
        
    def get_median_data(self, xtsys=None, ytsys=None):
        ''' Get optionally calibrated, optionally background subtracted
            data as median over polarization and antenna list in self.showants
        '''
        if xtsys is None:
            if self.docal:
                # Do calibration and optionally subtraction
                xtsys, ytsys = self.get_cal_data()
                if self.dosub:
                    xtsys, ytsys = self.get_bgsub_data(xtsys, ytsys)
            elif self.dosub:
                # No calibration, so select data and do subtraction
                xtsys = self.xdata[self.antlist,self.fidx[0]:self.fidx[1],self.tidx[0]:self.tidx[1]]
                ytsys = self.ydata[self.antlist,self.fidx[0]:self.fidx[1],self.tidx[0]:self.tidx[1]]
                xtsys, ytsys = self.get_bgsub_data(xtsys, ytsys)
        medxtsys = np.nanmedian(xtsys[self.showants,:,:],0)
        stdxtsys = np.nanstd(xtsys[self.showants,:,:],0)
        medytsys = np.nanmedian(ytsys[self.showants,:,:],0)
        stdytsys = np.nanstd(ytsys[self.showants,:,:],0)
        tsys = (medxtsys+medytsys)/2.
        stdtsys = np.sqrt(stdxtsys**2 + stdytsys**2)/2.
        return tsys, stdtsys
        
    def get_cal_data(self):
        # Select data
        xtsys = np.zeros((len(self.antlist), self.fidx[1] - self.fidx[0], self.tidx[1] - self.tidx[0]),dtype='float')
        ytsys = np.zeros((len(self.antlist), self.fidx[1] - self.fidx[0], self.tidx[1] - self.tidx[0]),dtype='float')
        # Apply calibration
        for i,j in enumerate(range(self.tidx[0], self.tidx[1])):
            xtsys[:,:,i] = ((self.xdata[self.antlist, self.fidx[0]:self.fidx[1], j] - 
                           self.offsun[self.antlist, 0, self.fidx[0]:self.fidx[1]])
                           *self.calfac[self.antlist, 0, self.fidx[0]:self.fidx[1]])
            ytsys[:,:,i] = ((self.ydata[self.antlist, self.fidx[0]:self.fidx[1], j] - 
                           self.offsun[self.antlist, 1, self.fidx[0]:self.fidx[1]])
                           *self.calfac[self.antlist, 1, self.fidx[0]:self.fidx[1]])
        return xtsys, ytsys

    def get_bgsub_data(self,xtsys=None, ytsys=None):
        ''' Get optionally calibrated data after background subtraction is applied.
        '''
        if xtsys is None:
            if self.docal:
                # Do calibration and optionally subtraction
                xtsys, ytsys = self.get_cal_data()
                if self.dosub:
                    xtsys, ytsys = self.get_bgsub_data(xtsys, ytsys)
            else:
                # No calibration, so select data and raw subtraction
                xtsys = self.xdata[self.antlist,self.fidx[0]:self.fidx[1],self.tidx[0]:self.tidx[1]]
                ytsys = self.ydata[self.antlist,self.fidx[0]:self.fidx[1],self.tidx[0]:self.tidx[1]]
            
        # Perform the background subtraction
        bgx, bgy = self.getbg(self.bidx, xtsys, ytsys)
        nt = self.tidx[1] - self.tidx[0]
        nf = self.fidx[1] - self.fidx[0]
        for i in range(nt):
            xtsys[:,0:nf,i] -= bgx
            ytsys[:,0:nf,i] -= bgy
        return xtsys, ytsys
        

    def get_data(self):
        ''' Get optionally calibrated, optionally background-subtracted data.
            If self.domedian is True, return the median of data over polarization
            and antenna list in self.showants.
        '''
        if self.docal:
            xtsys, ytsys = self.get_cal_data()
        else:
            xtsys = self.xdata[self.antlist,self.fidx[0]:self.fidx[1],self.tidx[0]:self.tidx[1]]
            ytsys = self.ydata[self.antlist,self.fidx[0]:self.fidx[1],self.tidx[0]:self.tidx[1]]

        if self.dosub:
            xtsys, ytsys = self.get_bgsub_data(xtsys, ytsys)
                
        if self.domedian:
            tsys, stdtsys = self.get_median_data(xtsys, ytsys)
        else:
            tsys = np.swapaxes(np.array([xtsys, ytsys]),1,0)
            stdtsys = None
            
        return tsys, stdtsys
        
    def getbg(self, bidx=None,  xtsys=None, ytsys=None):
        ''' Get background spectra for each antenna and polarization, applying
            calibration first if indicated by self.docal = True
        '''
        if bidx is None:
            bidx = self.bidx
        else:
            self.bidx = bidx
            
        if xtsys is None:
            # No data supplied, so generate it
            if self.docal:
                # Do calibration
                xtsys, ytsys = self.get_cal_data()
            else:
                # No calibration desired, so just select raw data
                xtsys = self.xdata[self.antlist,self.fidx[0]:self.fidx[1],self.tidx[0]:self.tidx[1]]
                ytsys = self.ydata[self.antlist,self.fidx[0]:self.fidx[1],self.tidx[0]:self.tidx[1]]

        # Generate median over supplied background indexes.  These are spectra for each antenna in self.antlist
        bgx = np.nanmedian(xtsys[:,:,bidx[0]:bidx[1]],2)
        bgy = np.nanmedian(ytsys[:,:,bidx[0]:bidx[1]],2)
        return bgx, bgy

    def explore(self):
        ''' Like show(), but provides a mouse-driven interface for exploring the plot
            after creation.  Only works for median data.
        '''
        import matplotlib.pylab as plt
        tsys, stdtsys = self.get_median_data()
        dlogtsys = stdtsys/tsys
        fig = plt.figure(figsize=(8,6))
        fig.suptitle('EOVSA Spectrogram '+self.time[self.tidx[0]].iso[:10],fontsize=25)
        spectrogram_ax = plt.axes([0.1,0.4,0.6,0.5])
        spectrogram_ax.set_ylabel('Frequency [GHz]')
        spectrum_ax = plt.axes([0.75,0.4,0.23,0.5])
        spectrum_ax.set_ylim(1.,tsys.max())
        try:
            spectrum_ax.set_yscale('log')
        except:
            pass
        spectrum_ax.set_xscale('log')
        spectrum_ax.set_xlim(1,18)
        spectrum_ax.set_ylabel('Flux Density [sfu]')
        spectrum_ax.set_xlabel('Frequency [GHz]')
        # Set initial spectruma and lightcurve to correspond to mid-range
        # time and frequency
        fghz = self.fghz[self.fidx[0]:self.fidx[1]]
        t = self.time[self.tidx[0]:self.tidx[1]].plot_date
        midf = len(fghz)/2
        midt = len(t)/2
        spec, = spectrum_ax.plot(fghz,tsys[:,midt],'.')
        p, ffit, sfit = tpfit(np.log(fghz),np.log(tsys[:,midt]),sigma=dlogtsys[:,midt])
        specfit, = spectrum_ax.plot(np.exp(ffit),np.exp(sfit))
        specpt, = spectrum_ax.plot(fghz[midf],tsys[midf,midt],'<',markersize=5,c='y')
        if self.cbar:
            lc_ax = plt.axes([0.1,0.1,0.48,0.25],sharex=spectrogram_ax)
        else:
            lc_ax = plt.axes([0.1,0.1,0.6,0.25],sharex=spectrogram_ax)
        lc, = lc_ax.plot_date(t,tsys[midf,:],'-')
        lcpt, = lc_ax.plot_date(t[midt],tsys[midf,midt],'^',markersize=5,c='y')
        lc_ax.set_ylim(1.,tsys.max())
        try:
            lc_ax.set_yscale('log')
        except:
            pass
        lc_ax.set_ylabel('Flux Density [sfu]')
        lc_ax.set_xlabel('Time [UT]')
        tstr = Time(t[midt],format='plot_date').iso[11:19]
        lctxt = lc_ax.text(0.02,0.9,'{:} UT, {:0.3f} GHz, {:0.3f} sfu'.format(tstr,fghz[midf],tsys[midf,midt]), transform=lc_ax.transAxes)
        sptxt = spectrum_ax.text(0.02,0.9,'{:} UT, {:0.3f} GHz, {:0.3f} sfu'.format(tstr,fghz[midf],tsys[midf,midt]), transform=spectrum_ax.transAxes)
        self.ax = spectrogram_ax
        self.show()

        def find_ij(x, y):
            return abs(utd - x).argmin(), abs(fghz - y).argmin()

        def onmove(event):
            if event.inaxes != spectrogram_ax: return
            i, j = abs(t - event.xdata).argmin(), abs(fghz - event.ydata).argmin()
            #print 'indexes are t=%f, f=%f'%(i, j)
            spec.set_data(fghz,tsys[:,i])
            p, ffit, sfit = tpfit(np.log(fghz),np.log(tsys[:,i]),sigma=dlogtsys[:,i])
            specfit.set_data(np.exp(ffit),np.exp(sfit))
            specpt.set_data(fghz[j],tsys[j,i])
            lc.set_data(t,tsys[j,:])
            lcpt.set_data(t[i],tsys[j,i])
            tstr = Time(t[i],format='plot_date').iso[11:19]
            lctxt.set_text('{:} UT, {:0.3f} GHz, {:0.3f} sfu'.format(tstr,fghz[j],tsys[j,i]))
            sptxt.set_text('{:} UT, {:0.3f} GHz, {:0.3f} sfu'.format(tstr,fghz[j],tsys[j,i]))
            fig.canvas.draw()
    
        def onclick(event):
            if event.inaxes != spectrogram_ax: return
            # print 'Turning on mouse-move events'
            fig.mid = fig.canvas.mpl_connect('motion_notify_event', onmove)
    
        def onrelease(event):
            if event.inaxes != spectrogram_ax: return
            # print 'Turning off mouse-move events'
            fig.canvas.mpl_disconnect(fig.mid)
    
        self.cid = fig.canvas.mpl_connect('button_press_event', onclick)
        self.rid = fig.canvas.mpl_connect('button_release_event', onrelease)

    def get_fit(self):
        ''' Returns an array (in this order) of peak frequency, peak flux density
            low-frequency index (slope), and high-frequency index (slope) for each
            time in self.tidx.  This is only valid for calibrated, background-
            subtracted burst emission.
        '''
        tsys, stdtsys = self.get_median_data()
        nf, nt = tsys.shape
        self.pfit = np.zeros((4,nt))
        dlogtsys = stdtsys/tsys
        for i in range(nt):
            p, ffit, sfit = tpfit(np.log(self.fghz[self.fidx[0]:self.fidx[1]]),np.log(tsys[:,i]),sigma=dlogtsys[:,i])
            # Calculate and return physical parameters S_pk, f_pk, alpha and beta
            # Note that if the fit does not work, S_pk and f_pk are determined from the
            # data, not the fit, and the slopes are nan.
            S_pk = np.exp(sfit).max()
            f_pk = np.exp(ffit[sfit.argmax()])
            lf_slope = p[1]
            hf_slope = p[1]-p[3]
            self.pfit[:,i] = [f_pk, S_pk, lf_slope, hf_slope]


def peval(x, a, b, c, d):
    ''' Function called by curve_fit() to evaluate fitting function.
    '''
    # Original Staehli expression S = exp(A)*f^alpha*[1-exp(-exp(B)*f^(-beta))]
    # can be written log S = A + alpha*log(f) + log(1-exp(-exp(B-beta*log(f))))
    # Hence, in the expression below:
    #   x = log(f) [base-e]
    #   a = A
    #   b = alpha
    #   c = B
    #   d = beta
    # Physical parameters are:
    #   Low-frequency slope = b
    #   High-frequency slope = alpha - beta = b - d
    #   Analytical expressions for S_pk and f_pk are not available, so
    #     the peak frequency and flux density are determined from the
    #     fit array
    return a+b*x+np.log(1-np.exp(-np.exp(c-d*x))) 

def tpfit(X,data,sigma=None): 
    ''' Given a set of log(frequencies) X and total power log(flux density) at each 
        frequency, fit the data with the Staehli function and return the fit parameters and 
        x,y arrays of smooth, fitted data ''' 
    from scipy.optimize import curve_fit
    
    # This complains if there are nans or infs in the data, so remove them first 
    X = X[np.isfinite(data)] 
    data = data[np.isfinite(data)]
    
    if sigma is not None:
        sigma = sigma[np.isfinite(data)]

    if len(data) < 4:
        return np.zeros(4),np.array([0,1]),np.array([0,1])

    y = data
    x = np.linspace(0,2.9,20000)

    try:
        params, pcov = curve_fit(peval, X, y, sigma=sigma)
        a, b, c, d = params
        y = peval(x, a, b, c, d)
    except:
        params = [np.nan]*4

    return params, x, y
