#!/usr/bin/env python
#
# Name:
#  calwidget.py -- A GUI (widget) interface to all of the calibration analysis procedures
#  
# History:
#  2017-Dec-27  DG
#    First began coding
#

import matplotlib
matplotlib.use('TkAgg')
import numpy as np
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg

# implement the default mpl key bindings
from matplotlib.backend_bases import key_press_handler
from matplotlib.pylab import subplots, close

import sys
if sys.version_info[0] < 3:
    import Tkinter as Tk
else:
    import tkinter as Tk
import ttk
from util import Time, nearest_val_idx, lobe
import cal_header as ch
from stateframe import extract
import refcal_anal as ra   #I'll try to eliminate this later...only needed for phase_diff()

class App():

    def __init__(self):

        self.root = Tk.Tk()
        self.root.protocol("WM_DELETE_WINDOW", self.quit)
        self.root.wm_title("Calibration Widget")
        
        tabsframe = Tk.Frame()
        tabsframe.pack(expand=1,fill=Tk.BOTH)
        # Add some tabs for different calibration types
        self.nb = ttk.Notebook(tabsframe)
        self.nb.pack(fill='both', expand='yes')
        fphacal = Tk.Frame()
        self.nb.add(fphacal, text='Phase Calibration')
        ftpcal = Tk.Frame()
        self.nb.add(ftpcal, text='Total Power Calibration')
        fgaincal = Tk.Frame()
        self.nb.add(fgaincal, text='Gain Calibration')

        # Fill in Phase Calibration Window
        pc_tlframe = Tk.Frame(fphacal)
        pc_tlframe.pack(side=Tk.LEFT, expand=True, fill=Tk.BOTH)
        pc_trframe = Tk.Frame(fphacal)
        pc_trframe.pack(side=Tk.LEFT, expand=True, fill=Tk.BOTH)
        pc_botframe = Tk.Frame(fphacal)
        pc_botframe.pack(side=Tk.BOTTOM)
        #   Date widget
        pc_dateframe = Tk.Frame(pc_tlframe)
        pc_dateframe.pack(side=Tk.TOP)
        Tk.Label(pc_dateframe, text='Date:', font='Helvetica 12').pack(side=Tk.LEFT)
        date_entry = Tk.Entry(pc_dateframe, width=12, font='Helvetica 12')
        date_entry.pack(side=Tk.LEFT)
        date_entry.insert(Tk.END, Time.now().iso[:10])
        date_entry.bind('<Return>',self.use_date)
        date_entry.bind('<Up>',self.date_prev)
        date_entry.bind('<Down>',self.date_next)
        #   Scan list widget
        pc_scanframe = Tk.Frame(pc_tlframe)
        pc_scanframe.pack(expand=False, fill=Tk.BOTH, side=Tk.TOP)
        self.pc_scanbox = Tk.Listbox(pc_scanframe, selectmode=Tk.SINGLE, width=35, height=10, font="Courier 10 bold")
        self.pc_scanscrl = Tk.Scrollbar(pc_scanframe,orient=Tk.VERTICAL)
        self.pc_scanscrl.pack(side=Tk.RIGHT, fill=Tk.Y)
        self.pc_scanbox.pack(side=Tk.LEFT, expand=True, fill=Tk.BOTH)
        self.pc_scanbox.config(yscrollcommand=self.pc_scanscrl.set)
        self.pc_scanscrl.config(command=self.pc_scanbox.yview)
        self.pc_scanbox.bind('<<ListboxSelect>>',self.scan_select)
        self.scan_selected = 'None'
        self.refcal_btn = Tk.Button(pc_tlframe, text='Analyze as Refcal', command=self.refcal_anal)
        self.refcal_btn.pack(side=Tk.TOP)
        self.refcal_btn.configure(state=Tk.DISABLED)
        self.refcalset_btn = Tk.Button(pc_tlframe, text='Set as Refcal', command=self.refcal_set)
        self.refcalset_btn.pack(side=Tk.TOP)
        self.refcalset_btn.configure(state=Tk.DISABLED)
        self.phacal_btn = Tk.Button(pc_tlframe, text='Analyze as Phasecal', command=self.phacal_anal)
        self.phacal_btn.pack(side=Tk.TOP)
        self.phacal_btn.configure(state=Tk.DISABLED)
        #   Results list widget
        pc_resultframe = Tk.Frame(pc_trframe)
        pc_resultframe.pack(side=Tk.TOP)
        self.resultvar = Tk.StringVar()
        Tk.Label(pc_resultframe, textvariable=self.resultvar, font='Helvetica 12').pack()
        self.resultvar.set('No Results Yet')
        self.pc_resultbox = Tk.Listbox(pc_trframe, selectmode=Tk.SINGLE, width=45, height=11, font="Courier 10 bold")
        self.pc_resultscrl = Tk.Scrollbar(pc_trframe,orient=Tk.VERTICAL)
        self.pc_resultscrl.pack(side=Tk.RIGHT, fill=Tk.Y)
        self.pc_resultbox.pack(side=Tk.LEFT, expand=True, fill=Tk.BOTH)
        self.pc_resultbox.config(yscrollcommand=self.pc_resultscrl.set)
        self.pc_resultscrl.config(command=self.pc_resultbox.yview)
        self.pc_resultbox.bind('<<ListboxSelect>>',self.band_select)
        #   Antenna Notebook
        self.nb_ant = ttk.Notebook(pc_botframe)
        self.nb_ant.pack(fill='both', expand='yes')
        fant = []      # Frame for each antenna
        self.fig_info = []  # Figure handle and axes for each antenna
        for i in range(15):
            fant.append(Tk.Frame())
            if i < 13:
                self.nb_ant.add(fant[i], text='Ant'+str(i+1))
                self.fig_info.append(subplots(1,2))
            elif i == 13:
                self.nb_ant.add(fant[i], text='Sum Amp')
                self.fig_info.append(subplots(2,13))
            else:
                self.nb_ant.add(fant[i], text='Sum Pha')
                self.fig_info.append(subplots(2,13))
            canvas = FigureCanvasTkAgg(self.fig_info[i][0], master=fant[i])
            canvas.show()
            canvas.get_tk_widget().pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)
            toolbar = NavigationToolbar2TkAgg(canvas, fant[i])
            toolbar.update()
            canvas._tkcanvas.pack(side=Tk.TOP, fill=Tk.BOTH, expand=True)

    def quit(self):
        ''' Have to explicitly call quit, in order to close all of the plot windows.
        '''
        close('all')
        exit()
        
    def use_date(self, event):
        ''' When user has selected a date (via <Return> in date box) this
            function verifies that the date is good, and if so, first
            checks for calibration results in the SQL database and lists
            dates and times of phase calibration observations in the scan 
            box.  For those scans with SQL calibrations, a notation '*' is
            made on the line.
        '''
        refcal_type = 8
        phacal_type = 9
        w = event.widget
        self.scan_selected = None
        try:
            mjd = Time(w.get()).mjd
        except:
            self.pc_scanbox.delete(0, Tk.END)
            self.pc_scanbox.insert(Tk.END, 'Error: Invalid Date.  Must be YYYY-MM=DD')
            return
        trange = Time([mjd+0.25,mjd+1.25],format='mjd')
        self.scan_dict = findscans(trange)
        sd = self.scan_dict
        self.pc_scanbox.delete(0, Tk.END)
        if sd['msg'] != 'Success':
            self.pc_scanbox.insert(Tk.END, sd['msg'])
            return
        self.pc_scanbox.insert(Tk.END, 'Time     Source   Duration [*]')
        self.pc_scanbox.insert(Tk.END, '-------- -------- -------- ---')
        self.pc_dictlist = []
        for i in range(len(sd['Timestamp'])):
            st_time = Time(sd['Timestamp'][i],format='lv')
            en_time = Time(sd['Timestamp'][i]+sd['duration'][i]*60.,format='lv')
            line = st_time.iso[11:19] + ' ' + sd['SourceID'][i] + '{:6.1f} m '.format(sd['duration'][i]) 
            # See if results exist in SQL database
            xml, buf = ch.read_cal(refcal_type, t=en_time)
            refcal_time = Time(extract(buf,xml['Timestamp']),format='lv')  # Mid-time of data
            dtr1 = st_time - refcal_time   # negative if in scan
            dtr2 = en_time - refcal_time   # positive if in scan
            if dtr1.jd < 0 and dtr2.jd > 0:
                line += ' R'
                x = extract(buf,xml['Refcal_Real']) + 1j*extract(buf,xml['Refcal_Imag'])
                sigma = extract(buf,xml['Refcal_Sigma'])
                flags = extract(buf,xml['Refcal_Flag'])
                fghz = extract(buf,xml['Fghz'])
                self.pc_dictlist.append({'fghz':fghz, 'sigma':sigma, 'x':x, 'flags':flags})
            else:
                xml, buf = ch.read_cal(phacal_type, t=en_time)
                phacal_time = Time(extract(buf,xml['Timestamp']),format='lv')  # Mid-time of data
                dtp1 = st_time - phacal_time
                dtp2 = en_time - phacal_time
                if dtp1.jd < 0 and dtp2.jd > 0:
                    line += ' P'
                    x = extract(buf,xml['Phacal_Amp'])*np.exp(1j*extract(buf,xml['Phacal_Pha']))
                    sigma = extract(buf,xml['Phacal_Sigma'])
                    flags = extract(buf,xml['Phacal_Flag'])
                    fghz = extract(buf,xml['Fghz'])
                    mbd = extract(buf,xml['MBD'])
                    mbd_flag = extract(buf,xml['Flag'])
                    self.pc_dictlist.append({'fghz':fghz, 'sigma':sigma, 'x':x, 'flags':flags, 'mbd':mbd,
                                             'mbd_flag':mbd_flag})
                else:
                    # Neither refcal nor phacal exists for this time, so set empty dictionary
                    self.pc_dictlist.append({})

            self.pc_scanbox.insert(Tk.END, line)
        
    def scan_select(self, event):
        ''' Get information on what scan has been selected. '''
        w = event.widget
        sel = map(int, w.curselection())
        if len(sel) == 1:
            if sel[0] < 2:
                # A header line was clicked, so ignore
                return
            line = w.get(sel[0])
            k = sel[0]-2
            self.scan_selected = k
            #self.
            self.pc_resultbox.delete(0, Tk.END)
            self.refcal_btn.configure(state=Tk.NORMAL)
            self.phacal_btn.configure(state=Tk.NORMAL)
            self.resultvar.set('Results for '+line[:19])
            fig1, ax1 = self.fig_info[13]
            fig2, ax2 = self.fig_info[14]
            if line[-1] == 'R':
                data = self.pc_dictlist[k]
                # This is a refcal so act accordingly
                nbad = sum(sum(data['flags'][:13],0),0)
                for i in range(34):
                    self.pc_resultbox.insert(Tk.END, '{:2d}'.format(int(nbad[i]))+' of 26 bad for band '+str(i+1))
                # Plot summary plots
                for i in range(13):
                    for j in range(2):
                        ax1[j,i].plot(data['fghz'],np.abs(data['x'][i,j]),'.')
                        ax2[j,i].plot(data['fghz'],np.unwrap(np.angle(data['x'][i,j])),'.')
                        ax1[j,i].set_ylim(0,1)
                        ax2[j,i].set_ylim(-20,20)
                for i in range(13):
                    ax1[0,i].set_title('Ant '+str(i+1))
                    ax2[0,i].set_title('Ant '+str(i+1))
                for j in range(2):
                    ax1[j,0].set_ylabel('Amplitude')
                    ax2[j,0].set_ylabel('Phase (rad)')
                self.refcalset_btn.configure(state=Tk.NORMAL)
            else:
                self.pc_resultbox.insert(Tk.END, 'Not yet analyzed.')
                # Clear summary plots
                for i in range(13):
                    for j in range(2):
                        ax1[j,i].cla()
                        ax2[j,i].cla()
                self.refcalset_btn.configure(state=Tk.DISABLED)
            fig1.canvas.draw()
            fig2.canvas.draw()
                    
    def band_select(self, event):
        ''' Get information on what band has been selected. '''
        w = event.widget
        sel = map(int, w.curselection())
        if len(sel) == 1:
            line = w.get(sel[0])
            band = int(line[-2:])
            # Update all plots for this band
            for i in range(13):
                fig, ax = self.fig_info[i]
                

            print band,'selected.'  

    def date_prev(self, event):
        w = event.widget
        date = Time(Time(w.get()).mjd - 1, format='mjd').iso[0:10]
        w.delete(0, Tk.END)
        w.insert(0, date)
        self.refcal_btn.configure(state=Tk.DISABLED)
        self.phacal_btn.configure(state=Tk.DISABLED)
        self.pc_resultbox.delete(0,Tk.END)
        self.pc_scanbox.delete(0,Tk.END)
        self.resultvar.set('No Scan Selected')
        
    def date_next(self, event):
        w = event.widget
        date = Time(Time(w.get()).mjd + 1, format='mjd').iso[0:10]
        w.delete(0,Tk.END)
        w.insert(0, date)
        self.refcal_btn.configure(state=Tk.DISABLED)
        self.phacal_btn.configure(state=Tk.DISABLED)
        self.pc_resultbox.delete(0,Tk.END)
        self.pc_scanbox.delete(0,Tk.END)
        self.resultvar.set('No Scan Selected')
        
    def refcal_set(self):
        ''' Indicates the unique refcal selected for use in phasecal
            analysis, by adding an asterisk on the line.
        '''
        lines = self.pc_scanbox.get(0,Tk.END)
        self.pc_scanbox.delete(0,Tk.END)
        i = self.scan_selected + 2
        # Clear any asterisks from lines with '*R'
        for k,line in enumerate(lines):
            if line[-2:] == '*R':
                line = line[:-2] + ' R'
            if k == i:
                line = line[:-2] + '*R'
            self.pc_scanbox.insert(Tk.END, line)

    def refcal_anal(self):
        # Do Reference Calibration analysis for currently selected line
        # Updates pc_dictlist with new refcal dictionary
        if self.scan_selected is None:
            return
        i = self.scan_selected
        file = self.scan_dict['filelist'][i]
        refcal = rd_refcal(file)
        out = refcal_anal(refcal)
        # Update the existing dictionary (may be empty) with the new one
        self.pc_dictlist[i].update(out)
        lines = self.pc_scanbox.get(0,Tk.END)
        self.pc_scanbox.delete(0,Tk.END)
        for k,line in enumerate(lines):
            if k == i+2:
                if line[-1] == 'R' or line[-1] == 'P':
                    line = line[:-2]+' R'
                else:
                    line += ' R'
            self.pc_scanbox.insert(Tk.END, line)
        print 'Analysis Complete.'
        
    def phacal_anal(self):
        # Do Phase Calibration analysis for currently selected line
        # Updates pc_dictlist with new phacal dictionary
        if self.scan_selected is None:
            return
        lines = self.pc_scanbox.get(0,Tk.END)
        # First identify refcal line (marked with '*R')
        refline = None
        for k,line in enumerate(lines):
            if line[-2:] == '*R':
                if refline is None:
                    refline = k - 2
                else:
                    print 'Error: Only one reference calibration allowed.'
                    return
        if refline is None:
            print 'Error: Please analyze and/or select a reference calibration.'
            return
        i = self.scan_selected
        file = self.scan_dict['filelist'][i]
        phacal = rd_refcal(file)
        out = refcal_anal(phacal)
        # Now calculate the phase difference wrt the appropriate refcal
        pcdif = ra.phase_diff(out,refcal=self.pc_dictlist[refline])
        # Update the existing dictionary (may be empty) with the new one
        out.update(pcdif)
        self.pc_dictlist[i].update(out)
        self.pc_scanbox.delete(0,Tk.END)
        for k,line in enumerate(lines):
            if k == i:
                if line[-1] == 'R' or line[-1] == 'P':
                    line = line[:-2]+' P'
                else:
                    line += ' P'
            self.pc_scanbox.insert(Tk.END, line)
        print 'Analysis Complete.'
        
def phase_diff(phacal, refcal):
    '''    'file': file, 
           'source': out['source'], 
           'vis': vis, 
           'bands': bds, 
           'fghz': fghz, 
           'times': out['time'], 
           'ha': out['ha'], 
           'dec': out['dec'], 
           'flag': np.zeros_like(vis, dtype=np.int), 
           'x':vis_median, [15,4,34] 
           'sigma':sigma, 
           'flags':flagavg}
    '''
    def mbdfunc0(fghz, mbd):
        # fghz: frequency in GHz
        # ph0 = 0: phase offset identically set to zero (not fitted)
        # mbd: multi-band delay associated with the phase_phacal - phase_refcal in ns
        return 2. * np.pi * fghz * mbd
        
    from scipy.optimize import curve_fit
    
    dpha = np.angle(phacal['x']) - np.angle(refcal['x'])
    flags = np.logical_or(phacal['flags'],refcal['flags']).astype(np.int)
    amp_pc = np.abs(phacal['x'])
    amp_rc = np.abs(refcal['x'])
    sigma = ((phacal['sigma']/amp_pc)**2. + (refcal['sigma']/amp_rc)**2)**0.5
    for ant in range(15):
        for pol in range(4):
            good, = np.where(flags[ant,pol] == 0)
            if len(good) > 3:
                x = fghz[good]
                y = np.unwrap(lobe(dpha[ant,pol,good]))
                p, pcov = curve_fit(mbdfunc0, x, y, p0=[0.], sigma=sigma[ant,pol,good], absolute_sigma=False)
                p = np.polyfit(x,y,1,w=w)
                dfit = np.std(y - np.polyval(p,x))
    
def findscans(trange):
    '''Identify phasecal scans
    '''
    import dbutil
    import dump_tsys
    tstart, tend = trange.lv.astype(int).astype(str)
    cursor = dbutil.get_cursor()
    verstr = dbutil.find_table_version(cursor, tstart, True)
    query = 'select Timestamp,Project,SourceID from hV'+verstr+'_vD1 where left(Project,8) = "PHASECAL" and Timestamp between '+tstart+' and '+tend+' order by Timestamp'
    projdict, msg = dbutil.do_query(cursor, query)
    if msg != 'Success':
        return {'msg':msg}
    if projdict == {}:
        return {'msg':'No PHASECAL scans for this day'}
    tsint = projdict['Timestamp'].astype(int)
    # Check UFDB file to get duration
    ufdb = dump_tsys.rd_ufdb(Time(int(tstart),format='lv'))
    ufdb2 = dump_tsys.rd_ufdb(Time(int(tstart)+86400.,format='lv'))
    for key in ufdb.keys():
        ufdb.update({key: np.append(ufdb[key], ufdb2[key])})
    ufdb_times = ufdb['ST_TS'].astype(float).astype(int)
    idx = nearest_val_idx(tsint,ufdb_times)
    fpath = '/data1/eovsa/fits/UDB/' + trange[0].iso[:4] + '/'
    dur = []
    file = []
    for i in idx:
        dur.append(((ufdb['EN_TS'].astype(float) - ufdb['ST_TS'].astype(float))[i])/60.)
        file.append(fpath+ufdb['FILE'][i])
    # Fix source ID to remove nulls
    srclist = np.array([str(i.replace('\x00','')) for i in projdict['SourceID']])
    return {'Timestamp': tsint, 'SourceID': srclist, 'duration': np.array(dur), 'filelist':np.array(file), 'msg': msg}
    
def rd_refcal(file, quackint=120., navg=3):
    ''' Reads a single UDB file representing a calibrator scan, and averages over the
        bands in the file
    '''
    from read_idb import read_idb, bl2ord
    out = read_idb([file], navg=navg, quackint=quackint)

    bds = np.unique(out['band'])
    nt = len(out['time'])
    nbd = len(bds)
    vis = np.zeros((15, 4, 34, nt), dtype=complex)
    fghz = np.zeros(34)
    # average over channels within each band
    o = out['x'][bl2ord[13,:13]]
    for bd in bds:
        idx = np.where(out['band'] == bd)[0]
        fghz[bd-1] = np.nanmean(out['fghz'][idx])
        vis[:13,:,bd-1] = np.mean(o[:, :, idx], axis=2)
    # Need to apply unrot to correct for feed rotation, before returning
    # *******
    return {'file': file, 'source': out['source'], 'vis': vis, 'bands': bds, 'fghz': fghz, 'times': out['time'], 'ha': out['ha'], 'dec': out['dec'], 'flag': np.zeros_like(vis, dtype=np.int)}
    
def refcal_anal(out):
    ''' Analyze the visibility data from rd_refcal and return time averaged visibility 
        values and flags.

    '''
    from copy import deepcopy
    vis = deepcopy(out['vis'])
    vis[np.where(out['flag'] == 1)] = np.nan
    times = out['times']
    nt = len(times)
    # compute standard deviation of the visibilities
    sigma = np.nanstd(vis, axis=3)
    amp = np.abs(vis)
    vis_median = np.nanmedian(vis,axis=3)
    amp_median = np.nanmedian(amp,axis=3)
    snr = amp_median / sigma
    flagavg = (snr < 1).astype(np.int)  # Will be unity where snr < 1
    flagavg[np.where(np.isnan(snr))] = 1
    out.update({'x':vis_median, 'sigma':sigma, 'flags':flagavg})
    return out

    
app = App()

Tk.mainloop()
