#!/usr/bin/env python
#
# Name:
#  calwidget.py -- A GUI (widget) interface to all of the calibration analysis procedures
#  
# History:
#  2017-Dec-27  DG
#    First began coding
#  2018-Jan-14  DG
#    Finally got a mostly complete and functional version running.
#  2018-Jan-16  DG
#    Fixed a couple of small bugs, esp. so we can analyze today's data
#  2018-Jan-23  DG
#    Change saving refcal to SQL, to ALWAYS ask about changing the reference time.
#    Also, another attempt to allow analyzing a date when a following date is missing.
#  2018-02-14  DG 
#    Added brute-force coarse delay calculation, and fixed phasecal to plot against band
#    instead of frequency.  Also write SQL time on each refcal line if found in SQL.
#  2018-02-17  DG
#    Replace antenna tabs with a single one with "reusable" plots, in an
#    attempt to speed things up.  Also finally solved the "keys do not work"
#    problem!
#  2019-01-01  DG
#    Plot phases using antenna 1 as reference (although plot baseline 1-14 for antenna 1)
#  2019-01-05  DG
#    Lots of new code and changes to allow combining LO and HI receiver calibrations for
#    a combined Refcal.  Also cleaned up the user interface.
#  2019-01-06  DG
#    Fixed a bug in writing Refcal back when read directly from SQL.
#  2019-06-22  DG
#    Start to make this work with new 52-band operations
#  2019-08-20  DG
#    Added ability to select two sets of time ranges to flag in each antenna/band
#  2019-08-23  DG
#    Fixed a bug in bands plotted.
#  2020-01-06  DG
#    Change line thickness for sigma map grid for every 5th ant and 10th band.  Also
#    add button for flagging all bands higher than a selected band.
#  2020-01-13 DG
#    Fixed a bug in fitting phase slopes with new 52-band data (just needed a nanmedian 
#    at one point in fix_time_drift())
#

import matplotlib
matplotlib.use('TkAgg')
import numpy as np
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg

# implement the default mpl key bindings
from matplotlib.backend_bases import key_press_handler
from matplotlib.pylab import subplots, close
from matplotlib.ticker import MaxNLocator

import sys
from time import sleep
if sys.version_info[0] < 3:
    import Tkinter as Tk
else:
    import tkinter as Tk
import ttk
from tkMessageBox import askyesno, showerror
from util import Time, nearest_val_idx, lobe, lin_phase_fit
import cal_header as ch
from stateframe import extract
import refcal_anal as ra   #I'll try to eliminate this later...only needed for phase_diff()

import tkSimpleDialog

class MyDialog(tkSimpleDialog.Dialog):

    def body(self, master):

        Tk.Label(master, text="Enter ISO-Format Time:").grid(row=0)

        self.e1 = Tk.Entry(master, width=24, font='Helvetica 12')
        self.e1.grid(row=0, column=1)
        try:
            t_ref = self.parent.t_ref
            self.e1.insert(Tk.END, t_ref.iso[:19])
        except:
            self.e1.insert(Tk.END, Time.now.iso[:19])
        return self.e1 # initial focus

    def apply(self):
        try:
            self.tout = Time(self.e1.get())
        except:
            self.tout = None

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
        self.pc_tlframe = Tk.Frame(fphacal)
        self.pc_tlframe.pack(side=Tk.LEFT, expand=True, fill=Tk.BOTH)
        pc_trframe = Tk.Frame(fphacal)
        pc_trframe.pack(side=Tk.LEFT, expand=True, fill=Tk.BOTH)
        pc_botframe = Tk.Frame(fphacal)
        pc_botframe.pack(side=Tk.BOTTOM)
        #   Date widget
        pc_dateframe = Tk.Frame(self.pc_tlframe)
        pc_dateframe.pack(side=Tk.TOP)
        Tk.Label(pc_dateframe, text='Date:', font='Helvetica 12').pack(side=Tk.LEFT)
        date_entry = Tk.Entry(pc_dateframe, width=12, font='Helvetica 12')
        date_entry.pack(side=Tk.LEFT)
        date_entry.insert(Tk.END, Time.now().iso[:10])
        date_entry.bind('<Return>',self.use_date)
        date_entry.bind('<Up>',self.date_prev)
        date_entry.bind('<Down>',self.date_next)
        #   Scan list widget
        pc_scanframe = Tk.Frame(self.pc_tlframe)
        pc_scanframe.pack(expand=False, fill=Tk.BOTH, side=Tk.TOP)
        self.pc_scanbox = Tk.Listbox(pc_scanframe, selectmode=Tk.SINGLE, width=39, height=30, font="Courier 10 bold")
        self.pc_scanscrl = Tk.Scrollbar(pc_scanframe,orient=Tk.VERTICAL)
        self.pc_scanscrl.pack(side=Tk.RIGHT, fill=Tk.Y)
        self.pc_scanbox.pack(side=Tk.LEFT, expand=True, fill=Tk.BOTH)
        self.pc_scanbox.config(yscrollcommand=self.pc_scanscrl.set)
        self.pc_scanscrl.config(command=self.pc_scanbox.yview)
        self.pc_scanbox.bind('<<ListboxSelect>>',self.scan_select)
        self.scan_selected = 'None'
        self.refcal_btn = Tk.Button(self.pc_tlframe, text='Analyze as Refcal', command=self.refcal_anal)
        self.refcal_btn.pack(side=Tk.TOP)
        self.refcal_btn.configure(state=Tk.DISABLED)
        rsframe = Tk.Frame(self.pc_tlframe)
        rsframe.pack(expand=False, fill=Tk.X, side=Tk.TOP)
        self.refcalset_btn = Tk.Button(rsframe, text='Set as Refcal', command=self.refcal_set)
        self.refcalset_btn.pack(side=Tk.LEFT, expand=1, fill=Tk.X)
        self.refcalset_btn.configure(state=Tk.DISABLED)
        self.extselect = Tk.BooleanVar()
        self.extselect_button = Tk.Checkbutton(rsframe, text="Extend Selection",
                variable=self.extselect, command=self.set_multi)
        self.extselect_button.configure(state=Tk.DISABLED)
        self.extselect_button.pack(side=Tk.LEFT, expand=1, fill=Tk.X)
        
        self.phacal_btn = Tk.Button(self.pc_tlframe, text='Analyze as Phasecal', command=self.phacal_anal)
        self.phacal_btn.pack(side=Tk.TOP)
        self.phacal_btn.configure(state=Tk.DISABLED)
        #self.user_time_btn = Tk.Button(self.pc_tlframe, text='Time request window', command=self.new_window)
        #self.user_time_btn.pack(side=Tk.TOP)
        self.fixdrift = Tk.BooleanVar()
        self.drift_button = Tk.Checkbutton(self.pc_tlframe, text="Fix Phase Drift vs. Time",
                variable=self.fixdrift)
        self.drift_button.configure(state=Tk.NORMAL)
        self.drift_button.pack(side=Tk.TOP)
        maxnbd = 52  # Default number of bands for initial interface

        #   Sigma map window
        pc_resultframe = Tk.Frame(pc_trframe)
        pc_resultframe.pack(side=Tk.TOP)
        self.resultvar = Tk.StringVar()
        Tk.Label(pc_resultframe, textvariable=self.resultvar, font='Helvetica 12').pack()
        self.resultvar.set('No Results Yet')
        self.ab_fig_info = subplots(1,1)
        self.ab_fig_info[0].set_size_inches(2.6,4.5,forward=True)
        ax = self.ab_fig_info[1]
        im = ax.pcolormesh(np.arange(14),np.arange(maxnbd+1),np.zeros((maxnbd,13)))
        for i in range(13):
            linewidth = 0.4 if i % 5 == 0 else 0.2 # make every 5th antenna line thicker
            ax.plot([i,i],[0,maxnbd],color='white',linewidth=linewidth)
#            ax.plot([i,i],[0,maxnbd],color='white',linewidth=0.2)
        for j in range(maxnbd):
            linewidth = 0.4 if j % 10 == 0 else 0.2 # make every 10th band line thicker
            ax.plot([0,13],[j,j],color='white',linewidth=linewidth)
#            ax.plot([0,13],[j,j],color='white',linewidth=0.2)
        self.ab_text = ax.text(2, maxnbd/2, 'No scan selected', color='white')
        ax.set_xlabel('Antenna Number')
        ax.set_ylabel('Band Number')
        bbox = ax.get_position().extents
        ax.set_position([bbox[0]+0.08,bbox[1],bbox[2]-bbox[0],bbox[3]-bbox[1]])
        #ax.set_title('No Results Yet')
        canvas = FigureCanvasTkAgg(self.ab_fig_info[0], master=pc_resultframe)
        canvas.mpl_connect('button_press_event',self.ab_select)
        canvas.draw()
        canvas.get_tk_widget().pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)

        self.allants = Tk.BooleanVar()
        Tk.Checkbutton(pc_resultframe, text="Apply to all antennas",
                variable=self.allants).pack(side=Tk.TOP, expand=0, fill=Tk.BOTH)
        self.allbands = Tk.BooleanVar()
        Tk.Checkbutton(pc_resultframe, text="Apply to all bands",
                variable=self.allbands).pack(side=Tk.TOP, expand=0, fill=Tk.BOTH)
        self.higherbands = Tk.BooleanVar()
        Tk.Checkbutton(pc_resultframe, text="Apply to all bands above selected one",
                variable=self.higherbands).pack(side=Tk.TOP, expand=0, fill=Tk.BOTH)
        self.apply_flags = Tk.Button(pc_resultframe, text='Apply Time Flagging', command=self.do_flags)
        self.apply_flags.pack(side=Tk.TOP, expand=0)
        self.save2sql = Tk.Button(pc_resultframe, text='Save to SQL', command=self.do_SQL)
        self.save2sql.pack(side=Tk.TOP, expand=0)
        
        #   Antenna Notebook
        self.nb_ant = ttk.Notebook(pc_botframe)
        self.nb_ant.pack(fill='both', expand='yes')
        self.nb_ant.bind("<<NotebookTabChanged>>", self.ant_tab_event)
        self.ant_selected = 0   # Currently selected antenna (0-based index)
        fant = []      # Frame for each antenna
        self.fig_info = []  # Figure handle and axes for each antenna
        fant.append(Tk.Frame())
        self.nb_ant.add(fant[0], text='Time History')
        self.fig_info.append(subplots(1,2))
        fant.append(Tk.Frame())
        self.nb_ant.add(fant[1], text='Sum Amp')
        self.fig_info.append(subplots(2,13))
        self.fig_info[-1][0].subplots_adjust(wspace=0, left=0.08, right=0.98)
        fant.append(Tk.Frame())
        self.nb_ant.add(fant[2], text='Sum Pha')
        self.fig_info.append(subplots(2,13))
        self.fig_info[-1][0].subplots_adjust(wspace=0, left=0.08, right=0.98)
        self.fig_info[-1][0].set_size_inches(9.0,5.2)
#            bbox = self.fig_info[-1][1].get_position().extents
#            self.fig_info[-1][1].set_position([bbox[0]-0.05,bbox[1],bbox[2]-bbox[0],bbox[3]-bbox[1]+0.1])
        for i in range(3):
            canvas = FigureCanvasTkAgg(self.fig_info[i][0], master=fant[i])
            if i == 0: canvas.mpl_connect('key_press_event',self.key_event)
            canvas.draw()
            canvas.get_tk_widget().pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)
            toolbar = NavigationToolbar2TkAgg(canvas, fant[i])
            toolbar.update()
            canvas._tkcanvas.pack(side=Tk.TOP, fill=Tk.BOTH, expand=True)
        self.status = Tk.Label(self.root, text="Status:", bd=1, relief=Tk.SUNKEN, anchor=Tk.W)
        self.status.pack(side=Tk.BOTTOM, fill=Tk.X)

    def quit(self):
        ''' Have to explicitly call quit, in order to close all of the plot windows.
        '''
        close('all')
        exit()
        
    def do_flags(self):
        ''' Apply the tflags.  All this is, is a call to refcal_anal for the selected scan
        '''
        self.status.config(text = 'Status: New time flags applied.')
        out = self.pc_dictlist[self.scan_selected]
        out = refcal_anal(out)
        self.pc_dictlist[self.scan_selected] = out
        self.scan_select()
        
    def do_SQL(self):
        #print 'Selected Save to SQL button.'
        k = self.scan_selected
        rk = self.ref_selected
        data = self.pc_dictlist[k]
        scan = self.pc_scanbox.get(k + 2)
        if scan[-1] == 'P':
            if self.saved[rk]:
                question = 'Save '+scan+' as a Phase Calibration?'
            else:
                showerror("Error",'The associated Reference Calibration must be saved first.')
                return
        elif scan[-1] == 'R':
            question = 'Save '+scan+' as a Reference Calibration scan?'
        else:
            showerror("Error",'This scan has not been analyzed.')
            return
        if self.saved[k]:
            question = 'This scan was read from SQL and has not been re-analyzed.  Save anyway?'
            
        if askyesno("Save to SQL",question):
            if scan[-1] == 'P':
                #Form correct dictionary for ch.phacal2sql()
                phacal = {'flag':data['flags'][:,:2], 'sigma':data['sigma'][:,:2], 'fghz':data['fghz'], 
                           'amp':np.abs(data['x'][:,:2]), 'pha':np.angle(data['x'][:,:2])}
                t_pha = Time(data['times'][0],format='jd')
                t_ed = Time(data['times'][-1],format='jd')
                phcal = {'phacal':phacal, 't_pha':t_pha, 't_bg':t_pha, 't_ed':t_ed, 
                           'poff':data['offsets'], 'pslope':data['mbd'], 'flag':data['mbd_flag']}
                ch.phacal2sql(phcal)
                self.status.config(text = 'Status: Phase Calibration saved to SQL Database at '+t_pha.iso)
                self.saved[k] = True
            elif scan[-1] == 'R':
                #Form correct dictionary for ch.refcal2sql()
                try:
                    t_ref = Time(data['times'][0],format='jd')
                    t_ed = Time(data['times'][-1],format='jd')
                    t_bg = t_ref
                except:
                    try:
                        t_bg = data['T_beg']
                        t_ed = data['T_end']
                        t_ref = t_bg
                    except:
                        showerror("Error",'Unknown time format.')
                        return
                rfcal = {'timestamp':t_ref, 't_bg':t_bg, 't_ed':t_ed, 'flag':data['flags'][:,:2],
                        'vis':data['x'][:,:2], 'sigma':data['sigma'][:,:2], 'fghz':data['fghz']}
                timestamp = t_ref
                question = 'Do you want to override SQL time '+t_ref.iso+'?'
                if askyesno("Override SQL Time",question):
                    self.root.t_ref = t_ref
                    d = MyDialog(self.root, title='Enter New SQL Time')
                    try:
                        if d.tout is None:
                            showerror("Error",'Unknown time format.  Please try Save to SQL button again.')
                            return
                        timestamp = Time(d.tout)
                    except:
                        self.status.config(text = 'Status: User canceled the dialog.')
                        return
                ch.refcal2sql(rfcal,timestamp)
                self.status.config(text = 'Status: Reference Calibration saved to SQL Database at '+timestamp.iso)
                self.saved[k] = True
                # Update the line with SQL time
                lines = self.pc_scanbox.get(0,Tk.END)
                self.pc_scanbox.delete(0,Tk.END)
                for i,line in enumerate(lines):
                    if k+2 == i:
                        line = line[:8] + timestamp.iso[10:19] + line[17:]
                    self.pc_scanbox.insert(Tk.END, line)
        
    def ant_tab_event(self, event):
        '''When user selects an antenna tab, this callback allows the newly
           exposed plot to be updated for current band and antenna.
           
           Antenna numbers are 0-based.
        '''
        tab = self.nb_ant.index(self.nb_ant.select())
        fig, ax = self.fig_info[tab]
        fig.canvas.draw()
        if tab == 0:
            # Update the plot on the antenna tab by "faking" an event
            try:
                event.xdata = self.ant_selected
                event.ydata = self.band_selected
                self.ab_select(event)
            except:
                pass
            
    def key_event(self, event):
        '''The user has pressed a key while the mouse is in an active
           plot window, so get information and act accordingly.
        '''
        #print 'Key',event.key,'at data coordinates',event.xdata,event.ydata
        if event.xdata is None or event.xdata < 2:
            # Indicates mouse is not in a window, or the window does not
            # contain valid times.
            self.status.config(text = 'Status: '+event.key+' ignored.  Not in window.')
            return
        if not 'tflags' in self.pc_dictlist[self.scan_selected].keys():
            self.status.config(text = 'Status: Selected scan has no time profiles (SQL scan?)')
            return
        key = event.key.upper()
        self.status.config(text = 'Status: '+key+' at data coordinates '+str(event.xdata)+' '+str(event.ydata))
        if key in ['A','B','X']:
            # This is a valid key, so act accordingly
            ant = self.ant_selected
            band = self.band_selected
            fig, ax = self.fig_info[0]
            if event.inaxes in ax:
                nlines = len(event.inaxes.lines)
                tflags = self.pc_dictlist[self.scan_selected]['tflags']
                # Check allants and allbands button states
                allants = self.allants.get()
                allbands = self.allbands.get()
                higherbands = self.higherbands.get()
                if key == 'A':
                    if nlines == 5:
                        # Erase last-drawn line to add a new one
                        ax[0].lines.pop()
                        ax[1].lines.pop()
                        nlines = len(event.inaxes.lines)
                    if nlines == 2 or nlines == 4:
                        # Okay to accept an "A" keystroke
                        ax[0].plot_date(event.xdata*np.ones(2),ax[0].get_ylim(),'g-')
                        ax[1].plot_date(event.xdata*np.ones(2),ax[1].get_ylim(),'g-')
                        fig.canvas.draw()
                        t1 = event.xdata
                        k = nlines/2 - 1  #  Either 0 or 1, depending on nlines
                        if allants and (allbands or higherbands):
                            if higherbands:
                                tflags[:,band:,0,k] = t1
                                tflags[:,band:,1,k] = 0
                            else:
                                tflags[:,:,0,k] = t1
                                tflags[:,:,1,k] = 0
                        elif allants:
                            tflags[:,band,0,k] = t1
                            tflags[:,band,1,k] = 0
                        elif allbands:
                            tflags[ant,:,0,k] = t1
                            tflags[ant,:,1,k] = 0
                        elif higherbands:
                            tflags[ant,band:,0,k] = t1
                            tflags[ant,band:,1,k] = 0
                        else:
                            tflags[ant,band,:,k] = [t1,0]
                elif key == 'B':
                    if nlines == 6:
                        # Erase last-drawn line to add a new one
                        ax[0].lines.pop()
                        ax[1].lines.pop()
                        nlines = len(event.inaxes.lines)
                    if nlines == 3 or nlines == 5:
                        # Okay to accept a "B" keystroke
                        ax[0].plot_date(event.xdata*np.ones(2),ax[0].get_ylim(),'r--')
                        ax[1].plot_date(event.xdata*np.ones(2),ax[1].get_ylim(),'r--')
                        self.last_key = key
                        fig.canvas.draw()
                        t2 = event.xdata
                        k = (nlines-1)/2 - 1  #  Either 0 or 1, depending on nlines
                        if allants and (allbands or higherbands):
                            if higherbands:
                                tflags[:,bands:,1,k] = t2
                            else:
                                tflags[:,:,1,k] = t2
                        elif allants:
                            tflags[:,band,1,k] = t2
                        elif allbands:
                            tflags[ant,:,1,k] = t2
                        elif higherbands:
                            tflags[ant,band:,1,k] = t2
                        else:
                            tflags[ant,band,1,k] = t2
                elif key == 'X':
                    ax[0].lines.pop()
                    ax[1].lines.pop()
                    fig.canvas.draw()
                    nlines = len(event.inaxes.lines)
                    if nlines == 2 or nlines == 4:
                        # Zero any time flags
                        k = nlines/2 - 1  #  Either 0 or 1, depending on nlines
                        if allants and (allbands or higherbands):
                            if higherbands:
                                tflags[:,:,band:,k] = 0
                            else:
                                tflags[:,:,:,k] = 0
                        elif allants:
                            tflags[:,band,:,k] = 0
                        elif allbands:
                            tflags[ant,:,:,k] = 0
                        elif higherbands:
                            tflags[ant,band:,:,k] = 0
                        else:
                            tflags[ant,band,:,k] = [0,0]
                    elif nlines == 3 or nlines == 5:
                        # Zero second time flags
                        k = (nlines-1)/2 - 1  #  Either 0 or 1, depending on nlines
                        if allants and (allbands or higherbands):
                            if higherbands:
                                tflags[:,band:,1,k] = 0
                            else:
                                tflags[:,:,1,k] = 0
                        elif allants:
                            tflags[:,band,1,k] = 0
                        elif allbands:
                            tflags[ant,:,1,k] = 0
                        elif higherbands:
                            tflags[ant,band:,1,k] = 0
                        else:
                            tflags[ant,band,1,k] = 0
                self.pc_dictlist[self.scan_selected]['tflags'] = tflags
                
    def use_date(self, event):
        ''' When user has selected a date (via <Return> in date box) this
            function verifies that the date is good, and if so, first
            checks for calibration results in the SQL database and lists
            dates and times of phase calibration observations in the scan 
            box.  For those scans with SQL calibrations, a notation '*' is
            made on the line.
        '''
        # Initialize the interface for a new date
        self.refcal_btn.configure(state=Tk.DISABLED)
        self.refcalset_btn.configure(state=Tk.DISABLED)
        self.phacal_btn.configure(state=Tk.DISABLED)
        #self.pc_resultbox.delete(0,Tk.END)
        self.pc_scanbox.delete(0,Tk.END)
        self.resultvar.set('No Scan Selected')
        self.ref_selected = None
        self.scan_selected = None
        self.ref2_selected = None
        self.scan2_selected = None
        self.band_selected = None
        sel = map(int, self.pc_scanbox.curselection())
        for s in sel:
            self.pc_scanbox.selection_clear(s)
        self.pc_scanbox.configure(selectmode=Tk.SINGLE)
        self.extselect_button.configure(state=Tk.DISABLED)
        self.extselect.set(0)
        self.refcalset_btn.configure(text='Set as Refcal')

        refcal_type = 8
        phacal_type = 9
        w = event.widget
        try:
            mjd = Time(w.get()).mjd
        except:
#            self.pc_scanbox.delete(0, Tk.END)
            self.pc_scanbox.insert(Tk.END, 'Error: Invalid Date.  Must be YYYY-MM=DD')
            return

        # At this point, check how many bands we should have (this changed on 2019-Feb-22)
        if mjd > 58536:
            self.maxnbd = 52
            from chan_util_52 import freq2bdname
        else:
            self.maxnbd = 34
            from chan_util_bc import freq2bdname
        self.freq2bdname = freq2bdname
        # Erase flags image
        fig, ax = self.ab_fig_info
        im = ax.pcolormesh(np.arange(14),np.arange(self.maxnbd+1),np.zeros((self.maxnbd,13)))
        self.ab_text.set_text('No scan selected')
        ax.set_title('')
        fig.suptitle('')
        self.ab_fig_info[0].canvas.draw()
        
        trange = Time([mjd+0.25,mjd+1.25],format='mjd')
        self.scan_dict = findscans(trange)
        sd = self.scan_dict
#        self.pc_scanbox.delete(0, Tk.END)
        if sd['msg'] != 'Success':
            self.pc_scanbox.insert(Tk.END, sd['msg'])
            return
        self.pc_scanbox.insert(Tk.END, 'Time     SQL Time Source   Duration [*]')
        self.pc_scanbox.insert(Tk.END, '-------- -------- -------- -------- ---')
        self.pc_dictlist = []
        self.saved = []
        for i in range(len(sd['Timestamp'])):
            st_time = Time(sd['Timestamp'][i],format='lv')
            en_time = Time(sd['Timestamp'][i]+sd['duration'][i]*60.,format='lv')
            line = st_time.iso[11:19] + '          ' + sd['SourceID'][i] + '{:6.1f} m  '.format(sd['duration'][i]) 
            # This scan is not a REFCAL unless proven otherwise
            not_a_refcal = True
            # This scan is not a PHACAL unless proven otherwise
            not_a_phacal = True
            # See if results exist in SQL database
            try:
                xml, buf = ch.read_cal(refcal_type, t=en_time)
                #refcal_time = Time(extract(buf,xml['Timestamp']),format='lv')  # Mid-time of data
                t_beg = Time(extract(buf,xml['T_beg']),format='lv')
                t_end = Time(extract(buf,xml['T_end']),format='lv')
                refcal_time = Time((t_beg.lv+t_end.lv)/2,format='lv')  # Mid-time of data
                mjd = refcal_time.mjd
                if mjd > 58536:
                    maxnbd = 52
                    from chan_util_52 import freq2bdname
                else:
                    maxnbd = 34
                    from chan_util_bc import freq2bdname
                dtr1 = st_time - refcal_time   # negative if in scan
                dtr2 = en_time - refcal_time   # positive if in scan
                if dtr1.jd < 0 and dtr2.jd > 0:
                    line += ' R'
                    SQL_time = Time(extract(buf,xml['SQL_timestamp']),format='lv').iso[10:19]
                    line = line.replace('          ',SQL_time+' ')
                    x = extract(buf,xml['Refcal_Real']) + 1j*extract(buf,xml['Refcal_Imag'])
                    sigma = extract(buf,xml['Refcal_Sigma'])
                    flags = extract(buf,xml['Refcal_Flag'])
                    fghz = extract(buf,xml['Fghz'])
                    bands = freq2bdname(fghz)
                    self.pc_dictlist.append({'refcal_time':refcal_time, 'T_beg': t_beg, 'T_end': t_end, 'fghz':fghz, 'sigma':sigma, 'x':x, 'flags':flags, 'bands':bands})
                    self.saved.append(True)
                    not_a_refcal = False
            except:
                pass
            if not_a_refcal:
                try:
                    xml, buf = ch.read_cal(phacal_type, t=en_time)
                    phacal_time = Time(extract(buf,xml['Timestamp']),format='lv')  # Mid-time of data
                    mjd = phacal_time.mjd
                    if mjd > 58536:
                        maxnbd = 52
                        from chan_util_52 import freq2bdname
                    else:
                        maxnbd = 34
                        from chan_util_bc import freq2bdname
                    dtp1 = st_time - phacal_time
                    dtp2 = en_time - phacal_time
                    if dtp1.jd < 0 and dtp2.jd > 0:
                        line += ' P'
                        SQL_time = Time(extract(buf,xml['SQL_timestamp']),format='lv').iso[10:19]
                        line = line.replace('          ',SQL_time+' ')
                        x = extract(buf,xml['Phacal_Amp'])*np.exp(1j*extract(buf,xml['Phacal_Pha']))
                        sigma = extract(buf,xml['Phacal_Sigma'])
                        flags = extract(buf,xml['Phacal_Flag'])
                        fghz = extract(buf,xml['Fghz'])
                        bands = freq2bdname(fghz)
                        mbd = extract(buf,xml['MBD'])
                        mbd_flag = extract(buf,xml['Flag'])
                        self.pc_dictlist.append({'fghz':fghz, 'sigma':sigma, 'x':x, 'flags':flags, 
                                             'mbd':mbd[:,:,1], 'offsets':mbd[:,:,0], 'mbd_flag':mbd_flag, 'bands':bands})
                        self.saved.append(True)
                        not_a_phacal = False
                except:
                    pass
            if not_a_refcal and not_a_phacal:
                # Neither refcal nor phacal exists for this time, so set empty dictionary
                self.pc_dictlist.append({})
                self.saved.append(False)

            nscans = len(self.pc_dictlist)
            self.pc_scanbox.insert(Tk.END, line)
        
    def set_multi(self):
        ''' Set scanbox selection mode to multiple if set '''
        if self.extselect.get():
            self.pc_scanbox.config(selectmode=Tk.MULTIPLE)
        else:
            self.pc_scanbox.config(selectmode=Tk.SINGLE)
            
    def scan_select(self, event=None):
        ''' Get information on what scan has been selected. '''
        w = self.pc_scanbox
        sel = map(int, w.curselection())
        if len(sel) == 3:
            # Third line selected--deselecting...
            for s in sel:
                if s != self.scan_selected+2 and s!= self.scan2_selected+2:
                    self.pc_scanbox.selection_clear(s)
                    return
        if len(sel) == 2:
            # Extended select means check if these are consistent with analysis as a single refcal
            # Lines sel[0] and sel[1] selected.
            if sel[1] < 2:
                # A header line was clicked, so clear selection and ignore
                if self.scan_selected == sel[0]-2:
                    self.pc_scanbox.selection_clear(sel[1])
                elif self.scan_selected == sel[1]-2:
                    self.pc_scanbox.selection_clear(sel[0])
                return
            # Check that both lines selected are analyzed REFCAL lines
            curscan = self.scan_selected
            if curscan == sel[0]-2:
                # sel[0] is the original line so set k to sel[1]
                k = sel[1]-2
            else:
                # sel[1] is the original line so set k to sel[0]
                k = sel[0]-2
            line0 = w.get(curscan+2)
            line1 = w.get(k+2)
            if line0[-1] != 'R':
                # Somehow original line is not a REFCAL!  Clear both and start over.
                self.pc_scanbox.selection_clear(curscan+2)
                self.pc_scanbox.selection_clear(k+2)             
            elif line1[-1] != 'R':
                print line1,'is not an already analyzed REFCAL scan.'
                self.pc_scanbox.selection_clear(k+2)
                return
            else:
                # Success!  Now do something useful...
                self.scan2_selected = k
                self.extselect_button.configure(state=Tk.DISABLED)
                self.refcalset_btn.configure(text='Set as Extended Refcal')
        elif len(sel) == 1:
            if sel[0] < 2:
                # A header line was clicked, so ignore
                return
            line = w.get(sel[0])
            k = sel[0]-2
            self.scan_selected = k
            self.scan2_selected = None
            self.band_selected = None
            #self.pc_resultbox.delete(0, Tk.END)
            self.refcal_btn.configure(state=Tk.NORMAL)
            if not self.ref_selected is None: self.phacal_btn.configure(state=Tk.NORMAL)
            self.resultvar.set('Sigma Map for '+line[:19])
            fig1, ax1 = self.fig_info[-2]
            fig2, ax2 = self.fig_info[-1]
            fig, ax = self.ab_fig_info
            ax.set_title(line[:19])
            fig.suptitle('Sigma Map')
            if line[-1] == 'R':
                data = self.pc_dictlist[k]
                # Convert frequency to band
                bands = data['bands'] #(data['fghz']*2 - 1).astype(np.int)
                good, = np.where(bands != -1)
                # This is a refcal so act accordingly
                flags = data['flags'][:13,:2]
                im = ax.pcolormesh(np.arange(14),np.arange(self.maxnbd+1),np.transpose(np.sum(flags,1)))
                #for i in range(13):
                #    ax.plot([i,i],[0,34],color='white',linewidth=0.2)
                #for j in range(34):
                #    ax.plot([0,13],[j,j],color='white',linewidth=0.2)
                self.ab_text.set_text('')
                # Plot summary plots
                for i in range(13):
                    for j in range(2):
                        ax1[j,i].cla()
                        ax2[j,i].cla()
                        ax1[j,i].plot(bands[good],np.abs(data['x'][i,j,good]),'.')
                        if i == 0:
                            # Special case for antenna 1
                            phz = np.unwrap(np.angle(data['x'][i,j,good]))
                            # Set 2pi wrap so that minimum of "U" (~ 7 GHz) is near 0
                            phz -= np.round(phz[14] / (2*np.pi)) * 2*np.pi
                        else:
                            #phz = np.unwrap(lobe(np.unwrap(np.angle(data['x'][i,j,good]) - np.angle(data['x'][0,j,good]))))
                            phz = lobe(np.unwrap(np.angle(data['x'][i,j,good]) - np.angle(data['x'][0,j,good])))
                        # Set 2pi wrap so that minimum of "U" (~ 7 GHz) is near 0
                        #phz = np.unwrap(np.angle(data['x'][i,j,good]))
                        #phz -= np.round(phz[14] / (2*np.pi)) * 2*np.pi
                        ax2[j,i].plot(bands[good],phz,'.')
                        ax1[j,i].set_ylim(0,0.5)
                        ax2[j,i].set_ylim(-2*np.pi,15)
                        if i != 0:
                            lab = ax1[j,i].get_yticklabels()
                            ax1[j,i].set_yticklabels(['']*len(lab))
                            ax2[j,i].set_yticklabels(['']*len(lab))
                for i in range(13):
                    ax1[0,i].set_title('Ant '+str(i+1), fontsize=10)
                    ax2[0,i].set_title('Ant '+str(i+1), fontsize=10)
                ax1[0,0].set_ylabel('XX Amplitude')
                ax2[0,0].set_ylabel('XX Phase (rad)')
                ax1[1,0].set_ylabel('YY Amplitude')
                ax2[1,0].set_ylabel('YY Phase (rad)')
                self.refcalset_btn.configure(state=Tk.NORMAL)
                    
            elif line[-1] == 'P':
                data = self.pc_dictlist[k]
                # Convert frequency to band
                bands = data['bands'] #(data['fghz']*2 - 1).astype(np.int)
                # This is a phacal so act accordingly
                flags = data['flags'][:13,:2]
                im = ax.pcolormesh(np.arange(14),np.arange(self.maxnbd+1),np.transpose(np.sum(flags,1)))
                #for i in range(13):
                #    ax.plot([i,i],[0,34],color='white',linewidth=0.2)
                #for j in range(34):
                #    ax.plot([0,13],[j,j],color='white',linewidth=0.2)
                self.ab_text.set_text('')
                #if not 'pdiff' in data.keys():
                if self.ref_selected:
                    data = phase_diff(data,self.pc_dictlist[self.ref_selected])
                # Plot summary plots
                for i in range(13):
                    for j in range(2):
                        ax1[j,i].cla()
                        ax2[j,i].cla()
                        good, = np.where(data['flags'][i,j] == 0)
                        if len(good) > 3:
                            try:
                                ax1[j,i].plot(bands[good],np.abs(data['x'][i,j,good]),'.')
                            except:
                                print 'ant',i+1,'pol',j+1,bands.shape,data['x'].shape,good
                            try:
                                phz = np.unwrap(data['pdiff'][i,j,good])
                                # unwrap starts with phz[0], so make sure it is in the lobe nearest to 0.
                                if phz[0] < -np.pi:
                                    phz += 2*np.pi
                                if phz[0] >= np.pi:
                                    phz -= 2*np.pi
                                ax2[j,i].plot(bands[good],phz,'.')
                            except:
                                pass
                            ax2[j,i].plot(bands,data['mbd'][i,j]*2*np.pi*data['fghz'])
                            ax1[j,i].set_ylim(0,0.5)
                            ax2[j,i].set_ylim(-8,8)
                for i in range(13):
                    ax1[0,i].set_title('Ant '+str(i+1))
                    ax2[0,i].set_title('Ant '+str(i+1))
                ax1[0,0].set_ylabel('XX Amplitude')
                ax2[0,0].set_ylabel('XX Phase Diff (rad)')
                ax1[1,0].set_ylabel('YY Amplitude')
                ax2[1,0].set_ylabel('YY Phase Diff (rad)')
            else:
                ax = self.ab_fig_info[1]
                im = ax.pcolormesh(np.arange(14),np.arange(self.maxnbd+1),np.zeros((self.maxnbd,13)))
                self.ab_text.set_text('Not yet analyzed')
#                self.pc_resultbox.insert(Tk.END, 'Not yet analyzed.')
                # Clear summary plots
                for i in range(13):
                    for j in range(2):
                        ax1[j,i].cla()
                        ax2[j,i].cla()
                self.refcalset_btn.configure(state=Tk.DISABLED)
            fig1.canvas.draw()
            fig2.canvas.draw()
            fig.canvas.draw()
                    
    def ab_select(self, event):
        ''' Selects antenna and band based on click of sigma image.
        '''
        if event.xdata is None:
            pass
        else:
            #print 'Ant =',np.floor(event.xdata), 'Band =',np.floor(event.ydata)
            if self.ab_text.get_text() != '':
                # No active band map, so do nothing
                pass
            else:
                ant  = int(np.floor(event.xdata))
                band = int(np.floor(event.ydata))
                self.band_selected = band  # 0-based index
                self.ant_selected = ant    # 0-based index
                k = self.scan_selected
                scan = self.pc_dictlist[k]
                vis = scan.get('vis',None)
                fig, ax = self.fig_info[0]
                self.nb_ant.select(0)
                ax[0].cla()
                ax[1].cla()
                if vis is None:
                    # This is from SQL.  No time profiles, so do nothing
                    fig.suptitle('No time history available until the selected scan is reanalyzed')
                else:
                    pdtimes = Time(scan['times'],format='jd').plot_date
                    # Update antenna plot for this band
                    fig.suptitle('Ant '+str(ant+1)+', Band '+str(band+1))
                    ax[0].plot_date(pdtimes,np.abs(vis[ant,0,band]),'.')
                    ax[0].plot_date(pdtimes,np.abs(vis[ant,1,band]),'.')
                    ax[1].plot_date(pdtimes,np.angle(vis[ant,0,band]),'.')
                    ax[1].plot_date(pdtimes,np.angle(vis[ant,1,band]),'.')
                    datamax = np.max(np.abs(vis[ant,:2,band]))
                    datamin = np.min(np.abs(vis[ant,:2,band]))
                    ax[0].set_ylim(0.001,max(1, datamax))
                    ax[0].set_yscale('log')
                    ax[1].set_ylim(-4,4)
                    ax[0].set_ylabel('Amplitude [arb units]')
                    ax[1].set_ylabel('Phase [rad]')
                    ax[0].set_xlabel('Time [UT]')
                    ax[1].set_xlabel('Time [UT]')
                    ax[0].xaxis.set_major_locator(MaxNLocator(3))
                    ax[0].xaxis.set_major_formatter(matplotlib.dates.DateFormatter("%H:%M"))
                    ax[1].xaxis.set_major_locator(MaxNLocator(3))
                    ax[1].xaxis.set_major_formatter(matplotlib.dates.DateFormatter("%H:%M"))
                    # Apply any existing time flags
                    if 'tflags' in self.pc_dictlist[self.scan_selected].keys():
                        for k in range(2):
                            tflags = self.pc_dictlist[self.scan_selected]['tflags'][ant,band,:,k]
                            if tflags[0] != 0:
                                ax[0].plot_date(tflags[0]*np.ones(2),ax[0].get_ylim(),'g-')
                                ax[1].plot_date(tflags[0]*np.ones(2),ax[1].get_ylim(),'g-')
                            if tflags[1] != 0:
                                ax[0].plot_date(tflags[1]*np.ones(2),ax[0].get_ylim(),'r--')
                                ax[1].plot_date(tflags[1]*np.ones(2),ax[1].get_ylim(),'r--')
                fig.canvas.draw()
                #fig.canvas.show()
                fig.canvas.get_tk_widget().focus_force()
                self.nb_ant.select(0)
                #print band,'selected.'
                #print self.pc_dictlist[self.scan_selected].keys()

    def date_prev(self, event):
        w = event.widget
        date = Time(Time(w.get()).mjd - 1, format='mjd').iso[0:10]
        w.delete(0, Tk.END)
        w.insert(0, date)
        
    def date_next(self, event):
        w = event.widget
        date = Time(Time(w.get()).mjd + 1, format='mjd').iso[0:10]
        w.delete(0,Tk.END)
        w.insert(0, date)
        
    def refcal_set(self):
        ''' Indicates the unique refcal selected for use in phasecal
            analysis, by adding an asterisk on the line.  Also sets
            self.ref_selected.
        '''
        text = self.refcalset_btn.config()['text'][-1]
        lines = self.pc_scanbox.get(0,Tk.END)
        self.pc_scanbox.delete(0,Tk.END)
        i = self.scan_selected
        j = None
        if text == 'Set as Extended Refcal':
            j = self.scan2_selected
        # Clear any asterisks from lines with '*R'
        for k,line in enumerate(lines):
            if line[-2:] == '*R':
                line = line[:-2] + ' R'
            if k == i+2:
                line = line[:-2] + '*R'
            if j:
                if k == j+2: line = line[:-2] + '*R'
            self.pc_scanbox.insert(Tk.END, line)
        self.ref_selected = i
        # Reset selection cleared by above insertion
        self.pc_scanbox.selection_set(i+2)
        if j: 
            self.ref2_selected = j
            self.pc_scanbox.selection_set(j+2)
            # Combine these two scans into a single calibration
            self.combine_refcal()
        self.scan_select()
        if text == 'Set as Refcal':
            self.extselect_button.configure(state=Tk.NORMAL)
        elif text == 'Set as Extended Refcal':
            self.extselect_button.configure(state=Tk.DISABLED)

    def refcal_anal(self):
        # Do Reference Calibration analysis for currently selected line
        # Updates pc_dictlist with new refcal dictionary
        if self.scan_selected is None:
            return
        self.status.config(text = 'Status: Analyzing Refcal -- please wait.')
        self.status.update()
        i = self.scan_selected
        file = self.scan_dict['filelist'][i]
        refcal = rd_refcal(file)
        if self.fixdrift.get():
            refcal = fix_time_drift(refcal)
        out = refcal_anal(refcal)
        # Update the existing dictionary (may be empty) with the new one
        self.pc_dictlist[i].update(out)
        self.saved[i] = False  # Mark as unsaved
        lines = self.pc_scanbox.get(0,Tk.END)
        self.pc_scanbox.delete(0,Tk.END)
        for k,line in enumerate(lines):
            if k == i+2:
                if line[-1] == 'R' or line[-1] == 'P':
                    line = line[:-2]+' R'
                else:
                    line += ' R'
            self.pc_scanbox.insert(Tk.END, line)
        # Reset selection cleared by above insertion
        self.pc_scanbox.selection_set(i+2)
        self.scan_select()
        self.status.config(text = 'Status: Analysis Complete.')
        
    def combine_refcal(self):
        # Combines a low- and high-frequency receiver refcal into a single one,
        # and writes the results to each of the scans separately (so either may
        # be written to SQL)
        from copy import deepcopy

        i = self.ref_selected
        j = self.ref2_selected
        lodict = None
        hidict = None
        print np.sum(np.array(self.pc_dictlist[i]['flags'][:13,:2]).astype(int)), np.sum(np.array(self.pc_dictlist[j]['flags'][:13,:2]).astype(int))
        if np.sum(np.array(self.pc_dictlist[i]['flags'][:13,:2]).astype(int)) > 500:
            lodict = self.pc_dictlist[i]
            loscan = i
        else:
            hidict = self.pc_dictlist[i]
            hiscan = i
        if np.sum(np.array(self.pc_dictlist[j]['flags'][:13,:2]).astype(int)) > 500:
            lodict = self.pc_dictlist[j]
            loscan = j
        else:
            hidict = self.pc_dictlist[j]
            hiscan = j
        if lodict is None or hidict is None:
            print 'Selected Refcal scans do not form a LO-HI pair.'
            return
        # The LO and HI receiver dicts have been identified.  Now determine phase slope of LO
        # relative to HI, and apply to correct the LO phases
        # print 'Writing out file /tmp/calwidget_out.txt'
        # f = open('/tmp/calwidget_out.txt','wb')
        # print lodict['fghz'].shape, lodict['fghz'].dtype
        # f.write(lodict['fghz'])
        # print lodict['x'].shape, lodict['x'].dtype
        # f.write(lodict['x'])
        # print lodict['flags'].shape, lodict['flags'].dtype
        # f.write(lodict['flags'])
        # print hidict['x'].shape, hidict['x'].dtype
        # f.write(hidict['x'])
        # print hidict['flags'].shape, hidict['flags'].dtype
        # f.write(hidict['flags'])
        # f.close()
        fghz = lodict['fghz']
        plo = np.angle(lodict['x'])
        phi = np.angle(hidict['x'])
        lobands, = np.where(fghz < 3.0)
        overlap, = np.where(np.logical_and(fghz > 3.0,fghz < 6.0))
        # Simply replace hidict values with lodict values for bands with frequency < 3 GHz, and copy flags. 
        hidict['x'][:,:,lobands] = lodict['x'][:,:,lobands]
        hidict['flags'][:,:,lobands] = lodict['flags'][:,:,lobands]
        #hidict['vis'][:,:,:4] = lodict['vis'][:,:,:4]
        #hidict['flag'][:,:,:4] = lodict['flag'][:,:,:4]
        # Now apply corrections (generally small) for the relevant ants 2-13 and pols XX and YY
        for i in range(12):
            for j in range(2):
                # Calculate phase slope of difference between LO and HI receiver for common frequencies
                pcal = lobe(plo[0,j,overlap] - phi[0,j,overlap])  # Ant 1 LO - HI phases (for calibration)
                ph = lobe(plo[i+1,j,overlap] - phi[i+1,j,overlap] - pcal)  # Calibrated difference LO - HI
                pout = lin_phase_fit(fghz[overlap],ph)  # Linear fit to phase difference
                # If standard deviation of fit is > 0.7 radians (40 degrees), skip the fit (no adjustment)
                if pout[2] < 0.7:
                    # Standard deviation is low enough to apply phase slope correction
                    p = pout[[1,0]]  # Swap pout order for consistency with polyval
                    pcor = np.polyval(p,fghz[lobands])   # These are the phases to subtract from current LO receiver phases
                    hidict['x'][i+1,j,lobands] = lodict['x'][i+1,j,lobands]*(np.cos(pcor)-1j*np.sin(pcor))
                    # Loop over common frequencies
                    for ibd in overlap:
                        # Check for bad (flagged) HI receiver bands, and replace with LO
                        pcor = np.polyval(p,fghz[ibd])
                        if hidict['flags'][i+1,j,ibd]:
                            hidict['x'][i+1,j,ibd] = lodict['x'][i+1,j,ibd]*(np.cos(pcor)-1j*np.sin(pcor))
        # The HI receiver values have been corrected, so now simply replace the lodict values so that both
        # scans have the same data.  This permits either one to be written to SQL
        self.pc_dictlist[loscan] = deepcopy(hidict)
        self.pc_scanbox.configure(selectmode=Tk.SINGLE)
        self.extselect_button.configure(state=Tk.DISABLED)
        self.extselect.set(0)
        self.refcalset_btn.configure(text='Set as Refcal')
        self.pc_scanbox.selection_set(hiscan+2)
        self.scan_selected = hiscan
        self.refcal_set()
        self.status.config(text = 'Status: Refcal Combine Complete.')
                
    def phacal_anal(self):
        # Do Phase Calibration analysis for currently selected line
        # Updates pc_dictlist with new phacal dictionary
        if self.scan_selected is None:
            return
        lines = self.pc_scanbox.get(2,Tk.END)
        if self.ref_selected is None:
            self.status.config(text = 'Status: Error: Please analyze and/or select a reference calibration.')
            return
        self.status.config(text = 'Status: Analyzing Phasecal -- please wait.')
        self.status.update()
        i = self.scan_selected
        file = self.scan_dict['filelist'][i]
        phacal = rd_refcal(file)
        if self.fixdrift.get():
            phacal = fix_time_drift(phacal)
        out = refcal_anal(phacal)
        # Now calculate the phase difference wrt the appropriate refcal
        pcout = phase_diff(out,self.pc_dictlist[self.ref_selected])
        self.pc_dictlist[i].update(out)
        self.saved[i] = False  # Mark as unsaved
        self.pc_scanbox.delete(2,Tk.END)
        for k,line in enumerate(lines):
            if k == i:
                if line[-1] == 'R' or line[-1] == 'P':
                    line = line[:-2]+' P'
                else:
                    line += ' P'
            self.pc_scanbox.insert(Tk.END, line)
        # Reset selection cleared by above insertion
        self.pc_scanbox.selection_set(i+2)
        self.scan_select()
        self.status.config(text = 'Status: Analysis Complete.')
        
def phase_diff(phacal, refcal):
    ''' Finds the delay slope (phase slope is 2*pi*fghz) of the difference
        between the input phase calibration and the input reference calibration.
        Adds some keywords to the phacal dict.  This does NOT fit for a phase
        offset, but still returns offsets of zero, in case this is needed in
        the future.  The sflags keyword is different from flags, because sflags
        can be set for either missing phase calibrations or missing reference
        calibrations.  The slope values are zero for entries flagged in sflags.
        
        2018-02-14  DG 
          Added brute-force coarse delay calculation
    '''
    def mbdfunc0(fghz, mbd):
        # fghz: frequency in GHz
        # ph0 = 0: phase offset identically set to zero (not fitted)
        # mbd: multi-band delay associated with the phase_phacal - phase_refcal in ns
        return 2. * np.pi * fghz * mbd
        
    def coarse_delay(fghz,phz):
        # Do a coarse search of delays corresponding to phase errors ranging from -0.1 to 0.1
        # Returns the delay value with the minimum sigma
        #tvals = np.arange(-0.1,0.1,0.01)
        #sigma = []
        #for t in tvals:
        #    sigma.append(np.std(lobe(phz - 2*np.pi*t*fghz)))
        #return tvals[np.argmin(np.array(sigma))]
        
        # Replace old coarse_delay (commented code above) with new lin_phase_fit function
        pout = lin_phase_fit(fghz,phz)
        return pout[1]/(2*np.pi)
        
    from scipy.optimize import curve_fit
    
    fghz = phacal['fghz']
    if len(fghz) != len(refcal['fghz']):
        self.status.config(text = 'Status: Phase and Reference calibrations have different frequencies.  No action taken.')
        return phacal
    dpha = np.angle(phacal['x'][:,:2]) - np.angle(refcal['x'][:,:2])
    flags = np.logical_or(phacal['flags'][:,:2],refcal['flags'][:,:2]).astype(np.int)
    amp_pc = np.abs(phacal['x'][:,:2])
    amp_rc = np.abs(refcal['x'][:,:2])
    sigma = ((phacal['sigma'][:,:2]/amp_pc)**2. + (refcal['sigma'][:,:2]/amp_rc)**2)**0.5
    slopes = np.zeros((15,2),np.float)
    offsets = np.zeros((15,2),np.float)
    flag = np.zeros((15,2),np.float)
    for ant in range(13):
        for pol in range(2):
            good, = np.where(flags[ant,pol] == 0)
            if len(good) > 3:
                x = fghz[good]
                t = coarse_delay(x,dpha[ant,pol,good])   # Get coarse delay 
                y = np.unwrap(lobe(dpha[ant,pol,good] - 2*np.pi*t*x))  # Correct for coarse delay
                p, pcov = curve_fit(mbdfunc0, x, y, p0=[0.], sigma=sigma[ant,pol,good], absolute_sigma=False)
                slopes[ant,pol] = p + t  # Add back coarse delay
                flag[ant,pol] = 0
            else:
                flag[ant,pol] = 1
    phacal.update({'mbd':slopes, 'mbd_flag':flag, 'flags': flags, 'offsets':offsets, 'pdiff':dpha})
    return phacal
    
def findscans(trange):
    '''Identify phasecal scans from UFDB files
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
    mjd0 = int(Time(int(tstart),format='lv').mjd)
    mjdnow = int(Time.now().mjd)
    if mjd0 < mjdnow:
        # The date is a previous day, so read a second ufdb file 
        # to ensure we have the whole local day
        try:
            ufdb2 = dump_tsys.rd_ufdb(Time(int(tstart)+86400.,format='lv'))
            for key in ufdb.keys():
                ufdb.update({key: np.append(ufdb[key], ufdb2[key])})
        except:
            # No previous day, so just skip it.
            pass
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
    from copy import deepcopy
    import dbutil as db
    
    out = read_idb([file], navg=navg, quackint=quackint)

    bds = np.unique(out['band'])
    nt = len(out['time'])
    nbd = len(bds)
    mjd = Time(out['time'][0],format='jd').mjd
    if mjd > 58536:
        maxnbd = 52
        from chan_util_52 import freq2bdname
    else:
        maxnbd = 34
        from chan_util_bc import freq2bdname
    vis = np.zeros((15, 4, maxnbd, nt), dtype=complex)
    fghz = np.zeros(maxnbd)
    # average over channels within each band
    o = out['x'][bl2ord[13,:13]]
    for bd in bds:
        idx = np.where(out['band'] == bd)[0]
        fghz[bd-1] = np.nanmean(out['fghz'][idx])
        vis[:13,:,bd-1] = np.mean(o[:, :, idx], axis=2)
    # Need to apply unrot to correct for feed rotation, before returning
    xml, buf = ch.read_cal(11, Time(out['time'][0],format='jd'))
    dph = extract(buf,xml['XYphase'])
    xi_rot = extract(buf,xml['Xi_Rot'])
    freq = extract(buf,xml['FGHz'])
    freq = freq[np.where(freq !=0)]
    
    band = np.array(freq2bdname(freq))
    bds, sidx = np.unique(band, return_index=True)
    nbd = len(bds)
    eidx = np.append(sidx[1:], len(band))
    dxy = np.zeros((14, maxnbd), dtype=np.float)
    xi = np.zeros(maxnbd, dtype=np.float)
    fghz = np.zeros(maxnbd)
    # average dph and xi_rot frequencies within each band, to convert to band representation
    for b, bd in enumerate(bds):
        fghz[bd - 1] = np.nanmean(freq[sidx[b]:eidx[b]])
        xi[bd - 1] = np.nanmean(xi_rot[sidx[b]:eidx[b]])
        for a in range(14):
            dxy[a, bd - 1] = np.angle(np.sum(np.exp(1j * dph[a, sidx[b]:eidx[b]])))
    bands = np.array(freq2bdname(fghz))
    # Read parallactic angles for this scan
    trange = Time(out['time'][[0,-1]],format='jd')
    times, chi = db.get_chi(trange)
    tchi = times.jd
    t = out['time']
    if len(t) > 0:
        vis2 = deepcopy(vis)
        idx = nearest_val_idx(t, tchi)
        pa = chi[idx]  # Parallactic angle for the times of this refcal.
        pa[:, [8, 9, 10, 12]] = 0.0
        nt = len(idx)  # Number of times in this refcal
        # Apply X-Y delay phase correction
        for a in range(13):
            a1 = lobe(dxy[a] - dxy[13])
            a2 = -dxy[13] - xi
            a3 = dxy[a] - xi + np.pi
            for j in range(nt):
                vis2[a, 1, :, j] *= np.exp(1j * a1)
                vis2[a, 2, :, j] *= np.exp(1j * a2)
                vis2[a, 3, :, j] *= np.exp(1j * a3)
        for j in range(nt):
            for a in range(13):
                vis[a, 0, :, j] = vis2[a, 0, :, j] * np.cos(pa[j, a]) + vis2[a, 3, :, j] * np.sin(pa[j, a])
                vis[a, 2, :, j] = vis2[a, 2, :, j] * np.cos(pa[j, a]) + vis2[a, 1, :, j] * np.sin(pa[j, a])
                vis[a, 3, :, j] = vis2[a, 3, :, j] * np.cos(pa[j, a]) - vis2[a, 0, :, j] * np.sin(pa[j, a])
                vis[a, 1, :, j] = vis2[a, 1, :, j] * np.cos(pa[j, a]) - vis2[a, 2, :, j] * np.sin(pa[j, a])
    # *******
    if fghz[1] < 1.:
        fghz[1] = 1.9290   # This band is missing, but no need to set its frequency to zero...
    return {'file': file, 'source': out['source'], 'vis': vis, 'bands': bands, 'fghz': fghz, 'times': out['time'], 'ha': out['ha'], 'dec': out['dec'], 'flag': np.zeros_like(vis, dtype=np.int)}
    
def refcal_anal(out):
    ''' Analyze the visibility data from rd_refcal and return time averaged visibility 
        values and flags.

    '''
    from copy import deepcopy
    vis = deepcopy(out['vis'])
    vis[np.where(out['flag'] == 1)] = np.nan
    times = out['times']
    mjd = Time(out['times'][0],format='jd').mjd
    if mjd > 58536:
        maxnbd = 52
    else:
        maxnbd = 34
    # Apply tflags
    if 'tflags' in out.keys():
        tflags = out['tflags']
        # Apply time flags
        for i in range(13):
            for j in range(maxnbd):
                for k in range(2):
                    if tflags[i,j,1,k] != 0.0:
                        jdrange = Time(tflags[i,j,:,k],format='plot_date').jd
                        if jdrange[0] > jdrange[1]:
                            jdrange = jdrange[::-1]
                        bad, = np.where(np.logical_and(times >= jdrange[0],times <= jdrange[1]))
                        if len(bad) > 0:
                            vis[i,:,j,bad] = np.nan
    else:
        # If there was no tflags key, add one of all zeros.
        out.update({'tflags':np.zeros((13,maxnbd,2,2),np.float)})
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
    
def fix_time_drift(out):
    nant, npol, nband, nt = out['vis'].shape
    for iant in range(nant):
        for ipol in range(2):
            slopes = []
            for iband in range(nband):
                phz = np.angle(out['vis'][iant,ipol,iband])
                #if out['flags'][iant,ipol,iband] == 0:
                p = lin_phase_fit(out['times'],phz)
                if p[2] < 0.7:
                    slopes.append(p[1]/out['fghz'][iband])
            if len(slopes) > 0:
                dpdt = np.nanmedian(slopes)  # Radians/GHz/Day
                for iband in range(nband):
                    pfit = dpdt*out['fghz'][iband]*(out['times']-out['times'][nt/2])
                    out['vis'][iant,ipol,iband] *= np.cos(pfit)-1j*np.sin(pfit)
    return out

if __name__ == "__main__":
    
    app = App()

    Tk.mainloop()
