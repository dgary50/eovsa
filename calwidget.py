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
from util import Time, nearest_val_idx, lobe
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
        self.pc_scanbox = Tk.Listbox(pc_scanframe, selectmode=Tk.SINGLE, width=35, height=10, font="Courier 10 bold")
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
        self.refcalset_btn = Tk.Button(self.pc_tlframe, text='Set as Refcal', command=self.refcal_set)
        self.refcalset_btn.pack(side=Tk.TOP)
        self.refcalset_btn.configure(state=Tk.DISABLED)
        self.phacal_btn = Tk.Button(self.pc_tlframe, text='Analyze as Phasecal', command=self.phacal_anal)
        self.phacal_btn.pack(side=Tk.TOP)
        self.phacal_btn.configure(state=Tk.DISABLED)
        #self.user_time_btn = Tk.Button(self.pc_tlframe, text='Time request window', command=self.new_window)
        #self.user_time_btn.pack(side=Tk.TOP)

        #   Sigma map window
        pc_resultframe = Tk.Frame(pc_trframe)
        pc_resultframe.pack(side=Tk.TOP)
        self.resultvar = Tk.StringVar()
        Tk.Label(pc_resultframe, textvariable=self.resultvar, font='Helvetica 12').pack()
        self.resultvar.set('No Results Yet')
        self.ab_fig_info = subplots(1,1)
        self.ab_fig_info[0].set_size_inches(2.6,4.5,forward=True)
        ax = self.ab_fig_info[1]
        im = ax.pcolormesh(np.arange(14),np.arange(35),np.zeros((34,13)))
        for i in range(13):
            ax.plot([i,i],[0,34],color='white',linewidth=0.2)
        for j in range(34):
            ax.plot([0,13],[j,j],color='white',linewidth=0.2)
        self.ab_text = ax.text(2, 17, 'No scan selected', color='white')
        ax.set_xlabel('Antenna Number')
        ax.set_ylabel('Band Number')
        bbox = ax.get_position().extents
        ax.set_position([bbox[0]+0.08,bbox[1],bbox[2]-bbox[0],bbox[3]-bbox[1]])
        #ax.set_title('No Results Yet')
        canvas = FigureCanvasTkAgg(self.ab_fig_info[0], master=pc_resultframe)
        canvas.mpl_connect('button_press_event',self.ab_select)
        canvas.show()
        canvas.get_tk_widget().pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)

        self.allants = Tk.BooleanVar()
        Tk.Checkbutton(pc_resultframe, text="Apply to all antennas",
                variable=self.allants).pack(side=Tk.TOP, expand=0, fill=Tk.BOTH)
        self.allbands = Tk.BooleanVar()
        Tk.Checkbutton(pc_resultframe, text="Apply to all bands",
                variable=self.allbands).pack(side=Tk.TOP, expand=0, fill=Tk.BOTH)
        self.apply_flags = Tk.Button(pc_resultframe, text='Apply Time Flagging', command=self.do_flags)
        self.apply_flags.pack(side=Tk.TOP, expand=0)
        self.save2sql = Tk.Button(pc_resultframe, text='Save to SQL', command=self.do_SQL)
        self.save2sql.pack(side=Tk.TOP, expand=0)
        
        #   Antenna Notebook
        self.nb_ant = ttk.Notebook(pc_botframe)
        self.nb_ant.pack(fill='both', expand='yes')
        self.nb_ant.bind("<<NotebookTabChanged>>", self.ant_tab_event)
        self.ant_tab = None   # Currently selected antenna tab
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
                self.fig_info[-1][0].subplots_adjust(wspace=0, left=0.08, right=0.98)
            else:
                self.nb_ant.add(fant[i], text='Sum Pha')
                self.fig_info.append(subplots(2,13))
                self.fig_info[-1][0].subplots_adjust(wspace=0, left=0.08, right=0.98)
            self.fig_info[-1][0].set_size_inches(9.0,5.2)
#            bbox = self.fig_info[-1][1].get_position().extents
#            self.fig_info[-1][1].set_position([bbox[0]-0.05,bbox[1],bbox[2]-bbox[0],bbox[3]-bbox[1]+0.1])
            canvas = FigureCanvasTkAgg(self.fig_info[i][0], master=fant[i])
            canvas.show()
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
        self.status.config(text = 'Status: Applying new time flags.  Select scan again to update screen.')
        out = self.pc_dictlist[self.scan_selected]
        out = refcal_anal(out)
        self.pc_dictlist[self.scan_selected] = out
        
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
                t_ref = Time(data['times'][0],format='jd')
                t_ed = Time(data['times'][-1],format='jd')
                rfcal = {'timestamp':t_ref, 't_bg':t_ref, 't_ed':t_ed, 'flag':data['flags'][:,:2],
                        'vis':data['x'][:,:2], 'sigma':data['sigma'][:,:2], 'fghz':data['fghz']}
                timestamp = t_ref
                if (t_ref.jd % 1) > 0.33:
                    question = 'This Reference Calibration is rather late in the day. Override SQL time?'
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
        
    def ant_tab_event(self, event):
        '''When user selects an antenna tab, this callback allows the newly
           exposed plots to be activated for keyboard events.
           
           Antenna numbers are 0-based.
        '''
        ant = self.nb_ant.index(self.nb_ant.select())
        antprev = self.ant_tab
        if antprev is None:
            pass
        else:
            fig, ax = self.fig_info[antprev]
            fig.canvas.mpl_disconnect(self.key_event)
            #print 'Disconnected from tab',antprev
            fig.canvas.show()
        self.ant_tab = None
        if ant < 13:
            # Enable the callback for keyboard events.
            fig, ax = self.fig_info[ant]
            fig.canvas.mpl_connect('key_press_event',self.key_event)
            self.ant_tab = ant
            #print 'Keyboard events enabled for tab',ant
            fig.canvas.show()
        # Update the plot on this tab by "faking" an event
        try:
            event.ydata = self.band_selected
            event.xdata = ant
            self.ab_select(event)
        except:
            pass
            
    def key_event(self, event):
        '''The user has pressed a key while the mouse is in an active
           plot window, so get information and act accordingly.
        '''
        #print 'Key',event.key,'at data coordinates',event.xdata,event.ydata
        if not 'tflags' in self.pc_dictlist[self.scan_selected].keys():
            self.status.config(text = 'Status: Selected scan has no time profiles (SQL scan?)')
            return
        if event.xdata is None or event.xdata < 2:
            # Indicates mouse is not in a window, or the window does not
            # contain valid times.
            self.status.config(text = 'Status: '+event.key+' ignored.  Not in window.')
            return
        key = event.key.upper()
        if key in ['A','B','X']:
            # This is a valid key, so act accordingly
            fig, ax = self.fig_info[self.ant_tab]
            if event.inaxes in ax:
                nlines = len(event.inaxes.lines)
                tflags = self.pc_dictlist[self.scan_selected]['tflags']
                # Check allants and allbands button states
                allants = self.allants.get()
                allbands = self.allbands.get()
                if key == 'A':
                    if nlines == 3:
                        # Erase last-drawn line to add a new one
                        ax[0].lines.pop()
                        ax[1].lines.pop()
                        nlines = len(event.inaxes.lines)
                    if nlines == 2:
                        # Okay to accept an "A" keystroke
                        ax[0].plot_date(event.xdata*np.ones(2),ax[0].get_ylim(),'g-')
                        ax[1].plot_date(event.xdata*np.ones(2),ax[1].get_ylim(),'g-')
                        fig.canvas.draw()
                        t1 = event.xdata
                        if allants and allbands:
                            tflags[:,:,0] = t1
                            tflags[:,:,1] = 0
                        elif allants:
                            tflags[:,self.band_selected,0] = t1
                            tflags[:,self.band_selected,1] = 0
                        elif allbands:
                            tflags[self.ant_tab,:,0] = t1
                            tflags[self.ant_tab,:,1] = 0
                        else:
                            tflags[self.ant_tab,self.band_selected,:] = [t1,0]
                elif key == 'B':
                    if nlines == 4:
                        # Erase last-drawn line to add a new one
                        ax[0].lines.pop()
                        ax[1].lines.pop()
                        nlines = len(event.inaxes.lines)
                    if nlines == 3:
                        # Okay to accept a "B" keystroke
                        ax[0].plot_date(event.xdata*np.ones(2),ax[0].get_ylim(),'r--')
                        ax[1].plot_date(event.xdata*np.ones(2),ax[1].get_ylim(),'r--')
                        self.last_key = key
                        fig.canvas.draw()
                        t2 = event.xdata
                        if allants and allbands:
                            tflags[:,:,1] = t2
                        elif allants:
                            tflags[:,self.band_selected,1] = t2
                        elif allbands:
                            tflags[self.ant_tab,:,1] = t2
                        else:
                            tflags[self.ant_tab,self.band_selected,1] = t2
                elif key == 'X':
                    ax[0].lines.pop()
                    ax[1].lines.pop()
                    fig.canvas.draw()
                    nlines = len(event.inaxes.lines)
                    if nlines == 2:
                        # Zero any time flags
                        if allants and allbands:
                            tflags[:,:,:] = 0
                        elif allants:
                            tflags[:,self.band_selected,:] = 0
                        elif allbands:
                            tflags[self.ant_tab,:,:] = 0
                        else:
                            tflags[self.ant_tab,self.band_selected,:] = [0,0]
                    elif nlines == 3:
                        # Zero second time flags
                        if allants and allbands:
                            tflags[:,:,1] = 0
                        elif allants:
                            tflags[:,self.band_selected,1] = 0
                        elif allbands:
                            tflags[self.ant_tab,:,1] = 0
                        else:
                            tflags[self.ant_tab,self.band_selected,1] = 0
                self.pc_dictlist[self.scan_selected]['tflags'] = tflags
                
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
        self.ref_selected = None
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
        self.saved = []
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
                self.saved.append(True)
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
                    self.pc_dictlist.append({'fghz':fghz, 'sigma':sigma, 'x':x, 'flags':flags, 
                                             'mbd':mbd[:,:,1], 'offsets':mbd[:,:,0], 'mbd_flag':mbd_flag})
                    self.saved.append(True)
                else:
                    # Neither refcal nor phacal exists for this time, so set empty dictionary
                    self.pc_dictlist.append({})
                    self.saved.append(False)

            nscans = len(self.pc_dictlist)
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
            self.band_selected = None

            #self.pc_resultbox.delete(0, Tk.END)
            self.refcal_btn.configure(state=Tk.NORMAL)
            if not self.ref_selected is None: self.phacal_btn.configure(state=Tk.NORMAL)
            self.resultvar.set('Sigma Map for '+line[:19])
            fig1, ax1 = self.fig_info[13]
            fig2, ax2 = self.fig_info[14]
            fig, ax = self.ab_fig_info
            ax.set_title(line[:19])
            fig.suptitle('Sigma Map')
            if line[-1] == 'R':
                data = self.pc_dictlist[k]
                # Convert frequency to band
                bands = (data['fghz']*2 - 1).astype(np.int)
                good, = np.where(bands != -1)
                # This is a refcal so act accordingly
                flags = data['flags'][:13,:2]
                im = ax.pcolormesh(np.arange(14),np.arange(35),np.transpose(np.sum(flags,1)))
                #for i in range(13):
                #    ax.plot([i,i],[0,34],color='white',linewidth=0.2)
                #for j in range(34):
                #    ax.plot([0,13],[j,j],color='white',linewidth=0.2)
                self.ab_text.set_text('')
                # Plot summary plots
                for i in range(13):
                    for j in range(2):
                        ax1[j,i].plot(bands[good],np.abs(data['x'][i,j,good]),'.')
                        # Set 2pi wrap so that minimum of "U" (~ 7 GHz) is near 0
                        phz = np.unwrap(np.angle(data['x'][i,j,good]))
                        phz -= np.round(phz[14] / (2*np.pi)) * 2*np.pi
                        ax2[j,i].plot(bands[good],phz,'.')
                        ax1[j,i].set_ylim(0,0.5)
                        ax2[j,i].set_ylim(-2*np.pi,20)
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
                bands = (data['fghz']*2 - 1).astype(np.int)
                # This is a phacal so act accordingly
                flags = data['flags'][:13,:2]
                im = ax.pcolormesh(np.arange(14),np.arange(35),np.transpose(np.sum(flags,1)))
                #for i in range(13):
                #    ax.plot([i,i],[0,34],color='white',linewidth=0.2)
                #for j in range(34):
                #    ax.plot([0,13],[j,j],color='white',linewidth=0.2)
                self.ab_text.set_text('')
                if not 'pdiff' in data.keys():
                    if self.ref_selected:
                        data = phase_diff(data,self.pc_dictlist[self.ref_selected])
                # Plot summary plots
                for i in range(13):
                    for j in range(2):
                        ax1[j,i].cla()
                        ax2[j,i].cla()
                        good, = np.where(data['flags'][i,j] == 0)
                        if len(good) > 3:
                            ax1[j,i].plot(bands,np.abs(data['x'][i,j]),'.')
                            try:
                                phz = np.unwrap(data['pdiff'][i,j,good])
                                if phz[0] < -np.pi:
                                    phz += 2*np.pi
                                ax2[j,i].plot(data['fghz'][good],phz,'.')
                            except:
                                pass
                            ax2[j,i].plot(data['fghz'],data['mbd'][i,j]*2*np.pi*data['fghz'])
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
                im = ax.pcolormesh(np.arange(14),np.arange(35),np.zeros((34,13)))
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
        ''' Selects antenna and band based on click of flags image.
        '''
        if event.xdata is None:
            pass
        else:
            #print 'Ant =',np.floor(event.xdata), 'Band =',np.floor(event.ydata)
        
            if self.ab_text.get_text() != '':
                # No active band map, so do nothing
                pass
            else:
                band = int(np.floor(event.ydata))
                self.band_selected = band
                k = self.scan_selected
                scan = self.pc_dictlist[k]
                vis = scan.get('vis',None)
                if vis is None:
                    # This is from SQL.  No time profiles, so do nothing
                    pass
                else:
                    pdtimes = Time(scan['times'],format='jd').plot_date
                    # Update all plots for this band
                    for i in range(13):
                       fig, ax = self.fig_info[i]
                       fig.suptitle('Ant '+str(i+1)+', Band '+str(band+1))
                       ax[0].cla()
                       ax[1].cla()
                       ax[0].plot_date(pdtimes,np.abs(vis[i,0,band]),'.')
                       ax[0].plot_date(pdtimes,np.abs(vis[i,1,band]),'.')
                       ax[1].plot_date(pdtimes,np.angle(vis[i,0,band]),'.')
                       ax[1].plot_date(pdtimes,np.angle(vis[i,1,band]),'.')
                       datamax = np.max(np.abs(vis[i,:2,band]))
                       datamin = np.min(np.abs(vis[i,:2,band]))
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
                           tflags = self.pc_dictlist[self.scan_selected]['tflags'][i,band,:]
                           if tflags[0] != 0:
                               ax[0].plot_date(tflags[0]*np.ones(2),ax[0].get_ylim(),'g-')
                               ax[1].plot_date(tflags[0]*np.ones(2),ax[1].get_ylim(),'g-')
                           if tflags[1] != 0:
                               ax[0].plot_date(tflags[1]*np.ones(2),ax[0].get_ylim(),'r--')
                               ax[1].plot_date(tflags[1]*np.ones(2),ax[1].get_ylim(),'r--')
                       fig.canvas.draw()
                self.nb_ant.select(int(np.floor(event.xdata)))
                #print band,'selected.'
                #print self.pc_dictlist[self.scan_selected].keys()

    def date_prev(self, event):
        w = event.widget
        date = Time(Time(w.get()).mjd - 1, format='mjd').iso[0:10]
        w.delete(0, Tk.END)
        w.insert(0, date)
        self.refcal_btn.configure(state=Tk.DISABLED)
        self.refcalset_btn.configure(state=Tk.DISABLED)
        self.phacal_btn.configure(state=Tk.DISABLED)
        #self.pc_resultbox.delete(0,Tk.END)
        # Erase flags image
        fig, ax = self.ab_fig_info
        im = ax.pcolormesh(np.arange(14),np.arange(35),np.zeros((34,13)))
        self.ab_text.set_text('No scan selected')
        ax.set_title('')
        fig.suptitle('')
        self.ab_fig_info[0].canvas.draw()
        self.pc_scanbox.delete(0,Tk.END)
        self.resultvar.set('No Scan Selected')
        self.ref_selected = None
        self.scan_selected = None
        self.band_selected = None
        
    def date_next(self, event):
        w = event.widget
        date = Time(Time(w.get()).mjd + 1, format='mjd').iso[0:10]
        w.delete(0,Tk.END)
        w.insert(0, date)
        self.refcal_btn.configure(state=Tk.DISABLED)
        self.refcalset_btn.configure(state=Tk.DISABLED)
        self.phacal_btn.configure(state=Tk.DISABLED)
        #self.pc_resultbox.delete(0,Tk.END)
        fig, ax = self.ab_fig_info
        im = ax.pcolormesh(np.arange(14),np.arange(35),np.zeros((34,13)))
        self.ab_text.set_text('No scan selected')
        ax.set_title('')
        fig.suptitle('')
        self.ab_fig_info[0].canvas.draw()
        self.pc_scanbox.delete(0,Tk.END)
        self.resultvar.set('No Scan Selected')
        self.ref_selected = None
        self.scan_selected = None
        self.band_selected = None
        
    def refcal_set(self):
        ''' Indicates the unique refcal selected for use in phasecal
            analysis, by adding an asterisk on the line.  Also sets
            self.ref_selected.
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
        self.ref_selected = i - 2

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
        self.status.config(text = 'Status: Analysis Complete.')
        
    def phacal_anal(self):
        # Do Phase Calibration analysis for currently selected line
        # Updates pc_dictlist with new phacal dictionary
        if self.scan_selected is None:
            return
        lines = self.pc_scanbox.get(2,Tk.END)
        if self.ref_selected is None:
            self.status.config(text = 'Status: Error: Please analyze and/or select a reference calibration.')
            return
        i = self.scan_selected
        file = self.scan_dict['filelist'][i]
        phacal = rd_refcal(file)
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
        self.status.config(text = 'Status: Analysis Complete.')
        
def phase_diff(phacal, refcal):
    ''' Finds the delay slope (phase slope is 2*pi*fghz) of the difference
        between the input phase calibration and the input reference calibration.
        Adds some keywords to the phacal dict.  This does NOT fit for a phase
        offset, but still returns offsets of zero, in case this is needed in
        the future.  The sflags keyword is different from flags, because sflags
        can be set for either missing phase calibrations or missing reference
        calibrations.  The slope values are zero for entries flagged in sflags.
    '''
    def mbdfunc0(fghz, mbd):
        # fghz: frequency in GHz
        # ph0 = 0: phase offset identically set to zero (not fitted)
        # mbd: multi-band delay associated with the phase_phacal - phase_refcal in ns
        return 2. * np.pi * fghz * mbd
        
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
                y = np.unwrap(lobe(dpha[ant,pol,good]))
                p, pcov = curve_fit(mbdfunc0, x, y, p0=[0.], sigma=sigma[ant,pol,good], absolute_sigma=False)
                slopes[ant,pol] = p
                flag[ant,pol] = 0
            else:
                flag[ant,pol] = 1
    phacal.update({'mbd':slopes, 'mbd_flag':flag, 'flags': flags, 'offsets':offsets, 'pdiff':dpha})
    return phacal
    
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
    mjd0 = int(Time(int(tstart),format='lv').mjd)
    mjdnow = int(Time.now().mjd)
    if mjd0 < mjdnow:
        # The date is a previous day, so read a second ufdb file 
        # to ensure we have the whole local day
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
    from copy import deepcopy
    import chan_util_bc as cu
    import dbutil as db
    
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
    xml, buf = ch.read_cal(11, Time(out['time'][0],format='jd'))
    dph = extract(buf,xml['XYphase'])
    xi_rot = extract(buf,xml['Xi_Rot'])
    freq = extract(buf,xml['FGHz'])
    freq = freq[np.where(freq != 0)]
    band = []
    for f in freq:
        band.append(cu.freq2bdname(f))
    bds, sidx = np.unique(band, return_index=True)
    nbd = len(bds)
    eidx = np.append(sidx[1:], len(band))
    dxy = np.zeros((14, 34), dtype=np.float)
    xi = np.zeros(34, dtype=np.float)
    fghz = np.zeros(34)
    # average dph and xi_rot frequencies within each band, to convert to 34-band representation
    for b, bd in enumerate(bds):
        fghz[bd - 1] = np.nanmean(freq[sidx[b]:eidx[b]])
        xi[bd - 1] = np.nanmean(xi_rot[sidx[b]:eidx[b]])
        for a in range(14):
            dxy[a, bd - 1] = np.angle(np.sum(np.exp(1j * dph[a, sidx[b]:eidx[b]])))
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
    return {'file': file, 'source': out['source'], 'vis': vis, 'bands': bds, 'fghz': fghz, 'times': out['time'], 'ha': out['ha'], 'dec': out['dec'], 'flag': np.zeros_like(vis, dtype=np.int)}
    
def refcal_anal(out):
    ''' Analyze the visibility data from rd_refcal and return time averaged visibility 
        values and flags.

    '''
    from copy import deepcopy
    vis = deepcopy(out['vis'])
    vis[np.where(out['flag'] == 1)] = np.nan
    times = out['times']
    # Apply tflags
    if 'tflags' in out.keys():
        tflags = out['tflags']
        # Apply time flags
        for i in range(13):
            for j in range(34):
                if tflags[i,j,1] != 0.0:
                    jdrange = Time(tflags[i,j],format='plot_date').jd
                    if jdrange[0] > jdrange[1]:
                        jdrange = jdrange[::-1]
                    bad, = np.where(np.logical_and(times >= jdrange[0],times <= jdrange[1]))
                    if len(bad) > 0:
                        vis[i,:,j,bad] = np.nan
    else:
        # If there was no tflags key, add one of all zeros.
        out.update({'tflags':np.zeros((13,34,2),np.float)})
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

if __name__ == "__main__":
    
    app = App()

    Tk.mainloop()
