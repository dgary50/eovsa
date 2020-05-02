#!/usr/bin/env python
#
# History:
#   2014-Dec-09  DG
#     Started this history log.  The PCapture window was slowing taking longer
#     and longer to refresh, as more plots were added.  Now explicitly clears
#     the plot before plotting a new one.0
#   2014-Dec-13  DG
#     Added text and highlight if ND is on.
#   2015-May-02  DG
#     The antenna MJD is suddenly glitching on occasion to very large values
#     and exceeding the datetime.timedelta limit in datime().  It now can get
#     no larger than today's MJD (in Communications section).  Also made a lot
#     of style edits to adhere to PEP8 style guide.
#   2015-May-29  DG
#     Converted from using datime() to using Time() based on astropy.
#   2015-Jun-25  DG
#     Now that Ant13's solar power station is online, changed SolPwr to display
#     data from both. 
#   2015-Jul-24  LK
#     Modified to support displaying stateframe data that are elements of arrays.
#   2015-Jul-25  DG
#     Changed ROACH Status output to show X and Y delays separately.
#   2015-Aug-27  DG
#     Changed order of antennas in Pointing section to allow separation of Az-El
#     and RA-Dec headings
#   2015-Sep-11  DG
#     Attempt to eliminate crashes due to bad values in tracking modes by using
#     np.clip() on PowerSwitch, RunControl, RunMode, DataMode.
#   2015-Sep-16  DG
#     Add code to make the "AT STOW" notification work for the old antennas (9-11 and 13).
#     Also changed "yellow" tracking warning limit to 0.005 degrees.
#   2015-Oct-13  DG
#     Added CryoRX tab to display cryoreceiver part of stateframe.
#   2015-Oct-14  DG
#     Adjustments to CryoRX tab.
#   2015-Oct-18  DG
#     Added FEMB to CryoRX tab, and cleaned up code a bit.
#   2015-Oct-28  DG
#     Added display of Local Sidereal Time, and fixed typo in Antenna tab red labels.
#   2015-Nov-21  DG
#     Sort saved-plot filenames.
#   2015-Nov-29  DG
#     Squashed some bugs where a timestamp or mjd of 0 was being converted to a Time()
#     object, which resulted in an annoying warning message.
#   2015-Dec-01  DG
#     Reduce font size of listboxes when the screen height is small.
#   2015-Dec-19  DG
#     Changed LNA output to print LNA name instead of number.
#   2015-Dec-30  DG
#     Changed Outlet heading to include outlet number.
#   2016-Jan-15  DG
#     Update expected STOW position for 27-m antennas to +20 Dec
#   2016-Mar-02  DG
#     Added code to read and display last CRIO command (and error)
#   2016-Mar-17  DG
#     Changed display of FEM voltage to FEM power
#   2017-Jan-18  DG
#     Added ant 14 tracking and receiver information to CryoRX display page
#   2017-Feb-09  DG
#     Changed antlist to remove ant 15, which is no longer planned to be used.
#     Also removed ant 15 (index 14) from definition of altantindex in update_display()
#     Also expanded space for FSeqFile display, to allow for longer FSEQ filenames
#   2018-Jan-10  DG
#     Added display of control room temperature, with red background if greater than 85 F
#   2018-Aug-25  DG
#     Added remaining antennas to temperature plot, and cleaned up the code.  Also added
#     Cryo-temperature (second stage) to front page, and fixed color coding to be red only
#     if temperature is out of range. Also changed startup page size and opened Temperature
#     tab on startup.
#   2018-Nov-17  DG
#     Fixed some deprecated function calls to call the replacement routines
#   2019-Jan-16  DG
#     Added indication of solar power not updating.
#   2019-Feb-23  DG
#     Fixed some annoying string display problems that were not there for earlier version
#     of Tkinter. Also finally killed the old "pcapture" tab, which had not been used in
#     forever.
#   2019-Nov-22  DG
#     Added red (error) color to LO1A Sweep Status
#   2020-Apr-23 OG
#     Added collapsable section for Antenna Last Command, increased the starting height 
#     of the opening window to 950 pixels, and set ROACH status to be expanded by default
#   2020-May-02  DG
#     Cleaned up some garbage in the "Task" string

from Tkinter import *
from ttk import *
from tkMessageBox import *
import tkFileDialog
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, \
                                              NavigationToolbar2Tk
import pylab as plt

import time
import sys
import os
import datetime
import subprocess
import socket
import struct
import copy
import numpy as np
import stateframe as stf
from util import *
from collections import deque
from math import fmod, pi
import antenna_control as ant_ctrl
import eovsa_lst as el


class App():

    def __init__(self):

        self.accini = stf.rd_ACCfile()
        if len(sys.argv) > 1:
            # Override host and port if redirected on local host
            # (and supplied as command line arguments)
            self.accini['host'] = sys.argv[1]
            self.accini['sfport'] = int(sys.argv[2])
        print 'Setting host to ', self.accini['host']
        print 'Setting port to ', self.accini['sfport']
        
        self.root = Tk()
        self.root.protocol("WM_DELETE_WINDOW", self.quit)
        self.root.wm_title('Stateframe Display')
        self.root.minsize(500,950)
        timeframe = Frame(self.root)
        timeframe.pack()
        toolbar = Frame(self.root)

        self.logsf = BooleanVar()
        self.CB = Checkbutton(toolbar, text="Log Stateframe",
                              variable=self.logsf, command=self.log_stateframe)
        self.accini['sf_file'] = None
        self.CB.pack(side=LEFT, expand=0, fill=BOTH)

        toolbar.pack(side=TOP, fill=X)
        
        style = Style()
        style.configure('BW.TLabel', foreground='black', background='yellow')
        style.configure('BG.TLabel', foreground='black', background='#8f8',
                        relief=RIDGE, width=13, anchor=CENTER, borderwidth=3)
        style.configure('BR.TLabel', foreground='black', background='#f88',
                        relief=RIDGE, width=13, anchor=CENTER, borderwidth=3)
        self.label = Label(timeframe, text='', style='BW.TLabel',
                           font="Helvetica 16 bold")
        self.label.pack(side=LEFT)
        self.lst_label = Label(timeframe, text='', font='Helvetica 11')
        self.lst_label.pack(side=RIGHT)
        textframe = Frame(self.root)
        textframe.pack(expand=1, fill=BOTH)

        # Calculate font size options (9-pt or 10-pt) based on screen size
        headerpix = 110  # Number of "overhead" pixels needed for window display
        # Number of 10-pixel font lines available on display
        maxlines = (self.root.winfo_screenheight() - headerpix)/20
        if maxlines < 46:
            font2use = "Courier 9 bold"
        else:
            font2use = "Courier 10 bold"
            
        # Attempt to add a tab
        self.nb = Notebook(textframe)
        self.nb.pack(fill='both', expand='yes')
        fmain = Frame()
        self.nb.add(fmain, text='Main')
        fantmain = Frame()
        self.nb.add(fantmain, text='Antennas')
        self.nb_ant = Notebook(fantmain)
        self.nb_ant.pack(fill='both', expand='yes')
        fant = []      # Frame for each antenna
        sfant = []     # Frame for each status area
        fantaz = []    # Frame for Az Status info for each antenna
        fantel = []    # Frame for El Status info for each antenna
        fantcn = []    # Frame for Central Status info for each antenna
        fanttrip = []  # Frame for text widgets for trip info for each antenna
        fantctrl = []  # Frame for text widgets for controller stateframe info
        self.antlist = np.array([1,2,3,4,5,6,7,8,9,10,11,12,13,14])
        nants = len(self.antlist)
        nlab = 20
        # 2-d array (of zeros) to hold labels
        self.Lbaz = [[0 for _ in range(nlab)] for _ in range(nants)]
        self.Lbel = [[0 for _ in range(nlab)] for _ in range(nants)]
        self.Lbcn = [[0 for _ in range(nlab)] for _ in range(nants)]
        self.Lbtrip = []   # Text widget for trip info
        self.Lbctrl1 = []  # Text widgets for controller stateframe info
        self.Lbctrl2 = []  # Text widgets for controller stateframe info
        self.Lbctrl3 = []  # Text widgets for controller stateframe info
        for i, iant in enumerate(self.antlist):
            fant.append(Frame())
            self.nb_ant.add(fant[i], text='Ant'+str(iant))
            sfant.append(Frame(fant[i]))
            sfant[i].pack(side=TOP)
            fantaz.append(Frame(sfant[i], relief=GROOVE, width=15, height=21))
            fantaz[i].pack(side=LEFT)
            fantel.append(Frame(sfant[i], relief=GROOVE, width=15, height=21))
            fantel[i].pack(side=LEFT)
            fantcn.append(Frame(sfant[i], relief=GROOVE, width=15, height=21))
            fantcn[i].pack(side=LEFT)
            fanttrip.append(Frame(sfant[i], width=60))
            fanttrip[i].pack(side=LEFT)
            fantctrl.append(Frame(fant[i]))
            fantctrl[i].pack(side=TOP)
            for j in range(nlab):
                # Empty labels for now
                self.Lbaz[i][j] = Label(fantaz[i], text='', style='BG.TLabel')
                self.Lbaz[i][j].pack(side=TOP)
                self.Lbel[i][j] = Label(fantel[i], text='', style='BG.TLabel')
                self.Lbel[i][j].pack(side=TOP)
                self.Lbcn[i][j] = Label(fantcn[i], text='', style='BG.TLabel')
                self.Lbcn[i][j].pack(side=TOP)
            self.Lbtrip.append(Listbox(fanttrip[i], selectmode=NONE, width=60,
                                       height=26))
            self.Lbtrip[i].pack(side=LEFT, fill=BOTH, expand=0)
            self.Lbctrl1.append(Listbox(fantctrl[i], selectmode=NONE, width=34,
                                        height=27, font=font2use))
            self.Lbctrl1[i].pack(side=LEFT, fill=BOTH, expand=0)
            self.Lbctrl2.append(Listbox(fantctrl[i], selectmode=NONE, width=34,
                                        height=27, font=font2use))
            self.Lbctrl2[i].pack(side=LEFT, fill=BOTH, expand=0)
            self.Lbctrl3.append(Listbox(fantctrl[i], selectmode=NONE, width=34,
                                        height=27, font=font2use))
            self.Lbctrl3[i].pack(side=LEFT, fill=BOTH, expand=0)

        fplot = Frame()
        self.nb.add(fplot, text='Temps')
        self.f2plot = Frame()
        self.nb.add(self.f2plot, text='Create')
#        fpng = Frame()
#        self.nb.add(fpng, text='PCapture')
        fcryo = Frame()
        self.nb.add(fcryo, text='CryoRX')

#        # pngfile tab--create a figure named 'pcapture' so we can refer to
#        # it later
#        self.pngtime = time.time()
#        self.pngplot = plt.figure('pcapture')
#        self.canvas1 = FigureCanvasTkAgg(self.pngplot, fpng)
#        print 'draw pcapture canvas'
#        self.canvas1.draw()
#        self.canvas1.get_tk_widget().pack(side=TOP, expand=1)
#        toolbar0 = NavigationToolbar2Tk(self.canvas1, fpng)
#        toolbar0.update()
#        self.canvas1._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)

        # plot tab
        self.prevtab = None
        self.plot1 = plt.figure(1)
        self.sub_plot1 = self.plot1.add_subplot(111)
        self.sub_plot1.grid()
        self.canvas = FigureCanvasTkAgg(self.plot1, fplot)
        print 'draw temperature canvas'
        self.canvas.draw()
        print 'draw canvas done'
        self.canvas.get_tk_widget().pack(side=TOP, expand=1)
        toolbar1 = NavigationToolbar2Tk(self.canvas, fplot)
        toolbar1.update()
        self.canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)

        # multi-plot tab
        self.plot2 = plt.figure(2)
        self.sub_plot2 = self.plot2.add_subplot(111)
        self.sub_plot2.grid()
        self.canvas2 = FigureCanvasTkAgg(self.plot2, self.f2plot)
        self.canvas2.draw()
        self.canvas2.get_tk_widget().pack(side=TOP, expand=1)
        toolbar2 = NavigationToolbar2Tk(self.canvas2, self.f2plot)
        toolbar2.update()
        self.canvas2._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)

        self.selectframe = Frame(self.f2plot)
        self.selectframe.pack(side=TOP, anchor=W)
        label = Label(self.selectframe, text='Select Item to Add:')
        label.pack(side=LEFT)
        btnframe = Frame(self.f2plot)
        btnframe.pack(side=TOP, anchor=W)
        saveframe = Frame(self.f2plot)
        saveframe.pack(side=TOP, anchor=E)

        # Add label for save name
        savelabel = Label(saveframe, text='Filename:')
        savelabel.pack(side=LEFT)

        # Add text box to enter the name of the new list
        self.content = StringVar()
        self.entry = Entry(saveframe, textvariable=self.content)
        self.entry.pack(side=LEFT)

        # Add label for save name
        savelabel = Label(saveframe, text='.txt')
        savelabel.pack(side=LEFT)

        # Button to save the current plot structure and it creates a new tab
        self.saveit = Button(saveframe, text='Save', command=self.save)
        self.saveit.pack(side=LEFT)

        # Add a clear button plot button
        clearbtn = Button(btnframe, text='Clear Graph', command=self.Clear)
        clearbtn.pack(side=RIGHT)
        
        # Add delete last plot button
        self.DeleteBtn = Button(btnframe, text='Delete Selected Item',
                                command=self.delete_selected)
        self.DeleteBtn.pack(side=RIGHT)

        # Add to Plot button
        self.Add2PlotBtn = Button(btnframe, text='Add Selected Item',
                                  command=self.add2plot)
        self.Add2PlotBtn.pack(side=RIGHT)

        self.miscplotlabel = []
        self.miscplotlocator = []

        # Initialize the state to have a single top-level dropdown list
        self.keyvars = []    # List of variables containing state of dropdown
        self.menu = []       # List of dropdown objects that currently exist
        # Add dropdown hierarchy starting with position 1 (entire hierarchy)
        self.add_dropdown(1)

        self.L3 = Listbox(fmain, selectmode=SINGLE, width=102, height=60,
                          font=font2use)
        self.L3.bind('<<ListboxSelect>>', self.toggle_heading)
        self.L3.pack(side=LEFT, fill=BOTH, expand=0)
        self.sectionDisplayState = [1, 0, 1, 0, 1, 1, 1, 0]
        self.colors = {'section0': '#757', 'section1': '#979',
                       'colhead': '#cfc', 'error': '#f88', 'warn': '#ff8',
                       'na': '#ddd', 'offsets': '#feb'}

        self.cryoLB = Listbox(fcryo, selectmode=SINGLE, width=102, height=46,
                          font=font2use)
        #self.cryoLB.bind('<<ListboxSelect>>', self.toggle_heading)
        self.cryoLB.pack(side=LEFT, fill=BOTH, expand=0)
        
        t = Time.now()
        self.label.configure(text=t.iso)
        self.lst_label.configure(text='  Local Sidereal Time:  '+str(el.eovsa_lst())[:8])

        # Previous non-blank "task" command from stateframe
        self.last_task = ''
        # Previous scan state (0 => not in scan, 1 => in scan)
        # A change from 0 to 1 means retrieve scan_header from acc
        self.last_scanstate = 0

        # Start logging stateframe data right away
        if socket.gethostname() == 'helios':
            second = False
            # Check whether a stateframe display is already running.
            # Only one can log data.
            mypid = str(os.getpid())
            out = subprocess.check_output(["pidof", "python"]).split()
            for i in range(len(out)):
                if mypid != out[i]:
                    if subprocess.check_output(["ps", "-lfp", out[i]]).find(
                                                'sf_display.py') != -1:
                        second = True
            os.environ['SF_LOGDIR'] = '/data/eovsa/stateframe_logs'
            if second:
                self.CB.configure(state=DISABLED)
            else:
                self.CB.invoke()
            pidlabel = Label(toolbar, text='My PID: '+mypid)
            pidlabel.pack(side=RIGHT, fill=X)

        # Make a "deque" (pronounced "deck") object and attach to self.
        # When a record is added to right side of que of length 3600, the
        # oldest record (on the left) will "age-off"
        self.que = deque([], 3600)  # Make maximum length equal to 1 h of data

        # Dictionary structure to keep the saved plots
        self.saved_dict = {}
        self.saved_dict_labels = {}

        '''read the saved list of plots'''
        self.path = 'saved_plots'

        plot_dict = {}
        # Loop over files in 'saved_plots'
        for filename in np.sort(os.listdir(self.path)):
            full_path = os.path.join(self.path, filename)
            if os.path.isfile(full_path):
                # Put contents of file (comma-separated list) into data dictionary
                # Key is filename "stem" (without .txt extension)
                key = filename[:-4]
                with open(full_path, 'r') as my_file:
                   plot_dict[key] = my_file.read()

        # Loop over files (key is filename)
        for key in plot_dict.iterkeys():
            temp = ''
            locator = []
            temp += plot_dict[key]
            # Split the data by the commas
            # This allows me to make them into arrays and use them in a much
            # easier way
            temp = temp.split(',')
            # Remove any leading or trailing white space or control characters,
            # which allows the file to be a bit more forgiving of formatting
            for i,item in enumerate(temp):
                temp[i] = item.strip()
            self.saved_dict_labels[key] = copy.deepcopy(temp)
            # Replace each list entry in temp with a list of keywords
            for i in range(len(temp)):
                temp[i] = temp[i].split(' ')
            # Loop over plot descriptor (multiple keywords)
            for l in range(len(temp)):
                sf = self.accini['sf']
                # Loop over keywords corresponding to plot descriptor
                for x in range(len(temp[l])):
                    try:
                        # Try to interpret as an integer 
                        # (and if it is, subtract 1 to make it zero-based)
                        temp[l][x] = int(temp[l][x])-1
                    except:
                        pass
                    sf = sf[temp[l][x]]
                locator.append(sf)
                self.saved_dict[key] = locator
            
        self.text2 = []
        
        #if there are saved plots, create new tabs and plotting devices
        if os.listdir(self.path):
            self.f3plot = Frame()
            self.nb.add(self.f3plot,text='SaveList')
            self.plot3 = plt.figure(3)
            self.sub_plot3 = self.plot3.add_subplot(111)
            self.sub_plot3.grid()
            self.canvas3 = FigureCanvasTkAgg(self.plot3, self.f3plot)
            self.canvas3.draw()
            self.canvas3.get_tk_widget().pack(side=TOP, expand=1)
            toolbar3 = NavigationToolbar2Tk(self.canvas3, self.f3plot)
            toolbar3.update()
            self.canvas3._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)

            self.saved_list = Listbox(self.f3plot,selectmode=SINGLE,width=35)
            self.Sy = Scrollbar(self.f3plot,orient=VERTICAL,command=self.yview)
            self.saved_list.config(yscrollcommand=self.Sy.set)
            self.Sy.pack(side=LEFT,fill=Y)
            self.saved_list.pack(side=LEFT)

            self.DeleteBtn = Button(self.f3plot, text = 'Create New Tab', command = self.new_tab)
            self.DeleteBtn.pack(side=RIGHT)
            for i in np.sort(plot_dict.keys()):
                self.saved_list.insert(END,i)

            self.plot_number = 3
            self.extra_plots = {}

            self.current_tab = None
            self.tab_change = self.nb.tab(self.nb.select(),'text')

        # Start the clock ticking
        self.root.after(1000 - int((t.datetime.microsecond)/1000.),self.inc_time)

    def quit(self):
        exit()

    #===============================
    def toggle_heading(self,evt):
        ''' When a user clicks on a section heading, this callback senses it and
            toggles whether the contents under the section heading is visible.
            Section headings are delineated by the background color, so this
            routine simply counts how many lines with section color are above the 
            clicked point
        '''
        w = evt.widget
        index = int(w.curselection()[0])
        res = w.itemcget(index,"bg")
        # Do this only if the clicked line is a section heading
        if res == self.colors['section0'] or res == self.colors['section1']:
            section = 0
            for i in range(index):
                if w.itemcget(i,"bg") == self.colors['section0'] or w.itemcget(i,"bg") == self.colors['section1']:
                    section += 1
            self.sectionDisplayState[section] = 1 - self.sectionDisplayState[section]

    #===============================
    '''button to create a new tab from a specific saved plots.
       this routine calls an outside program called plot_creator that
       takes in 2 arguments and returns the handles for a plot and canvas.
       I save these handles in a dictionary and use the name of the plots as the key
       and hadles for plot and canvas as values.
    '''
    def new_tab(self):
        if self.saved_list.curselection():
            self.plot_number += 1
            number = self.plot_number
            idx = self.saved_list.curselection()
            name = self.saved_list.get(idx)
            self.cur_plot = None
            temp = Frame()
            self.nb.add(temp,text=name)
            sub_plot, canvas = plot_creator(temp, name, number)
            self.extra_plots[name]=[sub_plot,canvas]
            canvas.show()
            canvas.get_tk_widget().pack(side=TOP, expand=1)
            canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)
            DeleteBtn = Button(temp, text = 'Delete this tab', command = self.forget)
            DeleteBtn.pack(side=RIGHT)            

    #===============================
    def add2plot(self):
        item = self.accini['sf']
        label = ''
        add_on = 0
        fmt_string = None
        size_of_elements = 0
        counter = 0

        for keyvar in self.keyvars:
            key = keyvar.get()
            label += key+' '
            try:
                # Try to interpret item as a number, in which case it is either an
                # index into an array or pointing to specific elements in an array.
                numkey = int(key)
                if len(item) == 3 and type(item[0]) is str and type(item[1]) is int:
                    item[2].reverse()
                    fmt_string = ''.join([i for i in item[0] if not i.isdigit()])
                    size_of_elements = struct.calcsize(fmt_string)
                    add_on *= item[2][counter]
                    add_on += numkey - 1
                    counter += 1
                    item[2].reverse()
                else:
                    item = item[numkey-1]
            except:
                # Not an index into an array, so go one level deeper
                item = item[key]

        if fmt_string is not None:
            item = [fmt_string, item[1] + (add_on * size_of_elements)]

        self.miscplotlabel.append(label)
        self.miscplotlocator.append(item)

    #===============================
    def add_dropdown(self, position):
        ''' Adds dropdowns to reflect current selection within stateframe.  Position is the
            number of the dropdown to add (1,2,3...)
        '''
        # Start with entire stateframe, and find the item associated with this position
        item = self.accini['sf']
        for i, keyvar in enumerate(self.keyvars[:position-1]):
            # Go into stateframe until we reach the current layer.  The resulting "item"
            # will be the last point in the chain, and position will give which dropdown
            # it refers to (numbered 1, 2, 3...).
            key = keyvar.get()
            try:
                # Try to interpret item as a number, in which case it is an
                # index into an array
                numkey = int(key)
                item = item[numkey-1]
            except:
                # Not an index into an array, so go one level deeper
                item = item[key]
        # We are at the current layer, so see if this is the end of the chain.
        thistype = type(item)
        thislen = len(item)
        if thistype is dict:
            # This is not the end of the chain, so get keys and enter into a new dropdown
            self.keyvars.append(StringVar(name='dropdown'+str(position)))
            sf_keys = item.keys()
            sf_keys.sort()
            sf_keys = [sf_keys[0]]+sf_keys
            self.menu.append(OptionMenu(self.selectframe, self.keyvars[-1], *sf_keys))
            # Set a trace on the variable so that we get a callback when this item is selected
            self.keyvars[-1].trace_variable('w',self.callback)
            self.menu[-1].pack(side = LEFT)
            # Check whether just-created dropdown requires another dropdown
            key = self.keyvars[-1].get()
            item = item[key]
            thistype = type(item)
            thislen = len(item)
        elif thislen > 3:
            position -= 1
        if thistype is list:
            if thislen == 2:
                # This is the end of the chain (plottable value)
                if type(item[0]) is str: 
                    return
                else:
                    position -= 1

            # This may be an array of values
            if thislen == 3:
                if type(item[0]) is str and type(item[1]) is int:
                    item[2].reverse()
                    for j in item[2]:
                        self.keyvars.append(StringVar())
                        sf_keys = []
                        for i in range(j):
                            sf_keys.append(str(i+1))
                        sf_keys = [sf_keys[0]]+sf_keys
                        # self.keyvars[-1].trace_variable('w',self.callback)
                        self.menu.append(OptionMenu(self.selectframe, self.keyvars[-1], *sf_keys))

                        self.menu[-1].pack(side = LEFT)
                    item[2].reverse()
                    return
                        # position += 2
                        # self.add_dropdown(position)

            # This is an array of values, so add a "numerical dropdown"
            self.keyvars.append(StringVar())
            sf_keys = []
            for i in range(thislen):
                sf_keys.append(str(i+1))
            sf_keys = [sf_keys[0]]+sf_keys
            self.menu.append(OptionMenu(self.selectframe, self.keyvars[-1], *sf_keys))
            # No callbacks on "numerical dropdown" since the hierarchy does not
            # change when selected
            # self.keyvars[-1].trace_variable('w',self.callback)
            self.menu[-1].pack(side = LEFT)
            position += 2
            self.add_dropdown(position)
        elif thistype is dict:
            # We do need to create another dropdown
            position += 1
            self.add_dropdown(position)

    #===============================
    def callback(self, name, index, mode):
        ''' Callback when a dropdown is changed
        '''
        key = self.root.globalgetvar(name)
        # Variable names are 'dropdownN' where N is an integer position in the hierarchy
        position = int(name[-1])
        n_rm = len(self.menu) - position
        # Remove all dropdowns higher than this one
        self.remove_dropdown(n_rm)
        # See if selected item is the end of the chain (plottable value)
        item = self.accini['sf']
        for keyvar in self.keyvars:
            # Go into stateframe until we reach the current layer.  The resulting "item"
            # will be the last point in the chain
            key = keyvar.get()
            try:
                # Try to interpret item as a number, in which case it is an
                # index into an array
                numkey = int(key)
                item = item[numkey-1]
            except:
                # Not an index into an array, so go one level deeper
                item = item[key]
        if type(item) is list:
            if type(item[0]) is dict:
                # Add the dropdown hierarchy starting with the selected position
                self.add_dropdown(position+1)
            elif len(item) == 3 and type(item[0]) is str:
                self.add_dropdown(position+1)
            else:
                pass
        else:
            # Add the dropdown hierarchy starting with the selected position
            self.add_dropdown(position+1)

    #===============================
    def remove_dropdown(self,n_rm):
        ''' Remove the number of dropdowns given by n_rm
        '''
        for i in range(n_rm):
            # Pop a menu and tell GUI to forget it
            menu = self.menu.pop()
            menu.pack_forget()
            menu.destroy()
            # Pop the corresponding Stringvar and delete it
            keyvar = self.keyvars.pop()

    #===============================
    '''Routine to delete the selected plot from the list
    '''
    def delete_selected(self):
        try:
            if len(self.miscplotlocator) and len(self.miscplotlabel):
                label = ''
                for keyvar in self.keyvars:
                    key = keyvar.get()
                    label += key+' '
            index = self.miscplotlabel.index(label)
            del self.miscplotlabel[index]
            del self.miscplotlocator[index]
        except:
            pass
        if len(self.miscplotlocator) == 0: self.Clear()
    
    #===============================
    def Clear(self):
        '''routine to clear the graph
        '''
        del self.miscplotlabel[:]
        del self.miscplotlocator[:]
        self.plot2.clf()
        self.sub_plot2 = self.plot2.add_subplot(111)
        self.sub_plot2.grid()
        self.canvas2.show()

    #===============================
    def forget(self):
        ''' This routine "forgets" the currently selected tab, removing it and its contents.
        '''
        self.nb.forget(self.nb.select())

    #===============================
    def save(self):
        '''this routine saves a set of plots. The routine writes the necessary
           information into a text file and opens a new tab to graph them.
           It creates a list to display the saved plots.
        '''
        cur_tab = self.nb.index(self.nb.select()) 
        total_tabs = self.nb.index('end')
        self.content = self.entry.get()
        self.newlabel = copy.deepcopy(self.miscplotlabel)
        temp = copy.deepcopy(self.miscplotlocator)
        if self.content and len(self.miscplotlabel) > 0:
            self.saved_dict[self.content]=temp
            self.saved_dict_labels[self.content]=self.newlabel
        #self.saved_keys.append(self.content)
        self.text3 = StringVar()
        self.text3.set(self.content)
        self.text2.append(self.text3)
        if (total_tabs - 1) == cur_tab and self.content and len(self.miscplotlabel) > 0:
            self.f3plot = Frame()
            self.nb.add(self.f3plot,text='SaveList')
            self.plot3 = plt.figure(3)
            self.sub_plot3 = self.plot3.add_subplot(111)
            self.sub_plot3.grid()
            self.canvas3 = FigureCanvasTkAgg(self.plot3, self.f3plot)
            self.canvas3.draw()
            self.canvas3.get_tk_widget().pack(side=TOP, expand=1)
            toolbar3 = NavigationToolbar2Tk(self.canvas3, self.f3plot)
            toolbar3.update()
            self.canvas3._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)

            self.saved_list = Listbox(self.f3plot,selectmode=SINGLE,width=35)
            self.Sy = Scrollbar(self.f3plot,orient=VERTICAL,command=self.yview)
            self.saved_list.config(yscrollcommand=self.Sy.set)
            self.Sy.pack(side=LEFT,fill=Y)
            self.saved_list.insert(END,self.content)
            self.saved_list.pack(side=LEFT)

            self.DeleteBtn = Button(self.f3plot, text = 'Create New Tab', command = self.new_tab)
            self.DeleteBtn.pack(side=RIGHT)
                                
        else:
            if self.content and len(self.miscplotlabel) > 0:
                self.saved_list.insert(END,self.content)

        if self.content and len(self.miscplotlabel) > 0:
            dirc = self.content + ".txt"
            var = 'saved_plots'

            f = open( os.path.join(var,self.content) + '.txt' , 'w')
            temp1 = ''
            temp2 = ''
            for i in range(len(self.miscplotlabel)):
                temp1 += self.saved_dict_labels[self.content][i].strip() + ', '
            temp1 = temp1[:-2]  # Remove comma-space from end of last label
            f.write(temp1)
            f.write(temp2)
            f.close()

    def yview(self, *args):
        self.saved_list.yview(*args)

# To open a new file i delete the contents of lines, self.mjd and
# the text widget. Then i proceed to populate them again
# similar to how it was done in the beginning.
    def Open(self):
        pass
        
# New option creates a new table with predetermined content.
    def New(self):
        pass
        
# The Save button will save the file as a text file, in a folder
# specified by the user. If the file exists, the program will ask
# the user if he wants to replace the file.
    def Save(self):
        pass

    def inc_time(self):

#   pos7 = self.S7.get()
#   pos8 = self.S8.get()
        data, msg = stf.get_stateframe(self.accini)

        if msg != 'No Error':
            print msg
            data = None
        elif data:
            # Compare data version number to accini version number
            version_change = False
            version = struct.unpack_from('d',data,8)[0]   # Get stateframe version from data
            if version > 0.0 and version != self.accini['version']:
                # The version number of the stateframe data has changed, so we need to reread
                # the ACC ini file (which will read a new stateframe.xml file and give us a new
                # sf dictionary.
                self.accini = stf.rd_ACCfile()
                version_change = True
            # Only log or further process non-zero stateframe data
            if version != 0.0:
                # This should be a good stateframe, so add it to the que for plotting
                self.que.append(data)
                # If there is an open log file, write raw data to it
                f = self.accini.get('sf_file',None)   
                if f:
                    date_change = int(time.time() / 86400) > int(os.path.getctime(f.name) / 86400)
                    if version_change or date_change:
                        # Looks like the version or date has changed, so open a new file
                        self.log_stateframe()
                        f = self.accini.get('sf_file')
                    f.write(data)

        curtab = self.nb.tab(self.nb.select(),'text')
        if data:
            if curtab[0:3] == 'Ant':
                curtab = self.nb_ant.tab(self.nb_ant.select(),'text')
                iant = int(curtab[3:])
                self.update_ant(data,iant)
            elif curtab == 'Temps':
                # If we have changed tabs, autoscale the plot
                if self.prevtab != curtab: self.sub_plot1.autoscale()
                self.plottemp()
            elif curtab == 'Create':
                # If we have changed tabs, autoscale the plot
                if self.prevtab != curtab: self.sub_plot2.autoscale()
                self.handle_tab()
            elif curtab == 'SaveList':
                # If we have changed tabs, autoscale the plot
                if self.prevtab != curtab: self.sub_plot3.autoscale()
                if(self.saved_list.curselection()):
                    self.handle_Extra()
            elif curtab == 'CryoRX':
                # Cryo-Receiver tab is selected
                self.cryo_display(data)
            elif curtab in self.extra_plots:
                # If we have changed tabs, autoscale the plot
                if self.prevtab != curtab: self.extra_plots[curtab][0].autoscale()
                self.cur_tab = curtab
                self.handle_Extra_plots()
#            elif curtab == 'PCapture':
#                # Case of PCapture tab showing
#                pngfile = '/common/tmp/dppcapture.png'
#                # Get creation/modification time of file, and display it
#                # if it is a new file (and the size of the file is > 100 kB).
#                (mode, ino, dev, nlink, uid, gid, size, atime, mtime, ctime) = os.stat(pngfile)
#                if ctime > self.pngtime and size > 100000:
#                    t = time.time()
#                    self.pngtime = ctime+5.  # Make sure pngtime is > ctime
#                    pngdata = plt.imread(pngfile)
#                    # Set the figure and clear it (otherwise the plots build up and take longer
#                    # to redraw)
#                    plt.figure('pcapture')
#                    plt.clf()
#                    self.pngplot.set_size_inches(10,10,forward=True)
#                    ax = plt.Axes(self.pngplot,[0,0,1,1])
#                    ax.set_axis_off()
#                    self.pngplot.add_axes(ax)
#                    ax.imshow(pngdata,origin='upper')
#                    self.canvas1.draw()
#                    print 'Update took',time.time()-t,'seconds.'
            else:
                if version != 0.0:
                    self.update_display(data)
                else:
                    self.L3.insert(0,'Stateframe is all zero--ACC late?  Message: '+msg)
                    self.L3.itemconfig(0,background="red",foreground='white') 
        else:
            try:
                self.L3.delete(0)
            except:
                pass
            self.L3.insert(0,'Stateframe could not be read--ACC is down?  Message: '+msg)
            self.L3.itemconfig(0,background="red",foreground='white') 

        self.prevtab = curtab

        t = Time.now()
        self.label.configure(text=t.iso)
        self.lst_label.configure(text='  Local Sidereal Time:  '+str(el.eovsa_lst())[:8])
        self.root.after(1000 - int((t.datetime.microsecond)/1000.),self.inc_time)

    def log_stateframe(self):
        '''Callback for when user clicks in the Log Stateframe checkbox.
           Also called when it is time to close a log file and open a new one.
        '''
        var = self.logsf.get()
        if var:
            # We are supposed to write a stateframe log, so get the directory
            try:
                logdir = os.environ['SF_LOGDIR']
            except:
                logdir = tkFileDialog.askdirectory()
                if logdir == '':
                    # User must have cancelled, so we have to do something.  Use generic location
                    if os.name == 'posix':
                        logdir = '/tmp'
                    elif os.name == 'nt':
                        logdir = os.environ['TEMP']
                    else:
                        logdir = os.getcwd()
                os.environ['SF_LOGDIR'] = logdir
            print 'Stateframe logs will be written to',logdir
            # Create file name from date
            # Get today's date
            t = Time.now()
            # Get the current version number from accini
            v = self.accini['version']
            filename = 'sf_'+t.iso[:10].replace('-','')+'_v'+str(v)+'.log'
            logfile = logdir+os.sep+filename
            if os.path.isfile(logfile):
                # Desired file name already exists, so simply append to it
                self.accini['sf_file'] = open(logfile,'a')
            else:
                # Need to open a new file, so first check if old file is open
                if self.accini['sf_file']:
                    self.accini['sf_file'].close()
                self.accini['sf_file'] = open(logfile,'wb')
        else:
            self.accini['sf_file'].close()
            self.accini['sf_file'] = None

    def update_ant(self,data,iant):
        # iant is the actual antenna number.  iant-1 is the index of the antenna
        # i is the ordinal number of the antenna in the list (self.antlist)
        sf = self.accini['sf']
        ant = sf['Antenna'][iant-1]['Controller']
        i = int(np.where(self.antlist == iant)[0])  # Behavior of numpy changed!
        #pos = self.S[i].get()
        self.Lbtrip[i].delete(0,END)

        # Status Bits
        az = stf.extract(data,ant['AzimuthMasterStatus'])
        az_stat = AzStatus(az)
        
        # Alert if drive is currently in a Trip state
        if az_stat.bits[0] == 1:
            # Azimuth drive is tripped, so display Trip info
            tripcode = stf.extract(data,ant['AzimuthTrip0'])
            tripinfo = ant_ctrl.get_trip(tripcode)
            text = ' Azimuth Trip: '+tripinfo['Trip']+' (code = '+str(tripinfo['Tripcode'])+')'
            self.Lbtrip[i].insert(END,text) 
            self.Lbtrip[i].itemconfig(END,bg=self.colors['error'])
            self.Lbtrip[i].insert(END,' Description:')
            self.Lbtrip[i].itemconfig(END,bg=self.colors['warn'])
            self.Lbtrip[i].insert(END,'    '+tripinfo['Description'])
            self.Lbtrip[i].insert(END,' Possible Reasons:')
            self.Lbtrip[i].itemconfig(END,bg=self.colors['warn'])
            for k in range(len(tripinfo['Clues'])):
                self.Lbtrip[i].insert(END,'    '+tripinfo['Clues'][k])            
        else:
            text = ' No Azimuth Trip'
            self.Lbtrip[i].insert(END,text)

        self.Lbtrip[i].insert(END,'')
 
        #self.Laz[i].insert(END,'=========Azimuth Status=========')
        j = 0
        for text, color in az_stat.get(range(24)):
            if text:
                if color == 'red':
                    self.Lbaz[i][j].configure(text=text,style='BR.TLabel')
                else:
                    self.Lbaz[i][j].configure(text=text,style='BG.TLabel')
                j += 1
        el = stf.extract(data,ant['ElevationStatus'])
        el_stat = ElStatus(el)

        # Alert if drive is currently in a Trip state
        if el_stat.bits[0] == 1:
            # Azimuth drive is tripped, so display Trip info
            tripcode = stf.extract(data,ant['ElevationTrip0'])
            tripinfo = ant_ctrl.get_trip(tripcode)
            text = ' Elevation Trip: '+tripinfo['Trip']+' (code = '+str(tripinfo['Tripcode'])+')'
            self.Lbtrip[i].insert(END,text) 
            self.Lbtrip[i].itemconfig(END,bg=self.colors['error'])
            self.Lbtrip[i].insert(END,' Description:')
            self.Lbtrip[i].itemconfig(END,bg=self.colors['warn'])
            self.Lbtrip[i].insert(END,'    '+tripinfo['Description'])
            self.Lbtrip[i].insert(END,' Possible Reasons:')
            self.Lbtrip[i].itemconfig(END,bg=self.colors['warn'])
            for k in range(len(tripinfo['Clues'])):
                self.Lbtrip[i].insert(END,'    '+tripinfo['Clues'][k])            
        else:
            text = ' No Elevation Trip'
            self.Lbtrip[i].insert(END,text)

        #self.L[i].insert(END,'=========Elevation Status=========')
        j = 0
        for text, color in el_stat.get(range(24)):
            if text:
                if color == 'red':
                    self.Lbel[i][j].configure(text=text,style='BR.TLabel')
                else:
                    self.Lbel[i][j].configure(text=text,style='BG.TLabel')
                j += 1
        cs = stf.extract(data,ant['CentralStatus'])
        cs_stat = CenStatus(cs)
        #self.L[i].insert(END,'=========Central Status=========')
        j = 0
        for text, color in cs_stat.get(range(32)):
            if text:
                if color == 'red':
                    self.Lbcn[i][j].configure(text=text,style='BR.TLabel')
                else:
                    self.Lbcn[i][j].configure(text=text,style='BG.TLabel')
                j += 1

        self.Lbctrl1[i].delete(0,END)
        self.Lbctrl2[i].delete(0,END)
        self.Lbctrl3[i].delete(0,END)
        keys = sorted(ant.keys())
        for j,key in enumerate(keys):
            value = str(stf.extract(data,ant[key]))
            nspaces = 33 - len(key) - len(value) - 1
            text = ' '+key+' '*nspaces+value
            if j < 27:
                self.Lbctrl1[i].insert(END,text)
            elif j < 54:
                self.Lbctrl2[i].insert(END,text)
            else:
                self.Lbctrl3[i].insert(END,text)

        # Fill in Trip text boxes if the antenna indicates a trip

    #for i in range(len(self.antlist)):
        #    self.L[i].yview_moveto(pos[0])

    def update_display(self,data):

        self.L3.delete(0,END)
        dtor = pi/180.

        sf = self.accini['sf']

        pswitch = ['OFF','ON ','-ERROR--']
        rmode = ['STOP    ','POSITION','VELOCITY','---NA---','TRACK   ','-ERROR--']
        rctrl = ['STANDBY','OPERATE','-ERROR--']
        dmode = ['AZ-EL ','RA-DEC','-ERROR--']
        t = Time(stf.extract(data,sf['Timestamp'])+0.000001,format='lv')
        mjd = t.mjd
        version = stf.extract(data,sf['Version'])
        line = 'SF v'+str(version)+' Time: '+t.iso
        task = stf.extract(data,sf['Schedule']['Task']).strip('\x00').replace('\t',' ').replace('\r\n','|')
        if task == '':
            task = self.last_task
        else:
            # Clean up the task string, which has some junk that does not need to be displayed.
            task = task[:-1]   # Remove trailing character (the '|' replacing the last \r\n)
            task = ' '.join(task.split()[1:])   # Remove the first "word" (a time within the current second)
            self.last_task = task

        line += ' Source: '+'  Task: '+task
        self.L3.insert(END,line)
        if version == 0:
            self.L3.itemconfig(END,bg=self.colors['error'])

        weather = sf['Schedule']['Data']['Weather']
        new_line = 'Weather>>'
        new_line += ' Wind: ' + str(int(stf.extract(data,weather['Wind']))) + ' mph ' + '<'+str(int(stf.extract(data,weather['AvgWind'])))+'>'
        dirs = ['N ','NE','E ','SE','S ','SW','W ','NW']
        direction = stf.extract(data,weather['WindDirection'])/dtor
        # Convert compass reading to cardinal direction
        idir = int(fmod(direction+22.5,360.)/45.)
        new_line += ' from ' + dirs[idir]
        new_line += '  Temp: ' + str(stf.extract(data,weather['Temperature']))[:4] + ' F'
        pressure = stf.extract(data,weather['Pressure'])
        new_line += '  Pressure: ' + str(pressure)[:6] + ' mbar'
        cr_temp = int(stf.extract(data,sf['Schedule']['Data']['Roach'][0]['Temp.ambient'])*90./5)/10. + 32
        new_line += '     Control Room Temp:'+str(cr_temp)+' F'
        self.L3.insert(END,new_line)
        if pressure == 0.0:
            self.L3.itemconfig(END,bg=self.colors['na'])
        if cr_temp > 85.:
            self.L3.itemconfig(END,bg=self.colors['error'])

        solpwr12 = sf['Schedule']['Data']['SolarPower'][0]
        solpwr13 = sf['Schedule']['Data']['SolarPower'][1]
        new_line = 'SolPwr12>>'
        sptime = Time(stf.extract(data, solpwr12['Timestamp']),format='lv')
        nsec = (t - sptime).value*86400
        if nsec > 120.0:
            new_line += 'Last read '+sptime.iso[:19]+', '+str(int(nsec))+'s'
            charge12 = 0.0
        else:
            charge12 = stf.extract(data,solpwr12['Charge'])
            new_line += ' Charge: ' + str(charge12)+'%  Volts:'+str(stf.extract(data,solpwr12['Volts']))[:4]+'VDC  Amps:'+str(stf.extract(data,solpwr12['Amps']))[:4]+'A'
        new_line += '  SolPwr13>>'
        sptime = Time(stf.extract(data, solpwr13['Timestamp']),format='lv')
        nsec = (t - sptime).value*86400
        if nsec > 120.0:
            new_line += 'Last read '+sptime.iso[:19]+', '+str(int(nsec))+'s'
            charge13 = 0.0
        else:
            charge13 = stf.extract(data,solpwr13['Charge'])
            new_line += ' Charge: ' + str(charge13)+'%  Volts:'+str(stf.extract(data,solpwr13['Volts']))[:4]+'VDC  Amps:'+str(stf.extract(data,solpwr13['Amps']))[:4]+'A'
        self.L3.insert(END,new_line)
        if charge12 == 0.0 or charge13 == 0.0:
            self.L3.itemconfig(END,bg=self.colors['error'])

        antn = [' 1 ',' 2 ',' 3 ',' 4 ',' 5 ',' 6 ',' 7 ',' 8 ',' 9 ','10 ','11 ','12 ','13 ',' A ',' B ']
        altantindex = np.array([0,1,2,3,4,5,6,7,11,8,9,10,12,13])
        # Expected STOW coordinates for the different types of antenna.  If/when Ant 12 is changed to
        # the old mount, its values must be changed.
        azstow = np.array([180,180,180,180,180,180,180,180,0,0,0,180,0,0,0])
        elstow = np.array([88,88,88,88,88,88,88,88,37,37,37,88,37,29,29])

        azel = stf.azel_from_stateframe(sf,data)
        antindex = self.antlist - 1

        sub1 = stf.extract(data,sf['LODM']['Subarray1'])
        sub2 = stf.extract(data,sf['LODM']['Subarray2'])
        subarray1 = []
        subarray2 = []
        for i in range(16):
            subarray1.append(sub1 & (1<<i) > 0)
            subarray2.append(sub1 & (1<<i) > 0)

        expand = 'Click to Expand_'
        collapse = 'Click to Collapse '
        # ================= Section 1: Antenna Tracking ===================
        nsubarray = 0
        ndbad = 0
        ncbad = 0
        for i in altantindex:
            if subarray1[i]:
                c = sf['Antenna'][i]['Controller']
                nsubarray += 1
                ra = RA_Angle(stf.extract(data,c['RAVirtualAxis'])/10000.,'degrees')
                dec = Dec_Angle(stf.extract(data,c['DecVirtualAxis'])/10000.,'degrees')
                if abs(azel['dAzimuth'][i]) > 0.02 or abs(azel['dElevation'][i]) > 0.02: 
                    ndbad += 1
                if ra.get() == 0.0 and dec.get() == 0.0: 
                    ncbad += 1
                    ndbad += 1
        heading = 'Antenna Tracking'            
        if not self.sectionDisplayState[0]:
            line = 'Cnx Good:{:2}  Bad:{:2}  Track Good:{:2}  Bad:{:2}'.format(nsubarray-ncbad,ncbad,nsubarray-ndbad,ndbad).center(100-len(heading)-len(expand),'_')
            headline = '_'+heading+line+expand
            self.L3.insert(END,headline)
            self.L3.itemconfig(END,bg=self.colors['section0'],fg=self.colors['na'])
        else:
            headline = ' '+heading+' '*(100-len(heading)-len(collapse))+collapse
            self.L3.insert(END,headline)
            self.L3.itemconfig(END,bg=self.colors['section0'],fg='white')

            head  = 'a# ---AZ--Act--EL---  ---AZ--Req--EL---  ---AZ--Err--EL---  Az-Off-El  RA-Off-Dc  ---RA---Req---Dec---'
            head2  = 'a# ---HA--Act--Dec--  ---HA--Req--Dec--  ---HA--Err--Dec--  Az-Off-El  RA-Off-Dc  ---RA---Req---Dec---'
            templ = '   xxx.xxxx sxx.xxxx  xxx.xxxx sxx.xxxx  sxx.xxxx sxx.xxxx  xx.x xx.x  xx.x xx.x  xx xx xx.x sxx xx xx'
            self.L3.insert(END,head)
            self.L3.itemconfig(END,bg=self.colors['colhead'])

            for i in altantindex:
                c = sf['Antenna'][i]['Controller']
                fe = sf['Antenna'][i]['Frontend']
                # If bit 9 of central status is set, dish is going to STOW
                tostow = stf.extract(data,c['CentralStatus']) & 256
                # If dish is in POSITION mode and the requested 
                # AZ, EL are 225, 15, dish is going to SERVICE
                toservice = (int(azel['RequestedAzimuth'][i] + 0.5) == 225 and 
                             int(azel['RequestedElevation'][i] + 0.5) == 15 and
                             rmode[stf.extract(data,c['RunMode'])] == 'POSITION')
                # If dish is at AZ, EL 180, 88, dish is AT STOW 
                atstow = (int(azel['ActualAzimuth'][i] + 0.5) == azstow[i] and 
                             int(azel['ActualElevation'][i] + 0.5) == elstow[i])
                # If dish is at AZ, EL 225, 15, dish is AT SERVICE
                atservice = (int(azel['ActualAzimuth'][i] + 0.5) == 225 and 
                             int(azel['ActualElevation'][i] + 0.5) == 15)
                line = antn[i]+'{:8.4f}'.format(azel['ActualAzimuth'][i])
                line += " "+'{:+8.4f}'.format(azel['ActualElevation'][i])
                bs = stf.extract(data,fe['BrightScram']['State'])
                ws = stf.extract(data,fe['WindScram']['State'])
                if bs:
                    line += "    BRIGHTSCRAM    "
                elif ws:
                    line += "     WINDSCRAM     "
                elif tostow:
                    line += "      TO STOW      "
                elif atstow:
                    line += "      AT STOW      "
                elif atservice:
                    line += "     AT SERVICE    "
                elif toservice:
                    line += "     TO SERVICE    "
                else:
                    line += "  "+'{:8.4f}'.format(azel['RequestedAzimuth'][i])
                    line += " "+'{:+8.4f}'.format(azel['RequestedElevation'][i])
                line += "  "+'{:+8.4f}'.format(azel['dAzimuth'][i])
                line += " "+'{:+8.4f}'.format(azel['dElevation'][i])
                azoff, eloff = stf.extract(data,c['AzOffset'])/10000., stf.extract(data,c['ElOffset'])/10000.
                raoff, decoff = stf.extract(data,c['RAOffset'])/10000.,stf.extract(data,c['DecOffset'])/10000.
                if i > 12:
                    line += " "+'{:+5.2f}{:+5.2f}'.format(azoff,eloff)
                    line += " "+'{:+5.2f}{:+5.2f}'.format(raoff,decoff)
                else:
                    line += "  "+'{:+4.1f} {:+4.1f}'.format(azoff,eloff)
                    line += "  "+'{:+4.1f} {:+4.1f}'.format(raoff,decoff)
                ra = RA_Angle(stf.extract(data,c['RAVirtualAxis'])/10000.,'degrees')
                line += "  "+ra.get('hms')[:10]
                dec = Dec_Angle(stf.extract(data,c['DecVirtualAxis'])/10000.,'degrees')
                line += " "+dec.get('dms')[:9]
                self.L3.insert(END,line)
                if abs(azoff) > 0 or abs(eloff) > 0 or abs(raoff) > 0 or abs(decoff) > 0:
                    # Antennas with non-zero offsets (color orange)
                    self.L3.itemconfig(END,bg=self.colors['offsets'])
                if abs(azel['dAzimuth'][i]) > 0.005 or abs(azel['dElevation'][i]) > 0.005:
                    # Antennas with tracking error > 0.005 (color yellow)
                    self.L3.itemconfig(END,bg=self.colors['warn'])
                if abs(azel['dAzimuth'][i]) > 0.02 or abs(azel['dElevation'][i]) > 0.02:
                    # Antennas with tracking error > 0.02 (color red)
                    self.L3.itemconfig(END,bg=self.colors['error'])
                if ws + bs:
                    # Antennas with wind or bright scram (color red)
                    self.L3.itemconfig(END,bg=self.colors['error'])
                if not subarray1[i]:
                    # Antennas out of service (color gray)
                    self.L3.itemconfig(END,bg=self.colors['na'])
                else:
                    # Antennas are in service, so check ra,dec coordinates for zeros
                    # (indicates cRIO not reporting)
                    if ra.get() == 0.0 and dec.get() == 0.0:
                        self.L3.itemconfig(END,bg=self.colors['error'])
                if i == 11:
                    self.L3.insert(END,head2)
                    self.L3.itemconfig(END,bg=self.colors['colhead'])


        # ================= Section 2: Communications ===================
        nsubarray = 0
        nabad = 0
        ncbad = 0
        for i in antindex:
            nsubarray += 1
            c = sf['Antenna'][i]['Controller']
            criomjd = int(mjd) + stf.extract(data,c['cRIOClockms'])/86400000.
            #criotimestr = Time(criomjd,format='mjd').iso
            # Controller MJD occasionally glitches to very large values, so limit it to no larger than today's MJD
            antmjd = min(stf.extract(data,c['SystemClockMJDay']),int(mjd)) + stf.extract(data,c['SystemClockms'])/86400000.
            #anttimestr = Time(antmjd,format='mjd').iso
            criotdif = (criomjd - mjd)*86400.
            anttdif = (antmjd - mjd)*86400.
            if abs(criotdif) > 1.0: ncbad +=1
            if abs(anttdif) > 1.0: nabad += 1

        heading = 'Communications'
        if not self.sectionDisplayState[1]:
            line = 'cRIO Good:{:2}  Bad:{:2}  Ant Good:{:2}  Bad:{:2}'.format(nsubarray-ncbad,ncbad,nsubarray-nabad,nabad).center(100-len(heading)-len(expand),'_')
            headline = '_'+heading+line+expand
            self.L3.insert(END,headline)
            self.L3.itemconfig(END,bg=self.colors['section1'],fg=self.colors['na'])
        else:
            headline = ' '+heading+' '*(100-len(heading)-len(collapse))+collapse
            self.L3.insert(END,headline)
            self.L3.itemconfig(END,bg=self.colors['section1'],fg='white')

            head = 'a# Pwr  Runctrl  Runmode  Datamode   cRIO_Time     Ant_Time     cRIO--TDIF--Ant'
            self.L3.insert(END,head)
            self.L3.itemconfig(END,bg=self.colors['colhead'])

            for i in antindex:
                c = sf['Antenna'][i]['Controller']
                criomjd = int(mjd) + stf.extract(data,c['cRIOClockms'])/86400000.
                if criomjd == 0.0:
                    criotimestr = '1904-01-01 00:00:00.000'
                else:
                    criotimestr = Time(criomjd,format='mjd').iso
                # Controller MJD occasionally glitches to very large values, so limit it to no larger than today's MJD
                antmjd = min(stf.extract(data,c['SystemClockMJDay']),int(mjd)) + stf.extract(data,c['SystemClockms'])/86400000.
                if antmjd == 0.0:
                    anttimestr = '1904-01-01 00:00:00.000'
                else:
                    anttimestr = Time(antmjd,format='mjd').iso
                criotdif = (criomjd - mjd)*86400.
                anttdif = (antmjd - mjd)*86400.
                line = antn[i]+pswitch[np.clip(stf.extract(data,c['PowerSwitch']),0,2)]
                line += "  "+rctrl[np.clip(stf.extract(data,c['RunControl']),0,2)]
                line += "  "+rmode[np.clip(stf.extract(data,c['RunMode']),0,5)]
                line += "  "+dmode[np.clip(stf.extract(data,c['DataMode']),0,2)]
                line += "  "+criotimestr[11:]+"  "+anttimestr[11:]+'{:8.3f}  {:8.3f}'.format(criotdif,anttdif)
                self.L3.insert(END,line)
                if abs(criotdif) > 1.0 or abs(anttdif) > 1.0 or line.find('-ERROR--') != -1:
                    self.L3.itemconfig(END,bg=self.colors['error'])
                if not subarray1[i]:
                    # Antennas out of service (color gray)
                    self.L3.itemconfig(END,bg=self.colors['na'])

        # ================= Section 3: Frequency Tuning ===================
        errstr = stf.extract(data,sf['LODM']['LO1A']['ERR']).split('\n')[0]
        errvals = errstr.split(',')
        if len(errvals) == 2:
            errno,errmsg = errvals
        else:
            errno = 0
            errmsg = errvals[0]
        stat = stf.extract(data,sf['LODM']['LO1A']['SweepStatus'])
        statdict = {0:'Stopped',8:'Sweeping',32:'Wait4Trg'}
        fseqfile = stf.extract(data,sf['LODM']['LO1A']['FSeqFile']).strip('\x00')
        try:
            status = statdict[stat]
        except KeyError:
            status = 'Unknown'

        heading = 'Frequency Tuning'
        if not self.sectionDisplayState[2]:
            line = '{:.8} <{:.19}> {:.30}'.format(status,fseqfile,errmsg).center(100-len(heading)-len(expand),'_')
            headline = '_'+heading+line+expand
            self.L3.insert(END,headline)
            self.L3.itemconfig(END,bg=self.colors['section0'],fg=self.colors['na'])
        else:
            headline = ' '+heading+' '*(100-len(heading)-len(collapse))+collapse
            self.L3.insert(END,headline)
            self.L3.itemconfig(END,bg=self.colors['section0'],fg='white')

            head = 'LO1A Sweep Status  FSeqFile             ESR STB  Err   ErrorMsg'
            self.L3.insert(END,head)
            self.L3.itemconfig(END,bg=self.colors['colhead'])

            esr = str(stf.extract(data,sf['LODM']['LO1A']['ESR']))
            stb = str(stf.extract(data,sf['LODM']['LO1A']['STB']))
            line = '    {:9.8}      {:20.19}  {:3} {:3} {:5} {:5}'.format(status,fseqfile,esr,stb,errno,errmsg)
            self.L3.insert(END,line)
            if errmsg.find('No error') != 1 :
                self.L3.itemconfig(END,bg=self.colors['error'])

        # ================= Section 4: Phase Tracking ===================
        c = sf['Schedule']['Data']
        scdtime = stf.extract(data,c['Timestamp'])
        if scdtime > 8640000000.:
            scdline = '**Timestamp contains garbage**'
        else:
            scdline = Time(scdtime+0.000001,format='lv').iso
        scan_state = stf.extract(data,c['ScanState'])
        if scan_state == -1:
            stateline = 'OFF'
        elif scan_state == 0:
            stateline = 'UNK'
        else:
            stateline = 'ON '
        if stf.extract(data,c['PhaseTracking']) == 0:
            trackline = 'False'
        else:
            trackline = 'True '
        heading = 'Phase Tracking'
        if not self.sectionDisplayState[3]:
            line = (scdline+'  State: '+stateline+'  Tracking: '+trackline).center(100-len(heading)-len(expand),'_')
            headline = '_'+heading+line+expand
            self.L3.insert(END,headline)
            self.L3.itemconfig(END,bg=self.colors['section1'],fg=self.colors['na'])
        else:
            headline = ' '+heading+' '*(100-len(heading)-len(collapse))+collapse
            self.L3.insert(END,headline)
            self.L3.itemconfig(END,bg=self.colors['section1'],fg='white')

            line = 'UVW Time: '+scdline
            line += '  Scan State: '+stateline
            line += '  Phase Tracking: '+trackline
            self.L3.insert(END,line)
            if scan_state == 0:
                self.L3.itemconfig(END,bg=self.colors['error'])
     
            head = 'a# --U[ns]-- --V[ns]-- --W[ns]-- DLA[ns]->----t---- ---t+1---'
            self.L3.insert(END,head)
            self.L3.itemconfig(END,bg=self.colors['colhead'])
            uvw = stf.extract(data,c['UVW'])
            delay = stf.extract(data,c['Delay'])
            for i in antindex:
                line = '{0:2d} {1:9.3f} {2:9.3f} {3:9.3f}          {4:9.3f} {5:9.3f}'.format(i+1,uvw[i,0],uvw[i,1],uvw[i,2],delay[0,i],delay[1,i])
                self.L3.insert(END,line)
                if np.all(uvw[1,0:3] == [0.0,0.0,0.0]):
                    # Gray out all lines if ant2 (or any ant other than 1) is all zeros 
                    self.L3.itemconfig(END,bg=self.colors['na'])

        # ================= Section 5: Power and Attenuation ===================
        nsubarray = 0
        nbbad = 0
        nfbad = 0
        for i in antindex:
            if subarray1[i]:
                nsubarray += 1
                fe = sf['Antenna'][i]['Frontend']['FEM']
                hpv = stf.extract(data,fe['HPol']['Power'])
                vpv = stf.extract(data,fe['VPol']['Power'])
                if hpv == 0.0 or vpv == 0.0: nfbad += 1
                be = sf['DCM'][i]
                hpv = stf.extract(data,be['HPol']['Voltage'])
                vpv = stf.extract(data,be['VPol']['Voltage'])
                if hpv == 0.0 or vpv == 0.0: nbbad += 1

        heading = 'Power and Attenuation'
        if not self.sectionDisplayState[4]:
            line = 'FE Good:{:2}  Bad:{:2}  BE Good:{:2}  Bad:{:2}'.format(nsubarray-nfbad,nfbad,nsubarray-nbbad,nbbad).center(100-len(heading)-len(expand),'_')
            headline = '_'+heading+line+expand
            self.L3.insert(END,headline)
            self.L3.itemconfig(END,bg=self.colors['section0'],fg=self.colors['na'])
        else:
            headline = ' '+heading+' '*(100-len(heading)-len(collapse))+collapse
            self.L3.insert(END,headline)
            self.L3.itemconfig(END,bg=self.colors['section0'],fg='white')

            head = '     Frontend:  --H-Channel---  --V-Channel---   Backend:  --H-Channel---  --V-Channel---'
            self.L3.insert(END,head)
            self.L3.itemconfig(END,bg=self.colors['colhead'])
            head = 'a#              --dBm-- -Attn-  --dBm-- -Attn-             -Volts- -Attn-  -Volts- -Attn-'
            self.L3.insert(END,head)
            self.L3.itemconfig(END,bg=self.colors['colhead'])

            for i in antindex:
                fe = sf['Antenna'][i]['Frontend']['FEM']
                hpv = stf.extract(data,fe['HPol']['Power'])
                vpv = stf.extract(data,fe['VPol']['Power'])
                ndon = ' '
                if stf.extract(data,fe['ND']) == 1:
                    ndon = '**ND-ON**'
                be = sf['DCM'][i]
                line = '{:2d}  {:>9}  {:7.4f} {:3d} {:3d} {:7.4f} {:3d} {:3d}            {:7.4f} {:5d}   {:7.4f} {:5d}'.format(i+1,ndon,
                    hpv,stf.extract(data,fe['HPol']['Attenuation']['First']),
                    stf.extract(data,fe['HPol']['Attenuation']['Second']),
                    vpv,stf.extract(data,fe['VPol']['Attenuation']['First']),
                    stf.extract(data,fe['VPol']['Attenuation']['Second']),
                    stf.extract(data,be['HPol']['Voltage']),stf.extract(data,be['HPol']['Attenuation']),
                    stf.extract(data,be['VPol']['Voltage']),stf.extract(data,be['VPol']['Attenuation']))
                self.L3.insert(END,line)
                if hpv == 0.0 or vpv == 0.0:
                    self.L3.itemconfig(END,bg=self.colors['error'])
                if ndon != ' ':
                    self.L3.itemconfig(END,bg=self.colors['warn'])
                if not subarray1[i]:
                    # Antennas out of service (color gray)
                    self.L3.itemconfig(END,bg=self.colors['na'])

        # ================= Section 6: Front End Temperature Controller ===================
        nsubarray = 0
        nlbad = 0
        nfbad = 0
        for i in antindex:
            if subarray1[i]:
                nsubarray += 1
                laird = sf['Antenna'][i]['Frontend']['TEC']
                lairdtemp = stf.extract(data,laird['Temperature'])
                if  lairdtemp < 24 or lairdtemp > 26: nlbad += 1
                fe = sf['Antenna'][i]['Frontend']['FEM']
                fetemp = stf.extract(data,fe['Temperature'])
                if fetemp < 20 or fetemp > 30: nfbad += 1

        heading = 'Front End Temperatures'
        if not self.sectionDisplayState[5]:
            line = 'FE Good:{:2}  Bad:{:2}  Laird Good:{:2}  Bad:{:2}'.format(nsubarray-nfbad,nfbad,nsubarray-nlbad,nlbad).center(100-len(heading)-len(expand),'_')
            headline = '_'+heading+line+expand
            self.L3.insert(END,headline)
            self.L3.itemconfig(END,bg=self.colors['section1'],fg=self.colors['na'])
        else:
            headline = ' '+heading+' '*(100-len(heading)-len(collapse))+collapse
            self.L3.insert(END,headline)
            self.L3.itemconfig(END,bg=self.colors['section1'],fg='white')

            head = 'a#  -FETemp (C)-  -Temp (C)- -TSet (C)-  --Dty%-   -Volts-  -Curr- -Alarm- -Error- -EHist-'
            self.L3.insert(END,head)
            self.L3.itemconfig(END,bg=self.colors['colhead'])
            for i in antindex:
                fe = sf['Antenna'][i]['Frontend']['FEM']
                fetemp = stf.extract(data,fe['Temperature'])
                laird = sf['Antenna'][i]['Frontend']['TEC']
                lairdtemp = stf.extract(data,laird['Temperature'])
                line = '{:2d}     {:6.3f}       {:6.3f}     {:6.3f}   {:8.3f}  {:7.3f}  {:6.3f}   {:3d}     {:3d}     {:3d}'.format(i+1,
                    fetemp,lairdtemp,stf.extract(data,laird['TReference']),
                    stf.extract(data,laird['DutyFactor']),stf.extract(data,laird['InputVoltage']),
                    stf.extract(data,laird['MainCurrent']),stf.extract(data,laird['Alarm']),
                    stf.extract(data,laird['Error']),stf.extract(data,laird['ErrorHistory']))
                if i == 13:
                    lairdtemp = stf.extract(data,sf['FEMA']['Thermal']['SecondStageTemp'])
                    line = '{:2d}     {:6.3f}       {:6.3f} K     --         --       --       --     ---     ---     ---'.format(i+1,fetemp,lairdtemp)
                self.L3.insert(END,line)
                if i == 13 and (lairdtemp < 10 or lairdtemp > 50):
                    self.L3.itemconfig(END,bg=self.colors['error'])
                if i < 13 and (lairdtemp < 24 or lairdtemp > 26):
                    self.L3.itemconfig(END,bg=self.colors['error'])
                if fetemp < 20 or fetemp > 30:
                    self.L3.itemconfig(END,bg=self.colors['error'])
                if not subarray1[i]:
                    self.L3.itemconfig(END,bg=self.colors['na'])

        # ================= Section 7: ROACH Status ===================
        nroaches = 8
        nrbad = 0
        r = sf['Schedule']['Data']['Roach']
        for i in range(nroaches):
            fpgafan = stf.extract(data,r[i]['Fan.fpga'])
            if fpgafan == 0: nrbad += 1
        
        heading = 'ROACH Status'
        if not self.sectionDisplayState[6]:
            line = 'Connected:{:2}  Not Connected:{:2}'.format(nroaches-nrbad,nrbad).center(100-len(heading)-len(expand),'_')
            headline = '_'+heading+line+expand
            self.L3.insert(END,headline)
            self.L3.itemconfig(END,bg=self.colors['section0'],fg=self.colors['na'])
        else:
            headline = ' '+heading+' '*(100-len(heading)-len(collapse))+collapse
            self.L3.insert(END,headline)
            self.L3.itemconfig(END,bg=self.colors['section0'],fg='white')

            head =  '   --------Voltage (VDC)------- ---------Current (A)------ ----Temperature--- ---------Delays--------'
            head2 = 'r# 1.0 1.5 1.8 2.5 3.3 5.0 12.0 1.0 1.5 1.8 2.5 3.3 5.0 12 AMB FPG PPC IN OUT  X--CHN1--Y  X--CHN2--Y'
            fmt = '{:2} {:3} {:3} {:3} {:3} {:3} {:3} {:4} {:3} {:3} {:3} {:3} {:3} {:3} {:3} {:2}  {:2}  {:2} {:2}  {:2} {:5d} {:5d} {:5d} {:5d}'
            self.L3.insert(END,head)
            self.L3.itemconfig(END,bg=self.colors['colhead'])
            self.L3.insert(END,head2)
            self.L3.itemconfig(END,bg=self.colors['colhead'])
      
            for i in range(nroaches):
                v1 = str(stf.extract(data,r[i]['Voltage.1v'])+0.05)[:3]
                v15 = str(stf.extract(data,r[i]['Voltage.1v5'])+0.05)[:3]
                v18 = str(stf.extract(data,r[i]['Voltage.1v8'])+0.05)[:3]
                v25 = str(stf.extract(data,r[i]['Voltage.2v5'])+0.05)[:3]
                v33 = str(stf.extract(data,r[i]['Voltage.3v3'])+0.05)[:3]
                v5 = str(stf.extract(data,r[i]['Voltage.5v'])+0.05)[:3]
                v12 = str(stf.extract(data,r[i]['Voltage.12v'])*3+0.05)[:4]
                c1 = str(stf.extract(data,r[i]['Current.1v'])+0.05)[:3]
                c15 = str(stf.extract(data,r[i]['Current.1v5'])+0.05)[:3]
                c18 = str(stf.extract(data,r[i]['Current.1v8'])+0.05)[:3]
                c25 = str(stf.extract(data,r[i]['Current.2v5'])+0.05)[:3]
                c33 = str(stf.extract(data,r[i]['Current.3v3'])+0.05)[:3]
                c5 = str(stf.extract(data,r[i]['Current.5v'])+0.05)[:3]
                c12 = str(stf.extract(data,r[i]['Current.12v'])+0.05)[:3]
                tamb = str(stf.extract(data,r[i]['Temp.ambient'])+0.5).split('.')[0]
                tfpga = str(stf.extract(data,r[i]['Temp.fpga'])+0.5).split('.')[0]
                tin = str(stf.extract(data,r[i]['Temp.inlet'])+0.5).split('.')[0]
                tout = str(stf.extract(data,r[i]['Temp.outlet'])+0.5).split('.')[0]
                tppc = str(stf.extract(data,r[i]['Temp.ppc'])+0.5).split('.')[0]
                dla0x = stf.extract(data,r[i]['Delay0x'])
                dla0y = stf.extract(data,r[i]['Delay0y'])
                dla1x = stf.extract(data,r[i]['Delay1x'])
                dla1y = stf.extract(data,r[i]['Delay1y'])
                line = fmt.format(i+1,v1,v15,v18,v25,v33,v5,v12,c1,c15,c18,c25,c33,c5,c12,
                                  tamb,tfpga,tppc,tin,tout,dla0x,dla0y,dla1x,dla1y)
                self.L3.insert(END,line)
                fpgafan = stf.extract(data,r[i]['Fan.fpga'])
                if fpgafan == 0:
                    self.L3.itemconfig(END,bg=self.colors['na'])
                status = hex(stf.extract(data,r[i]['Status']))
                if status != '0xc00000' and status != '0x0':
                    # For some reason, a good ROACH reports '0xc00000' for status
                    self.L3.itemconfig(END,bg=self.colors['error'])

        # ================= Section 8: Last Antenna Command ===================
        heading = 'Antenna Last Command'
        
        if not self.sectionDisplayState[7]:
            headline = ' '+heading+' '*(100-len(heading)-len(expand))+expand
            self.L3.insert(END,headline)
            self.L3.itemconfig(END,bg=self.colors['section0'],fg=self.colors['na'])
        else:
            headline = ' '+heading+' '*(100-len(heading)-len(collapse))+collapse
            self.L3.insert(END,headline)
            self.L3.itemconfig(END,bg=self.colors['section0'],fg='white')
            for i in antindex:
                parser = sf['Antenna'][i]['Parser']
                an='  '+'{:2d}'.format(i)+":   "
                self.L3.insert(END,an+stf.extract(data,parser['Command']).strip('\x00')+' '+str(stf.extract(data,parser['CommErr'])))

    def cryo_display(self,data):
        ''' Creates the CryoRX page to display Cryo Receiver information
        '''
        onoff = ['  OFF  ','  ON   ']
        self.cryoLB.delete(0,END)
        sf = self.accini['sf']
        cryos = [sf['FEMA']] #, sf['FEMB']]
        cryostr = ['FEMA','FEMB']
        for k,cryo in enumerate(cryos):
            ps = cryo['PowerStrip']
            lvtstamp = stf.extract(data,cryo['Timestamp'])
            if int(lvtstamp) == 0:
                line = 'Time: 1904-01-01 00:00:00.000   Version: '+'{:4.1f}'.format(stf.extract(data,cryo['Version']))
            else:
                line = 'Time: '+Time(lvtstamp,format='lv').iso+'   Version: '+'{:4.1f}'.format(stf.extract(data,cryo['Version']))
            self.cryoLB.insert(END,line)
            # Section on outlet status
            line = cryostr[k]+' Outlets               Status'
            self.cryoLB.insert(END,line)
            line = 'Voltage [V]   Current [A]  Computer  Brick   LNA5VBB  LNA12V   2ndAmp   OptRabit RFSwitch NoiseDiode'
            self.cryoLB.insert(END,line)
            self.cryoLB.itemconfig(END,bg=self.colors['colhead'])
            line = '--1-----2--  --1-----2--   ---4----  --3---  ---6---  ---5---  ---7---  ---2---- ---1---- ----8-----'
            self.cryoLB.insert(END,line)
            self.cryoLB.itemconfig(END,bg=self.colors['colhead'])
            cur = stf.extract(data,ps['Current'])
            vol = stf.extract(data,ps['Volts'])
            line = '{:>5.0f} {:>5.0f}  {:>5.1f} {:>5.1f}  '.format(vol[0],vol[1],cur[0],cur[1])
            line += '  '+onoff[stf.extract(data,ps['ComputerStatus'])]+'  '+onoff[stf.extract(data,ps['DeltaTauBrickStatus'])]
            line += '  '+onoff[stf.extract(data,ps['LNA5VBiasBBStatus'])]+'  '+onoff[stf.extract(data,ps['LNA12VBiasStatus'])]
            line += '  '+onoff[stf.extract(data,ps['SecondAmpsStatus'])]+'  '+onoff[stf.extract(data,ps['OpticalTxRabbitStatus'])]
            line += '  '+onoff[stf.extract(data,ps['RFSwitchStatus'])]+'  '+onoff[stf.extract(data,ps['NoiseDiodeStatus'])]
            self.cryoLB.insert(END,line)
            # Section on Receiver
            rx = cryo['Receiver']
            hifreq = 'Disabled'
            lofreq = 'Disabled'
            nd = 'OFF'
            if stf.extract(data,rx['HiFreqEnabled']):
                hifreq = 'Enabled '
            if stf.extract(data,rx['LoFreqEnabled']):
                lofreq = 'Enabled '
            if stf.extract(data,rx['NoiseDiodeEnabled']):
                nd = 'ON '
            line = cryostr[k]+' Receiver     HighFrequency: '+hifreq+'  LowFrequency: '+lofreq+'  NoiseDiode: '+nd
            self.cryoLB.insert(END,line)
            line = 'LNA  DrainVoltage  DrainCurrent  GateAVoltage  GateACurrent  GateBVoltage  GateBCurrent'
            self.cryoLB.insert(END,line)
            self.cryoLB.itemconfig(END,bg=self.colors['colhead'])
            line = '----  ---Volts----  -----mA-----  ---Volts----  -----mA-----  ---Volts----  -----mA-----'
            self.cryoLB.insert(END,line)
            self.cryoLB.itemconfig(END,bg=self.colors['colhead'])
            lnas = {'Hi H':0,'Lo H':1,'Lo V':2,'Hi V':3}
            lna_order = ['Lo H','Lo V','Hi H','Hi V']
            for x in lna_order:
                lna = rx['LNAs'][lnas[x]]
                line = '{:4s}   {:8.2f}      {:8.2f}      {:8.2f}      {:8.2f}      {:8.2f}      {:8.2f}'.format(x,
                    stf.extract(data,lna['DrainVoltage']),stf.extract(data,lna['DrainCurrent']),
                    stf.extract(data,lna['GateAVoltage']),stf.extract(data,lna['GateACurrent']),
                    stf.extract(data,lna['GateBVoltage']),stf.extract(data,lna['GateBCurrent']))
                self.cryoLB.insert(END,line)
            # Section on Temperatures
            rt = cryo['Thermal']
            fe14temp = stf.extract(data,sf['Antenna'][13]['Frontend']['FEM']['Temperature'])

            line = 'Temps:  FocusBox  RadShield  2ndStage  1stStage  HiFrqPlate  LoFrqHorn  HiFrqHorn  LoFrqLNA  HiFrqLNA'
            self.cryoLB.insert(END,line)
            self.cryoLB.itemconfig(END,bg=self.colors['colhead'])
            line = '        ---C----  ----K----  ---K----  ---K----  ----K-----  ----K----  ----K----  ---K----  ---K----'
            self.cryoLB.insert(END,line)
            self.cryoLB.itemconfig(END,bg=self.colors['colhead'])
            line = '        '
            line += '  {:5.1f}  '.format(fe14temp)
            line += '   {:5.1f}  '.format(stf.extract(data,rt['RadiationShieldTemp']))
            line += '   {:5.1f}  '.format(stf.extract(data,rt['SecondStageTemp']))
            line += '   {:5.1f}   '.format(stf.extract(data,rt['FirstStageTemp']))
            line += '   {:5.1f}    '.format(stf.extract(data,rt['HiFreq15KPlateTemp']))
            line += '   {:5.1f}   '.format(stf.extract(data,rt['LowFreqFeedhornTemp']))
            line += '   {:5.1f}   '.format(stf.extract(data,rt['HiFreqFeedhornTemp']))
            line += '   {:5.1f}  '.format(stf.extract(data,rt['LowFreqLNATemp']))
            line += '   {:5.1f}  '.format(stf.extract(data,rt['HiFreqLNATemp']))
            self.cryoLB.insert(END,line)
            # Section on Motors
            serv = cryo['FRMServo']
            selrx = '***Unset***'
            if stf.extract(data,serv['SelectedRx']) == 1:
                selrx = 'Low Freq RX '
            elif stf.extract(data,serv['SelectedRx']) == 2:
                selrx = 'High Freq RX'
            homed = '       '
            if stf.extract(data,serv['Homed']) == 1: homed = '*HOMED*'
            line = ' Motor      Fault    Limit  Current  Position  PosErr  PosOff  Home Indicator:   Selected RX: '
            self.cryoLB.insert(END,line)
            self.cryoLB.itemconfig(END,bg=self.colors['colhead'])
            line = '--------  ---------  -----  -------  --------  ------  ------      '+homed+'      '+selrx
            self.cryoLB.insert(END,line)
            self.cryoLB.itemconfig(END,bg=self.colors['colhead'])
            motors = [serv['PositionAngle'],serv['RxSelect'],serv['ZFocus']]
            motstr = ['Pangle     ','RXSelect   ','ZFocus     ']
            for j,motor in enumerate(motors):
                fault = 'No Fault '
                if stf.extract(data,motor['AmplifierFault']):
                    fault = 'Amp Fault'
                limit = 'None '
                if stf.extract(data,motor['MinusLimit']):
                    limit = 'LoLim'
                if stf.extract(data,motor['PlusLimit']):
                    limit = 'HiLim'
                line = motstr[j]+fault+'  '+limit+'  {:5.1f}    {:5.1f}   {:5.1f}   {:5.1f}'.format(
                  stf.extract(data,motor['MotorCurrent']),stf.extract(data,motor['Position']),
                  stf.extract(data,motor['PositionError']),stf.extract(data,motor['PositionOffset']))
                self.cryoLB.insert(END,line)
            line = ''
            self.cryoLB.insert(END,line)
            self.cryoLB.insert(END,line)

        pswitch = ['OFF','ON ','-ERROR--']
        rmode = ['STOP    ','POSITION','VELOCITY','---NA---','TRACK   ','-ERROR--']
        rctrl = ['STANDBY','OPERATE','-ERROR--']
        dmode = ['AZ-EL ','RA-DEC','-ERROR--']

        antn = [' A ']#,' B ']
        altantindex = np.array([13])#,14])
        # Expected STOW coordinates for the different types of antenna.  If/when Ant 12 is changed to
        # the old mount, its values must be changed.

        azel = stf.azel_from_stateframe(sf,data)
        antindex = altantindex

        sub1 = stf.extract(data,sf['LODM']['Subarray1'])
        sub2 = stf.extract(data,sf['LODM']['Subarray2'])
        subarray1 = []
        subarray2 = []
        for i in range(16):
            subarray1.append(sub1 & (1<<i) > 0)
            subarray2.append(sub1 & (1<<i) > 0)

        self.cryoLB.insert(END,'====== Ant 14 Tracking ======')
        head  = 'a# ---HA--Act--Dec--  ---HA--Req--Dec--  ---HA--Err--Dec--  Az-Off-El  RA-Off-Dc  ---RA---Req---Dec---'
        templ = '   xxx.xxxx sxx.xxxx  xxx.xxxx sxx.xxxx  sxx.xxxx sxx.xxxx  xx.x xx.x  xx.x xx.x  xx xx xx.x sxx xx xx'
        self.cryoLB.insert(END,head)
        self.cryoLB.itemconfig(END,bg=self.colors['colhead'])
        for k,i in enumerate(altantindex):
            c = sf['Antenna'][i]['Controller']
            fe = sf['Antenna'][i]['Frontend']
            # If bit 9 of central status is set, dish is going to STOW
            tostow = stf.extract(data,c['CentralStatus']) & 256
            # If dish is in POSITION mode and the requested 
            # AZ, EL are 225, 15, dish is going to SERVICE
            toservice = (int(azel['RequestedAzimuth'][i] + 0.5) == 225 and 
                         int(azel['RequestedElevation'][i] + 0.5) == 15 and
                         rmode[stf.extract(data,c['RunMode'])] == 'POSITION')
            # If dish is at AZ, EL 180, 88, dish is AT STOW 
            atstow = (int(azel['ActualAzimuth'][i] + 0.5) == 0 and 
                         int(azel['ActualElevation'][i] + 0.5) == 29)
            # If dish is at AZ, EL 225, 15, dish is AT SERVICE
            atservice = (int(azel['ActualAzimuth'][i] + 0.5) == 225 and 
                         int(azel['ActualElevation'][i] + 0.5) == 15)
            line = antn[k]+'{:8.4f}'.format(azel['ActualAzimuth'][i])
            line += " "+'{:+8.4f}'.format(azel['ActualElevation'][i])
            bs = stf.extract(data,fe['BrightScram']['State'])
            ws = stf.extract(data,fe['WindScram']['State'])
            if bs:
                line += "    BRIGHTSCRAM    "
            elif ws:
                line += "     WINDSCRAM     "
            elif tostow:
                line += "      TO STOW      "
            elif atstow:
                line += "      AT STOW      "
            elif atservice:
                line += "     AT SERVICE    "
            elif toservice:
                line += "     TO SERVICE    "
            else:
                line += "  "+'{:8.4f}'.format(azel['RequestedAzimuth'][i])
                line += " "+'{:+8.4f}'.format(azel['RequestedElevation'][i])
            line += "  "+'{:+8.4f}'.format(azel['dAzimuth'][i])
            line += " "+'{:+8.4f}'.format(azel['dElevation'][i])
            azoff, eloff = stf.extract(data,c['AzOffset'])/10000., stf.extract(data,c['ElOffset'])/10000.
            raoff, decoff = stf.extract(data,c['RAOffset'])/10000.,stf.extract(data,c['DecOffset'])/10000.
            line += " "+'{:+5.2f}{:+5.2f}'.format(azoff,eloff)
            line += " "+'{:+5.2f}{:+5.2f}'.format(raoff,decoff)
            ra = RA_Angle(stf.extract(data,c['RAVirtualAxis'])/10000.,'degrees')
            line += "  "+ra.get('hms')[:10]
            dec = Dec_Angle(stf.extract(data,c['DecVirtualAxis'])/10000.,'degrees')
            line += " "+dec.get('dms')[:9]
            self.cryoLB.insert(END,line)
            if abs(azoff) > 0 or abs(eloff) > 0 or abs(raoff) > 0 or abs(decoff) > 0:
                # Antennas with non-zero offsets (color orange)
                self.cryoLB.itemconfig(END,bg=self.colors['offsets'])
            if abs(azel['dAzimuth'][i]) > 0.005 or abs(azel['dElevation'][i]) > 0.005:
                # Antennas with tracking error > 0.005 (color yellow)
                self.cryoLB.itemconfig(END,bg=self.colors['warn'])
            if abs(azel['dAzimuth'][i]) > 0.02 or abs(azel['dElevation'][i]) > 0.02:
                # Antennas with tracking error > 0.02 (color red)
                self.cryoLB.itemconfig(END,bg=self.colors['error'])
            if ws + bs:
                # Antennas with wind or bright scram (color red)
                self.cryoLB.itemconfig(END,bg=self.colors['error'])
            if not subarray1[i]:
                # Antennas out of service (color gray)
                self.cryoLB.itemconfig(END,bg=self.colors['na'])
            else:
                # Antennas are in service, so check ra,dec coordinates for zeros
                # (indicates cRIO not reporting)
                if ra.get() == 0.0 and dec.get() == 0.0:
                    self.cryoLB.itemconfig(END,bg=self.colors['error'])

        self.cryoLB.insert(END,' ') # Blank line
        self.cryoLB.insert(END,'====== Ant 14 Communication ======')
        head = 'a# Pwr  Runctrl  Runmode  Datamode   cRIO_Time     Ant_Time     cRIO--TDIF--Ant'
        self.cryoLB.insert(END,head)
        self.cryoLB.itemconfig(END,bg=self.colors['colhead'])
        t = Time(stf.extract(data,sf['Timestamp'])+0.000001,format='lv')
        mjd = t.mjd

        for k,i in enumerate(antindex):
            c = sf['Antenna'][i]['Controller']
            criomjd = int(mjd) + stf.extract(data,c['cRIOClockms'])/86400000.
            if criomjd == 0.0:
                criotimestr = '1904-01-01 00:00:00.000'
            else:
                criotimestr = Time(criomjd,format='mjd').iso
            # Controller MJD occasionally glitches to very large values, so limit it to no larger than today's MJD
            antmjd = min(stf.extract(data,c['SystemClockMJDay']),int(mjd)) + stf.extract(data,c['SystemClockms'])/86400000.
            if antmjd == 0.0:
                anttimestr = '1904-01-01 00:00:00.000'
            else:
                anttimestr = Time(antmjd,format='mjd').iso
            criotdif = (criomjd - mjd)*86400.
            anttdif = (antmjd - mjd)*86400.
            line = antn[k]+pswitch[np.clip(stf.extract(data,c['PowerSwitch']),0,2)]
            line += "  "+rctrl[np.clip(stf.extract(data,c['RunControl']),0,2)]
            line += "  "+rmode[np.clip(stf.extract(data,c['RunMode']),0,5)]
            line += "  "+dmode[np.clip(stf.extract(data,c['DataMode']),0,2)]
            line += "  "+criotimestr[11:]+"  "+anttimestr[11:]+'{:8.3f}  {:8.3f}'.format(criotdif,anttdif)
            self.cryoLB.insert(END,line)
            if abs(criotdif) > 1.0 or abs(anttdif) > 1.0 or line.find('-ERROR--') != -1:
                self.cryoLB.itemconfig(END,bg=self.colors['error'])
            if not subarray1[i]:
                # Antennas out of service (color gray)
                self.cryoLB.itemconfig(END,bg=self.colors['na'])

        self.cryoLB.insert(END,' ') # Blank line
        self.cryoLB.insert(END,'====== Ant 14 Power and Attenuation ======')
        head = '     Frontend:  --H-Channel---  --V-Channel---   Backend:  --H-Channel---  --V-Channel---'
        self.cryoLB.insert(END,head)
        self.cryoLB.itemconfig(END,bg=self.colors['colhead'])
        head = 'a#              --dBm-- -Attn-  --dBm-- -Attn-             -Volts- -Attn-  -Volts- -Attn-'
        self.cryoLB.insert(END,head)
        self.cryoLB.itemconfig(END,bg=self.colors['colhead'])
        for k,i in enumerate(antindex):
            fe = sf['Antenna'][i]['Frontend']['FEM']
            hpv = stf.extract(data,fe['HPol']['Power'])
            vpv = stf.extract(data,fe['VPol']['Power'])
            ndon = ' '
            if stf.extract(data,fe['ND']) == 1:
                ndon = '**ND-ON**'
            be = sf['DCM'][i]
            line = '{:2d}  {:>9}  {:7.4f} {:3d} {:3d} {:7.4f} {:3d} {:3d}            {:7.4f} {:5d}   {:7.4f} {:5d}'.format(i+1,ndon,
                hpv,stf.extract(data,fe['HPol']['Attenuation']['First']),
                stf.extract(data,fe['HPol']['Attenuation']['Second']),
                vpv,stf.extract(data,fe['VPol']['Attenuation']['First']),
                stf.extract(data,fe['VPol']['Attenuation']['Second']),
                stf.extract(data,be['HPol']['Voltage']),stf.extract(data,be['HPol']['Attenuation']),
                stf.extract(data,be['VPol']['Voltage']),stf.extract(data,be['VPol']['Attenuation']))
            self.cryoLB.insert(END,line)
            if hpv == 0.0 or vpv == 0.0:
                self.cryoLB.itemconfig(END,bg=self.colors['error'])
            if ndon != ' ':
                self.cryoLB.itemconfig(END,bg=self.colors['warn'])
            if not subarray1[i]:
                # Antennas out of service (color gray)
                self.cryoLB.itemconfig(END,bg=self.colors['na'])
                
    def plottemp(self):
        ''' Just a simple routine to see if I can plot something once a second...
        '''
        sf = self.accini['sf']
        laird = []
        fe = []
        for i in range(13):
            laird.append(sf['Antenna'][i]['Frontend']['TEC'])
            fe.append(sf['Antenna'][i]['Frontend']['FEM'])
        fe.append(sf['Antenna'][13]['Frontend']['FEM'])  # Ant 14 ambient temperature
        laird.append(sf['FEMA']['Thermal']['SecondStageTemp'])   # Ant 14 cyro temperature
        windkey = sf['Schedule']['Data']['Weather']['AvgWind']
        npts = len(self.que)
#        daydif = (datetime.datetime(1901,1,1)-datetime.datetime(1,1,1)).days
        temps = np.zeros((npts,14,2),dtype=np.float)
        tm = np.zeros(npts,dtype=np.float)   # Timestamp
        amb = np.zeros(npts,dtype=np.float)
        wind = np.zeros(npts,dtype=np.float)
        for i in range(npts):
            amb[i] = stf.extract(self.que[i],sf['Schedule']['Data']['Weather']['Temperature'])
            tm[i] = stf.extract(self.que[i],sf['Timestamp'])
            for j in range(13):
                temps[i,j,0] = stf.extract(self.que[i],laird[j]['Temperature'])
                temps[i,j,1] = stf.extract(self.que[i],fe[j]['Temperature'])
            temps[i,13,0] = stf.extract(self.que[i],laird[13])
            temps[i,13,1] = stf.extract(self.que[i],fe[13]['Temperature'])
            wind[i] = stf.extract(self.que[i],windkey)
        # Exactly zero temperatures mean missing data (except ambient on rare cold days!), so flag with NaN
        amb[np.where(amb == 0.0)[0]] = np.NaN
        amb = (amb-32)*5./9   # Convert ambient to Celcius temperature
        temps[np.where(temps == 0.0)] = np.NaN
        wind[np.where(wind == 0.0)[0]] = np.NaN
        # Get current axis extent
        extent = self.sub_plot1.axis()
        auto = self.sub_plot1.get_autoscale_on()
        # Clear the existing "axes" to start anew with the plots
        self.sub_plot1.cla()
        # Add all of the temperatures to the plot
        tim = Time(tm,format='lv').plot_date
        for i in range(13):
            self.sub_plot1.plot_date(tim,temps[:,i,0],label='Laird A'+str(i+1),marker=None,linestyle='-')
            self.sub_plot1.plot_date(tim,temps[:,i,1],label='FEM A'+str(i+1),marker='.',markersize=0.5)
        # Plot Ant 14 temperatures in RED to highlight them somewhat
        self.sub_plot1.plot_date(tim,temps[:,13,0],label='A14 Cryo [K]',marker=None,color='red',linestyle='-')
        self.sub_plot1.plot_date(tim,temps[:,13,1],label='A14 FE Box'+str(i+1),marker='.',color='red',markersize=0.5)
        self.sub_plot1.plot_date(tim,amb,label='Ambient',marker=None,color='magenta',linestyle='-',linewidth=2)
        self.sub_plot1.plot_date(tim,wind,label='AvgWind',marker=None,color='gray',linestyle='-',linewidth=2)
        # Set appropriate titles, legend, etc.
        self.sub_plot1.set_title('Laird Controller and FEM Temp')
        self.sub_plot1.set_ylabel('Temperature [C]')
        self.sub_plot1.set_xlabel('Time [s]')
        self.sub_plot1.legend(loc='lower left',fontsize='x-small')
        self.sub_plot1.grid()
        # Set up date formatting
        fmt = plt.DateFormatter('%H:%M:%S')
        self.sub_plot1.xaxis.set_major_formatter(fmt)
        if not auto:
           self.sub_plot1.axis(extent)
        # Make it visible
        self.canvas.draw()

    #=========================
    def handle_tab(self):
        ''' This tab permits selection of any part of stateframe for plotting, and optionally
            plots any desired selections.  All we have to do here is check if there is anything
            to plot.
        '''
        plots = []
        nplots = len(self.miscplotlocator)
#        daydif = (datetime.datetime(1901,1,1)-datetime.datetime(1,1,1)).days
        if nplots == 0:
            pass
#            self.sub_plot2.cla()
#            self.canvas.draw()
        elif nplots > 0:
            # Get current axis extent
            extent = self.sub_plot2.axis()
            auto = self.sub_plot2.get_autoscale_on()

            # Clear the existing "axes" to start anew with the plots
            self.sub_plot2.cla()
            npts = len(self.que)
            tm = np.zeros(npts,'float')
            for i in range(npts):
                tm[i] = stf.extract(self.que[i],['d',0])
            tim = Time(tm,format='lv').plot_date
            
            colors = ['blue','red','green','orange','cyan','black','magenta','0.75']
            for locator in self.miscplotlocator:
                plots.append(np.zeros(npts,'float'))
                for i in range(npts):
                    item = stf.extract(self.que[i],locator)
                    if type(item) is list:
                        item = item[0]
                    plots[-1][i] = item
                plots[-1][np.where(plots[-1] == 0.0)[0]] = np.NaN
            for i in range(nplots):
                # Add all of the plots
                self.sub_plot2.plot_date(tim,plots[i],label=self.miscplotlabel[i],marker=None,color=colors[i % 8],linestyle='-')
            self.sub_plot2.legend(loc='lower left',fontsize='x-small')
            self.sub_plot2.grid()
            # Set up date formatting
            fmt = plt.DateFormatter('%H:%M:%S')
            self.sub_plot2.xaxis.set_major_formatter(fmt)
            if not auto:
                self.sub_plot2.axis(extent)
            # Make it visible
            self.canvas2.draw()

            
#=========================
    def handle_Extra(self):
        '''This tab is for plotting previously saved graph structures'''

        idx = self.saved_list.curselection()
        name = self.saved_list.get(idx)

        # Get current axis extent
        extent = self.sub_plot3.axis()
        auto = self.sub_plot3.get_autoscale_on()

        self.sub_plot3.cla()
        plots = []
        npts = len(self.que)
        tm = np.zeros(npts,'float')
        for i in range(npts):
            tm[i] = stf.extract(self.que[i],['d',0])
        tim = Time(tm,format='lv').plot_date

        colors = ['blue','red','green','orange','cyan','black','magenta','0.75']
        temp_plots2 = self.saved_dict_labels[name]
        temp_plots = self.saved_dict[name]
        n_plots = len(temp_plots)

        for locator in temp_plots:
            plots.append(np.zeros(npts,'float'))
            for i in range(npts):
                item = stf.extract(self.que[i],locator)
                if type(item) is list:
                    item = item[0]
                plots[-1][i] = item
            plots[-1][np.where(plots[-1] == 0.0)[0]] = np.NaN
        
        for i in range(n_plots):
            self.sub_plot3.plot_date(tim,plots[i],label=temp_plots2[i],marker=None,color=colors[i % 8],linestyle='-')

        self.sub_plot3.legend(loc='lower left',fontsize='x-small')      
        self.sub_plot3.grid()
        # Set up date formatting
        fmt = plt.DateFormatter('%H:%M:%S')
        self.sub_plot3.xaxis.set_major_formatter(fmt)
        if not auto:
           self.sub_plot3.axis(extent)
        # Make it visible
        self.canvas3.draw()

#=========================
    def handle_Extra_plots(self):

        #self.extra_plots[name]=[sub_plot,canvas]
        name = self.cur_tabself.AzStatus
        sub_plot = self.extra_plots[name][0]
        canvas = self.extra_plots[name][1]

        # Get current axis extent
        extent = sub_plot.axis()
        auto = sub_plot.get_autoscale_on()

        sub_plot.cla()
        plots = []
        npts = len(self.que)
        tm = np.zeros(npts,'float')
        for i in range(npts):
            tm[i] = stf.extract(self.que[i],['d',0])
        tim = Time(tm,format='lv').plot_date

        colors = ['blue','red','green','orange','cyan','black','magenta','0.75']
        temp_plots2 = self.saved_dict_labels[name]
        temp_plots = self.saved_dict[name]
        n_plots = len(temp_plots)

        for locator in temp_plots:
            plots.append(np.zeros(npts,'float'))
            for i in range(npts):
                plots[-1][i] = stf.extract(self.que[i],locator)
            plots[-1][np.where(plots[-1] == 0.0)[0]] = np.NaN
        
        for i in range(n_plots):
            sub_plot.plot_date(tim,plots[i],label=temp_plots2[i],marker=None,color=colors[i % 8],linestyle='-')

        sub_plot.legend(loc='lower left',fontsize='x-small')      
        sub_plot.grid()
        # Set up date formatting
        fmt = plt.DateFormatter('%H:%M:%S')
        sub_plot.xaxis.set_major_formatter(fmt)
        if not auto:
           sub_plot.axis(extent)
        # Make it visible
        canvas.draw()
            
        
#=================================
class Status:
    def __init__(self, value, sizebytes):
        self.bits = []
        for i in range(sizebytes*8):
            self.bits.append( (value & (1 << i)) >> i )
        self.on = ['']*sizebytes*8
        self.off = ['']*sizebytes*8
    def set(self, idx, onlist, offlist):
        # Given a list of indexes, ON strings and OFF strings
        # set the values of hte ON and OFF strings
        for i, k in enumerate(idx):
            self.on[k] = onlist[i]
            self.off[k] = offlist[i]
    def get(self, idx):
        out = []
        for i in idx:
            if self.bits[i]:
                out.append((self.on[i],'red'))
            else:
                out.append((self.off[i],'green'))
        return out

class AzStatus(Status):
    def __init__(self, value, sizebytes=3):
        self.bits = []
        for i in range(sizebytes*8):
            self.bits.append( (value & (1 << i)) >> i )
        self.on = ['']*sizebytes*8
        self.off = ['']*sizebytes*8
        onlist = ['Tripped','Inactive','Local','Brake','Lo Soft Lim','Hi Soft Lim','Lo Hard Lim','Hi Hard Lim']
        offlist = ['Healthy','Active','Remote','Brake','Lo Soft Lim','Hi Soft Lim','Lo Hard Lim','Hi Hard Lim']
        idx = range(8)
        self.set(idx,onlist,offlist)
        onlist = ['Position','Speed','Pos Offset','Vel Offset','Axis Lock','Lo DemandLim','Hi DemandLim']
        offlist = ['Position','Speed','Pos Offset','Vel Offset','Axis Lock','Lo Demand Ok','Hi Demand Ok']
        idx = range(8,15)
        self.set(idx,onlist,offlist)
        onlist = ['Drive Ena','Permit','SpdDemandLim','Brake Alarm','Brake Disable']
        offlist = ['Drive Ena','Permit','SpdDemandLim','Brake Okay','Brake Enable']
        idx = [16,17,19,21,22]
        self.set(idx,onlist,offlist)
        
class ElStatus(Status):
    def __init__(self, value, sizebytes=3):
        self.bits = []
        for i in range(sizebytes*8):
            self.bits.append( (value & (1 << i)) >> i )
        self.on = ['']*sizebytes*8
        self.off = ['']*sizebytes*8
        onlist = ['Tripped','Inactive','Local','Brake','Lo Soft Lim','Hi Soft Lim','Lo Hard Lim','Hi Hard Lim']
        offlist = ['Healthy','Active','Remote','Brake','Lo Soft Lim','Hi Soft Lim','Lo Hard Lim','Hi Hard Lim']
        idx = range(8)
        self.set(idx,onlist,offlist)
        onlist = ['Position','Speed','Pos Offset','Vel Offset','Axis Lock','Lo DemandLim','Hi DemandLim']
        offlist = ['Position','Speed','Pos Offset','Vel Offset','Axis Lock','Lo Demand Ok','Hi Demand Ok']
        idx = range(8,15)
        self.set(idx,onlist,offlist)
        onlist = ['Drive Ena','Permit','SpdDemandLim','Brake Alarm','Brake Disable']
        offlist = ['Drive Ena','Permit','SpdDemandLim','Brake Okay','Brake Enable']
        idx = [16,17,19,21,22]
        self.set(idx,onlist,offlist)

class CenStatus(Status):
    def __init__(self, value, sizebytes=4):
        self.bits = []
        for i in range(sizebytes*8):
            self.bits.append( (value & (1 << i)) >> i )
        # Swap order of second and third bits since that matches better with the other status words
        tmp = self.bits[2]
        self.bits[2] = self.bits[1]
        self.bits[1] = tmp

        self.on = ['']*sizebytes*8
        self.off = ['']*sizebytes*8
        onlist = ['PwrConAux','Standby','Local','Clock Init','SNTP Dead','30s Timeout','<10 TrkPts','STOW']
        offlist = ['PwrConAux','Operate','Remote','Clock Okay','SNTP Okay','No Timeout','>10 TrkPts','STOW']
        idx = range(8)
        self.set(idx,onlist,offlist)
        onlist = ['Stowing','RA-Dec Mode','Track Mode','TimOutDis']
        offlist = ['Not Stowing','Az-El mode','Setpnt Mode','TimOutEna']
        idx = [8,9,10,12]
        self.set(idx,onlist,offlist)
        # Special case of 2-bit track mode
        trkmode = self.bits[29]*2+self.bits[28]
        onlist = ['TrkArrayInit','El Online','Az Online','Track']
        offlist = ['TrkArrayOk','El Offline','Az Offline']
        if trkmode == 0:
            offlist.append('Stop')
        elif trkmode == 1:
            offlist.append('Position')
            self.bits[28] = 0   # Was 1, so set to zero
        elif trkmode == 2:
            offlist.append('Velocity')
        else:
            offlist.append('')
        idx = [16,18,20,22]
        self.set(idx,onlist,offlist)
        onlist = ['Offset On','RAOffset On','AzOffset On']
        offlist = ['Offset Stop','NoRAOffset','NoAzOffset','Refr+Pnt Cor']
        # Special case of 2-bit corrections
        cormode = self.bits[23]*2+self.bits[22]
        if cormode == 1:
            onlist.append('Refr Cor Off')
        elif cormode == 2:
            onlist.append('Pnt Cor Off')
            self.bits[22] = 1   # Was zero, so set to 1
        elif cormode == 3:
            onlist.append('All Cor Off')
        else:
            onlist.append('')
        idx = [24,26,27,28]
        self.set(idx,onlist,offlist)

def plot_creator(frame,name,number):
    name = plt.figure(number)
    sub_plot = name.add_subplot(111)
    sub_plot.grid()
    canvas = FigureCanvasTkAgg(name, frame)
    #canvas.show()
    #canvas.get_tk_widget().pack(side=TOP, expand=1)
    toolbar = NavigationToolbar2Tk(canvas, frame)
    #toolbar.update()
    #canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)

    return sub_plot, canvas

app = App()

mainloop()
