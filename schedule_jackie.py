#!/usr/bin/env python
'''
   Main application for the EOVSA Schedule, which sends commands to the
   Array Control Computer according to a timed schedule and runs various
   routines locally to generate source coordinates, track tables, uvw
   and delay information.'''
#  ____      _             _        _              
# / ___| ___| |__  ___  __| |_   _ | | ___  
# \___ \/ __| '_ \/ _ \/ _` | | | || |/ _ \
#  ___) |(__| | | | __/ (_| | |_| || |  __/
# |____/\___|_| |_\___|\__,_|\__,_||_|\___|
# 
import os, signal
os.chdir('/home/sched/Dropbox/PythonCode/Current')
from Tkinter import *
import ttk
from tkMessageBox import *
from tkFileDialog import *
from ftplib import FTP
import urllib2
import util
import threading
import subprocess
import roach
from eovsa_tracktable import *
from eovsa_array import *
from eovsa_lst import eovsa_ha
from math import pi
from readvla import *
from scan_header import *
from gen_schedule_sf import *
import stateframe
from aipy.phs import PointingError
import corr, time, numpy, socket, struct, sys
import ephem
import eovsa_cat
from eovsa_visibility import scan_visible

# Ensure that output to "terminal" goes to log file.
sys.stdout = open('/tmp/schedule.log','w')

# Show the scanheader dictionary as undefined
global sh_dict, sf_dict
sh_dict = {}
sf_dict = {}

#============================
class FuncThread(threading.Thread):
    ''' Defines a class for passing arguments to a function that will run on a
        separate thread.  Was used for call to execute_cmds(), but was
        unsuccessful so NOT CURRENTLY IN USE.
    '''
    def __init__(self, target, *args):
        self._target = target
        self._args = args
        threading.Thread.__init__(self)
    def run(self):
        self._target(*self._args)

#============================
def init_scanheader_dict(version=0.0):
    ''' Create the initial scan header dictionary, with reasonable defaults.
        Entries will be overridden before creating the scan_header file using
        scan_header.py
    '''
    global sh_dict

    d = util.datime()
    mjd_ = d.get()                    # Get mjd of datime() object
    timestamp = (mjd_ - 16480)*86400.  # Convert mjd to LabVIEW timestamp

    aa = eovsa_array()
    aa.date = str(aa.date)[:11]+'00:00'  # Set date to 0 UT
    mjd0 = aa.date + 15019.5   # Convert from ephem date to mjd

    # Read delay center file from ACC parm directory
    try:
        dlafile = urllib2.urlopen('ftp://acc.solar.pvt/parm/delay_centers.txt',timeout=1)
        dcen = []
        for line in dlafile.readlines():
            if line[0] != '#':
                # Skip comment lines and take second number as delay center [ns]
                dcen.append(float(line.split()[1]))
    except:
        print util.datime().get('str'),'ACC connection for delay centers timed out.'
        dcen = [0]*16
    
    sh_dict = {'timestamp': timestamp,
               'Version': version,
               'project':'Normal Observing',
               'operator':'Kjell Nelin',
               'comments':'None',
               'version':['1.0.0','1.0.0'],
               'nants':16,
               'antlist':numpy.array([1,0,0,0,0,0,2,4,0,0,0,0,0,0,0,3]),  # 1,7,16,8 for prototype
               'antpos': aa,
               'pbfwhm': [1.22*30*180*3600./210./pi]*13 + [1.22*30*180*3600./2700./pi]*2 + [0.0],
               'mount': [1]*13 + [2]*2 + [0],
               'scan_type':'test',
               'source_id':None,
               'track_mode':'PLANET',
               'epoch':'DATE',
               'dra': 0.0, 'ddec': 0.0, 'ha': 0.0, 'dha': 0.0,
               'sun_info': mjd0,    # Solar P-angle, B0-angle, and radius will be calculated for this date
               'pol': [-1,-2,0,0],  # Default RR, LL polarization
               'max_file_size':100, # MB, max IDB file size
               'intval': 20,
               'nintval': 50,
               'dlacen': numpy.array(dcen),
               'time_at_acc0':util.datime(),
               'roach_brd_clk':[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],
               'katadc':[{},{},{},{},{},{},{},{}]}   # Empty dictionaries (one for each ROACH) for katadc sensor data

#============================
def init_sched_dict():
    '''Create the dictionary that represents the schedule portion of the
       stateframe.  These entries are overwritten once per second as long
       as the schedule is running.
    '''
    global sf_dict
    d = util.datime()
    mjd_ = d.get()                    # Get mjd of datime() object
    timestamp = (mjd_ - 16480)*86400.  # Convert mjd to LabVIEW timestamp

    # Create schedule dictionary
    sf_dict = {'timestamp': timestamp,
               'timestamp1': timestamp,
               'scan_state': -1,
               'phase_tracking': False,
               'uvw': numpy.array([[0.0,0.0,0.0]]*16),
               'uvw1': numpy.array([[0.0,0.0,0.0]]*16),
               'delay': numpy.zeros(16),
               'delay1': numpy.zeros(16),
               'SolPwr':[{},{}],
               'delays':[{},{},{},{},{},{},{},{}],    # Empty dictionaries (one for each ROACH) for ROACH delays
               'sensors':[{},{},{},{},{},{},{},{}]}   # Empty dictionaries (one for each ROACH) for ROACH sensor data

#============================
def set_uvw(aa,d,srcname):
    '''Set u,v,w coordinates and delays (=-w) [ns] for antenna array aa
       at time specified by datime() object d (1 s in advance of current time).  
       The source name must exist in the source catalog (aa.cat) attached to aa.
    '''
    global sf_dict, sh_dict
         
    if srcname is None:
        init_sched_dict()
        return
    mjd_ = d.get()                         # Get mjd of datime() object (t+1)
    aa.date = mjd_ - 15019.5               # Convert mjd to ephem time
#    try:
#        # Newly calculate coordinates for geosat and add it to source list
#	sat=sf_dict['geosat']
#	sat.compute(aa)
#	geosat=aipy.amp.RadioFixedBody(sat.ra,sat.dec,name=sat.name)
#        # This will either add a new source or replace the existing one
#	aa.cat.add_srcs(geosat)
#    except:
#        # Probably there is no geosat defined.  Just continue.
#        pass
 
    # Fill in scan header dictionary, so that it will be ready at any time for writing
    # at the start of a scan.  These should be for the current time rather than 1 s in
    # the future, since we have not yet "computed" for new time.
    src = aa.cat[srcname]
    sh_dict['ra'] = src.ra
    sh_dict['dec'] = src.dec
    sh_dict['ha'] = eovsa_ha(src)

    cat = aa.cat                           # Make a copy of catalog
    cat.compute(aa)                        # Compute for time t+1
    sf_dict['uvw'] = sf_dict['uvw1'].copy()       # Transfer previously calculated uvw (time t)
    sf_dict['delay'] = sf_dict['delay1'].copy()   # Transfer previously calculated delay (time t)
    sf_dict['timestamp'] = sf_dict['timestamp1']  # Transfer previous timestamp (time t)
    timestamp = (mjd_ - 16480)*86400.      # Convert mjd to LabVIEW timestamp (t+1)
    sf_dict['timestamp1'] = timestamp
    uvw = sf_dict['uvw1']                  # Prepare for next uvw (t+1)
    # Loop over antennas to generate uvw of each relative to antenna 1 (index 0)
    try:
        for i in range(16):
            uvw[i] = aa.gen_uvw(0,i,src=cat[srcname]).transpose()[0,0]
        sf_dict['phase_tracking'] = True
    except PointingError:
        # Source is below horizon
        uvw = numpy.array([[0.0,0.0,0.0]]*16)
        sf_dict['phase_tracking'] = False
    # Store result for t+1 stateframe dictionary
    sf_dict['uvw1'] = uvw
    # Store corresponding delays (= -w coordinate)
    sf_dict['delay1'] = -uvw[:,2]

#============================
def mjd(line=None):
    # Returns the MJD for the time in a line of the schedule,
    # or the current time if no line is supplied.
    d = util.datime()
    if line:
        d.set(line[:19],'str')
    return d.get()

#============================
class App():
    '''Main application
    '''
    def __init__(self):
        '''Create the user interface and initialize it
        '''
        self.root = Tk()
        self.root.geometry("+100+0")

        # Check whether a schedule is already running.  Only one is allowed.
        self.mypid = str(os.getpid())
        out = subprocess.check_output(["pidof","python"]).split()
        for i in range(len(out)):
            if self.mypid != out[i]: 
                if subprocess.check_output(["ps","-lfp",out[i]]).find('schedule.py') != -1:
                    showerror('Error','Another schedule is already running (PID '+out[i]+')')
                    exit()
        # Set function to handle signal alarm, if it should go off.  It is set for 5 s in inc_time().
        signal.signal(signal.SIGALRM, self.wake_up)

        # Read ACC.ini file to get values for global variables
        # binsize, xmlpath, scdport and sfport
        try:
            self.accini = stateframe.rd_ACCfile()
        except urllib2.URLError:
            showerror('Error','ACC unreachable, or ACC.ini file does not exist\n'
                                     +'Cannot continue.')
            exit()

        # Set window title from command-line argument and
        # set global variable for checking it.  The window
        # title will be used to define which subarray.
        
        global title
        if len(sys.argv) > 1:
           title = sys.argv[1]
        else:
           # Default to Subarray1
           title = 'Subarray1'
        self.root.wm_title('Schedule for '+title)
        timeframe = Frame(self.root)
        timeframe.pack()

        self.menu = Menu(self.root)

        filemenu = Menu(self.menu, tearoff = 0)
        self.menu.add_cascade(label = 'File', menu = filemenu)
        filemenu.add_command(label = 'New', command = self.New)
        filemenu.add_command(label = 'Save', command = self.Save)
        filemenu.add_command(label = 'Open', command = self.Open)

        self.root.config(menu = self.menu)

        self.error = ''   # Error string to be added to source on front panel

        pageframe = Frame(self.root)
        pageframe.pack(expand = 1, fill = BOTH)
        # Attempt to add a tab
        self.nb = ttk.Notebook(pageframe)
        self.nb.pack(fill='both',expand='yes')

        fmain = Frame()
        self.nb.add(fmain,text='Main')

        toolbar = Frame(fmain)
        
        self.var = IntVar()
        R0 = Radiobutton(toolbar, text='D', variable=self.var, value=0)
        R0.pack(side = LEFT)

        R1 = Radiobutton(toolbar, text='H', variable=self.var, value=1)
        R1.pack(side = LEFT)

        R2 = Radiobutton(toolbar, text='M', variable=self.var, value=2)
        R2.pack(side = LEFT)

        R3 = Radiobutton(toolbar, text='S', variable=self.var, value=3)
        R3.pack(side = LEFT)
        
##        self.cmd_not = BooleanVar()
##        self.CB = Checkbutton(toolbar, text="Disable Commands", variable=self.cmd_not)
##        self.CB.pack(side=LEFT, expand=0, fill=BOTH)

        # Widget for current source and phase tracking state
        self.source = Label(toolbar,text='            Source: None        '+'Phase Tracking: False')
        self.source.pack(side=LEFT)
        toolbar.pack(side=TOP, fill=X)

        self.label = Label(timeframe,text='',bg="yellow",font="Helvetica 16 bold")
        self.label.pack()
        #textframe = Frame(self.root)
        textframe = Frame(fmain)
        textframe.pack(expand = True, fill = BOTH)

        # Main Listbox widget for Macro commands
        Lframe = Frame(textframe)
        Lframe.pack(expand = False, fill=Y,side=LEFT)
        self.L = Listbox(Lframe,selectmode=EXTENDED,width=45)
        self.Sx = Scrollbar(Lframe,orient=HORIZONTAL,command=self.xview)
        self.L.config( xscrollcommand = self.Sx.set)
        self.L.bind('<<ListboxSelect>>',self.display_ctl)

        # Associated Listbox widget for status (Waiting..., Running..., or Done)
        self.status = Listbox(Lframe,width=8)
        self.Sx.pack(side=BOTTOM, fill=X)
        self.L.pack(side=LEFT, expand=True, fill = Y)

        self.S = Scrollbar(Lframe)
        self.S.pack( side = LEFT, fill = BOTH)
        self.status.pack(side = LEFT, fill = Y, expand = True)
        self.S.config (command = self.yview)
        self.L.config( yscrollcommand = self.S.set)
        self.status.config( yscrollcommand = self.S.set)
        
        # Atomic Command Listbox
        Rframe = Frame(textframe)
        Rframe.pack(expand = True, fill=BOTH,side=LEFT)
        self.L2 = Listbox(Rframe,selectmode=SINGLE,width=25)
        self.L2.pack(side=LEFT, fill = BOTH, expand = True)
        self.L.atomlist = self.L2   # Associate main listbox (L) with L2

        self.S2 = Scrollbar(Rframe)
        self.S2.pack(side = RIGHT, fill = BOTH)
        self.S2.config( command = self.L2.yview)
        self.L2.config( yscrollcommand = self.S2.set)

        self.downbutton = Button(fmain, text = '- 1', command=self.Decrease_cmd)
        self.downbutton.pack(side=LEFT)

        self.upbutton = Button(fmain, text = '+ 1', command=self.Increase_cmd)
        self.upbutton.pack(side=LEFT)

        self.Insert = Button(fmain,text="Insert", command=self.Insert)
        self.Insert.pack(side = LEFT)

        # Toggle controls which state the button is, Go/Stop
        self.Toggle = 1
        #Stop/Go Button
        self.B2 = Button(fmain,text='GO', width = 8, command=self.toggle_state, background='GREEN')
        self.B2.pack(side = LEFT)

        self.TodayBtn = Button(fmain,text = 'Today', command = self.Today)
        self.TodayBtn.pack(side = LEFT)

        self.ClearBtn = Button(fmain,text = 'Clear', width = 5, command=self.Clear)
        self.ClearBtn.pack(side=LEFT)

        # Entry widget, for entering a new line into the schedule
##        content = StringVar()
##        self.E1 = Entry(textvariable = content)
##        self.E1.pack( side = RIGHT)

#        self.B1 = Button(text="Create new event", command=self.write)
#        self.B1.pack( side = RIGHT)

        #self.B3 = Button(text="Raw Command", command=self.send_cmd)
        #self.B3.pack( side = RIGHT)
        #command = StringVar()
        #self.E3 = Entry(textvariable = command)
        rawlabel = Label(fmain,text='Raw Command:')
        rawlabel.pack(side=LEFT)
        E3 = Entry(fmain)
        E3.pack( side = LEFT, expand=True, fill=X)
        E3.bind('<Return>',self.send_cmd)

        d = util.datime()
        self.label.configure(text=d.get('str'))

        # Setup Project Tab
        fproj = Frame()
        fproj.pack(side=LEFT)
        self.nb.add(fproj,text='Project')
        fline1 = Frame(fproj)
        fline1.pack(side=TOP)
        tmplabel = Label(fline1,text='Project:')
        tmplabel.pack(side=LEFT)
        self.wproj = Entry(fline1,width=32)
        self.wproj.pack( side = LEFT, expand=False, fill=X)
        self.wproj.bind('<Return>',self.handle_entry)
        
        tmplabel = Label(fline1,text='Operator:')
        tmplabel.pack(side=LEFT)
        self.woper = Entry(fline1,width=16)
        self.woper.pack( side = LEFT, expand=False, fill=X)
        self.woper.bind('<Return>',self.handle_entry)
        self.woper.config(state=DISABLED)
        
        fline2 = Frame(fproj)
        fline2.pack(side=TOP)
        tmplabel = Label(fline2,text='Comment:')
        tmplabel.pack(side=LEFT)
        self.wcomm = Entry(fline2,width=50)
        self.wcomm.pack( side = LEFT, expand=False, fill=X)
        self.wcomm.bind('<Return>',self.handle_entry)

        # Make sure window can never be smaller than initially created size
        self.root.update_idletasks()
        self.root.minsize(self.root.winfo_width(),self.root.winfo_height())

        # Generate the Antenna Array object, used to calculate source coordinates,
        # uvw, delays, etc.
        self.aa = eovsa_cat.eovsa_array_with_cat()
        global sf_dict, sh_dict
        if sh_dict == {}:
            init_scanheader_dict(self.accini['version'])
        if sf_dict == {}:
            init_sched_dict()
            # Generate schedule item XML file on init
            junk = gen_schedule_sf(sf_dict)#,mk_xml=True)
         
        self.Open('solar.scd')
        
        # Establish connect to ROACHes.
        self.roaches =[]
        self.connect2roach()

        # Initialize $WAIT settings
        self.waitmode = False  # Not in $WAIT mode
        self.nextctlline = 0
        self.wait = 0
        self.solpwr = [{},{}]  # Empty solar power dictionary pair
        self.sensors = [{},{},{},{},{},{},{},{}]  # Empty ROACH sensor dictionaries
        self.delays = [{},{},{},{},{},{},{},{}]  # Empty ROACH delays dictionaries
        self.w = {} # Empty weather dictionary

        # Start the clock ticking
        self.prev = time.time()
        self.root.after(1,self.inc_time)

    #============================
    def xview(self, *args):
        #Definition for the scrollbar to control two windows and the same time.
        self.L.xview(*args)
        self.status.xview(*args)

    #============================
    def yview(self, *args):
        #Definition for the scrollbar to control two windows and the same time.
        self.L.yview(*args)
        self.status.yview(*args)

    #============================
    def connect2roach(self):
        # This will eventually be a loop over all active ROACHes, and must
        # tolerate a missing ROACH
        roach_ips = ('roach1.solar.pvt','roach2.solar.pvt','roach3.solar.pvt')
        if len(self.roaches) != 0:
            # Some roaches are already connected, so stop them and reconnect
            for r in self.roaches:
                if r.fpga: r.fpga.stop()
        self.roaches = []
        for roach_ip in roach_ips:
            # Make connection to ROACHes
            rnum = int(roach_ip[5:6])-1
            self.roaches.append(roach.Roach(roach_ip, self.accini['boffile']))
            if self.roaches[-1].msg == 'Success':
                print roach_ip,'is reachable'
                try:
                    self.roaches[-1].brd_clk = self.roaches[-1].fpga.est_brd_clk()
                    print roach_ip,'clock is',self.roaches[-1].brd_clk
                    if self.roaches[-1].brd_clk < 199 or self.roaches[-1].brd_clk > 201:
                        print roach_ip,'clock NOT 200 MHz.'
                except:
                    print roach_ip,'could NOT read FPGA clock speed. Will mark unreachable'
                    self.roaches[-1].brd_clk = 0.0
                    self.roaches[-1].fpga.stop()
                    self.roaches[-1].fpga = None
                sh_dict['roach_brd_clk'][rnum] = self.roaches[-1].brd_clk
            else:
                print roach_ip,'is unreachable!',self.roaches[-1].msg
                self.roaches.pop()

    #============================
    def wake_up(self):
        # This is called whenever the 5-second alarm goes off, indicating the
        # process is stuck in sk_wait.  We simply send ourselves a SIGINT (ctrl-C),
        # which should hopefully do it, but we should also log the fact by setting
        # a flag in the self object.
        self.error = 'The 5-s-alarm went off!'
        # Try to reestablish connection to the ROACHes, and set self.fpga accordingly
        # This will keep dla2roach() from hanging.
        self.connect2roach()
        os.kill(self.mypid, signal.SIGINT)

    #============================
    def send_cmd(self,event):
        command = event.widget.get()
        s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        try:
            s.connect((self.accini['host'],self.accini['scdport']))
            s.send(command)
            time.sleep(0.01)
            s.close()
        except:
            pass

    #============================
    def handle_entry(self,event):
        w = event.widget
        command = w.get()
        if w == self.wproj:
            print 'Project is:',command
        elif w == self.woper:
            print 'Observer is:',command
        elif w == self.comm:
            print 'Comment is:',command
        else:
            print 'unknown widget'

    #============================
    def display_ctl(self,event):
        '''Callback for when user clicks in the Macro command window
        '''
        w = event.widget
        if self.Toggle == 0:
            # Schedule is running, so clear selection and do nothing
            w.selection_clear(0,END)
        else:
            sel = map(int, w.curselection())
            if len(sel) == 1:
                line = w.get(sel[0])
                cmds = line[20:].split(' ')
                f2 = open(cmds[0].rstrip()+'.ctl')
                w.atomlist.delete(0,END)
                for ctlline in f2.readlines():
                    w.atomlist.insert(END,ctlline.rstrip('\n'))
                f2.close()

    #============================
    def Clear(self):
        #Button to clear the status and the highlight. 
        self.L.selection_clear(0,END)
        self.status.configure( state = NORMAL)
        self.status.delete(0,END)
        self.state=['']*self.lastline
        for x in self.state:
            self.status.insert(END,x)        
        self.L.itemconfig(self.curline,background="white")
        self.curline = 0
        self.status.configure( state = DISABLED)

    #============================
    def Open(self, filename=None):
        ''' To open a new file, delete the contents of lines, and
            the text widget. Then proceed to populate them again,
            checking any PHASECAL lines to add them to the catalog.
        '''
        if self.Toggle == 0:
            # Schedule is running, so do nothing
            return
        self.status.configure(state = NORMAL)
        if filename is None:
            init_dir = os.getcwd()
            f = askopenfile(initialdir = init_dir, mode = 'r',
                            filetypes = [('SCD Files','*.scd'),('all files','*.*')])
            if f is None:
                # User cancelled, so do nothing.
                self.status.configure(state = DISABLED)
                return        
        else:
            f = open(filename,'r')
        try:    #takes care of an empty line, if there is one, in
                #the file being read.
            lines = f.readlines()
        except AttributeError:
            pass
        else:
            self.L.delete(0, END)
            self.curline = 0
            self.lastline = len(lines)
            for i,line in enumerate(lines):
                self.L.insert(END, line.rstrip('\n'))
                if line[20:28] == 'PHASECAL':
                    name = line[29:37]
                    # Get start and stop time for this calibrator, as ephem-compatible strings
                    dtstart = util.datime()
                    dtstart.set(mjd(line))
                    dtstop = util.datime()
                    if i == self.lastline-1:
                        # Should never happen, but just in case the PHASECAL
                        # is the last line of the schedule, assume 15 minutes duration
                        dtstop.set(mjd(line) + 15.*60./86400.)
                    else:
                        line2 = lines[i+1]
                        dtstop.set(mjd(line2))
                    tstop_str = dtstop.get('str').replace('-','/')
                    # Check visibility for source
                    try:
                        src=self.aa.cat[name]
                        visible = scan_visible(src,self.aa,dtstart,dtstop)
                        if not visible:
                            self.error = 'Warning, source '+name+' not visible at scheduled time: Schedule line '+str(i+1)
                    except:
                        self.error = 'Err: source '+name+' not found in catalog.  Not added to schedule!'

        self.state=['']*len(lines)
        self.status.configure(state = DISABLED)
        if filename is None:
            filenamelist = f.name.split('/')
            self.filename = filenamelist[len(filenamelist)-1:][0]
        else:
            self.filename = filename
        
    #============================
    def New(self):
        # New option creates a new table with predetermined content.
        if self.Toggle == 0:
            # Schedule is running, so do nothing
            return
        self.status.configure( state = NORMAL)
        d = util.datime()
        time = d.get('str')
        self.L.delete(0, END)
        self.L.insert(0,time[:19] + ' ' + 'SUN')
        self.L.insert(1,time[:19] + ' ' + 'REWIND')
        self.curline = 0
        self.lastline = 2
        self.status.configure( state = DISABLED)
        
    #============================
    def Save(self):
        ''' The Save button will save the file as a text file, in a folder
            specified by the user. If the file exists, the program will ask
            the user if he wants to replace the file.
        '''
        if self.Toggle == 0:
            # Schedule is running, so do nothing
            return
        try:
            #The Exception is to allow the user to cancel the save.
            init_dir = os.getcwd()
            fileout = asksaveasfile(initialdir = init_dir, 
                                    initialfile = self.filename, mode = 'w',
                                    filetypes = [('SCD Files','*.scd'),
                                    ('all files','*.*')])
            for i in range(self.lastline):
                fileout.write(self.L.get(i)+'\n')
            fileout.close()
        except AttributeError:
            pass

    #============================
    def adjust_selection(self,sel,delt):
        if len(sel) == 0:
            sel = range(self.lastline)
        d = util.datime()
        for i in sel:
            line = self.L.get(i)
            time_ = mjd(line) + delt
            d.set(time_)
            tstr = d.get('str')
            temp1 = tstr[:19] + ' ' + line[20:]
            self.L.delete(i)
            self.L.insert(i,temp1)
        for i in sel:
            # Selection is getting unset, so reset it.
            self.L.selection_set(i)

    #============================
    def Decrease_cmd(self):
        ''' The decrease button will decrease the selected time by one. The options
            are Hours, Minutes and Seconds. If no selection is made
            it will create an error.
        '''
        sel = map(int, self.L.curselection())  # list of line indexes selected
        if self.var.get() == 0:
            one_day = 1.
            self.adjust_selection(sel,-one_day)
        elif self.var.get() == 1:
            one_hour = 1./24
            self.adjust_selection(sel,-one_hour)
        elif self.var.get() == 2:
            one_minute = 1./1440
            self.adjust_selection(sel,-one_minute)
        elif self.var.get() == 3:
            one_second = 1./86400
            self.adjust_selection(sel,-one_second)

    #============================
    def Increase_cmd(self):
        ''' The increase button will increase the selected time by one. The options
            are Hours, Minutes and Seconds. If no selection is made
            it will create an error.
        '''
        sel = map(int, self.L.curselection())  # list of line indexes selected
        if self.var.get() == 0:
            one_day = 1.
            self.adjust_selection(sel,one_day)
        elif self.var.get() == 1:
            one_hour = 1./24
            self.adjust_selection(sel,one_hour)
        elif self.var.get() == 2:
            one_minute = 1./1440
            self.adjust_selection(sel,one_minute)
        elif self.var.get() == 3:
            one_second = 1./86400
            self.adjust_selection(sel,one_second)


    #============================
    def Today(self,d=None):
        ''' The button Today will change the schedule to start on the current date.
            This will NOT change times. Lines that start a day after the start line
            of the schedule will start on today + 1.
        '''
        if d is None:
            d = util.datime()
        # Determine how many days from date of first line to today
        line = self.L.get(0)
        days = int(d.get()) - int(mjd(line))
        print 'Adding ',days,'days.'
        for i in range(self.lastline):
            line = self.L.get(i)
            linemjd = mjd(line) + days
            d.set(linemjd)
            time_ = d.get('str')[:10]
            line = time_ + line[10:]
            self.L.delete(i)
            self.L.insert(i,line)
#            self.logf.write(str(i)+' '+line+' '+str(linemjd)+'\n')
        self.curline = 0
        # If this is the standard solar.scd file, do an auto-generate for date of first line
        if self.filename == 'solar.scd':
            line = self.L.get(0)
            linemjd = mjd(line)
            d.set(linemjd)
            self.autogen(d)

    #============================
    def autogen(self,d):
        # Auto-generate the standard solar schedule
        dtemp = util.datime()
        # Determine sunrise, sunset times for this day
        mjdrise, mjdset = suntimes(d)
        # First solar line starts at sunrise
        dtemp.set(mjdrise)
        line = self.L.get(1)
        line = dtemp.get('str')[:19] + line[19:]
        self.L.delete(1)
        self.L.insert(1,line)
        # First calibrator line starts 15 min earlier
        dtemp.set(mjdrise - 15.*60./86400.)
        line = self.L.get(0)
        line = dtemp.get('str')[:19] + line[19:]
        self.L.delete(0)
        self.L.insert(0,line)
        # Last calibrator line starts at sunset
        dtemp.set(mjdset)
        line = self.L.get(self.lastline-2)
        line = dtemp.get('str')[:19] + line[19:]
        self.L.delete(self.lastline-2)
        self.L.insert(self.lastline-2,line)
        # Last (REWIND) line starts 15 min later
        dtemp.set(mjdset + 15.*60./86400.)
        line = self.L.get(END)
        line = dtemp.get('str')[:19] + line[19:]
        self.L.delete(END)
        self.L.insert(END,line)

        # Read calibrator database
        cal = readvlacaldb()
        # Find sources within 15-35 deg of Sun (also returns antenna array aa)
        srclistnarrow, aa = findcal(cal,dt=d,dtheta=[15,35])
        # Early and late in day, need a wider search window, 15-55 degrees of the Sun
        srclistwide, aa = findcal(cal,dt=d,dtheta=[15,60])
        # Sort by flux density (both narrow and wide)
        fluxnarrow = []
        for src in srclistnarrow:
            fluxnarrow.append(src.mag)
        fsortnarrow = sorted(fluxnarrow,reverse=True)
        fluxwide = []
        for src in srclistwide:
            fluxwide.append(src.mag)
        fsortwide = sorted(fluxwide,reverse=True)

        # Now go through calibrator lines one by one and
        # select appropriate source in VLA cal database
        for i in range(self.lastline):
            line = self.L.get(i)
            if line[20:28] == 'PHASECAL':
                dtstart = util.datime()
                dtstart.set(mjd(line))
                dtstop = util.datime()
                if i == self.lastline-1:
                    # Should never happen, but just in case the PHASECAL
                    # is the last line of the schedule, assume 15 minutes duration
                    dtstop.set(mjd(line) + 15.*60./86400.)
                else:
                    line2 = self.L.get(i+1)
                    dtstop.set(mjd(line2))
                # Loop over calibrators, highest flux first
                for f in fsortnarrow:
                    # These are 15-min observations, so make sure calibrator
                    # is visible and will remain visible
                    idx = fluxnarrow.index(f)
                    jys = fluxnarrow[idx]
                    src = srclistnarrow[idx]
                    visible = scan_visible(src,self.aa,dtstart,dtstop)
                    if visible:
                        # Take first visible source, since it will be the one
                        # with the highest flux in sorted list
                        break
                if visible:
                    line = line[:29] + src.name + line[37:]
                    # If this source is not already in the source list, add it
                    try:
                        blah = self.aa.cat[src.name]
                    except KeyError:
                        # This source is not already in the list, so add it.
                        # Since src is an ephem.FixedBody, it must be converted
                        # to an aipy.phs.RadioFixedBody
                        radiosrc = aipy.amp.RadioFixedBody(src.a_ra,src.a_dec,name=src.name,jys=jys,mfreq=1.0)
                        radiosrc.compute(aa)
                        radiosrc.name = src.name
                        self.aa.cat.add_srcs(radiosrc)
                else:
                    # No source found for narrow window, so try the wide one
                    for f in fsortwide:
                        # These are 15-min observations, so make sure calibrator
                        # is visible and will remain visible
                        idx = fluxwide.index(f)
                        jys = fluxwide[idx]
                        src = srclistwide[idx]
                        visible = scan_visible(src,self.aa,dtstart,dtstop)
                        if visible:
                            # Take first visible source, since it will be the one
                            # with the highest flux in sorted list
                            break
                    if visible:
                        line = line[:29] + src.name + line[37:]
                        # If this source is not already in the source list, add it
                        try:
                            blah = self.aa.cat[src.name]
                        except KeyError:
                            # This source is not already in the list, so add it.
                            # Since src is an ephem.FixedBody, it must be converted
                            # to an aipy.phs.RadioFixedBody
                            radiosrc = aipy.amp.RadioFixedBody(src.a_ra,src.a_dec,name=src.name,jys=jys,mfreq=1.0)
                            radiosrc.compute(aa)
                            radiosrc.name = src.name
                            self.aa.cat.add_srcs(radiosrc)
                    else:
                        # No source found for wide, so mark line as a failure
                        line = line[:29] + 'No Src!!' + line[37:]
                        print 'No source after searching the following sources: '
                        for f in fsortwide:
                            idx = fluxwide.index(f)
                            jys = fluxwide[idx]
                            src = srclistwide[idx]
                            visible = scan_visible(src,self.aa,dtstart,dtstop)
                            print src.name, jys, visible, src.ra, src.dec, '     ', src.az, src.alt
                self.L.delete(i)
                self.L.insert(i,line)

    #============================
    def Insert(self):
        ''' Insert button will take the text written in the entry widget and
            insert it at the insertion indicator
        '''
        self.status.configure( state = NORMAL )
        sel = map(int, self.L.curselection())
        if len(sel) == 1:
            index = sel[0]
            self.content = self.E1.get().upper()
            if self.content:
                line = self.L.get(index)
                line = line[:20] + self.content
                self.L.insert(index,line)
                self.lastline += self.lastline
        else:
            # If there is no string do not do anything.
            pass
        self.status.insert(END, '')
        self.status.configure( state = DISABLED)

    #============================
    def toggle_state(self):
        now = mjd()
        if self.Toggle  == 0:    #Stop was pressed
            self.Toggle = 1
            self.B2.configure(text = 'Go')
            self.B2.configure(background = 'GREEN')
            self.downbutton.configure(state = NORMAL)
            self.upbutton.configure(state = NORMAL)
##            self.B1.configure(state = NORMAL)
            self.Insert.configure(state = NORMAL)
            self.ClearBtn.configure(state = NORMAL)
            self.TodayBtn.configure(state = NORMAL)
            self.L.configure( state = NORMAL)
            

        else:    #Go was Pressed
            self.Toggle = 0
            self.status.configure( state = NORMAL)
            self.L.selection_clear(0,END)
            self.B2.configure(text = 'STOP')
            self.B2.configure(background = 'RED')
            self.curline = self.lastline-1  # Will be overridden if another line is determined to be current line

            # Go through the status lines to find the expired ones.
            for i in range(self.lastline):
##                if self.mjd[i] >= now:
                line = self.L.get(i)
                if mjd(line) >= now:
                    # Lines in the future
                    if i == 0:
                        # If the first line, mark it Waiting
                        self.status.delete(i)
                        self.status.insert(i,'Waiting...')
                        self.curline = i
                        self.L.itemconfig(i,background="orange")
                    else:
                        # If not the first line, mark earlier line for starting.
                        self.curline = i-1
                        self.status.delete(self.curline)
                        self.status.insert(self.curline,'Started...')
                        self.L.itemconfig(self.curline,background="orange")
                    break
                else:
                    # Mark lines in the past as Skipped.
                    self.status.delete(i)
                    self.status.insert(i,'Skipped')
            # Clear remaining lines
            for i in range(min(self.curline+1,self.lastline),self.lastline):
                self.status.delete(i)
                self.status.insert(i,'')
            if self.curline == (self.lastline-1):
                self.status.delete(self.curline)
                self.status.insert(self.curline,'Started...')

            # Find the file associated with the Macro command on the current 
            # line and fill in the L2 Listbox
            line = self.L.get(self.curline)
            cmds = line[20:].split(' ')
            f2 = open(cmds[0].rstrip()+'.ctl')
            self.L2.delete(0,END)
            lines = f2.readlines()
            for ctlline in lines:
                # Check for hash mark (#) in line other than first character
                # (hash mark in first character means a comment)
                if '#' in ctlline[1:]:
                    # We have a substitution to do
                    ihash = ctlline[1:].find('#')+1
                    i = int(ctlline[ihash+1:ihash+2])
                    ctlline = ctlline[:ihash]+cmds[i]+ctlline[ihash+2:]
                self.L2.insert(END,ctlline.rstrip('\n'))
            f2.close()

            self.L.see(min(self.curline+5,END))
            self.status.see(min(self.curline+5,END))            
            self.downbutton.configure(state = DISABLED)
            self.upbutton.configure(state = DISABLED)
##            self.B1.configure(state = DISABLED)
            self.Insert.configure(state = DISABLED)
            self.ClearBtn.configure(state = DISABLED)
            self.TodayBtn.configure(state = DISABLED)
            self.status.configure(state = DISABLED)
                    

    #============================
    def inc_time(self):
        global sf_dict, sh_dict
        self.status.configure(state=NORMAL)
        # First set the timer to wake us on the next second
        d = util.datime()
        tdif = int((d.dt.microsecond)/1000.)
        self.root.after(1000 - tdif, self.inc_time)
        # Update the clock
        self.label.configure(text=d.get('str')[:19])

        # Set an alarm for 5 seconds.  If the process hangs for more than that, we will send ourselves
        # a SIGINT via the Callback self.wake_up(), which should recover from sk_wait hang up
        signal.alarm(5)  # This will reset alarm if it has not gone off yet
        tnow = time.time()
        self.telapsed = tnow - self.prev
        self.prev = tnow
        telapsed = int(self.telapsed*1000)
        if telapsed == 999 or telapsed == 1000:
            pass
        else:
            # If elapsed time is not nominal (e.g. 999 or 1000), write it to log file.
            print util.datime().get('str'),str(int(self.telapsed*1000))
            sys.stdout.flush()  # Flush stdout (/tmp/schedule.log) once per second so we can see the output.

        # Update phase tracking (u,v,w and delays)
        srcname = sh_dict['source_id']
        try:
            # Generate a datime() object at exactly the next upcoming second (time t+1)
            dsec = util.datime()
            dsec.set(dsec.get()  + (1 - dsec.dt.microsecond/1000000.)/86400.)
            src = self.aa.cat[srcname]        # This causes KeyError if source is not found
        except KeyError:
            # The current scan header source ID is not in the source catalog
            srcname = None
        # Debug info, simply logs that we have started this procedure
        #sys.stdout.write('+')
        #sys.stdout.flush()  # Flush stdout (/tmp/schedule.log) once per second so we can see the output.
        set_uvw(self.aa,dsec,srcname)
        #sys.stdout.write('-')
        #sys.stdout.flush()  # Flush stdout (/tmp/schedule.log) once per second so we can see the output.
        self.source.configure(text='    Source: '+(str(srcname)+'            ')[:12]
                                  +'Phase Tracking: '
                                  +str(bool(sf_dict['phase_tracking'])) + '    '+self.error)
        self.error = ''
        # Send integer delays to ROACHes
        self.dla2roach()

        # Update weather information in sf_dict (reads from OVRO weather station)
        self.w = stateframe.weather()
        sf_dict.update(self.w)

        # Once per minute, update the information from the Solar Power station(s)
        if d.dt.second == 0:
            self.solpwr = stateframe.rd_solpwr()
        sf_dict.update({'SolPwr':self.solpwr})

        # Read ROACH sensor data, but only one each minute, staggered over different times
        # since for all 8 ROACHes this can take more than 0.5 s
        for i in range(len(self.roaches)):
            if d.dt.second == 5*i+5:
                r = self.roaches[i]
                rnum = int(r.roach_ip[5:6])-1
                if r.fpga:
                    r.get_sensor_dict()
                    if r.msg == 'Success':
                        self.sensors[rnum].update(r.sensors)
                    else:
                        self.sensors[rnum] = {}
                else:
                    self.sensors[rnum] = {}

        # Read ROACH delay values
        for i in range(len(self.roaches)):
            r = self.roaches[i]
            rnum = int(r.roach_ip[5:6])-1
            if r.fpga:
                r.get_delays()
                if r.msg == 'Success':
                    delays = dict(zip(['dx0','dy0','dx1','dy1'],r.delays))
                else:
                    delays = dict(zip(['dx0','dy0','dx1','dy1'],[0,0,0,0]))
            else:
                delays = dict(zip(['dx0','dy0','dx1','dy1'],[0,0,0,0]))
            self.delays[rnum].update(delays)        

        for i in range(8):
            sf_dict['sensors'][i].update(self.sensors[i])
            sf_dict['delays'][i].update(self.delays[i])

        # Get current stateframe (from ACC) and update sf_dict with Azimuth, Elevation, TrackFlag 
        # and parallactic angle information from it (all in degrees!)
        data, msg = stateframe.get_stateframe(self.accini)
        if msg == 'No Error':
            version = struct.unpack_from('d',data,8)[0]   # Get stateframe version from data
            if version > 0.0 and version != self.accini['version']:
                # The version number of the stateframe data has changed, so we need to reread
                # the ACC ini file (which will read a new stateframe.xml file and give us a new
                # sf dictionary.
                self.accini = stateframe.rd_ACCfile()
            sf = self.accini['sf']
            sf_dict.update(stateframe.azel_from_stateframe(sf,data))
            # Flag unused antennas as not tracking
            sf_dict['TrackFlag'] = (sf_dict['TrackFlag']) & (sh_dict['antlist'] != 0)
        else:
            self.error = msg

        # Create schedule part of stateframe from sf_dict
        fmt, buf, sched_xmlfile = gen_schedule_sf(sf_dict)

        # Open socket to ACC
        s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)

        try:
            # Try to connect and send schedule items of stateframe to ACC
            # Uses "schedule" command and encloses data buffer in square brackets
            time.sleep(0.01)
            s.connect((self.accini['host'],self.accini['scdsfport']))
            s.settimeout(0.5)
            s.sendall(buf)
            time.sleep(0.02)
            s.close()
        except socket.timeout: 
            print util.datime().get('str'),'Socket time-out when writing sched stateframe to ACC'
            s.close()
        except:
            self.error = 'Err: Cannot write sched stateframe to ACC'
            
        if self.Toggle == 0:
            # Schedule is in the GO state.
            # First check if the current line needs to be started or stopped
            line = self.L.get(self.curline)
            status = self.status.get(self.curline)
            now = mjd()
            if (mjd(line) <= now and status == 'Waiting...') or status == 'Started...':
                # This line has not been started, so do so now
                self.status.delete(self.curline)
                self.status.insert(self.curline,'Running...')
                self.L.itemconfig(max(self.curline-1,0),background="white")
                self.L.itemconfig(self.curline,background="orange")
                self.L.see(min(self.curline+5,END))
                self.status.see(min(self.curline+5,END))
#******
                # Change to spawn this task as non-blocking function
                # but make sure it returns, or there is some semaphore
                # behavior with error checking
                self.execute_cmds()
                #t1 = FuncThread(execute_cmds,self)
            elif status == 'Running...':
                if self.waitmode:
                    # If $WAIT is currently in force, decrement self.wait
                    # When self.wait = 0, continue executing commands starting with
                    # line self.nextctlline, which should be line following $WAIT
                    self.wait -= 1
                    if self.wait == 0:
                        self.execute_cmds()
                        # After executing, turn off waitmode
                        self.waitmode = False
                        self.nextctlline = 0
                nextline = self.L.get(self.curline+1)
                if mjd(nextline) <= now:
                    # Next line should be running
                    self.status.delete(self.curline)
                    self.status.insert(self.curline,'Done')
                    self.curline += 1
                    self.status.delete(self.curline)
                    self.status.insert(self.curline,'Running...')
                    self.L.itemconfig(self.curline-1,background="white")
                    self.L.itemconfig(self.curline,background="orange")
                    self.L.see(min(self.curline+5,END))
                    self.status.see(min(self.curline+5,END))
#******
                    # Change to spawn this task as non-blocking function
                    # but make sure it returns, or there is some semaphore
                    # behavior with error checking
                    self.execute_cmds()
                    #t1 = FuncThread(execute_cmds,self)
        self.status.configure(state=DISABLED)
        # Debug info, simply logs that we have exited this procedure
        #sys.stdout.write('-')
        #sys.stdout.flush()  # Flush stdout (/tmp/schedule.log) so we can see this '-'.

    #============================
    def dla2roach(self):
        '''Set integer delays and send to all ROACHes for next second, to be
           ready for next 1 PPS.  Ant 1 delay is fixed at 5000 steps (delay
           depends on step size), and all other antennas are set to
           5000 + t_cen[i] + t_geom[i]/step_size, the latter being calculated
           from sf_dict['delay1'][i] (nsec).  The step_size is 1/f_ADC,
           where f_ADC is ADC clock frequency, in GHz.
        '''
        global sh_dict, sf_dict
        # Delay centers.  These are read from file 'delay_centers.txt' in the
        # ACC parm directory by init_scan_dict().
        # We assume for now that the X channel and Y channel delays are the same.
        # The 0.800 is because current clock frequency is 800 MHz.  This should be a
        # parameter. The 5000 is to center delays mid-range (ROACHes have a range
        # of 10,000), and provide integer delay.
        dlax = ((sh_dict['dlacen'] + sf_dict['delay1'])*0.800 + 5000).astype(int)
        dlay = dlax

        for r in self.roaches:
            if r.fpga:
                a1,a2 = r.ants
                r.set_delays([dlax[a1-1],dlay[a1-1],dlax[a2-1],dlay[a2-1]])
                if r.msg != 'Success':
                    self.error = r.msg+' '+r.roach_ip

    #============================
    def sequence2roach(self,sequence):
        '''Set frequency sequence on the ROACH boards, so that application of band-dependent
           coefficients is properly applied.  The sequence numbers are 0-based band numbers
        '''
        # Convert from comma-separated variables to zero-based band numbers
        bands = numpy.array(sequence.split(',')).astype('int')-1
        for r in self.roaches:
            if r.fpga:
                r.set_sequence(bands)
                if r.msg != 'Success':
                    self.error = r.msg+' '+r.roach_ip
    
    #============================
    def execute_cmds(self):
        '''Execute the atomic commands associated with the current line
           of the schedule.  First read the commands from the associated
           file and enter them into the L2 Listbox.  Then read them one
           at a time from the Listbox and execute them.
        '''
        global sf_dict, sh_dict, title
        # Get time range of this Macro command
        mjd1 = mjd(self.L.get(self.curline))
        mjd2 = mjd(self.L.get(self.curline+1))
        # Find the file associated with the Macro command on the current 
        # line and fill in the L2 Listbox
        line = self.L.get(self.curline)
        cmds = line[20:].split(' ')
        f2 = open(cmds[0].rstrip()+'.ctl')
        self.L2.delete(0,END)
        # Current options for source ID
        if cmds[0].upper() == 'SUN':
            sh_dict['source_id'] = 'Sun'
            sh_dict['track_mode'] = 'PLANET'
        elif cmds[0].upper() == 'PHASECAL':
            sh_dict['source_id'] = cmds[1]
            sh_dict['track_mode'] = 'RADEC '
        elif cmds[0].upper() == 'STARBURST':
            sh_dict['source_id'] = cmds[1]
            sh_dict['track_mode'] = 'RADEC '
            print 'Source is',cmds[1]
        elif cmds[0].upper() == 'GEOSAT':
            sh_dict['source_id'] = cmds[1].replace('_',' ')
            # These are geostationary satellites so far.  If/when we add
            # moving satellite capability, track_mode for those should be 'SATELL'
            sh_dict['track_mode'] = 'FIXED '
            try:
                f = urllib2.urlopen('http://www.celestrak.com/NORAD/elements/geo.txt',timeout=2)
                lines = f.readlines()
                for i,line in enumerate(lines):
                     if line.find(sh_dict['source_id']) == 0:
                         break
                if i < len(lines):
                    sat=ephem.readtle(lines[i],lines[i+1],lines[i+2])
                    sf_dict['geosat']=sat
                    sat.compute(self.aa)
                    geosat=aipy.amp.RadioFixedBody(sat.ra,sat.dec,name=sat.name)
                    self.aa.cat.add_srcs(geosat)
                else:
                    print 'Geosat named ',sh_dict['source_id'],'not found!'
                    sh_dict['source_id']=None
            except:
                print util.datime().get('str'),'Connection to Celestrak timed out.'
                sh_dict['source_id']=None
        else:
            sh_dict['source_id'] = None
        lines = f2.readlines()
        for ctlline in lines:
            # Check for hash mark (#) in line other than first character
            # (hash mark in first character means a comment)
            if '#' in ctlline[1:]:
                # We have a substitution to do
                ihash = ctlline[1:].find('#')+1
                i = int(ctlline[ihash+1:ihash+2])
                ctlline = ctlline[:ihash]+cmds[i]+ctlline[ihash+2:]
            self.L2.insert(END,ctlline.rstrip('\n'))
        f2.close()
        # Now read the atomic commands one at a time, executing locally
        # those starting with $, and sending the others to the ACC.
        if self.waitmode:
            # If $WAIT is in effect, start with next following line
            sline = self.nextctlline
        else:
            sline = 0
        for i in range(sline,len(lines)):
            ctlline = self.L2.get(i)
            # Send line to ACC (lines that start with '$' will be entered
            # into stateframe without execution by ACC.
            s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
            # Skip comments
            if ctlline[0] == '#':
                pass
            else:
                try:
                    s.connect((self.accini['host'],self.accini['scdport']))
                    s.send(ctlline)
                    time.sleep(0.01)
                    s.close()
                except:
                    pass
                if ctlline[0] == '$':
                    # This line starts with '$', so execute locally
                    #==== MK_TABLES ====
                    if ctlline.split(' ')[0].upper() == '$MK_TABLES':
                        cmd, fname, src = ctlline.split(' ')
                        tbl = make_tracktable(src,self.aa,mjd1,mjd2)
                        # Write out to file with .radec extension
                        fname = fname+'.radec'
                        f = open('/tmp/'+fname,'w')
                        f.write(tbl)
                        f.close()
                        time.sleep(0.01)
                        # Send tracktable file to acc
                        f = open('/tmp/'+fname,'r')
                        acc = FTP(self.accini['host'])
                        acc.login()
                        acc.cwd('parm')
                        acc.storlines('STOR '+fname,f)
                        acc.close()
                        f.close()
                        f = urllib2.urlopen('ftp://acc.solar.pvt/parm/'+fname)
                        tbl_echo = ''
                        for line in f.readlines():
                            tbl_echo += line.rstrip()+'\n'
                        if tbl != tbl_echo:
                            print 'Error: Transfer of track table',fname,'failed!'
                            print tbl,'not equal\n',tbl_echo
                    #==== SCAN-START ====
                    elif ctlline.split(' ')[0].upper() == '$SCAN-START':
                        # Do any tasks here that are required to start a new scan
                        # Set scan state to on, and create scan header
                        sf_dict['scan_state'] = 1
                        sys.stdout.write('Started new scan\n')
                        # We need an initial call to set_uvw() in order to set RA, Dec and HA
                        # coordinates in scan header dictionary.
                        srcname = sh_dict['source_id']
                        try:
                            # Generate a datime() object at exactly the next upcoming second (time t+1)
                            dsec = util.datime()
                            dsec.set(dsec.get()  + (1 - dsec.dt.microsecond/1000000.)/86400.)
                            src = self.aa.cat[srcname]        # This causes KeyError if source is not found
                        except KeyError:
                            # The current scan header source ID is not in the source catalog
                            srcname = None
                        set_uvw(self.aa,dsec,srcname)
                        
                        print 'Current RA, Dec, HA:',sh_dict['ra'],sh_dict['dec'],sh_dict['ha']
                        sys.stdout.flush()
                        #scan_header(sh_dict)
                        # Update roach katadc sensor data
                        sys.stdout.write('There are '+str(len(self.roaches))+'active ROACHes\n')
                        sys.stdout.flush()
                        for r in self.roaches:
                            rnum = int(r.roach_ip[5:6]) - 1
                            sys.stdout.write('Reading KATADC for '+r.roach_ip+'...')
                            sys.stdout.flush()
                            if r.fpga:
                                r.get_katadc_dict()
                                if r.msg == 'Success':
                                    sh_dict['katadc'][rnum].update(r.katadc)                            
                                    sys.stdout.write(r.msg+'\n')
                                    sys.stdout.flush()
                                else:
                                    # In case of failure, set to empty dictionary
                                    sh_dict['katadc'][rnum] = {}                            
                                    sys.stdout.write('Failed:'+r.msg+'\n')
                                    sys.stdout.flush()
                                # This fails, for some reason--probably just takes too long
                                #sys.stdout.write('Reading clock...')
                                #sys.stdout.flush()
                                #r.brd_clk = r.fpga.est_brd_clk()
                                #if r.msg == 'Success':
                                #    sh_dict['roach_brd_clk'][rnum].update(r.brd_clk)
                                #    sys.stdout.write(r.msg+'\n')
                                #    sys.stdout.flush()
                                #else:
                                #    # In case of failure, set to empty dictionary
                                #    sh_dict['katadc'][rnum] = {}                            
                                #    sys.stdout.write('Failed:'+r.msg+'\n')
                                #    sh_dict['roach_brd_clk'][rnum] = 0                            
                            else:
                                # In case of no communication, set to empty dictionary
                                sh_dict['katadc'][rnum] = {}
                                sh_dict['roach_brd_clk'][rnum] = 0                            
                                sys.stdout.write('FPGA communication failed\n')
                                sys.stdout.flush()
                        # Add information to scan header and store on ACC
                        scan_header(sh_dict)
                    #==== SCAN-STOP ====
                    elif ctlline.split(' ')[0].upper() == '$SCAN-STOP':
                        sf_dict['scan_state'] = -1
                    #==== WAIT ====
                    elif ctlline.split(' ')[0].upper() == '$WAIT':
                        # Need to wait for given number of seconds, so set self.waitmode to True,
                        # set self.nextctlline to point to next following line, and record duration
                        try:
                            dur = int(ctlline.split(' ')[1])
                        except:
                            print 'Could not interpret duration on $WAIT command--defaulting to 10 s'
                            dur = 10
                        self.waitmode = True
                        self.nextctlline = i+1
                        self.wait = dur
                    #==== KATADC_GET ====
                    elif ctlline.split(' ')[0].upper() == '$KATADC_GET':
                        # Get the standard deviation for each KatADC (assumes frequency tuning is static)
                        # Note, takes about 0.2 s for each active ROACH
                        for r in self.roaches:
                            r.get_attn()
                            rnum = int(r.roach_ip[5:6])
                            if r.msg == 'Success':
                                sdev = dict(zip(['sdev.adc0.h','sdev.adc0.v','sdev.adc1.h','sdev.adc1.v'],r.sdev))
                                sh_dict['katadc'][rnum].update(sdev)
                            else:
                                # In case of failure, set to empty dictionary
                                sh_dict['katadc'][rnum] = {}
                    #==== REWIND ====
                    elif ctlline.split(' ')[0].upper() == '$REWIND':
                        # Get date of first line of current schedule
                        # and increment by 1 day, then autogenerate a
                        # new schedule
                        self.toggle_state()  # Turn off schedule
                        # Get time of first line in schedule and add a day
                        mjd1 = mjd(self.L.get(0))
                        d = util.datime()
                        d.set(mjd1+1)
                        self.Today(d)
                        self.Clear()
                        self.toggle_state()  # Turn schedule back on
                else:
                    cmds = ctlline.split(' ')
                    if cmds[0].upper() == 'FSEQ-FILE':
                        # This is an FSEQ-FILE command, so find and set frequency sequence
                        # Just FTP sequence file from ACC
                        fseqfile = urllib2.urlopen('ftp://acc.solar.pvt/parm/'+cmds[1])
                        nrpt = None     # Initially not defined
                        fsequence = ''  # Initially empty
                        for line in fseqfile.readlines():
                            # Find DWELL line (contains 35 dwell times, in s, one for each
                            # defined band) and split into 35 "repeat" numbers
                            keywd = 'LIST:DWELL'
                            if line.find(keywd) == 0:
                                dwellseq = line[len(keywd):].split(',')
                                if len(dwellseq) != 35:
                                    print 'FSEQ file',cmds[1],'DWELL line must have 35 entries.'
                                    break
                                # Find nearest-integer number of 0.02 s periods
                                nrpt = (numpy.array(dwellseq).astype('float')/0.02 + 0.5).astype('int')
                            keywd = 'LIST:SEQUENCE'
                            if line.find(keywd) == 0:
                                if nrpt is None:
                                    print 'FSEQ file',cmds[1],'DWELL line must come before SEQUENCE line.'
                                    break
                                bands = numpy.array(line[len(keywd):].split(',')).astype('int')
                                # Step through bands in fsequence and repeat them according to
                                # nrpt in order to form a 50-element sequence
                                for band in bands:
                                    for i in range(nrpt[band]):
                                        fsequence += str(band)+','
                                break
                        if fsequence == '':
                            print 'FSEQ file',cmds[1],'not successfully interpreted.'
                        else:
                            sh_dict.update({'fsequence':fsequence[:-1]})   # -1 removes trailing ','
                            self.sequence2roach(fsequence[:-1])

app = App()

mainloop()
