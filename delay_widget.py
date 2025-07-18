#!/usr/bin/env python
'''
   Interactive tool to adjust measured phasecal phases and determine
   delays.  Requires a saved NPZ file containing the phasecal data.
'''
# History:
#   2017-Mar-05  DG
#      First complete version
#   2017-Mar-07  DG
#      Corrected sign error in Y-X delay, and now write to database
#      and ACC without specifying a time (uses current time).  That
#      means the current time should not be too different from the
#      time in the data, but I leave that up to the user.
#   **Breaking News**
#      The new 300 MHz correlator is working!  Delays determined based
#      on data from that correlator have flipped signs!  So I made the
#      change on saving them.
#   2017-Apr-07  DG
#      Cleaned up the "show" button output, and set the default directory
#      for npz files to the phasecal web directory (/common/webplots/phasecal)
#   2017-Jul-08  DG
#      Added new menu item for setting the Ant 14 low-frequency receiver
#      delays.  These delays are written into the Ant 15 slot by the
#      cal_header.dla_update2sql() routine.  Note that delays for other
#      antennas are not updated, even if they are non-zero (a warning is given).
#   2018-Jun-08  DG
#      Relatively important change to display antenna phases with respect to
#      ant 1, except for ant 1 itself, which displays baseline 1-14. This
#      turns out to be a relatively minor change.  Also add a checkbox to
#      allow the user to mark missing antennas.
#   2020-May-02  DG
#      Added some much-needed legends to the plots
#   2023-Oct-07  DG
#      Add refant checkbox and ability to change the reference antenna.  This is
#      for the display only.  The delay center reference remains with Ant 1.
#   2025-May-20  DG
#      Significant changes to work with 16 antennas, with Ant A as Ant 16.  Note
#      that this code is only for updating current delays, so does not have
#      to work for past data.

from Tkinter import *
from tkFileDialog import askopenfile
from tkMessageBox import askyesno
import os
import numpy as np
from util import Time, lobe, common_val_idx
import read_idb as ri
import matplotlib.pylab as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, \
                                              NavigationToolbar2Tk

#============================
class App():
    '''Main application
    '''
    def __init__(self):
        '''Create the user interface and initialize it
        '''
        self.root = Tk()
        #self.root.geometry("+100+0")
        # Set window title
        self.root.wm_title('Delay Widget')
        
        self.menu = Menu(self.root)

        filemenu = Menu(self.menu, tearoff = 0)
        self.menu.add_cascade(label = 'File', menu = filemenu)

        filemenu.add_command(label = 'Open npz File...', command = self.Opennpz)
        filemenu.add_command(label = 'Save Delays to ACC', command = self.Save)
        filemenu.add_command(label = 'Save LoRX Delays to ACC', command = self.SaveLoRX)

        self.root.config(menu = self.menu)

        fmain = Frame(self.root)
        fmain.pack()
        line1 = Frame(fmain)
        line1.pack()
        Label(line1, text='Ant A Y-X Delay', font="Helvetica 14").pack(side=LEFT,expand=1)
        self.dlaAvar = StringVar()
        self.dlaA = Entry(line1, textvariable=self.dlaAvar, font="Helvetica 14", width=5)
        self.dlaAvar.set('0.0')
        self.dlaA.bind("<Return>", self.fetch)
        self.dlaA.pack(side=LEFT)
        dlaAupdown = Frame(line1)
        dlaAupdown.pack(side=LEFT)
        self.dlaAupbtn = Button(dlaAupdown, text=u'\u25B2', command=self.upAdla, borderwidth=0, pady=0)
        self.dlaAupbtn.pack(padx=1,pady=0)
        self.dlaAdnbtn = Button(dlaAupdown, text=u'\u25BC', command=self.downAdla, borderwidth=0, pady=0)
        self.dlaAdnbtn.pack(padx=1,pady=0)

        line2 = Frame(fmain)
        line2.pack()
        fleft = Frame(line2)
        fleft.pack(side=LEFT)
        leftant = Frame(fleft)
        leftant.pack(side=TOP)
        Label(leftant, text='Ant', font="Helvetica 14").pack(side=LEFT)
        self.antvar = StringVar()
        self.ant = Entry(leftant, textvariable=self.antvar, font="Helvetica 14", width=3)
        self.antvar.set('1')
        #self.ant.bind("<Return>", self.fetch)
        self.ant.pack(side=LEFT)
        #antbtns = Frame(leftant)
        #antbtns.pack()
        antupdown = Frame(leftant)
        antupdown.pack()
        self.antupbtn = Button(antupdown, text=u'\u25B2', command=self.up, borderwidth=0, pady=0)
        self.antupbtn.pack(padx=1,pady=0)
        self.antdnbtn = Button(antupdown, text=u'\u25BC', command=self.down, borderwidth=0, pady=0)
        self.antdnbtn.pack(padx=1,pady=0)
        
        leftdla = Frame(fleft)
        leftdla.pack()
        Label(leftdla, text='X Delay', font="Helvetica 14").pack(side=LEFT)
        self.delays = np.zeros(15,np.float)
        self.dlavar = StringVar()
        self.dla = Entry(leftdla, textvariable=self.dlavar, font="Helvetica 14", width=5)
        self.dlavar.set(str(self.delays[0]))
        self.dla.bind("<Return>", self.fetch)
        self.dla.pack(side=LEFT)
        #dlabtns = Frame(leftdla)
        dlaupdown = Frame(leftdla)
        dlaupdown.pack()
        self.dlaupbtn = Button(dlaupdown, text=u'\u25B2', command=self.updla, borderwidth=0, pady=0)
        self.dlaupbtn.pack(padx=1,pady=0)
        self.dladnbtn = Button(dlaupdown, text=u'\u25BC', command=self.downdla, borderwidth=0, pady=0)
        self.dladnbtn.pack(padx=1,pady=0)

        leftxydla = Frame(fleft)
        leftxydla.pack()
        Label(leftxydla, text='Y-X Delay', font="Helvetica 14").pack(side=LEFT)
        self.xydelays = np.zeros(15,np.float)
        self.xydlavar = StringVar()
        self.xydla = Entry(leftxydla, textvariable=self.xydlavar, font="Helvetica 14", width=5)
        self.xydlavar.set(str(self.xydelays[0]))
        self.xydla.bind("<Return>", self.fetch)
        self.xydla.pack(side=LEFT)
        #dlabtns = Frame(leftdla)
        xydlaupdown = Frame(leftxydla)
        xydlaupdown.pack()
        self.xydlaupbtn = Button(xydlaupdown, text=u'\u25B2', command=self.xyupdla, borderwidth=0, pady=0)
        self.xydlaupbtn.pack(padx=1,pady=0)
        self.xydladnbtn = Button(xydlaupdown, text=u'\u25BC', command=self.xydowndla, borderwidth=0, pady=0)
        self.xydladnbtn.pack(padx=1,pady=0)

        var = IntVar()
        leftckbox = Frame(fleft)
        leftckbox.pack()
        self.chkbox = Checkbutton(leftckbox, text="Mark Ant as Missing", variable=var, command=self.cb)
        self.chkbox.var = var
        self.chkbox.pack(side=LEFT)
        self.missing = np.zeros(15,dtype=int) # None of the antennas are missing

        refvar = IntVar()
        leftckbox = Frame(fleft)
        leftckbox.pack()
        self.refchkbox = Checkbutton(leftckbox, text="Mark Ant as Refant", variable=refvar, command=self.refcb)
        self.refchkbox.refvar = refvar
        self.refchkbox.pack(side=LEFT)
        self.refant = np.array([1]+[0]*14)  # Ant 1 is the default refant
        self.refchkbox.select()

        fright = Frame(line2)
        fright.pack()
        
        self.fig = plt.figure(1,(10,5))
        self.ax = []
        loc = [231,232,234,235,133]
        for i in loc:
            self.ax.append(self.fig.add_subplot(i))
            #self.ax[i].grid()
        self.canvas = FigureCanvasTkAgg(self.fig, fright)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(side=TOP, expand=1)
        toolbar1 = NavigationToolbar2Tk(self.canvas, fright)
        toolbar1.update()
        self.canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)

        savebtn = Button(fmain, text='Show', command=self.show)
        savebtn.pack(side=LEFT,padx=3,pady=3)
        quitbtn = Button(fmain, text='Exit', command=self.root.quit)
        quitbtn.pack(side=LEFT,padx=3,pady=3)
        self.label1 = Label(fmain,text='',font='Courier 10')
        self.label1.pack()
        self.label2 = Label(fmain,text='',font='Courier 10')
        self.label2.pack()
        self.label3 = Label(fmain,text='',font='Courier 10')
        self.label3.pack()

        # Some test data to play with
        self.data_source = 'Simulation'
        self.fghz = np.array([3.4, 3.5, 3.6, 4.4, 4.6, 4.8])
        self.pol = np.array([0,2,3,1])
        nf = len(self.fghz)
        taux = np.random.randn(15)
        tauy = taux + np.random.randn(15)*0.1
        tauxA = 0.3
        tauyA = -1.1
        df = np.random.randn(15,4,nf)*0.1
        self.ph = np.zeros((15,4,nf),float)
        for k,fghz in enumerate(self.fghz):
            self.ph[:,0,k] = 2*np.pi*fghz*(tauxA - taux) + df[:,0,k] # XX
            self.ph[:,1,k] = 2*np.pi*fghz*(tauyA - tauy) + df[:,1,k] # YY
            self.ph[:,2,k] = 2*np.pi*fghz*(tauyA - taux) + df[:,2,k] # XY
            self.ph[:,3,k] = 2*np.pi*fghz*(tauxA - tauy) + df[:,3,k] # YX
        self.doplot(ant=1)

    def cb(self):
        # Handle the checkbox widget
        ant_str = self.ant.get()  # Current antenna showing
        ant = int(ant_str)
        self.missing[ant-1] = self.chkbox.var.get()  # List indicating missing ant (if 1)
        
    def refcb(self):
        # Handle the refant checkbox widget
        ant_str = self.ant.get()  # Current antenna showing
        ant = int(ant_str)
        #print 'Ant',ant,'selected.', self.chkbox.var.get(), 'is the state.'
        self.refant = np.zeros(15,np.int)
        self.refant[ant-1] = 1  # List indicating reference antenna (if 1)
        #print 'Refant selection going in is',self.refant
        refant, = np.where(self.refant == 1)
        if len(refant) == 0:
            # If no refant selected, set refant to ant 1
            self.refant = np.array([1]+[0]*14)  # Ant 1 is the default refant
            if ant == 1: self.refchkbox.select()
        #print 'Refant selection going out is',self.refant

    def fetch(self, e):
        ant_str = self.ant.get()
        ant = int(ant_str)
        if e.widget == self.dla:
            dla_str = self.dla.get()
            dla = float(dla_str)
            self.delays[ant-1] = dla
        elif e.widget == self.xydla:
            xydla_str = self.xydla.get()
            xydla = float(xydla_str)
            self.xydelays[ant-1] = xydla
        self.doplot(ant)

    def show(self):
        delays_str = '   '
        xydelays_str = '   '
        delays = np.append(self.delays[0] - self.delays,self.delays[0])
        # Do not change delays for a missing antenna
        bad, = np.where(self.missing == 1)
        delays[bad] = 0.0
        self.label1.configure(text='Ant    1,    2,    3,    4,    5,    6,    7,    8,    9,   10,   11,   12,   13,   14,   15,   16')
        fmt = '{:5.1f},'*16
        delays_str = fmt.format(*delays)
        self.label2.configure(text='   '+delays_str[:-1])
        xydelays = np.append(self.xydelays,-float(self.dlaA.get()))
        xydelays_str = fmt.format(*xydelays)
        self.label3.configure(text='   '+xydelays_str[:-1])

    def up(self):
        ant_str = self.ant.get()
        ant = int(ant_str)
        nsolant = 15
        if ant < nsolant:
            ant += 1
        else:
            ant = nsolant
        self.antvar.set(str(ant))
        # Show current missing status
        if self.missing[ant-1]:
            self.chkbox.select()
        else:
            self.chkbox.deselect()
        if self.refant[ant-1]:
            self.refchkbox.select()
        else:
            self.refchkbox.deselect()
        print 'Up: The current refant selection is',self.refant
        dla = self.delays[ant-1]
        self.dlavar.set(str(dla))
        dla = self.xydelays[ant-1]
        self.xydlavar.set(str(dla))
        self.doplot(ant)

    def down(self):
        ant_str = self.ant.get()
        ant = int(ant_str)
        if ant > 1:
            ant -= 1
        else:
            ant = 1
        self.antvar.set(str(ant))
        # Show current missing status
        if self.missing[ant-1]:
            self.chkbox.select()
        else:
            self.chkbox.deselect()
        if self.refant[ant-1]:
            self.refchkbox.select()
        else:
            self.refchkbox.deselect()
        print 'Down: The current refant selection is',self.refant
        dla = self.delays[ant-1]
        self.dlavar.set(str(dla))
        dla = self.xydelays[ant-1]
        self.xydlavar.set(str(dla))
        self.doplot(ant)

    def upAdla(self):
        ant_str = self.ant.get()
        ant = int(ant_str)
        dla_str = self.dlaA.get()
        dla = float(dla_str)
        dla += 0.1
        self.dlaAvar.set(str(dla))
        self.doplot(ant)

    def downAdla(self):
        ant_str = self.ant.get()
        ant = int(ant_str)
        dla_str = self.dlaA.get()
        dla = float(dla_str)
        dla -= 0.1
        self.dlaAvar.set(str(dla))
        self.doplot(ant)

    def updla(self):
        ant_str = self.ant.get()
        ant = int(ant_str)
        dla_str = self.dla.get()
        dla = float(dla_str)
        dla += 0.1
        self.delays[ant-1] = dla
        self.dlavar.set(str(dla))
        self.doplot(ant)

    def downdla(self):
        ant_str = self.ant.get()
        ant = int(ant_str)
        dla_str = self.dla.get()
        dla = float(dla_str)
        dla -= 0.1
        self.delays[ant-1] = dla
        self.dlavar.set(str(dla))
        self.doplot(ant)

    def xyupdla(self):
        ant_str = self.ant.get()
        ant = int(ant_str)
        dla_str = self.xydla.get()
        dla = float(dla_str)
        dla += 0.1
        self.xydelays[ant-1] = dla
        self.xydlavar.set(str(dla))
        self.doplot(ant)

    def xydowndla(self):
        ant_str = self.ant.get()
        ant = int(ant_str)
        dla_str = self.xydla.get()
        dla = float(dla_str)
        dla -= 0.1
        self.xydelays[ant-1] = dla
        self.xydlavar.set(str(dla))
        self.doplot(ant)

    def doplot(self,ant=1):
        polstr = ['XX','XY','YX','YY']
        dla = self.delays[ant-1]
        ydla = self.xydelays[ant-1]
        # Also need delay settings for ant 1
        dla1 = self.delays[0]
        ydla1 = self.xydelays[0]
        ydlaA = np.float(self.dlaA.get())
        for i,ax in enumerate(self.ax[:4]):
            if self.pol[i] == 0:
                # XX => use only the ant X delay 
                tau = dla
                tau1 = dla1
            elif self.pol[i] == 1:
                # YY => use the ant (X delay + Y-X delay) + Ant A Y-X delay
                tau = dla + ydla + ydlaA
                tau1 = dla1 + ydla1 + ydlaA
            elif self.pol[i] == 2:
                # XY => use the ant X delay + Ant A Y-X delay
                tau = dla + ydlaA
                tau1 = dla1 + ydlaA
            else:
                #YX => use the ant (X delay + Y-X delay)
                tau = dla + ydla
                tau1 = dla1 + ydla1
            ax.cla()
            refant, = np.where(self.refant==1)
            refant = refant[0]
            if ant == refant+1:
                ax.plot(self.fghz,lobe(self.ph[ant-1,self.pol[i]] - 2*np.pi*self.fghz*tau),'.',label=polstr[i])
                if self.pol[i] == 0:
                    pxx = lobe(self.ph[ant-1,self.pol[i]] - 2*np.pi*self.fghz*tau)
                if self.pol[i] == 1:
                    pyy = lobe(self.ph[ant-1,self.pol[i]] - 2*np.pi*self.fghz*tau)
                if self.pol[i] == 2:
                    pxy = lobe(self.ph[ant-1,self.pol[i]] - 2*np.pi*self.fghz*tau)
                if self.pol[i] == 3:
                    pyx = lobe(self.ph[ant-1,self.pol[i]] - 2*np.pi*self.fghz*tau)
            else:
                ax.plot(self.fghz,lobe(self.ph[ant-1,self.pol[i]] - self.ph[refant,self.pol[i]] - 2*np.pi*self.fghz*(tau-tau1)),'.',label=polstr[i])
                if self.pol[i] == 0:
                    pxx = lobe(self.ph[ant-1,self.pol[i]] - 2*np.pi*self.fghz*tau)
                if self.pol[i] == 1:
                    pyy = lobe(self.ph[ant-1,self.pol[i]] - 2*np.pi*self.fghz*tau)
                if self.pol[i] == 2:
                    pxy = lobe(self.ph[ant-1,self.pol[i]] - 2*np.pi*self.fghz*tau)
                if self.pol[i] == 3:
                    pyx = lobe(self.ph[ant-1,self.pol[i]] - 2*np.pi*self.fghz*tau)
            ax.set_ylim(-4,4)
            ax.legend(fontsize=9,loc='lower right')
        self.ax[4].cla()
        self.ax[4].plot(self.fghz,lobe(pyy - pxx),'.')
        self.ax[4].plot(self.fghz,lobe(pyx - pxy),'.')
        self.ax[4].set_title('YY - XX Phase')
        #self.ax[4].plot(self.fghz,lobe(pxy - pyx),'.')
        self.ax[4].set_ylim(-4,4)
        self.canvas.draw()

    def Opennpz(self):
        init_dir = '/common/webplots/phasecal/'
        f = askopenfile(initialdir = init_dir, mode = 'r',
                        filetypes = [('NPZ files','*.npz'),('all files','*.*')])
        if f is None:
            # User cancelled, so do nothing.
            return
        data = np.load(f, allow_pickle=True)
        self.root.wm_title('Delay Widget '+f.name)
        k = data.keys()
        out = data[k[0]].item()
        self.fghz = out['fghz']
        # Set time for delays as end time of data (for lack of a better choice)
        self.time = Time(out['time'][-1],format='jd')
        nsolant = 15
        self.ph = np.angle(np.sum(out['x'][ri.bl2ord[0:nsolant,nsolant]],3))
        self.data_source = 'Data'
        self.doplot()
                
    def Save(self):
        # Send saved delays to the SQL database and the ACC (if the data source is not 'Simulation')
        if self.data_source == 'Data':
            if (Time.now().mjd - self.time.mjd) > 1.0:
                question = "Warning: Data more than a day old. Are you sure you want to save delays to SQL and ACC?"
            else:
                question = "Save delays to SQL and ACC?"
            import cal_header as ch
            # Calculate delays relative to Ant 1 and tack Ant1-A delay at end
            delays = self.delays[0] - self.delays
            delays = np.append(delays,self.delays[0])
            # Have to change the sign of Ant A Y-X delay, hence the minus sign
            xydelays = np.append(self.xydelays,-float(self.dlaA.get()))
            # Do not change delays for a missing antenna
            bad, = np.where(self.missing == 1)
            delays[bad] = 0.0
            if askyesno("Write Delays",question):
                # All Y-X delays need a sign flip, hence the minus sign
                # ************ This block commented out due to loss of SQL **************
                ch.dla_update2sql(delays,-xydelays)
                #ch.dla_update2sql(-delays,xydelays)  # 300 MHz design uses flipped signs!
                ch.dla_censql2table()
                # Replaced by
                #ch.dla_update2table(delays,-xydelays)


    def SaveLoRX(self):
        # Send saved Ant A delays as Ant 15 to the SQL database and the ACC (if the data source is not 'Simulation')
        # Note, ONLY Ant A delay is changed, and it is ascribed to Ant 15.  No other delay changes are made.
        if self.data_source == 'Data':
            if (Time.now().mjd - self.time.mjd) > 1.0:
                question = "Warning: Data more than a day old. Are you sure you want to save delays to SQL and ACC?"
            else:
                question = "Save delays to SQL and ACC?"
            import cal_header as ch
            
            delays = self.delays[0] - self.delays
            delays = np.append(delays,self.delays[0])
            # Have to change the sign of Ant A Y-X delay, hence the minus sign
            xydelays = np.append(self.xydelays,-float(self.dlaA.get()))
            # Do not change delays for a missing antenna
            bad, = np.where(self.missing == 1)
            delays[bad] = 0.0
            # Check that the delays to all antennas except Ant A are zero, or give warning if not
            for i in range(15):
                if delays[i] != 0.0:
                    question = "Some Ant1-15 delays are not 0, but will NOT be updated.  Save anyway?"
                    break
                if xydelays[i] != 0.0:
                    question = "Some Ant1-15 delays are not 0, but will NOT be updated.  Save anyway?"
                    break
            
            if askyesno("Write Delays",question):
                # All Y-X delays need a sign flip, hence the minus sign
                # ************ This block commented out due to loss of SQL **************
                ch.dla_update2sql(delays,-xydelays,lorx=True)
                #ch.dla_update2sql(-delays,xydelays,lorx=True)  # 300 MHz design uses flipped signs!
                ch.dla_censql2table()
                # Replaced by
                #ch.dla_update2table(delays,-xydelays,lorx=True)
                
app = App()

mainloop()
        
