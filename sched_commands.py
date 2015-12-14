#!/usr/bin/env python
'''
          Minimal widget for sending commands to ACC.
'''
#
# History:
#   2015-Jun-11  RG
#     First Written
#   2015-Nov-29  DG
#     Added $LNA-INIT command, for convenience
#   2015-Dec-01  DG
#     Changed default width of text entry box.
#

import os, time
import urllib2
from Tkinter import *
import stateframe
import socket

class App():
    '''Main application
    '''
    def __init__(self):
        '''Create the user interface and initialize it
        '''
        self.root = Tk()
        self.cmd_list = []
        self.ptr = -1
        global title
        title = 'Command box'
        self.root.wm_title(title)
        fmain = Frame(self.root)
        fmain.pack()

#        fmain = Frame()
#        self.nb.add(fmain,text='Main')

        #self.B3 = Button(text="Raw Command", command=self.send_cmd)
        #self.B3.pack( side = RIGHT)
        #command = StringVar()
        #self.E3 = Entry(textvariable = command)
        rawlabel = Label(fmain,text='Raw Command:')
        rawlabel.pack(side=LEFT)
        E3 = Entry(fmain, width=35)
        E3.pack( side = LEFT, expand=True, fill=X)
        E3.bind('<Return>',self.send_cmd)
        E3.bind('<Up>',self.up)
        E3.bind('<Down>',self.down)
        self.accini = stateframe.rd_ACCfile()

    def send_cmd(self,event):
        '''Send a raw command.'''
        command = event.widget.get()
        if command != '':
            if command.split(' ')[0].upper() == '$LNA-INIT':
                # Get LNA_settings.txt file from ACC and send the series of
                # commands needed to set the LNA voltages
                userpass = 'admin:observer@'
                lnafile = urllib2.urlopen('ftp://'+userpass+'acc.solar.pvt/parm/LNA_settings.txt',timeout=1)
                lines = lnafile.readlines()
                lnafile.close()
                lnas = {0:'hh',1:'lh',2:'lv',3:'hv'}
                lnas_a = [{},{},{},{}]
                lnas_b = [{},{},{},{}]
                try:
                    for i,line in enumerate(lines):
                        if line.find('[FEMA]') == 0:
                            # Found FEMA section, so take lines i+2 through i+5 as data lines
                            for k in range(4):
                                lna,fstr,polstr,model,sn,vdrain,vg1,vg2,idrain =lines[i+2+k].split()
                                lnas_a[int(lna)] = {'vd':float(vdrain),'vg1':float(vg1),'vg2':float(vg2)}
                        if line.find('[FEMB]') == 0:
                            # Found FEMB section, so take lines i+2 through i+5 as data lines
                            for k in range(4):
                                lna,fstr,polstr,model,sn,vdrain,vg1,vg2,idrain =lines[i+2+k].split()
                                lnas_b[int(lna)] = {'vd':float(vdrain),'vg1':float(vg1),'vg2':float(vg2)}
                except:
                    print 'Error reading/parsing LNA_settings.txt file from ACC'
                        
                try:
                    for i in range(4):
                        cmdstr = 'LNA-ENABLE '+lnas[i]+' on ANT14'
                        self.execute_ctlline(cmdstr)
                        cmdstr = 'LNA-DRAIN '+lnas[i]+' '+str(lnas_a[i]['vd'])+' ANT14'
                        self.execute_ctlline(cmdstr)
                        cmdstr = 'LNA-GATE1 '+lnas[i]+' '+str(lnas_a[i]['vg1'])+' ANT14'
                        self.execute_ctlline(cmdstr)
                        cmdstr = 'LNA-GATE2 '+lnas[i]+' '+str(lnas_a[i]['vg2'])+' ANT14'
                        self.execute_ctlline(cmdstr)
                        cmdstr = 'LNA-ENABLE '+lnas[i]+' on ANT15'
                        self.execute_ctlline(cmdstr)
                        cmdstr = 'LNA-DRAIN '+lnas[i]+' '+str(lnas_b[i]['vd'])+' ANT15'
                        self.execute_ctlline(cmdstr)
                        cmdstr = 'LNA-GATE1 '+lnas[i]+' '+str(lnas_b[i]['vg1'])+' ANT15'
                        self.execute_ctlline(cmdstr)
                        cmdstr = 'LNA-GATE2 '+lnas[i]+' '+str(lnas_b[i]['vg2'])+' ANT15'
                        self.execute_ctlline(cmdstr)
                except:
                    print 'Error sending LNA_settings to ACC'        
        
        
            self.execute_ctlline(command)
        self.cmdlist(command) 
        event.widget.delete(0,END)             

    def up(self,event):
        if len(self.cmd_list) <= -self.ptr:
            pass
        else:
            self.ptr -= 1
            self.use_ptr(event)
 
    def down(self,event):
        if self.ptr > -1:
            pass
        else:
            if self.ptr == -1:
                event.widget.delete(0,END)
                self.ptr = 0
            else:
                self.ptr += 1
                self.use_ptr(event)
        
    def use_ptr(self,event):
        event.widget.delete(0,END) 
        event.widget.insert(0,self.cmd_list[self.ptr])
        
    def execute_ctlline(self,ctlline,mjd1=None,mjd2=None):
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
                
    def cmdlist(self,command):
        if self.cmd_list == [] and command != '':
            self.cmd_list.append(command)
        elif command == '' or self.cmd_list[-1] == command:
            pass
        else:
            if len(self.cmd_list) == 20:
                self.cmd_list.pop(0)
            self.cmd_list.append(command)
        self.ptr = 0
        
app = App()

mainloop()
