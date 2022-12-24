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
# History:
#   2014-Nov-15  DG
#      Started this history log.  Added $TRIPS command and added antenna
#      diagnostics update every 5 minutes.  Extended height of window
#      (controlled by Macro list box) to 15 lines.
#   2014-Nov-17  DG
#      Changed antenna diagnostics update to occur at half-minute mark
#      of 5 minute updates, to avoid collision with other commands.
#   2014-Nov-28  DG
#      Fixed the default ROACH channel assignments 'Antlist' in the scan_header
#      to be [1,2,3,4,0,0,0,0,0,0,0,0,0,0,0,0].  Also, no longer flag antennas that
#      are not in 'Antlist' as not tracking.
#   2014-Nov-29  DG
#      Added information to sh_dict['project'] identifying the type of
#      observation, mainly to allow finding of SOLPNTCAL right now, but
#      likely to be of general use.  Also improved visibility of errors
#      when writing to the SQL Server.  Had to put 'Antlist' back to (wrong)
#      original value due to bug in Fortran software--for now...
#   2014-Dec-07  DG
#      Add default project/source ID in cases where none is specifically
#      defined/hard-coded.  Note that default Track Mode is 'FIXED' in this
#      case.  Also tried to get $SCAN-START to work as a raw command.
#   2014-Dec-08  DG
#      Finally fixed problem with scan_state getting set to -1. It was
#      happening in get_uvw() when srcname was None.  It should be possible
#      (e.g. for total power scans) to take data with no source name... Hmm,
#      still not working.  Also changed default observation Project from
#      "Normal Observing" to "NormalObserving".
#   2015-Jan-19  DG
#      Attempt to run 4 ROACH boards (2 separate systems)
#      Updated: antlist, roach_ips
#   2015-Feb-11  DG
#      Implement $DLASWEEP command, which sweeps delay one per second
#      on given antenna from a start value to a stop value relative to
#      current delay.
#   2015-Feb-14  DG
#      Rather major change to interpret a "Raw Command" just like any
#      other atomic command in a .ctl file.  This should make all
#      legitimate atomic commands, even those interpreted by the schedule
#      (that start with $) runnable as a "Raw Command." 
#   2015-Feb-15  DG
#      Several changes to get the code working for Geosats again.
#   2015-Mar-02  DG
#      Change the 5000 step delay offset to 5000.5 to make astype(int)
#      act as a round().
#   2015-Mar-29  DG
#      On starting a new scan, reads DCM_master_table.txt and creates
#      dcm.txt according to the frequency sequence, then transfers it to
#      the ACC.
#   2015-Apr-02  JV
#      Modifications to enable running a second subarray: use commands
#      > python schedule.py Subarray2
#      for OVSA, or
#      > python schedule.py Starburst
#      for Starburst.  In either case, a master schedule (Subarray1), initiated
#      by running schedule.py without args, must be running, because only the
#      master schedule updates the antenna diagnostics in the ACC stateframe,
#      and only the master schedule writes stateframe data to the SQL database.
#      The second subarray schedule writes stateframe data to the ACC using a
#      different port than the master schedule.
#   2015-Apr-03  JV
#      - Modified App.get_subarray_pid() (used to check if a schedule with same
#        name is already running) to be case insensitive.
#      - Added a '$SUBARRAY' command to execute_ctlline() which runs SUBARRAY1
#        when run by schedule 1 and SUBARRAY2 when run by schedule 2.  $SUBARRAY
#        command can be invoked in two ways (shown using examples):
#         $SUBARRAY ant1-8 ant15
#         $SUBARRAY default.antlist starburst
#        where the 3rd argument in the 2nd example is the name of an antlist specified
#        in the file default.antlist.
#      - Wrote function get_antlist(antlistname,antlistfile) to read antlistfile and
#        return the antlist defined for name antlistname.  This way we can modify
#        default.antlist (or other antlistfile) instead of modifying all the .ctl files
#        when we want to move ants in and out of the array.
#   2015-Apr-03  DG
#      Add s = socket.socket(socket.AF_INET, socket.SOCK_STREAM) ahead of s.connect()
#      in places where it was missing.
#   2015-May-29  DG
#      Converted from using datime() to using Time() based on astropy.
#   2015-Jun-12  DG
#      Added command history code written by Rob Gelosa
#   2015-Jun-16  DG
#      FTP to ACC now requires a username and password
#   2015-Jun-19  DG
#      Now gets ROACH antenna assignments from eovsa_corr.ini, so that this
#      can be changed in a single place and propagate correctly to the rest
#      of the system.
#   2015-Jun-25  DG
#      Now that Ant13's solar power station is online, changed rd_solpwr() to 
#      read from either power station depending on supplied url.
#   2015-Jun-26  DG
#      Finally found and fixed the bug that was preventing automatic update of 
#      choice of calibrators.
#   2015-Jul-07  DG
#      Added subbw for 400 MHz to scan_header (not used in DPP?)
#   2015-Jul-25  DG
#      Changes to allow for different X and Y delay centers in file 
#      acc:/parm/delay_centers.txt.  Added dlaceny to sh_dict.
#   2015-Jul-27  DG
#      Added optional polarization to $DLASWEEP command (X, Y, or omitted=>both)
#   2015-Jul-29  DG
#      I found that the DCMTABLE command was not being issued to send DCM.TXT
#      to the ACC.  This is now automatically done just before any DCMAUTO-ON
#      command is sent.
#   2015-Aug-09  DG
#      Added provision to sweep all delays in $DLASWEEP command if ant = 0 is 
#      specified.  This is useful for total power polarization measurements, 
#      if only X or Y is swept.
#   2015-Aug-17  DG
#      Changed to allow setting different adc clock frequency in one place
#      [in connect2roach(), i.e. self.brd_clk_freq]
#   2015-Oct-13  DG
#      Added $LNA-INIT command to read ACC file LNA_settings.txt and send the
#      corresponding commands to the ACC.
#   2015-Oct-22  DG
#      Finally had a chance to debug $LNA-INIT command.  I ended up having to
#      create a new procedure sendctlline() that creates the socket, sends, and
#      closes the socket for each line.  Then it was possible to send the
#      multiple lines needed for $LNA-INIT.
#   2015-Oct-24  DG
#      Added PLANET macro command to track a PLANET other than the Sun.
#   2015-Oct-27  DG
#      Expand to 6 ROACHes (three pairs) [change to antlist, and to list of roach IPs]
#      --one ROACH (#6) seems messed up, so going back to 4 ROACHes for now.
#   2015-Nov-29  DG
#      Updated $LNA-INIT command to use new command names and syntax.
#   2016-Jan-16  DG
#      Added code for $PCYCLE command, to power-cycle a device in the field
#      (antenna, crio, or fem) using the new Viking relay controllers
#   2016-Jan-19  DG
#      Added attempt to read pwr_cycle queue
#   2016-Feb-27  DG
#      This time I really did add sending of DCM.TXT before DCMAUTO-ON
#   2016-Feb-28  DG
#      Changed printing of elapsed time to print only if more than 10 ms
#      early or late (should mean shorter log file...).  Also fixed a not 
#      very important bug in dlasweep.
#   2016-Feb-29  DG
#      Added $SCAN-RESTART command, to just turn on the scan state
#   2016-Mar-08  DG
#      Changed delay offset from 5000 to 8000, to reflect new range of
#      16-ant correlator (0-16000 steps).
#   2016-Mar-15  JRV
#      Added function check_27m_sun, which is called once per second (by
#      inc_time) to make sure 27-m's do not go w/in 10 degrees of Sun 
#   2016-Mar-17  DG
#      Changed delay offset from 8000 to 9000, to reflect new range of
#      16-ant correlator (0-16000 steps).
#   2016-Mar-30  DG
#      I discovered that I was writing the wrong polarizations.
#   2016-May-04  DG
#      Implemented writing of DCM table to sql server as well as to ACC.
#      Also added $FEM-INIT command to reset the FEM attenuations to
#      their optimal value (for power level 3 dBm)
#   2016-May-20  DG
#      Starting a new scan now takes much longer (not sure why), so changed
#      the wake_up() timer from 10 s to 15 s
#   2016-Aug-17  DG
#      Change to connect2roach() to minimize number of FTP accesses to 
#      files on ACC.
#   2016-Sep-07  DG
#      Aborted the preceding, since it never worked, and the ACC is fixed...
#   2016-Oct-26  DG
#      Added $CAPTURE-1S command, to capture 1 s of data packets on dpp.
#      Data recording must first be stopped, using $SCAN-STOP.
#   2016-Nov-04  DG
#      Add $SCAN-START NODATA option, so that a scan is set up for running
#      but no data is recorded (needed prior to $CAPTURE-1S).
#   2016-Nov-15  DG
#      Implement two new commands, $PA-TRACK and $PA-STOP, to initiate and
#      abort automatic tracking of the position angle setting of the 27-m
#      focus rotation mechanism.  $PA-TRACK requires an antenna number, 
#      which is the antenna from which to read the value of parallactic angle.
#   2016-Nov-22  DG
#      Changed the code to read delay centers from the database at the
#      start of each scan (so code does not have to be reloaded every time).
#      Note that it reads from the SQL database, not the ACC file
#      /parm/delay_centers.txt, so proper procedures are needed to
#      ensure that these are the same!  The dppxmp program uses that file.
#   2016-Nov-25  DG
#      Delays on ant 13 are still going negative.  Change Ant1 delay to
#      11000, to counteract it.  I hope this does not drive other delays
#      past 16000
#   2016-Nov-26  DG
#      Added CALPNTCAL command handling, to do 27-m pointing measurement on
#      calibrators.
#   2016-Dec-10  DG
#      Added $PA-SWEEP command, which rotates the FRM from -PA to PA at
#      a given rate.  Uses the same PA_thread variable so that $PA-STOP
#      works to stop it, and it cannot run if $PA-TRACK is running (and
#      vice versa).  Also fixed a number of problems where I was referring
#      to np instead of numpy--there is no np in this module!
#   2016-Dec-12  DG
#      Changed autogen() to expect the first line to be an ACQUIRE line,
#      second to be a PHASECAL, and third to be SUN.
#   2016-Dec-18  DG
#      Discovered that there is no difference between make_tracktable and 
#      make_geosattable, so I eliminated the latter.  Also changed split(' ')
#      to just split() everywhere, to remove pointless requirement of only
#      one space between elements of schedule commands.
#   2017-Jan-05  DG
#      It turns out that there were some subtle but important differences
#      in make_geosattable(), so I reinstated it.
#   2017-Jan-06  DG
#      Discovered that X and Y delays were being set independently, which
#      made the difference in delay between X and Y change by 1 step.  This
#      unwanted behavior explains the random changes in X vs. Y (and XX vs. YY)
#      delay that we see in the data.  Now only X controls when delays are
#      changed, and X and Y are changed together.
#   2017-Feb-06  DG
#      Major change to set all source coordinates as J2000 TOPOCENTRIC.  This
#      means setting the EPOCH string in the scan_header to '2000', and converting
#      RA and Dec coordinates from TOPOCENTRIC of DATE as follows:
#         ra_j2000 = (src.ra - src.g_ra) + src.a_ra
#         dec_j2000 = (src.dec - src.g_dec) + src.a_dec
#      where src.ra and src.dec are the old topocentric coordinates of date,
#      src.g_ra, src.g_dec are the geocentric coordinates of date, and
#      src.a_ra, src.a_dec are the astrometric (J2000, geocentric) coordinates.
#   2017-Feb-08  DG
#      The above strategy did not work.  I determined that the call to check the 27-m 
#      position was calling aa.set_jultime(), which has a side-effect of setting the 
#      epoch to date.  I replaced that with a change to aa.date.  I also verified
#      that uvw coordinates are calculated correctly with the current scheme.  Still
#      not good enough.  The epoch is getting reset somewhere else.  I finally gave
#      up and just set the epoch to J2000 once per second!
#    2017-Mar-05  DG
#      The 300 MHz correlator is working!  However, that led me to discover that
#      a 200 MHz board clock speed was hard-coded here.  I have changed it to
#      check the ADC clock speed in the ROACH config file (attached to the
#      roach objects when connecting to the ROACHes), and set it to 1/4 of that.
#      If for some reason that fails, the clock speed defaults to 300 MHz.  I
#      also increased the delay offset to 16000.  Also added reading of new ACC
#      file with the ROACH sync time, in init_scanheader_dict(), which replaces 
#      the erroneous value we have been using in the scan_header.  I also changed
#      time reference time of the uvw to correspond to the current second.
#    2017-Mar-06 DG
#      Changed delay offset to scale with adc_clk frequency, so that it works for
#      either 200 or 300 MHz design.
#    2017-Apr-22  DG
#      Updated $CAPTURE-1S handling to allow it to work when no <stem> argument is given.
#    2017-May-18  DG
#      Commented out lines relating STARBURST
#    2017-Jun-28  DG
#      Added "crossed" keyword handling for $PA_ADJUST command
#    2017-Jul-08  DG
#      Changes to allow detection of Ant 14 receiver position (based on stateframe
#      information from FEMA, and remembered via self.lorx), and setting of Ant 14
#      delays based on that, to those in slot for Ant 15.
#    2017-Jul-10  DG
#      After some confusing results, I realized that the ACC delay_centers.txt
#      file also had to be changed, because the DPP is using that.  Therefore,
#      the whole scheme was updated to read from SQL at a $SCAN-START, do the 
#      swap at that point, if necessary, and then create and FTP the file. This
#      solves another problem that the ACC file could in principle deviate from
#      the SQL table.  Now it cannot (if all goes smoothly).
#    2017-Sep-01  DG
#      Added update_status() to create a status file with the currently running
#      schedule.
#    2018-Apr-09  DG
#      My earlier change of split(' ') to just split() everywhere seems to have
#      vanished, so I reedited this change.
#    2018-Jun-08  DG
#      Added a distinct PHASECAL_LO command, in order to selectively change the
#      receiver commands for a scan using the low-frequency receiver.  At the
#      moment, a different set of DCM attenuations are needed for that receiver,
#      but this ability is probably a good idea in general.
#    2018-Sep-18  DG
#      Added new functionality for New menu item!  Automatically make today's
#      solar schedule.
#    2018-Oct-24  DG
#      Rearrange code to check if another schedule is running BEFORE opening a
#      new schedule.log file.  This avoids overwriting the schedule.log file of
#      a running schedule.
#    2019-Feb-22  DG
#      Import chan_util from new chan_util_52, which defines things for new
#      IF filters, e.g. 52 channels of 325 MHz bandwidth.
#    2019-Nov-24  DG
#      Added $WSCRAM-LIMIT command to set the default windscram limit for 27-m and
#      code to automatically set the windscram limit to 0 if the weather station
#      information is "stale" (older than 5 minutes).  One side effect is that the
#      windscram limit will be set to the default 17 mph on restarting the schedule.
#    2020-Mar-08  DG
#      Added $PLUSDELAY command to add a given delay to the current ANT 14 LO-FRQ 
#      RCVR delays (both X and Y).  This does write a new delay-center record to
#      the SQL database, so use cal_header.delete_cal() to remove it if desired.
#    2020-Aug-09  DG
#      Changed the way the low-frequency receiver delays for Ant 14 get set.  It
#      was just checking for the receiver position, but that failed whenever there
#      was a glitch in the stateframe at the moment of a $SCAN-START command.
#      Now it sets self.lorx to True whenever there is a command RX-SELECT LO.
#    2021-Jan-31  DG
#      Rename existing /tmp/schedule.log file to /tmp/schedule_old.log for debug
#      purposes.
#    2021-May-09  DG
#      Major change to DCM.txt.  Today I learned that the DCM attenuations are
#      switching correctly, but with a lag of 7 slots (140 ms), which means that
#      the desired attenuation settings were not getting applied.  This
#      rather badly messed up the ADC levels, but it was never possible to know
#      this until recently when we created a means to synchronously measure
#      the ADC levels.  A change was made to sequence2dcmtable() to simply
#      rotate the list of attenuation settings used by the system (file 
#      DCM_table.txt) by 7 lines to compensate.  If the lag amount ever changes,
#      or the code sending the DCM attenuations is fixed, simply change the
#      value of "lag" accordingly (0 for no lag).
#    2021-Aug-11  DG
#      After years of occasional schedule hangs, I now realize that the wake_up()
#      call handler requires two other arguments!  Hopefully this will not crash
#      anymore.
#    2021-Aug-18  DG
#      Fix a bug related to non-existent /tmp/schedule_old.log on Helios reboot.
#      Also print what antennas are associated with each ROACH in connect2roach().
#      Also change old log file name to /tmp/schedule<yyyymmdd_hhmm>.log so that
#      it does not get overwritten
#    2021-Aug-20  DG
#      Avoid bug in execute_ctlline() where an empty line in a .ctl file caused a crash.
#    2021-Sep-29  DG
#      A wake_up() call did not restart the schedule, so now I have added code in 
#      wake_up() to cancel the inc_time timer and restart it.
#    2022-Jan-12  DG
#      Added a FLARE* project ID, to trigger the new fast recording mode.
#    2022-Mar-07  DG
#      Multiple changes due to loss of SQL.  Blocks changes are preceded by
#            # ************ This block commented out due to loss of SQL **************
#    2022-Mar-08  DG
#      Added log_stateframe() routine to make logging to file work.
#    2022-Mar-14  DG
#      Changes to use new chan_info_52.py code to define a Chan_Info object.  The main
#      purpose is to enable a fast FLARE mode by specifying a DWELL mode for
#      one band specified in a dwellXX.fsq file.
#    2022-Apr-02  DG
#      Restored the SQL code (for now, but still log the scanheader and stateframe
#      just in case)
#    2022-Apr-10  DG
#      Fixed bug in writing scanheader to SQL.
#    2022-May-14  DG
#      Added a "No 27m" checkbox, which makes the Today button create a solar
#      schedule that skips all 27-m calibrations.  The 27-m will move for the SKYCALTEST
#      if it is in the subarray (remove it by editing the default.antlist file).
#    2022-Dec-11  DG
#      Explicitly set implicit_transactions to OFF, since it seems to be getting turned
#      on somehow, resulting in a spurious error on writing stateframe to SQL.
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
import threading, pwr_cycle
import subprocess
import roach
from eovsa_tracktable import *
from eovsa_array import *
from eovsa_lst import eovsa_ha
from math import pi
from readvla import *
import chan_info_52 as ci
from scan_header import scan_header
from gen_schedule_sf import *
import stateframe, stateframedef
from aipy.phs import PointingError
import corr, time, numpy, socket, struct, sys
import ephem
import eovsa_cat
from eovsa_visibility import scan_visible
from whenup import whenup
#import starburst
from matplotlib.mlab import find
import cal_header
import adc_cal2
import pcapture2
from whenup import make_sched, remove_cal

# Determine whether this is the master schedule (Subarray1) or controlling a second subarray
# To run the master schedule, just type > python schedule.py
#                                    or > python schedule.py Subarray1
# To control a second subarray,    type > python schedule.py Subarray2
#                                    or > python schedule.py Starburst
if len(sys.argv) < 2: # master schedule (Subarray1)
    subarray_name = 'Subarray1'
else:
    subarray_name = sys.argv[1] # this option is for a 2nd OVSA subarray - Subarray2 is suggested name,
                                # but any name should work other than Starburst or Subarray1

mypid = str(os.getpid())

#============================
def get_subarray_pid(subarray_name):
    # check whether a subarray with subarray_name ('Subarray1' for master, or 'Subarray2' or 'Starburst')
    # is running in another instance of schedule.py (case-insensitive).
    # return PID if it is running, -1 if it is not
    pidlist = subprocess.check_output(["pidof","python"]).split() # list of PIDs for all python processes
    for pid in pidlist:
        if mypid != pid:
            ps_out = subprocess.check_output(["ps","-lfp",pid])
            ind = ps_out.find('schedule.py') # position in string at which 'schedule.py' can be found (-1 if not found)
            if ind != -1:  # if this PID is running schedule.py
                if subarray_name=='Subarray1': # if we are checking if the master schedule is running (Subarray1)
                    if len(ps_out)==ind+len('schedule.py\n'): # check if this PID is for a master schedule (no extra args after schedule.py)
                        return pid
                elif ps_out.upper().find(subarray_name.upper()) != -1: # if the name of the subarray we are checking was used for this PID
                    return pid
    return -1

# Check whether a schedule is already running.  Only one of each subarray is allowed.
match_pid = get_subarray_pid(subarray_name)
if match_pid != -1:
    showerror('Error','Another '+subarray_name+' schedule is already running (PID '+ match_pid + ')')
    exit()

# Ensure that output to "terminal" goes to log file.
if len(sys.argv)<2: # for master schedule write to schedule.log
    datstr = util.Time.now().iso[:16].replace('-','').replace(':','').replace(' ','_')
    try:
        # Rename existing log file for debug purposes
        os.rename('/tmp/schedule.log','/tmp/schedule'+datstr+'.log')
    except:
        print 'Could not rename /tmp/schedule.log to /tmp/schedule'+datstr+'.log.  Perhaps a Helios reboot?'
    sys.stdout = open('/tmp/schedule.log','w')
else: # use a different log file name for the 2nd subarray (schedule_Subarray2.log or schedule_Starburst.log)
    sys.stdout = open('/tmp/schedule_'+sys.argv[1]+'.log','w')

# Show the scanheader dictionary as undefined
global sh_dict, sf_dict
userpass = 'admin:observer@'
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
def init_scanheader_dict(version=37.0):
    ''' Create the initial scan header dictionary, with reasonable defaults.
        Entries will be overridden before creating the scan_header file using
        scan_header.py
    '''
    global sh_dict

    #userpass = 'admin:observer@'
    t = util.Time.now()
    mjd_ = t.mjd      # Get mjd of Time() object
    timestamp = t.lv  # Get LabVIEW timestamp

    aa = eovsa_array()
    aa.date = str(aa.date)[:11]+'00:00'  # Set date to 0 UT
    #print t.iso,aa.epoch
    mjd0 = aa.date + 15019.5   # Convert from ephem date to mjd

    try:
        f = urllib2.urlopen('ftp://'+userpass+'acc.solar.pvt/parm/acc0time.txt',timeout=1)
        mjdacc0 = np.double(f.readline().split()[0])
        f.close()
    except:
        print t.iso,'ACC connection for acc0 time timed out.  Reading from /tmp/acc0time.txt'
        f = open('/tmp/acc0time.txt','r')
        mjdacc0 = np.double(f.readline().split()[0])
        f.close()
      
    # ************ This block commented out due to loss of SQL **************
    try:
        xml, buf = cal_header.read_cal(4)
        dcenters = stateframe.extract(buf,xml['Delaycen_ns'])
        dcen  = dcenters[:,0]
        dceny = dcenters[:,1]
    except:
        print t.iso,'SQL connection for delay centers failed.'
        dcen = [0]*16
        dceny = [0]*16
    # ************ End of block *********
    # Replaced by:
    # delaydict = cal_header.ACCdlatable2dict()
    # if delaydict != {}:
        # dcen  = delaydict['Delaycen_ns'][:,0]
        # dceny = delaydict['Delaycen_ns'][:,1]
    # else:
        # print t.iso,'ACC transfer of delay centers failed.  Delay center not updated'
        # dcen = [0]*16
        # dceny = [0]*16
    
    try:
        # Read eovsa_corr.ini file from ACC and get ROACH antenna assignments. 
        inifile = urllib2.urlopen('ftp://'+userpass+'acc.solar.pvt/parm/eovsa_corr.ini',timeout=1)
        antasgn = ''
        for line in inifile.readlines():
            if line.find('antasgn') == 0:
                # Keep appending antennas to the string (final will end in ',').
                antasgn += line.strip().split('=')[1]+','
        # Convert string to numpy array
        antlist = numpy.array(antasgn[:-1].split(',')).astype('int')
    except:
        print t.iso,'ACC connection for eovsa_corr.ini (ROACH antenna assignments) timed out.'
        antlist = numpy.arange(16)   # Assume [bad assumption!] that antennas are assigned in order
    print t.iso, 'Antlist is:', antlist
    
    sh_dict = {'timestamp': timestamp,
               'Version': version,
               'project':'NormalObserving',
               'operator':'Kjell Nelin',
               'comments':'None',
               'version':['1.0.0','1.0.0'],
               'chinfo':ci.Chan_Info(),  # Handle to Chan_Info object
               'nants':16,
               'antlist':antlist,
               #'antlist':numpy.array([1,2,3,4,0,0,0,0,0,0,0,0,0,0,0,0]),  # 1,2,3,4 for prototype
               'antpos': aa,
               'pbfwhm': [1.22*30*180*3600./210./pi]*13 + [1.22*30*180*3600./2700./pi]*2 + [0.0],
               'mount': [1]*13 + [2]*2 + [0],
               'scan_type':'test',
               'source_id':'None',
               'track_mode':'PLANET',
               'epoch':'2000',
               'dra': 0.0, 'ddec': 0.0, 'ha': 0.0, 'dha': 0.0,
               'sun_info': mjd0,    # Solar P-angle, B0-angle, and radius will be calculated for this date
               'pol': [-5,-6,-7,-8],  # Default XX, YY, XY, YX polarization
               'max_file_size':100, # MB, max IDB file size
               'intval': 20,
               'nintval': 50,
               'subbw': 0.4/4096,   # Case of 400 MHz IF bandwidth (not used?)
               'dlacen': numpy.array(dcen),
               'dlaceny': numpy.array(dceny),
               'time_at_acc0': Time(mjdacc0,format='mjd'),
               'roach_brd_clk': [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],
               'katadc':[{},{},{},{},{},{},{},{}]}   # Empty dictionaries (one for each ROACH) for katadc sensor data

#============================
def init_sched_dict():
    '''Create the dictionary that represents the schedule portion of the
       stateframe.  These entries are overwritten once per second as long
       as the schedule is running.
    '''
    global sf_dict
    t = util.Time.now()
    mjd_ = t.mjd      # Get mjd of Time() object
    timestamp = t.lv  # Get LabVIEW timestamp

    # Create schedule dictionary
    sf_dict = {'timestamp': timestamp,
               'timestamp1': timestamp,
               'scan_state': -1,
               'phase_tracking': False,
               'uvw': numpy.array([[0.0,0.0,0.0]]*16),
               'uvw1': numpy.array([[0.0,0.0,0.0]]*16),
               'delay': numpy.zeros(16),
               'delay1': numpy.zeros(16),
               'geosat': None,  # Initally, no geosats defined
               'SolPwr':[{},{}],
               'delays':[{},{},{},{},{},{},{},{}],    # Empty dictionaries (one for each ROACH) for ROACH delays
               'sensors':[{},{},{},{},{},{},{},{}]}   # Empty dictionaries (one for each ROACH) for ROACH sensor data

#============================
def set_uvw(aa,t,srcname):
    '''Set u,v,w coordinates and delays (=-w) [ns] for antenna array aa
       at time specified by datime() object d (1 s in advance of current time).  
       The source name must exist in the source catalog (aa.cat) attached to aa.
    '''
    global sf_dict, sh_dict
         
    if srcname is None:
        if sf_dict['scan_state'] != 1:
            init_sched_dict()    
        return
    mjd_ = t.mjd                 # Get mjd of Time() object (t+1)
    aa.date = mjd_ - 15019.5     # Convert mjd to ephem timebase
    aa.epoch = '2000/1/1 12:00'  # Seems silly to set epoch to J2000 every second, but it keeps getting set to date...
    if sf_dict['geosat']:
        try:
            # Newly calculate coordinates for geosat (an ephem.EarthSatellite object) 
            # and add it to source catalog.  If it already exists, it will be overwritten.
            sat = sf_dict['geosat']
            sat.compute(aa)
            geosat=aipy.amp.RadioFixedBody(sat.ra,sat.dec,name=sat.name)
            # This will either add a new source or replace the existing one
            aa.cat.add_srcs([geosat,geosat])
        except:
            # Probably there is no geosat defined.  Just continue.
            pass

    # Fill in scan header dictionary, so that it will be ready at any time for writing
    # at the start of a scan.  These should be for the current time rather than 1 s in
    # the future, since we have not yet "computed" for new time.
    src = aa.cat[srcname]
    src.compute(aa)
    # Save coordinates as TOPOCENTRIC J2000.  These should agree with JPL Horizons astrometric
    # coordinates for the OVRO site.
    sh_dict['ra'] = (src.ra - src.g_ra) + src.a_ra
    sh_dict['dec'] = (src.dec - src.g_dec) + src.a_dec
    sh_dict['ha'] = eovsa_ha(src)

    cat = aa.cat                           # Make a copy of catalog
    cat.compute(aa)                        # Compute for time t+1
    #print t.iso,aa.epoch
    sf_dict['uvw'] = sf_dict['uvw1'].copy()       # Transfer previously calculated uvw (time t)
    sf_dict['delay'] = sf_dict['delay1'].copy()   # Transfer previously calculated delay (time t)
    sf_dict['timestamp'] = sf_dict['timestamp1']  # Transfer previous timestamp (time t)
    timestamp = t.lv      # Get LabVIEW timestamp (t+1)
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
    t = util.Time.now()
    if line:
        t = util.Time(line[:19])
    return t.mjd

#============================
def get_antlist(key='sun',filename='default.antlist'):
    # antlist = get_antlist(key='sun',filename='default.antlist')
    #
    # Load filename into a dictionary of format {antlistname: antlist}
    # then return the antlist for the key.
    # Supports spaces, commas, etc in antlist. key is not case sensitive,
    # but filename is case sensitive.
    # If key is not found in filename, returns an empty string.
    f = open(filename)
    antlist_dict = {}
    for line in f:
        if len(line.split())<1 or line[0]=='#':  # skip blank or commented lines
            pass
        else:
            linekey = line.split()[0].lower()    # name of antlist
            l = len(linekey)
            antlist = line[l:].strip()           # antlist
            antlist_dict[linekey] = antlist
    f.close()
    try:
        return antlist_dict[key.lower()]
    except:
        print 'Warning: get_antlist() could not find an antlist of name', key,
        'in file', filename + '. In this case get_antlist() returns an empty string.'
        return ''

#============================
class App():
    '''Main application
    '''
    def __init__(self):
        '''Create the user interface and initialize it
        '''
        self.root = Tk()
        self.root.geometry("+100+0")
        
        self.subarray_name = subarray_name # Subarray1 for master schedule, sys.argv[1] for 2nd schedule
        
        # Define self.mypid, since this will be used by get_subarray_pid() and wake_up()
        self.mypid = mypid
        
        # Set function to handle signal alarm, if it should go off.  It is set for 15 s in inc_time().
        signal.signal(signal.SIGALRM, self.wake_up)

        # Read ACC.ini file to get values for global variables
        # binsize, xmlpath, scdport and sfport
        try:
            self.accini = stateframe.rd_ACCfile()
        except urllib2.URLError:
            showerror('Error','ACC unreachable, or ACC.ini file does not exist\n'
                                     +'Cannot continue.')
            exit()
        
        global sf_dict, sh_dict
        sh_folder = '/common/Tables/scanheader/'
        if self.subarray_name == 'Subarray1':
            self.sh_datfile = sh_folder+'scan_header.dat'
        else:
            self.sh_datfile = shfolder+'scan_header_' + self.subarray_name + '.dat'
        if sh_dict == {}:
            init_scanheader_dict()
#            if self.subarray_name == 'Starburst':
#                # write over some defaults (such as for dlacen) with Starburst defaults and add Starburst-specific entries
#                sh_dict_starburst = starburst.init_sh_dict()
#                sh_dict.update(sh_dict_starburst)
#            else:
            scan_header(sh_dict,self.sh_datfile)
        if sf_dict == {}:
            init_sched_dict()
            # Generate schedule item XML file on init
            junk = gen_schedule_sf(sf_dict)#,mk_xml=True)
        
        if self.subarray_name == 'Subarray1':
            # Connect to SQL Server database, or set to None if cannot connect
            # Only do this for the master schedule (Subarray1)
            connstr = "DRIVER={FreeTDS};SERVER=192.168.24.106,1433; \
                                        DATABASE=eOVSA06;UID=aaa;PWD=I@bsbn2w;"
            try:
                sh, shver = stateframe.xml_ptrs(sh_folder+'scan_header.xml')
            except:
                sh, shver = stateframe.xml_ptrs(None)  # Reads from ACC                
            try:
                cnxn = stateframedef.pyodbc.connect(connstr)
                cursor = cnxn.cursor()
                sfbrange, outlist = stateframedef.sfdef(self.accini['sf'])
                shbrange, outlist = stateframedef.sfdef(sh)
            except:
                showerror('Warning','Cannot connect to SQL server\n')
                cnxn = None
                cursor = None
                sfbrange = None
                shbrange = None
            self.sql = {'cnxn':cnxn,'cursor':cursor,'sfbrange':sfbrange,'shbrange':shbrange}
        else:
            self.sql = {'cnxn':None,'cursor':None,'sfbrange':None,'shbrange':None}
        
        # Set window title from command-line argument, which is saved in self.subarray_name
        self.root.wm_title('Schedule for '+self.subarray_name)
        
        timeframe = Frame(self.root)
        timeframe.pack()
        self.no27m = BooleanVar()
        self.CB = Checkbutton(timeframe, text="No 27m",
                              variable=self.no27m)
        self.CB.pack(side=LEFT, expand=0, fill=BOTH)

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
        self.L = Listbox(Lframe,selectmode=EXTENDED,width=45,height=30)
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
        E3.bind('<Up>',self.up)
        E3.bind('<Down>',self.down)
        self.cmd_history = []
        self.ptr = 0

        t = util.Time.now()
        self.label.configure(text=t.iso)
        
        # Set default Ant 14 receiver position to the High-Frequency setting, i.e. lorx = False
        self.lorx = False

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

        # Listbox widget for calibrator times
        CTframe = Frame(fproj)
        CTframe.pack(expand = False, fill=Y,side=LEFT)
        self.CalTimeLB = Listbox(CTframe,selectmode=EXTENDED,width=45,height=10,font="Courier 10 bold")
        out = whenup()
        for line in out['lines']:
            self.CalTimeLB.insert(END,line)
        self.CalTimeLB.pack(side=LEFT)
            
        # Make sure window can never be smaller than initially created size
        self.root.update_idletasks()
        self.root.minsize(self.root.winfo_width(),self.root.winfo_height())

        # Generate the Antenna Array object, used to calculate source coordinates,
        # uvw, delays, etc.
        sys.stdout.flush()
        self.aa = eovsa_cat.eovsa_array_with_cat()
        #print t.iso,self.aa.epoch,' Initial epoch'

        self.aa.epoch = '2000/1/1 12:00'      # Get coordinates in J2000 epoch
        sys.stdout.flush()
         
        #self.Open('solar.scd')
        self.New()  # Generate a solar schedule for today
        self.filename = 'solar.scd'
        
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
        self.PAthread = None
        self.wlimit = 17  # Default wind limit (mph) for 27-m antenna
        self.stale = True # Default status of weather station information (will be immediately set to False if not stale)
        
        # *************** New code due to loss of SQL ***************
        # Log both stateframe and scanheader data to files instead of SQL
        self.accini['sf_file'] = None
        self.accini['sh_file'] = None
        self.log_stateframe()

        # Start the clock ticking
        self.prev = time.time()
        self.tmr = self.root.after(1,self.inc_time)

    #============================
    
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
        # Connect to all ROACHs used by the array.  For Starburst subarray, connect to Starburst ROACHs
        #  instead of OVSA ROACHs.
        
        if self.subarray_name == 'Subarray1': # MASTER SCHEDULE ONLY: connect to OVSA ROACHs
            roachModule = roach
            roach_ips = ('roach1.solar.pvt','roach2.solar.pvt','roach3.solar.pvt','roach4.solar.pvt',
                         'roach5.solar.pvt','roach6.solar.pvt','roach7.solar.pvt','roach8.solar.pvt')
            boffile_name = self.accini['boffile']
            self.brd_clk_freq = None   # Start with no clock defined
            #self.brd_clk_freq = 200
#        elif self.subarray_name == 'Starburst': # STARBURST ONLY: connect to Starburst ROACHs
#            roachModule = starburst.roach
#            roach_ips = starburst.roach.get_roach_ips()
#            boffile_name = starburst.roach.get_boffile_name()
#            self.brd_clk_freq = starburst.roach.get_brd_clk_freq()
        else:  # OVSA SCHEDULE 2: do not connect to any ROACHs (self.roaches will be an empty list)
            roach_ips = ()
        
        if len(self.roaches) != 0:
            # Some roaches are already connected, so stop them and reconnect
            for r in self.roaches:
                if r.fpga: r.fpga.stop()
        self.roaches = []
        
        # This will eventually be a loop over all active ROACHes, and must
        # tolerate a missing ROACH
        for roach_ip in roach_ips:
            # Make connection to ROACHes
            rnum = int(roach_ip[5:6])-1
            #if len(self.roaches) > 0:
                #cfg = self.roaches[0].cfg
                #for line in cfg:
                    #if line.find('adc clock') != -1:
                        #self.brd_clk_freq = int(line.strip().split('=')[1])/4
                        #break
            #else:
                #cfg = None
            #if self.brd_clk_freq is None:
                #self.brd_clk_freq = 300
            self.roaches.append(roachModule.Roach(roach_ip, boffile_name))#, cfg))
            if self.roaches[-1].msg == 'Success':
                print roach_ip,'serving ants',self.roaches[-1].ants,'is reachable'
                self.roaches[-1].dlasweep = None
                try:
                    self.roaches[-1].brd_clk = self.roaches[-1].fpga.est_brd_clk()
                    print roach_ip,'clock is',self.roaches[-1].brd_clk
                    if self.brd_clk_freq is None:
                        # Board clock has not been set yet, so set it once (from first ROACH)
                        self.brd_clk_freq = int(self.roaches[-1].brd_clk)
                    if abs(self.roaches[-1].brd_clk-self.brd_clk_freq)>1:
                        print roach_ip,'clock NOT', self.brd_clk_freq, 'MHz.'
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
    def wake_up(self, signum, frame):
        # This is called whenever the 15-second alarm goes off, indicating the
        # process is stuck in sk_wait.  We simply send ourselves a SIGINT (ctrl-C),
        # which should hopefully do it, but we should also log the fact by setting
        # a flag in the self object.
        self.error = 'The 15-s-alarm went off!'
        print self.error,'Signal:',signum,'at frame',frame
        sys.stdout.flush()
        # Try to reestablish connection to the ROACHes, and set self.fpga accordingly
        # This will keep dla2roach() from hanging.
        self.connect2roach()
        os.kill(int(self.mypid), signal.SIGINT)
        # Kill the inc_time timer and restart it.
        self.root.after_cancel(self.tmr)
        self.prev = time.time()
        self.tmr = self.root.after(1,self.inc_time)


    #============================
    # Raw command history routines

    def send_cmd(self,event):
        '''Send a raw command to ACC.'''
        global sf_dict, sh_dict
        command = event.widget.get()
        if command != '':
            self.execute_ctlline(command)
        self.add2cmdlist(command) 
        event.widget.delete(0,END)

    def add2cmdlist(self,command):
        ''' Adds a command to a 20-element command history
        '''
        if self.cmd_history == [] and command != '':
            # First non-blank command, so append it
            self.cmd_history.append(command)
        elif command == '' or self.cmd_history[-1] == command:
            # Blank command, or same as previous so skip appending
            pass
        else:
            # Unique non-blank command
            if len(self.cmd_history) == 20:
                # More than 20 commands in history, so discard first
                self.cmd_history.pop(0)
            self.cmd_history.append(command)
        self.ptr = 0    # Set pointer to bottom of command history

    def up(self,event):
        # User used up-arrow in Raw Command box
        if len(self.cmd_history) <= -self.ptr:
            # Already at top of command history
            pass
        else:
            # Decrement pointer and insert command string into widget
            self.ptr -= 1
            event.widget.delete(0,END) 
            event.widget.insert(0,self.cmd_history[self.ptr])

 
    def down(self,event):
        # User used down-arrow in Raw Command box
        if self.ptr == 0:
            # Already at bottom of command history
            pass
        else:
            if self.ptr == -1:
                # Pointer was previously at last command, so enter blank
                event.widget.delete(0,END)
                self.ptr = 0
            else:
                # Increment pointer and insert command string into widget
                self.ptr += 1
                event.widget.delete(0,END) 
                event.widget.insert(0,self.cmd_history[self.ptr])

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
                cmds = line[20:].split()
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
            checking any PHASECAL lines to confirm that they will be visible.
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
                ctl_cmd = line.split()[2]
                if ctl_cmd == 'PHASECAL' or ctl_cmd == 'PHASECAL_LO' or ctl_cmd == 'STARBURST' or ctl_cmd == 'CALPNTCAL':
                    name = line.split()[3]
                    # Get start and stop time for this calibrator, as ephem-compatible strings
                    if i == self.lastline-1:
                        # Should never happen, but just in case the PHASECAL
                        # is the last line of the schedule, assume 15 minutes duration
                        trange = util.Time([mjd(line),mjd(line) + 15.*60./86400.],format='mjd')
                    else:
                        line2 = lines[i+1]
                        trange = util.Time([mjd(line),mjd(line2)],format='mjd')
                    if ctl_cmd == 'STARBURST':
                        check_27m = True
                        check_2m = False
                    elif ctl_cmd == 'PHASECAL' or ctl_cmd == 'PHASECAL_LO' or ctl_cmd == 'CALPNTCAL':
                        check_27m = True
                        check_2m = True
                    # Check visibility for source
                    try:
                        src=self.aa.cat[name]
                        visible = scan_visible(src,self.aa,trange,check_27m,check_2m)
                        if not visible:
                            self.error = 'Warning, source '+name+' not visible at scheduled time: Schedule line '+str(i+1)
                            print self.error
                    except:
                        self.error = 'Err: source '+name+' not found in catalog.  Not added to schedule!'

        self.state=['']*len(lines)
        self.status.configure(state = DISABLED)
        if filename is None:
            filenamelist = f.name.split('/')
            self.filename = filenamelist[len(filenamelist)-1:][0]
        else:
            self.filename = filename
        # Update the status file in /common/webplots for display on the status web page
        self.update_status()
        
    #============================
    def New(self):
        # New option creates a new table with predetermined content.
        if self.Toggle == 0:
            # Schedule is running, so do nothing
            return
        self.status.configure( state = NORMAL)
        t = util.Time.now()
        scd = make_sched(t=t)
        if self.no27m.get() == 1:
            scd = remove_cal(scd)
        self.L.delete(0, END)
        for line in scd:
            self.L.insert(END,line)
        self.curline = 0
        self.lastline = len(scd)
        self.status.configure( state = DISABLED)
        self.filename = 'solar.scd'
        
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
            t = util.Time(mjd(line) + delt,format='mjd')
            self.L.delete(i)
            self.L.insert(i,t.iso[:19] + ' ' + line[20:])
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
    def Today(self,t=None):
        ''' The button Today will change the schedule to start on the current date.
            This will NOT change times. Lines that start a day after the start line
            of the schedule will start on today + 1.
        '''
        if t is None:
            t = util.Time.now()
        # If this is the standard solar.scd file, do an auto-generate for today
        if self.filename == 'solar.scd':
            self.New()
        else:
            # Determine how many days from date of first line to today
            line = self.L.get(0)
            days = int(t.mjd) - int(mjd(line))
            print 'Adding ',days,'days.'
            for i in range(self.lastline):
                line = self.L.get(i)
                linemjd = mjd(line) + days
                t2 = util.Time(linemjd,format='mjd')
                line = t2.iso[:10] + line[10:]
                self.L.delete(i)
                self.L.insert(i,line)
            self.curline = 0

    #============================
    def autogen(self,t):
        # Auto-generate the standard solar schedule
        # Determine sunrise, sunset times for this day
        mjdrise, mjdset = suntimes(t)
        # First solar line starts at sunrise
        trise = util.Time(mjdrise,format='mjd')
        line = self.L.get(2)
        line = trise.iso[:19] + line[19:]
        self.L.delete(2)
        self.L.insert(2,line)
        # First calibrator line starts 15 min earlier
        trise = util.Time(mjdrise - 15.*60./86400.,format='mjd')
        line = self.L.get(1)
        line = trise.iso[:19] + line[19:]
        self.L.delete(1)
        self.L.insert(1,line)
        # First calibrator ACQUIRE line starts 20 min earlier
        trise = util.Time(mjdrise - 20.*60./86400.,format='mjd')
        line = self.L.get(0)
        line = trise.iso[:19] + line[19:]
        self.L.delete(0)
        self.L.insert(0,line)
        # Last calibrator line starts at sunset
        tset = util.Time(mjdset,format='mjd')
        line = self.L.get(self.lastline-2)
        line = tset.iso[:19] + line[19:]
        self.L.delete(self.lastline-2)
        self.L.insert(self.lastline-2,line)
        # Last (REWIND) line starts 15 min later
        tset = util.Time(mjdset + 15.*60./86400.,format='mjd')
        line = self.L.get(END)
        line = tset.iso[:19] + line[19:]
        self.L.delete(END)
        self.L.insert(END,line)

        # Read calibrator database
        cal = readvlacaldb()
        # Find sources within 15-35 deg of Sun (also returns antenna array aa, not used)
        srclistnarrow, aa = findcal(cal,t=t,dtheta=[15,35])
        # Early and late in day, need a wider search window, 15-55 degrees of the Sun
        srclistwide, aa = findcal(cal,t=t,dtheta=[15,60])
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
                if i == self.lastline-1:
                    # Should never happen, but just in case the PHASECAL
                    # is the last line of the schedule, assume 15 minutes duration
                    trange = util.Time([mjd(line),mjd(line) + 15.*60./86400.],format='mjd')
                else:
                    line2 = self.L.get(i+1)
                    trange = util.Time([mjd(line),mjd(line2)],format='mjd')
                # Loop over calibrators, highest flux first
                for f in fsortnarrow:
                    # These are 15-min observations, so make sure calibrator
                    # is visible and will remain visible
                    idx = fluxnarrow.index(f)
                    jys = fluxnarrow[idx]
                    src = srclistnarrow[idx]
                    visible = scan_visible(src,self.aa,trange,True,True)
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
                        visible = scan_visible(src,self.aa,trange,True,True)
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
                            visible = scan_visible(src,self.aa,trange)
                            print src.name, jys, visible, src.ra, src.dec, '     ', src.az, src.alt
                self.L.delete(i)
                self.L.insert(i,line)
                if i == 1:
                    # Check if earlier line is ACQUIRE, and if so, make sure it goes to the same
                    # calibrator!
                    acline = self.L.get(0)
                    if acline[20:27] == 'ACQUIRE':
                        acline = acline[:28] + src.name
                        self.L.delete(0)
                        self.L.insert(0,acline)

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
            cmds = line[20:].split()
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
    def check_27m_sun(self,sf,data):
        pass
    def check_27m_sun2(self,sf,data):    # Disable for now by renaming (2017-09-05  DG)
        # For each 27m, check whether it is too close to Sun (min_dist: 10 degrees);
        # if it is, use 'position' command to send to stow position (HA = 0, dec = 29)
        #
        # If ant is in runmode TRACK or POSITION, then trigger if Requested position
        #   is within min_dist of Sun OR if Actual position is within min_dist of Sun
        #   for 3 seconds (it should take < 2 sec to slew across Sun, which is okay)
        #   - timer restarts if this function sends a position command for the antenna
        # If ant is in runmode STOP or VELOCITY, then trigger only if Actual position
        #   is within min_dist of Sun for 3 seconds (do not trigger based on Requested
        #   position because Requested position is not accurate in these modes)
        #
        # If all Controller monitor data from antenna is zero, does not make check
        #
        # Stow position (0 29) can be 6 degrees from Sun in summer.  Code does not trigger
        #   position command if antenna is already within 0.5 degrees of stow position.
        #
        # Input params: sf is accini['sf'], data is from get_stateframe (which was just run by inc_time)
        
        global sf_dict
        global sun_timer
        global trigger_timer
        
        try:
            k = sun_timer.keys()
        except:
            sun_timer = {}
        try:
            k = trigger_timer.keys()
        except:
            trigger_timer = {}
        
        min_dist = 10.    # minimum_distance from Sun
        max_time = 2.     # maximum time near Sun (to allow slewing)
        safe_pos = [0,29] # safe position to send antenna to - stow position
        trigger_cadence = 20 # minimum time between sending command in response to trigger
        test_mode = False  # prints a bunch of diagnostic info to log if this is True

        # define dict to convert runmode ID# from controller to meaning
        runmode_map = {0: 'stop', 1: 'position', 2: 'velocity', 4: 'track'}
        
        # get coords of Sun
        sun = self.aa.cat['Sun']
#        sun = self.aa.cat['AMC-8 (GE-8)'] # for testing purposes can put name of geosat to avoid
        self.aa.date = util.Time.now().mjd - 15019.5  # set date to present time
        sun.compute(self.aa)
        sun_coords = (sun.ra,sun.dec)
        lst = self.aa.sidereal_time() # current LST
        
        # calc RA,dec of safe pos (to avoid triggering when within 0.5 deg of safe pos)
        safe_coords = (lst-safe_pos[0]*pi/180.,safe_pos[1]*pi/180.)
        
        for antnum in [14]:
            # augment trigger_timer and skip this antenna if it's less than trigger_cadence since last trigger
            try:
                trigger_timer[antnum] = trigger_timer[antnum] + 1
            except KeyError:
                trigger_timer[antnum] = 1000   # start high so that it won't wait 20s before the first trigger
            if trigger_timer[antnum] < trigger_cadence:
                if test_mode: print 'Skipping ant '+ str(antnum) + ' because <' +str(trigger_cadence)+' sec since last trigger'
                continue
            
            # skip this antenna if all crio monitor data for this ant is zero
            c = sf['Antenna'][antnum-1]['Controller']
            sflist = array([stateframe.extract(data,c[k]) for k in c.keys()])
            if len(find(sflist != 0)) == 0:
                if test_mode: print 'Skipping ant '+ str(antnum)+' because no stateframe data from controller'
                continue
            
            trigger = False  # set this to True if too close to Sun
            
            # get runmode from stateframe data passed to me
            rm = stateframe.extract(data,c['RunMode'])
            runmode = runmode_map[rm]
            
            if test_mode:
                print '-----------'
                print 'LST:', lst, '- Sun coords:', sun_coords
                print 'Antnum:', antnum
                print 'Runmode:', rm, runmode
            
            # determine whether to trigger moving to safe position
            
            # all runmodes: trigger if Actual pos near Sun for > max_time (2 sec)
            # check distance of Actual position from Sun
            HA_actual = sf_dict['ActualAzimuth'][antnum-1] * pi/180.
            dec_actual = sf_dict['ActualElevation'][antnum-1] * pi/180.
            actual_coords = (lst-HA_actual,dec_actual)
            actual_dist = ephem.separation(actual_coords,sun_coords)*180./pi
            if actual_dist > min_dist:
                # if Actual position is not too near Sun, reset timer to zero
                sun_timer[antnum] = 0.
            else:
                try:
                    t = sun_timer[antnum]
                except KeyError:
                    if test_mode: print 'KeyError! setting suntimer to 0'
                    sun_timer[antnum] = 0.
                # if Actual position is too close to Sun, increment timer by one second
                sun_timer[antnum] = sun_timer[antnum] + 1.
            if sun_timer[antnum] > max_time:
                # trigger if Actual position has been close to Sun for more than max_time
                # (and not within 0.5 degrees of safe_pos)
                safe_actual_dist = ephem.separation(actual_coords,safe_coords)
                if test_mode: print 'Dist between stow and actual position:', safe_actual_dist
                if safe_actual_dist > 0.5:
                    trigger = True

            if test_mode:
                print 'Actual coords:', actual_coords, '- Actual dist:', actual_dist
                print 'Sun timer:', sun_timer[antnum]
                print 'Trigger based on Sun timer:', trigger
                
            # runmodes POSITION and TRACK: trigger if Requested pos too near Sun
            if runmode in ['track','position']:
                HA_requested = sf_dict['RequestedAzimuth'][antnum-1] * pi/180.
                dec_requested = sf_dict['RequestedElevation'][antnum-1] * pi/180.
                requested_coords = (lst-HA_requested,dec_requested)
                requested_dist = ephem.separation(requested_coords,sun_coords)*180./pi
                if test_mode: print 'Requested coords:', requested_coords, '- Requested dist:', requested_dist
                if requested_dist < min_dist:
                    # trigger if Requested position is too close to Sun and not within 0.5 degrees of stow
                    safe_requested_dist = ephem.separation(requested_coords,safe_coords)
                    if test_mode: print 'Dist between stow and requested position:', safe_requested_dist
                    if safe_requested_dist > 0.5:
                        trigger = True
            if test_mode: print 'Trigger:', trigger
            
            if trigger:
                trigger_timer[antnum] = 0. # reset timer so antenna has time to process command and slew off Sun
                # position command causes antenna to go to safe_pos, enter 'position'
                #  runmode, but it will not lose its tracktable
                safe_str = str(safe_pos[0]) + ' ' + str(safe_pos[1])
                stop_cmd = 'stop ant' + str(antnum)
                pos_cmd = 'position ' + safe_str + ' ant' + str(antnum)
                self.sendctlline(stop_cmd)
                self.sendctlline(pos_cmd)
                print 'Antenna ' + str(antnum) + ' actual or requested position too close to Sun, sending to safe position: ' + safe_str
                print stop_cmd
                print pos_cmd
            
        sys.stdout.flush()
                
    #============================
    def inc_time(self):
        global sf_dict, sh_dict

        self.status.configure(state=NORMAL)

        # First set the timer to wake us on the next second
        t = util.Time.now()
        tdif = int((t.datetime.microsecond)/1000.)
        self.root.after(1000 - tdif, self.inc_time)
        # Update the clock
        self.label.configure(text=t.iso[:19])

        # Set an alarm for 15 seconds.  If the process hangs for more than that, we will send ourselves
        # a SIGINT via the Callback self.wake_up(), which should recover from sk_wait hang up
        signal.alarm(15)  # This will reset alarm if it has not gone off yet
        tnow = time.time()
        self.telapsed = tnow - self.prev
        self.prev = tnow
        telapsed = int(self.telapsed*1000)
        if telapsed > 990 and telapsed < 1010:
            pass
        else:
            # If elapsed time is not nominal (e.g. 990 or 1010), write it to log file.
            print t.iso,str(int(self.telapsed*1000))
            sys.stdout.flush() # Flush stdout (/tmp/schedule.log or /tmp/schedule_[self.subarray_name].log) so we can see this '-'.

        # Attempt to read from spawned task pwr_cycle.ant_toggle() queue.  Reads up to 10
        # items at a time unless queue is empty.
        for i in range(10):
            try:
                msg = pwr_cycle.q.get_nowait()
                print t.iso,msg
            except:
                break
        # Attempt to read from spawned task pcapture2.capture_1s() queue.
        try:
            msg = pcapture2.q.get_nowait()
            print t.iso,msg
        except:
            pass
            
        # If this schedule is running the second subarray, confirm that Subarray1 is running; if it is not,
        # print a warning message to the log file.
        if self.subarray_name != 'Subarray1':
            subarray1_pid = self.get_subarray_pid('Subarray1')
            if subarray1_pid == -1:
                print util.datime().get('str'), \
                      'Warning: The master schedule (Subarray1) is not running.  This means that antenna diagnostic ' + \
                      'information will not be updated in the ACC stateframe and no data will be written to the SQL database.'
                sys.stdout.flush() # Flush stdout (/tmp/schedule.log or /tmp/schedule_[self.subarray_name].log) once per second so we can see the output.

        # Update phase tracking (u,v,w and delays)
        srcname = sh_dict['source_id']
        try:
            # Generate a Time() object at exactly the next upcoming second (time t+1)
            t2 = util.Time.now()
            tsec = util.Time(t2.mjd  + (1 - t2.datetime.microsecond/1000000.)/86400.,format='mjd')
            src = self.aa.cat[srcname]        # This causes KeyError if source is not found
        except KeyError:
            # The current scan header source ID is not in the source catalog
            srcname = None
        # Debug info, simply logs that we have started this procedure
        #sys.stdout.write('+')
        #sys.stdout.flush() # Flush stdout (/tmp/schedule.log or /tmp/schedule_[self.subarray_name].log) so we can see this '-'.
        set_uvw(self.aa,tsec,srcname)
        #sys.stdout.write('-')
        #sys.stdout.flush() # Flush stdout (/tmp/schedule.log or /tmp/schedule_[self.subarray_name].log) so we can see this '-'.
        self.source.configure(text='    Source: '+(str(srcname)+'            ')[:12]
                                  +'Phase Tracking: '
                                  +str(bool(sf_dict['phase_tracking'])) + '    '+self.error)
        self.error = ''
       
        # Send integer delays to ROACHs - DIFFERENT FOR STARBURST
#        if self.subarray_name == 'Starburst':
#            starburst.roach.dla2roach(self,sh_dict,sf_dict)
#        else:
        self.dla2roach()

        # Update weather information in sf_dict (reads from OVRO weather station)
        self.w = stateframe.weather()
        sf_dict.update(self.w)
        # If weather information is "stale" (older than 5 minutes), set wscram-limit to 0 to force
        # Ant 14 to be kept stowed.
        try:
            tdifw = t - Time(self.w['mtSampTime'].replace('/','-'))
            if tdifw.value > 300./86400.:
                if self.stale is False:
                    self.stale = True
                    # Open socket to ACC
                    s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
                    try:
                        # Send commands to update antenna trip information
                        s.connect((self.accini['host'],self.accini['scdport']))
                        s.send('WSCRAM-LIMIT 0 ANT14')
                        time.sleep(0.01)
                        s.close()
                    except:
                        pass
            else:
                if self.stale is True:
                    self.stale = False
                    # Open socket to ACC
                    s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
                    try:
                        # Send command to update windscram limit for Ant 14
                        s.connect((self.accini['host'],self.accini['scdport']))
                        s.send('WSCRAM-LIMIT '+str(self.wlimit)+' ANT14')
                        time.sleep(0.01)
                        s.close()
                    except:
                        pass
        except:
            # The above calculation of tdif failed--probably a glitch in reading the weather, so leave state as is
            pass

        # Once per minute, update the information from the Solar Power station(s)
        if t.datetime.second == 0:
            # Updates first solar power station on the minute
            self.solpwr[0] = stateframe.rd_solpwr('http://data.magnumenergy.com/MW5127')
        if t.datetime.second == 1:
            # Updates second solar power station one second later
            self.solpwr[1] = stateframe.rd_solpwr('http://data.magnumenergy.com/MW5241')
        sf_dict.update({'SolPwr':self.solpwr})

        # Read ROACH sensor data, but only one each minute, staggered over different times
        # since for all 8 ROACHes this can take more than 0.5 s
        for i in range(len(self.roaches)):
            if t.datetime.second == 5*i+5:
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
        self.cr_temp = stateframe.control_room_temp()
        self.sensors[0]['temp.ambient'] = self.cr_temp

        # MASTER SCHEDULE ONLY: Update antenna diagnostics, but only once every 5 minutes (300 s)
        if self.subarray_name == 'Subarray1':
            if int(t.mjd * 86400.) % 300 == 30:
                # Open socket to ACC
                s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
                try:
                    # Send commands to update antenna trip information
                    s.connect((self.accini['host'],self.accini['scdport']))
                    s.send('UPDATEAZIMUTHDIAGNOSTICS 1')
                    time.sleep(0.01)
                    s.send('UPDATEELEVATIONDIAGNOSTICS 1')
                    time.sleep(0.01)
                    s.close()
                except:
                    pass
        
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

#        # STARBURST ONLY: update sf_dict with Starburst-specific monitor data
#        if self.subarray_name == 'Starburst':
#            sf_dict_starburst = starburst.get_sf_dict(self)
#            sf_dict.update(sf_dict_starburst) # make sure that entries in sf_dict_starburst have unique keys so that we don't overwrite sf_dict entries

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
                # MASTER SCHEDULE ONLY: If we are connected to the SQL database, we will need to create and send
                # a new stateframe definition
                if self.subarray_name == 'Subarray1':
                    result = stateframedef.load_deftable(sdict=self.accini['sf'],version=version)
                    if not result:
                        sys.stdout.write('Error loading new stateframe definition')
                        sys.stdout.flush()
            sf = self.accini['sf']
            sf_dict.update(stateframe.azel_from_stateframe(sf,data))
            # Flag unused antennas as not tracking (no!  this is just the ROACH assignments)
            # sf_dict['TrackFlag'] = (sf_dict['TrackFlag']) & (sh_dict['antlist'] != 0)
            
            # Check that 27-m antennas are not too close to Sun, and if they are, send them to stow position (using 'position' command)
            self.check_27m_sun(sf,data)
        else:
            self.error = msg
                
        # ************ This block commented out due to loss of SQL **************
        # If we are connected to the SQL database, send converted stateframe (only master schedule is connected)
        if msg == 'No Error' and self.sql['cnxn']:
            self.sql['cursor'].execute('set implicit_transactions off')
            bufout = stateframedef.transmogrify(data, self.sql['sfbrange'])
            try:
                self.sql['cursor'].execute('insert into fBin (Bin) values (?)', 
                                           stateframedef.pyodbc.Binary(bufout))
                #sys.stdout.write('*')
                #sys.stdout.flush()
                self.sql['cnxn'].commit()
            except Exception as e:
                print e
                # An exception could be an error, or just that the entry was already inserted
                self.error = 'Err: Cannot write stateframe to SQL'
        # ********** End of block **************
        f = self.accini.get('sf_file',None)   
        if f:
            date_change = int(time.time() / 86400) > int(os.path.getctime(f.name) / 86400)
            if date_change:
                # Looks like the date has changed, so open a new file
                self.log_stateframe()
                f = self.accini.get('sf_file')
        try:
            if msg == 'No Error': 
                f.write(data)
        except:
            print Time.now().iso+' Error writing stateframe to log file'

        # Create schedule part of stateframe from sf_dict
        # Subarray1 writes Weather, SolarPower, Roach whereas Subarray2/Starburst don't
        if self.subarray_name == 'Subarray1':
            fmt, buf, sched_xmlfile = gen_schedule_sf(sf_dict)
#        else:
#            # SUBARRAY2 (OVSA) AND STARBURST: use starburst module's gen_schedule2_sf to create binary buffer to write to ACC
#            # if it is Subarray2 (OVSA), starburst.gen_schedule2_sf writes default values for the Starburst-specific data
#            fmt, buf, sched_xmlfile = starburst.gen_schedule2_sf(sf_dict) # using mk_xml=False, so make sure schedule2_stateframe.xml already exists

        # Open socket to ACC - PORT DEPENDS ON SUBARRAY
        if self.subarray_name == 'Subarray1':
            portkey = 'scdsfport'
        else:
            portkey = 'scd2sfport'
        s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        try:
            # Try to connect and send schedule items of stateframe to ACC
            # Uses "schedule" command and encloses data buffer in square brackets
            time.sleep(0.01)
            s.connect((self.accini['host'],self.accini[portkey]))
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
                    print 'Waiting...',self.wait
                    sys.stdout.flush()
                    if self.wait == 0:
                        self.execute_cmds()
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
        #sys.stdout.flush()  # Flush stdout (/tmp/schedule.log or /tmp/schedule_[self.subarray_name].log) so we can see this '-'.

    def log_stateframe(self):
        '''Called on init, or when it is time to close a log file and open 
           a new one.  This logs both stateframe and scanheader.  This is
           only needed if the SQL server is down.
        '''
        global sf_dict, sh_dict
        # Create file name from date
        # Get today's date
        t = Time.now()
        v = self.accini['version']
        logfile = '/mnt/data1/Tables/stateframe/sf_'+t.iso[:10].replace('-','')+'_v'+str(v)+'.log'
        if os.path.isfile(logfile):
            # Desired file name already exists, so simply append to it
            self.accini['sf_file'] = open(logfile,'a')
        else:
            # Need to open a new file, so first check if old file is open
            if self.accini['sf_file']:
                self.accini['sf_file'].close()
            self.accini['sf_file'] = open(logfile,'wb')
        # Create scanheader log filename from date
        # Get the current version number from accini
        v = sh_dict['Version']
        logfile = '/common/Tables/scanheader/sh_'+t.iso[:10].replace('-','')+'_v'+str(v)+'.log'
        if os.path.isfile(logfile):
            # Desired file name already exists, so simply append to it
            self.accini['sh_file'] = open(logfile,'a')
        else:
            # Need to open a new file, so first check if old file is open
            if self.accini['sh_file']:
                self.accini['sh_file'].close()
            self.accini['sh_file'] = open(logfile,'wb')


    #============================
    def dla2roach(self):
        '''Set integer delays and send to all ROACHes for next second, to be
           ready for next 1 PPS.  Ant 1 delay is fixed at 9000 steps (delay
           depends on step size), and all other antennas are set to
           9000 + t_cen[i] + t_geom[i]/step_size, the latter being calculated
           from sf_dict['delay1'][i] (nsec).  The step_size is 1/f_ADC,
           where f_ADC is ADC clock frequency, in GHz.
        '''
        global sh_dict, sf_dict
        # Delay centers.  These are read from the SQL database whenever a scan starts.
        # The delay center offset (dlaoff) is now scaled to adc_clk, so that it is
        # 11,000 for 800 MHz, and 16,500 for 1200 MHz. 
        # ****** Send delay0, which is actually the delay for the next upcoming second ******
        adc_clk = self.brd_clk_freq*4./1000.
        dlaoff = int(16500.*adc_clk/1.2)

        # This swap is now done at $SCAN-START (and written to ACC)
        #dcenidx = numpy.arange(16)
        #if self.lorx:
        #    # If the low-frequency receiver is in place (i.e. an RX-SELECT LO ANT14 command
        #    # was sent), replace the Ant 14 delay centers with those in the Ant 15 slot.
        #    dcenidx[13:15] = [14,13]
        #dlax = numpy.round((sh_dict['dlacen'][dcenidx] - sf_dict['delay'])*adc_clk + dlaoff)
        #dlay = dlax + (sh_dict['dlaceny'] - sh_dict['dlacen'])[dcenidx]
        dlax = numpy.round((sh_dict['dlacen'] - sf_dict['delay'])*adc_clk + dlaoff)
        dlay = dlax + (sh_dict['dlaceny'] - sh_dict['dlacen'])
        
        for r in self.roaches:
            if r.fpga:
                a1,a2 = r.ants
                # Handle the case of a ROACH that has been set to sweep delays
                # The dlasweep dictionary is normally None, but can be {'ant':n,'dla':a,'dlastop':b}
                # which will sweep the delay by 1 step per second for delays from a to b relative
                # to the current nominal delay.
                if r.dlasweep is None:
                    r.set_delays([dlax[a1-1],dlay[a1-1],dlax[a2-1],dlay[a2-1]])
                else:
                    if r.dlasweep['ant'] == a1:
                        # Increment delay by one
                        r.dlasweep['dla'] += 1
                        dla = r.dlasweep['dla']
                        if r.dlasweep['dla'] > r.dlasweep['dlastop']:
                            # The delay is greater than the stop delay, so cancel sweep
                            # by setting dictionary to None
                            dla = 0
                            r.dlasweep = None
                        else:
                            if r.dlasweep['pol'] == 'X':
                                r.set_delays([dlax[a1-1]+dla,dlay[a1-1],dlax[a2-1],dlay[a2-1]])
                                print 'DLASWEEP Ant',a1,'X delay',dla
                            elif r.dlasweep['pol'] == 'Y':
                                r.set_delays([dlax[a1-1],dlay[a1-1]+dla,dlax[a2-1],dlay[a2-1]])
                                print 'DLASWEEP Ant',a1,'Y delay',dla
                            else:
                                r.set_delays([dlax[a1-1]+dla,dlay[a1-1]+dla,dlax[a2-1],dlay[a2-1]])
                                print 'DLASWEEP Ant',a1,'X and Y delay',dla
                    elif r.dlasweep['ant'] == a2:
                        # Increment delay by one
                        r.dlasweep['dla'] += 1
                        dla = r.dlasweep['dla']
                        if r.dlasweep['dla'] > r.dlasweep['dlastop']:
                            # The delay is greater than the stop delay, so cancel sweep
                            # by setting dictionary to None
                            dla = 0
                            r.dlasweep = None
                        if r.dlasweep['pol'] == 'X':
                            r.set_delays([dlax[a1-1],dlay[a1-1],dlax[a2-1]+dla,dlay[a2-1]])
                            print 'DLASWEEP Ant',a2,'X delay',dla
                        elif r.dlasweep['pol'] == 'Y':
                            r.set_delays([dlax[a1-1],dlay[a1-1],dlax[a2-1],dlay[a2-1]+dla])
                            print 'DLASWEEP Ant',a2,'Y delay',dla
                        else:
                            r.set_delays([dlax[a1-1],dlay[a1-1],dlax[a2-1]+dla,dlay[a2-1]+dla])
                            print 'DLASWEEP Ant',a2,'X and Y delay',dla
                    elif r.dlasweep['ant'] == 0:
                        # If ant is 0, sweep delays for both antennas 
                        #             (used for total power polarization tests)
                        # Increment delay by one
                        r.dlasweep['dla'] += 1
                        dla = r.dlasweep['dla']
                        if r.dlasweep['dla'] > r.dlasweep['dlastop']:
                            # The delay is greater than the stop delay, so cancel sweep
                            # by setting dictionary to None
                            dla = 0
                            r.dlasweep = None
                        if r.dlasweep['pol'] == 'X':
                            r.set_delays([dlax[a1-1]+dla,dlay[a1-1],dlax[a2-1]+dla,dlay[a2-1]])
                            print 'DLASWEEP Ant',a1,'and',a2,'X delay',dla
                        elif r.dlasweep['pol'] == 'Y':
                            r.set_delays([dlax[a1-1],dlay[a1-1]+dla,dlax[a2-1],dlay[a2-1]+dla])
                            print 'DLASWEEP Ant',a1,'and',a2,'Y delay',dla
                        else:
                            r.set_delays([dlax[a1-1]+dla,dlay[a1-1]+dla,dlax[a2-1]+dla,dlay[a2-1]+dla])
                            print 'DLASWEEP Ant',a1,'and',a2,'X and Y delay',dla
                            
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
    def sequence2dcmtable(self,sequence):
        '''Use frequency sequence to set dcmtable.txt and send to ACC
           The sequence numbers are 0-based band numbers
        '''
        # Convert from comma-separated variables to zero-based band numbers
        bands = numpy.array(sequence.split(',')).astype('int')-1
        # ************ This block commented out due to loss of SQL **************
        # Read current DCM_Master_Table
        dcm, buf = cal_header.read_cal(2)
        dcm_m_attn = stateframe.extract(buf,dcm['Attenuation'])
        dcm_attn = dcm_m_attn[bands]
        # ************ End of block ***********
        # Replaced by
        # dcm_dict = cal_header.ACC_DCMtable2dict()
        # dcm_attn = dcm_dict['DCMattn'][bands]

        # Rest is unchanged
        lines = []
        for line in dcm_attn:
            l = ' '.join(map(str,line))
            lines.append(l)
        # On 2021 May 09, it was determined that the attenuation table is only
        # applied in the DCM after some lag, which at present is 0 slots.
        lag = 0   # Change this if needed to account for a different lag (lag=0 for no lag)
        # This rotates the lines list by the number of entries given by lag
        lines = lines[-lag:] + lines[:-lag]
#        for i in range(lag):
#            lines.append(lines.pop(0))
        g = open('DCM_table.txt','w')
        for line in lines:
            g.write(line+'\n')
        g.close()
        # ************ This line commented out due to loss of SQL **************
        cal_header.dcm_table2sql(lines)         # Note that I don't think this caltype is used anyway
        
        # Connect to ACC /parm directory and transfer dcm.txt file
        try:
            g = open('DCM_table.txt','r')
            acc = FTP('acc.solar.pvt')
            acc.login('admin','observer')
            acc.cwd('parm')
            # Send DCM table lines to ACC
            print acc.storlines('STOR dcm.txt',g)
            g.close()
            print 'Successfully wrote dcm.txt to ACC'
        except:
            print 'Cannot FTP dcm.txt to ACC'
    
    #============================
    def execute_cmds(self):
        '''Execute the atomic commands associated with the current line
           of the schedule.  First read the commands from the associated
           file and enter them into the L2 Listbox.  Then read them one
           at a time from the Listbox and execute them.
        '''
        global sf_dict, sh_dict
        # Update the status file in /common/webplots for display on the status web page
        self.update_status()
        # Get time range of this Macro command
        mjd1 = mjd(self.L.get(self.curline))
        mjd2 = mjd(self.L.get(self.curline+1))
        # Find the file associated with the Macro command on the current 
        # line and fill in the L2 Listbox
        line = self.L.get(self.curline)
        cmds = line[20:].split()
        f2 = open(cmds[0].rstrip()+'.ctl')
        self.L2.delete(0,END)
        # Current options for source ID
        if cmds[0].upper() == 'SUN': 
            sh_dict['project'] = 'NormalObserving'
            sh_dict['source_id'] = 'Sun'
            sh_dict['track_mode'] = 'PLANET'
        elif cmds[0].upper() == 'SOLPNTCAL':
            sh_dict['project'] = 'SOLPNTCAL'
            sh_dict['source_id'] = 'Sun'
            sh_dict['track_mode'] = 'PLANET'
        elif cmds[0].upper()[:5] == 'FLARE':
            sh_dict['project'] = cmds[0].upper()
            sh_dict['source_id'] = 'Sun'
            sh_dict['track_mode'] = 'PLANET'
        elif cmds[0].upper() == 'PLANET':
            sh_dict['project'] = 'PLANET'
            sh_dict['source_id'] = cmds[1]
            sh_dict['track_mode'] = 'PLANET'
        elif cmds[0].upper() == 'FEATTNTEST':
            sh_dict['project'] = 'FEATTNTEST'
            sh_dict['source_id'] = 'Sun'
            sh_dict['track_mode'] = 'PLANET'
        elif cmds[0].upper() == 'PHASECAL':
            sh_dict['project'] = 'PHASECAL'
            sh_dict['source_id'] = cmds[1]
            sh_dict['track_mode'] = 'RADEC '
        elif cmds[0].upper() == 'PHASECAL_LO':
            sh_dict['project'] = 'PHASECAL'
            sh_dict['source_id'] = cmds[1]
            sh_dict['track_mode'] = 'RADEC '
        elif cmds[0][:5].upper() == 'PACAL':
            sh_dict['project'] = 'PHASECAL'
            sh_dict['source_id'] = cmds[1]
            sh_dict['track_mode'] = 'RADEC '
        elif cmds[0].upper() == 'CALPNTCAL':
            sh_dict['project'] = 'CALPNTCAL'
            sh_dict['source_id'] = cmds[1]
            sh_dict['track_mode'] = 'RADEC '
        elif cmds[0].upper() == 'STARBURST':
            sh_dict['project'] = 'STARBURST'
            sh_dict['source_id'] = cmds[1]
            sh_dict['track_mode'] = 'RADEC '
            print 'Source is',cmds[1]
        elif cmds[0].upper() == 'GEOSAT' or cmds[0].upper() == 'DELAYCAL':
            sh_dict['project'] = 'GEOSAT'
            sh_dict['source_id'] = cmds[1].replace('_',' ')
            # These are geostationary satellites so far.  If/when we add
            # moving satellite capability, track_mode for those should be 'SATELL'
            sh_dict['track_mode'] = 'FIXED '
            try:
                f = urllib2.urlopen('http://www.celestrak.com/NORAD/elements/geo.txt',timeout=20)
                lines = f.readlines()
            except:
                print util.Time.now().iso,'Connection to Celestrak timed out.'
                sh_dict['source_id']='None'
                lines = ['']
            for i,line in enumerate(lines):
                 if line.find(sh_dict['source_id']) == 0:
                     break
            if i < len(lines):
                # This creates an ephem.EarthSatellite object, which does the
                # right thing in calculating coordinates when the time in aa is updated
                sat=ephem.readtle(lines[i],lines[i+1],lines[i+2])
                sf_dict['geosat']=sat
                sat.compute(self.aa)
                # Unfortunately, aipy cannot deal with an ephem.EarthSatellite object,
                # so this creates a fake RadioFixedBody for the current RA,Dec of the 
                # satellite, to be added to the source catalog. This has to be updated 
                # once per second in set_uvw()
                geosat=aipy.amp.RadioFixedBody(sat.ra,sat.dec,name=sat.name)
                self.aa.cat.add_srcs([geosat,geosat])
            else:
                print 'Geosat named ',sh_dict['source_id'],'not found!'
                sh_dict['source_id']='None'
        else:
            # Default project is just the first command on line (truncate to 32 chars)
            sh_dict['project'] = cmds[0][:32]
            print 'Default project:',cmds[0][:32]
            if len(cmds) == 1:
                # Case of only one string on command line
                sh_dict['source_id'] = 'None'
            else:
                # Default source ID is second string on command line (truncate to 12 chars)
                sh_dict['source_id'] = cmds[1][:12]
            print 'Default source:',sh_dict['source_id']
            sh_dict['track_mode'] = 'FIXED '        
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
            # After executing, turn off waitmode
            self.waitmode = False
#            self.nextctlline = 0
        else:
            sline = 0
        for i in range(sline,len(lines)):
            ctlline = self.L2.get(i)
            self.execute_ctlline(ctlline,mjd1,mjd2)
            if ctlline.split()[0].upper() == '$WAIT':
                self.nextctlline = i+1
                break

    def sendctlline(self,ctlline):
        try:
            s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
            s.connect((self.accini['host'],self.accini['scdport']))
            s.send(ctlline)
            time.sleep(0.01)
            s.close()
        except:
            pass
            
    def interpret_pcycle(self, ctlline):
        ''' Interprets the control line containing a $PCYCLE command, of the
            form $PCYCLE <device> <antlist>, where <cevice> is one of 'ant',
            'fem', 'frontend', 'crio', and <antlist> is the usual antenna list.
            This is all NOT case-sensitive.  Note that only the first three
            characters of the device name are examined.  If the device is 'ant',
            then the device name is optional, i.e.
               $PCYCLE ant1 ant2-4 ant7
            will work.
            
            The device and antenna list are returned.  If line cannot be
            interpreted, device is None
        '''
        valid = {'ANT':'ANT','FEM':'FRONTEND','FRO':'FRONTEND','CRI':'CRIO','OTH':'OTHER'}
        tokens = ctlline.upper().split()
        # First token is just the command.  Check whether the second token is valid:
        if len(tokens) == 1:
            return None, None
        try:
            # Check that first three characters of second token matches valid
            # devices, and set device to full spelling.
            device = valid[tokens[1][:3]]
        except KeyError:
            return None, None
        if device != 'ANT': tokens = tokens[1:]
        # Now check remaining tokens for antenna number
        ants = []
        for token in tokens[1:]:
            try:
                # Try to interpret ANTn where n is an integer
                junk, num = token.split('ANT')
                if num != '': ants.append(int(num))
            except:
                try:
                    # Try to interpret ANTm-n, where m and n are integers
                    a1, a2 = num.split('-')
                    for i in range(int(a1),int(a2)+1):
                        ants.append(i)
                except:
                    return None, None
        return device, ants
        
    #============================
    def execute_ctlline(self,ctlline,mjd1=None,mjd2=None):
        # Send line to ACC.  Lines that start with '$' will be entered
        # into stateframe without execution by ACC.  Lines that start
        # with '$*' are Starburst commands.
        # Skip comments
        if len(ctlline) == 0:
            print util.Time.now().iso[:19],'Empty line in .ctl file?'
            pass
        elif ctlline[0] == '#':
            pass
        else:
            if ctlline.strip().upper() == 'DCMAUTO-ON':
                self.sendctlline('DCMTABLE DCM.TXT')
            self.sendctlline(ctlline)
#            try:
#                s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
#                s.connect((self.accini['host'],self.accini['scdport']))
#                s.send(ctlline)
#                time.sleep(0.01)
#                s.close()
#            except:
#                pass
            if ctlline[0] == '$':
                # This line starts with '$', so execute locally
#                if ctlline[1] == '*':
#                    # This is a Starburst-specific command, pass it to starburst module to handle
#                    # Pass it sh_dict and sf_dict in case it needs to modify values - dicts are mutable so should
#                    #   be able to do this
#                    starburst.execute_ctlline(self,ctlline,sh_dict,sf_dict,mjd1,mjd2)
                #==== MK_TABLES ====
                if ctlline.split()[0].upper() == '$MK_TABLES':
                    cmd, fname, src = ctlline.split()
                    if mjd1 is None:
                        d = util.datime()
                        mjd1 = d.get()
                        mjd2 = mjd1+1
                    if fname.upper() == 'GEOSAT_TAB':
                        tbl = make_geosattable(sf_dict['geosat'],self.aa,mjd1,mjd2)
                    else:
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
                    acc.login('admin','observer')
                    acc.cwd('parm')
                    acc.storlines('STOR '+fname,f)
                    acc.close()
                    f.close()
                    #userpass = 'admin:observer@'
                    f = urllib2.urlopen('ftp://'+userpass+'acc.solar.pvt/parm/'+fname)
                    tbl_echo = ''
                    for line in f.readlines():
                        tbl_echo += line.rstrip()+'\n'
                    if tbl != tbl_echo:
                        print 'Error: Transfer of track table',fname,'failed!'
                        print tbl,'not equal\n',tbl_echo
                #==== FEM-INIT ====
                elif ctlline.split()[0].upper() == '$FEM-INIT':
                    ant_str = 'ant1-13'
                    t = threading.Thread(target=adc_cal2.set_fem_attn, kwargs={'ant_str':ant_str})
                    t.daemon = True
                    t.start()
                #==== PLUSDELAY ====
                elif ctlline.split()[0].upper() == '$PLUSDELAY':
                    # Really specific command that adds given delay argument (in nsec) to the Ant 14 LO-FRQ RCVR
                    # delays (both X and Y) in the current delay center table.  This writes the changed record as
                    # a new record in the SQL database.  It can be deleted later using cal_header.delete_cal().
                    try:
                        cmd, dla = ctlline.strip().split()
                        dla = numpy.float(dla)
                        # Read current delay center table, which creates output table in /tmp
                        cal_header.dla_censql2table()
                        time.sleep(0.1)
                        # Read the table from the disk
                        f = open('/tmp/delay_centers.txt','r')
                        lines = f.readlines()
                        f.close()
                        time.sleep(0.1)
                        # Add the given delay to both X and Y delays in line 18 (Ant 15 line, which is really Ant 14 LO_FRQ RCVR)
                        vals = array(map(float,lines[18].strip().split())) + [0,dla,dla]
                        lines[18] = '  {:2d}    {:9.3f}    {:9.3f}\n'.format(int(vals[0]),vals[1],vals[2])
                        # Write the table back to the disk
                        f = open('/tmp/delay_centers.txt','w')
                        for line in lines:
                            f.write(line)
                        f.close()
                        time.sleep(0.1)
                        # Write the table to the SQL database
                        cal_header.dla_centable2sql(filename='/tmp/delay_centers.txt')
                    except:
                        print util.Time.now().iso,'Could not interpret PLUSDELAY arguments.'
                #==== CAPTURE-1S ====
                elif ctlline.split()[0].upper() == '$CAPTURE-1S':
                    # Use $CAPTURE-1S <stem> where <stem> is a string to add to the end of the
                    # capture filename.  The capture is done on the dpp.  This will take a few
                    # seconds to complete.
                    try:
                        cmd, stem = ctlline.strip().split()
                    except:
                        # Must be no stem given, so use ''
                        stem = ''
                    # Capture 1 s of data on dpp
                    t = threading.Thread(target=pcapture2.capture_1s, kwargs={'stem':stem})
                    t.daemon = True
                    t.start()
                elif ctlline.split()[0].upper() == '$WSCRAM-LIMIT':
                    # Update the 27-m windscram limit
                    try:
                        wlimit = int(ctlline.split()[1])
                        self.wlimit = wlimit
                        self.stale = True   # This will force sending of new limit at next 1-sec tick
                    except:
                        pass
                #==== SCAN-START ====
                elif ctlline.split()[0].upper() == '$SCAN-START':
                    # Command by itself is normal scan start, while $SCAN-START NODATA
                    # means set up scan, but do not take data.  Used for some calibrations.
                    nodata = ctlline.strip().split()[-1].upper()
                    # Do any tasks here that are required to start a new scan
                    sys.stdout.write('Started new scan\n')

                    # This block commented out 2020-08-09 due to too-likely failure.  The setting of
                    #   self.lorx is now done on sending an RX-SELECT LO command
                    # # Check for Ant 14 low-frequency receiver status
                    # self.lorx = False   # Default (normal) position is high frequency receiver
                    # # Check Ant 14 Receiver Position Status
                    # data, msg = stateframe.get_stateframe(self.accini)
                    # FEMA = self.accini['sf']['FEMA']
                    # if stateframe.extract(data,FEMA['Timestamp']) != 0:
                        # # This is a valid record, so proceed
                        # if stateframe.extract(data,FEMA['PowerStrip']['RFSwitchStatus']) == 0:
                            # # The switch position is right for LoRX
                            # RX_pos = stateframe.extract(data,FEMA['FRMServo']['RxSelect']['Position']) + stateframe.extract(data,FEMA['FRMServo']['RxSelect']['PositionError'])
                            # if RX_pos < 150.:
                                # # Consistent with LoRX being in position, or heading there, so set as True
                                # self.lorx = True
                                # print 'Ant 14 delays will be set for LO-Frequency Receiver'
                            # else:
                                # print 'Ant 14 outlet set for LO-Frequency Receiver, but RxSelect position is wrong.'
                                # print 'Ant 14 delays will be set for HI-Frequency Receiver.'
                        # else:
                            # print 'Ant 14 delays will be set for HI-Frequency Receiver.'
                    # else:
                        # print 'LO-Frequency Receiver check failed due to bad (0) stateframe.'

                    # ************ This block commented out due to loss of SQL **************
                    xml, buf = cal_header.read_cal(4)
                    try:
                        xml, buf = cal_header.read_cal(4)
                        dcenters = stateframe.extract(buf,xml['Delaycen_ns'])
                        if self.lorx:
                            # If the LO-frequency receiver is active, put delays in slot for Ant 15 into Ant 14
                            dcenters[13] = dcenters[14]
                        timestr = Time(int(stateframe.extract(buf, xml['Timestamp'])), format='lv').iso
                        f = open('/tmp/delay_centers.txt', 'w')
                        f.write('# Antenna delay centers, in nsec, relative to Ant 1\n')
                        f.write('#     Date: ' + timestr + '\n')
                        f.write('# Note: For historical reasons, dppxmp needs four header lines\n')
                        f.write('# Ant  X Delay[ns]  Y Delay[ns]\n')
                        fmt = '{:4d}   {:10.3f}   {:10.3f}\n'
                        for i in range(16):
                            f.write(fmt.format(i + 1, *dcenters[i]))
                        f.close()
                        time.sleep(0.1)  # Make sure file has time to be closed.
                        f = open('/tmp/delay_centers.txt', 'r')
                        acc = FTP('acc.solar.pvt')
                        acc.login('admin', 'observer')
                        acc.cwd('parm')
                        # Send DCM table lines to ACC
                        print acc.storlines('STOR delay_centers.txt', f)
                        f.close()
                        print 'Successfully wrote delay_centers.txt to ACC'
                        
                        sh_dict['dlacen']  = dcenters[:,0]
                        sh_dict['dlaceny'] = dcenters[:,1]
                    except:
                        print util.Time.now().iso,'SQL connection for delay centers failed.  Delay center not updated'
                    # ************* End of block ****************
                    # Replaced by:
                    # delaydict = cal_header.ACCdlatable2dict()
                    # if delaydict == {}:
                        # print util.Time.now().iso,'ACC transfer of delay centers failed.  Delay center not updated'
                    # else:
                        # sh_dict['dlacen']  = delaydict['Delaycen_ns'][:,0]
                        # sh_dict['dlaceny'] = delaydict['Delaycen_ns'][:,1]
                    # Fetch current delay centers from SQL database, and write them to
                    # the ACC file /parm/delay_centers.txt, which is used by the dppxmp program

                    if self.subarray_name == 'Subarray1':
                        try:
                            s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
                            # Send commands to update antenna trip information
                            s.connect((self.accini['host'],self.accini['scdport']))
                            s.send('UPDATEAZIMUTHDIAGNOSTICS 1')
                            time.sleep(0.01)
                            s.send('UPDATEELEVATIONDIAGNOSTICS 1')
                            time.sleep(0.01)
                            s.close()
                        except:
                            pass

                    # We need an initial call to set_uvw() in order to set RA, Dec and HA
                    # coordinates in scan header dictionary.
                    srcname = sh_dict['source_id']
                    try:
                        # Generate a Time() object at exactly the next upcoming second (time t+1)
                        t2 = util.Time.now()
                        tsec = util.Time(t2.mjd  + (1 - t2.datetime.microsecond/1000000.)/86400.,format='mjd')
                        src = self.aa.cat[srcname]        # This causes KeyError if source is not found
                    except KeyError:
                        # The current scan header source ID is not in the source catalog
                        srcname = None
                    sh_dict['timestamp'] = tsec.lv
                    if srcname is not None:
                        set_uvw(self.aa,tsec,srcname)
                        print 'Current RA, Dec, HA:',sh_dict['ra'],sh_dict['dec'],sh_dict['ha']
                        sys.stdout.flush()
                    # Read KATADC status registers.  This can take a long time...
                    sys.stdout.write('There are '+str(len(self.roaches))+' active ROACHes\n')
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
                            
                    # Read acc0time.txt file from ACC and update scan header
                    try:
                        f = urllib2.urlopen('ftp://'+userpass+'acc.solar.pvt/parm/acc0time.txt',timeout=1)
                        mjdacc0 = np.double(f.readline().split()[0])
                        f.close()
                    except:
                        print t.iso,'ACC connection for acc0 time timed out.  Reading from /tmp/acc0time.txt'
                        f = open('/tmp/acc0time.txt','r')
                        mjdacc0 = np.double(f.readline().split()[0])
                        f.close()
                    sh_dict['time_at_acc0'] = Time(mjdacc0,format='mjd')

#                    if self.subarray_name == 'Starburst':
#                        # get a dictionary with any Starburst-specific scan header data and add it to sh_dict
#                        sh_dict_starburst = starburst.get_sh_dict(self,ctlline)
#                        sh_dict.update(sh_dict_starburst)
#                        # make a copy of scan_header file including Starburst-specific data and copy it to Starburst server
#                        starburst.write_scan_header(sh_dict,self.sh_datfile)
#                    else: # write OVSA scan header and store on ACC
                    scan_header(sh_dict,self.sh_datfile)
                    
                    # ************ This block commented out due to loss of SQL **************
                    # If we are connected to the SQL database, send converted scan header
                    if self.sql['cnxn']:
                        f = open(self.sh_datfile)
                        data = f.read()
                        f.close()
                        bufout = stateframedef.transmogrify(data, self.sql['shbrange'])
                        try:
                            self.sql['cursor'].execute('insert into hBin (Bin) values (?)', 
                        stateframedef.pyodbc.Binary(bufout))
                            self.sql['cnxn'].commit()
                            sys.stdout.write('Scan Header Record successfully written to SQL Server\n')
                            sys.stdout.flush()
                        except:
                            # An exception could be an error, or just that the entry was already inserted
                            sys.stdout.write('Writing Scan Header record to SQL Server FAILED\n')
                            sys.stdout.flush()
                            self.error = 'Err: Cannot write scan header to SQL'
                    # ************* End of block ***************
                    # Also write scan header data to log file
                    f2 = open(self.sh_datfile)
                    data = f2.read()
                    f2.close()
                    f = self.accini.get('sh_file')
                    try:
                        f.write(data)
                    except:
                        print Time.now().iso+' Error writing scan header to log file'

                    
                    if nodata == 'NODATA':
                        pass
                    else:
                        # Set scan state to on
                        sf_dict['scan_state'] = 1
                #==== SCAN-RESTART ====
                elif ctlline.split()[0].upper() == '$SCAN-RESTART':
                    # This command is for restarting a scan with the same setup as the
                    # previously running scan, where only the scan state must be turned on
                    # Set scan state to on
                    sf_dict['scan_state'] = 1
                #==== SCAN-STOP ====
                elif ctlline.split()[0].upper() == '$SCAN-STOP':
                    sf_dict['scan_state'] = -1
                #==== PA-SWEEP ====
                elif ctlline.split()[0].upper() == '$PA-SWEEP':
                    # Rotate 27-m focus rotation mechanism to sweep through a given angle 
                    # range (does nothing if PA adjustment routine is already running)
                    #   Usage: $PA-SWEEP PA rate, where angle is swept from -PA to PA at
                    #                             rate of 1-degree-per-rate [s]
                    print 'Got '+ctlline.split()[0].upper()+' command.'
                    if self.PAthread is None or not self.PAthread.is_alive():
                        # Thread is not already running, so it is safe to proceed
                        try:
                            PA,rate = map(numpy.int,ctlline.strip().split()[-2:])
                        except:
                            # Reading arguments failed, so use defaults
                            PA = 80
                            rate = 3
                        print 'PA and rate are ',PA,rate
                        try:
                            # Spawn the stateframe.PA_sweep() routine to update PA once/rate
                            self.PAthread = threading.Thread(target=stateframe.PA_sweep,kwargs={'PA':PA,'rate':rate})
                            self.PAthread.daemon = True
                            self.PAthread.start()
                            print 'PAthread started.'
                        except:
                            # Something went wrong
                            print 'Error spawning PA_sweep task'
                            pass
                #==== PA-TRACK ====
                elif ctlline.split()[0].upper() == '$PA-TRACK':
                    # Track 27-m focus rotation mechanism to correct for parallactic angle 
                    # of given antenna (does nothing if antenna not correctly specified)
                    #   Usage: $PA-TRACK ant4 <CROSSED>
                    print 'Got '+ctlline.split()[0].upper()+' command.'
                    if self.PAthread is None or not self.PAthread.is_alive():
                        # Thread is not already running, so it is safe to proceed
                        antstr = ctlline.strip().split()[1].upper()
                        crossed = False
                        if len(ctlline.strip().split()) == 3:
                            if ctlline.strip().split()[2].upper() == 'CROSSED':
                                crossed = True
                        print 'Given antenna is '+antstr
                        try:
                            # Spawn the stateframe.PA_adjust() routine to update PA once/minute
                            antn = pcapture2.ant_str2list(antstr)[0]
                            print 'Antenna index is',antn
                            self.PAthread = threading.Thread(target=stateframe.PA_adjust,kwargs={'ant':antn,'crossed':crossed})
                            self.PAthread.daemon = True
                            self.PAthread.start()
                            print 'PAthread started.'
                        except:
                            # Antenna not correctly specified, so do not spawn routine
                            print 'Antenna specification no good?'
                            pass
                #==== PA-STOP ====
                elif ctlline.split()[0].upper() == '$PA-STOP':
                    # Send Abort string to stateframe.PA_adjust() routine.  Note that
                    # abort may not be acted upon until up to 1 s later.
                    if self.PAthread and self.PAthread.is_alive():
                        stateframe.q.put_nowait('Abort')
                #==== TRIPS ====
                elif ctlline.split()[0].upper() == '$TRIPS':
                    try:
                        # Send commands to update antenna trip information
                        s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
                        s.connect((self.accini['host'],self.accini['scdport']))
                        s.send('UPDATEAZIMUTHDIAGNOSTICS 1')
                        time.sleep(0.01)
                        s.send('UPDATEELEVATIONDIAGNOSTICS 1')
                        time.sleep(0.01)
                        s.close()
                    except:
                        pass
                #==== DLASWEEP ====
                elif ctlline.split()[0].upper() == '$DLASWEEP':
                    try:
                        vals = ctlline.split()
                        print vals
                        if len(vals) == 4:
                            junk,ant,dla,dlastop = vals
                            pol = None
                        elif len(vals) == 5:
                            junk,ant,dla,dlastop,pol = vals
                        for r in self.roaches:
                            a1,a2 = r.ants
                            if int(ant) == a1 or int(ant) == a2 or int(ant) == 0:
                                r.dlasweep = {'ant':int(ant),'dla':int(dla),'dlastop':int(dlastop),'pol':pol}
                    except:
                        print 'Could not interpret $DLASWEEP command'
                #==== WAIT ====
                elif ctlline.split()[0].upper() == '$WAIT':
                    # Need to wait for given number of seconds, so set self.waitmode to True,
                    # set self.nextctlline to point to next following line, and record duration
                    try:
                        dur = int(ctlline.split()[1])
                    except:
                        print 'Could not interpret duration on $WAIT command--defaulting to 10 s'
                        dur = 10
                    self.waitmode = True
                    #self.nextctlline = i+1
                    self.wait = dur
                    print 'Initializing wait for',dur,'seconds'
                    #break
                #==== PCYCLE ====
                elif ctlline.split()[0].upper() == '$PCYCLE':
                    # Cycle the power of some device (antenna controller, fem, or crio)
                    # for a given antenna.
                    device, ants = self.interpret_pcycle(ctlline)
                    if device is None:
                        print 'Error interpreting $PCYCLE command',ctlline
                    else:
                        # Since device is not None, interpreting tokens succeeded.
                        if device == 'ANT':
                            # Turn off all antennas that are to be power cycled
                            antstr = ''
                            for antnum in ants:
                                antstr += ' ant'+str(antnum)
                            self.sendctlline('powerswitch 0 '+antstr)
                        # Spawn tasks to perform power cycle (each takes awhile)
                        for antnum in ants:
                            t = threading.Thread(target=pwr_cycle.ant_toggle, args=(antnum, device))
                            t.daemon = True
                            t.start()
                #==== KATADC_GET ====
                elif ctlline.split()[0].upper() == '$KATADC_GET':
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
                elif ctlline.split()[0].upper() == '$REWIND':
                    # Get date of first line of current schedule
                    # and increment by 1 day, then autogenerate a
                    # new schedule
                    self.toggle_state()  # Turn off schedule
                    # Get time of first line in schedule and add a day
                    mjd1 = mjd(self.L.get(0))
                    t = util.Time(mjd1+1,format='mjd')
                    self.Today(t)
                    self.Clear()
                    self.toggle_state()  # Turn schedule back on
                #==== LNA_INIT ====
                elif ctlline.split()[0].upper() == '$LNA-INIT':
                    # Get LNA_settings.txt file from ACC and send the series of
                    # commands needed to set the LNA voltages
                    #userpass = 'admin:observer@'
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
                            self.sendctlline(cmdstr)
                            cmdstr = 'LNA-DRAIN '+lnas[i]+' '+str(lnas_a[i]['vd'])+' ANT14'
                            self.sendctlline(cmdstr)
                            cmdstr = 'LNA-GATE1 '+lnas[i]+' '+str(lnas_a[i]['vg1'])+' ANT14'
                            self.sendctlline(cmdstr)
                            cmdstr = 'LNA-GATE2 '+lnas[i]+' '+str(lnas_a[i]['vg2'])+' ANT14'
                            self.sendctlline(cmdstr)
                            cmdstr = 'LNA-ENABLE '+lnas[i]+' on ANT15'
                            self.sendctlline(cmdstr)
                            cmdstr = 'LNA-DRAIN '+lnas[i]+' '+str(lnas_b[i]['vd'])+' ANT15'
                            self.sendctlline(cmdstr)
                            cmdstr = 'LNA-GATE1 '+lnas[i]+' '+str(lnas_b[i]['vg1'])+' ANT15'
                            self.sendctlline(cmdstr)
                            cmdstr = 'LNA-GATE2 '+lnas[i]+' '+str(lnas_b[i]['vg2'])+' ANT15'
                            self.sendctlline(cmdstr)
                    except:
                        print 'Error sending LNA_settings to ACC'

                #==== SUBARRAY ====
                elif ctlline.split()[0].upper() == '$SUBARRAY':
                    # run the SUBARRRAY1 command if this is the master schedule,
                    # otherwise run the SUBARRAY2 command
                    print '$SUBARRAY line is:',ctlline
                    if ctlline.find('.antlist') == -1:
                        # there is no .antlist file in this line --> the antlist
                        # should be directly specified in the line, e.g.:
                        # $SUBARRAY ant1 ant7-8,ant10
                        l = len('$SUBARRAY ')
                        antlist = ctlline[l:]
                    else:
                        # a .antlist file is specified --> read antlist from
                        # the specified .antlist file
                        antlistfile = ctlline.split()[1]
                        try:
                            antlistname = ctlline.split()[2]
                            antlist = get_antlist(antlistname,antlistfile)
                        except:
                            antlist = ''
                        if antlist == '':
                            self.error = '$SUBARRAY: antlist name not in ' + antlistfile
                            return
                    if self.subarray_name == 'Subarray1':
                        N = 1
                    else:
                        N = 2
                    cmd = 'SUBARRAY' + str(N) + ' ' + antlist
                    try: # send appropriate command to ACC
                        s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
                        s.connect((self.accini['host'],self.accini['scdport']))
                        s.send(cmd)
                        time.sleep(0.01)
                        s.close()
                        print 'ctl cmd \'' + ctlline + '\' sent to ACC as \'' + cmd + '\''
                        sys.stdout.flush()
                    except:
                        print 'ctl cmd \'' + cmd + '\' not succesfully sent to ACC'
                        sys.stdout.flush()
                        pass
            else:
                cmds = ctlline.split()
                if cmds[0].upper() == 'FSEQ-FILE':
                    # This is an FSEQ-FILE command, so find and set frequency sequence
                    # First set up channel info for this sequence
                    sh_dict['chinfo'].fseq2nsavg(cmds[1])
                    # Then FTP sequence file from ACC
                    fseqfile = urllib2.urlopen('ftp://'+userpass+'acc.solar.pvt/parm/'+cmds[1])
                    nrpt = None     # Initially not defined
                    fsequence = ''  # Initially empty
                    for line in fseqfile.readlines():
                        # Find DWELL line (contains 35 dwell times, in s, one for each
                        # defined band) and split into 35 "repeat" numbers
                        keywd = 'LIST:DWELL'
                        if line.find(keywd) == 0:
                            dwellseq = line[len(keywd):].split(',')
#                            if len(dwellseq) != 35:
#                                print 'FSEQ file',cmds[1],'DWELL line must have 35 entries.'
#                                break
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
                                for i in range(nrpt[band-1]):
                                    fsequence += str(band)+','
                            break
                    if fsequence == '':
                        print 'FSEQ file',cmds[1],'not successfully interpreted.'
                        # Default to allowing all channels in RFI mask
                        sh_dict.update({'chanmask':numpy.array([1]*204800,'byte')})
                    else:
                        sh_dict.update({'fsequence':fsequence[:-1]})   # -1 removes trailing ','
                        chanmask = ci.get_chanmask(fsequence[:-1])
                        sh_dict.update({'chanmask': chanmask})
                        # Get nominal Chan2Wide assignment, then multiply by chanmask and update it
                        fseqlist = fsequence[:-1].rsplit(',')
                        item = []
                        for band in fseqlist:
                            ch = sh_dict['chinfo'].chan_asmt(int(band))
                            item += ch
                        sh_dict.update({'chan2wide':numpy.array(item)*sh_dict['chanmask']})
                        self.sequence2roach(fsequence[:-1])
                        self.sequence2dcmtable(fsequence[:-1])
                elif cmds[0].upper() == 'RX-SELECT':
                    if cmds[1].upper() == 'LO':
                        self.lorx = True # Next $SCAN-START will use low-frequency receiver delays for Ant14
                    else:
                        self.lorx = False # Next $SCAN-START will use high-frequency receiver delays for Ant14

    def update_status(self):
        # Read the current schedule from the list window and write to output file
        fileout = open('/common/webplots/status.txt','w')
        for i in range(self.lastline):
            line = self.L.get(i)
            if i == self.curline:
                line = '* '+line
            else:
                line = '  '+line
            fileout.write(line+'\n')
        fileout.close()


app = App()

mainloop()
