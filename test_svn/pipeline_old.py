# Note mirexec, numpy are found in /usr/lib/python2.7/dist-packages/
# Added idbfits procedure in udb_process, jmm, 23-oct-2015
from numpy import *
from mirexec import *
import time
import fdb
# for deleting old datasets
import os
import shutil
# for IDB fits processing
import dump_tsys_ext
from util import Time
import calibration

# defines a varplt class for miriad python:
class TaskVarplt(TaskBase): 
    _keywords = ['vis','device','log','xaxis','yaxis','nxy','xrange','yrange'] 
    _options = ['dtime','compress','overlay','unwrap']

global dppdatadir
dppdatadir = '/dppdata1/'
#Interim database files
global idbdir
idbdir = '/dppdata1/IDB/'
#FDB files for IDB
global fdbdir
fdbdir = '/dppdata1/FDB/'
#Utility data files
global udbdir
udbdir = '/data1/UDB/'
#IFDB files have 1 entry for each of the IDB files
global ifdbdir
ifdbdir = '/data1/IFDB/'
#UFDB files have 1 entry for each of the UDB files
global ufdbdir
ufdbdir = '/data1/UFDB/'
#Text files to be processed into fits files
global udbtxtdir
udbtxtdir = '/data1/UDBTXT/'
global idbfinaldir
idbfinaldir = '/data1/eovsa/fits/IDB/'

def udb_init_fdbfiles(ndays=7,fdbtype='FDB'):
    '''Reads in the available FDB files, for the last N days'''
    # check for what type of filedb you're looking for
    # currently FDB, IFDB
    if fdbtype == 'FDB':
        xdir = fdbdir
    elif fdbtype == 'IFDB':
        xdir = ifdbdir
    elif fdbtype == 'UFDB':
        xdir = ufdbdir
    else:
        xdir = fdbdir

    #local time in seconds from 1970
    ptime = time.mktime(time.localtime())
    #create a list of FDB files, today, yesterday, the day before, etc....
    fdbfiles = []
    for j in range(ndays):
    #j days ago
        p3 = ptime-j*86400.0
    #change format for input into string time
        p3tuple = time.localtime(p3)
        ppp = time.strftime('%Y%m%d', p3tuple)
        if fdbtype == 'FDB':
            fdbfiles = fdbfiles+[xdir+fdbtype+ppp+'.txt']
        else:
            yyyy = time.strftime('%Y', p3tuple)+'/'
            fdbfiles = fdbfiles+[xdir+yyyy+fdbtype+ppp+'.txt']

    #reverse order
    fdbfiles.reverse()
    return fdbfiles
#End of udb_init_fdbfiles

def udb_fb2process(ndays=7):
    '''Reads in the last N days of FDB files, and IFDB files, finds
    the IDB files that need processing, and returns the fdb ifdb
    classes for those files'''

    #get FDB, IFDB filenames
    fdbfiles = udb_init_fdbfiles(ndays=ndays)
    if len(fdbfiles) == 0:
        print 'UDB_FB2PROCESS Error: No FDB files:'
        ifboops = []
        fboops = []
        return fboops, ifboops

    ifdbfiles = udb_init_fdbfiles(ndays=ndays,fdbtype='IFDB')
    if len(ifdbfiles) == 0:
        # Not necessarily an error
        print 'UDB_FB2PROCESS: No IFDB files:'

    #read in the files
    fb=[]
    for j in range(len(fdbfiles)):
        fbj = fdb.fdb_read(fdbfiles[j])
        if len(fbj) > 0:
            fb = fb+fbj

    ifb=[]
    if len(ifdbfiles) > 0:
        for j in range(len(ifdbfiles)):
            ifbj = fdb.pfdb_read(ifdbfiles[j])
            if len(ifbj) > 0:
                ifb = ifb+ifbj

    #check for non-empty fb's first
    if len(fb) == 0:
        print "no data to process"
        otp = []
        return fb, otp

    #if I am here, I have IDB files from the FDB. Do I have PFDB entries?
    if len(ifb) == 0:
        print "no IFDB entries, process all FDB entries"
        #For each FDB entry, create a IFDB entry, with pstatus = 0
        for j in range(len(fb)):
            ifbj = fdb.pfiledb(fb[j].fileid,fb[j].scanid,fb[j].sourceid,fb[j].projectid,fb[j].st_ts,fb[j].en_ts,0)
            ifb.append(ifbj)

        return fb, ifb

    #if I am here, I have both FDB and IFDB entries, for each FDB,
    #check to see if there is a corresponding IFDB. If not, create one
    #with pstatus = 0
    iflist = fdb.fdb_list_fileid(ifb)
    iidarray = array(iflist)
    ifb_out = []
    for j in range(len(fb)):
        fileidj = fb[j].fileid
        ifbtemp = extract(iidarray == fileidj, ifb)
        if(len(ifbtemp) > 0):
            ifb_outj = ifbtemp[0]
        else:
            ifb_outj = fdb.pfiledb(fb[j].fileid,fb[j].scanid,fb[j].sourceid,fb[j].projectid,fb[j].st_ts,fb[j].en_ts,0)
            
        ifb_out.append(ifb_outj)

    return fb, ifb_out
#End of udb_fb2process

def udb_scan2process(ifb):
    '''Given a set of IFDB entries, pick out the files that need
    processing, the first scan that has any file with pstatus = 0 is to
    be processed'''
    #check the input first
    ifb_out = []
    ss_ifb_out = []
    if len(ifb) == 0:
        print "No ifb entries input, No processing"
        return ifb_out, ss_ifb_out
    
    #here we will go scan by scan, first get the unique scan ids'
    scan_ids = fdb.fdb_uniq_scan(ifb)
    nscan = len(scan_ids)
    #for each unique scan id, pick out the ifb entries
    redo_scan = 0
    for j in range(nscan):
        if redo_scan == 0:
            ifb_test, ss_ifb_test = fdb.fdb_extract_1scan(ifb, scan_ids[j])
            #check for partially completed scans
            for i in range(len(ifb_test)):
                if ifb_test[i].pstatus==0:
                    redo_scan = 1
                #endif 
            #endfor 
            #if I have a partially completed scan,
            #keep non-zero pstatuses, jmm, 2103-10-23
            if redo_scan == 1:
                for i in range(len(ifb_test)):
                    #ifb_test[i].pstatus=0
                    #if we have a scan to redo, then populate ifb_out,
                    #which os a list, while ifb_test is an array
                    ifb_out.append(ifb_test[i])
                #endfor 
                ss_ifb_out = ss_ifb_test.tolist()
            #endif 
        #endif
    #endfor
    return ifb_out, ss_ifb_out
# end of udb_scan2process
        
def udb_process1scan(ifb_1):
    '''Given the set of ifb entries, process
    them into a single UDB Miriad dataset for this scan'''

        
    # defines a varplt class for miriad python: class
    # TaskVarplt(TaskBase): _keywords =
    # ['vis','device','log','xaxis','yaxis','nxy','xrange','yrange']
    # _options = ['dtime','compress','overlay','unwrap']

    # a call to the varplt task TaskVarplt(vis='temp_pytest',
    # xaxis='ut', yaxis='ytsys', device='/xs',
    # xrange='0.975,0.977',options='overlay').run()
    
    # Uvaver is the class we'll call here: vis is a set of files, out
    # is the output file, line='channel,50,1,10,10' selects 50
    # channels, starting with channel 1, summing over 10 channel, then
    # skipping to every 10th channel, interval = '0.666667' specifies
    # 4 second intervals

    # TaskUVAver(vis='/dppdata1/IDB/IDB20131220232416,/dppdata1/IDB/IDB20131220232516',
    # out='temp_pytest1',
    # line='channel,50,1,10,10',interval='0.066667').run()
    # switched to 1 second intervals, jmm, 2014-06-17
    #first get the list of input files
    n1 = len(ifb_1)
    filelist = fdb.fdb_list_fileid(ifb_1)
    print filelist
    #Next, name the output files, this should be based on the start
    #time of the first file, changed from using scan_id, 2014-06-17,
    #jmm
    scan_id = ifb_1[0].scanid
    source_id = ifb_1[0].sourceid
    project_id = ifb_1[0].projectid
    st_ts = ifb_1[0].st_ts
    en_ts = ifb_1[n1-1].en_ts

    #Take care of two digit year here
    #    ufileid = 'UDB'+'20'+scan_id
    ufileid = 'U'+ifb_1[0].fileid[1:]
    yyyy = ufileid[3:7]
    ufilename = udbdir+yyyy+'/'+ufileid
    ufb = fdb.pfiledb(ufileid, scan_id, source_id, project_id, st_ts, en_ts, 0)

    #call taskuvaver on the file list, which looks like it needs to be
    #a comma-separated string
    filelist_str = ''
    filelist_str = idbdir+filelist[0]
    if n1 > 1:
        #Not just for testing anymore, something will crash if there
        #are too many files, it looks like the most files that can be
        #processed might be 63, but we'll stick with 60 -- jmm,
        #2014-08-23
        if n1 > 60:
            print "UDB_PROCESS1SCAN: Too many Files N1 = ", n1
            print "UDB_PROCESS1SCAN: Reset to 60"
            n1 = 60
        #end of not just for testing anymore
        for j in range(n1-1):
            filelist_str = filelist_str+','+idbdir+filelist[j+1]

    #If the file exists, you need to delete it before uvavering
    if os.path.isdir(ufilename) == True:
        print "UDB_PROCESS1SCAN: dataset: ", ufilename, " will be deleted"
        shutil.rmtree(ufilename)

    #Dropped 4 second averaging, 2014-06-17, jmm
    TaskUVAver(vis=filelist_str,out=ufilename,line='channel,50,1,10,10').run()

    return [ufb]
# End of udb_process1scan

def udb_process(ndays=7):
    ''' This is the main program, that reads the FDB files, decides
    what to do, if anything, prcesses a scan, updates the DB and
    finishes'''

    #initialize output
    ufb_1 = []
    #first read the FDB files from DPP, and determine if there are any
    #files that need to be processed. Here fb are from the DPP, and
    #will not be touched. The ifb are the ones that we will work with
    #and should have the same number of elements as the fb; both are
    #lists'''
    fb, ifb = udb_fb2process(ndays=ndays)
    #Check for success
    if len(ifb) == 0:
        print 'UDB_PROCESS Error: No filedbs input, No processing'
    else:
        #Now process the first scan with any files with pstatus = 0
        ifb_1, ss_ifb_1 = udb_scan2process(ifb)
        #Only process if you have to
        if len(ifb_1) == 0:
            print 'UDB_PROCESS: No Scans to process, finished successfully'
        else:
            print 'UDB_PROCESS: Processing Scan: '+ifb_1[0].scanid
            ufb_1 = udb_process1scan(ifb_1)
            #Full database management
            if len(ufb_1) == 0:
                print 'UDB_PROCESS Error: No output from UDB_PROCESS1SCAN'
            else:
                #Update ifdb.pstatus for the processed files, output
                #idbfits and copy IDB files if necessary
                n1 = len(ifb_1)
                for j in range(n1):
                    if ifb_1[j].pstatus == 0:
                        print 'UDB_PROCESS: Processing IDB: '+ifb_1[j].fileid
                        fj, info = process_1ifb(ifb_1[j])
                    #end if
                    ifb_1[j].pstatus = 1
                    ifb[ss_ifb_1[j]].pstatus = 1
#                    print 'UDB_PROCESS: Reset status for:'+ifb[ss_ifb_1[j]].fileid
#                    print 'UDB_PROCESS: Reset status for:'+ifb_1[j].fileid
                #End j loop
                #Write out the new ifb entries
                ifb_files_out = update_ifdb(ifb)

                #Write out ufdb entry, ufb_1 is a 1 element list
                ufb_files_out = update_ufdb(ufb_1, ndays=ndays)
                        
                print 'UDB_PROCESS: Processed Scan: '+ifb_1[0].scanid
                print 'UDB_PROCESS: finished successfully'

            #Endelse
        #Endelse
    #Endelse
    return ufb_1
# End of udb_process

def process_1ifb(ifb0):
    '''Creates an IDB fits file for the input ifb entry, and copies
    the IDB Miriad file to the local disk. Added antennalist,
    source_id and project_id, output 2015-01-06, jmm'''

    #string may be a list
    if isinstance(ifb0, list):
        ifb = ifb0[0]
    else:
        ifb = ifb0
    #end else
    filelist = [idbdir+ifb.fileid]
    #Check validity
    file_test = dump_tsys_ext.valid_miriad_dataset(filelist)
    if len(file_test) == 0 or file_test[0] == False :
        print "process_1ifb: Bad Miriad Dataset:"
        print filelist[0]
        return '', ''
    #endif
    #If you are here you have a good IDB dataset, error check anyway
    xdat = dump_tsys_ext.rd_miriad_tsys_file(filelist)
    if xdat == None:
        print "process_1ifb: Bad Miriad Dataset (2):"
        print filelist[0]
        return '', ''
    #endif
    #Hold onto antennalist for output
    antlist = xdat['antennalist']
    #Calibrate, but only for normal observing
    sid = ifb.sourceid[:3]
    prid = ifb.projectid[:15]
    if sid == "Sun" and prid == "NormalObserving":
        t1 = Time(xdat['ut_mjd'][0],format='mjd')
        fghz, calfac, offsun = calibration.sp_read_calfac(t1)
        if fghz != None:
            xdat = calibration.sp_apply_cal(xdat, fghz, calfac, offsun)
            calflag = True
        else:
            calflag = False
        #endif
    else:
        calflag = False
    #endif
    #write the file
    print "CALFLAG", calflag
    file_out = dump_tsys_ext.tsys_writetofits(xdat, calflag)
    if len(file_out) == 0:
        print "process_1ifb: write failed:"
        print filelist[0]
        return file_out, ''
    #endif
    #copy the IDB dataset from /dppdata1/IDB to /data1/eovsa/fits/IDB/yyyymmdd
    filename = filelist[0]
    yyyymmdd = filename[len(filename)-14:len(filename)-6]
    outdir = idbfinaldir+yyyymmdd
    if os.path.isdir(outdir) != True:
        os.mkdir(outdir)
    #end if
    full_filename = outdir+'/'+ifb.fileid
    #If the file exists, you need to delete it before using shutil.copytree
    if os.path.isdir(full_filename) == True:
        print "process_1ifb: dataset: ", full_filename, " will be overwritten"
        shutil.rmtree(full_filename)
    #end if
    shutil.copytree(filename, full_filename)
#the second output is passed into udb_fitsfile
    return file_out, antlist
#End of process_1ifb

def udb_process_spec(ndays=7):
    '''This process reads the UDB filedb and creates spectral FITS
    files for each scan as necessary '''

    #Read in UFDB files and check status.
    ufdbfiles = udb_init_fdbfiles(ndays=ndays,fdbtype='UFDB')
    if len(ufdbfiles) == 0:
        print 'UDB_PROCESS_SPEC Error: No UFDB files:'
        ufboops = []
        return ufboops

    ufb = []
    for j in range(len(ufdbfiles)):
        ufbj = fdb.pfdb_read(ufdbfiles[j])
        if len(ufbj) > 0:
            ufb = ufb+ufbj

    #check for non-empty ufb's
    if len(ufb) == 0:
        print "UDB_PROCESS_SPEC: no data to process"
        return ufb

    #If I am here, I have a valid ufb, but do I have status = 0, If
    #so, get the first one and process it.
    do_this_file = -1
    for j in range(len(ufb)):
        if do_this_file == -1:
            if ufb[j].pstatus == 0:
                do_this_file = j

    #the variable do_this_file now contains the first filename, unless
    #it's -1, then no processing happens
    if do_this_file == -1:
        print "UDB_PROCESS_SPEC: no data to process"
        print 'UDB_PROCESS_SPEC: finished successfully'
        ufb_out = []
        idlfile = "/home/user/workdir/udb_process_spec.pro"
        f = open(idlfile, 'w')
        s = "pr_path_lib, 'eovsa_write_pwrfits', /multi\n"
        f.write(s)
        s = "message, /info, 'No data to process'\n"
        f.write(s)
        s = "end\n"
        f.write(s)
        f.close
        return ufb_out

    ufileid = ufb[do_this_file].fileid
    uprojectid = ufb[do_this_file].projectid
    yyyy = ufileid[3:7]
    print 'Processing: ', ufileid
    # Check for YYYY directories
    if os.path.isdir(udbdir+yyyy) != True:
        os.mkdir(udbdir+yyyy)
    if os.path.isdir(udbtxtdir+yyyy) != True:
        os.mkdir(udbtxtdir+yyyy)
    ufilename = udbdir+yyyy+'/'+ufileid
    # Write text files for xtsys, ytsys, sfreq, sdf  and header for this file
    xtsysfile = udbtxtdir+yyyy+'/'+ufileid+'_xtsys.txt'
    ytsysfile = udbtxtdir+yyyy+'/'+ufileid+'_ytsys.txt'
    sfreqfile = udbtxtdir+yyyy+'/'+ufileid+'_sfreq.txt'
    sdffile = udbtxtdir+yyyy+'/'+ufileid+'_sdf.txt'
    TaskVarplt(vis=ufilename,xaxis='time',yaxis='xtsys',log=xtsysfile).run()
    TaskVarplt(vis=ufilename,xaxis='time',yaxis='ytsys',log=ytsysfile).run()
    TaskVarplt(vis=ufilename,xaxis='time',yaxis='sfreq',log=sfreqfile).run()
    TaskVarplt(vis=ufilename,xaxis='time',yaxis='sdf',log=sdffile).run()
    hdrfile = udbtxtdir+yyyy+'/'+ufileid+'_hdr.txt'
    TaskPrintHead(in_=ufilename,log=hdrfile).run()

    # Here figure out how to call the IDL process from here, write a
    # batch file and run ssw_batch using the subprocess command
    # 2014-06-26. Instead, just create the batch file and we will run
    # ssw_batch directly from the shell script

    # 2014-12-15 added projectid to eovsa_write_pwrfits.pro call
 
    idlfile = "/home/user/workdir/udb_process_spec.pro"
    f = open(idlfile, 'w')
    s = "pr_path_lib, 'eovsa_write_pwrfits', /multi\n"
    f.write(s)
    s = "fileid = '"+ufileid+"'\n"
    f.write(s)
    s = "hdrfile = '"+hdrfile+"'\n"
    f.write(s)
    s = "sfreqfile = '"+sfreqfile+"'\n"
    f.write(s)
    s = "sdffile = '"+sdffile+"'\n"
    f.write(s)
    s = "xtsysfile = '"+xtsysfile+"'\n"
    f.write(s)
    s = "ytsysfile = '"+ytsysfile+"'\n"
    f.write(s)
    s = "project_id = '"+uprojectid+"'\n"
    f.write(s)
    s = "eovsa_write_pwrfits, fileid, hdrfile, sfreqfile, sdffile, xtsysfile, ytsysfile, outfil=outfil, project_id = project_id\n"
    f.write(s)
    s = "eovsa_plot_tpower_fits, outfil\n"
    f.write(s)
    s = "end\n"
    f.write(s)
    f.close()

#subprocess is not dealing correctly with ssw_batch arguments, so write a shell script to spawn
#    shellfile = "/home/user/workdir/udb_process_specpro.csh"
#    f = open(shellfile, 'w')
#    f.write("#! /bin/csh\n")
#    f.write(" \n")
#    f.write("ssw_batch "+idlfile+"\n")
#    f.close

#hard coded, BECAUSE PYTHON SUCKS
#    subprocess.call(["/bin/csh", "/home/user/workdir/udb_process_specpro.sh"])

#    os.spawnl("ssw_batch udb_process_spec.pro")
#    subprocess.call(["ssw_batch", "udb_process_spec.pro"])
#    subprocess.call(["ssw_batch", "/home/user/workdir/udb_process_spec.pro"])

#    subprocess.call(["rm", "-f", idlfile])


    ufb[do_this_file].pstatus = 1
    ufb_out = [ufb[do_this_file]]
    #Write out ufdb entry, ufb_out is a 1 element list
    ufb_files_out = update_ufdb(ufb_out, ndays=ndays)
    print 'UDB_PROCESS_SPEC: Processed Spectra for: '+ufileid
    print 'UDB_PROCESS_SPEC: finished successfully'

    return ufb_out
# End of udb_process_spec

def update_ifdb(ifb):
    ''' Updates ifdb files, Here we just write out the files for the
    previous ndays days'''

    files_out = []
    #Get days for each ifb entry
    ifdays = fdb.fdb_list_fileid(ifb)
    # 'FDByyyymmddhhmmss' extract subscript range 3 to 11 to get the date
    # 'IFDByyyymmddhhmmss' extract subscript range 4 to 12 to get the date
    for j in range(len(ifdays)):
        y = str.find(ifdays[j], 'DB')
        x0 = y+2
        x1 = y+10
        ifdays[j] = ifdays[j][x0:x1]

    #Get the filenames for the files,
    days = fdb.fdb_uniq_day(ifb)
    ndays_out = len(days)
    if ndays_out == 0:
        print 'UPDATE_IFDB Error: No days to process'
    else:
        #need an array of days for extract commands
        idyarr = array(ifdays)
        for j in range(ndays_out):
            #Find files with dates matching each day here, build a
            #list to output
            ifbj = extract(idyarr==days[j], ifb)
            ifbj_out = ifbj.tolist()
            if len(ifbj_out)==0:
                print 'UPDATE_IFDB Error: No output for: '+days[j]
            else:
                #filename
                yyyy = days[j]
                yyyy = yyyy[:4]
                filenamej = ifdbdir+'/'+yyyy+'/'+'IFDB'+days[j]+'.txt'
                #check for yyyy directory
                if os.path.isdir(ifdbdir+yyyy) != True:
                    os.mkdir(ifdbdir+yyyy)
                fdb.pfdb_write(ifbj_out, filenamej)
                print 'Wrote: '+filenamej
                files_out = files_out+[filenamej]

    return files_out
    
# End of update_ifdb

def update_ufdb(ufb, ndays=7):
    '''Updates ufdb files, need to read the last ndays files, check to
    see if you're overwriting an entry, if not then append the new
    scan to the list and output. If so, replace the old entry, then
    output. Note that ufb is a 1-element list'''

    if len(ufb) == 0:
        print 'UPDATE_UFDB ERROR: No UFB input'

    files_out = []
    #get filenames
    ufdb_files = udb_init_fdbfiles(ndays=ndays,fdbtype='UFDB')
    
    ufb_in = []
    if len(ufdb_files) == 0:
        #maybe there are none
        print 'UPDATE_UFDB: No UFDB files in range:'
    else:
        for j in range(len(ufdb_files)):
            ufbj = fdb.pfdb_read(ufdb_files[j])
            if len(ufbj) > 0:
                ufb_in = ufb_in+ufbj

    # If there is no ufb_in, then it's easy
    if len(ufb_in) == 0:
        ufb_out = ufb
    else:
        # Check to see if there is an entry for the current ufb
        uflist_in = fdb.fdb_list_fileid(ufb_in)
        ufarray = array(uflist_in)
        # return the subscript
        ss_ufb_in = extract(ufarray==ufb[0].fileid, array(range(len(ufb_in))))
        #I think ss_ufb_in should be a scalar..
        if len(ss_ufb_in) == 0:
            ufb_out = ufb_in+ufb
        else:
            # Get all of the entries that don't have this fileid
            ss_ufb_not = extract(ufarray!=ufb[0].fileid, array(range(len(ufb_in))))
            if len(ss_ufb_not) == 0:
                ufb_out = ufb
            else:
                ufb_out = []
                for j in range(len(ss_ufb_not)):
                    ufb_out = ufb_out+[ufb_in[ss_ufb_not[j]]]
                ufb_out = ufb_out+ufb
        
    #Here, now output the ufb_out entries
    files_out = []
    #Get days for each ufb entry
    ufdays = fdb.fdb_list_fileid(ufb_out)
    for j in range(len(ufdays)):
        y = str.find(ufdays[j], 'DB')
        x0 = y+2
        x1 = y+10
        ufdays[j] = ufdays[j][x0:x1]

    #Get the filenames for the files,
    days_out = fdb.fdb_uniq_day(ufb_out)
    #Adjustment due to lack 
    ndays_out = len(days_out)
    if ndays_out == 0:
        print 'UPDATE_UFDB Error: No days to process'
    else:
        #need an array of days for extract commands
        udyarr = array(ufdays)
        for j in range(ndays_out):
            #Find files with dates matching each day here, build a
            #list to output
            ufbj = extract(udyarr==days_out[j], ufb_out)
            ufbj_out = ufbj.tolist()
            if len(ufbj_out)==0:
                print 'UPDATE_UFDB Error: No output for: '+days[j]
            else:
                #filename
                yyyy = days_out[j]
                yyyy = yyyy[:4]
                filenamej = ufdbdir+'/'+yyyy+'/'+'UFDB'+days_out[j]+'.txt'
                #check for yyyy
                if os.path.isdir(ufdbdir+yyyy) != True:
                    os.mkdir(ufdbdir+yyyy)
                fdb.pfdb_write(ufbj_out, filenamej)
                print 'Wrote: '+filenamej
                files_out = files_out+[filenamej]

    return files_out
# End of ufdb_update
