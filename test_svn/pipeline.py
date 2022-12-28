# Note mirexec, numpy are found in /usr/lib/python2.7/dist-packages/
# Added idbfits procedure in udb_process, jmm, 23-oct-2015 changed
# spec processing to udb_process_1spec, jmm, 12-jan-2016 Dropped
# mirexec, uses aipy, creates UDB files from scratch using aipy, jmm,
# 1-mar-2016. New version uses udb_util to create new-style UDB files
# with no frequency averaging, but time averaged over 60 seconds, jmm,
# 8-jan-2017.
#
# 13-Jan-2021 dg
#   Commented out creation of FITS files and plots (long overdue!)

from numpy import *
import time
import fdb
# for deleting old datasets
import os, sys
import shutil
# for IDB fits processing, UDB file creation, dump_tsys_ext is old,
# but still needed for FITS file processing; udb_util is new
import dump_tsys_ext
import udb_util

from util import Time

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
#udbdir = '/data1/UDB/' Moved UDB files, 2016-07-08, jmm
udbdir = '/data1/eovsa/fits/UDB/'
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
global one_day
one_day = 86400.0

def udb_init_ifdbfiles_1day(day_in,fdbtype='IFDB'):
    '''Reads in the available IFDB files, for the given day, formatted YYYYMMDD'''
    # check for what type of filedb you're looking for
    # currently FDB, IFDB
    if fdbtype == 'IFDB':
        xdir = ifdbdir
    elif fdbtype == 'FDB':
        xdir = fdbdir
    elif fdbtype == 'UFDB':
        xdir = ufdbdir
    else:
        xdir = fdbdir
    #endif
    ppp = day_in
    if fdbtype == 'FDB':
        fdbfile = xdir+fdbtype+ppp+'.txt'
    else:
        yyyy = day_in[0:4]+'/'
        fdbfile = xdir+yyyy+fdbtype+ppp+'.txt'
    #endif

    #Now read in the ifdb entries from the file
    ifb = fdb.pfdb_read(fdbfile)
    return ifb

#End of udb_init_ifdbfiles_1day

def udb_init_fdbfiles(ndays=5,fdbtype='FDB',day0=None):
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
    #endif
    if day0 != None:
        #Ok time in seconds from 1970
        tform = '%Y-%m-%d %H:%M:%S'
        ptime = time.mktime(time.strptime(day0, tform))
    else:
        #local time in seconds from 1970
        ptime = time.mktime(time.localtime())
    #endelse
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
        #endif
    #endfor

    #reverse order
    fdbfiles.reverse()
    return fdbfiles
#End of udb_init_fdbfiles

def udb_fb2process(ndays=5,day0=None):
    '''Reads in the last N days of FDB files, and IFDB files, finds
    the IDB files that need processing, and returns the fdb ifdb
    classes for those files'''

    #get FDB, IFDB filenames
    fdbfiles = udb_init_fdbfiles(ndays=ndays,day0=day0)
    if len(fdbfiles) == 0:
        print 'UDB_FB2PROCESS Error: No FDB files:'
        sys.stdout.flush()
        ifboops = []
        fboops = []
        return fboops, ifboops
    #endif
    ifdbfiles = udb_init_fdbfiles(ndays=ndays,fdbtype='IFDB',day0=day0)
    if len(ifdbfiles) == 0:
        # Not necessarily an error
        print 'UDB_FB2PROCESS: No IFDB files:'
        sys.stdout.flush()
    #endif

    #read in the files
    fb=[]
    for j in range(len(fdbfiles)):
        fbj = fdb.fdb_read(fdbfiles[j])
        if len(fbj) > 0:
            fb = fb+fbj
        #endif
    #endfor
    ifb=[]
    if len(ifdbfiles) > 0:
        for j in range(len(ifdbfiles)):
            ifbj = fdb.pfdb_read(ifdbfiles[j])
            if len(ifbj) > 0:
                ifb = ifb+ifbj
            #endif
        #endfor
    #endif

    #check for non-empty fb's first
    if len(fb) == 0:
        print "no data to process"
        sys.stdout.flush()
        otp = []
        return fb, otp
    #endif

    #if I am here, I have IDB files from the FDB. Do I have PFDB entries?
    if len(ifb) == 0:
        print "no IFDB entries, process all FDB entries"
        sys.stdout.flush()
        #For each FDB entry, create a IFDB entry, with pstatus = 0
        for j in range(len(fb)):
            ifbj = fdb.pfiledb(fb[j].fileid,fb[j].scanid,fb[j].sourceid,fb[j].projectid,fb[j].st_ts,fb[j].en_ts,0)
            ifb.append(ifbj)
        #endfor
        return fb, ifb
    #endif

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
        #endif
        ifb_out.append(ifb_outj)
    #endfor
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
        sys.stdout.flush()
        return ifb_out, ss_ifb_out
    #endif
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
            #keep non-zero pstatuses, jmm, 2013-10-23
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
        
def udb_process1scan(ifb_1,reprocess=False):
    '''Given the set of ifb entries, process them into a single UDB
    Miriad dataset for this scan'''

    n1 = len(ifb_1)
    filelist = fdb.fdb_list_fileid(ifb_1)
    print filelist
    sys.stdout.flush()
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

    # Call udbfile_create for this filelist
    filelist_full = []
    for j in range(len(filelist)):
        #do not process files with pstatus = 666
        if ifb_1[j].pstatus != 666:
            if reprocess == True:
                filelistj = filelist[j]
                yyyymmdd = filelistj[len(filelistj)-14:len(filelistj)-6]
                filelist_full.append(idbfinaldir+yyyymmdd+'/'+filelist[j])
            else:
                filelist_full.append(idbdir+filelist[j])
            #endelse
        #endif
    #endfor
    n1 = len(filelist_full)
    if n1 == 0:
        print 'UDB_PROCESS1SCAN: No good files to process'
        sys.stdout.flush()
        return [], []
    #endif
    if n1 > 50:
        print "UDB_PROCESS1SCAN: Too many Files N1 = ", n1
        print "UDB_PROCESS1SCAN: Reset to 50"
        sys.stdout.flush()
        n1 = 50
        filelist_full = filelist_full[0:n1-1]
    #endif

    #If the file exists, you need to delete it
    if os.path.isdir(ufilename) == True:
        print "UDB_PROCESS1SCAN: dataset: ", ufilename, " will be deleted"
        sys.stdout.flush()
        shutil.rmtree(ufilename)
    #endif
    ufile_out, bad_filename = udb_util.udbfile_create(filelist_full, ufilename)
    if(len(bad_filename) > 0):
        print type(ufile_out)
        print "UDB_PROCESS1SCAN: Error creating: ", ufilename
        print "UDB_PROCESS1SCAN: Bad Filename: ", bad_filename
        sys.stdout.flush()
        #strip the path out of the bad_filename here
        bbb = bad_filename[len(bad_filename)-17:]
        return ufile_out, bbb
    #endif
    if(len(ufile_out) == 0):
        print "UDB_PROCESS1SCAN: Error creating: ", ufilename
        sys.stdout.flush()
        return [], []
    else:
        return [ufb], []
    #endif

# End of udb_process1scan

def udb_process(ndays=5,day0=None,reprocess=False):
    ''' This is the main program, that reads the FDB files, decides
    what to do, if anything, prcesses a scan, updates the DB and
    finishes'''

    #initialize output
    ufb_1 = []
    #first read the FDB files from DPP, and determine if there are any
    #files that need to be processed. Here fb are from the DPP, and
    #will not be touched. The ifb are the ones that we will work with
    #and should have the same number of elements as the fb; both are
    #lists
    fb, ifb = udb_fb2process(ndays=ndays,day0=day0)
    #Check for success
    if len(ifb) == 0:
        print 'UDB_PROCESS Error: No filedbs input, No processing'
        sys.stdout.flush()
    else:
        #Now process the first scan with any files with pstatus = 0
        ifb_1, ss_ifb_1 = udb_scan2process(ifb)
        #Only process if you have to
        if len(ifb_1) == 0:
            print 'UDB_PROCESS: No Scans to process, finished successfully'
            sys.stdout.flush()
        else:
            print 'UDB_PROCESS: Processing Scan: '+ifb_1[0].scanid
            sys.stdout.flush()
            ufb_1, bad_filename = udb_process1scan(ifb_1,reprocess=reprocess)
            if len(bad_filename) > 0:
                ifb_1_list = fdb.fdb_list_fileid(ifb_1)
                n1 = len(ifb_1) #change status for bad file to 666
                print ifb_1_list
                print bad_filename
                for j in range(n1):
                    if ifb_1_list[j] == bad_filename:
                        ifb_1[j].pstatus = 666
                        ifb[ss_ifb_1[j]].pstatus = 666            
                    #endif
                #endfor
                #Update IFDB
                ifb_files_out = update_ifdb(ifb)
                print "Updated IFDB for Bad filename: "
                sys.stdout.flush()
                return []   # ufb_1 #This will probably segmentation fault, but maybe not -- yes, it DOES. Now return empty list
            #endif bad filename
            #Full database management
            if len(ufb_1) == 0:
                print 'UDB_PROCESS Error: No output from UDB_PROCESS1SCAN'
            else:
                #Update ifdb.pstatus for the processed files, output
                #idbfits and copy IDB files if necessary
                n1 = len(ifb_1)
#                have_antlist = False
                for j in range(n1):
                    if ifb_1[j].pstatus != 666: #check as bad files will crash
#                        if reprocess == False:
                        if ifb_1[j].pstatus == 0:
                            print 'UDB_PROCESS: Processing IDB: '+ifb_1[j].fileid
                            process_1ifb(ifb_1[j], reprocess=reprocess)
#                            fj, antlist = process_1ifb(ifb_1[j], reprocess=reprocess)  # process_1ifb() no longer returns anything
#                            have_antlist = True
#                        else: #only do this once for antlist
#                            if ifb_1[j].pstatus == 0: # and have_antlist == False:
#                                print 'UDB_PROCESS: Processing IDB: '+ifb_1[j].fileid
#                                fj, antlist = process_1ifb(ifb_1[j], reprocess=reprocess)
#                                have_antlist = True
                        ifb_1[j].pstatus = 1
                        ifb[ss_ifb_1[j]].pstatus = 1
#                    else:
#                        antlist = '1 2 3 4 5 6 7 8 9 10 11 12 13 14 15, 16'
                    #endelse
#                    print 'UDB_PROCESS: Reset status for:'+ifb[ss_ifb_1[j]].fileid
#                    print 'UDB_PROCESS: Reset status for:'+ifb_1[j].fileid
                #End j loop
                #Write out the new ifb entries
                ifb_files_out = update_ifdb(ifb)
                
                #here process the UDB spec fits files by spawning an IDL process
#                ufb_1a = udb_process_1spec(ufb_1, antlist, reprocess=reprocess)

#                if len(ufb_1a) == 0:
#                    print 'UDB_PROCESS: Problem with Scan to FITS: '+ifb_1[0].scanid
#                    sys.stdout.flush()
                #Write out ufdb entry, ufb_1 is a 1 element list
#                    ufb_files_out = update_ufdb(ufb_1, ndays=ndays, day0=day0)
#                else:
                ufb_files_out = update_ufdb(ufb_1, ndays=ndays, day0=day0)
#                ufb_files_out = update_ufdb(ufb_1a, ndays=ndays, day0=day0)
                print 'UDB_PROCESS: Processed Scan: '+ifb_1[0].scanid
                print 'UDB_PROCESS: finished successfully'
                sys.stdout.flush()
                #Endelse
            #Endelse
        #Endelse
    #Endelse
    return ufb_1
# End of udb_process

def process_1ifb(ifb0, reprocess=False):
    # Creates an IDB fits file for the input ifb entry (no longer does this)
    '''
    Copies the IDB Miriad file to the local disk. Added antennalist output
    2015-01-06, jmm'''

    #ifb0 may be a list
    if isinstance(ifb0, list):
        ifb = ifb0[0]
    else:
        ifb = ifb0
    #end else
    filelist = [idbdir+ifb.fileid]
    #Check validity
    file_test = udb_util.valid_miriad_dataset(filelist)
    if len(file_test) == 0 or file_test[0] == False :
        print "process_1ifb: Bad Miriad Dataset:"
        print filelist[0]
        sys.stdout.flush()
        return #'', ''
    #endif
    #If you are here you have a good IDB dataset, error check anyway
#    xdat = dump_tsys_ext.rd_miriad_tsys_file(filelist)
#    if xdat == None:
#        print "process_1ifb: Bad Miriad Dataset (2):"
#        print filelist[0]
#        sys.stdout.flush()
#        return '', ''
    #endif
    #Hold onto antennalist for output
#    antlist = xdat['antennalist']
    #Calibrate, but only for normal observing, No calibration
#    sid = ifb.sourceid[:3]
#    prid = ifb.projectid[:15]
#    calflag = False
#    if sid == "Sun" and prid == "NormalObserving":
#        t1 = Time(xdat['ut_mjd'][0],format='mjd')
#        fghz, calfac, offsun = calibration.sp_read_calfac(t1)
#        if fghz != None:
#            xdat = calibration.sp_apply_cal(xdat, fghz, calfac, offsun)
#            calflag = True
#        else:
#            calflag = False
#        #endif
#    else:
    #endif
    #write the file
#    print "CALFLAG", calflag
#    file_out = dump_tsys_ext.tsys_writetofits(xdat, calflag)
#    if len(file_out) == 0:
#        print "process_1ifb: write failed:"
#        print filelist[0]
#        sys.stdout.flush()
#        return file_out, ''
    #endif
    if reprocess == True:
        #you are done
        return #file_out, antlist
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
        sys.stdout.flush()
        shutil.rmtree(full_filename)
    #end if
    shutil.copytree(filename, full_filename)
#the second output is passed into udb_fitsfile
    return #file_out, antlist
#End of process_1ifb

def udb_process_1spec(ufb_1_in, antlist, reprocess=False):
    '''This process inputs a single UFDB entry, creates a FITS file,
    then spawns an IDL process to make summary plots'''

   #string may be a list                                                                                                           
    if isinstance(ufb_1_in, list):
        ufb = ufb_1_in[0]
    else:
        ufb = ufb_1_in
    #end else                                                                                                                      

    #If I am here, I have a valid ufb, but do I have status = 0, If
    #so, process it.
    if ufb.pstatus == 0:
        do_this_file = 1
    else:
        do_this_file = -1
    #end else                                                                                                                      
    #if do_this_file is -1, then no processing happens
    #check for reprocess
    if reprocess == True:
        idlfile = "/home/user/workdir/udb_reprocess_spec.pro"
    else:
        idlfile = "/home/user/workdir/udb_process_spec.pro"
    #endelse
    if do_this_file == -1:
        print "UDB_PROCESS_1SPEC: no data to process"
        print 'UDB_PROCESS_1SPEC: finished successfully'
        sys.stdout.flush()
        ufb_out = []
        f = open(idlfile, 'w')
        s = "pr_path_lib, 'eovsa_write_pwrfits', /multi\n"
        f.write(s)
        s = "message, /info, 'No data to process'\n"
        f.write(s)
        s = "end\n"
        f.write(s)
        f.close
        return ufb_out
    #endif

    ufileid = ufb.fileid
    uprojectid = ufb.projectid
    uantlist = antlist
    print 'UDB_PROCESS_1SPEC - Antlist: ', uantlist
    yyyy = ufileid[3:7]
    print 'UDB_PROCESS_1SPEC - Processing: ', ufileid
    sys.stdout.flush()
    # Check for YYYY directories
    if os.path.isdir(udbdir+yyyy) != True:
        os.mkdir(udbdir+yyyy)
    if os.path.isdir(udbtxtdir+yyyy) != True:
        os.mkdir(udbtxtdir+yyyy)
    ufilename = udbdir+yyyy+'/'+ufileid
    # write a fits file, ignoring tpcal for the moment
    print ufilename
    sys.stdout.flush()
    udat = dump_tsys_ext.rd_miriad_tsys_file(ufilename)
    if udat == None:
        print "UDB_PROCESS_1SPEC: Bad Miriad Dataset (2):"
        print ufilename
        sys.stdout.flush()
        ufb_out = []
        return ufb_out
    #endif
    calflag = False
    ufile_out = dump_tsys_ext.tsys_writeudbfits(udat, calflag)
    if len(ufile_out) == 0:
        print "UDB_PROCESS_1SPEC: write failed:"
        print filelist[0]
        sys.stdout.flush()
        return []
    #endif

    # Here figure out how to call the IDL process from here, write a
    # batch file and run ssw_batch using the subprocess command
    # 2014-06-26. Instead, just create the batch file and we will run
    # ssw_batch directly from the shell script

    # 2014-12-15 added projectid to eovsa_write_pwrfits.pro call
    # 2016-04-19 The fits file already has been written, so just plot it
 
    if reprocess == True:
        idlfile = "/home/user/workdir/udb_reprocess_spec.pro"
    else:
        idlfile = "/home/user/workdir/udb_process_spec.pro"
    #endelse
    f = open(idlfile, 'w')
    s = "pr_path_lib, 'eovsa_write_pwrfits', /multi\n"
    f.write(s)
    s = "fitsfile = '"+ufile_out+"'\n"
    f.write(s)
    s = "eovsa_plot_tpower_fits, fitsfile\n"
    f.write(s)
    s = "end\n"
    f.write(s)
    f.close()

    ufb.pstatus = 1
    ufb_out = [ufb]

    print 'UDB_PROCESS_1SPEC: Processed Spectra for: '+ufileid
    print 'UDB_PROCESS_1SPEC: finished successfully'
    sys.stdout.flush()

    return ufb_out
# End of udb_process_1spec

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
    #endfor
    #Get the filenames for the files,
    days = fdb.fdb_uniq_day(ifb)
    ndays_out = len(days)
    if ndays_out == 0:
        print 'UPDATE_IFDB Error: No days to process'
        sys.stdout.flush()
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
                sys.stdout.flush()
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
                sys.stdout.flush()
                files_out = files_out+[filenamej]
            #endif
        #endfor
    #endif
    return files_out
    
# End of update_ifdb

def update_ufdb(ufb, ndays=5, day0=None):
    '''Updates ufdb files, need to read the last ndays files, check to
    see if you're overwriting an entry, if not then append the new
    scan to the list and output. If so, replace the old entry, then
    output. Note that ufb is a 1-element list'''

    if len(ufb) == 0:
        print 'UPDATE_UFDB ERROR: No UFB input'
        sys.stdout.flush()

    files_out = []
    #get filenames
    ufdb_files = udb_init_fdbfiles(ndays=ndays,fdbtype='UFDB',day0=day0)
    
    ufb_in = []
    if len(ufdb_files) == 0:
        #maybe there are none
        print 'UPDATE_UFDB: No UFDB files in range:'
        sys.stdout.flush()
    else:
        for j in range(len(ufdb_files)):
            ufbj = fdb.pfdb_read(ufdb_files[j])
            if len(ufbj) > 0:
                ufb_in = ufb_in+ufbj
            #endif
        #endfor
    #endif
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
                #endfor
                ufb_out = ufb_out+ufb
            #endif
        #endif
    #endif
    #Here, now output the ufb_out entries
    files_out = []
    #Get days for each ufb entry
    ufdays = fdb.fdb_list_fileid(ufb_out)
    for j in range(len(ufdays)):
        y = str.find(ufdays[j], 'DB')
        x0 = y+2
        x1 = y+10
        ufdays[j] = ufdays[j][x0:x1]
    #endfor

    #Get the filenames for the files,
    days_out = fdb.fdb_uniq_day(ufb_out)
    #Adjustment due to lack 
    ndays_out = len(days_out)
    if ndays_out == 0:
        print 'UPDATE_UFDB Error: No days to process'
        sys.stdout.flush()
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
                sys.stdout.flush()
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
                sys.stdout.flush()
                files_out = files_out+[filenamej]
            #endif
        #endfor
    #endif

    return files_out
# End of ufdb_update

def udb_reprocess_days(start_day, end_day):
    '''Reprocess UDB data from the start to end day. Reads in the IFDB files 
    for the time interval, and the previous and next days, sets status flags 
    to zero for all of the IFDB entries for the given day, and reprocesses 
    until the files are finished. Start and end days are passed in as time 
    strings YYYY-MM-DD 00:00:00 (or %Y-%m-%d HH:MM:SS). Days are inclusive, 
    so that if start_day = end_day, then that 1 day is processed.'''

    tform = '%Y-%m-%d %H:%M:%S'
    time0 = time.mktime(time.strptime(start_day, tform))
    time1 = time.mktime(time.strptime(end_day, tform))
    #number of days to load
    ndays0 = int((time1-time0)/one_day)+1
    
    from util import get_idbdir, Time
    global idbfinaldir
    tobj = Time(start_day)
    idbfinaldir = get_idbdir(tobj)

    #get ifdb files for this time range, set status to 0 for all and rewrite
    ifdb_files = udb_init_fdbfiles(ndays=ndays0, fdbtype='IFDB', day0=end_day)
    if len(ifdb_files) == 0:
        print "UDB_REPROCESS_DAYS: No files to process"
        print "Start_day = ", start_day
        print "End_day = ", end_day
        sys.stdout.flush()
        return []
    #endif
    ifb=[]
    for j in range(len(ifdb_files)):
        ifbj = fdb.pfdb_read(ifdb_files[j])
        if len(ifbj) > 0:
            ifb = ifb+ifbj
        #endif
    #endfor
    if len(ifb) == 0:
        print "UDB_REPROCESS_DAYS: No files to process ???"
        print "Start_day = ", start_day
        print "End_day = ", end_day
        sys.stdout.flush()
        return []
    #endif
    #set status to zero
    for j in range(len(ifb)):
        ifb[j].pstatus = 0
        print "Reset Pstatus: ", ifb[j].fileid
        sys.stdout.flush()
    #endfor
    #And rewrite the IFDB files
    ifdbfilesout = update_ifdb(ifb)
    #Now reprocess the ifdbs, adding days insures that we get scans
    #that overlap day boundaries
    ndays2 = 2+ndays0
    #day02 will be day0 input
    day02 = time.strftime(tform, time.localtime(time1+one_day))

    ufb_out = ['dummy']
    while len(ufb_out) > 0:
        ufb_out = udb_process(ndays=ndays2,day0=day02,reprocess=True)
        if len(ufb_out) > 0:
            print "UDB_REPROCESS_DAYS: Processed ", ufb_out
        else:
            print "UDB_REPROCESS_DAYS: No Processing "
        sys.stdout.flush()
        #endelse
    #endwhile

    print "UDB_REPROCESS_DAYS: Finished Processing "
    print "Start_day = ", start_day
    print "End_day = ", end_day
    sys.stdout.flush()
    return ['OK']
#end of udb_reprocess_days

def udb_reprocess(ndays=5):
    '''Reads in a file /data1/processing/last_day_reprocessed.txt 
    and reprocess the previous ndays days'''

    #Check for the input file:
    filename = '/data1/processing/last_day_reprocessed.txt'
    # check for ok file
    if os.path.isfile(filename) == False:
        print "no file: ", filename
        sys.stdout.flush()
        fb = []
        return fb
    else:
        print "file found", filename
        sys.stdout.flush()

    f = open(filename)
    lines = f.readlines()
    f.close()
    temp_l = lines[0]
    len_l = len(temp_l)
    time1 = temp_l[:len_l-1]
    tform = '%Y-%m-%d %H:%M:%S'
    time1 = time.mktime(time.strptime(time1, tform))
    end_day = time.strftime(tform, time.localtime(time1-one_day))
    start_day =  time.strftime(tform, time.localtime(time1-ndays*one_day))

    print "UDB_REPROCESS: PROCESSING "+start_day+" TO  "+end_day
    sys.stdout.flush()

    isok = udb_reprocess_days(start_day, end_day)

    #If processing didn't fail, write out the new file
    if len(isok) != 0:
        f = open(filename, 'w')
        s = start_day+'\n'
        f.write(s)
        f.close()
    #endif
#end of udb_reprocess
