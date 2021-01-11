def udb_process2(ndays=7):
    ''' This is the main program, that reads the FDB files, decides
    what to do, if anything, processes a scan, updates the DB and
    finishes, IDB and UDB fits processing is now done here, jmm,
    2016-01-04'''

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
                antlist = '' #just in case it's not defined, but it always should be
                for j in range(n1):
                    if ifb_1[j].pstatus == 0:
                        print 'UDB_PROCESS: Processing IDB: '+ifb_1[j].fileid
                        fj, antlist = process_1ifb(ifb_1[j])
                    #end if
                    ifb_1[j].pstatus = 1
                    ifb[ss_ifb_1[j]].pstatus = 1
#                    print 'UDB_PROCESS: Reset status for:'+ifb[ss_ifb_1[j]].fileid
#                    print 'UDB_PROCESS: Reset status for:'+ifb_1[j].fileid
                #End j loop

                udbfitsfile = udb_fitsfile(ufb_1, antlist)

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
# End of udb_process2

def udb_fitsfile(ufb0, antlist):
    '''Creates a UDB fits file for the input ufb entry 2015-01-06, jmm'''

    #string may be a list
    if isinstance(ufb0, list):
        ufb = ufb0[0]
    else:
        ufb = ufb0
    #end else

    #extract year, to get full file path
    ufileid = ufb.fileid
    yyyy = ufileid[3:7]
    ufilename = udbdir+yyyy+'/'+ufileid
    filelist = [ufileid]
    #Check validity
    file_test = dump_tsys_ext.valid_miriad_dataset(filelist)
    if len(file_test) == 0 or file_test[0] == False :
        print "udb_fitsfile: Bad Miriad Dataset:"
        print filelist[0]
        return ''
    #endif
    #If you are here you have a good IDB dataset, error check anyway
    xdat = dump_tsys_ext.rd_miriad_tsys_file(filelist)
    if xdat == None:
        print "udb_fitsfile: Bad Miriad Dataset (2):"
        print filelist[0]
        return ''
    #endif

    xdat['antennalist'] = antlist
    xdat['scanid'] = ufb.scanid
    xdat['proj'] = ufb.projectid

    #Calibrate, but only for normal observing
    sid = ufb.sourceid[:3]
    prid = ufb.projectid[:15]
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
    file_out = dump_tsys_ext.tsys_writeudbfits(xdat, calflag)
    if len(file_out) == 0:
        print "udb_fitsfile: write failed:"
        print filelist[0]
        return file_out, ''
    #endif

    return file_out
#End of udb_fitsfile

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
