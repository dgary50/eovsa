#!/usr/bin/python

import os
import pipeline

#First check for a lock file, if none, then process
lockfilename = '/data1/processing/udb_process_spec.lock'
if os.path.isfile(lockfilename) == True:
    print lockfilename, ' Exists, returning'
else:
    print ' Starting UDB_PROCESS_SPEC'
    f = open(lockfilename, 'w')
    f.write("UDB process in progress")
    f.close
    
    # process 1 scan
    pipeline.udb_process_spec(ndays=14)

    print 'Removing ', lockfilename
    os.remove(lockfilename)
    print ' Ending UDB_PROCESS_SPEC'

#That's it for now

    

