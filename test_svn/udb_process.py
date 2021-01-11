#!/usr/bin/python

import os, sys
import pipeline
import time

#First check for a lock file, if none, then process
lockfilename = '/data1/processing/udb_process.lock'
if os.path.isfile(lockfilename) == True:
    t = time.localtime()
    print time.asctime(),lockfilename, ' Exists, returning'
    sys.stdout.flush()
else:
    print ' Starting UDB_PROCESS'
    sys.stdout.flush()
    f = open(lockfilename, 'w')
    f.write("UDB process in progress")
    f.close
    
    # process 1 scan
    pipeline.udb_process()

    print 'Removing ', lockfilename
    os.remove(lockfilename)
    print ' Ending UDB_PROCESS'
    sys.stdout.flush()

#That's it for now

    

