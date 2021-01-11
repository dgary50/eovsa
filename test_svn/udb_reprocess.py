#!/usr/bin/python

import os
import pipeline

#First check for a lock file, if none, then process
lockfilename = '/data1/processing/udb_reprocess.lock'
if os.path.isfile(lockfilename) == True:
    print lockfilename, ' Exists, returning'
else:
    print ' Starting UDB_REPROCESS'
    f = open(lockfilename, 'w')
    f.write("UDB process in progress")
    f.close
    
    # process 5 days
    pipeline.udb_reprocess(ndays=5)

    print 'Removing ', lockfilename
    os.remove(lockfilename)
    print ' Ending UDB_REPROCESS'

#That's it for now

    

