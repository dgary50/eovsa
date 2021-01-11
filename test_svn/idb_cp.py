#!/usr/bin/python

#python script to copy old IDB files from /dppdata1/IDB to
#/data1/eovsa/fits/IDB/YYYYMMDD

import os
import shutil
import mirtest

idbdir = '/dppdata1/IDB/'
idbfinaldir = '/data1/eovsa/fits/IDB/'

filelist = os.listdir(idbdir)
print filelist[0:10]
print len(filelist)

#filelist has all of the files
for j in range(len(filelist)):
#for j in range(10):
    filename = filelist[j]
    init_filename = idbdir+filename
    ok, ok_list, bad_list = mirtest.valid_miriad_dset0([init_filename])
    if ok[0] == True:
        yyyymmdd = filename[len(filename)-14:len(filename)-6]
        yyyy = yyyymmdd[0:4]
        if(yyyy == '2014' or yyyy == '2015'):
#            print "idb_cp: ", filename, " will be copied"
            outdir = idbfinaldir+yyyymmdd
            if os.path.isdir(outdir) != True:
                os.mkdir(outdir)
            #end if
            full_filename = outdir+'/'+filename
            #If the file exists, you need to delete it before using shutil.copytree
            if os.path.isdir(full_filename) == True:
#                print "idb_cp: dataset: ", full_filename, " will be overwritten"
                shutil.rmtree(full_filename)
            #end if
            shutil.copytree(init_filename, full_filename)
        #end if
    #end if
#    else:
#        print "idb_cp: ", filename, " will not be copied"
    #end else
#end for

#End of idb_cp
