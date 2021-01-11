# Pretty much all of these filedb programs have no error checking,
# except for the read programs which return a [-1] if files are not
# found. Error checking is to be done in the pipeline.py module

#for array manipulation
from numpy import *
#for file existence
import os

class ifiledb:

    #file_id, scan_id, source_id, start and end accumulations, start
    #and end second and Musecond, start and end stateframe timestamps
    
    def __init__(self, fileid, scanid, sourceid, projectid, st_acc, en_acc, st_sec, st_usec, en_sec, en_usec, st_ts, en_ts):
        self.fileid=fileid
        self.scanid=scanid
        self.sourceid=sourceid
        self.projectid=projectid
        self.st_acc=int(st_acc)
        self.en_acc=int(en_acc)
        self.st_sec=int(st_sec)
        self.st_usec=int(st_usec)
        self.en_sec=int(en_sec)
        self.en_usec=int(en_usec)
        self.st_ts=float(st_ts)
        self.en_ts=float(en_ts)
# Ends this class definition        

class pfiledb:
    # info for each IDB file or PDB file (pdb being generic for
    # pipeline database), and a status flag

    def __init__(self, fileid, scanid, sourceid, projectid, st_ts, en_ts, pstatus):
        self.fileid=fileid
        self.scanid=scanid
        self.sourceid=sourceid
        self.projectid=projectid
        self.st_ts=float(st_ts)
        self.en_ts=float(en_ts)
        self.pstatus=int(pstatus)
# Ends this class definition        

def fdb_read(filename):
    '''Reads in an FDB file and returns a list of ifiledb
    classes. Returns an empty list if there is no file'''

    # check for ok file
    if os.path.isfile(filename) == False:
        print "no file: ", filename
        fb = []
        return fb
    else:
        print "file found", filename

    f = open(filename)
    lines = f.readlines()
    f.close()

    #strip the header line
    lines = lines[1:]
    #remove the \n, carefully
    nlines = len(lines)
    for j in range(nlines):
        temp_l = lines[j]
        len_l = len(temp_l)
        temp_l = temp_l[:len_l-1]
        lines[j] = temp_l

    #build up a list of ifiledbs, each record has two lines for
    #readability, but if st_sec is zero or en_sec is zero, something
    #is up, and we will not include this
    fb = []
    nfb = nlines/2
    for j in range(nfb):
        l1 = 2*j
        l2 = l1+1
        line1 = str.split(lines[l1])
        line2 = str.split(lines[l2])
#        if((len(line1) != 3 and len(line1) != 4) or len(line2) != 8): temporary, 2014-12-08
        if(len(line2) != 8):
            print 'FDB_READ: Bad file entry: ', line1[0]
        else:
            if(len(line1) == 3):
                fbj = ifiledb(line1[0], line1[1], line1[2], 'NONE', line2[0], line2[1], line2[2], line2[3], line2[4], line2[5], line2[6], line2[7])
            else:
                fbj = ifiledb(line1[0], line1[1], line1[2], line1[3], line2[0], line2[1], line2[2], line2[3], line2[4], line2[5], line2[6], line2[7])
                
            if(fbj.st_ts == 0 or fbj.en_ts == 0):
                print 'FDB_READ: Bad file entry: ', fbj.fileid
            elif(fbj.st_sec == 0 or fbj.en_sec == 0):
                print 'FDB_READ: Bad file entry: ', fbj.fileid
            elif(fbj.en_ts < fbj.st_ts):
                print 'FDB_READ: Bad file entry: ', fbj.fileid
            elif(fbj.en_sec < fbj.st_sec):
                print 'FDB_READ: Bad file entry: ', fbj.fileid
            else:
                fb.append(fbj)

    return fb
#End of fdb_read

def pfdb_read(filename):
    '''Reads in a pipeline FDB file and returns a list of pfiledb
    classes. Returns an empty list if the file doesn't exist.'''

    # check for ok file
    if os.path.isfile(filename) == False:
        print "no file: ", filename
        fb = []
        return fb

    f = open(filename)
    lines = f.readlines()
    f.close()

    #strip the header line
    lines = lines[1:]
    #remove the \n, carefully
    nlines = len(lines)
    for j in range(nlines):
        temp_l = lines[j]
        len_l = len(temp_l)
        temp_l = temp_l[:len_l-1]
        lines[j] = temp_l

    #build up a list of filedbs
    fb = []
    for j in range(nlines):
        line1 = str.split(lines[j])
        #careful to set pstatus to be an integer
        if(len(line1) == 6):
            fbj = pfiledb(line1[0], line1[1], line1[2], 'NONE', line1[3], line1[4], line1[5])
        elif(len(line1) == 7):
            fbj = pfiledb(line1[0], line1[1], line1[2], line1[3], line1[4], line1[5], line1[6])
        else:
            print 'PFDB_READ: Bad file entry: ', fbj.fileid

        if(fbj.st_ts == 0 or fbj.en_ts == 0):
            print 'PFDB_READ: Bad file entry: ', fbj.fileid
        else:
            fb.append(fbj)

    return fb
#End of pfdb_read

def fdb_write(fb, filename):
    '''Write an ifiledb list to a file'''
    #sort by st_ts first
    fb = sorted(fb,key = lambda filedb: filedb.st_ts)
    f = open(filename, 'w')
    s = "FILE: SCANID: SOURCEID: PROJECTID: ST_ACC: EN_ACC: ST_SEC: ST_USEC: EN_SEC: EN_USEC: ST_TS: EN_TS:\n"
    f.write(s)
    nlines = len(fb)
    for j in range(nlines):
        s1 = fb[j].fileid+' '+fb[j].scanid+' '+fb[j].sourceid+' '+fb[j].projectid+'\n'
        f.write(s1)
        s2 = str(fb[j].st_acc)+' '+str(fb[j].en_acc)+' '+str(fb[j].st_sec)+' '+str(fb[j].st_usec)+' '+str(fb[j].en_sec)+' '+str(fb[j].en_usec)+' '+str(fb[j].st_ts)+' '+str(fb[j].en_ts)+'\n'
        f.write(s2)

    f.close
#end of fdb_write

def pfdb_write(fb, filename):
    '''Write a pfiledb list to a file'''
    #sort by st_ts first
    fb = sorted(fb,key = lambda pfiledb: pfiledb.st_ts)
    f = open(filename, 'w')
    s = "FILE: SCANID: SOURCEID: PROJECTID ST_TS: EN_TS: PSTATUS:\n"
    f.write(s)
    nlines = len(fb)
    for j in range(nlines):
        s1 = fb[j].fileid+' '+fb[j].scanid+' '+fb[j].sourceid+' '+fb[j].projectid+' '+str(fb[j].st_ts)+' '+str(fb[j].en_ts)+' '+str(fb[j].pstatus)+'\n'
        f.write(s1)

    f.close
#end of pfdb_write

def fdb_list_fileid(fb, path=''):
    '''Extract a list of fileid values from a list of filedbs'''
    nlines = len(fb)
    olist = []
    for j in range(nlines):
        if (len(path) == 0):
            olist = olist+[fb[j].fileid]
        else:
            olist = olist+[path+'/'+fb[j].fileid]

    return olist
#end of fdb_list_fileid

def fdb_list_scanid(fb):
    '''Extract a list of scanid values from a list of filedbs'''
    nlines = len(fb)
    olist = []
    for j in range(nlines):
        olist = olist+[fb[j].scanid]

    return olist
#end of fdb_list_scanid

def fdb_uniq_scan(fb):
    '''From an input list of fdb classes, extract the unique values of scanid'''
    flist = fdb_list_scanid(fb)
    scarray = array(flist)
    uel = unique(scarray)
    #uel is an array, and I need a list
    nuel = len(uel)
    u = []
    for j in range(nuel):
        u = u+[uel[j]]

    return u
#end of fdb_uniq_scan

def fdb_uniq_day(fb):
    '''From an input list of fdb classes, extract the unique values of the date'''
    flist = fdb_list_fileid(fb)

    # you need to find the unique days in the list of fileids:
    # 'FDByyyymmddhhmmss' means you need to extract subscript range 3
    # to 11 to get the date, but it may be 4 to 12
    for j in range(len(flist)):
        y = str.find(flist[j], 'DB')
        x0 = y+2
        x1 = y+10
        flist[j]=flist[j][x0:x1]

    dyarray = array(flist)
    uel = unique(dyarray)
    #uel is an array, and I need a list
    nuel = len(uel)
    u = []
    for j in range(nuel):
        u = u+[uel[j]]

    return u
#end of fdb_uniq_day

def fdb_extract_1scan(fb,scanid_in):
    '''From an input list of fdb classes, extract all of the elements
    with scanid equal to the input value, Also return the list of fileids'''
    
    flist = fdb_list_scanid(fb)
    #make an array so that numpy extract can be used
    scarray = array(flist)
    #extract the elements where the scanid is the same as the input one
    fb_temp = extract(scarray == scanid_in, fb)
    #extract the subscripts of the list, which is what you need for output
    ss_fb_out = extract(scarray == scanid_in, array(range(len(scarray))))
    #fb_out needs to be a list, and not an array
    fb_out = fb_temp.tolist()
    #make a fileid list too
    #flist_out = fdb_list_fileid(fb_out)
    return fb_out, ss_fb_out
#end of fdb_extract_1scan

def fdb_extract_trange(fb, tstart, tend):
    '''From an input list of fdb classes, extract all of the elements
    with start times Ge to trange0 and Lt to trange1, Also return the
    list of fileids'''
    

    #First get the start times
    nlines = len(fb)
    times = []
    for j in range(nlines):
        times = times+[fb[j].st_ts]

    #make an array so that numpy extract can be used
    tim_arr = array(times)

    #but I can't use And for arrays
    x = tim_arr >= tstart
    y = tim_arr < tend
    for j in range(len(x)):
        x[j] = x[j] and y[j]

    #extract the elements where the scanid is the same as the input one
    fb_temp = extract(x, fb)
    #extract the subscripts of the list, which is what you need for output
    ss_fb_out = extract(x, array(range(len(tim_arr))))
    #fb_out needs to be a list, and not an array
    fb_out = fb_temp.tolist()
    #make a fileid list too
    #flist_out = fdb_list_fileid(fb_out)
    return fb_out, ss_fb_out
#end of fdb_extract_1scan
    
     
