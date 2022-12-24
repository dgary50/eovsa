def make_fsq(band):

    f = open('band'+str(band)+'.fsq','w')
    print 'band'+str(band)+'.fsq'
    f.write('LIST:DWELL '+'1,'*51+'1\n')
    f.write('LIST:SEQUENCE '+str(band)+'\n')
    f.close()

def make_flare_fsq(dwellband, bandlist=None, acc=True, check=True, verbose=False):
    ''' Make a "fast-mode" frequency sequence that dwells on a
        a given band and monitors another set of bands (logarithmically
        spaced over the remaining frequencies unless specified by 
        bandlist).  The fseq file is printed to the screen and saved
        in /tmp/.
        
        dwellband    integer band number, from 1 to 52 (but do not use 2 or 3)
        bandlist     optional list of integer band numbers 1-52.  If you want
                       only the dwell band, give a bandlist that is only the
                       dwellband (but enclose in []).  Default bandlist is
                       [5,6,8,11,14,18,24,31,40,52]
        acc          optional boolean.  If True (default), send the fseq file 
                       to the ACC.
        check        optional boolean.  If True (default), perform a check
                       that the sequence does not exceed the maximum 512
                       channels that dppxmp can handle.
                       the terminal.
        verbose      optional boolean.  If True, print the final fseq file to
        Examples:
            make_flare_fsq(13,[13])     dwells full time on band 13
            make_flare_fsq(29)          dwells for 0.8 on band 29, plus
                                            takes snapshot on 
    '''
    snapshot_time = ['0.01975']*52   # Default duration for all bands is 20 ms
    if bandlist:
        nbands = len(bandlist)
    else:
        nbands = 10
        bandlist = [5,6,8,11,14,18,24,31,40,52]  # Standard bands
    if dwellband in bandlist:
        # Avoid repeated band
        nbands -= 1
    else:
        # Insert dwellband into bandlist
        bandlist.insert(0,dwellband)
        bandlist.sort()
    dwell_time = 1.0 - nbands*0.02
    snapshot_time[dwellband-1] = str(dwell_time)
    str1 = 'LIST:DWELL '
    for s in snapshot_time:
        str1 += s+','
    str1 = str1[:-1]
    if verbose: print str1
    str2 = 'LIST:SEQUENCE '
    for b in bandlist:
        str2 += str(b)+', '
    str2 = str2[:-2]
    if verbose: print str2
    if nbands > 0:
        fsqfile = 'dwell'+str(dwellband)+'+'+str(nbands)+'.fsq'
    else:
        fsqfile = 'dwell'+str(dwellband)+'.fsq'
    f = open('/tmp/'+fsqfile,'w')
    if verbose: print fsqfile
    f.write(str1+'\n')
    f.write(str2+'\n')
    f.close()
    if acc:
        from ftplib import FTP
        acc = FTP('acc.solar.pvt')
        acc.login('admin','observer')
        acc.cwd('parm')
        # Send file to ACC
        f = open('/tmp/'+fsqfile,'r')
        acc.storlines('STOR '+fsqfile,f)
        f.close()
    if check and not acc:
        print 'Must use acc=True if check=True!'
        return None
    if check:
        import chan_info_52 as ci
        print ''
        chinfo = ci.Chan_Info()
        chinfo.fseq2nsavg(fsqfile)
        if chinfo.fast_nsavg in chinfo.nsavg:
            print 'Fseq file',fsqfile,'is valid and results in',chinfo.tot_scichan(),'science channels.'
            return fsqfile
        else:
            print 'Fseq file',fsqfile,'is NOT valid and will be ignored. It has > 512 science channels.'
            return None
    else:
        return fsqfile
