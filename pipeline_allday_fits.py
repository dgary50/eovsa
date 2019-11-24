if __name__ == '__main__':
    import pipeline_cal as pc
    import eovsa_fits as ef
    import glob
    from util import Time
    import sys, os
    from __future__ import print_function

    print(sys.argv)
    try:
        argv = sys.argv[3:]
        year = argv[0]
        month = argv[1]
        day = argv[2]
        t = Time([year+'-'+month+'-'+day+' 20:00:00'])
        if '--clearcache' in argv:
            clearcache = True
        else:
            clearcache = False
    except:
        print('Error interpreting command line arguments--will analyze data from yesterday.')
        # No arguments (or no arguments given), so default to yesterday's data to analyze
        mjdnow = Time.now().mjd
        t = Time(mjdnow-1,format='mjd')
        year = t.iso.split('-')[0]
        month = t.iso.split('-')[1]
        day = t.iso.split('-')[2].split(' ')[0]
        t = Time([t.iso])
        clearcache = True
    # Change to standard working directory and delete any existing IDB files there
    outpath = '/data1/dgary/HSO/'
    os.chdir(outpath)
    os.system('rm -rf IDB*')
    # Run first (and lengthy!) task to create corrected IDB files for the entire day
    pc.allday_udb_corr(t, outpath=outpath)
    # Process the entire day's IDB files to create fits files
    pc.allday_process(path=outpath)
    files = glob.glob(outpath+year+'/'+month+'/'+day+'/*_TP_*.fts')
    files = files.sort()
    spec = ef.eovsa_combinefits(files, freqgaps=True, outpath=outpath, ac_corr=True, doplot=False)
    files = glob.glob(outpath+year+'/'+month+'/'+day+'/*_X_*.fts')
    files = files.sort()
    spec = ef.eovsa_combinefits(files, freqgaps=True, outpath=outpath, ac_corr=True, doplot=False)
    if clearcache:
        os.system('rm -rf IDB*')
