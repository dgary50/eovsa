def html_movie(t,flare=False):
    ''' After the quick-look images are created for a given date, this routine will
        write the movie.html file that allows them to be viewed as a movie.  Just
        call this with a Time() object containing the desired date.
    '''
    import glob
    from util import Time
    datstr = t.iso[:10]
    tstr = t.iso[11:16]
    if flare:
        dir = '/common/webplots/qlook_flares/'+datstr.replace('-','/')+'/'+tstr+'/png/'
    else:
        dir = '/common/webplots/qlookimg_10m/'+datstr.replace('-','/')+'/'
    files = glob.glob(dir+'/*.png')
    files.sort()
    nfiles = len(files)
    f = open('/common/webplots/qlookimg_10m/2017/08/21/eclipse_20170821.html','r')
    lines = f.readlines()
    nlines = len(lines)
    f.close()
    skiplines = []
    for i,line in enumerate(lines):
        k = line.find('var imax')
        if k != -1:
            lines[i] = line[:10]+'{:3d}'.format(nfiles)+line[13:]
        k = line.find('urls[')
        if k != -1: 
            skiplines.append(i)
    #print skiplines
    f = open(dir+'/movie_'+datstr.replace('-','')+'.html','w')
    for i in range(skiplines[1]-1):
        f.write(lines[i])
    for i in range(nfiles):
        if flare:
            f.write('urls[{:d}]=url_path+"'.format(i)+files[i][50:]+'";\n')
        else:
            f.write('urls[{:d}]=url_path+"'.format(i)+files[i][40:]+'";\n')
    for i in range(skiplines[-1]+1,nlines):
        f.write(lines[i])
    f.close()
    