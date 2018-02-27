def html_movie(t, dir=None, imgprefix='', htmlname=None):
    ''' After the quick-look images are created for a given date, this routine will
        write the movie.html file that allows them to be viewed as a movie.  Just
        call this with a Time() object containing the desired date.
    '''
    import glob
    from util import Time
    datstr = t.iso[:10]
    if not dir:
        dir = '/common/webplots/qlookimg_10m/' + datstr.replace('-', '/') + '/'
    files = glob.glob(dir + '/{}*.png'.format(imgprefix))
    files.sort()
    nfiles = len(files)
    f = open('/common/webplots/qlookimg_10m/2017/08/21/eclipse_20170821.html', 'r')
    lines = f.readlines()
    nlines = len(lines)
    f.close()
    skiplines = []
    for i, line in enumerate(lines):
        k = line.find('var imax')
        if k != -1:
            lines[i] = line[:10] + '{:3d}'.format(nfiles) + line[13:]
        k = line.find('urls[')
        if k != -1:
            skiplines.append(i)
    #print skiplines
    if not htmlname:
        htmlname = dir + '/movie_' + datstr.replace('-', '') + '.html'
    f = open(htmlname, 'w')
    for i in range(skiplines[1] - 1):
        f.write(lines[i])
    for i in range(nfiles):
        f.write('urls[{:d}]=url_path+"'.format(i) + files[i][40:] + '";\n')
    for i in range(skiplines[-1] + 1, nlines):
        f.write(lines[i])
    f.close()
    print('html saved to {}.'.format(htmlname))
