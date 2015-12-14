import urllib2
import numpy as np

def sat_info(name):
    namestr = 'http://www.lyngsat.com/'+name+'.html'
    try:
        f = urllib2.urlopen(namestr)
    except:
        print '404 Error'
        return None
    lines = f.readlines()
    f.close()
    satname = lines[2].split('<title>')[1].split(' at ')[0]
    try:
        satloc = lines[2].split(' at ')[1].split('\xb0')[0]
        freq = []
        poln = []
        
        for line in lines:
            if line.find('align="center"><font face="Verdana"><font size=2><b>') != -1:
                freq.append(line.split('<b>')[1].split('&nbsp')[0])
                poln.append(line.split('&nbsp;')[1][0])
        if poln == []:
            print 'No frequencies listed'
            return None
        print 'Success!'
    except:
        # Probably this satellite is not active, so return default None values
        print 'Satellite not active'
        return None
    out = {'name':satname, 'loc':satloc, 'freqlist':np.array(freq).astype('int'), 'pollist':np.array(poln)}
    return out

def get_sat_info():
    f = urllib2.urlopen('http://www.lyngsat.com/tracker/america.html')
    lines = f.readlines()
    f.close()
    names = []
    for line in lines:
        if line.find('<font face="Arial"><font size=2><a href="http://www.lyngsat.com/tracker/') != -1:
            names.append(line.split('http://www.lyngsat.com/tracker/')[1].split('.html')[0])
    out = []
    for name in names:
        print 'Reading information for satellite',name,
        outi = sat_info(name)
        if outi is not None:
            out.append(sat_info(name))
    return out

def print_sat_names(satlist, Freq1, Freq2):
    for sat in satlist:
        f = sat['freqlist']
        idx, = np.where (np.logical_and(f > Freq1,f < Freq2))    
        if len(idx) == 0:
            pass
        else: 
            print sat['name'], sat['loc'], f[idx]

