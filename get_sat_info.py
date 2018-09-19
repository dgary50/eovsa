#
# Gets communication satellite information from the Lyngsat web pages,
# including name, longitude, frequencies, and polarizations
#
# Another important web site is http://www.satbeams.com/satellites,
# which has links to each satellite with beacon and beam information.
#
# History:
#   2016-Mar-22  DG
#      Started this history log.  Added "names" keyword to get_sat_info(),
#      to return information only for the satellites matching the names
#      list.  Names are converted to upper case for comparison, so case
#      is not significant.
#   2016-Aug-05  DG
#      Update get_sat_info() for plotting if the doplot keyword is True.
#
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
    try:
        satname = lines[2].split('<title>')[1].split(' at ')[0]
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
    # Find unique frequencies, and use indexes to get corresponding poln
    freq, idx = np.unique(np.array(freq).astype('int'),return_index=True)
    poln = np.array(poln)[idx]
    out = {'name':satname, 'loc':satloc, 'freqlist':freq, 'pollist':poln}
    return out

def get_sat_info(names=None,doplot=False):
    import matplotlib.pylab as plt
    from util import Time
    f = urllib2.urlopen('http://www.lyngsat.com/tracker/america.html')
    lines = f.readlines()
    f.close()
    found_names = []
    # Convert names list (if any) to upper case
    if not names is None:
        names = np.array([name.upper() for name in names])
    for line in lines:
        if line.find('<font face="Verdana"><font size=2><a href="https://www.lyngsat.com/tracker/') != -1:
            if line.find('bgcolor=#ffffff') == -1:
                name = line.split('https://www.lyngsat.com/tracker/')[1].split('.html')[0]
                if names is None:
                    # If no name list given, mark all satellites as found
                    found = True
                else:
                    # Satellite is found if name is in names
                    found = len(np.where(name.upper() == names)[0]) == 1
                if found:
                    found_names.append(name)
    out = []
    for name in found_names:
        print 'Reading information for satellite',name,
        outi = sat_info(name)
        if outi is not None:
            out.append(outi)
    if doplot:
        for i,sat in enumerate(out):
            nf = len(sat['freqlist'])
            plt.plot((float(sat['loc'])-118)*np.ones(nf),sat['freqlist'],'.')
            plt.text(float(sat['loc'])-118,sat['freqlist'][0],str(i),ha='center',va='top')
            plt.xlabel('HA [deg]')
            plt.ylabel('Frequency [MHz]')
            plt.title('Geosat Information for '+Time.now().iso)
    return out

def print_sat_names(satlist, Freq1, Freq2):
    for sat in satlist:
        f = sat['freqlist']
        idx, = np.where (np.logical_and(f > Freq1,f < Freq2))    
        if len(idx) == 0:
            pass
        else: 
            print sat['name'], sat['loc'], f[idx]

