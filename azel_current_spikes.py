import dbutil as db
from util import Time
import numpy as np

def show_spikes(az,el,filename=None):
    naz = len(az)
    nel = len(el)
    n = naz
    if nel > n:
        n = nel
    
    if filename == None:
        print '                Azimuth               |                  Elevation'
        print 'Ant   Universal Time    Current (mA)  |  Ant   Universal Time    Current (mA)'
    
        for i in range(n):
            if i < naz:
                line = ' %2d %19s %12.3f  |' % (az[i][0], Time(az[i][1],format='lv').iso[:19], az[i][2])
            else:
                line = '                                      |'
        
            if i<nel:
                line = line + '   %2d %19s %12.3f' % (el[i][0], Time(el[i][1],format='lv').iso[:19], el[i][2])
            else:
                line = line + '                                      ' 
    
            print line
    else:
        file1 = open(filename, "w")                               
        file1.write('                Azimuth               |                  Elevation\n')
        file1.write('Ant   Universal Time    Current (mA)  |  Ant   Universal Time    Current (mA)\n')
        for i in range(n):
            if i < naz:
                line = ' %2d %19s %12.3f  |' % (az[i][0], Time(az[i][1],format='lv').iso[:19], az[i][2])
            else:
                line = '                                      |'
        
            if i<nel:
                line = line + '   %2d %19s %12.3f\n' % (el[i][0], Time(el[i][1],format='lv').iso[:19], el[i][2])
            else:
                line = line + '                                      \n' 
    
            file1.write(line)
        
        file1.close()

def findmotorcurrentspikes(trange,antlist=None,thresh=100,filename=None):
    if antlist == None:
        antlist = range(1,14)
        
    startday = int(trange[0].mjd)
    endday = int(trange[1].mjd+1)
    az=[]
    el=[]
    
    for d in range(startday,endday):
        for i,a in enumerate(antlist):
            if d == startday:
                t1 = trange[0]
            else:
                t1 = Time(d,format='mjd')
            
            if d == endday:
                t2 = trange[1]
            else:
                t2 = Time(d+1,format='mjd')
            
            fld = ['Ante_Cont_AzimuthMotorCurrent','Ante_Cont_ElevationMotorCurrent']  
            data, msg = db.loadsfdata(fld,[t1,t2],a)
            
            for f in data.keys():
                if f != 'Timestamp':
                    din=np.abs(np.array(data[f]))
                    elements, = np.where((din>thresh) & (din < 10000))
                    if elements.size > 0:
                        for j in elements:
                            if "Azimuth" in f:
                                az.append([a, data['Timestamp'][j],data[f][j]])
                            else:
                                el.append([a, data['Timestamp'][j],data[f][j]])
    
    show_spikes(az,el)
    if filename != None:
        show_spikes(az,el,filename)
            
    return az,el
    
    
    
