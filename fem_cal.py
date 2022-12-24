#
# Routine to read FEM calibration measurement files and fit
# log(voltage) vs. power with a quartic, in order to generate
# coefficients to be entered into the crio.ini files.
#
#   2016-Mar-19  DG
#     First written.
#   2021-Jul-28 OG
#     Modified so that ant can be either antenna number or 
#     list of 2 files [HPOL,VPOL] 

import numpy as np
import matplotlib.pylab as plt
import glob

def fem_cal(ant):
    ''' Reads HPOL and VPOL fem calibration files from plots 
        the data and a quartic fit, and prints the coefficients
        of the fit.
        If ant is a number, it plots the most recent data.
        If ant is a list of two files, it will plot the data
        in these files.
    '''
    
    tp=-1
    if  isinstance(ant,int):
        hfiles = glob.glob('Antenna '+str(ant)+' H*.txt')
        if len(hfiles) > 1: hfiles = np.sort(hfiles)
        print(hfiles[-1])
        f = open(hfiles[-1],'r')
        tp=0
    elif isinstance(ant,list):
        if isinstance(ant[0],str) and isinstance(ant[1],str):
            tp=1
            print 'openning file ' + ant[0]
            f = open(ant[0],'r')
    
    if tp==-1:
        print "ant must be an integer from 1 to 13 or a list of 2 files, HPOL first and VPOL second"
        return
    
    lines = f.readlines()
    f.close()
    hp = []
    hv = []
    for i,line in enumerate(lines):
        if i > 0:
            vals = np.array(line.strip().split()).astype('float')
            hp.append(vals[0])
            hv.append(vals[4])
    hp = np.array(hp)
    hv = np.array(hv)
    pfit = np.polyfit(np.log(hv),hp,4)
    lhv = np.linspace(np.log(hv).min(),np.log(hv).max(),100)
    plt.plot(np.log(hv),hp,'o')
    plt.plot(lhv,np.polyval(pfit,lhv))
    for i,p in enumerate(pfit[::-1]):
        print 'HPOL.c{:d} = {:10.7f}'.format(i,p)
    
    if tp==0: 
        vfiles = glob.glob('Antenna '+str(ant)+' V*.txt')
        if len(vfiles) > 1: vfiles = np.sort(vfiles)
        print(vfiles[-1])
        f = open(vfiles[-1],'r')
    else:
        print 'openning file ' + ant[1]
        f = open(ant[1],'r')
    
    lines = f.readlines()
    f.close()
    vp = []
    vv = []
    for i,line in enumerate(lines):
        if i > 0:
            vals = np.array(line.strip().split()).astype('float')
            vp.append(vals[0])
            vv.append(vals[-1])
    vp = np.array(vp)
    vv = np.array(vv)
    pfit = np.polyfit(np.log(vv),vp,4)
    lvv = np.linspace(np.log(vv).min(),np.log(vv).max(),100)
    plt.plot(np.log(vv),vp,'*')
    plt.plot(lvv,np.polyval(pfit,lvv))
    plt.ylabel('Power [dBm]')
    plt.xlabel('Log Voltage [V]')
    if tp==0:
        plt.title('Antenna '+str(ant))
    else:
        pos = ant[0].rfind("/")
        print ant[0]
        print pos
        plt.title(ant[0][pos+1:pos+10].rstrip())
    for i,p in enumerate(pfit[::-1]):
        print 'VPOL.c{:d} = {:10.7f}'.format(i,p)
    plt.show()

def striplines(infile,pol):
    f = open(infile,'r')
    inlines = f.readlines()
    f.close()
    
    vals=[]
    vv = []
    pp = []
    for i,l in enumerate(inlines[1:]):
        l=l.strip()
        if l != "":
            v = np.array(l.split()).astype('float')
            vals.append(v)
            
    if pol == 'h':
        vals = [v for v in vals if v[2] != 31 and v[4] > 0.0]
        
    elif pol=='v':
        vals = [v for v in vals if v[2] != 31 and v[-1] > 0.0]
    
    tv = []
    for i,v in enumerate(vals):
        if i == 0:
            tv=[v]
        else:
            if v[2] == tv[-1][2] and v[3] == tv[-1][3]:
                tv.append(v)
            else:
                p = len(tv) // 2
                pp.append(np.median(np.array(tv)[:,4]))
                if pol == 'h':
                    vv.append(tv[p][4])
                elif pol == 'v':
                    vv.append(tv[p][-1])
                
                tv=[v]
    
    p = len(tv) // 2
    #pp.append(np.median(tv[0:len(tv)][0]))
    pp.append(np.median(np.array(tv)[:,4]))
    if pol == 'h':
        vv.append(tv[p][4])
    elif pol == 'v':
        vv.append(tv[p][-1])
    
    return np.array(pp),np.array(vv)            
        
def process_fem_cal(hfile,vfile,plottitle=None):
    #get HPOL FEM Calibration
    hp, hv = striplines(hfile,'h')
    #get VPOL FEM Calibratio
    vp, vv = striplines(vfile,'v')
    
    #Display the data
    print "HPOL DATA"
    for i in range(hp.size):
        print str(hp[i])," ",str(hv[i])
    
    print "VPOL DATA"
    for i in range(vp.size):
        print str(vp[i])," ",str(vv[i])
    
    #HPOL fit
    pfit = np.polyfit(np.log(hv),hp,4)
    lhv = np.linspace(np.log(hv).min(),np.log(hv).max(),100)
    plt.plot(np.log(hv),hp,'o')
    plt.plot(lhv,np.polyval(pfit,lhv))
    for i,p in enumerate(pfit[::-1]):
        print 'HPOL.c{:d} = {:10.7f}'.format(i,p)
    
    #VPOL fit
    pfit = np.polyfit(np.log(vv),vp,4)
    lvv = np.linspace(np.log(vv).min(),np.log(vv).max(),100)
    plt.plot(np.log(vv),vp,'*')
    plt.plot(lvv,np.polyval(pfit,lvv))
    for i,p in enumerate(pfit[::-1]):
        print 'VPOL.c{:d} = {:10.7f}'.format(i,p)
    
    #Add  plot labels
    if plottitle is not None:
         plt.title(plottitle)
    plt.ylabel('Power [dBm]')
    plt.xlabel('Log Voltage [V]')
    plt.show()
