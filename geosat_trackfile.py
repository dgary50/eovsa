from eovsa_lst import *

def geosat_trackfile(name, hadeg, decdeg):
    # Generates a 1-hour track file for a geosynchronous satellite
    # at location hadeg, decdeg (degrees)
    d = datime()
    mjd = d.get()
    hadeg = -3.80
    decdeg = -5.895
    ha = RA_Angle(hadeg*0.0174533)
    lines = []
    for i in range(60):
        mjd1 = mjd+i*60/86400.
        d.set(mjd1)
        lst = eovsa_lst(d)
        ra = lst - ha
        lines.append(str(int(ra.get()*10000/0.0174533))+' '+str(int(decdeg*10000))+' '
                    +str(int(d.get()))+' '+str(int((d.get() % 1)*86400000)))

    f = open(name+'_tab.radec','w')
    for line in lines:
        f.write(line+'\n')
    f.close()
