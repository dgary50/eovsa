import eovsa_cat as ec
import copy
satnames = array(['O3B FM']*11)
sats = []
for i in range(2,13):
    sats.append(satnames[i-2]+str(i))
aa = ec.eovsa_array_with_cat()
t = aa.get_jultime()
ot = copy.copy(t)
for satname in sats:
    t = copy.copy(ot)
    aa.set_jultime(ot)
    for i in range(144):
        t += 300./86400.
        aa.set_jultime(t)
        sat = aa.cat[satname]
        sat.compute(aa)
        if sat.alt > 10.*pi/180.:
            print aa.date,satname,sat.ra,sat.dec

