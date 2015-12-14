from util import Time
import dbutil as db
import numpy as np

hcurve = []
vcurve = []
t = Time(['2015-10-28 6:17','2015-10-28 6:20','2015-10-28 6:23','2015-10-28 6:30','2015-10-28 6:33','2015-10-28 6:36','2015-10-28 6:39']).lv.astype('int')
for lv in t:
    query = 'select Timestamp,Ante_Fron_FEM_HPol_Voltage from fv61_vD15 where (I15 % 15) = 13 and Timestamp between '+str(lv-90)+' and '+str(lv+90)+' order by Timestamp'
    hc, msg = db.do_query(cursor,query)
    hcurve.append(hc)
    query = 'select Timestamp,Ante_Fron_FEM_VPol_Voltage from fv61_vD15 where (I15 % 15) = 13 and Timestamp between '+str(lv-90)+' and '+str(lv+90)+' order by Timestamp'
    vc, msg = db.do_query(cursor,query)
    vcurve.append(vc)

f,ax = subplots(2,1)
hlabel = ['175 mm','150 mm','125 mm','100 mm','75 mm','50 mm','25 mm']
for i,h in enumerate(hcurve):
    x = (h['Timestamp']-t[i])*np.cos(13*pi/180)
    ax[0].plot(x,h['Ante_Fron_FEM_HPol_Voltage'],label=hlabel[i])
for i,v in enumerate(vcurve):
    x = (v['Timestamp']-t[i])*np.cos(13*pi/180)
    ax[1].plot(x,v['Ante_Fron_FEM_VPol_Voltage'],label=hlabel[i])
ax[0].set_ylim(1.2,2.2)
ax[0].legend(fontsize=10)
ax[0].set_xlim(-50,50)
ax[1].set_ylim(1.2,2.2)
ax[1].set_xlim(-50,50)
ax[1].legend(fontsize=10)
suptitle('27-m Drift Scans on Moon vs. Focus',fontsize='18')
ax[0].set_ylabel('HPol Voltage')
ax[1].set_ylabel('VPol Voltage')
ax[1].set_xlabel('Arcmin from Nominal Center')
ax[0].plot([-16.5,-16.5],[1.2,2.2],color='blue')
ax[0].plot([16.5,16.5],[1.2,2.2],color='blue')
ax[1].plot([-16.5,-16.5],[1.2,2.2],color='blue')
ax[1].plot([16.5,16.5],[1.2,2.2],color='blue')
