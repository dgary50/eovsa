import stateframe
from read_xml import *
import copy as cp
import os
from pylab import *

def sf_log_temp():
    accini = stateframe.rd_ACCfile()
    f = open('/tmp/stateframe.log','rb')
    nrecs = os.path.getsize('/tmp/stateframe.log')/accini['binsize']/10
#    t1fem = zeros(nrecs)
#    t1tec = zeros(nrecs)
#    t7fem = zeros(nrecs)
#    t7tec = zeros(nrecs)
#    t8fem = zeros(nrecs)
#    t8tec = zeros(nrecs)
#    tair = zeros(nrecs)
    for i in range(nrecs):
        try:
            data = f.read(accini['binsize']*10)
        except:
            break
        keys, template, fmt = cp.deepcopy(accini['xml'])
        sf = decode_xml_data(data, keys, template, fmt)

        plot(i,sf['Antenna'][0]['Frontend']['FEM']['Temperature'],'b.')
        plot(i,sf['Antenna'][0]['Frontend']['TEC']['Temperature'],'r.')
        plot(i,sf['Antenna'][6]['Frontend']['FEM']['Temperature'],'g.')
        plot(i,sf['Antenna'][6]['Frontend']['TEC']['Temperature'],'.',color='orange')
#        plot(i,sf['Antenna'][7]['Frontend']['FEM']['Temperature'],'.',color='yellow')
        plot(i,sf['Antenna'][7]['Frontend']['TEC']['Temperature'],'.',color='yellow')
        plot(i,(sf['Schedule']['Data']['Weather']['Temperature'] - 32)*5./9,'k.')

#    plot(t1fem)
#    plot(t1tec)
#    plot(t7fem)
#    plot(t7tec)
#    plot(t8fem)
#    plot(t8tec)
#    plot(tair)
