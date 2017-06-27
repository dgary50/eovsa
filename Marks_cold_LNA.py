#!/usr/bin/env python

import sys
sys.path.insert(0,"/home/antctl/starburst/femTools/src")
bb = __import__('bb_worker')

bbw = bb.BBWorker()

# values for 300K

# amp SN560D
bbw.execute('LNA-ENABLE LV ON'.split())
bbw.execute('LNA-DRAIN LV 1.2'.split())
bbw.execute('LNA-GATE1 LV 0.63'.split())
bbw.execute('LNA-GATE2 LV 0.63'.split())
# amp SN568D
bbw.execute('LNA-ENABLE LH ON'.split())
bbw.execute('LNA-DRAIN LH 1.2'.split())
bbw.execute('LNA-GATE1 LH 1.06'.split())
bbw.execute('LNA-GATE2 LH 0.83'.split())
# amp SN40A158
bbw.execute('LNA-ENABLE HV ON'.split())
bbw.execute('LNA-DRAIN HV 1.2'.split())
bbw.execute('LNA-GATE1 HV -0.70'.split())
bbw.execute('LNA-GATE2 HV -0.75'.split())
#amp SN40A179
bbw.execute('LNA-ENABLE HH ON'.split())
bbw.execute('LNA-DRAIN HH 1.2'.split())
bbw.execute('LNA-GATE1 HH -0.6'.split())
bbw.execute('LNA-GATE2 HH -0.3'.split())

sys.exit()

