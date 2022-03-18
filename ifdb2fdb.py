import sys
from util import Time
import os

#This program converts IFDB files back into FDB files.
#You need to run this from the folder that contains the IFDB file.
#It  will create a directory called FDB that contains the FDB files.

t = Time(sys.argv[1])
endt = Time(sys.argv[2])

if not os.path.isdir('FDB'):
    print "Creating FDB Folder"
    os.mkdir('FDB')
    
while t.jd < endt.jd:
    fname = "IFDB" + t.iso[:10].replace('-','') + '.txt'
    if os.path.isfile(fname):
        f1 = open(fname,'r')
        lines = []
        for line in f1:
            lines.append(line)
        f1.close()
        
        fname = 'FDB/' + fname[1:]
        f2 = open(fname,'w')
        f2.write(' FILE: SCANID: SOURCEID: PROJECTID: ST_ACC: EN_ACC: ST_SEC: ST_USEC: EN_SEC: EN_USEC: ST_TS: EN_TS: \n')
        for line in lines[1:]:
            f2.write(line[:18]+'  '+line[18:76]+'\n')
            f2.write('           0           0  ' + str(int(float(line[77:89])-2082844800)) + '      0  ' + str(int(float(line[90:102])-2082844800)) + '      0  ' + line[77:88] + '  ' + line[90:101] + '\n')
        f2.close()
    t = Time(t.jd+1.0, format='jd')
    

