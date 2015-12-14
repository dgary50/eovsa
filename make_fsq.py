def make_fsq(band):

    f = open('band'+str(band)+'.fsq','w')
    print 'band'+str(band)+'.fsq'
    f.write('LIST:DWELL '+'1,'*34+'1\n')
    f.write('LIST:SEQUENCE '+str(band)+'\n')
    f.close()
