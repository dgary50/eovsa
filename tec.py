#   2015-Jun-16  DG
#      FTP to cRIOs now requires a username and password
#
import urllib2
import numpy as np

def rd_teclog(crio=1):
    # List of expected commands, if all goes well...
    expected_commands = np.array(['> $R?\n']*97)
    for i in range(97): expected_commands[i] = expected_commands[i][0:4]+str(i)+expected_commands[i][4:]
    # List of expected responses to the commands
    expected_values = np.array(['+2.500000e+01\n', '+8.250000e+00\n', '+1.059999e-01\n',
                             '+4.289999e+02\n', '+2.000000e+00\n', '+3.000000e+00\n',
                             '+1.000000e+02\n', '+5.000000e-01\n', '+1.000000e+02\n',
                             '+5.000000e-02\n', '+5.500000e+00\n', '+2.500000e+00\n',
                             '+1.000000e-01\n', '134\n', '+5.000000e+00\n', '+5.000000e+00\n',
                             '1\n', '+2.000000e+01\n', '+8.000000e+00\n', '+4.000000e+00\n',
                             '+2.000000e+00\n', '+1.199999e+01\n', '+4.000000e+00\n', '1\n',
                             '+2.000000e+01\n', '+8.000000e+00\n', '+4.000000e+00\n', '+2.000000e+00\n',
                             '+1.199999e+01\n', '+4.799999e+01\n', '+0.000000e+00\n',
                             '+0.000000e+00\n', '+1.000000e+00\n', '+0.000000e+00\n',
                             '+1.000000e+00\n', '+1.000000e+00\n', '+0.000000e+00\n',
                             '+1.000000e+00\n', '+0.000000e+00\n', '+1.000000e+00\n',
                             '+0.000000e+00\n', '+1.000000e+00\n', '+1.059999e-01\n', '226\n',
                             '245\n', '+3.000000e+01\n', '+1.000000e+01\n', '+1.500000e+01\n',
                             '+1.000000e-01\n', '+2.000000e+00\n', '+1.000000e-01\n',
                             '+2.000000e+00\n', '+1.000000e-01\n', '+1.299999e+01\n', '+7.000000e+00\n',
                             '12\n', '4\n', '4\n', '4\n', '+1.396917e-03\n', '+2.378257e-04\n',
                             '+9.372652e-08\n', '+1.396917e-03\n', '+2.378257e-04\n',
                             '+9.372652e-08\n', '+1.396917e-03\n', '+2.378257e-04\n',
                             '+9.372652e-08\n', '+6.843508e-04\n', '+2.898552e-04\n',
                             '-8.177021e-13\n', '+8.000000e+01\n', '-4.000000e+01\n',
                             '+5.000000e+01\n', '-1.000000e+01\n', '+5.000000e+01\n',
                             '-1.000000e+01\n', '+6.000000e+01\n', '-1.000000e+01\n', 
                             '+7.593999e+02\n', '+3.057699e+03\n', '+2.987579e+04\n', 
                             '+7.593999e+02\n', '+3.057699e+03\n', '+0.000000e+00\n', 
                             '+7.593999e+02\n', '+0.000000e+00\n', '+2.987579e+04\n', 
                             '+2.965139e+03\n', '+2.883676e+04\n', '+7.821895e+04\n', '351\n', '67\n',
                             '+8.000000e+00\n', '300\n', '200\n', '21006\n'])
    try:
        userpass = 'admin:observer@'
        ftpadr = 'ftp://'+userpass+'crio'+str(crio)+'.solar.pvt/tec.txt'
    except:
        print 'Cannot create FTP address from crio:',crio
        return {}
    try:
        f = urllib2.urlopen(ftpadr,timeout=0.5)
        lines = np.array(f.readlines())
        f.close()
    except:
        print 'Cannot retrieve TEC log file from',ftpadr
        return {}
    idx = np.array([],'int')
    for i in range(len(lines)):
        if len(lines[i].split()) == 2: 
            idx = np.append(idx,i)
    if len(idx) == 0:
        print 'No correct commands found.  Possibly may have to issue TEC$A command to turn off broadcast?'
    rd_cmds = lines[idx[0:96]]
    responses = lines[idx[0:96]+1]
    # Compare responses with expected values (output is an array of True/False)
    ok = np.in1d(responses,expected_values)
    bad = np.where(~ok)[0]
    for i in bad:
        print rd_cmds[i].strip(),'gave bad response',responses[i].strip(),'Expected:',expected_commands[i].strip(),'=',expected_values[i].strip()
    return {}
