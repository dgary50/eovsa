import subprocess, glob
from util import Time
def dpp_status():
    datstr = Time.now().iso[:10].replace('-','')
    out1 = sort(glob.glob('/data1/IDB/IDB'+datstr+'*'))[-5:]
    command = 'python dpp_eth_mon.py'
    out2 = subprocess.check_output(command.split(),stderr=subprocess.STDOUT).split('\n')[0]
    command = 'ps -C dppxmp4'
    ps = subprocess.Popen(command.split(),stdout=subprocess.PIPE)
    out3 = ps.communicate()[0].split('\n')[1]
    return out1, out2, out3