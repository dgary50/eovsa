#!/usr/bin/env python

import subprocess, time

def rd_packetnum(source):
    res = subprocess.Popen(['ifconfig',source],stdout=subprocess.PIPE)
    for line in res.stdout.readlines():
        if line.find('RX packets') != -1:
            return int(line.split(':')[1].split(' ')[0])

def get_packetrates():
    start_pkt2 = rd_packetnum('eth2')
    start_pkt3 = rd_packetnum('eth3')
    time.sleep(1)
    end_pkt2 = rd_packetnum('eth2')
    end_pkt3 = rd_packetnum('eth3')
    return [end_pkt2 - start_pkt2, end_pkt3 - start_pkt3]

if __name__ == "__main__":
    print get_packetrates()
