"""
    STARBURST ACC/FEANTA Server Runner
    Author: Lokbondo Kung
    Email: lkkung@caltech.edu
"""

import os
import sys
import time
import signal
import feanta_server
import pdu_worker
import brick_worker
import bb_worker
import cryostat_worker
from daemon import runner

def instantiate(pid_file, log_file):
    # Instantiate workers.
    pdu = pdu_worker.PDUWorker()
    brick = brick_worker.BrickWorker()
    bb = bb_worker.BBWorker()
    cryo = cryostat_worker.CryoWorker()

    # Instantiate server.
    server = feanta_server.ServerDaemon(pid_file)

    # Setup log file.
    server.set_log_file(log_file)

    # Link workers.
    server.link_worker(pdu)
    server.link_worker(brick)
    server.link_worker(bb)
    server.link_worker(cryo)

    return server, None

def start(server):
    # Start server.
    server_runner = runner.DaemonRunner(server)
    server_runner.do_action()

def __kill(pidfile):
    pid = 0
    with open(pidfile, 'r') as pid_file:
        pid = int(pid_file.readline())

    # Try killing the daemon process
    try:
        i = 0
        while 1:
            os.kill(pid, signal.SIGTERM)
            time.sleep(0.1)
            i = i + 1
            if i % 10 == 0:
                os.kill(pid, signal.SIGHUP)
    except OSError, err:
        err = str(err)
        if err.find("No such process") > 0:
            if os.path.exists(pidfile):
                os.remove(pidfile)
        else:
            print str(err)
            sys.exit(1)

def stop(pid_file):
    __kill(pid_file)

