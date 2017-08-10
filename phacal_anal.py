''' Routine to diagnose all phacal observations data based on the provided refcal.
'''
#  2017-06-26 DG
#    Added phacal_anal() routine to plot phase difference wrt refcal, and return
#    the fitted phase slopes and offsets
#  2017-06-28 SJ

import read_idb as ri
from util import Time, ant_str2list
import numpy as np
import matplotlib.pyplot as plt
from glob import glob
import matplotlib.dates as mdates
import os
# import pcapture2 as p
import pdb
import refcal_anal as ra
import cal_header as ch


def phacal_anal(refcal, minsnr=0.7, doplot_refcal=False, verbose=False):
    '''Fit the phase difference between a phase calibration (or another refcal)
       and the reference calibration.  Returns the phase slopes and, optionally,
       the offsets, along with the relevant times, as a dictionary.

       Returns:
       pfit    A python dictionary with the following keys:
           't_pha' : Central time of the phase calibration data
           't_ref' : Central time of the reference calibration data
           'poff'  : Phase offsets (zero if fitoffsets=False) for each antenna, polarization
           'pslope': Phase slope for each antenna, polarization
           'flag'  : Indication of bad data for each antenna, polarization (0 = good, 1 = bad)
           'phacal': Actual instance of phase calibration dictionary, containing frequency-dependent
                       amplitudes, phases, etc.
    '''
    if isinstance(refcal, Time):
        refcal = ra.sql2refcal(refcal)
    elif isinstance(refcal, dict):
        pass
    else:
        raise ValueError('The refcal must be a refcal dict or timestamp')
    timestamp = refcal['timestamp']
    dhr = timestamp.LocalTime.utcoffset().total_seconds() / 60 / 60 / 24
    btime = Time(np.fix(timestamp.mjd + dhr) - dhr, format='mjd')
    etime = Time(btime.mjd + 1, format='mjd')
    trange = Time([btime, etime])
    out = ra.rd_refcal(trange, projid='PHASECAL')
    out = ra.unrot_refcal(out)
    scanidx = []
    for ll in range(len(out['scanlist'])):
        if out['tstlist'][ll].mjd <= timestamp.mjd <= out['tedlist'][ll].mjd:
            pass
        else:
            scanidx.append(ll)
    # ra.graph(out, bandplt=[5, 7, 9, 21])
    nscanidx = len(scanidx)
    scanidx_sql = []
    ct_accept = 0
    print '{} PHASECAL scan(s) found...'.format(nscanidx)
    for idx, ll in enumerate(scanidx):
        print 'processing PHASECAL {}/{} ...'.format(idx + 1, nscanidx)
        phacal = ra.phase_diff(ra.refcal_anal(out, scanidx=[ll], minsnr=minsnr, doplot=doplot_refcal), refcal=refcal)
        prompt = ''
        while not (prompt.lower() in ['y', 'n']):
            prompt = raw_input('Do you want to accept this phacal results? [y/n]')
        if prompt.lower() == 'n':
            print 'PHASECAL {}/{} abort ...'.format(idx + 1, nscanidx)
            scanidx_sql.append({'id': ll, 'status': 'rejected'})
        elif prompt.lower() == 'y':
            print 'PHASECAL {}/{} sending to SQL datebase'.format(idx + 1, nscanidx)
            scanidx_sql.append({'id': ll, 'status': 'accepted'})
            ct_accept += 1
            ch.phacal2sql(phacal)
        plt.close('all')
    print '{} out of {} PHASECAL results were sent to SQL database.'.format(ct_accept, nscanidx)
    print 'PHASECAL: time/            source/          status'
    for idx, ll in enumerate(scanidx_sql):
        print '{0}: {1}~{2}/ {3}/ {4}'.format(ll['id'], out['tstlist'][ll['id']].isot[:-4],
                                              out['tedlist'][ll['id']].isot[:-4], out['srclist'][ll['id']],
                                              ll['status'])
