'''
   Module for plotting the median of the front-end RF detector voltages
   from the stateframe SQL database, as a crude flare monitor'''
#
# History:
#   2024-Mar-17  DG
#      Entirely new replacement for old flare_monitor.py routine.  This one
#      reads the realtime flaretest files for a given date and finds flares,
#      indicating the flare times with a black step curve.
#
import numpy as np
from util import Time

def get_projects(t, nosql=False):
    ''' Read all projects from SQL for the current date and return a summary
        as a dictionary with keys Timestamp, Project, and EOS (another timestamp)
    '''
    if nosql == True:
        return get_projects_nosql(t)
    import dbutil
    # timerange is 12 UT to 12 UT on next day, relative to the day in Time() object t
    trange = Time([int(t.mjd) + 12./24,int(t.mjd) + 36./24],format='mjd')
    tstart, tend = trange.lv.astype('str')
    cursor = dbutil.get_cursor()
    mjd = t.mjd
    # Get the project IDs for scans during the period
    verstrh = dbutil.find_table_version(cursor,trange[0].lv,True)
    if verstrh is None:
        print 'No scan_header table found for given time.'
        return {}
    query = 'select Timestamp,Project from hV'+verstrh+'_vD1 where Timestamp between '+tstart+' and '+tend+' order by Timestamp'
    projdict, msg = dbutil.do_query(cursor, query)
    if msg != 'Success':
        print msg
        return {}
    elif len(projdict) == 0:
        # No Project ID found, so return data and empty projdict dictionary
        print 'SQL Query was valid, but no Project data were found.'
        return {}
    projdict['Timestamp'] = projdict['Timestamp'].astype('float')  # Convert timestamps from string to float
    for i in range(len(projdict['Project'])): projdict['Project'][i] = projdict['Project'][i].replace('\x00','')
    projdict.update({'EOS':projdict['Timestamp'][1:]})
    projdict.update({'Timestamp':projdict['Timestamp'][:-1]})
    projdict.update({'Project':projdict['Project'][:-1]})
    cursor.close()
    return projdict

def get_projects_nosql(t):
    ''' Read all projects from FDB file for the current date and return a summary
        as a dictionary with keys Timestamp, Project, and EOS (another timestamp)
    '''
    import dump_tsys as dt
    # timerange is 12 UT to 12 UT on next day, relative to the day in Time() object t
    trange = Time([int(t.mjd) + 12./24,int(t.mjd) + 36./24],format='mjd')
    tstart = t.iso[2:10].replace('-','')+'120000'
    t2 = Time(int(t.mjd)+1, format='mjd')
    tend = t2.iso[2:10].replace('-','')+'120000'
    fdb = dt.rd_fdb(t)
    fdb2 = dt.rd_fdb(t2)
    if fdb == {}:
        # No FDB file found, so return empty project dictionary
        print 'No Project data [FDB file] found for the given date.'
        return {}
    if fdb2 == {}:
        pass
    else:
        #  Concatenate the two dicts into one
        if fdb['FILE'][-1] == '':
            fdb = dict([(k, np.concatenate((fdb.get(k,[])[:-1],fdb2.get(k,[])))) for k in set(fdb)|set(fdb2)])
        else:
            fdb = dict([(k, np.concatenate((fdb.get(k,[]),fdb2.get(k,[])))) for k in set(fdb)|set(fdb2)])  
    # Get "good" indexes for times between 12 UT on date and 12 UT on next date
    gidx, = np.where(np.logical_and(fdb['SCANID']>tstart,fdb['SCANID']<tend))        
    scanid,idx = np.unique(fdb['SCANID'][gidx],return_index=True)
    sidx = gidx[idx]   # Indexes into fdb for the start of each scan
    eidx = np.concatenate((sidx[1:]-1,np.array([gidx[-1]])))   # Indexes into fdb for the end of each scan
    # Get the project IDs for scans during the period
    projdict = {'Timestamp':fdb['ST_TS'][sidx].astype(float),
                'Project':fdb['PROJECTID'][sidx],
                'EOS':fdb['EN_TS'][eidx].astype(float)}
    return projdict
    
def rd_RT(datstr=None, nsigma=3, nflare=10, nbgnd=60, debug=False):
    ''' Reads all flaretest real-time files available for the given date
        provides as datstr (in the form 'yyyy-mm-dd'), and returns a
        dictionary of background-subtracted data, times (in plot_date format)
        and a flare flag.
    '''
    import glob
    from copy import copy
    
    if datstr is None:
        datstr = Time.now().iso
    yyyy = datstr[:4]
    mm = datstr[5:7]
    dd = datstr[8:10]
    mjd = int(Time(datstr).mjd)+0.5
    mjd2 = mjd + 0.6
    datstr2 = Time(mjd2,format='mjd').iso[:10]
    yyyy2 = datstr2[:4]
    mm2 = datstr2[5:7]
    dd2 = datstr2[8:10]
    fstr = '/data1/eovsa/fits/FTST/'+yyyy+'/'+mm+'/flaretest_'+yyyy[2:]+mm+dd+'*.txt'
    files = sorted(glob.glob(fstr))
    fstr2 = '/data1/eovsa/fits/FTST/'+yyyy2+'/'+mm2+'/flaretest_'+yyyy2[2:]+mm2+dd2+'*.txt'
    files += sorted(glob.glob(fstr2))
    if len(files) == 0:
        print('No files for this date')
        return {}
    filelist = []
    for file in files:
        s = file[-16:-4]
        mjdfile = Time('20'+s[:2]+'-'+s[2:4]+'-'+s[4:6]+' '+s[6:8]+':'+s[8:10]+':'+s[10:]).mjd
        if mjdfile > mjd and mjdfile < mjd2:
            filelist.append(file)
    if len(filelist) == 0:
        print('No files during the solar day for this date')
        return {}
            
    data = []
    t = []
    bgnd = [[],[],[]]
    flareflag = []
    
    if debug:
        bk_vals = []
        sigma_vals = []
        numvotes = []
        inflr = []
        bgnd_evol = []
        
    inflare = 0
#    nflare = 10
#    nsigma = 3
    for file in filelist:
        fh = open(file,'r')
        # Read the most recent lines of the file in a loop:
        lines = fh.readlines()
        fh.close()
        source = lines[0].split()[1].replace('\x00','')
        project = lines[1].split()[1].replace('\x00','')
        scanid = lines[2].split()[1].replace('\x00','')
        if project == 'NormalObserving':
            #import pdb; pdb.set_trace()
            for k,line in enumerate(lines[8:]):
                cols = line.split()
                datstr3 = cols[0][:4]+'-'+cols[0][4:6]+'-'+cols[0][6:8]
                timstr = cols[1][0:2]+':'+cols[1][2:4]+':'+cols[1][4:]
                t.append(Time(datstr3+' '+timstr).plot_date)
                # Skip if any of the data points are zero
                if np.float64(cols[3]) != 0 and np.float64(cols[4]) != 0 and np.float64(cols[5]) != 0:
                    if len(bgnd[0]) < nbgnd:
                        # Add this point to background
                        bgnd[0].append(np.float64(cols[3]))
                        bgnd[1].append(np.float64(cols[4]))
                        bgnd[2].append(np.float64(cols[5]))
                        flareflag.append(0)
                    else:
                        # Do this only if the background has at least nbgnd pts
                        # Get votes for this sample
                        nvotes = 0   # Number of "votes" for being in a flare
                        for i in range(3):
                            bk = np.mean(np.array(bgnd[i]))
                            sigma = np.std(np.array(bgnd[i]))
                            if debug:
                                bk_vals.append(bk)
                                sigma_vals.append(sigma)
                            # Determine if we should be in the flare state
                            if (np.float64(cols[3+i]) - bk) > nsigma*sigma:
                                nvotes += 1
                        if debug:
                            numvotes.append(nvotes)
                        if k < 120:
                            # Ignore votes within 2 min of the start of a scan
                            nvotes = 0
                        if nvotes <=1 :
                            # This point is not consistent with being in a flare, so decrement inflare (or zero)
                            inflare = max(inflare-1,0)
                            #if inflare == 0:
                            #    # Definitely not in a flare, so add this point to the background
                            for i in range(3):
                                # Remove first background point and add new point to end
                                bgnd[i].pop(0)
                                bgnd[i].append(np.float64(cols[3+i]))
                            flareflag.append(0)
                            #else:
                            #    # This sample is not consistent with being in a flare, but the flare window is not yet expired
                            #    # so go ahead and add it
                            #    flareflag.append(1)
                        else:
                            # This point is consistent with being in a flare, so increment inflare (max = nflare)
                            # and save data for level check
                            if inflare == 0:
                                # This is the first flare point, so initialize flaredat
                                flaredat = [np.float64(cols[3:6])]
                                printflare = True
                            flaredat.append(np.float64(cols[3:6]))
                            inflare = min(inflare+1,nflare)
                            if inflare == nflare:
                                if len(flaredat) == nbgnd:
                                    # The algorithm thinks it has found a flare, but check that it is not just a level change
                                    print(Time(t[-1],format='plot_date').iso,'flaredat has',len(flaredat),'points')
                                    fvotes = 0
                                    for i in range(3):
                                        # Compare background standard deviation with standard deviation of nbgnd pts of flare data.
                                        # We need at least two votes of f_std greater than 2.5*b_std
                                        b_std = np.std(np.array(bgnd[i]))
                                        f_std = np.std(np.array(flaredat)[nflare:,i])
                                        if f_std > 2.5*b_std:
                                            fvotes += 1
                                        print('b_mean:',np.mean(np.array(bgnd[i])), 'f_mean:', np.mean(np.array(flaredat)[nflare:,i]),
                                              'b_std:', np.std(np.array(bgnd[i])),  'f_std:',  np.std(np.array(flaredat)[nflare:,i]))
                                    print('Flare votes:',fvotes)
                                    if fvotes <= 1:
                                        # Not a flare, so remove earlier flare flags and set background to "flare" points
                                        for i in range(nbgnd):
                                            flareflag[-i] = 0
                                        for i in range(3):
                                            bgnd[i] = np.array(flaredat)[:,i].tolist()
                                        inflare = 0
                                        flareflag.append(0)
                                    else:
                                        # Definitely in a flare, so set flare flag
                                        flareflag.append(nvotes)
                                    printflare = False 
                                else:
                                    # Definitely in a flare, so set flare flag
                                    flareflag.append(nvotes)
                            else:
                                if flareflag[-1] == 1 and inflare > nflare/2:
                                    # We were in a flare state and half the flare window is not yet expired, 
                                    # so even though this point may not be in the flare, go ahead and add it
                                    flareflag.append(nvotes)
                                else:
                                    # Half the flare window is expired, so consider this point as not in a flare
                                    flareflag.append(0)
                    if debug:
                        inflr.append(inflare)
                    for i in range(3):
                        bk = np.mean(np.array(bgnd[i]))
                        data.append(np.float64(cols[3+i])-bk)
                else:
                    # One or more of the datapoints was 0, so set all to 0, but remain agnostic about flareflag
                    data.append(0)
                    data.append(0)
                    data.append(0)
                    if debug:
                        bk_vals.append(0)
                        bk_vals.append(0)
                        bk_vals.append(0)
                        sigma_vals.append(0)
                        sigma_vals.append(0)
                        sigma_vals.append(0)
                        numvotes.append(0)
                        inflr.append(0)
                    if len(flareflag) == 0:
                        flareflag.append(0)
                    else:
                        flareflag.append(flareflag[-1])
                if debug:
                    bgnd_evol.append(copy(bgnd[0]))

    if len(data) == 0:
        print('No solar scan files for this date')
        return {}
    # Remove any flare flags of less than nbgnd s duration
    f_on = np.clip(np.array(flareflag),0,1)
    transitions_up, = np.where(f_on - np.roll(f_on,1) == 1)
    transitions_down, = np.where(f_on - np.roll(f_on,1) == -1)
    tup = []
    tdown = []
    for i, tran in enumerate(transitions_up):
        if np.sum(f_on[tran:tran+nbgnd-nflare]) == nbgnd-nflare:
            tup.append(transitions_up[i])
            tdown.append(transitions_down[i])
    fgood = np.zeros_like(f_on)
    for i in range(len(tup)):
        fgood[tup[i]:tdown[i]] = 1

    data = np.array(data)
    nt = len(t)
    data.shape = (nt,3)
    # Get indexes where data are zero (indicates dropouts)
    bad, = np.where(data[:,0] == 0)
    for i in range(3):
        data[bad,i] = np.nan
    
    # fig, ax = plt.subplots(1,1)
    # fig.suptitle('EOVSA data for '+datstr)
    # fig.set_figheight(3)
    # fig.set_figwidth(14)
    # ax.plot_date(t,data[:,0],'-',color='C0')
    # ax.plot_date(t,data[:,1],'-',color='C1')
    # ax.plot_date(t,data[:,2],'-',color='C2')
    # maxdata = np.nanmax(data[:,0])
    # ax.set_ylim(-100,max(maxdata,320000))
    # ax.plot_date(t,np.array(flareflag)*fgood*100000,'-',color='k')
    if debug:
        return {'time':np.array(t), 'data':data, 'flareflag':np.array(flareflag)*fgood, 'bk_vals':bk_vals, 'sigma_vals':sigma_vals, 'inflare':inflr, 'nvotes':numvotes, 'bgnd_evol':bgnd_evol}
    else:
        return {'time':np.array(t), 'data':data, 'flareflag':np.array(flareflag), 'finalflag':np.array(flareflag)*fgood}

if __name__ == "__main__":
    ''' For non-interactive use, use a backend that does not require a display
        Usage python /common/python/current/flare_monitor.py "2014-12-20" <skip>
        If optional argument skip is given, the time-consuming creation of the
        xdata spectrum (XSP file) is skipped.
    '''
    import glob, shutil
    import matplotlib, sys, util
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    known = ['GAIN','PHAS','SOLP']  # known calibration types (first 4 letters)
    t = Time.now()
    if len(sys.argv) >= 2:
        try:
            t = Time(sys.argv[1])
        except:
            print 'Cannot interpret',sys.argv[1],'as a valid date/time string.'
            exit()
    print t.iso[:19],': ',

    if t.mjd % 1 < 0.5:
        t = Time(t.mjd - 0.5, format='mjd')
    projdict = get_projects_nosql(t)
    ut = [Time(projdict['Timestamp'],format='lv').plot_date[0]]*2
    out = rd_RT(t.iso)
    f, ax = plt.subplots(1,1)
    f.set_size_inches(10,3)
    ax.set_ylim(-100,320000)
    if out != {}:
        ax.plot_date(out['time'],out['finalflag']*100000,'-',color='k')
        ax.plot_date(out['time'],out['data'][:,0],'-',color='C0')
        ax.plot_date(out['time'],out['data'][:,1],'-',color='C1')
        ax.plot_date(out['time'],out['data'][:,2],'-',color='C2')
        maxdata = np.nanmax(out['data'][:,0])
        ax.set_ylim(-100,max(maxdata,320000))
    ax.set_xlabel('Time [UT]')
    ax.set_ylabel('Flare Detector [arb. units]')
    ax.set_title('EOVSA Flare Monitor for '+t.iso[:10])
    ax.set_xlim(int(ut[0])+13/24.,int(ut[0])+26/24.)  # Time plot ranges from 13 UT to 02 UT
    if projdict == {}:
        print 'No annotation can be added to plot for',t.iso[:10]
    else:
        nscans = len(projdict['Project'])
        SOS = Time(projdict['Timestamp'],format='lv').plot_date
        EOS = Time(projdict['EOS'],format='lv').plot_date
        yran = np.array(ax.get_ylim())
        for i in range(nscans):
            uti = SOS[i]*np.array([1.,1.])
            plt.plot_date(uti,yran,'g',lw=0.5)
            if projdict['Project'][i] == 'NormalObserving' or projdict['Project'][i] == 'Normal Observing':
                ax.text(uti[0],yran[1]*0.935,'SUN',fontsize=8, clip_on=True)
            elif projdict['Project'][i] == 'None':
                ax.text(uti[0],yran[1]*0.975,'IDLE',fontsize=8, clip_on=True)
            elif projdict['Project'][i][:4] == 'GAIN':
                ax.text(uti[0],yran[1]*0.955,'GCAL',fontsize=8, clip_on=True)
            elif projdict['Project'][i] == 'SOLPNTCAL':
                ax.text(uti[0],yran[1]*0.955,'TPCAL',fontsize=8, clip_on=True)
            elif projdict['Project'][i] == 'PHASECAL':
                ax.text(uti[0],yran[1]*0.955,'PCAL',fontsize=8, clip_on=True)
            else:
                ax.text(uti[0],yran[1]*0.975,projdict['Project'][i],fontsize=8)
        if len(projdict['EOS']) == nscans:
            for i in range(nscans):
                uti = EOS[i]*np.array([1.,1.])
                plt.plot_date(uti,yran,'r--',lw=0.5)
                uti = np.array([SOS[i],EOS[i]])
                if projdict['Project'][i] == 'NormalObserving':
                    plt.plot_date(uti,yran[1]*np.array([0.93,0.93]),ls='-',marker='None',color='#aaffaa',lw=2,solid_capstyle='butt')
                elif projdict['Project'][i][:4] in known:
                    plt.plot_date(uti,yran[1]*np.array([0.95,0.95]),ls='-',marker='None',color='#aaaaff',lw=2,solid_capstyle='butt')
                else:
                    plt.plot_date(uti,yran[1]*np.array([0.97,0.97]),ls='-',marker='None',color='#ffaaaa',lw=2,solid_capstyle='butt')
    datstr = t.iso[:10].replace('-','')
    plt.savefig('/common/webplots/flaremon/FLM'+datstr+'.png',bbox_inches='tight')
    plt.close(f)
    print 'Plot written to /common/webplots/flaremon/FLM'+datstr+'.png'
    # Copy the most recent two files to fixed names so that the web page can find them.
    flist = np.sort(glob.glob('/common/webplots/flaremon/XSP20*'))
    shutil.copy(flist[-1],'/common/webplots/flaremon/XSP_latest.png')
    shutil.copy(flist[-2],'/common/webplots/flaremon/XSP_later.png')

