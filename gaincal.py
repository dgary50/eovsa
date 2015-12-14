# 
# RG
# This program is designed to help display and analyze data
#  collected on the attenuators during gain calibration
#  (FEATTNTEST) scans so that the attenuators may be 
#  adjusted correctly, and to learn the impact of the 
#  of the attenuators and there gain on data that 
#  has been collected. 
#
# RG
# June 30, 2015
# More updates, can now accomodate gaincal data with errors. 
#
# RG
# July 3, 2015
# Now greatly improved. show_dB_ratio has been completely rewritten. 
#  It now creates tables of attenuations.
#
# RG
# July 5, 2015
# Fixed some problems in code regarding how it handles an 
#  attenuation of 0, which is often read when an error 
#  briefly occurs. Also, find_gaincal now returns the 
#  timerange or timeranges of gaincal scans on a given
#  day. 
#
# RG
# July 15, 2015
# I have many many additions and changes since last update. The program
#  now has feature to find and graph the "tsys" with corrections for the n
#  and attenuation states. There is a new schedule item (GAINCALTEST.scd)
#  that is used by the program as well.
#
# RG
# July 21, 2015
# The program now finds and corrects for the attenuations of both the first 
#  and second front end attenuators. The attenuations for the first front
#  end attenuator are found by running, and correctly using the data, from
#  a FEATTNTEST2, as opposed to the FEATTNTEST, which is used to find the 
#  attenuations of the second front end attenuator. 
#
# *note: All attenuations, calculation, plots, etc. found in this program are
#  predicated on the assumption that the "0" attenuator is actually 0, or 
#  extremely close to it*

import spectrogram_fit as sp
from util import Time
import util
import numpy as np
import matplotlib.pylab as plt
import subprocess
import dbutil
import stateframedef
import astropy.table as tbl 
from matplotlib.font_manager import FontProperties
    
def show_dB_ratio(trange, test='FEATTNTEST'):
    # This routine creates a list of the change in decibels 
    #  for the change between the attenuatorsfor each of 
    #  the four cycles
    if test == 'FEATTNTEST':
        cycles = 4
        attenuator = 1
    else:
        if test == 'FEATTNTEST2':
            cycles = 1
            attenuator = 0
        else:
            print 'Please input valid attenuation test'
    s = sp.Spectrogram(trange)
    s.docal = False
    s.dosub = False
    s.domedian = False
    tsys, std = s.get_data()
    stateidx, res = get_state_idx(trange, cycles=cycles, attenuator=attenuator)

    all_tsys_vals = []
    for a in range(8):
        a_vals = []
        for p in range(2):
            p_vals = []
            for f in range(tsys.shape[2]):
                f_vals = []
                for i in range(6):
                    tsys_vals = []
                    for j in range(stateidx.shape[1]):
                        _vals = []
                        for vals in stateidx[i][j]:
                            _vals.append(tsys[a, p, f, vals])
                        tsys_vals.append(_vals)
                    f_vals.append(tsys_vals)
                p_vals.append(f_vals)
            a_vals.append(p_vals)                       
        all_tsys_vals.append(a_vals)
    all_tsys_vals = np.array(all_tsys_vals)

    all_avgs = []
    for a in range(8):
        a_avgs = []
        for p in range(2):
            p_avgs = []
            for f in range(tsys.shape[2]):
                f_avgs = []
                for i in range(6):
                    tsys_avgs = []
                    for j in range(stateidx.shape[1]):
                        tsys_avgs.append(sum(all_tsys_vals[a, p, f, i, j])/len(all_tsys_vals[a, p, f, i, j]))
                    f_avgs.append(tsys_avgs)
                p_avgs.append(f_avgs)
            a_avgs.append(p_avgs)                       
        all_avgs.append(a_avgs)
    all_avgs = np.array(all_avgs)
    
    noise_level = []
    for ant in range(tsys.shape[0]):
	pol_noise = []
        state, = np.where(np.logical_and(res['Ante_Fron_FEM_HPol_Atte_Second'][ant::15].astype('int') == 31,res['Ante_Fron_FEM_HPol_Atte_First'][ant::15].astype('int') == 31))
        for pol in range(tsys.shape[1]):
            freq_noise = []	    
	    for freq in range(tsys.shape[2]):
	        avgs_noise = []
                for cycles in range(stateidx.shape[1]):
                    avg_noise = []
		    for index in state:
                        if index <= stateidx[:, cycles][-1][-1]:
                            try:
		                if np.logical_and(tsys[ant, pol, freq, index] <= 0.005, index < tsys.shape[3]):
		                    avg_noise.append(tsys[ant, pol, freq, index])
                            except:
                                pass
                    avgs_noise.append(np.average(avg_noise))
		freq_noise.append(avgs_noise)
	    pol_noise.append(freq_noise)
	noise_level.append(pol_noise)
    noise_level = np.array(noise_level)

    all_dBratios = []
    for a in range(8):
        a_dBratios = []
        for p in range(2):
            p_dBratios = []
            for f in range(tsys.shape[2]):
                f_dBratios = []
                for j in range(stateidx.shape[1]):               
                    tsys_dBratios = []
                    for i in range(5):
                        tsys_dBratio = 10*np.log10((all_avgs[a, p, f, i, j] - noise_level[a, p, f, j])/(all_avgs[a, p, f, i+1, j] - noise_level[a, p, f, j]))
                        tsys_dBratios.append(tsys_dBratio)             
                    f_dBratios.append(tsys_dBratios)
                p_dBratios.append(f_dBratios)
            a_dBratios.append(p_dBratios)                       
        all_dBratios.append(a_dBratios)
    all_dBratios = np.array(all_dBratios)
    return all_dBratios, noise_level



def show_attenuation(trange, test='FEATTNTEST'):  
    dBlist, x = show_dB_ratio(trange, test)
    antennaattn = []
    for a in range(8):
        polarizationattn = []
        for p in range(2):
            cycleattn = []
            for c in range(dBlist.shape[3]):           
		attenuation = []
                for attenuator in range(5):
                    dBs = dBlist[a, p, :, c, :(attenuator+1)]
		    attns = np.sum(dBs[:, :], 1)
		    attenuation.append(attns)
		cycleattn.append(attenuation)		
            polarizationattn.append(cycleattn)
        antennaattn.append(polarizationattn)
    allattns = np.array(antennaattn)
    return allattns
   
def get_all_attns(trange, test='FEATTNTEST', from_corrected_attns=False):
    # This returns an array of all 32 possible attenuation
    #  states. All attenuationsare some combination of attenuations 
    #  0, 1, 2, 4, 8, and 16.
    #  It is called in attn_noises and get_reverseattn.    
    if from_corrected_attns == False:
        attns = show_attenuation(trange, test)
        zeros_ = np.zeros((attns.shape[0], attns.shape[1], attns.shape[2], attns.shape[4]))
        a = attns[:, :, :, 0, :]
        b = attns[:, :, :, 1, :]
        c = attns[:, :, :, 2, :]
        d = attns[:, :, :, 3, :]
        e = attns[:, :, :, 4, :]
    else:
        if from_corrected_attns == True:
            attens, bla, delta = get_attn_states(trange, test)
            attns = attens[:, :, :, :] * (1 + delta[:, :, :, :])
            zeros_ = np.zeros((attns.shape[0], attns.shape[1], attns.shape[2]))
            a = attns[:, :, :, 1]
            b = attns[:, :, :, 2]
            c = attns[:, :, :, 3]
            d = attns[:, :, :, 4]
            e = attns[:, :, :, 5]           
    f = np.array(a + b)
    g = np.array(a + c)
    h = np.array(b + c)
    i = np.array(f + c)
    j = np.array(a + d)
    k = np.array(b + d)
    l = np.array(f + d)
    m = np.array(c + d)
    n = np.array(g + d)
    o = np.array(h + d)
    p = np.array(f + m)
    aa = np.array(a + e)
    bb = np.array(b + e)
    cc = np.array(c + e)
    dd = np.array(d + e)
    ff = np.array(f+e)
    gg = np.array(g + e)
    hh = np.array(h + e)
    ii = np.array(i + e)
    jj = np.array(j + e)
    kk = np.array(k + e)
    ll = np.array(l + e)
    mm = np.array(m + e)
    nn = np.array(n + e)
    oo = np.array(o + e)
    pp = np.array(p + e)
    fullattnlist = [zeros_, a, b, f, c, g, h, i, d, j, k, l, m, n, o, p, e, aa, bb, ff, cc, gg, hh, ii, dd, jj, kk, ll, mm, nn, oo, pp]
    fullattnlist = np.array(fullattnlist)
    return fullattnlist  

def get_all_avg_attns(trange, test='FEATTNTEST', from_corrected_attns=False):
        all_attns_avg = []
        all_attns = get_all_attns(trange, test, from_corrected_attns)
        if from_corrected_attns == False:
	    for attn in range(all_attns.shape[0]):
	        ant_avgs = []
	        for ant in range(all_attns.shape[1]):
		    pol_avgs = []
		    for pol in range(all_attns.shape[2]):
		        freq_avgs = []
		        for freq in range(all_attns.shape[4]):
		            freq_avgs.append(np.average(all_attns[attn, ant, pol, :, freq]))
		        pol_avgs.append(freq_avgs)
		    ant_avgs.append(pol_avgs)
	        all_attns_avg.append(ant_avgs)
	    all_attns_avg = np.array(all_attns_avg)
        else:
            all_attns_avg = all_attns
                   
        return all_attns_avg    
    
def show_attn_comparison(trange):
    #This routine finds the difference between 
    #  actual attenuation levels with the noise
    #  on and off. It does it for both off sun 
    #  and on sun.
    attns = get_all_attns(trange)
    antenna = []
    for a in range(8):
        polarization = []
        for p in range(2):
            attenuation = []
            for attenuator in range(31):
                NDdifference = []
                NDdifference_offsun = attns[attenuator, a, p, 0, :] - attns[attenuator, a, p, 1, :]
                NDdifference_onsun = attns[attenuator, a, p, 2, :] - attns[attenuator, a, p, 3, :]
	        NDdifference.append(NDdifference_offsun)
                NDdifference.append(NDdifference_onsun)
                attenuation.append(NDdifference)
	    polarization.append(attenuation)
        antenna.append(polarization)
    attn_comparison = np.array(antenna)
    return attn_comparison  

def show_attn_plt(trange, test='FEATTNTEST'):
    #This creates useful plots showing how far off the actual
    #  attenuations are from the theoretical attenuations,
    #  for attenuations 1, 2, 4, 8, and 16. It shows this 
    #  difference for all frequencies. 
    attns = get_all_attns(trange, test)
    idx = np.array([0, 1, 2, 4, 8, 16])
    if test == 'FEATTNTEST':
        figure_names = ['Second Attenuator: ND On, Off Sun \n Attenuation(dB) vs Attenuator' , 'Second Attenuator: ND Off, Off Sun \n Attenuation(dB) vs Attenuator' , 'Second Attenuator: ND On, On Sun \n Attenuation(dB) vs Attenuator' , 'Second Attenuator: ND Off, On Sun \n Attenuation(dB) vs Attenuator']
    else:
        if test == 'FEATTNTEST2':
            figure_names = ['First Attenuator \n Attenuation(dB) vs Attenuator' ]
    antennas = np.arange(1, 9)
    polarizations = [' x', ' y']
    for j in range(attns.shape[3]):
        f, ax = plt.subplots(4, 4, sharex=True, sharey=True)
        for i in range(16):
            ax[i/4, i % 4].plot(idx,attns[idx, i/2, i % 2, j, :],'.') 
            ax[i/4, i % 4].plot(idx, idx)
            ax[i/4, i % 4].text(1, 17.5, 'Antenna ' + str(antennas[i/2]) + polarizations[i % 2])
            ax[i/4, i % 4].set_ylim([0, 22])
        plt.suptitle(figure_names[j])
        plt.show()

def get_state_idx(trange, cycles=4, attenuator=1):
    #This program creates an array of shape (6, 4) which contains the 
    #  times in which each attenuator is in each state. 6 attenuators,
    #  4 cycles. 
    firstorsecond = ['First', 'Second']
    s = sp.Spectrogram(trange)
    s.docal = False
    s.dosub = False
    s.domedian = False
    cursor = dbutil.get_cursor()
    res, msg = dbutil.do_query(cursor,'select Timestamp,Ante_Fron_FEM_HPol_Atte_First,Ante_Fron_FEM_HPol_Atte_Second from fV54_vD15 where Timestamp between '+str(trange[0].lv)+' and '+str(trange[1].lv))
    cursor.close()
    if msg == 'Success':
        antlist = []
        for i in [0, 1, 2, 4, 8, 16]:
            statelist = []
            for j in range(15):
                state, = np.where(np.logical_and(res['Ante_Fron_FEM_HPol_Atte_' + firstorsecond[attenuator]][j::15].astype('int') == i,res['Ante_Fron_FEM_HPol_Atte_' + firstorsecond[attenuator-1]][j::15].astype('int') != 0))
                statelist.append(state)
            statelist = np.array(statelist)
            antlist.append(statelist)
        states = np.array(antlist)
        states = np.rollaxis(states, 1)
        for i in range(15):
            for j in range(6):
                states[i, j] = res['Timestamp'][i::15][states[i, j]]
    else:
        print 'failure'
        return None
    time_array = (s.time.lv+0.001).astype('int')
    time_list = list(time_array)
    attns = ['0', '1', '2', '4', '8', '16']
    common_list = []
    for j in range(6):
        # Antenna 1 is used as the reference antenna here. 
        #  Earlier versions had only indices which were shared 
        #  for attenuations AND antennas, but because of a
        #  small timing error that occur between antennas
        #  during the scan itself, this older version would
        #  fail sometimes. 
        i1, i2 = util.common_val_idx(time_array,states[0,j])
        if i1.shape == i2.shape:
            common_ant_list = i2
        else: 
            print 'There is a problem with antenna '+str(i)+' at attenuation '+attns[j]
        common_list.append(common_ant_list)
    
    final_indices = []
    final_indices1 = []
    for i in range(6):
        index_list = []
        for indxs in common_list[i]:
            try:
                index_list.append(time_list.index(states[0,i][indxs]))
            except:
                pass
        final_indices1.append(index_list)
    for i in range(6):
        indices_array = np.array(final_indices1[i])
        final_indices.append(indices_array)
    final_indices = np.array(final_indices) 

    rolled_indices = []
    for i in range(6):
        rolled = np.roll(final_indices[i], -1)
        rolled_indices.append(rolled)
    subtracted_list = []
    for j in range(6):
        subtracted_list.append(rolled_indices[j] - final_indices[j])
    break_lists = []
    for k in range(6):
        break_list = []
        for indx in range(subtracted_list[k].shape[0]):
            if np.absolute(subtracted_list[k][indx]) <= 2:
                break_list.append(indx)
            else: 
                break_list.append(-1)
        break_lists.append(break_list)
    for i in range(6):
        for indx in range(int(len(break_lists[i]))-1):
            try:
                if break_lists[i][indx] == break_lists[i][indx-1]:
                    break_lists[i].pop(indx)
            except:
                pass
    break_list = []
    for j in range(6):
        breaklist = np.array(break_lists[j])
        break_list.append(breaklist)
    break_spots = []
    for i in range(6):
        try:
            break_spot = []
            for indx in range(len(break_list[i])):
                if break_list[i][indx] == -1:
                    break_spot.append(indx)
            break_spots.append(break_spot)
        except:
            pass
    split_lists = []
    for k in range(6):
        steps_list = [break_list[k][0:break_spots[k][0]]]
        for j in range(cycles-1):
            try:
                steps_list.append(break_list[k][1 + break_spots[k][j]:break_spots[k][j+1]])
            except:
                pass            
        split_lists.append(steps_list)
    split_lists = np.array(split_lists)  
    final_grouped_indices = []
    for i in range(6):
        grouped_indices = []
        for j in range(cycles):
            try:
                indices_ = []
                for indxs in split_lists[i][j]:
                    indices_.append(rolled_indices[i][indxs])
                grouped_indices.append(indices_)
            except:
                pass
        final_grouped_indices.append(grouped_indices)
    final_grouped_indices = np.array(final_grouped_indices)
    for i in range(6):
        for j in range(cycles):
            try:
                for k in range(1,int(len(final_grouped_indices[i][j]))-1):
                    try:
                        for m in range(len(final_ped_indices[i][j])):
                            if (final_grouped_indices[i][j][k-1] + 3) <= final_grouped_indices[i][j][k]:
                                final_grouped_indices[i][j].pop(k-1) 
                            if (final_grouped_indices[i][j][k+1]-3) >= final_grouped_indices[i][j][k]:
                                final_grouped_indices[i][j].pop(k+1)       
                    except:
                        pass
            except:
                pass
    return final_grouped_indices, res

def find_gaincal(t=None, scan_length=6, findwhat='FEATTNTEST'):
    # This will find the project on the day of "t" or, going backwards, the nearest
    #  day to "t" 
    #  This routine returns a timerange that may be used as the timerange in any 
    #  of the above programs. It will return the appropriate timerange for a 
    #  FEATTNTEST scan on the date for which t is given to this function, or a
    #  list of appropriate time range. If the length of the scan is not 6 minute,
    #  please set scan_length to the appropriate length. 
    ''' Makes an SQL query to find the FEATTNTEST or GAINCALTESTscans for the date given in
        the Time() object t.  A list of timestamps is returned, along with 
        the timestamp of the object provided.
    '''
    loop_ = 1
    if t is None: 
        # Get today's date
        t = util.Time.now()
    timestamp = int(t.lv)
    while loop_ == 1:
        stimestamp = timestamp - (timestamp % 86400)  # Start of day
        etimestamp = stimestamp + 86399               # End of day
        # Open handle to SQL database
        cursor = dbutil.get_cursor()
        # Try to find a scan header with project SOLPNTCAL (only works after 2014 Nov. 30)
        verstr = dbutil.find_table_version(cursor,timestamp,True)
        if verstr is None:
            print 'No scan_header table found for given time.'
            return [], timestamp
        # First retrieve the Project from all scan headers for the day
        cursor.execute('select timestamp,Project from hV'+verstr+'_vD1 where timestamp between '+str(stimestamp)+' and '+str(etimestamp)+' order by timestamp')
	data = np.transpose(np.array(cursor.fetchall()))
        names = stateframedef.numpy.array(cursor.description)[:,0]
        cursor.close()
        if len(data) == 0:
            # No FEATTNTEST found, so return empty list (and timestamp)
            return [], timestamp
        else:
            projdict = dict(zip(names,data))
            projdict['timestamp'] = projdict['timestamp'].astype('float')  # Convert timestamps from string to float
        good = np.where(projdict['Project'] == findwhat)[0]
        if len(good) != 0:
            if len(good) == 1:
                loop_ = 0
                tgc = [projdict['timestamp'][good]], timestamp
                start_= Time(tgc[0][0], format = 'lv').iso
                end_ = Time(tgc[0][0]+60*scan_length, format = 'lv').iso
                trange = Time([start_[0] , end_[0]])
                return trange
            else:
                loop_ = 0
                tgc = [projdict['timestamp'][good]], timestamp
                start_trange = Time(tgc[0][0], format = 'lv').iso
                end_trange = Time(tgc[0][0]+60*scan_length, format = 'lv').iso
                tranges = []
                for i in range(start_trange.shape[0]):
                    trange = Time([start_trange[i], end_trange[i]])
                    tranges.append(trange)
                return tranges
        else:
            timestamp = timestamp - 60*60*24
            loop_ = 1
        
def get_attn_states(trange, test='FEATTNTEST'):
    allattns = get_all_attns(trange, test, from_corrected_attns=False)
    main_attns = []
    for attn in [0, 1, 2, 4, 8, 16]:
        main_attns.append(allattns[attn, :, :, :, :])
    main_attns = np.array(main_attns)
    targets = [0, 1, 2, 4, 8, 16]
    polars = ['x', 'y']
    mainattns = []
    for a in range(allattns.shape[1]):
        pattns = []
        for p in range(allattns.shape[2]):
            fattns = []
            for f in range(allattns.shape[4]):
                attns_ = []
                attns_avg = []
                attns_.append(np.average(main_attns[0, a, p, :, f]))              
                for attns in range(1,6):
                    attns_avg = []
                    steps_ = main_attns[attns, a, p, :, f]
                    stepdif = abs(main_attns[attns, a, p, :, f] - main_attns[attns-1, a, p, :, f])
                    attns_avg.append(attns_[attns-1]+stepdif)
                    attns_.append(np.average(attns_avg))
                fattns.append(attns_)
            pattns.append(fattns)
        mainattns.append(pattns)  
    mainattns = np.array(mainattns) 

    attenuations_ = ['0', '1', '2', '4', '8', '16']
    for a in range(mainattns.shape[0]):
        for p in range(mainattns.shape[1]):
            for f in range(mainattns.shape[2]):
                for attns in range(mainattns.shape[3]):
                    if (targets[attns] - 10) <= mainattns[a, p, f, attns] <= (targets[attns] + 10):
                        pass
                    else:
                        mainattns[a, p, f, attns] = targets[attns]                             
    avg_per_attn = []
    for a in range(mainattns.shape[0]):
        avg_for_pols = []
        for p in range(mainattns.shape[1]):
            avg_for_attn = []
            for attn in range(mainattns.shape[3]):
                avg_for_attn.append(np.average(mainattns[a, p, :, attn]))    
            avg_for_pols.append(avg_for_attn)
        avg_per_attn.append(avg_for_pols)
    avg_per_attn = np.array(avg_per_attn)
    diff_per_freq = []
    for a in range(mainattns.shape[0]):
        per_pol = []
        for p in range(mainattns.shape[1]):
            per_attn = []
            for attn in range(mainattns.shape[3]):
                per_freq = []
                for freq in range(mainattns.shape[2]):
                    per_freq.append(avg_per_attn[a, p, attn] - mainattns[a, p, freq, attn])
                per_attn.append(per_freq)
            per_pol.append(per_attn)
        diff_per_freq.append(per_pol)
    diff_per_freq = np.array(diff_per_freq)
    freq_percent = []
    for a in range(mainattns.shape[0]):
        per_pol = []
        for p in range(mainattns.shape[1]):
            per_freq = []
            for freq in range(mainattns.shape[2]):
                per_attn = []
                for attn in range(mainattns.shape[3]):
                    per_attn.append(diff_per_freq[a, p, attn, freq]/avg_per_attn[a, p, attn])
                per_freq.append(per_attn)
            per_pol.append(per_freq)
        freq_percent.append(per_pol)
    freq_percent = np.array(freq_percent) 
    delta_percent = np.nan_to_num(freq_percent)       
    return mainattns, avg_per_attn, delta_percent
    
def make_attntable(trange, type='attenuation', test='FEATTNTEST'):
    #The type can be 'attenuation' or 'delta'. delta is the percent difference at each frequency between     
    #the average measured attenuation and the measured attenuation at that frequency. 
    polars = ['x', 'y']
    mainattns, avg_per_attn, delta_percent = get_attn_states(trange, test, tsys, s)
    if type == 'attenuation':
        rowlabels = []
        for i in range(16):
            label_ = 'Ant' + str((i/2)+1) + ',' + polars[i % 2]
            rowlabels.append(label_)
        cols = []
        for j in range(5):
            rows = []
            for i in range(16):
                row = tbl.Column(data=avg_per_attn[i/2, i % 2, j], name=(str(j)+'dB'))
                rows.append(row)
            cols.append(rows)
        cols.insert(0, rowlabels)
        cols = np.array(np.transpose(cols))
        t = tbl.Table(cols, names=['Ant, Pol' , '1 dB', '2 dB', '4 dB', '8 dB', '16 dB'])
        print t
    else:
        if type == 'delta':
            columnlabels = ['freq']
            for i in range(16):
                label_ = str('Ant ' + str((i/2)+1) + ' ' + polars[i % 2])
                columnlabels.append(label_)
            cols = []    
            for i in range(163):
                rows = [i+1]
                for j in range(16):
                    row = delta_percent[j/2, j % 2, i]
                    rows.append(row)
                cols.append(rows)
            #cols.insert(0, rowlabels)
            cols = np.array(cols) 
            t = tbl.Table(cols, names = columnlabels)  
            print t
        else:
            print 'Please select either attenuation or delta for the table'
            t = None
    return t

def make_allattndicts(trange, per_ant=False):
    # This routine makes dictionaries of tables to help sort the attenuation 
    #  information from a gaincal scan. It can sort it in two ways:
    #  A dictionary of frequencies, with the attenuation for each
    #  antenna and polarization in each item, or a dictionary of
    #  antennas and polarization, with the attenuation information
    #  for each frequency. To easily view a large table here, 
    #  just use dictionary['key'].more() this will bring up an 
    #  interactive guide from astropy to help go through the data. 
    if per_ant == False:
        polars = ['x', 'y']
        mainattns, avg_per_attn, delta_percent = get_attn_states(trange)
        all_tables = []
        for f in range(163):
            rowlabels = []
            for i in range(16):
                label_ = 'Ant' + str((i/2)+1) + ',' + polars[i % 2]
                rowlabels.append(label_)
            cols = []
            for j in range(5):
                rows = []
                for i in range(16):
                    row = tbl.Column(data=mainattns[i/2, i % 2, f, j], name=(str(j)+'dB'))
                    rows.append(row)
                cols.append(rows)
            cols.insert(0, rowlabels)
            cols = np.array(np.transpose(cols))
            t = tbl.Table(cols, names=['Ant, Pol' , '1 dB', '2 dB', '4 dB', '8 dB', '16 dB'])
            all_tables.append(t)
        freqs = np.arange(1,164).tolist()
        for freq in range(len(freqs)):
            freqs[freq] = 'Frequency channel ' + str(freqs[freq])
        freq_attndict = {}
        for i in range(163):
            freq_attndict[freqs[i]] = all_tables[i]
        return freq_attndict
    if per_ant == True:
        polars = ['x', 'y']
        mainattns, avg_per_attn, delta_percent = get_attn_states(trange)
        all_tables = []
        for i in range(16):
            cols = [np.arange(1, 164).tolist()]
            for j in range(5):
                rows = []
                for f in range(163):
                    row = tbl.Column(data=mainattns[i/2, i % 2, f, j])
                    rows.append(row)
                cols.append(rows)   
            t = tbl.Table(cols, names=['Freqs', '1 dB', '2 dB', '4 dB', '8 dB', '16 dB']) 
            all_tables.append(t)
        ant_attndict = {}
        columnlabels = []
        for i in range(16):
            label_ = 'Antenna ' + str((i/2)+1) + ',' + polars[i % 2] + ' polarization'
            columnlabels.append(label_)
        for j in range(16):
            ant_attndict[columnlabels[j]] = all_tables[j]        
        return ant_attndict

def attn_noises(trange_gaincal):
        # This function is used with "get_reverseattn." It find and returns background noise, and returns the "res"
        #  from dbutil.do_query and returns the "tsys" from spectrogram_fit and get_data().
        #  The first parameter it takes should be a GAINCALTEST trange, which may be located with find_gaincal, 
        #  The second parameter it takes should be a FEATTNTEST trange, which may be located with find_gaincal as
        #  well. 
        #
        #  PLEASE NOTE: ANY trange with data may be used as "trange_gaincal", use a trange from a GAINCALTEST and the other 
        #  file as "trange_other" if you want the noise to be calibrated from the GAINCALTEST file, which will most likely
        #  be more recent than the FEATTNTEST file it would otherwise take the noise from. 
        s = sp.Spectrogram(trange_gaincal)
	s.docal = False
	s.dosub = False
	s.domedian = False        
	tsys, std = s.get_data()
        trange_feattncal = find_gaincal()
        if type(trange_feattncal) == list:
            trange_feattncal = trange_feattncal[-1]
        else:
            pass
        trange_feattncal2 = find_gaincal(t = Time('2015-07-21 00:00'), scan_length=5, findwhat='FEATTNTEST2')
        if type(trange_feattncal2) == list:
            trange_feattncal2 = trange_feattncal2[-1]
        else:
            pass
        ratios, calfilenoise = show_dB_ratio(trange_feattncal)
        ratios1, calfilenoise1 = show_dB_ratio(trange_feattncal2, test='FEATTNTEST2')

	cursor = dbutil.get_cursor()
	res, msg = dbutil.do_query(cursor,'select Timestamp,Ante_Fron_FEM_HPol_Atte_First,Ante_Fron_FEM_HPol_Atte_Second from fV54_vD15 where Timestamp between '+str(trange_gaincal[0].lv)+' and '+str(trange_gaincal[1].lv))
	cursor.close()

        idx1, idx2 = util.common_val_idx(res['Timestamp'][0::15].astype('int'), (s.time.lv+0.5).astype('int'))       
        idx3, idx4 = util.common_val_idx(res['Timestamp'].astype('int'), (s.time.lv+0.5).astype('int'))
        marker = -1
        while idx1[-1] > idx2[-1]:
            idx1 = np.delete(idx1, -1)
            marker += 1
        tsys = tsys[:, :, :, idx1]

        calfilenoise_ = []
        for ant in range(calfilenoise.shape[0]):
            calfilenoisepol = []
            for pol in range(calfilenoise.shape[1]):
                calfilenoisefreq = []
                for freq in range(calfilenoise.shape[2]):
                    calfilenoisefreq.append(np.average(calfilenoise[ant, pol, freq, :]))
                calfilenoisepol.append(calfilenoisefreq)
            calfilenoise_.append(calfilenoisepol)
        calfilenoise = np.array(calfilenoise_)

	noise_level = []
	for ant in range(tsys.shape[0]):
	    pol_noise = []
	    for pol in range(tsys.shape[1]):
		freq_noise = []
		state, = np.where(np.logical_and(res['Ante_Fron_FEM_HPol_Atte_Second'][ant::15].astype('int') == 31,res['Ante_Fron_FEM_HPol_Atte_First'][ant::15].astype('int') == 31))
		for freq in range(tsys.shape[2]):
		    avg_noise = []
		    for index in state:
                        try:
		            if np.logical_and(tsys[ant, pol, freq, index] <= 0.005, index < tsys.shape[3]):
		                avg_noise.append(tsys[ant, pol, freq, index])
                        except:
                            pass
		    freq_noise.append(np.average(avg_noise))
		pol_noise.append(freq_noise)
	    noise_level.append(pol_noise)
	noise_level = np.array(noise_level)

        for ant in range(tsys.shape[0]):
	    for pol in range(tsys.shape[1]):
		for freq in range(tsys.shape[2]):                      
		    if np.isnan(noise_level[ant, pol, freq]) == False:
                        pass
                    else:
                        if np.isnan(noise_level[ant, pol, freq]) == True:                           
                            try:
                                noise_level[ant, pol, freq] = calfilenoise[ant, pol, freq]
                            except:
                                pass

        return tsys, res, noise_level, idx1, idx3, marker, trange_feattncal, trange_feattncal2
    

def get_reverseattn(trange_gaincal, trange_other=None, first_attn_base=5, second_attn_base=3, corrected_attns=False):
        # This function finds and return the tsys, the noise corrected tsys, and the noise and attenuation
        #  tsys. The first parameter it takes should be a GAINCALTEST trange, which may be located with find_gaincal, 
        #  The second parameter it takes should be a FEATTNTEST trange, which may be located with find_gaincal as
        #  well. It may or may not take a third parameter. Any form of file that was recorded by the system from 
        #  from the antennas may be inserted here, whether a flare or just quiet sun data.
        #
        #  PLEASE NOTE: ANY trange with data may be used as "trange_gaincal", use a trange from a GAINCALTEST and the other 
        #  file as "trange_other" if you want the noise to be calibrated from the GAINCALTEST file, which will most likely
        #  be more recent than the FEATTNTEST file it would otherwise take the noise from. 
        tsys1, res1, noise_level, idx1a, idx3a, marker1, trange_feattncal, trange_feattncal2  = attn_noises(trange_gaincal)
        if corrected_attns == True:
            all_attns_avg = get_all_attns(trange_feattncal, test='FEATTNTEST', from_corrected_attns=True)
            all_attns_avg1 = get_all_attns(trange_feattncal2, test='FEATTNTEST2', from_corrected_attns=True)
        else:
            if corrected_attns == False:
                all_attns_avg = get_all_avg_attns(trange_feattncal, test='FEATTNTEST')
                all_attns_avg1 = get_all_avg_attns(trange_feattncal2, test='FEATTNTEST2')
        
        if trange_other == None:
            tsys = tsys1
            res = res1
            idx1 = idx1a
            idx3 = idx3a
            marker = marker1
        else:
            s = sp.Spectrogram(trange_other)
      	    s.docal = False
	    s.dosub = False
	    s.domedian = False        
	    tsys, std = s.get_data()
            cursor = dbutil.get_cursor()
	    res, msg = dbutil.do_query(cursor,'select Timestamp,Ante_Fron_FEM_HPol_Atte_First,Ante_Fron_FEM_HPol_Atte_Second from fV54_vD15 where Timestamp between '+str(trange_other[0].lv)+' and '+str(trange_other[1].lv))
	    cursor.close()
            idx1, idx2 = util.common_val_idx(res['Timestamp'][0::15].astype('int'), (s.time.lv+0.5).astype('int'))
            marker = -1
            while idx1[-1] > idx2[-1]:
                idx1 = np.delete(idx1, -1)
                marker += 1
            tsys = tsys[:, :, :, idx1]
            idx3, idx4 = util.common_val_idx(res['Timestamp'].astype('int'), (s.time.lv+0.5).astype('int'))
            res['Timestamp'] = res['Timestamp'][idx3]

        if noise_level.shape[2] < tsys.shape[2]:
            freqs_for_range = noise_level.shape[2]
        else:
            freqs_for_range = tsys.shape[2]
                            
	tsys_noise_corrected = []
	for ant in range(tsys.shape[0]):
	    pol_corrected = []
	    for pol in range(tsys.shape[1]):
		freq_corrected = []
		for freq in range(freqs_for_range):
		    index_corrected = []
		    for index in range(tsys.shape[3]):
		        index_corrected.append(tsys[ant, pol, freq, index] - noise_level[ant, pol, freq])
		    freq_corrected.append(index_corrected)
		pol_corrected.append(freq_corrected)
	    tsys_noise_corrected.append(pol_corrected)
	tsys_noise_corrected= np.array(tsys_noise_corrected)
        
        freqloopslist = [tsys_noise_corrected.shape[2], all_attns_avg.shape[3], all_attns_avg1.shape[3]]
        freqloops = min(freqloopslist)

        if tsys.shape[3] < len(res['Ante_Fron_FEM_HPol_Atte_Second'][0::15]):
            indexloops = tsys.shape[3]
        else:
            if tsys.shape[3] >= len(res['Ante_Fron_FEM_HPol_Atte_Second'][0::15]):
                indexloops = len(res['Ante_Fron_FEM_HPol_Atte_Second'][0::15])-1  
        if tsys.shape[3] < len(res['Ante_Fron_FEM_HPol_Atte_First'][0::15]):
            indexloops1 = tsys.shape[3]
        else:
            if tsys.shape[3] >= len(res['Ante_Fron_FEM_HPol_Atte_First'][0::15]):
                indexloops1 = len(res['Ante_Fron_FEM_HPol_Atte_First'][0::15])-1   
        idxstart = marker + (15-8)
        xory = ['x' , 'y']
	ant_postcorrected = []
	for ant in range(tsys.shape[0]):
	    pol_postcorrected = []
	    for pol in range(tsys.shape[1]):
		freq_postcorrected = []
		for freq in range(freqloops):
                     
                    indices_postcorrected = []
                    for indx in range(indexloops): 
                        testlevel = res['Ante_Fron_FEM_HPol_Atte_Second'][ant::15][indx+idxstart] 
                        
                        if 0 <= testlevel <= 31:
                            pass
                        else:
                            print 'Problem with the attenuation of antenna ' + str(ant) + xory[pol] + ' at frequency channel ' + str(freq) + ' and time index '  + str(indx) + '. The attenuation is showing: ' + str(testlevel)       
                            testlevel = 0       
		        indices_postcorrected.append(10**((all_attns_avg[testlevel, ant, pol, freq]-all_attns_avg[second_attn_base, ant, pol, freq])/10)*tsys_noise_corrected[ant, pol, freq, indx])
                    indices_postcorrected1 = []
                    for indx in range(indexloops1): 
                        testlevel = res['Ante_Fron_FEM_HPol_Atte_First'][ant::15][indx+idxstart]                        
                        if 0 <= testlevel <= 31:
                            pass
                        else:
                            print 'Problem with the attenuation of antenna ' + str(ant) + xory[pol] + ' at frequency channel ' + str(freq) + ' and time index '  + str(indx) + '. The attenuation is showing: ' + str(testlevel)       
                            testlevel = 0      
		        indices_postcorrected1.append(10**((all_attns_avg1[testlevel, ant, pol, freq]-all_attns_avg1[first_attn_base, ant, pol, freq])/10)*indices_postcorrected[indx])                            
		    freq_postcorrected.append(indices_postcorrected1)
		pol_postcorrected.append(freq_postcorrected)
	    ant_postcorrected.append(pol_postcorrected)
	tsys_attn_noise_corrected = np.array(ant_postcorrected)
        
        return tsys_attn_noise_corrected, tsys_noise_corrected, tsys

def show_reverseattn(trange_gaincal, trange_other=None, first_attn_base=5, second_attn_base=3, corrected_attns=False):
        # This routine returns the same things as get_reverseattn, but it also graphs the three different "tsys"
        #  variables. It interactively asks the user which antenna, which polarization, and which frequency they 
        #  wish to view. Once it creates a graph, it will ask the user whether they wish to see more graphs or not.
        #  The first parameter it takes should be a GAINCALTEST trange, which may be located with find_gaincal, 
        #  The second parameter it takes should be a FEATTNTEST trange, which may be located with find_gaincal as
        #  well. It may or may not take a third parameter. Any form of file that was recorded by the system from 
        #  from the antennas may be inserted here, whether a flare or just quiet sun data. 
        #
        #  PLEASE NOTE: ANY trange with data may be used as "trange_gaincal", use a trange from a GAINCALTEST and the other 
        #  file as "trange_other" if you want the noise to be calibrated from the GAINCALTEST file, which will most likely
        #  be more recent than the FEATTNTEST file it would otherwise take the noise from. 
        tsys_attn_noise_corrected, tsys_noise_corrected, tsys = get_reverseattn(trange_gaincal, trange_other, first_attn_base, second_attn_base, corrected_attns)
        plotting_ = True
        while plotting_ == True:
            print ' '
            antenna = input('Which antenna would you like to see? You can choose from 1 to 8 : ')
            print ' '
            polarization = input('Which polarization would you like to see? Enter 0 for x, and 1 for y: ') 
            print ' '
            channel = input('Which channel would you like to see? You can choose from channel 0 to ' +str(tsys_attn_noise_corrected.shape[2]) + ' : ')

            antenna = int(antenna)-1
            polarization = int(polarization)
            channel = int(channel)

            plt.figure()
            tsys_ = plt.plot(tsys[antenna, polarization, channel, :], label='tsys')
            tsys_noise_corrected_ = plt.plot(tsys_noise_corrected[antenna, polarization, channel, :], label='noise corrected tsys')
            tsys_attn_noise_corrected_ = plt.plot(tsys_attn_noise_corrected[antenna, polarization, channel, :], label='noise and attn corrected tsys')
            plt.yscale('log')
            fontP = FontProperties()
            fontP.set_size('small')
            plt.legend(prop = fontP)
            plt.show()

            print ' '
            yesorno = input('Would you like to see another plot? Please answer 0 for yes or 1 for no. ')
            if yesorno == 0:
                pass
            else:
                print ' '
                print 'See you next time!'
                plotting_ = False
                return tsys_attn_noise_corrected, tsys_noise_corrected, tsys


















               
                

        

      
            


    





