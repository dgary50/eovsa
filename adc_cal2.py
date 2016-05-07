#
# Routines to set up system for ADC calibration/DCM attenuation setting
#
#   2016-Feb-20  DG
#     First written.
#   2016-Mar-19  DG
#     Converted to adc_cal2.py (for now), and added set_fem_attn() routine.
#   2016-May-05  DG
#     Added a new routine DCM_master_attn_cal() to do the quickest possible
#     calibration--takes only 5 minutes.  It returns DCM_lines, ready to be
#     inserted into the SQL table by
#        import cal_header 
#        cal_header.dcm_master_table2sql(DCM_lines)
#

import time
import numpy as np
import roach as r
import stateframe as stf

def acc_tune(band,acc):
    if type(band) is int:
        fsqfile = 'BAND'+str(band)+'.FSQ'
    elif type(band) is str:
        if band.lower() == 'solar.fsq' or band.lower() == 'pcal.fsq':
            fsqfile = band.lower()
    else:
        print 'Error: Unknown band',band
        return
    cmds = ['FSEQ-OFF','FSEQ-INIT','WAIT','FSEQ-FILE '+fsqfile.lower(), 'FSEQ-ON']
    send_cmds(cmds,acc)

def send_cmds(cmds,acc):
    ''' Sends a series of commands to ACC.  The sequence of commands
        is not checked for validity!
        
        cmds   a list of strings, each of which must be a valid command
    '''
    import socket, stateframe

    for cmd in cmds:
        #print 'Command:',cmd
        s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        try:
            s.connect((acc['host'],acc['scdport']))
            s.send(cmd)
            time.sleep(0.01)
            s.close()
        except:
            print 'Error: Could not send command',cmd,' to ACC.'
    return

def ant_str2list(ant_str):
    ant_list = []
    try:
        grps = ant_str.split()
        for grp in grps:
            antrange = grp[3:].split('-')
            if len(antrange) == 1:
                if antrange != '':
                    ant_list.append(int(antrange[0])-1)
            elif len(antrange) == 2:
                ant_list += range(int(antrange[0])-1,int(antrange[1]))
    except:
        print 'Error: cannot interpret ant_str',ant_str
        return None
    return np.array(ant_list)
                    
def set_fem_attn(level=3,ant_str='ant1-15'):
    ''' Read current power and attenuation values in FEMs, and set the attenuation
        such that the power level will be as close to "level" as possible.
    '''
    accini = stf.rd_ACCfile()
    acc = {'host': accini['host'], 'scdport':accini['scdport']}
    
    ant_list = ant_str2list(ant_str)
    if ant_list is None:
        return 'Bad antenna list'
    accini = stf.rd_ACCfile()
    hatn1 = np.zeros((10,15),dtype='int')
    vatn1 = np.zeros((10,15),dtype='int')
    hatn2 = np.zeros((10,15),dtype='int')
    vatn2 = np.zeros((10,15),dtype='int')
    hpwr = np.zeros((10,15),dtype='float')
    vpwr = np.zeros((10,15),dtype='float')
    # Read 10 instances of attenuation and power, and take the median to avoid
    # glitches
    for i in range(10):
        # Read the frontend attenuations and powers for each antenna
        data, msg = stf.get_stateframe(accini)
        for iant in range(15):
            fem = accini['sf']['Antenna'][iant]['Frontend']['FEM']
            hatn1[i,iant] = stf.extract(data,fem['HPol']['Attenuation']['First'])
            vatn1[i,iant] = stf.extract(data,fem['VPol']['Attenuation']['First'])
            hatn2[i,iant] = stf.extract(data,fem['HPol']['Attenuation']['Second'])
            vatn2[i,iant] = stf.extract(data,fem['VPol']['Attenuation']['Second'])
            hpwr[i,iant] = stf.extract(data,fem['HPol']['Power'])
            vpwr[i,iant] = stf.extract(data,fem['VPol']['Power'])
        time.sleep(1)
    hatn1 = np.median(hatn1,0).astype('int')
    vatn1 = np.median(vatn1,0).astype('int')
    hatn2 = np.median(hatn2,0).astype('int')
    vatn2 = np.median(vatn2,0).astype('int')
    hpwr = np.median(hpwr,0)
    vpwr = np.median(vpwr,0)
    hatn2 = np.clip(hatn2 - (level - hpwr).astype('int'),0,31)
    vatn2 = np.clip(vatn2 - (level - vpwr).astype('int'),0,31)
    # Send attenuation commands to each antenna in ant_list
    for iant in ant_list:
        hatn = str(hatn1[iant])+' '+str(hatn2[iant])+' ant'+str(iant+1)
        vatn = str(vatn1[iant])+' '+str(vatn2[iant])+' ant'+str(iant+1)
        send_cmds(['HATTN '+hatn,'VATTN '+vatn],acc)
    return 'Success'
    
def set_dcm_attn(roach_list,fem_level=5,nd_state='on',adc_nom=25,ant_list='ant1-13',do_plot=False):
    ''' Set the FEM attenuation to the given value, switch ND state to the given state
        and cycle through the IF bands, and find the optimum DCM attenuation needed to 
        attain the given ADC level.
        
        Returns the corresponding DCM table
    '''
    import copy
    accini = stf.rd_ACCfile()
    acc = {'host': accini['host'], 'scdport':accini['scdport']}
    
    # Switch noise diode to requested state
    if nd_state == 'on':
        send_cmds(['ND-ON '+ant_list],acc)
    else:
        send_cmds(['ND-OFF '+ant_list],acc)
    # Set FEM power level to requested level
    if set_fem_attn(fem_level,ant_list) == 'Failure':
        return
    time.sleep(5)
    dcm_table = np.zeros((34,32), dtype='int')
    # Set DCM state to standard values
    send_cmds(['DCMAUTO-OFF '+ant_list,'DCMATTN 12 12 '+ant_list],acc)
    # Cycle through bands to get ADC levels
    for band in range(34):
        # Start with nominal DCM attenuation
        send_cmds(['DCMATTN 12 12 '+ant_list],acc)
        time.sleep(0.3)
        print 'Band:',band+1
        acc_tune(band+1,acc)
        time.sleep(3)
        # Go through ROACH list
        for i,ro in enumerate(roach_list):
            # Get ADC levels at nominal setting
            dcm_base = np.array([12,12,12,12],dtype='int')
            r.adc_levels([ro])
            # Calculate new attenuation to achieve nominal level
            ch_atn = np.clip(((20*np.log10(ro.adc_levels/adc_nom)+dcm_base + 1)/2).astype('int')*2,0,30)
            # Set new attenuation levels and check result
            send_cmds(['DCMATTN '+str(ch_atn[0])+' '+str(ch_atn[1])+' ant'+str(ro.ants[0])],acc)
            time.sleep(1)
            send_cmds(['DCMATTN '+str(ch_atn[2])+' '+str(ch_atn[3])+' ant'+str(ro.ants[1])],acc)
            time.sleep(1)
            r.adc_levels([ro])
            ch_atn2 = np.clip(((20*np.log10(ro.adc_levels/adc_nom)+ch_atn + 1)/2).astype('int')*2,0,30)
            print '  ',ro.roach_ip[:6],'Attn:',ch_atn,'Check:',ch_atn2
            dcm_table[band,np.array(((ro.ants[0]-1)*2,(ro.ants[0]-1)*2+1,(ro.ants[1]-1)*2,(ro.ants[1]-1)*2+1))] = copy.copy(ch_atn2)
    return dcm_table
    
# Insert 62 dB into FEMs, cycle through bands, get ADC levels (optionally plot results)
# Set FEM power between 2 and 3 dBm with ND-ON, cycle through bands, get ADC levels (optionally plot results)
# Use ADC level tests to "guess" best DCM attenuation settings

def adc_cal(roach_list,ant_list='ant1-15',do_plot=False):
    ''' Perform a sequence of FEM settings, using ADC levels to 
        deduce optimum DCM attenuation settings for all 34 bands.
        This can also reveal problems in FEM or DCM hardware.
        TAKES ABOUT 17 MINUTES TO COMPLETE
        
        roach_list  a set of ROACH objects created with roach.py
        ant_list    a list of antennas in the form of a string,
                      e.g. "ant1-5 ant7" on which to adjust FEMs
                      Default is all antennas, and an empty string
                      means all antennas in current subarray.
        do_plot     if True, makes a summary plot of results
        
        Returns numpy arrays :
                adc_nosig[34, nroach, 4] (no-signal ADC levels)
                adc_ndoff[34, nroach, 4] (ADC levels for ND-OFF)
                adc_ndon [34, nroach, 4] (ADC levels for ND-ON)
    '''
    accini = stf.rd_ACCfile()
    acc = {'host': accini['host'], 'scdport':accini['scdport']}

    n = len(roach_list)
    adc_nosig = np.zeros((34,n,4),dtype='float')
    adc_ndoff = np.zeros((34,n,4),dtype='float')
    adc_ndon = np.zeros((34,n,4),dtype='float')
    # Set DCM state to standard values
    send_cmds(['DCMAUTO-OFF '+ant_list,'DCMATTN 12 12 '+ant_list],acc)
    # Set FEM attenuations to maximum
    send_cmds(['FEMATTN 15 '+ant_list],acc)
    # Cycle through bands to get "zero-input" ADC levels
    for band in range(34):
        acc_tune(band+1,acc)
        time.sleep(3)
        # Go through roaches and find optimum ADC levels
        for i,ro in enumerate(roach_list):
            r.adc_levels([ro])
            
            adc_nosig[band,i] = ro.adc_levels
    # Set FEM attenuations to nominal
    send_cmds(['FEMATTN 0 '+ant_list],acc)
    # Cycle through bands to get "nd-on" ADC levels
    send_cmds(['ND-ON '+ant_list],acc)
    for band in range(34):
        acc_tune(band+1,acc)
        time.sleep(3)
        r.adc_levels(roach_list)
        for i,ro in enumerate(roach_list):
            adc_ndon[band,i] = ro.adc_levels
    # Cycle through bands to get "nd-off" ADC levels
    send_cmds(['ND-OFF '+ant_list],acc)
    for band in range(34):
        acc_tune(band+1,acc)
        time.sleep(3)
        r.adc_levels(roach_list)
        for i,ro in enumerate(roach_list):
            adc_ndoff[band,i] = ro.adc_levels
    if do_plot: 
        plot_adc_cal(roach_list, adc_nosig, adc_ndoff, adc_ndon)
    return adc_nosig, adc_ndoff, adc_ndon

def DCM_master_attn_cal(roach_list,ant_list='ant1-15'):
    ''' Perform a sequence of FEM settings, using ADC levels to 
        deduce optimum DCM attenuation settings for all 34 bands
        and return the master table DCM attenuation lines.
        TAKES ABOUT 5 MINUTES TO COMPLETE

        roach_list  a set of ROACH objects created with roach.py
        ant_list    a list of antennas in the form of a string,
                      e.g. "ant1-5 ant7" on which to adjust FEMs
                      Default is all antennas, and an empty string
                      means all antennas in current subarray.
        
        Returns DCM_lines list of strings, one for each band
    '''
    accini = stf.rd_ACCfile()
    acc = {'host': accini['host'], 'scdport':accini['scdport']}

    n = len(roach_list)
    adc_ndon = np.zeros((34,n,4),dtype='float')
    # Set DCM state to standard values
    send_cmds(['DCMAUTO-OFF '+ant_list,'DCMATTN 12 12 '+ant_list],acc)
    # Set FEM attenuations to nominal
    send_cmds(['FEMATTN 0 '+ant_list],acc)
    # Cycle through bands to get "nd-on" ADC levels
    send_cmds(['ND-ON '+ant_list],acc)
    for band in range(34):
        acc_tune(band+1,acc)
        time.sleep(3)
        r.adc_levels(roach_list)
        for i,ro in enumerate(roach_list):
            adc_ndon[band,i] = ro.adc_levels
    send_cmds(['ND-OFF '+ant_list],acc)
    DCM_lines = make_DCM_table(roach_list,adc_ndon,dcm_base=12,adc_nom=30)
    return DCM_lines
    
def plot_adc_cal(roach_list,adc_nosig,adc_ndoff,adc_ndon):
    import matplotlib.pylab as plt
    n = len(roach_list)
    chans = ['X','Y']
    f, ax = plt.subplots(n,4)
    f.set_size_inches(10,2.5*n, forward=True)
    for i in range(n):
        rstr = 'Roach'+roach_list[i].roach_ip[5:6]
        for j in range(4):
            ant = roach_list[i].ants[j / 2]
            chan = chans[j % 2]
            astr = ' Ant '+str(ant)+chan+':'
            ax[i,j].plot(adc_nosig[:,i,j],'.')
            ax[i,j].plot(adc_ndoff[:,i,j],'.')
            ax[i,j].plot(adc_ndon[:,i,j],'.')
            ax[i,j].set_ylim(0, 60)
            ax[i,j].text(5,50,rstr+astr,fontsize=10)
    plt.show()

def make_DCM_table(roach_list,adc_ndon,dcm_base=12,adc_nom=30):

    DCMlines = []
    DCMlines.append('#      Ant1  Ant2  Ant3  Ant4  Ant5  Ant6  Ant7  Ant8  Ant9 Ant10 Ant11 Ant12 Ant13 Ant14 Ant15')
    DCMlines.append('#      X  Y  X  Y  X  Y  X  Y  X  Y  X  Y  X  Y  X  Y  X  Y  X  Y  X  Y  X  Y  X  Y  X  Y  X  Y')
    DCMlines.append('#     ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- -----')
    for band in range(1,35):
        out = np.zeros(32,dtype='int') + 12  # Default to 12 dB if not present
        for i in range(len(roach_list)):
            # Calculate DCM attenuation for the 4 channels on this roach at this band
            # The target standard deviation is adc_nom (default is 30), and the base
            # attenuation at which the observations were made is dcm_base (default is 12 dB).
            # This uses the ratio of standard deviations to determine the factor in dB
            # needed to change it to the target standard deviation.  The division by two,
            # conversion to integer, and multiplication by 2 is because the attenuation steps
            # are in units of 2 dB.  The result is clipped to be between 0 and 30 dB.
            ch_atn = np.clip(((10*np.log(adc_ndon[band-1,i,:]/adc_nom)+dcm_base + 1)/2).astype('int')*2,0,30)
            # Determine the two antennas on this roach (-1 converts to 0-based index)
            ant1,ant2 = roach_list[i].ants - 1
            # Use indexes to assign the 4 channels to the right place in the array
            out[np.array((ant1*2,ant1*2+1,ant2*2,ant2*2+1))] = ch_atn
        DCMlines.append('{:2} :  {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2}'.format(band,*out[:30]))
    return DCMlines
        
def adc_check(roach_list,dcmlines,ant_list='ant1-15',do_plot=False):
    ''' Perform a sequence of FEM settings, using ADC levels to 
        deduce optimum DCM attenuation settings for all 34 bands.
        This can also reveal problems in FEM or DCM hardware.
        TAKES ABOUT 17 MINUTES TO COMPLETE
        
        roach_list  a set of ROACH objects created with roach.py
        ant_list    a list of antennas in the form of a string,
                      e.g. "ant1-5 ant7" on which to adjust FEMs
                      Default is all antennas, and an empty string
                      means all antennas in current subarray.
        do_plot     if True, makes a summary plot of results
        
        Returns numpy arrays :
                adc_nosig[34, nroach, 4] (no-signal ADC levels)
                adc_ndoff[34, nroach, 4] (ADC levels for ND-OFF)
                adc_ndon [34, nroach, 4] (ADC levels for ND-ON)
    '''
    accini = stf.rd_ACCfile()
    acc = {'host': accini['host'], 'scdport':accini['scdport']}

    n = len(roach_list)
    adc_nosig = np.zeros((34,n,4),dtype='float')
    adc_ndoff = np.zeros((34,n,4),dtype='float')
    adc_ndon = np.zeros((34,n,4),dtype='float')
    # Set DCM state to standard values
    send_cmds(['DCMAUTO-OFF '+ant_list,'DCMATTN 12 12 '+ant_list],acc)
    # Set FEM attenuations to maximum
    send_cmds(['FEMATTN 15 '+ant_list],acc)
    # Cycle through bands to get "zero-input" ADC levels
    for band in range(34):
        acc_tune(band+1,acc)
        line = dcmlines[band+3]
        for ant in range(1,16):
            send_cmds(['DCMATTN'+line[ant*6-1:(ant+1)*6-1]+' ant'+str(ant)],acc)
            time.sleep(1)
        r.adc_levels(roach_list)
        for i,ro in enumerate(roach_list):
            adc_nosig[band,i] = ro.adc_levels
    # Set FEM attenuations to nominal
    send_cmds(['FEMATTN 0 '+ant_list],acc)
    # Cycle through bands to get "nd-on" ADC levels
    send_cmds(['ND-ON '+ant_list],acc)
    for band in range(34):
        acc_tune(band+1,acc)
        line = dcmlines[band+3]
        for ant in range(1,16):
            send_cmds(['DCMATTN'+line[ant*6-1:(ant+1)*6-1]+' ant'+str(ant)],acc)
            time.sleep(1)
        r.adc_levels(roach_list)
        for i,ro in enumerate(roach_list):
            adc_ndon[band,i] = ro.adc_levels
    # Cycle through bands to get "nd-off" ADC levels
    send_cmds(['ND-OFF '+ant_list],acc)
    for band in range(34):
        acc_tune(band+1,acc)
        line = dcmlines[band+3]
        for ant in range(1,16):
            send_cmds(['DCMATTN'+line[ant*6-1:(ant+1)*6-1]+' ant'+str(ant)],acc)
            time.sleep(1)
        r.adc_levels(roach_list)
        for i,ro in enumerate(roach_list):
            adc_ndoff[band,i] = ro.adc_levels
    if do_plot: 
        plot_adc_cal(roach_list, adc_nosig, adc_ndoff, adc_ndon)
    return adc_nosig, adc_ndoff, adc_ndon
