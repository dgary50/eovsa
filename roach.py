#
#   Routines that talk to the ROACH boards and KatADCs.
# 
# History:
#   2015-Jan-20  DG
#      Started this history log.  Really important change!  We have to
#      configure the KatADCs immediately after loading the bof file, to
#      avoid having the FPGA clock come up at twice the nominal speed.
#      This has now been added to the load() routine.
#   2015-Feb-17  DG
#      Change mcount start to 700 (trying to synchronize X packets).
#      Also made changes to configuration to send ROACH1 and ROACH2
#      packets to eth2, and ROACH3 and ROACH4 packets to eth3.  It was
#      odd-numbered ROACHes to eth2 and even-numbered to eth3.
#   2015-Apr-04  DG
#      Added adc_levels, which makes a plot of the ADC levels for
#      a list of roach objects on a single, nicely-formatted plot, and
#      labels the plot with standard deviation.
#   2015-Jun-16  DG
#      FTP to ACC now requires a username and password
#   2015-Jul-25  DG
#      I was storing the X and Y delays in reverse order!  It did not
#      matter when they were both the same, but now that they can 
#      differ, they have to be stored correctly.  Yikes.  I was also
#      getting them in reverse order to write into stateframe...
#   2015-Sep-02  DG
#      Slight update to print the correct expected mcount value (for
#      a 200 MHz clock).
#   2015-Oct-27  DG
#      Attempt to get three pairs of ROACHes working. Slight change to
#      print boffile status during initial connection.
#   2016-Jan-02  DG
#      Cleaned up some bugs that became apparent when testing the new
#      16-antenna correlator.  It was basically associated with IP
#      address assignments and also ARP table errors.
#   2016-Jan-03  DG
#      Added setup of multiplicative "equalization" coefficients for
#      registers eq_x0_coeff, in 16-antenna design (skipped if nbds
#      is not 8).
#   2016-Jan-22  DG
#      Changed "pcycle" wait time to 90 s if more than 6 roaches are
#      being rebooted (60 s was not quite long enough).
#   2016-Feb-10  DG
#      Added import of qdr, and calibration of QDR memory.  The design
#      uses three of them, so this adds about 2 minutes to the startup
#      time.  Probably I should spawn tasks so that the calibration step
#      will run in parallel.
#   2016-Feb-14  DG
#      Added tvg_state() function to turn test-vector generator on/off.
#   2016-Feb-20  DG
#      Added do_plot keyword to adc_levels(), and return of adc_levels
#      (standard deviation) for 4 channels in each roach.
#   2016-Feb-26  DG
#      Fix adc_levels() plotting bug, when only one roach is used.  Also
#      added "helper" routine to ROACH object, self.set_eq() to set
#      the equalizer coefficients.
#

import corr, qdr, struct, numpy, time, copy, sys
import urllib2, subprocess
from katcp import Message
from pwr_cycle import pwr_cycle, pwr_off

# Start of ROACH object definition
# ========================
class Roach():
    def __init__(self, roach_ip=None, boffile=None):

        # Some initializations to be added to ROACH object
        self.sdev = [0]*4
        self.sdevold = [0]*4
        self.db = [0]*4
        self.dbnew = [19]*4
        self.dbold = [0]*4
        self.roach_ip = None
        self.boffile = None
        self.fpga = None
        self.msg = None

        if roach_ip is None:
            roach_ip = 'roach1.solar.pvt'
        try:
            tok = roach_ip.split('.')
        except:
            # Roach IP address is malformed
            self.msg = 'Bad IP address'
            self.fpga = None
            return
        if tok[0][:-1] != 'roach':
            # Roach IP address not rational--must start with 'roach<n>'
            self.msg = 'Bad IP address--must start with "roach<n>"'
            self.fpga = None
            return
        elif len(tok) == 1:
            # Probably 'solar.pvt' was not included, so add it
            roach_ip += '.solar.pvt'

        # This is the handle for communicating with ROACH
        fpga = corr.katcp_wrapper.FpgaClient(roach_ip,timeout=3)
        # Allow time for connection
        if fpga.wait_connected(1):
            self.fpga = fpga
            self.roach_ip = roach_ip
        else:
            self.msg = 'Could not connect to ROACH'
            self.fpga = None
            return

        if boffile is None:
            # Read boffile name from ACC's ini file
            try:
                userpass = 'admin:observer@'
                f = urllib2.urlopen('ftp://'+userpass+'acc.solar.pvt/ni-rt/startup/acc.ini',timeout=0.5)
                lines = f.readlines()
                f.close()
                for line in lines:
                    if line.find('boffile = ') == 0:
                        boffile = line[10:].strip()
            except:
                msg = 'FTP of ftp://acc.solar.pvt/parm/acc.ini timed out.'
                boffile = None

        self.boffile = boffile

        # See if ROACH design is loaded
        try:
            devs = self.fpga.listdev()
        except:
            self.msg = 'ROACH design is not yet loaded.'
            return

        # All is good, so find out what antennas are connected
        self.ants,self.msg = self.get_ants(self.boffile)

        # Get initial sensor state
        self.get_sensor_dict()      # Defines self.sensors

#======= Roach.get_ants ========
    def get_ants(self,boffile=None):
        ''' Read configuration file based on boffile, and return the antennas
            assigned to this ROACH.  Also attaches contents of config file to
            self.cfg
        '''
        self.boffile = boffile
        ants = None
        cfg = None
        if boffile:
            # Read configuration file to find out what antennas are connected
            cfgfile = boffile[:-3]+'ini'
            try:
                userpass = 'admin:observer@'
                c = urllib2.urlopen('ftp://'+userpass+'acc.solar.pvt/parm/'+cfgfile,timeout=0.5)
                cfg = c.readlines()
                c.close()
            except:
                msg = 'FTP of '+cfgfile+' timed out.'
            if cfg:
                # Find out which ROACH we are
                rnum = int(self.roach_ip[5:6])
                roachants,msg = rd_ini(cfg,'antasgn',[8,2])
                if msg is 'int':
                    ants = roachants[rnum-1]
                    msg = 'Success'
        else:
            msg = 'No boffile defined'
        self.cfg = cfg
        return ants,msg

#======= Roach.load ========
    def load(self,boffile=None,config=True):
        ''' Load the ROACH with bitstream given by boffile.  If boffile
            is omitted, the user is prompted to enter the boffile from
            a list of boffiles.  Normally this will also call self.config()
            to configure the design, unless config is set to False.
        '''
        if self.fpga.is_connected():
            bofs = self.fpga.listbof()
            bofs.sort()
            if boffile:
                pass
            else:
                for i,bof in enumerate(bofs):
                    if bof[:10] == 'eovsa_corr': print i,':',bof
                bofn = raw_input('Enter number of boffile to load:')
                try:
                    boffile = bofs[int(bofn)]
                except:
                    print 'Entered value not a string'
                    self.msg = 'Count not load boffile--error in input.'
                    self.boffile = None
                    return
            try:
                bofn = bofs.index(boffile)
            except:
                self.msg = 'Count not load boffile--boffile not found.'
                self.boffile = None
                return
            try:
                status = 'Progdev failed.'
                status = self.fpga.progdev(boffile)
                self.boffile = boffile
                if status is 'ok':
                    self.msg = 'Success'
                else:
                    self.msg = 'Could not load boffile--'+status
            except:
                self.msg = 'Could not load boffile--'+status
        else:
            self.msg = 'Could not load boffile--client not connected'
            self.boffile = None
        if self.boffile:
            # Configure KATADC registers
            addr = [0x0000, 0x0001, 0x0002, 0x0003, 0x0009, 0x000A, 0x000B, 0x000E, 0x000F]
            val  = [0x7FFF, 0xBAFF, 0x007F, 0x807F, 0x03FF, 0x007F, 0x807F, 0x00FF, 0x007F]
            #if interleaved: val[4] = 0x23FF # Uncomment this line for interleaved mode
            for i in range(len(addr)):
                print('Setting ADC register %04Xh to 0x%04X' % (addr[i], val[i]))
                # Program both ZDOKs (this could be made smarter if needed).
                corr.katadc.spi_write_register(self.fpga, 0, addr[i], val[i])
                corr.katadc.spi_write_register(self.fpga, 1, addr[i], val[i])

        if self.boffile and config:
            self.config()

#======= Roach.config ========
    def config(self,dbg=False,print_arp=False):
        ''' Read config file and configure the firmware for operation.
            Leaves swreg_rst set to 1 (i.e. 10 GbE cores held in reset).
            Config files are in file '<boffile_name>.cfg'
        '''
        self.ants,self.msg = self.get_ants(self.boffile)
        if self.msg != 'Success':
            return
            
        bpow = 2**numpy.array([24,16,8,0])
        cfg = self.cfg
#        try:
#            c = urllib2.urlopen('ftp://acc.solar.pvt/parm/'+cfgfile,timeout=0.5)
#            cfg = c.readlines()
#            c.close()
#        except:
#            self.msg = 'Not configured--FTP of',cfgfile,'timed out.'
#            return
        # Find out which ROACH we are
        rnum = int(self.fpga.host[5:6])
        
        # Get ADC clock frequency [MHz]
        clk, self.msg = rd_ini(cfg,'adc clock')
        if clk is None:
            return
        clkHz = clk*1000000
        # Check that FPGA clock is 1/4 of clk
        fpga_clk = int(self.fpga.est_brd_clk() + 0.5)
        if abs(clk - fpga_clk*4) > 4.0:
            self.msg = 'FPGA clock '+str(fpga_clk)+' is not compatible with ADC clock '+str(clk)
        #    return

        # Set nominal KatADC attenuation
        attn, self.msg = rd_ini(cfg,'adcattn',[8,4])
        if attn is None:
            return
        db = (attn[rnum-1]*2).astype('int')
        if dbg:
            print 'ADC Atten',db,'packed =',hex((db*bpow).sum())
        else:
            self.fpga.write_int('swreg_adc_atten',(db*bpow).sum())

        # Set nominal KatADC enable
        en, self.msg = rd_ini(cfg,'adcen',[8,4])
        if en is None:
            return
        ena = en[rnum-1].astype('int')
        if dbg:
            print 'ADC Enable',ena,'packed =',hex((ena*bpow).sum())
        else:
            self.fpga.write_int('swreg_adc_en',(ena*bpow).sum())

        # Set Phase Switching Pattern
        psw, self.msg = rd_ini(cfg,'pspattern',[16,16])
        if psw is None:
            return
        psw1 = psw[[self.ants[0]-1,self.ants[1]-1]]
        psw1.shape = 32
        pow32 = 2**numpy.arange(31,-1,-1) 
        if dbg:
            print 'Phase Switching',psw1,'packed =',hex((psw1*pow32).sum())
        else:
            self.fpga.write_int('swreg_phase_switch',(psw1*pow32).sum())
        
        # Set Delays
        dla, self.msg = rd_ini(cfg,'delay',[16,1])
        if dla is None:
            return
        # Convert to integer delays in units of ADC clock
        dlas = 5000 + (dla[[self.ants[0]-1,self.ants[1]-1]]*1000./clk).astype('int')
        if dbg:
            print 'Delay0',dlas[0],'packed =',hex(dlas[0]*(2**16)+dlas[0])
            print 'Delay1',dlas[1],'packed =',hex(dlas[1]*(2**16)+dlas[1])
        else:
            self.fpga.write_int('swreg_d0',dlas[0]*(2**16)+dlas[0])
            self.fpga.write_int('swreg_d1',dlas[1]*(2**16)+dlas[1])

        # Set packet delay and fft shift
        fftshift, self.msg = rd_ini(cfg,'fft shift')
        if self.msg != 'int':
            self.msg = 'Incorrect format for fftshift'
            return
        pktdla, self.msg = rd_ini(cfg,'packet delay')
        if self.msg != 'int':
            self.msg = 'Incorrect format for packet delay'
            return
        if dbg:
            print 'Packet Delay',pktdla
            print 'FFT Shift',fftshift,'packed =',hex((2**16)*pktdla+fftshift)
        else:
            self.fpga.write_int('swreg_pkt_fft',(2**16)*pktdla+fftshift)

        # Set frequency sequence
        fseq, self.msg = rd_ini(cfg,'fseq',[1,50])
        if fseq is None:
            return
        fseq.shape = 50
        self.set_sequence(fseq)

        # Determine number of boards being run (if 8, means 16-ant correlator)
        nbds, self.msg = rd_ini(cfg,'n boards')
        if self.msg != 'int':
            self.msg = 'Incorrect format for "n boards"'
            return
        self.nboards = nbds

        # Set accumulation start and length (calculate from ADC clock)
        acc_start = int(0.001 * clkHz / 8192.) + 54   # ~1 ms delay
        acc_len   = int(0.019 * clkHz / 8192.) - 63   # ~19 ms duration
        if nbds == 8:
            # 16-ant correlator uses different units for x_acc_len (number of 256-sample blocks)
            # This is 7 (7*256 = 1792), but because of 0-based scheme it is set to 6
            x_acc_len = 6
        else:
            x_acc_len = acc_len

        if dbg:
            print 'Acc Length',acc_len
            print 'Acc Start',acc_start,'packed =',hex((2**16)*acc_len+acc_start)
        else:
            self.fpga.write_int('swreg_acc_len',(2**16)*acc_len + acc_start)

        # Set up board information
        act_bds, self.msg = rd_ini(cfg,'actual boards')
        if self.msg != 'int':
            self.msg = 'Incorrect format for "actual boards"'
            return
        if dbg:
            print 'Board ID',rnum-1
            print 'Actual Boards',act_bds,'packed =',hex((2**16)*(rnum-1) + act_bds)
        else:
            self.fpga.write_int('swreg_board_info',(2**16)*(rnum-1) + act_bds)

        # Load ROACH address and port information
        fx, self.msg = rd_ini(cfg,'fx',[2,8])
        if fx is None:
            return
        if dbg:
#            for i in range(nbds):
#                print 'Roach'+str(i)+' FX0 IP address',fx[0,i],'packed =',ip2int(fx[0,i])
#                print 'Roach'+str(i)+' FX1 IP address',fx[1,i],'packed =',ip2int(fx[1,i])
            if nbds == 2:
                # Case of 4-antenna prototype
                j = int(self.roach_ip[5:6])-1
                # If j is odd, use next lower index -- this only works for pairs of ROACHes
                j = (j/2)*2
                print 'Roach'+str(j)+' FX0 IP address','swreg_ip0_roach0',fx[0,j],'packed =',ip2int(fx[0,j])
                print 'Roach'+str(j)+' FX1 IP address','swreg_ip1_roach0',fx[1,j],'packed =',ip2int(fx[1,j])
                print 'Roach'+str(j+1)+' FX0 IP address','swreg_ip0_roach1',fx[0,j+1],'packed =',ip2int(fx[0,j+1])
                print 'Roach'+str(j+1)+' FX1 IP address','swreg_ip1_roach1',fx[1,j+1],'packed =',ip2int(fx[1,j+1])
            else:
                for j in range(8):
                    print 'Roach'+str(j)+' FX0 IP address','swreg_ip0_roach'+str(j),fx[0,j],'packed =',ip2int(fx[0,j])
                    print 'Roach'+str(j)+' FX1 IP address','swreg_ip1_roach'+str(j),fx[1,j],'packed =',ip2int(fx[1,j])
        else:
            if nbds == 2:
                # Case of 4-antenna prototype
                j = int(self.roach_ip[5:6])-1
                # If j is odd, use next lower index -- this only works for pairs of ROACHes
                j = (j/2)*2
                self.fpga.write_int('swreg_ip0_roach0',ip2int(fx[0,j]))
                self.fpga.write_int('swreg_ip1_roach0',ip2int(fx[1,j]))
                self.fpga.write_int('swreg_ip0_roach1',ip2int(fx[0,j+1]))
                self.fpga.write_int('swreg_ip1_roach1',ip2int(fx[1,j+1]))
            else:
                for j in range(8):
                    self.fpga.write_int('swreg_ip0_roach'+str(j),ip2int(fx[0,j]))
                    self.fpga.write_int('swreg_ip1_roach'+str(j),ip2int(fx[1,j]))
                
        fx_port, self.msg = rd_ini(cfg,'fx port')
        if self.msg != 'int':
            self.msg = 'Incorrect format for FX port'
            return
        if dbg:
            print 'FX Port',fx_port
        else:
            self.fpga.write_int('swreg_fx_udp_port',fx_port)

        # Load multipurpose control register
        pps, self.msg = rd_ini(cfg,'pps',[1,8])
        if self.msg is None:
            return
        mypps = pps[rnum-1]
        mypps_sel = 0
        if mypps == -1:
            # A value of -1 for mypps means use synthetic pps
            mypps = 0
            mypps_sel = 1
        poln, self.msg = rd_ini(cfg,'poln',[1,8])
        if self.msg is None:
            return
        mypoln = poln[rnum-1]
        adcprot, self.msg = rd_ini(cfg,'adcprot',[1,8])
        if self.msg is None:
            return
        myadcprot = adcprot[rnum-1]
        vgen, self.msg = rd_ini(cfg,'vgen',[1,8])
        if self.msg is None:
            return
        myvgen = vgen[rnum-1]
        pows = 2**numpy.array([31,30,29,28,23])
        vect = numpy.array([mypps_sel,mypps,mypoln,myadcprot,myvgen]).astype('int')
        if dbg:
            print 'Control Register vector [PPS_sel,PPS,Poln,ADCProt,VGen]',vect,'packed',hex((vect*pows).sum())
        else:
            self.fpga.write_int('swreg_ctrl',(vect*pows).sum())

        # Load X control register
        pows = 2**numpy.array([16,8,4,0])
        vect = numpy.array([x_acc_len,pktdla,rnum-1,mypoln])
        if dbg:
            print 'X Control Register vector [Acc Length,Pktdla,Board,Poln]',vect,'packed',hex((vect*pows).sum())
        else:
            self.fpga.write_int('swreg_ctrl_x',(vect*pows).sum())

        # Load X-engine output sequence
        xop, self.msg = rd_ini(cfg,'xop',[1,11])
        if dbg:
            print 'X-engine output sequence',xop
        else:
            for i,op in enumerate(xop):
                # Toggle in values and addresses by setting lsb to 1, then 0
                self.fpga.write_int('swreg_x_op_select',(2**16)*i + (2**4)*op + 1)
                self.fpga.write_int('swreg_x_op_select',0)

        # Read X and P interface assignments
        dpp_P, self.msg = rd_ini(cfg,'dpp_P',[1,8])
        if dpp_P is None:
            return
        dpp_X, self.msg = rd_ini(cfg,'dpp_X',[1,8])
        if dpp_X is None:
            return
        mac0, self.msg = rd_ini(cfg,'mac0')
        mac, self.msg = rd_ini(cfg,'mac',[2,1])
        if mac is None:
            return
        # Determine the subnet this ROACH is going to use for X and P packets
        subnet = int(dpp_P[rnum-1].split('.')[-2])
        if int(dpp_X[rnum-1].split('.')[-2]) != subnet:
            self.msg = 'Subnet for X and P packets for this ROACH do not agree'
            return

        # Load dpp address and port number
        dppn, self.msg = rd_ini(cfg,'dpp',[2,1])
        if dppn is None:
            return
        if dbg:
            print 'DPP IP Address',dppn[subnet-1],'packed =',ip2int(dppn[subnet-1])
        else:
            self.fpga.write_int('swreg_dpp_ip',ip2int(dppn[subnet-1]))
        dpp_port, self.msg = rd_ini(cfg,'dpp port')
        if self.msg != 'int':
            self.msg = 'Incorrect format for DPP port'
            return
        if dbg:
            print 'DPP Port',dpp_port
        else:
            self.fpga.write_int('swreg_dpp_port',dpp_port)

        # Manually fix Tx_P and Tx_X ARP table (for some reason tap_start does not work right)
        # Configuring 10GbE core Tx_P (transmitting P data)
        arp = self.fpga.read('Tx_P', 256 * 8, 0x3000)
        arp_tab = numpy.array(struct.unpack('>256Q', arp))
        for roach_id in range(8):
            idxP = int(dpp_P[roach_id].split('.')[-1])
            idxX = int(dpp_X[roach_id].split('.')[-1])
            if int(dpp_P[roach_id].split('.')[-2]) == subnet:
                arp_tab[idxP] = mac2int(mac0) + ip2int(dpp_P[roach_id])  #  MAC address for 10.0.x.n1 (ROACHn)
                arp_tab[idxX] = mac2int(mac0) + ip2int(dpp_X[roach_id])  #  MAC address for 10.0.x.n2 (ROACHn)
        idx = int(dppn[subnet-1].split('.')[-1])
        arp_tab[idx] = mac2int(mac[subnet-1])   # MAC address for 10.0.1.100 or 10.0.2.100 (dpp eth2 or eth3)
        self.fpga.config_10gbe_core('Tx_P', mac2int(mac0) + ip2int(dpp_P[rnum-1]),
                               ip2int(dpp_P[rnum-1]), dpp_port, arp_tab.tolist())
        if print_arp:
            self.fpga.print_10gbe_core_details( 'Tx_P', arp = True )


        # Configuring 10GbE core Tx_X (transmitting X data)
        arp = self.fpga.read('Tx_X', 256 * 8, 0x3000)
        arp_tab = numpy.array( struct.unpack('>256Q', arp))
        for roach_id in range(8):
            idxP = int(dpp_P[roach_id].split('.')[-1])
            idxX = int(dpp_X[roach_id].split('.')[-1])
            if int(dpp_P[roach_id].split('.')[-2]) == subnet:
                arp_tab[idxP] = mac2int(mac0) + ip2int(dpp_P[roach_id])  #  MAC address for 10.0.x.n1 (ROACHn)
                arp_tab[idxX] = mac2int(mac0) + ip2int(dpp_X[roach_id])  #  MAC address for 10.0.x.n2 (ROACHn)
        idx = int(dppn[subnet-1].split('.')[-1])
        arp_tab[idx] = mac2int(mac[subnet-1])   # MAC address for 10.0.1.100 or 10.0.2.100 (dpp eth2 or eth3)
        self.fpga.config_10gbe_core('Tx_X', mac2int(mac0) + ip2int(dpp_X[rnum-1]),
                               ip2int(dpp_X[rnum-1]), dpp_port, arp_tab.tolist())
        if print_arp:
            self.fpga.print_10gbe_core_details( 'Tx_X', arp = True )

        # Configuring 10GbE core Rx_Tx_FX0 (data transmission between F and X engines)
        arp = self.fpga.read( 'Rx_Tx_FX0', 256 * 8, 0x3000 )
        arp_tab = numpy.array( struct.unpack( '>256Q', arp ) )
        for roach_id in range(8):
            idx1 = int(fx[0,roach_id].split('.')[-1])
            idx2 = int(fx[1,roach_id].split('.')[-1])
            arp_tab[idx1] = mac2int(mac0) + ip2int(fx[0,roach_id])  #  MAC address for 10.0.0.x ROACH_id
            arp_tab[idx2] = mac2int(mac0) + ip2int(fx[1,roach_id])  #  MAC address for 10.0.0.x ROACH_id
        self.fpga.config_10gbe_core('Rx_Tx_FX0', mac2int(mac0) + ip2int(fx[0,rnum-1]),
                               ip2int(fx[0,rnum-1]), fx_port, arp_tab.tolist())
        if print_arp:
            self.fpga.print_10gbe_core_details( 'Rx_Tx_FX0', arp = True )

        # Configuring 10GbE core Rx_Tx_FX1 (data transmission between F and X engines)
        arp = self.fpga.read( 'Rx_Tx_FX1', 256 * 8, 0x3000 )
        arp_tab = numpy.array( struct.unpack( '>256Q', arp ) )
        for roach_id in range(8):
            idx1 = int(fx[0,roach_id].split('.')[-1])
            idx2 = int(fx[1,roach_id].split('.')[-1])
            arp_tab[idx1] = mac2int(mac0) + ip2int(fx[0,roach_id])  #  MAC address for 10.0.0.x ROACH_id
            arp_tab[idx2] = mac2int(mac0) + ip2int(fx[1,roach_id])  #  MAC address for 10.0.0.x ROACH_id
        self.fpga.config_10gbe_core('Rx_Tx_FX1', mac2int(mac0) + ip2int(fx[1,rnum-1]),
                               ip2int(fx[1,rnum-1]), fx_port, arp_tab.tolist())
        if print_arp:
            self.fpga.print_10gbe_core_details( 'Rx_Tx_FX1', arp = True )

        # Reset 10 GbE cores
        if not dbg:
            self.fpga.write_int( 'swreg_rst', 1 )
#            self.fpga.write_int( 'swreg_rst', 0 )  # Commented out to hold 10  GbE cores in reset state

        # Set nominal values for multiplicative coefficients, in the case of the 16-antenna correlator
        if nbds == 8:
            for i in range(4): self.set_eq(xn=i,coeff=10)
#            coefficients = numpy.ones(2**13, dtype='>I') << (10 + 6)
#            for i in range(4): self.fpga.write('eq_x'+str(i)+'_coeffs',coefficients.tostring())
            # Calibrate qdrs
            devs = self.fpga.listdev()
            # Look for string 'qdr' in list of devs
            qdrs = [s for s in devs if 'qdr' in s]
            if qdrs != []:
                for qdrstr in qdrs:
                    if qdrstr[5:] == 'ctrl':
                        q = qdr.Qdr(self.fpga,qdrstr[:4])
                        print 'QDR',qdrstr[:4],'calibration status:',q.qdr_cal()
                

        self.msg = 'Success'

#======= Roach.set_delays ========
    def set_delays(self,dlas):
        ''' Sets integer delays (called as often as once each second)
            Order of delays is [dlax[0],dlay[0],dlax[1],dlay[1]]
        '''
        if self.fpga.is_connected():
            self.fpga.write_int('swreg_d0',dlas[0]*(2**16)+dlas[1])
            self.fpga.write_int('swreg_d1',dlas[2]*(2**16)+dlas[3])
            self.msg = 'Success'
        else:
            self.msg = 'Could not set delays--client not connected'

#======= Roach.get_delays ========
    def get_delays(self):
        ''' Gets integer delays from ROACH (called as often as once each second)
            Result is attached to the ROACH object (self.delays)
            Order of delays returned is [dlax[0],dlay[0],dlax[1],dlay[1]]
        '''
        if self.fpga.is_connected():
            d0 = self.fpga.read_int('swreg_d0')
            d1 = self.fpga.read_int('swreg_d1')
            self.delays = [d0 >> 16, d0 & 0xffff, d1 >> 16, d1 & 0xffff]
            self.msg = 'Success'
        else:
            self.msg = 'Could not get delays--client not connected'

#======= Roach.set_sequence ========
    def set_sequence(self,bands):
        ''' Sets frequency sequence (called once each scan change)
            The bands is an array of 0-based band numbers
        '''
        if self.fpga.is_connected():
            for i in range(10):
                idx = i * 5
                sw_reg_name = 'swreg_sky_freq' + str(i)
                sw_reg_value = bands[idx:idx+5]*(2**numpy.array([0,6,12,18,24]))
                self.fpga.write_int(sw_reg_name, sw_reg_value.sum())
            self.msg = 'Success'
        else:
            self.msg = 'Could not set frequency sequence--client not connected'

#======= Roach.get_katadc_dict ========
    def get_katadc_dict(self):
        ''' Create or update dictionary of temperatures of KatADCs.  This
            was part of sensors, but it takes too long to call very
            often, so it is split out to call once per scan.
        '''
        keys = []
        values = []
        if self.fpga.is_connected():
            keys.append('temp.ambient0')
            val = corr.katadc.get_ambient_temp(self.fpga,0)
            values.append(val)
            keys.append('temp.ambient0.status')
            if val < 10 or val > 40:
                values.append('error')
            else:
                values.append('nominal')
            keys.append('temp.adc0')
            val = corr.katadc.get_adc_temp(self.fpga,0)
            values.append(val)
            keys.append('temp.adc0.status')
            if val < 30 or val > 99:
                values.append('error')
            else:
                values.append('nominal')
            keys.append('temp.ambient1')
            val = corr.katadc.get_ambient_temp(self.fpga,1)
            values.append(val)
            keys.append('temp.ambient1.status')
            if val < 10 or val > 40:
                values.append('error')
            else:
                values.append('nominal')
            keys.append('temp.adc1')
            val = corr.katadc.get_adc_temp(self.fpga,1)
            values.append(val)
            keys.append('temp.adc1.status')
            if val < 30 or val > 99:
                values.append('error')
            else:
                values.append('nominal')
            self.katadc = dict(zip(keys,values))
            self.msg = 'Success'
        else:
            self.msg = 'Could not update katadc sensors--client not connected'

#======= Roach.get_sensor_dict ========
    def get_sensor_dict(self):
        ''' Create or update sensor dictionary
        '''
        if self.fpga.is_connected():
            # Multiplicative factor to apply to each.  This assumes specific
            # ROACH2 sensor sampling--this could change if ROACH
            # tcpborphserver3 is updated
            factor = [0]+[0.001]*3+[1]*4+[0.001]*18
            keys = []
            values = []
            reply, sensors = self.fpga.blocking_request(Message.request('sensor-list'))
            if reply.arguments[0] == 'ok':
                # Got sensor list okay
                n = int(reply.arguments[1])
                reply, vals = self.fpga.blocking_request(Message.request('sensor-value'))
                self.sensor_dict = {}
                if reply.arguments[0] == 'ok':
                    # Got sensor values okay
                    if n == int(reply.arguments[1]):
                        # The numbers of sensors and values agree!
                        # Skip first "sensor" which is the mode
                        for i in range(1,n):
                            name = sensors[i].arguments[0][4:]
                            keys.append(name)
                            values.append(int(vals[i].arguments[-1])*factor[i])
                            keys.append(name+'.status')
                            values.append(vals[i].arguments[-2])
                    else:
                        self.sensors = {}
                        self.msg = 'Could not init sensors--names,values are different lengths'
                        return
                else:
                    self.sensors = {}
                    self.msg = 'Could not init sensors:',reply.arguments[0]
                    return
            else:
                self.sensors = {}
                self.msg = 'Could not init sensors:',reply.arguments[0]
                return
            
            self.sensors = dict(zip(keys,values))
            self.msg = 'Success'
        else:
            self.msg = 'Could not update sensors--client not connected'

#======= Roach.get_xin_data ========
    def get_xin_data(self):
        impol1 = numpy.zeros((4,1024),dtype='byte')
        repol1 = numpy.zeros((4,1024),dtype='byte')
        impol0 = numpy.zeros((4,1024),dtype='byte')
        repol0 = numpy.zeros((4,1024),dtype='byte')
        raw = numpy.zeros((4,1024),dtype='int')
        for i in range(4):
            self.fpga.write_int('swreg_snap_select',i * (2**10))
            self.fpga.write_int('x_in_ctrl',0)
            # self.fpga.write_int('x_in_ctrl',7)
            # Per Dave M, can't write 7 here, has to be one
            self.fpga.write_int('x_in_ctrl',1)
            time.sleep(0.02)
            data = self.fpga.read('x_in_bram',4096,0)
            raw[i,:] = numpy.array(struct.unpack('>1024I',data))
            for j,val in enumerate(raw[i,:]):
                impol1[i,j] = val >> (28-16)
                repol1[i,j] = val >> (24-16) & 15
                impol0[i,j] = val >> (20-16) & 15
                repol0[i,j] = val >> (16-16) & 15
        self.xin = {'raw':raw,'impol1':impol1,'repol1':repol1,'impol0':impol0,'repol0':repol0}

#======= Roach.get_attn ========
    def get_attn(self, sdev_target=30.0, grab=None):
        ''' Get the attenuation settings for both channels of both ADCs
            associated with this board, calculate the standard deviation
            of the digitized signal on each, and use it to calculate a
            better attenuation to achieve sdev_target. Optionally grab
            an ADC data snapshot using grab=chan (0,1,2 or 3).
        '''
        for chan in range(4):
            # Set up to snap from ADC snap block, and apply signal (positive-
            # going pulse from 0 to 1) to trigger at next clock cycle.
            self.fpga.write_int('swreg_snap_select',chan)
            self.fpga.write_int('adc_data_adc_ctrl',0)
            self.fpga.write_int('adc_data_adc_ctrl',7)
            time.sleep(0.02)
            # Read data (8192 samples) from BRAM snap block
            data = self.fpga.read('adc_data_adc_bram',2048*4,0)
            # Unpack data to numpy array and get standard deviation
            blah = numpy.array(struct.unpack('>8192b',data))
            if grab == chan:
                self.data =  blah
            sdev1 = blah.std()
            # Ensure that a malfunctioning or disconnected ADC does not force
            # all attenuation to be removed.
            if sdev1 < 1:
                sdev1 = 3
            self.sdev[chan] = sdev1
            # Calculate the attenuation change needed to reach sdev_target
            # Note that these are attenuator units of 0.5 dB steps
            dbinc = 20*numpy.log(sdev1/sdev_target)
            # Read current attenuation from ADC attenuation register and
            # extract it from its packed form
            atnreg = self.fpga.read_int('swreg_adc_atten',0)
            shft = 8*(3-chan)
            db = (atnreg & (0xff<<shft))>>shft
            self.db[chan] = db/2.0
            # Determine new attenuation value (in 0.5-dB step units), but limit
            # to range 10-63 (5 dB to 31.5 dB)
            dbnew = int(db + dbinc)
            if dbnew > 63:
                dbnew = 63  # 31.5 dB
            elif dbnew < 10:
                dbnew = 10  # 5 dB
            self.dbnew[chan] = dbnew/2.0

#======= Roach.set_attn ========
    def set_attn(self, sdev_target=30.0, update=False):
        ''' Set the attenuation settings for both channels of both ADCs
            associated with this board.  If update is False or not supplied,
            use the self.dbnew settings.
        '''
        # Shift bits according to channel
        #       Chan 0  Chan 1  Chan 2  Chan 3
        shft = [24,     16,     8,      0]
        if update:
            self.get_attn(sdev_target)
        try:
            dbnew = self.dbnew
        except:
            # Looks like dbnew is not yet set, so return an error
            self.msg = 'DBnew not set--use get_attn() or set_attn(update=True)'
            return self
        # Set attenuation register bits (note the units are 0.5 dB steps)
        atnreg = 0
        for chan in range(4):
            atnreg += (int(dbnew[chan]*2)<<shft[chan])
        # Write it back to the ROACH
        self.fpga.write_int('swreg_adc_atten',atnreg)
        # Wait for things to settle, then grab a new sample and report back
        time.sleep(0.1)
        # Save old values for comparison
        self.sdevold = copy.copy(self.sdev)
        self.dbold = copy.copy(self.db)
        # Read new values resulting from attenuation change
        self.get_attn(sdev_target)

#======= Roach.get_xdata ========
    def get_xdata(self):
        ''' Grab X-in BRAM data, for diagnostic purposes.
        '''
        xdata = numpy.zeros([4,4,1024],'float')
        for chan in range(4):
            # Set up to snap from x-engine snap block, and apply signal
            # (positive-going pulse from 0 to 1) to trigger at next clock cycle.
            self.fpga.write_int('swreg_snap_select',chan<<20)
            self.fpga.write_int('x_in_ctrl',0)
            self.fpga.write_int('x_in_ctrl',1)
            time.sleep(0.02)
            # Read data (4096 samples) from BRAM snap block
            data = self.fpga.read('adc_data_adc_bram',4096,0)
            # Unpack data to bit-packed numpy array
            blah = numpy.array(struct.unpack('>1024I',data))
            # Each value is in lower 2-bytes, packed in 4-bit parts
            for i,val in enumerate(blah):
                xdata[chan,0][i] = ((val & 0x0000f000) >> 12)
                xdata[chan,1][i] = ((val & 0x00000f00) >> 8)
                xdata[chan,2][i] = ((val & 0x000000f0) >> 4)
                xdata[chan,3][i] = ((val & 0x0000000f))
        neg = numpy.where(xdata > 8.0)
        xdata[neg] = xdata[neg] - 15
        xdata = xdata/8.
        self.xdata = xdata

#======= Roach.get_xdata ========
    def tvg_state(self,state='on',verbose=False):
        ''' Turn test-vector generator on or off.
        '''
        current_val = self.fpga.read_uint('swreg_ctrl')
        if verbose: print 'Initial state :',binary_repr(2**32+current_val)[1:]
        if state == 'on':
            # Set test-vector-generator bit to 1
            current_val = current_val | (2**32 >> 9)
        else:
            # Set test-vector-generator bit to 0
            current_val -= current_val & (2**32 >> 9)
        self.fpga.write_int('swreg_ctrl', current_val)
        if verbose: print 'After tvg set :',binary_repr(2**32+current_val)[1:]
        # Toggle tvg enable, ensuring that it starts at 0
        current_val -= current_val & (2**32 >> 12)  # Zeros the tvg enable bit
        self.fpga.write_int('swreg_ctrl', current_val)
        if verbose: print 'TVG enable off:',binary_repr(2**32+current_val)[1:]
        time.sleep(0.1)
        current_val |= (2**32 >> 12)
        self.fpga.write_int('swreg_ctrl', current_val)
        if verbose: print 'TVG enable on :',binary_repr(2**32+current_val)[1:]
        
#======= Equalizer Gain ========
    def set_eq(self,file=None,xn=0,ifb=None,coeff=1.0+0j,update=False):
        ''' Given a list of roach objects, and optionally a source of gains,
            set the equalizer gains.
            
            file    If not None, is the name of a file that contains 
                      34 x 128 x 4 complex gains, as one set of 128 values
                      for each of the 34 IF bands, for each of x0, x1, x2, x3.  
                      If not None, the other keywords, if present, are
                      ignored.
            xn      The equalizer to set, either 0, 1, 2 or 3
            ifb     The IF band to set, ranging from 1 to 34. If None,
                       all IF bands are set to the same values
            coeff   The complex value(s) to set the coefficients to.  If an array,
                       it should have length 128.  If not, all 128 values
                       for target IF band are set to the same value.
            update  A boolean value--if True, multiplies coefficient to
                       what is there, otherwise overwrites (default)
        '''

        if file is None:
            eqname = 'eq_x'+str(xn)+'_coeffs'
            # Read equalizer coefficients buffer from ROACH
            buf = self.fpga.read(eqname,8192*4)
            # Unpack coefficients to an array of 16-bit integers.
            vals = numpy.array(struct.unpack('>16384H',buf),dtype='>H')
            # Complex representation as numpy complex array (length 8192)
            cvals = vals[::2] + 1j*vals[1::2]
            # Convert supplied coeff to complex value(s) scaled by 2**6
            coeff *= 64+0j
            if type(coeff) == complex:
                # Single value given, so set all 128 values for the band
                # to the same value
                newvals = numpy.ones(128)*(int(numpy.real(coeff)) + 1j*int(numpy.imag(coeff)))
            elif type(coeff) == numpy.ndarray:
                if len(coeff) == 128:
                    newvals = numpy.real(coeff).astype('int') + 1j*numpy.imag(coeff).astype('int')
                else:
                    print 'Error in length of coeff parameter. Must be scalar or length 128.'
                    print 'No update performed.'
                    self.msg = 'Error in length of coeff parameter'
                    return
            else:
                    print 'Error in type of coeff parameter. Must be scalar numpy.ndarray.'
                    print 'No update performed.'
                    self.msg = 'Error in type of coeff parameter'
                    return
            # At this point, I have current coefficients cvals, and new coefficients newvals,
            # both in complex form.
            if ifb is None:
                # Apply to all IF bands
                if update:
                    # Interpret supplied coeff as multiplicative factor
                    for i in range(64):
                        cvals[i*128:(i+1)*128] *= newvals
                else:
                    # Interpret supplied coeff as new values to replace existing
                    for i in range(64):
                        cvals[i*128:(i+1)*128] = newvals
            else:
                # Apply only to given IF band
                if update:
                    # Interpret supplied coeff as multiplicative factor
                    cvals[(ifb-1)*128:ifb*128] *= newvals
                else:
                    # Interpret supplied coeff as new values to replace existing
                    cvals[(ifb-1)*128:ifb*128] = newvals
        else:
            print 'File not yet implemented'
            self.msg = 'File not yet implemented'
            return
        # Convert complex values in cvals to packed buffer, and write to ROACH
        coefficients = numpy.zeros(16384,dtype='>H')
        coefficients[::2] = numpy.real(cvals).astype('>H')
        coefficients[1::2] = numpy.imag(cvals).astype('>H')        
        self.fpga.write(eqname,coefficients.tostring())
        self.msg = 'Success'
   
#======= arm ========
def arm(roach_list=None):
        ''' Arm all ROACH boards as simultaneously as possible, just
            after current second tick.  Arming an already armed board has
            no effect.
        '''
        # Do some sanity checks on roach list
        if roach_list is None:
            print 'Must provide a list of connected ROACH objects'
            return
        if not isinstance(roach_list,list):
            print 'Must provide a list of connected ROACH objects'
            return
        if not isinstance(roach_list[0],Roach):
            print 'List is not ROACH object list'
            return

        # Check connection to each board prior to arming
        ngood = 0
        for roach in roach_list:
            if roach.fpga.is_connected():
                ngood += 1
            else:
                print 'Roach',roach.roach_ip,'not connected!'        
        if ngood != len(roach_list):
            print 'At least one ROACH not connected, so none were armed'
            return
        
        # Wait for turn of second (+ 200 ms) to ensure arming well
        # before next 1 pps tick
        time.sleep(1 - (time.time() % 1) + 0.1)
        for roach in roach_list:
            roach.fpga.write_int('swreg_arm',3)
        time.sleep(0.1)
        for roach in roach_list:
            roach.fpga.write_int('swreg_arm',0)

        # Should be armed.  Wait for next 1 PPS pulse and check result
        time.sleep(2.0)
        
        # Check if synchronized
        for roach in roach_list:
            sync = roach.fpga.read_int('sync')
            if sync:
                print roach.roach_ip,'synced =',True
            else:
                print roach.roach_ip,'synced =',False

def reload(roach_list=None,pcycle=False):
    ''' This routine reloads and arms a list of boards given by either
        their ip addresses, or as a list of ROACH objects, with option
        to power-cycle the boards. It assumes that the boffile to use
        is the one in the acc.ini file, or in the roach objects themselves.
    '''
    # Do some sanity checks on roach list
    ips = False
    if roach_list is None:
        print 'Must provide a list of ip addresses or ROACH objects'
        return
    if not isinstance(roach_list,list):
        print 'Must provide a list of ip addresses or ROACH objects'
        return
    if not isinstance(roach_list[0],Roach):
        # This is presumably a list of ip addresses
        if not isinstance(roach_list[0],str):
            print 'Roach list is neither ROACH objects nor strings.'
            return
        else:
            ips = True

    # If pcycle is True, loop over roach boards and cycle their power
    if pcycle:
        print 'Power cycling',roach_list
        for roach in roach_list:
            if ips:
                rnum = roach[5:6]
            else:
                rnum = roach.roach_ip[5:6]
            result = pwr_cycle('pduroach.solar.pvt',rnum)
            if not result:
                print 'Connection timed out. Try again.'
                return
        # All ROACH boards have been power cycled.
        if len(roach_list) < 7:
            # Wait 60 s for bootup before proceeding
            print 'Sleeping 60 seconds while ROACH(es) bootup completes'
            time.sleep(60)
        else:
            # Wait 90 s for bootup before proceeding
            print 'Sleeping 90 seconds while ROACH(es) bootup completes'
            time.sleep(90)

    # If ips is True, create ROACH objects from ip addresses
    if ips:
        # Loop over roach list and replace ip address with ROACH object
        for i,roach in enumerate(roach_list):
            roach_list[i] = Roach(roach)
            print roach_list[i].roach_ip,'boffile status',roach_list[i].msg

    # Now loop over ROACH object list to (re)load boffile
    for i, roach in enumerate(roach_list):
        roach.load(boffile=roach_list[i].boffile)
        print roach.roach_ip,'status',roach.msg
        
    # All set, so arm the boards
    arm(roach_list)
    # Sleep for two seconds, then arm X-engine
    #time.sleep(2)
    # Take the 10 GbE cores out of reset
    for roach in roach_list:
        roach.fpga.write_int( 'swreg_rst', 0 )

    for roach in roach_list:
        ctrl = roach.fpga.read_int('swreg_ctrl')
        print 'SWREG_CTRL is originally:',ctrl
        roach.fpga.write_int('swreg_ctrl',ctrl+1)
        ctrl = roach.fpga.read_int('swreg_ctrl')
        print 'SWREG_CTRL is now:',ctrl

    # Sleep for 2 s, then get current mcount on first roach
    time.sleep(2)
    mc = roach_list[0].fpga.read_int('rx_mcount_fx0',0)
    mcstart = 700
    if len(roach_list) == 8:
        mcstart = 350
    print 'MCount is:',mc,'Future MCount should be:',mc+mcstart
    # Set mcount to desired start value
    mc = mcstart
    # Set mcount_start in all roaches for 16 s from now
    for roach in roach_list:
        #roach.fpga.write_int('vacc_target_mcnt',mc)
        roach.fpga.write_int('swreg_mcount_start',mc)
    time.sleep(16)
    for i,roach in enumerate(roach_list):
        print 'MCount for roach',i,roach.fpga.read_int('rx_mcount_fx0',0)

# Helper routines that are not part of the ROACH object
#======= rd_ini ========
def rd_ini(cfg,key,shape=None):
    ''' Given the lines of a ROACH configuration file, find and return
        the entity requested in the key string e.g. 'fft shift'. The shape
        item says how many lines to read, and how many entries on each line,
        e.g. [8,4] means read 8 lines with 4 values each (adding 1-8 to
        end of key for keyword).  All entities are returned as lists.
    '''
    # See if key exists in the file
    line = next((x for x in cfg if x.find(key) != -1),None)
    if line is None:
        return None,'Key not found '+key
    idx = cfg.index(line)
    
    def single_val(line):
        # Just a single number or string is expected
        line = line.strip()
        chars = line[line.find('=')+1:].strip()
        try:
            # See if the value can be interpreted as an integer
            val = int(chars)
            return chars,'int'         # Return integer value
        except:
            try:
                # See if the value can be interpreted as a float
                val = float(chars)
                return chars,'float'     # Return float value
            except:
                return chars,'string'   # Return as string

    if shape is None:
        # Just a single number or string is expected
        val, vartype = single_val(line)
        if vartype == 'string':
            return val,vartype
        elif vartype == 'int':
            return int(val),vartype
        else:
            return float(val),vartype
    items = []  # Will just append linearly and set shape at end
    for i in range(shape[0]):
        # This is one or more lines
        if shape[0] == 1:
            newkey = key
        else:
            newkey = key+str(i+1)
        # See if newkey exists in the file
        line = next((x for x in cfg if x.find(newkey) != -1),None)
        if line is None:
            return None,'Incorrect shape--too few lines '+newkey
        # Got this line, so parse for ',' separators
        tok = line.strip().split(',')
        if len(tok) == shape[1]:
            # Get first number in line, and variable type
            val,vartype = single_val(tok[0])
            items.append(val.strip(' '))
            for item in tok[1:]:
                items.append(item.strip(' '))
        else:
            return None,'Incorrect shape--too few values in line '+newkey
    values = numpy.array(items)
    # Reduce to remove unit dimensions, e.g. [8,1] => 8, or [1,8] => 8
    if shape[0] == 1:
        values.shape = shape[1]
    elif shape[1] == 1:
        values.shape = shape[0]
    else:
        values.shape = shape
    if vartype == 'string':
        return values,vartype
    else:
        return values.astype(vartype),vartype

#======= ip2int ========
def ip2int(ip_string):
    ''' Convert a string like 192.168.24.20 to an integer
    '''
    bpow = 2**numpy.array([24,16,8,0])
    tok = numpy.array(ip_string.split('.')).astype('int')
    return (bpow*tok).sum()
    
#======= mac2int ========
def mac2int(mac_string):
    ''' Convert a string like 00:60:dd:44:b3:f4 to an integer
    '''
    bpow = 2**numpy.array([40,32,24,16,8,0])
    tok = mac_string.split(':')
    intar = numpy.array([int(a,16) for a in tok])
    return (bpow*intar).sum()

#======= ADC_Levels ========
def adc_levels(ro,do_plot=False):
    ''' Given a list of roach objects, grab and display the current
        ADC level on each, optionally plotting the result.
    '''
    import matplotlib.pylab as plt
    n = len(ro)
    if do_plot:
        f, ax = plt.subplots(n,4)
        f.set_size_inches(10,2.5*n, forward=True)
    chans = ['X','Y']
    # Loop over roaches in list
    for j,rn in enumerate(ro):
        rn.adc_levels = []  # Clear any existing ADC level for this ROACH
        # Loop over the four channels on each roach
        for i in range(4):
            ant = rn.ants[i / 2]
            chan = chans[i % 2]
            rn.get_attn(grab=i)
            rn.adc_levels.append(int((rn.data).std()*10000. + 0.5)/10000.) # Append ADC level to list, with 0.0001 resolution
            if do_plot:
                if n == 1:
                    ax[i].cla()
                    ax[i].plot(rn.data,'.')
                    ax[i].tick_params(axis='both',which='major',labelsize=10)
                    ax[i].set_ylim([-128,128])
                    ax[i].set_xlim([0,4096])
                    ax[i].xaxis.set_ticks([0, 2048, 4096])
                    ax[i].text(100,100,'Roach'+rn.roach_ip[5:6]+' Ant '+str(ant)+chan+':',fontsize=10)
                    ax[i].text(4000,100,str((rn.data).std())[:5],fontsize=10,horizontalalignment='right')
                else:
                    ax[j,i].cla()
                    ax[j,i].plot(rn.data,'.')
                    ax[j,i].tick_params(axis='both',which='major',labelsize=10)
                    ax[j,i].set_ylim([-128,128])
                    ax[j,i].set_xlim([0,4096])
                    ax[j,i].xaxis.set_ticks([0, 2048, 4096])
                    ax[j,i].text(100,100,'Roach'+rn.roach_ip[5:6]+' Ant '+str(ant)+chan+':',fontsize=10)
                    ax[j,i].text(4000,100,str((rn.data).std())[:5],fontsize=10,horizontalalignment='right')
        rn.adc_levels = numpy.array(rn.adc_levels)  # Convert ADC levels to numpy array
    if do_plot: 
        plt.show()
        plt.tight_layout()
            

 
    

#======= rmon =======
def rmon():
    ''' Read contents of DPP_PACKET_TEST file, which is output once per
        minute on the dpp, as well as the Myricom interface packet counts.
    '''
    # DPP_PACKET_TEST file
    p = subprocess.Popen(['ssh','user@dpp','cat','DPP_PACKET_TEST.txt'],stdout=subprocess.PIPE)
    lines = p.stdout.readlines()
    vals = [0,0]
    if len(lines) > 10:
        vals = numpy.array(lines[2].split()).astype('int')
    elif len(lines) > 3:
        vals = numpy.array(lines[3].split()).astype('int')
        
    # Myricom interface board packet counts (runs dpp_eth_mon.py on dpp)
    p = subprocess.Popen(['ssh','user@dpp','python','dpp_eth_mon.py'],stdout=subprocess.PIPE)
    lines = p.stdout.readlines()
    return numpy.append(vals,numpy.array(lines[0].strip().strip('[]').split(',')).astype('int'))
    
