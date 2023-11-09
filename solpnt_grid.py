import dump_tsys
from calibration import skycal_anal
import dbutil
import stateframe
import urllib2
from util import Time, common_val_idx

cursor = dbutil.get_cursor()
solpntdict = dbutil.get_dbrecs(cursor,dimension=15,timestamp=Time('2023-10-13 17:40:22'),nrecs=300)
solpntdict = dbutil.get_dbrecs(cursor,dimension=15,timestamp=Time('2023-10-13 17:40:22'),nrecs=1700)
antlist = np.array([0,1,2,3,4,5,6,7,8,9,10,11,12])
ra = (solpntdict['Ante_Cont_RAVirtualAxis'][:,antlist]*np.pi/10000./180.).astype('float')
dec = (solpntdict['Ante_Cont_DecVirtualAxis'][:,antlist]*np.pi/10000./180.).astype('float')
rao = (solpntdict['Ante_Cont_RAOffset'][:,antlist]).astype('float')
deco = (solpntdict['Ante_Cont_DecOffset'][:,antlist]).astype('float')
times = solpntdict['Timestamp'][:,0].astype('int64').astype('float')
outdict = stateframe.azel_from_sqldict(solpntdict)
trk = np.logical_and(outdict['dAzimuth'][:,antlist]<0.0020,outdict['dElevation'][:,antlist]<0.0020)
pnt = {'Timestamp':times[0],'tstamps':times,'antlist':antlist,'trjfile':'SOLPNT_GRID.TRJ',
             'ra':ra,'dec':dec,'rao':rao,'deco':deco,'trk':trk}
userpass = 'admin:observer@'
trjfile = urllib2.urlopen('ftp://'+userpass+'acc.solar.pvt/parm/'+pnt['trjfile'])
trjlines = trjfile.readlines()
trjfile.close()
trjrao = []
trjdeco = []
for line in trjlines:
    rao, deco, junk = line.split(' ')
    trjrao.append(int(rao))
    trjdeco.append(int(deco))
antidx = range(len(antlist))
mask = np.zeros((1700,len(antlist),len(trjlines)),dtype='bool')
for i in range(1700):
    # Loop over each RAO,DECO position and calculate "good data" mask for each
    # desired position
    for j in range(len(trjlines)):
        # True if RAO and DECO matches this line
        m = np.logical_and(pnt['rao'][i,antidx] == trjrao[j],pnt['deco'][i,antidx] == trjdeco[j])
        # True if above AND is tracking
        mask[i,:,j] = np.logical_and(m,pnt['trk'][i,antidx])
ra0 = np.median(pnt['ra'])
dec0 = np.median(pnt['dec'])
proc = {'Timestamp':pnt['Timestamp'],'tstamps':pnt['tstamps'],'antlist':antlist,
            'ra':pnt['ra'][:,antidx],'dec':pnt['dec'][:,antidx],'ra0':ra0,'dec0':dec0,
            'rao':np.array(trjrao)*np.cos(dec0), 'deco':np.array(trjdeco), 
            'mask':mask}
trange = Time([pnt['Timestamp'],pnt['Timestamp']+1700.],format='lv')
skycal = skycal_anal(trange[0], do_plot=False, last=True, desat=True)
otp = dump_tsys.rd_miriad_tsys_16(trange, tref=trange[0], skycal=skycal, desat=True)

nf = 451
tsecsql = proc['tstamps']
# The otp times, if from miriad, are slightly off from exact times, so round to the ms
tsecdata = (Time(otp['ut_mjd'],format='mjd').lv + 0.001).astype('int64').astype('float')
sqlidx, dataidx = common_val_idx(tsecsql, tsecdata)
offsets = [-50000,-20000,-10000,-5000,-2000,-1000,0,1000,2000,5000,10000,20000,50000]
xidx = np.zeros((169),dtype=int)
yidx = np.zeros((169),dtype=int)
for i,o in enumerate(offsets):
    blah, = where(abs((proc['rao']-o)) < 500)
    xidx[blah] = i
    blah, = where(abs((proc['deco']-o)) < 500)
    yidx[blah] = i
tsys0 = otp['tsys'][:13,0,:,dataidx]               # Size 13, 451, nt
tsys1 = otp['tsys'][:13,1,:,dataidx]
m = proc['mask'][sqlidx]   # Size nt, 13, 169
images = np.zeros((13,13,13,2,451),dtype=float64)
for i in range(169):
    for iant in range(13):
        good, = np.where(m[:,iant,i])
        for j in range(nf):
            images[yidx[i],xidx[i],iant,0,j] = np.nanmedian(tsys0[good,iant,j])
            images[yidx[i],xidx[i],iant,1,j] = np.nanmedian(tsys1[good,iant,j])
x = proc['rao'][:13]
y = proc['deco'][::-13]
im = pcolormesh(x,y,images[:,:,10,1,200])