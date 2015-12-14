#!/usr/bin/env python

from mirexec import *
import os,sys,shutil

if len(sys.argv) < 2:
    print 'Usage: ./gen_full_disk.py stem unna, e.g.'
    print './gen_full_disk.py 2.000GHZ_2.40_B4_RCP un'
    print 'for MODEL_2.000GHZ_2.40_B4_RCP uniform weighting'
    print 'Default is natural (robust = 0.5'
    exit()

stem = sys.argv[1]
robust = 0.5
wlab = 'NA'
if sys.argv[2] == "un":
    robust = -0.5
    wlab = 'UN'

out = stem+'_'+wlab
model = 'MODEL_'+stem
disk = 'DISK_'+stem[:13]
freq = stem[:4]

# generate visibility file with just noise, baseunit converts to meters
shutil.rmtree('uvgen.tmp',ignore_errors=True)
TaskUVGen(source='empty', ant='ant.tmp', baseunit=-3.33564, corr='1,1,0,100',
          freq=freq+',0.0', radec='0.0,10.0', harange='-6.5,6.5,.083333',
          ellim=10, lat=37, systemp=1500, out='uvgen.tmp').run()

TaskUVPlot(vis='uvgen.tmp', axis='uc,vc', device='2/xs', nobase=True, 
          equal=True, nxy='1,1').run()

# Convert MODEL FITS file to Miriad format
shutil.rmtree(model,ignore_errors=True)
TaskFits(op='xyin',in_=model+'.FITS',out=model).run()

# Convert DISK FITS file to Miriad format
shutil.rmtree(disk,ignore_errors=True)
TaskFits(op='xyin',in_=disk+'.FITS',out=disk).run()

# add in the model - original had max of .017, for TB=2.e6 cell=2.6018 should
# have 243.5 Jy/pixel at 5 GHz, so rescaled image
shutil.rmtree('uvgen.sun',ignore_errors=True)
TaskUVModel(vis='uvgen.tmp', model=model, add=True, line='channel,1,1',
      out='uvgen.sun').run()

# now map and clean
# subtract disk for clean to be compatible with MAXEN
shutil.rmtree('uvgen.sub',ignore_errors=True)
TaskUVModel(vis='uvgen.sun', model=disk, subtract=True, line='channel,1,1',
      out='uvgen.sub').run()

shutil.rmtree('sub.beam',ignore_errors=True)
shutil.rmtree('sub.map',ignore_errors=True)
TaskInvert(vis='uvgen.sub', map='sub.map', beam='sub.beam',line='channel,1', 
      sup='200,200', robust=robust, double=True, imsize='1024,1024', 
      cell='2.4,2.4').run() 

TaskCgDisp(device='1/xs',in_='sub.beam', type='p', region='abspixel,box(769,769,1280,1280)').run()
shutil.rmtree('sub.clean',ignore_errors=True)
TaskClean(map='sub.map', beam='sub.beam', niters=5000, gain=.05, phat=0.3,
          clip=0.8, out='sub.clean', speed=-2, minpatch=127).run()
shutil.rmtree('sub.cm',ignore_errors=True)
TaskRestore(map='sub.map', beam='sub.beam', model='sub.clean', 
            out='sub.cm').run()

# restore disk in the clean (Jy/pixel) plane
# Unfortunately, mirexec does not have the ability to get a value from
# a command, but we can directly call gethd and capture the output by
# piping STDOUT.
import subprocess
proc = subprocess.Popen(['gethd','in=sub.cm/bmaj'],stdout=subprocess.PIPE)
bmaj = float(proc.stdout.read())*206265
proc = subprocess.Popen(['gethd','in=sub.cm/bmin'],stdout=subprocess.PIPE)
bmin = float(proc.stdout.read())*206265
proc = subprocess.Popen(['gethd','in=sub.cm/bpa'],stdout=subprocess.PIPE)
bpa = float(proc.stdout.read())
# convol does not seem to be trustworthy - calculate factor myself
sfac = 1.333*bmaj*bmin/2.4/2.4
print "Convolution rescaling is ",sfac

shutil.rmtree('sun.clean',ignore_errors=True)
TaskMaths(exp='(sub.clean+'+disk+')', out='sun.clean').run()
shutil.rmtree('sun.map',ignore_errors=True)
shutil.rmtree('sun.beam',ignore_errors=True)
TaskInvert(vis='uvgen.sun', map='sun.map', beam='sun.beam', line='channel,1',
       sup='200,200', robust=robust, double=True, imsize='1024,1024', 
       cell='2.4,2.4').run()
shutil.rmtree('sun.cm',ignore_errors=True)
TaskRestore(map='sun.map', beam='sun.beam', model='sun.clean', out='sun.cm',
       fwhm=str(bmaj)+','+str(bmin), pa=bpa).run()
TaskCgDisp(device='1/xs', in_='sun.cm').run()

# maximum entropy deconv: need to clean out bright stuff down to 50K first
# calculate cutoff for 50 K and this beam size
# turns out we want about 1000 as cutoff: 200K apparently
cutoff = 4*0.8993*bmaj*bmin
print "Cutoff is ",cutoff
shutil.rmtree('tmp.clean',ignore_errors=True)
# use earlier disk-subtracted maps, we don't care
TaskClean(map='sub.map', beam='sub.beam', niters=5000, gain=.05, phat=0.3, 
          clip=0.8, out='tmp.clean', speed=-2, minpatch=127, 
          cutoff=cutoff).run()
shutil.rmtree('uvgen.sub',ignore_errors=True)
TaskUVModel(vis='uvgen.sun', model='tmp.clean', subtract=True,
            line='channel,1,1', out='uvgen.sub').run()

shutil.rmtree('sub.map',ignore_errors=True)
shutil.rmtree('sub.beam',ignore_errors=True)
shutil.rmtree('sub.maxen',ignore_errors=True)
TaskInvert(vis='uvgen.sub', map='sub.map', beam='sub.beam', 
          line='channel,1', sup='200,200', robust=robust, double=True, 
          imsize='1024,1024', cell='2.4,2.4').run()
# No TaskMaxen, so we have to "roll our own"
# Note that only the used keywords are defined.  If you want to use others,
# add them to the _keywords list.
class TaskMaxen(TaskBase):
    _keywords = ['map','beam','out','niters','measure','rms','tol','default']

TaskMaxen(map='sub.map', beam='sub.beam', out='sub.maxen', niters=400,
         measure='gull', rms=0.6, tol=0.005, default=disk).run()
# combine models assuming they match original map+beam
shutil.rmtree('sun.maxen',ignore_errors=True)
TaskMaths(exp='(sub.maxen+tmp.clean)', out='sun.maxen').run()
# restore using map/beam with no subtraction from earlier
shutil.rmtree('sun.maxen_cnvl',ignore_errors=True)
TaskRestore(map='sun.map', beam='sun.beam', model='sun.maxen',
            out='sun.maxen_cnvl', fwhm=str(bmaj)+','+str(bmin), pa=bpa).run()
TaskCgDisp(device='2/xs', in_='sun.maxen_cnvl').run()
# No TaskConvol, so we have to "roll our own"
# Note that only the used keywords are defined.  If you want to use others,
# add them to the _keywords list.
class TaskConvol(TaskBase):
    _keywords = ['map','out','fwhm','pa','scale']

shutil.rmtree('sun.model',ignore_errors=True)
TaskConvol(map=model, out='sun.model', fwhm=str(bmaj)+','+str(bmin), 
           pa=bpa).run()
       # pa=bpa, scale=sfac

shutil.rmtree(out+'.uv',ignore_errors=True)
shutil.move('uvgen.sun',out+'.uv')
TaskFits(op='xyout', in_='sun.beam', out=out+'.BEAM.fits').run()
TaskFits(op='xyout', in_='sun.model', out=out+'.MODEL.fits').run()
TaskFits(op='xyout', in_='sun.maxen_cnvl', out=out+'.MAXEN.fits').run()
TaskFits(op='xyout', in_='sun.cm', out=out+'.CLEAN.fits').run()
shutil.rmtree('sun.beam',ignore_errors=True)
shutil.rmtree('sun.clean',ignore_errors=True)
shutil.rmtree('sun.cm',ignore_errors=True)
shutil.rmtree('sun.map',ignore_errors=True)
shutil.rmtree('sun.maxen',ignore_errors=True)
shutil.rmtree('sun.maxen_cnvl',ignore_errors=True)
shutil.rmtree('sun.model',ignore_errors=True)
shutil.rmtree('sub.beam',ignore_errors=True)
shutil.rmtree('sub.clean',ignore_errors=True)
shutil.rmtree('sub.cm',ignore_errors=True)
shutil.rmtree('sub.map',ignore_errors=True)
shutil.rmtree('sub.maxen',ignore_errors=True)

