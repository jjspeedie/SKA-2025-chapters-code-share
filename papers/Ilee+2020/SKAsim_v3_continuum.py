
# Simple imaging simulation script

# Need to add 
#
# measures.observatory.directory: (/path/to/this/folder)/data
#
# to your ~/.casarc file


import numpy as np
import math
import time
import sys

try:
    MOD_NAME = sys.argv[3]
    NEIT = float(sys.argv[4])
    BMAJ = sys.argv[5]    
except:
    print("Please specify a valid input model name (Usage: > casa -c SKAsim_v3_continuum.py [modname.fits] [NEIT] [BMAJ]")
    sys.exit(1)


DO_SIMULATE = True                  # Do a new simulated MS (rather than using existing one)
OUT_ROOT = 'sim_'+MOD_NAME[:-5]
DUMP_TIME = 0.14                    # integration time in seconds
SCAN_TIME = (DUMP_TIME*1.)/3600.    # make each scan last a few integrations, units of hours
NUM_SCANS = 31                      # accumulate this many scans (use odd number)
HA_RANGE = 8.0                      # total range of HA coverage in hours (use > 0.0)

FREQ = 12.5                         # centre frequency in GHz
FREQ_FBW = 0.3                      # frequency fractional bandwidth per channel
NUM_CHAN = 1                        # number of spectral channels
NUM_SPW = 1                         # number of spectral windows (use odd number)
FREQ_SPW = 0.01                     # offset between SPWs in GHz

NEBW = 5.0                          # Noise effective total bandwidth in GHz
#NEIT = 1000.                        # Noise effective total integration time in hours
NPOL = 2.                           # Noise effective number of polarisations

RA = '16h26m30.0s'
DEC = '-24d24m00.0s'

TSKY_FILE = 'data/SKA_Tsky.txt';         # File with columns freq, tsky10, tsky50, tsky90, by default uses tsky10
TAU_FILE = 'data/SKA_tau.txt';           # File with columns freq, tau10, tau50, tau90, by default uses tau10
LOW_FILE = 'data/SKALA4_Specs.txt';      # File with columns freq in GHz, Low station A/T in m2/K

#MOD_NAME = 'continuum_12GHz.fits'   # model sky image
#MOD_NAME = 'new_RT.fits'

MOD_SCALE = 10.0                    # Gives F(2.5cm) ~ 2mJy

NPIX = 1000                        # size of output image (2x padding is used for gridding and FFT)
#CELL = '0.0015arcsec'               # pixel size
CELL = '0.00714286arcsec'              # pixel size
DO_NAT = False                      # Make use of NATURAL weighting 
DO_UNIF = True                      # Make use of UNIFORM weighting followed by Gauusian tapering
FOV = '10arcsec'                    # FoV for sidelobe minimisation during UNIFORM weighting
#BMAJ = '0.034arcsec'                  # Gaussian taper after UNIFORM weighting
DO_ROB = False                      # Make use of ROBUST weighting 
ROBUST = 0.                         # Briggs robustness parameter


# ----------------------------------------------------------------------------

FREQ_RES = str(FREQ*FREQ_FBW*1000.) +'MHz'                  # channel freq resolution

fracbw = FREQ_FBW*float(NUM_CHAN)
#output_file = OUT_ROOT + '_freq_' + str(FREQ) + '_bm_' + str(BMAJ) + str(fracbw) + 'fbw' + str(MOD_SCALE) + 'modscale'

output_file = OUT_ROOT + '_freq_' + str(FREQ) + '_bm_' + str(BMAJ) + str(NEIT) + 'neit' + str(MOD_SCALE) + 'modscale'

xfreq = str(FREQ) + 'GHz'
output_ms = output_file + '.MS'
if (DO_UNIF):
  output_im = output_file + '.fw' + FOV
if (DO_NAT):
  output_im = output_file + '.na'
if (DO_ROB):
  output_im = output_file + '.ro' + str(ROBUST)
fieldDir = me.direction( 'J2000', RA, DEC )
REF_TIME = me.epoch('UTC','today');                        

if (DO_SIMULATE):

  os.system( 'rm -rf skymodel \n' )
  os.system( 'rm -rf temp \n' )
  os.system( 'rm -rf ' + output_ms + '\n' )
  os.system( 'rm -rf ' + output_file + '-vp.tab' + '\n' )
  sm.open( output_ms );

  vp.reset()
  if (FREQ > 0.05 and FREQ < 0.35) :
    TEL_NAME = 'SKA1-LOW';
    TEL_FILE = 'data/ska1low.cfg';            
    vp.setpbairy( telescope = TEL_NAME, dopb = True, dishdiam = '38.0m', blockagediam = '0.0m', maxrad = '12.0deg', reffreq = '300.0MHz')
    vp.saveastable( tablename = output_file + '-vp.tab' )
  if (FREQ > 0.58 and FREQ < 3.05) :  
    TEL_NAME = 'SKA1-MID';
    TEL_FILE = 'data/ska1mid.cfg';            
    vp.setpbairy( telescope = TEL_NAME, dopb = True, dishdiam = '15.0m', blockagediam = '0.0m', maxrad = '12.0deg', reffreq = '300.0MHz')
    vp.saveastable( tablename = output_file + '-vp.tab' )
  if ((FREQ > 0.35 and FREQ < 0.58) or (FREQ > 3.05 and FREQ < 50)) :
    TEL_NAME = 'SKA1-MID';
    TEL_FILE = 'data/ska1mid-m.cfg';			# 133 dishes
    #TEL_FILE = 'data/ska1mid_64skaonly.cfg';
    vp.setpbairy( telescope = TEL_NAME, dopb = True, dishdiam = '15.0m', blockagediam = '0.0m', maxrad = '12.0deg', reffreq = '300.0MHz')
    vp.saveastable( tablename = output_file + '-vp.tab' )
  # read configuration file
  f = open( TEL_FILE, 'r' )
  xx = []; yy = []; zz = []; dd = []; ss = [];
  while True:
    line = f.readline()
    if not line: break
    items = line.split()
    if ( items[0] != '#' ):
      xx.append( float( items[0] ) )
      yy.append( float( items[1] ) )
      zz.append( float( items[2] ) )
      dd.append( float( items[3] ) )
      ss.append( str( items[4] ) )
  f.close()
  # make a primary beam pattern:
  # hpbw = str( FOV / 60. ) + 'deg'
  # vp.setpbgauss( telescope = TEL_NAME, dopb = True, halfwidth = hpbw, reffreq = xfreq)
    # need appropriate entry in Observatories table
  posska = me.observatory( TEL_NAME);
  #print posska
  sm.setconfig(telescopename = TEL_NAME, x = xx, y = yy, z = zz, dishdiameter = dd, antname = ss, mount = 'ALT-AZ', 
             coordsystem = 'global', referencelocation = posska )
  jvar = NUM_SPW +1
  for j in range(1, jvar):
    # define the spectral windows
    ifnum = (j - 1)
    xspwname = 'IF' + str(ifnum)
    spwfreq = str(FREQ + float( j - ((NUM_SPW + 1) / 2) )*FREQ_SPW - float((NUM_CHAN / 2))*FREQ*FREQ_FBW) + 'GHz'
    sm.setspwindow( spwname = xspwname, freq = spwfreq, deltafreq = FREQ_RES, freqresolution = FREQ_RES, nchannels = NUM_CHAN, stokes = 'XX' ) 
    print (xspwname, spwfreq)
  sm.setfeed( mode = 'perfect X Y', pol = [''] )
  sm.setfield( sourcename = 'ptcentre', sourcedirection = fieldDir )
  sm.setlimits( shadowlimit = 0.001, elevationlimit = '15.0deg' )
  sm.setauto( autocorrwt = 0.0 )
  sm.settimes( integrationtime = str( DUMP_TIME ) + 's', usehourangle = True, referencetime = REF_TIME )
  ivar = NUM_SCANS +1
  for i in range(1, ivar):
    # set the start and end times.
    if (NUM_SCANS > 1): 
      begin = float( i - ((NUM_SCANS + 1) / 2) )/float(NUM_SCANS-1)*HA_RANGE  -  SCAN_TIME / 2.
    else:
      begin =  - SCAN_TIME / 2.    
    end = begin +  SCAN_TIME 
    # convert to string.
    begin = str( begin ) + 'h'
    end = str( end ) + 'h'
    for j in range(1, jvar):
      ifnum = (j - 1)
      xspwname = 'IF' + str(ifnum)
      sm.observe( sourcename= 'ptcentre', spwname= xspwname, starttime = begin, stoptime = end )      
      
  # get model image from fits file
  os.system('rm -rf temp')
  os.system('rm -rf skymodel')
  ia.fromfits( infile = MOD_NAME )
  print ia.shape()
  cs = ia.coordsys()
  print cs.axiscoordinatetypes()  
  ia.adddegaxes( outfile = 'temp', spectral = False, stokes = 'I')   
  ia.close()
#  im3 = ia.subimage(outfile = 'temp', overwrite = True)  
#  im3.done()
#  ia.close()                                                      
  ia.open('temp')
  pixstring = 'temp * ' + str(MOD_SCALE)
  myim=ia.imagecalc( outfile = 'skymodel', pixels=pixstring, overwrite = True )
  myim.done()
  ia.close()
  ia.open('skymodel')
  cs = ia.coordsys()
  cs.setdirection(refval=RA+' '+DEC,refpix=[float(NPIX)/2,float(NPIX)/2])
  cs.setreferencevalue( xfreq, 'spectral' )
  cs.setincrement( FREQ_RES, 'spectral' )
  print cs.axiscoordinatetypes()  
  ia.setcoordsys(cs.torecord())
  sm.setvp( dovp = True, usedefaultvp = False, vptable = output_file + '-vp.tab', dosquint = False )
  # generate the predicted visibilities.
  sm.predict( imagename = 'skymodel', incremental = True)
  #ia.close()


  # work out the noise model appropriate for SKA1-Low

  # read data file
  f = open( LOW_FILE, 'r' )
  frqz = []; aontz = []; 
  while True:
    line = f.readline()
    if not line: break
    items = line.split()
    if ( items[0] != '#' ):
      frqz.append( float( items[0] ) )
      aontz.append( float( items[1] ) )
  f.close()

  skad = 1.
  skan = 0.

  if (FREQ > 0.05 and FREQ < 0.35) :  
    skad = 38.
    skan = 512.
    aont_tot = 0.
    for i in range(1, NUM_CHAN):
      ifreq = FREQ - float((NUM_CHAN/2 - i + 1))*FREQ*FREQ_FBW 
      iaontlow = skan*np.interp(ifreq,frqz,aontz)
      aont_tot = aont_tot + iaontlow
    aontlow = aont_tot/NUM_CHAN
  else:
    aontlow = 0.

  # work out the noise model appropriate for SKA1-Mid

  # read tsky file
  f = open( TSKY_FILE, 'r' )
  frqx = []; tsky10x = []; tsky50x = []; tsky90x = []; 
  while True:
    line = f.readline()
    if not line: break
    items = line.split()
    if ( items[0] != '#' ):
      frqx.append( float( items[0] ) )
      tsky10x.append( float( items[1] ) )
      tsky50x.append( float( items[2] ) )
      tsky90x.append( float( items[3] ) )
  f.close()
  # read tau file
  f = open( TAU_FILE, 'r' )
  frqy = []; tau10y = []; tau50y = []; tau90y = []; 
  while True:
    line = f.readline()
    if not line: break
    items = line.split()
    if ( items[0] != '#' ):
      frqy.append( float( items[0] ) )
      tau10y.append( float( items[1] ) )
      tau50y.append( float( items[2] ) )
      tau90y.append( float( items[3] ) )
  f.close()

  tsky10_tot = 0.
  tau10_tot = 0.
  if (NUM_CHAN > 1):
    for i in range(1, NUM_CHAN):
      ifreq = FREQ - float((NUM_CHAN/2 - i + 1))*FREQ*FREQ_FBW 
      itsky10 = np.interp(ifreq,frqx,tsky10x)
      itau10 = np.interp(ifreq,frqy,tau10y)
      tsky10_tot = tsky10_tot + itsky10
      tau10_tot = tau10_tot + itau10
    tsky10 = tsky10_tot/NUM_CHAN
    tau10 = tau10_tot/NUM_CHAN
  else:
    tsky10 = np.interp(FREQ,frqx,tsky10x)
    tau10 = np.interp(FREQ,frqy,tau10y)

  if (FREQ > 0.35 and FREQ < 50) :
    skad = 15.
    skan = 133.
  if (FREQ > 0.35 and FREQ < 0.95) :
    trcv = 15.+30.*(FREQ-0.75)**2
  if (FREQ >= 0.95 and FREQ <= 1.76) :
    trcv = 7.5
  if (FREQ > 1.76 and FREQ < 3.05) :
    trcv = 7.5
  if (FREQ > 3.05 and FREQ < 5.18) :
    trcv = 7.5
  if (FREQ > 5.18 and FREQ < 50) :
    trcv = 4.4+0.69*FREQ
  etaa0 = 0.92 - 0.04*abs(log10(FREQ))
  lamb = 3.e8/(FREQ*1.e9)
  etad = 1.-20.*(lamb/skad)**1.5
  errp = 280.e-6
  errs = 154.e-6
  a_p = 0.89
  a_s = 0.98
  Del = 2.*(a_p*errp**2+a_s*errs**2)**0.5
  delph = 2.*pi*Del/lamb
  etaph = exp(-delph**2)
  etaa = etaa0 * etad * etaph
  if (FREQ <= 0.35 or FREQ >= 50) :
    etaa = 1.e-20
  tspl = 3.

  if (FREQ > 0.58 and FREQ < 3.05) :
    mkd = 13.5
    mkn = 64.
  else :
    mkn = 0.
    mkd = 1.e20
  mtrcv = 1.e20
  if (FREQ > 0.58 and FREQ < 0.9) :
    mtrcv = 11.-4.5*(FREQ-0.58)
  if (FREQ > 0.9 and FREQ < 1.67) :
    mtrcv = 7.5+6.8*(abs(FREQ-1.65))**1.5
  if (FREQ > 1.67 and FREQ < 3.05) :
    mtrcv = 7.5
  metaa0 = 0.8 - 0.04*abs(log10(FREQ))
  metad = 1.-20.*(lamb/mkd)**1.5
  merrp = 480.e-6
  merrs = 265.e-6
  ma_p = 0.89
  ma_s = 0.98
  mDel = 2.*(ma_p*merrp**2+ma_s*merrs**2)**0.5
  mdelph = 2.*pi*mDel/lamb
  metaph = exp(-mdelph**2)
  metaa = metaa0 * metad * metaph
  if (FREQ <= 0.58 or FREQ >= 3.05) :
    metaa = 1.e-20
  mtspl = 4.

  kb = 1.38e-16/1.e4/1.e-23
  pi = 3.1415926535
  hpm = 6.63e-34
  kbm = 1.38e-23
  fnu = FREQ*1.e9

  tsysp = trcv+tsky10+tspl
  if (FREQ > 0.35 and FREQ < 50) :  
    tsysx = (hpm*fnu/kbm)/(exp(hpm*fnu/(kbm*tsysp))-1.)
  else:
    tsysx = 1.e20 
  tsys = tsysx/exp(-tau10)

  aeff = skan*etaa*pi*skad**2/4.
  aont = aeff/tsys
  
  mtsysp = mtrcv+tsky10+mtspl
  if (FREQ > 0.58 and FREQ < 3.05) :  
    mtsysx = (hpm*fnu/kbm)/(exp(hpm*fnu/(kbm*mtsysp))-1.)
  else:
    mtsysx = 1.e20 
    mkn = 0.
  mtsys = mtsysx/exp(-tau10)

  maeff = mkn*metaa*pi*mkd**2/4.
  maont = maeff/mtsys

  caont = aont+maont+aontlow
  csefd = 2.*kb/caont                  # full array SEFD in Jy
  nvis = (skan + mkn)*(skan+mkn-1)/2.  # total number of instantaneous correlated baselines
  vsefd = csefd*(nvis)**0.5            # average SEFD per visibility in Jy
  eta_q = 0.9                          # total system efficiency
  dnueff = NEBW*1.e9/float(NUM_SPW)/float(NUM_CHAN) # effective bandwidth for noise calculation
  dps = SCAN_TIME*3600./DUMP_TIME                   # dumps per scan
  dtaueff = NEIT*3600./float(NUM_SCANS)/dps         # effective integration time for noise calculation
  vrmsjy = vsefd/(eta_q*(NPOL*dnueff*dtaueff)**0.5) # effective visibility noise
  
  # Double the noise for the smallest beam
  if BMAJ == '0.034arcsec':
    sm.setnoise(mode='simplenoise', simplenoise=str(2*vrmsjy) + 'Jy')
    #print output_im, 'vrmsjy = ', 2*vrmsjy*1.0e3
  else: 
    sm.setnoise(mode='simplenoise', simplenoise=str(vrmsjy) + 'Jy')
    #print output_im, 'vrmsjy = ', vrmsjy*1.0e3

  
  sm.corrupt()

  sm.close()
  sm.done()
  ia.close()

spwmax = NUM_SPW
# prepare input measurement set.
#fixvis(vis=output_ms, outputvis=output_ms+'.rightdir', field='', phasecenter='J2000 00h00m00.s -30d00m00.0s')
#output_ms = output_ms+'.rightdir'
im.open( output_ms )
#im.fixvis(phasedirs='J2000 00h00m00.s -30d00m00.0s')
#im.done()
#im.open( output_ms )
im.defineimage( nx = NPIX, ny = NPIX, cellx = CELL, celly = CELL, stokes = 'XX', spw = range(0,spwmax), phasecenter = fieldDir, mode = 'mfs' )
im.setoptions( padding = 2.0 )
im.selectvis( spw = range(0,spwmax) )
if (DO_NAT):
  im.weight( type = 'natural')
if (DO_UNIF):
  im.weight( type = 'uniform', fieldofview = FOV)
  im.filter( type = 'gaussian', bmaj = BMAJ, bmin = BMAJ) 
if (DO_ROB):
  im.weight( type = 'briggs', rmode = 'norm', robust = ROBUST)
# overwrite existing images:
os.system('rm -rf ' + output_im + '.map')
os.system('rm -rf ' + output_im + '.map.fits')
os.system('rm -rf ' + output_im + '.psf')
os.system('rm -rf ' + output_im + '.psf.fits')




im.setvp( dovp = True, usedefaultvp = False, vptable = output_file + '-vp.tab', dosquint = False )
im.makeimage( type = 'observed', image=output_im + '.map' )



im.makeimage( type = 'psf', image=output_im + '.psf' )
params = im.fitpsf( output_im + '.psf' )
im.makeimage( type = 'pb', image=output_im + '.pb' )

im.makeimage( type = 'model', image=output_im + '.model' )
im.makeimage( type = 'corrected', image=output_im + '.corrected' )

#im.setbeam(bmaj='1.18arcsec', bmin='1.18arcsec', bpa='76.6deg')
#im.restore(model=output_im + '.map' , complist='', image='toto.restored' , residual='toto.residual' )

im.done()

ia.open(output_im + '.map')
stats = ia.statistics()
print stats["rms"]
ia.setrestoringbeam(major=BMAJ,minor=BMAJ,pa='0deg')
ia.close()


exportfits(imagename = output_im + '.map', fitsimage=output_im + '.map.fits')
#exportfits(imagename = output_im + '.psf', fitsimage=output_im + '.psf.fits')
#exportfits(imagename = output_im + '.pb', fitsimage=output_im + '.pb.fits')
#exportfits(imagename = output_im + '.corrected', fitsimage=output_im + '.corrected.fits')


rmDir = "rm -rf `ls -1 -d "+output_file+"*/`"
#os.system(rmDir)
rmDir = "rm -rf skymodel"
os.system(rmDir)
rmDir = "rm -rf temp"
os.system(rmDir)

