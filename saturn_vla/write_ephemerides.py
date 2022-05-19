#! /usr/bin/env python3
from datetime import datetime, timedelta
from urllib.request import urlopen, urlretrieve
from numpy import array, radians
from common import is_number 
from astropy.io import fits
import sys, os, argparse, re

# adapted from nirc2_reduce/get_ephem.py
def naif_lookup(target):
  '''looks up NAIF ID for target in list at naif_id_table.txt'''
  target = target.upper().strip(', \n')
  with open('../naif_id_table.txt','r') as f:
    for line in f:
      l = line.split(',')
      if is_number(target):
        #print(l[0].strip(', \'\n') == target)
        if l[0].strip(', \'\n') == target:
          code = target
      else:
        if l[1].strip(', \'\n') == target:
          code = l[0].strip(', \n')
    try:
      code
    except:
      msg = 'WARN (naif_lookup): NAIF code not in lookup table.'
      msg += 'If code fails, ensure target can be queried in Horizons.'
      print(msg)
      code = target
                          
  if len(code) == 7: #minor body
    if code[0] == '2': #asteroid
      return code[1:]+';'
    elif code[0] =='1': #comet
      msg = 'ERROR (naif_lookup): Comets cannot be looked up by NAIF ID; '
      msg += 'Horizons generates multiple codes for every comet. '
      msg += 'Try your search in the Horizons Web tool, '
      msg += 'select the target body you want, '
      msg += 'and then copy the exact string into this code and it *may* work.'
      sys.exit(msg)
    return code
  return code

# adapted from nirc2_reduce/get_ephem.py
def get_ephemerides(code, obs_code, tstart, tend, stepsize, quantities = None) :
  """
  input NAIF target code, e.g. 501 for Io, and date in the format:
  'YYYY-MM-DD HH:MM'
  For example: data=get_ephem.get_ephemerides('501','2017-06-09 08:24')
  Returns a list containing (in string format):
  UTdate and time,sun,moon,RA (J2000),DEC (J2000),dra,ddec,azimuth,elevation,
  Airmass,Extinction,APmag,s-brt,Ang-Diam("),ang-sep("),visibility,Ob-lon,Ob-lat,NP.ang,NP.dist
  """
  if quantities == None:
    quantities = '1,3,4,8,9,12,13,14,15,17,19,20'

  tstart_obj = datetime.strptime(tstart,'%Y-%m-%d %H:%M')
  tend_obj = datetime.strptime(tend,'%Y-%m-%d %H:%M')
  if tend_obj - tstart_obj <= timedelta(minutes = 1):
    print('End time before start time. Setting end time to one minute after start time.')
    tend_obj = tstart_obj + timedelta(minutes = 1)
  
  tstart_UT = datetime.strftime(tstart_obj,"'%Y-%m-%d %H:%M'")
  tend_UT = datetime.strftime(tend_obj,"'%Y-%m-%d %H:%M'")

  http = "https://ssd.jpl.nasa.gov/horizons_batch.cgi?batch=1"
  make_ephem = "&MAKE_EPHEM='YES'&TABLE_TYPE='OBSERVER'"
  command    = "&COMMAND=" + str(code)
  center     = "&CENTER="+str(obs_code)  #568 is Mauna Kea, 662 is Lick, etc
  t_start    = "&START_TIME=" + tstart_UT
  t_stop     = "&STOP_TIME=" + tend_UT
  t_step     = "&STEP_SIZE='" + stepsize + "'"
  quantities = "&QUANTITIES='%s'"%quantities
  csv        = "&CSV_FORMAT='YES'"

  url = http+make_ephem+command+center+t_start+t_stop+t_step+quantities+csv
  url = url.replace(" ", "%20")
  url = url.replace("'", "%27")
  try:
    with urlopen(url) as response:
      ephem = response.readlines()  
  except:
    msg = 'ERROR (get_ephem): Could not retrieve query from Horizons. '
    msg += 'Check Internet connection and input URL'
    sys.exit(msg)

  inephem = False
  data = []
  for i in range(len(ephem)) :
    if type(ephem[i]) != str:
      ephem[i] = ephem[i].decode('UTF-8')
    if inephem == False:
      #get observatory lat, lon, alt for later
      if ephem[i].startswith('Center geodetic :'):
        l = ephem[i].split(':')[1]
        latlonalt = l.split()[0]
        [lon,lat,alt] = [float(s.strip(', \n')) for s in latlonalt.split(',')]
        observatory_coords = [lat,lon,alt]
      if ephem[i].startswith('$$SOE') :
        inephem=True
        #data = [ephem[i+1].decode('UTF-8').split(',')]
    elif inephem == True:
      if ephem[i].startswith("$$EOE") :
        inephem=False
      else:
        data.append(ephem[i].split(','))
  if len(data) == 0:
    msg = 'ERROR (get_ephem): Ephemeris data not found. '
    msg += 'Check that the target has valid ephemeris data for the specified time range.'
    raise ValueError(msg)
  else:
    out = array(data)[:,:-1]
  return out, observatory_coords

if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument('-d', '--dir',
    default = '.',
    help = 'directory of the fits file'
    )
  parser.add_argument('-f', '--fits',
    default = 'tbd',
    help = 'base name of the fits file without extension'
    )
  parser.add_argument('-b', '--obs',
    default = 'tbd',
    help = 'observation facility'
    )
  parser.add_argument('-s', '--shift',
    default = '0',
    help = 'shift oblat by'
    )
  args = vars(parser.parse_args())

  print('Running write_ephemerides ...')

  hdu = fits.open('%s/%s.fits' % (args['dir'], args['fits']))[0]
  date = hdu.header['DATE-OBS']
  planet = hdu.header['OBJECT']
  del1 = hdu.header['CDELT1']
  del2 = hdu.header['CDELT2']

  naif = naif_lookup(planet)

  if args['obs'] == 'vla':
    obscode = -5
  elif args['obs'] == 'alma':
    obscode = -7
  elif args['obs'] == 'keck':
    obscode = 568
  elif args['obs'] == 'lick':
    obscode = 662
  else:
    raise ValueError('Observation Code Not Found.')
  tstart = datetime.strptime(date, '%Y-%m-%dT%H:%M:%S.%f')
  tend = tstart + timedelta(minutes = 60)

  tstart = datetime.strftime(tstart, '%Y-%m-%d %H:%M')
  tend = datetime.strftime(tend, '%Y-%m-%d %H:%M')
  ephem = get_ephemerides(naif, obscode, tstart, tend, '1 minutes')[0][0]
  ephem = [val.strip(' ') for val in ephem]
  #print(ephem)

  time = ephem[0]
  ra, dec = ephem[3], ephem[4]
  dra, ddec = float(ephem[5]), float(ephem[6])
  ob_lon, ob_lat = float(ephem[16]), float(ephem[17])
  print(ob_lat)
  ob_lat += float(args['shift'])
  print(ob_lat)
  np_ang = float(ephem[20])
  dist = float(ephem[24])*1.496e8 # au -> km
  pix1 = dist*abs(radians(del1))
  pix2 = dist*abs(radians(del2))

  os.system('mkdir -p ../dat')
  with open('../dat/%s.info' % args['fits'], 'w') as file:
    file.write('# EPHEMERIDES\n')
    file.write('FITS_FILE         = %s/%s.fits\n' % (args['dir'],args['fits']))
    file.write('PLANET            = %s\n' % planet)
    file.write('NAIF_ID_CODE      = %s\n' % naif)
    file.write('FACILITY          = %s\n' % args['obs'])
    file.write('FACILITY_CODE     = %d\n' % obscode)
    file.write('TIME_START        = %s\n' % tstart.replace(' ', ','))
    file.write('TIME_END          = %s\n' % tend.replace(' ', ','))
    file.write('RA                = %s\n' % ra.replace(' ', ','))
    file.write('DEC               = %s\n' % dec.replace(' ' , ','))
    file.write('DRA               = %s\n' % dra)
    file.write('DDEC              = %s\n' % ddec)
    file.write('OBLON (W)         = %s\n' % ob_lon)
    file.write('OBLAT (GRAPHIC)   = %s\n' % ob_lat)
    file.write('NPANG (DEG)       = %s\n' % np_ang)
    file.write('DIST (KM)         = %s\n' % dist)
    file.write('DEL1 (DEG)        = %s\n' % del1)
    file.write('DEL2 (DEG)        = %s\n' % del2)
    file.write('PIXSCALE1 (KM)    = %s\n' % pix1)
    file.write('PIXSCALE2 (KM)    = %s\n' % pix2)
