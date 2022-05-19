#! /usr/bin/env python3
from numpy import sort, zeros, nan

def is_number(s):
  try:
    float(s)
    return True
  except ValueError:
    return False

def read_info(fname):
  info = {}
  with open(fname, 'r') as file:
    for line in file.readlines():
      if line[0] == '#':
        continue
      field, value = line.split('=')
      field = field.strip(' \n')
      value = value.strip(' \n')
      #print(field)
      if 'FITS_FILE' == field:
        info['fits'] = value
      if 'PLANET' == field:
        info['target'] = value
      if 'TIME_START' in field:
        info['date'] = value.split(',')[0]
      if 'OBLON' in field:
        info['oblon'] = float(value)
      if 'OBLAT' in field:
        info['oblat'] = float(value)
      if 'PIXSCALE1' in field:
        info['scale1'] = float(value)
      if 'PIXSCALE2' in field:
        info['scale2'] = float(value)
      if 'PIXELS' in field:
        info['xsize'] = float(value.split('x')[0])
        info['ysize'] = float(value.split('x')[1])
      if 'CENTER' in field:
        info['xcenter'] = float(value.split(',')[0])
        info['ycenter'] = float(value.split(',')[1])
      if 'RE' in field:
        info['re'] = float(value)
      if 'RP' in field:
        info['rp'] = float(value)
  return info

def regrid(lon, lat, tb, mu):
  longrid = sort(list(set(lon)))
  latgrid = sort(list(set(lat)))

  tbgrid = zeros((len(longrid), len(latgrid)))
  mugrid = zeros((len(longrid), len(latgrid)))
  for i,v in enumerate(tb):
    ilon = list(longrid).index(lon[i])
    ilat = list(latgrid).index(lat[i])
    tbgrid[ilon,ilat] = v
    mugrid[ilon,ilat] = mu[i]
  tbgrid[tbgrid == 0.] = nan
  mugrid[mugrid == 0.] = nan

  return longrid, latgrid, tbgrid, mugrid
