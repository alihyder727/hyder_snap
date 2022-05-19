#! /usr/bin/env python3
from pylab import *

bands = ['s','c','x','u','k','q']
nbands = len(bands)

data = genfromtxt('../dat/vla_saturn_u.txt')
lat1 = data[:,0]
lat2 = data[:,1]
tb_avg = data[:,2]
tb_std = data[:,3]
mu = data[:,4]

fname = '../dat/vla_saturn_csxukq.txt'
with open(fname, 'w') as file:
  file.write('# VLA smeared map of saturn, 0.2 degree resolution, 2 - 48.5 GHz\n')
  file.write('# S Band, 13.33 - 15.04 cm, 2 - 2.25 GHz\n')
  file.write('# C Band, 7.08 - 7.53 cm, 4 - 4.25 GHz\n')
  file.write('# X Band, 2.5 - 3.76 cm, 8 - 12 GHz\n')
  file.write('# U Band, 1.66 - 2.5 cm, 12 - 18 GHz\n')
  file.write('# K Band, 1.15 - 1.67 cm, 18 - 26 GHz\n')
  file.write('# Q Band, 0.62 - 0.75 cm, 40 - 48.4 GHz\n')
  file.write('# Brightness temperature (K)\n')
  print('Writing brightness temperature ...')
  file.write('#%5s%10s%10s%10s%10s%10s%10s%10s%10s\n' % ('IDX', 'LATG1', 'LATG2', 'S', 'C', 'X',
      'U', 'K', 'Q'))
  tb_bands = zeros(nbands)
  for i in range(len(lat1)):
    for j in range(nbands):
      data = genfromtxt('../dat/vla_saturn_%s.txt' % bands[j])
      i1 = searchsorted(data[:,0], lat1[i], side = 'right')-1
      i2 = searchsorted(data[:,0], lat2[i], side = 'left')
      tb = sum(data[i1:i2,2]*data[i1:i2,6])
      num = sum(data[i1:i2,6])
      tb_bands[j] = tb/num
    file.write('%6d%10.1f%10.1f' % (i, lat1[i], lat2[i]))
    for j in range(nbands):
      file.write('%10.1f' % tb_bands[j])
    file.write('\n')

  stb_bands = zeros(nbands)
  file.write('# Brightness temperature standard deviation (K)\n')
  print('Writing brightness temperature standard deviation ...')
  file.write('#%5s%10s%10s%10s%10s%10s%10s%10s%10s\n' % ('IDX', 'LATG1', 'LATG2', 'S', 'C', 'X',
      'U', 'K', 'Q'))
  for i in range(len(lat1)):
    for j in range(nbands):
      data = genfromtxt('../dat/vla_saturn_%s.txt' % bands[j])
      i1 = searchsorted(data[:,0], lat1[i], side = 'right')-1
      i2 = searchsorted(data[:,0], lat2[i], side = 'left')
      stb = sum(data[i1:i2,3]*data[i1:i2,3]*data[i1:i2,6])
      num = sum(data[i1:i2,6])
      stb_bands[j] = sqrt(stb/num)
    file.write('%6d%10.2f%10.2f' % (i, lat1[i], lat2[i]))
    for j in range(nbands):
      file.write('%10.2f' % stb_bands[j])
    file.write('\n')

  mu_bands = zeros(nbands)
  file.write('# cosine emission angle (mu = cos theta)\n')
  print('Writing cosine emission angle ...')
  file.write('#%5s%10s%10s%10s%10s%10s%10s%10s%10s\n' % ('IDX', 'LATG1', 'LATG2', 'S', 'C', 'X',
      'U', 'K', 'Q'))
  for i in range(len(lat1)):
    for j in range(nbands):
      data = genfromtxt('../dat/vla_saturn_%s.txt' % bands[j])
      i1 = searchsorted(data[:,0], lat1[i], side = 'right')-1
      i2 = searchsorted(data[:,0], lat2[i], side = 'left')
      mu = sum(data[i1:i2,4]*data[i1:i2,6])
      num = sum(data[i1:i2,6])
      mu_bands[j] = mu/num
    file.write('%6d%10.2f%10.2f' % (i, lat1[i], lat2[i]))
    for j in range(nbands):
      file.write('%10.2f' % mu_bands[j])
    file.write('\n')

  smu_bands = zeros(nbands)
  file.write('# cosine emission angle standard deviation (mu sigma)\n')
  print('Writing cosine emission angle standard deviation ...')
  file.write('#%5s%10s%10s%10s%10s%10s%10s%10s%10s\n' % ('IDX', 'LATG1', 'LATG2', 'S', 'C', 'X',
      'U', 'K', 'Q'))
  for i in range(len(lat1)):
    for j in range(nbands):
      data = genfromtxt('../dat/vla_saturn_%s.txt' % bands[j])
      i1 = searchsorted(data[:,0], lat1[i], side = 'right')-1
      i2 = searchsorted(data[:,0], lat2[i], side = 'left')
      smu = sum(data[i1:i2,5]*data[i1:i2,5]*data[i1:i2,6])
      num = sum(data[i1:i2,6])
      smu_bands[j] = sqrt(smu/num)
    file.write('%6d%10.3f%10.3f' % (i, lat1[i], lat2[i]))
    for j in range(nbands):
      file.write('%10.3f' % smu_bands[j])
    file.write('\n')

print('output file written to %s' % fname)
