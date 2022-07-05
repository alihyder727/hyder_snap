#! /usr/bin/env python3
from pylab import *
from glob import glob
from scipy.interpolate import interp1d
from h5py import File
import argparse
import matplotlib, re

parser = argparse.ArgumentParser()
parser.add_argument('--ch',
        choices = ['1','2','3','4','5','6'],
        default = '1',
        help='select mwr channel')
parser.add_argument('--lmin',
        default = '-60.',
        help='minimum latitude')
parser.add_argument('--lmax',
        default = '60.',
        help='maximum latitude')
args = vars(parser.parse_args())

latmin, latmax = float(args['lmin']), float(args['lmax'])
mwr_channel = int(args['ch'])
pjs = [1,3,4,5,6,7,8,9,12]
fig, axs = subplots(2, 1, figsize = (12, 8), sharex = True)

# read mwr original data
tb, ld = [], []
ff = File('./fao_ess_2020_data.h5', 'r')
for pj in pjs:
  data = ff['channel %d' % mwr_channel]['perijove %02d' % pj]
  coeff = data['coeff'][:]
  lat = data['latitude'][:]
  R45 = data['R45'][:]
  ix = where((lat > latmin) & (lat < latmax))[0]
  lat = lat[ix]
  tb.append(coeff[0,ix])
  ld.append(R45[ix])
tb = array(tb).T
ld = array(ld).T

# read gravity correction data
data = genfromtxt('data/gravity_correction_factor.txt')
rows = len(data)
data = data.reshape(6, rows//6, -1)
latx = data[0,:,0]
gtb = interp1d(latx, data[mwr_channel-1,:,1])
gld = interp1d(latx, data[mwr_channel-1,:,2])

# gravity correction
for i in range(len(lat)):
  tb[i,:] *= gtb(abs(lat[i]))
  ld[i,:] *= gld(abs(lat[i]))

# Tb filtering
eps_tb = zeros((len(lat), len(pjs)))
for i in range(len(lat)):
  a, x = list(tb[i]), list(range(len(pjs)))
  while len(a) > 1:
    eps_max = 0.
    for j,v in enumerate(x):
      eps_tb[i,v] = tb[i,v] - (sum(a) - tb[i,v])/(len(a) - 1)
      if abs(eps_tb[i,v]) >= eps_max:
        k = j
        eps_max = abs(eps_tb[i,v])
    a.pop(k)
    x.pop(k)

# Ld filtering
eps_ld = zeros((len(lat), len(pjs)))
for i in range(len(lat)):
  a, x = list(ld[i]), list(range(len(pjs)))
  while len(a) > 1:
    eps_max = 0.
    for j,v in enumerate(x):
      eps_ld[i,v] = ld[i,v] - (sum(a) - ld[i,v])/(len(a) - 1)
      if abs(eps_ld[i,v]) >= eps_max:
        k = j
        eps_max = abs(eps_ld[i,v])
    a.pop(k)
    x.pop(k)

# save filtered data
with open('data/ch%d_gfilter.txt' % mwr_channel, 'w') as file:
  file.write('# CH1 filtered data\n')
  file.write('# Tb\n')
  for i in range(len(lat)):
    file.write('%10.4f' % lat[i])
    for j in range(len(pjs)):
      file.write('%10.4f' % tb[i,j])
    file.write('\n')
  file.write('# eps Tb\n')
  for i in range(len(lat)):
    file.write('%10.4f' % lat[i])
    for j in range(len(pjs)):
      file.write('%10.4f' % eps_tb[i,j])
    file.write('\n')
  file.write('# ld\n')
  for i in range(len(lat)):
    file.write('%10.4f' % lat[i])
    for j in range(len(pjs)):
      file.write('%10.4f' % ld[i,j])
    file.write('\n')
  file.write('# eps ld\n')
  for i in range(len(lat)):
    file.write('%10.4f' % lat[i])
    for j in range(len(pjs)):
      file.write('%10.4f' % eps_ld[i,j])
    file.write('\n')
