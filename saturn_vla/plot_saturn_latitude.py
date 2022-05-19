#! /usr/bin/env python3
from pylab import *

#bands = ['s', 'c', 'x', 'u', 'k', 'q']
bands = ['q','k','u','x','c','s']
nbands = len(bands)

# disk Tb
coeffs = genfromtxt('../dat/limb_disk.txt', delimiter = ',')

fig, axs = subplots(nbands, 1, figsize = (10, 10), sharex = True)
subplots_adjust(hspace = 0.12, wspace = 0.08)

for i in range(nbands):
  if bands[i] == 's':
    row = 0
  elif bands[i] == 'c':
    row = 1
  elif bands[i] == 'x':
    row = 2
  elif bands[i] == 'u':
    row = 4
  elif bands[i] == 'k':
    row = 6
  elif bands[i] == 'q':
    row = 7
  else :
    pass
  coeff = coeffs[row,1:-1]
  ax = axs[i]
  data = genfromtxt('../dat/vla_saturn_%s.txt' % bands[i])
  lat1 = data[:,0]
  lat2 = data[:,1]
  tb_avg = data[:,2]
  tb_std = data[:,3]
  mu = data[:,4]
  ax.errorbar((lat1+lat2)/2., tb_avg, yerr = tb_std,
    fmt = '.', color = '0.7', alpha = 0.5)
  ax.step(lat1, tb_avg, 'C1', where = 'post', label = '%s band' % bands[i].upper())
  tb_disk = (coeff[0]*pow(mu, coeff[1]) + coeff[5]*pow(mu, coeff[6]))/2.
  ax.step(lat1, tb_disk, 'C0', where = 'post', label = '%s band disk' % bands[i].upper())
  ax.set_xlim([0., 90.])
  ax.set_ylabel('T$_b$ (K)', fontsize = 10)
  if i == nbands - 1:
    ax.set_xlabel('Planetographic Latitude (degree)', fontsize = 10)
  ax.legend()

savefig('../figs/saturn_tb_latitude.png', bbox_inches = 'tight')
