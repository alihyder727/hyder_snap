#! /usr/bin/env python3
from pylab import *
from netCDF4 import Dataset

clat_grid = genfromtxt('data/jup_latitude_grid.txt', dtype = str)
grav_grid = genfromtxt('data/jup_gravity_grid.txt', dtype = str)

lats = zeros(len(clat_grid))
tb0 = zeros((len(clat_grid), 6))
ld = zeros((len(clat_grid), 6))
for i,clat in enumerate(clat_grid):
  lats[i] = float(clat)
  data = Dataset('data/juno_mwr-modelA-%s.nc' % clat, 'r')
  tb = data['radiance'][0,:,0,0]
  tb0[i,:] = tb[::4]
  tb45 = tb[3::4]
  ld[i,:] = (tb0[i,:] - tb45)/tb0[i,:]*100.

fig, axs = subplots(6, 2, figsize = (12, 8), sharex = True)
for i in range(6):
  # tb
  ax = axs[5-i,0]
  ax.plot(lats, tb0[:,i], 'k')
  ax.plot(-lats, tb0[:,i], 'k')
  ax.set_ylabel('T$_b$ (K)')
  if i == 0:
    ax.set_xlabel('Planetocentric latitude (degree)')

  # ld
  ax = axs[5-i,1]
  ax.plot(lats, ld[:,i], 'k')
  ax.plot(-lats, ld[:,i], 'k')
  ax.set_ylabel('R$_45$ (%)')
  if i == 0:
    ax.set_xlabel('Planetocentric latitude (degree)')

savefig('figs/fig_modelA_lat.png', bbox_inches = 'tight')
