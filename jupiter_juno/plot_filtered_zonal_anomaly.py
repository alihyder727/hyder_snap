#! /usr/bin/env python3
import numpy.ma as ma
from pylab import *

pjs = ['PJ1', 'PJ3', 'PJ4', 'PJ5', 'PJ6', 'PJ7', 'PJ8', 'PJ9', 'PJ12']
data = genfromtxt('data/mwr_zonal_anomaly_flag.txt')
nrows, ncols = data.shape
# tb/ld, ch, lat, pj
data = data.reshape(6, 2, nrows//12, ncols)
lat = data[0,0,:,0]

X, Y = meshgrid(range(len(pjs)), lat)

fig, axs = subplots(2, 3, figsize = (12, 6), sharex = True, sharey = True)
subplots_adjust(hspace = 0.08, wspace = 0.08)
vmins = array([-20, -20, -15, -15, -10, -5])

# Tb
for i in range(6):
  tb = ma.masked_array(data[i,0,:,1:], mask = data[i,0,:,1:] == 0.)

  ax = axs[i//3,i%3]
  h = ax.pcolormesh(X, Y, tb, shading = 'nearest', cmap = 'plasma',
    vmin = vmins[i], vmax = -vmins[i])
  ax.set_ylim((-60, 60))
  ax.set_xticks(range(len(pjs)))
  ax.set_xticklabels(pjs)
  colorbar(h, ax = ax)
  if i%3 == 0:
    ax.set_ylabel('Planetocentric latitude (degree)')

savefig('figs/fig_filtered_zonal_tb_anomaly.png', bbox_inches = 'tight')
close()

fig, axs = subplots(2, 3, figsize = (12, 6), sharex = True, sharey = True)
subplots_adjust(hspace = 0.08, wspace = 0.08)
vmins = array([-1, -1, -1, -1, -1, -1])
# Ld
for i in range(6):
  ld = ma.masked_array(data[i,1,:,1:], mask = data[i,1,:,1:] == 0.)

  ax = axs[i//3,i%3]
  h = ax.pcolormesh(X, Y, ld, shading = 'nearest', cmap = 'plasma',
    vmin = vmins[i], vmax = -vmins[i])
  ax.set_ylim((-60, 60))
  ax.set_xticks(range(len(pjs)))
  ax.set_xticklabels(pjs)
  colorbar(h, ax = ax)
  if i%3 == 0:
    ax.set_ylabel('Planetocentric latitude (degree)')

savefig('figs/fig_filtered_zonal_ld_anomaly.png', bbox_inches = 'tight')
close()
