#! /usr/bin/env python3
from pylab import *
from netCDF4 import Dataset
import argparse
import matplotlib, re

latmin, latmax = -60., 60.
pjs = [1,3,4,5,6,7,8,9,12]
freqs = [0.6, 1.25, 2.6, 5.2, 10., 22.]

fig, axs = subplots(6, 1, figsize = (12, 8), sharex = True)

# plot observed data from PJ1 - PJ12
for j in range(6):
  # read filtered data
  data = genfromtxt('data/ch%d_gfilter.txt' % (j+1,))
  nrows, ncols = data.shape
  data = data.reshape((4, nrows//4, ncols))
  lat = data[0,:,0]
  tb = data[0,:,1:]
  eps_tb = data[1,:,1:]
  ld = data[2,:,1:]
  eps_ld = data[3,:,1:]

  ax = axs[5-j]
  for i, pj in enumerate(pjs):
    ix = (abs(eps_tb[:,i]) < 5.) & (abs(eps_ld[:,i]) < 0.5)
    #ax.step(lat[ix], tb[ix,i], where = 'mid', label = 'PJ%d' % pj)
    ax.step(lat[ix], ld[ix,i], where = 'mid', label = 'PJ%d' % pj)

# plot Model A data
data = Dataset('data/juno_mwr-modelA.nc', 'r')
tb = data['radiance'][0,:,0,0]
tb0 = tb[::4]
tb45 = tb[3::4]
ld = (tb0 - tb45)/tb0*100.

for j in range(6):
  ax = axs[5-j]
  #ax.plot([latmin, latmax], [tb0[j], tb0[j]], '0.7', linewidth = 4)
  ax.plot([latmin, latmax], [ld[j], ld[j]], '0.7', linewidth = 4)
  ax.set_ylabel('%s GHz' % freqs[j])

ax = axs[5]
ax.set_xlabel('Planetocentric Latitude (degree)', fontsize = 12)
ax.set_xlim([latmax, latmin])
ax = axs[4]
ax.legend(ncol = 5, fontsize = 8)

#show()
savefig('figs/fig_ld_gfilter.png', bbox_inches = 'tight')
