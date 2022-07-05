#! /usr/bin/env python3
from pylab import *
from glob import glob
from netCDF4 import Dataset

# ModelA-lat
clat_grid = genfromtxt('data/jup_latitude_grid.txt', dtype = str)

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


output = 'data/gravity_correction_factor.txt'
with open(output, 'w') as file:
  file.write('# Juno MWR Tb gravity correction factor\n')
  file.write('# Planetocentric latitude, Tb correction factor to equator, Ld correction factor to equator\n')
  for i in range(6):
    file.write('# Channel %d\n' % (1+i,))
    for j in range(len(lats)):
      file.write('%8.2f %10.4f %10.4f\n' % (lats[j], tb0[0,i]/tb0[j,i],
        ld[0,i]/ld[j,i]))
print('Output file written to %s' % output)
