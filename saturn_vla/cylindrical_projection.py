#! /usr/bin/env python3
import matplotlib
matplotlib.use('Agg')
from pylab import *
from mpl_toolkits.axes_grid1 import make_axes_locatable
from common import read_info, regrid
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-f', '--fits',
  default = 'tbd',
  help = 'base name of the fits file without extension'
  )
parser.add_argument('-n', '--nadir',
  action = 'store_true',
  default = False,
  help = 'base name of the fits file without extension'
  )
parser.add_argument('--vmin',
  default = 'none',
  help = 'base name of the fits file without extension'
  )
parser.add_argument('--vmax',
  default = 'none',
  help = 'base name of the fits file without extension'
  )
parser.add_argument('-g',
  action = 'store_true',
  help = 'whether to regrid'
  )
args = vars(parser.parse_args())

print('Running cylindrical projection ...')

if args['nadir']:
  data = genfromtxt('../dat/%s-nadir.dat' % args['fits'])
  tb = data[:,6]
else:
  data = genfromtxt('../dat/%s.dat' % args['fits'])
  tb = data[:,5]

lat = data[:,1]
lon = data[:,3]
mu = data[:,4]
if args['g']:
  lon, lat, tb, mu = regrid(lon, lat, tb, mu)
  lon1 = hstack(([2*lon[0] - lon[1]], lon, [2*lon[-1] - lon[-2]]))
  lat1 = hstack(([2*lat[0] - lat[1]], lat, [2*lat[-1] - lat[-2]]))
  lon1 = (lon1[1:] + lon1[:-1])/2.
  lat1 = (lat1[1:] + lat1[:-1])/2.

ang = arccos(mu)/pi*180.

if args['vmax'] == 'none':
  vmax = tb.max()
else:
  vmax = float(args['vmax'])

if args['vmin'] == 'none':
  vmin = tb.max() - 20.
else:
  vmin = float(args['vmin'])

# observation point
try:
  info = read_info('../dat/%s.info' % args['fits'])
  ob_lat, ob_lon = info['oblat'], info['oblon']
  plot_ob = True
except:
  plot_ob = False

fig, axs = subplots(2, 1, figsize = (12, 8), sharex = True, sharey = True)
subplots_adjust(wspace = 0.08, hspace = 0.08)

ax = axs[0]
if args['g']:
  h1 = ax.pcolormesh(lon1, lat1, ang.T)
else:
  h1 = ax.scatter(lon, lat, c = ang, s = 0.2)
ax.set_xlim([0., 360.])
ax.set_ylabel('Planetographic Latitude', fontsize = 15)
if plot_ob:
  if lon.min() < 0.:
    ax.scatter(180. - ob_lon, ob_lat, c = 'r', s = 20)
  else:
    ax.scatter(ob_lon, ob_lat, c = 'r', s = 20)

# add colorbar
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size = "5%", pad = 0., aspect = 20)
colorbar(h1, cax = cax, label = 'Emission Angle')

ax = axs[1]
if args['g']:
  h2 = ax.pcolormesh(lon1, lat1, tb.T, cmap = 'gist_heat')
else:
  h2 = ax.scatter(lon, lat, c = tb, s = 0.2,
    cmap = 'gist_heat', vmin = vmin, vmax = vmax)
if lon.min() < 0.:
  ax.set_xlim([-180., 180.])
  ax.set_xlabel('East Longitude', fontsize = 15)
  if plot_ob:
    ax.scatter(180. - ob_lon, ob_lat, c = 'b', s = 20)
else:
  ax.set_xlim([360., 0.])
  ax.set_xlabel('West Longitude', fontsize = 15)
  if plot_ob:
    ax.scatter(ob_lon, ob_lat, c = 'b', s = 20)
ax.set_ylabel('Planetographic Latitude', fontsize = 15)
ax.grid('on')

# add colorbar
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size = "5%", pad = 0., aspect = 20)
colorbar(h2, cax = cax, label = 'Tb')

if args['nadir']:
  savefig('../img/%s-cylin-nadir.png' % args['fits'], dpi = 800,
    bbox_inches = 'tight')
else:
  savefig('../img/%s-cylin.png' % args['fits'], dpi = 800,
    bbox_inches = 'tight')
