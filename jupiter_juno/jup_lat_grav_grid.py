#! /usr/bin/env python3
from snapy.planet_gravity import jupiter_gravity
from pylab import *

clat = linspace(0., 85., 18)
grav = ['%.2f' % x for x in list(map(jupiter_gravity, clat))]
grav_grid = ' '.join(grav)
clat_grid = ' '.join(['%.2f' % x for x in clat])

with open('data/jup_latitude_grid.txt', 'w') as file:
  file.write(clat_grid)

with open('data/jup_gravity_grid.txt', 'w') as file:
  file.write(grav_grid)
