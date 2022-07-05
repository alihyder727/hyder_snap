#! /usr/bin/env python3
from pylab import *
import matplotlib
#matplotlib.pyplot.style.use('dark_background')

data = genfromtxt('data/gravity_correction_factor.txt')
rows = len(data)
data = data.reshape(6, rows//6, -1)
lat = data[0,:,0]

fig, axs = subplots(2, 3, figsize = (18, 10), sharex = True)
subplots_adjust(hspace = 0.08)

# CH1,2
ax = axs[0,0]
for i in [0,1]:
  ax.plot(lat, data[i,:,1], 'o-', label = 'CH%d' % (i+1))
ax.set_ylabel('T$_b$ correction factor')
ax.legend()

ax = axs[1,0]
for i in [0,1]:
  ax.plot(lat, data[i,:,2], 'o-', label = 'CH%d' % (i+1))
ax.legend()
ax.set_xlabel('Planetocentric latitude (degree)')
ax.set_ylabel('R$_{45}$ correction factor')

# CH3,4
ax = axs[0,1]
for i in [2,3]:
  ax.plot(lat, data[i,:,1], 'o-', label = 'CH%d' % (i+1))
ax.legend()

ax = axs[1,1]
for i in [2,3]:
  ax.plot(lat, data[i,:,2], 'o-', label = 'CH%d' % (i+1))
ax.legend()
ax.set_xlabel('Planetocentric latitude (degree)')

# CH5,6
ax = axs[0,2]
for i in [4,5]:
  ax.plot(lat, data[i,:,1], 'o-', label = 'CH%d' % (i+1))
#ax.set_ylabel('T$_b$ correction factor')
ax.legend()

ax = axs[1,2]
for i in [4,5]:
  ax.plot(lat, data[i,:,2], 'o-', label = 'CH%d' % (i+1))
ax.legend()
ax.set_xlabel('Planetocentric latitude (degree)')

#show()
savefig('figs/fig_gravity_correction.png', bbox_inches = 'tight')
