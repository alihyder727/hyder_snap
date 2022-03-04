#! /usr/bin/env python3
from matplotlib.gridspec import GridSpec
from pylab import *
from astropy.io import fits
from snapy.harp.utils import get_sample_pressure

fname = 'vla_saturn-test'

# read parameters in input file
pres = get_sample_pressure(fname + '.inp')

hdul = fits.open('%s.fits' % fname)
par = hdul[0].data*1.E3 # kg/kg -> g/kg
val = hdul[1].data
lnp = hdul[2].data
msk = hdul[3].data

nstep, nwalker, ndim = par.shape
std_par, avg_par = zeros(ndim), zeros(ndim)
for i in range(ndim):
  std_par[i] = std(par[:,:,i].flatten())
  avg_par[i] = mean(par[:,:,i].flatten())
print('mean = ', avg_par)
print('std = ', std_par)

fig = figure(figsize = (12, 10))
gs = GridSpec(ndim + 1, ndim - 1, figure = fig)
subplots_adjust(hspace = 0.08, wspace = 0.08)

# correlation plot
for i in range(ndim - 1):
  ax = fig.add_subplot(gs[0, i])
  ax.plot([0, 0], [-5, 5], '0.7', linewidth = 2)
  ax.plot([-5, 5], [0, 0], '0.7', linewidth = 2)
  for j in range(nwalker):
    ax.scatter((par[:,j,i] - avg_par[i])/std_par[i],
      (par[:,j,i+1] - avg_par[i+1])/std_par[i+1], s = 1, alpha = 0.5)
  ax.set_xlabel('%.1f bar ($\sigma$)' % pres[i])
  ax.xaxis.tick_top()
  ax.xaxis.set_label_position('top')
  ax.set_xlim([-5, 5])
  ax.set_ylim([-5, 5])
  if i == 0:
    ax.set_ylabel('%.1f bar ($\sigma$)' % pres[i+1])
  elif i == ndim - 2:
    ax.set_ylabel('%.1f bar ($\sigma$)' % pres[i+1])
    ax.yaxis.tick_right()
    ax.yaxis.set_label_position('right')
  else:
    ax.set_yticklabels([])

# mcmc step plot
for i in range(ndim):
  ax = fig.add_subplot(gs[-(i+1), :-1])
  for j in range(nwalker):
    ax.plot(range(len(par)), par[:,j,i])
  ax.set_ylabel('%.1f bar (g/kg)' % pres[i])
  ax.set_xlim([0, nstep])
  if i == 0:
    ax.set_xlabel('MCMC Steps')
  else:
    ax.xaxis.set_ticklabels([])

# histgram
for i in range(ndim):
  #hist, edges = histogram(par[:,:,i].flatten(), bins = 10)
  ax = fig.add_subplot(gs[-(i+1), -1])
  ax.hist(par[:,:,i].flatten(), orientation = 'horizontal', bins = 20)
  ax.yaxis.tick_right()
  if i != 0:
    ax.set_xticklabels([])

savefig('%s-mcmc.png' % fname, bbox_inches = 'tight')
