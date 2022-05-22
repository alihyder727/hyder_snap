#! /usr/bin/env python3
#import matplotlib
#matplotlib.use('Agg')
from pylab import *
from netCDF4 import Dataset
from astropy.io import fits
from matplotlib import gridspec
from matplotlib.ticker import NullFormatter
from snapy.harp.utils import get_sample_pressure, get_inversion_vars
import argparse, re, h5py, glob, os

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input',
    required = True,
    help = 'case file'
    )
parser.add_argument('--mode',
    choices = ['init', 'best'],
    default = 'init',
    help = 'select plotting mode'
    )
parser.add_argument('--pmax',
    default = '50.',
    help = 'maximum pressure'
    )
parser.add_argument('--pmin',
    default = '0.2',
    help = 'minimum pressure'
    )
parser.add_argument('--nburn',
    default = '500',
    help = 'burn-in steps'
    )
args = vars(parser.parse_args())

if __name__ == '__main__':
  pmin, pmax, p1bar = float(args['pmin']), float(args['pmax']), 1.
  colors = ['C0','C1', 'C2', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9']

  try:
    data = Dataset('%s-mcmc.nc' % args['input'])
  except FileNotFoundError:
    data = Dataset('%s-main.nc' % args['input'])
    if args['mode'] == 'best':
      raise ValueError('This input file does not support best mode')

  temp_ad = data['temp'][0,:,3,0]
  nh3_ad = data['vapor2'][0,:,3,0]*1.E3 # kg/kg -> g/kg

  temp_base = data['temp'][0,:,0,0]
  nh3_base = data['temp'][0,:,0,0]*1.E3 # kg/kg -> g/kg

  if args['mode'] == 'init':
    temp = data['temp'][0,:,0,0]
    temp_propose = data['temp'][0,:,1,0]
    nh3 = data['vapor2'][0,:,0,0]*1.E3 # kg/kg -> g/kg
    pres = data['press'][0,:,0,0]/1.E5  # pa -> bar
    ivar = []
  elif args['mode'] == 'best':
    psample = get_sample_pressure(args['input'] + '.inp')
    ivar = get_inversion_vars(args['input'] + '.inp')
    print(ivar)

    hdul = fits.open('%s.fits' % args['input'])
    par = hdul[0].data
    val = hdul[1].data
    lnp = hdul[2].data
    nstep, nwalker, ndim = par.shape
    istep, iwalker = unravel_index(argmax(lnp), lnp.shape)
    temp = data['temp'][istep,:,3,iwalker]
    nh3 = data['vapor2'][istep,:,3,iwalker]*1.E3 # kg/kg -> g/kg
    pres = data['press'][istep,:,3,iwalker]/1.E5  # pa -> bar
    nburn = int(args['nburn'])

    print(par[istep,iwalker,:])
    print('Best lnp = ', lnp[istep,iwalker])

  fig, axs = subplots(2, 4, figsize = (12, 8),
    gridspec_kw = {'height_ratios': [2,5], 'width_ratios': [4,1,4,1]})
  subplots_adjust(hspace = 0.08, wspace = 0.08)

# contribution function
  ff = h5py.File('contribution_function.h5', 'r')
  cf2 = ff['bellotti']['contribution_function_channel_2'][:]
  cf3 = ff['bellotti']['contribution_function_channel_3'][:]
  cf4 = ff['bellotti']['contribution_function_channel_4'][:]
  cf5 = ff['bellotti']['contribution_function_channel_5'][:]
  cf6 = ff['bellotti']['contribution_function_channel_6'][:]
  pcon = ff['bellotti']['pressure'][:]

# 1. temperature, 1 bar to pmin
  ax = axs[0,0]
  ax1 = ax.twiny()
  ax1.plot(cf5, pcon, 'k--', color = '0.7', linewidth = 2)
  ax1.plot(cf6, pcon, 'k--', color = '0.7', linewidth = 2)
  ax1.set_ylim([p1bar, pmin])
  ax1.xaxis.set_major_formatter(NullFormatter())

  ax.plot(temp_ad, pres, linewidth = 1, color = 'k')
  ax.plot(temp, pres, linewidth = 2, color = 'C3')
  ax.xaxis.tick_top()
  ax.xaxis.set_label_position('top')
  ax.set_xlim([100., 180.])
  ax.set_ylim([p1bar, pmin])
  ax.set_ylabel('Pressure (bar)')
  ax.set_xlabel('Temperature (K)')


# 2. temperature, pmax to 1 bar
  ax = axs[1,0]
  ax.plot(temp_ad, pres, linewidth = 1, color = 'k')
  ax.plot(temp, pres, linewidth = 2, color = 'C3')
  ax.set_xlim([100., 400.])
  ax.set_yscale('log')
  ax.set_ylim([pmax, p1bar])
  ax.set_xlabel('Temperature (K)')
  ax.set_ylabel('Pressure (bar)')

# 3. temperature anomaly, 1 bar to pmin
  ax = axs[0,1]
  #ax.plot(temp_propose - temp_ad, pres, linewidth = 1, color = 'C3', linestyle = '--')
  ax.plot(temp - temp_ad, pres, linewidth = 2, color = 'C3')
  ax.plot([0., 0.], [p1bar, pmin], 'k--')
  if 0 in ivar:
    xbest = par[istep,iwalker,:ndim//2]
    xerror = std(par[nburn:,:,:ndim//2], axis = (0,1))
    ax.errorbar(xbest, psample, xerr = xerror,
      capsize = 10, fmt = 'o', color = 'C3')
  ax.set_ylim([p1bar, pmin])
  ax.yaxis.set_major_formatter(NullFormatter())
  ax.xaxis.set_major_formatter(NullFormatter())

# 4. temperature anomaly, pmax to 1 bar
  ax = axs[1,1]
  #ax.plot(temp_propose - temp_ad, pres, linewidth = 1, color = 'C3', linestyle = '--')
  ax.plot(temp - temp_ad, pres, linewidth = 2, color = 'C3')
  ax.plot([0., 0.], [pmax, p1bar], 'k--')
  if 0 in ivar:
    xbest = par[istep,iwalker,:ndim//2]
    xerror = std(par[nburn:,:,:ndim//2], axis = (0,1))
    ax.errorbar(xbest, psample, xerr = xerror,
      capsize = 10, fmt = 'o', color = 'C3')
  ax.set_ylim([p1bar, pmin])
  ax.set_yscale('log')
  ax.set_ylim([pmax, p1bar])
  ax.yaxis.set_major_formatter(NullFormatter())
  ax.set_xlabel('T - T$_{ad}$ (K)')

# 5. NH3, H2O mixing ratio, 1 bar to pmin
  ax = axs[0,2]
  ax.plot(nh3_ad, pres, linewidth = 1, color = 'k')
  ax.plot(nh3, pres, linewidth = 2, color = 'C2')

  ax.set_xlabel('NH$3$ (g/kg)')
  ax.set_xscale('log')
  ax.set_xlim([1E-3, 10.])
  ax.set_ylim([p1bar, pmin])
  ax.xaxis.tick_top()
  ax.xaxis.set_label_position("top")
  #ax.set_xlabel('q (g/kg)')
  ax.yaxis.set_major_formatter(NullFormatter())

# 6. NH3, H2O mixing ratio, pmax to 1 bar
  ax = axs[1,2]
  ax.plot(nh3_ad, pres, linewidth = 1, color = 'k')
  ax.plot(nh3, pres, linewidth = 2, color = 'C2')

  ax.set_yscale('log')
  ax.set_xlim([2.0, 2.8])
  ax.set_ylim([pmax, p1bar])
  ax.set_xlabel('NH$3$ (g/kg)')
  ax.yaxis.set_major_formatter(NullFormatter())

  ax = axs[1,0].twiny()
  ax.plot(cf2, pcon, 'k--', color = '0.7', linewidth = 2)
  ax.plot(cf3, pcon, 'k--', color = '0.7', linewidth = 2)
  ax.plot(cf4, pcon, 'k--', color = '0.7', linewidth = 2)
  ax.plot(cf5, pcon, 'k--', color = '0.7', linewidth = 2)
  ax.plot(cf6, pcon, 'k--', color = '0.7', linewidth = 2)
  ax.set_ylim([pmax, p1bar])
  ax.xaxis.set_major_formatter(NullFormatter())

# 7. NH3 mixing ratio anomaly, 1 bar to pmin
  ax = axs[0,3]
  ax.plot(nh3 - nh3_ad, pres, linewidth = 2, color = 'C2')
  ax.plot([0., 0.], [p1bar, pmin], 'k--')
  if 2 in ivar:
    if 0 in ivar:
      xbest = par[istep,iwalker,ndim//2:]*1.E3
      xerror = std(par[nburn:,:,ndim//2:], axis = (0,1))*1.E3
      ax.errorbar(xbest, psample, xerr = xerror,
        capsize = 10, fmt = 'o', color = 'C2')

  ax.set_ylim([p1bar, pmin])
  ax.yaxis.set_major_formatter(NullFormatter())
  ax.xaxis.set_major_formatter(NullFormatter())

# 8. NH3 mixing ratio anomaly, pmax to 1 bar
  ax = axs[1,3]
  ax.plot(nh3 - nh3_ad, pres, linewidth = 2, color = 'C2')
  ax.plot([0., 0.], [pmax, p1bar], 'k--')
  if 2 in ivar:
    if 0 in ivar:
      xbest = par[istep,iwalker,ndim//2:]*1.E3
      xerror = std(par[nburn:,:,ndim//2:], axis = (0,1))*1.E3
      ax.errorbar(xbest, psample, xerr = xerror, 
        capsize = 10, fmt = 'o', color = 'C2')

  ax.set_yscale('log')
  #ax.set_xlim([2.0, 3.])
  ax.set_xlabel('q - q$_{ad}$ (g/kg)')
  ax.set_ylim([pmax, p1bar])
  ax.yaxis.set_major_formatter(NullFormatter())

  savefig('%s_profiles.png' % os.path.basename(args['input']), bbox_inches = 'tight')
  #show()
