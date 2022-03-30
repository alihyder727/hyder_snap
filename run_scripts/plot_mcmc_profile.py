#! /usr/bin/env python3
import argparse, glob, os
from pylab import *
from netCDF4 import Dataset
from snapy.harp.utils import *
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib
matplotlib.use('Agg')

parser = argparse.ArgumentParser()
parser.add_argument('-d',
    action = 'store_true',
    help = 'plot differential'
    )
parser.add_argument('-i', '--input',
    required = True,
    help = 'true default atmospheric profiles'
    )
parser.add_argument('--truth',
    default = 'none',
    help = 'true atmospheric profiles'
    )
parser.add_argument('--dir',
    default = '.',
    help = 'save directory'
    )
parser.add_argument('--var',
    default = 'nh3',
    help = 'which variable to plot'
    )
parser.add_argument('--pmax',
    default = '100.',
    help = 'maximum pressure'
    )
parser.add_argument('--pmin',
    default = '0.2',
    help = 'minimum pressure'
    )
parser.add_argument('--nburn',
    default= '100',
    help = 'number of burn-in steps'
    )
args = vars(parser.parse_args())
pmin, pmax = float(args['pmin']), float(args['pmax'])

def read_variable_simulate(name, case, out = 'out2'):
  try :
    data = Dataset('%s-mcmc.nc' % case)
    pres = data['press'][0,:,0,0]/1.E5  # par -> bar
  except FileNotFoundError :
    data = Dataset('%s.out2.nc' % case)
    pres = data['press'][0,:,0,0]/1.E5  # par -> bar
    if name == 'tem':
      data = Dataset('%s.out3.nc' % case)
    else :
      data = Dataset('%s.out2.nc' % case)

  # var : (ntime, nlevel, nmodel, nwalker)
  if name == 'tem' :
    var = data['temp'][:,:,3,:]
    var_ad = data['temp'][0,:,3,0]
    var_base = data['temp'][0,:,0,0]
  elif name == 'nh3' : # kg/kg -> g/kg
    var = data['vapor2'][:,:,3,:]*1.E3
    var_ad = data['vapor2'][0,:,3,0]*1.E3
    var_base = data['vapor2'][0,:,0,0]*1.E3
  else :
    raise ValueError("Variable %s not found" % name)

  return var_ad, var_base, var, pres

def read_variable_truth(name, case, out = 'out2'):
  try :
    data = Dataset('%s-main.nc' % case)
  except FileNotFoundError :
    data = Dataset('%s.%s.nc' % (case, out))

  if name == 'tem':
    var = data['temp'][0,:,0,0]
  elif name == 'nh3':
    var = data['vapor2'][0,:,0,0]*1.E3
  else :
    raise ValueError("Variable %s not found" % name)

  return var

def read_tbld_truth(case, nfreq, i45, out = 'out4'):
  try :
    data = Dataset('%s-main.nc' % case)
  except FileNotFoundError :
    data = Dataset('%s.%s.nc' % (case, out))

  tb, ld = [], []
  # tb_truth is the anomaly with respect to the baseline
  for i in range(nfreq):
    tb.append(data['b%dtoa' % (i+1,)][0,0,0,0])
    ld.append((data['b%dtoa' % (i+1,)][0,0,0,0]
      - data['b%dtoa' % (i+1,)][0,i45,0,0])/tb[i]*100.)
  return tb, ld

def read_tbld_simulate(case, nfreq, i45, out = 'out4'):
  try :
    data = Dataset('%s-mcmc.nc' % case)
  except FileNotFoundError :
    data = Dataset('%s.%s.nc' % (case, out))

  # tb_base is the baseline model
  tb_base, ld_base = [], []
  for i in range(nfreq):
    tb_base.append(data['b%dtoa' % (i+1,)][0,0,0,0])
    ld_base.append((data['b%dtoa' % (i+1,)][0,0,0,0] 
        - data['b%dtoa' %(i+1,)][0,i45,0,0])/tb_base[i]*100.)
  tb_base, ld_base = array(tb_base), array(ld_base)

  # tb_ad is the adiabatic model
  tb_ad, ld_ad = [], []
  for i in range(nfreq):
    tb_ad.append(data['b%dtoa' % (i+1,)][0,0,3,0])
    ld_ad.append((data['b%dtoa' % (i+1,)][0,0,3,0] 
        - data['b%dtoa' %(i+1,)][0,i45,3,0])/tb_ad[i]*100.)
  tb_ad, ld_ad = array(tb_ad), array(ld_ad)

  tb, ld = [], []
  for i in range(nfreq):
    tb.append(data['b%dtoa' % (i+1,)][:,0,3,:])
    ld.append((data['b%dtoa' % (i+1,)][:,0,3,:] 
      - data['b%dtoa' % (i+1,)][:,i45,3,:])/tb[i]*100.)
  tb, ld = array(tb), array(ld)

  return tb_ad, tb_base, tb, ld_ad, ld_base, ld

def plot_mcmc_profile(name, nburn):
  # read atmospheric profiles
  var_ad, var_base, var, pres = read_variable_simulate(name, args['input'])

  # read auxiliary information from the input file
  radio_bands = get_rt_bands('%s.inp' % args['input'])
  amu, aphi = get_ray_out('%s.inp' % args['input'])
  freq = radio_bands[:,0]
  nfreq = len(freq)
  i45 = list(map(int, amu)).index(45)
  
  # read mwr tbld
  tb_ad, tb_base, tb, ld_ad, ld_base, ld = read_tbld_simulate(args['input'], nfreq, i45)

  # true ammonia
  if args['truth'] != 'none':
    # read truth variable
    var_truth = read_variable_truth(name, args['truth'])
    # read truth tbld
    tb_truth, ld_truth = read_tbld_truth(args['truth'], nfreq, i45)

  obsfile = athinput('%s.inp' % args['input'])['inversion']['obsfile']
  if obsfile != 'none':
    dirname = os.path.dirname(args['input'])
    if dirname == '': dirname = '.'
    data = genfromtxt('%s/%s' % (dirname, obsfile), max_rows = 13)
    tb_truth = data[1::2]
    ld_truth = data[2::2]

  # make plots
  fig, axs = subplots(2, 2, figsize = (12, 10),
    gridspec_kw = {'height_ratios':[1,4], 'width_ratios':[4,1]})
  subplots_adjust(hspace = 0.04, wspace = 0.04)

  # preprocess
  nstep, nlevel, nwalker = var.shape

  # tb -> tb anomaly, ld -> ld anomaly
  tb -= tb_base.reshape(nfreq,1,1)
  ld -= ld_base.reshape(nfreq,1,1)
  if args['truth'] != 'none':
    tb_truth -= tb_base
    ld_truth -= ld_base
  if obsfile != 'none':
    diff = athinput('%s.inp' % args['input'])['inversion']['differential']
    if diff == 'false':
      tb_truth -= tb_base
      ld_truth -= ld_base

  # plot differential
  if args['d']:
    var -= var_base.reshape(1,nlevel,1)
    if args['truth'] != 'none':
      var_truth -= var_base.reshape(1,nlevel,1)

  # brightness temperature anomaly
  ax = axs[0,0]
  # average over walkers
  tb_avg = mean(tb, axis = 2)
  ax.plot(range(nstep), zeros(nstep), '0.7', linewidth = 2)
  for i in range(nfreq):
    ax.plot(range(nstep), tb_avg[i], label = '%.1f GHz' % freq[i], color = 'C%d' % (i+1))
  ax.set_xlim([0, nstep-1])
  ax.set_ylabel("Tb' (K)", fontsize = 12)
  ax.xaxis.tick_top()
  ax.xaxis.set_label_position('top')
  ax.set_xlabel('MCMC steps')
  ax.legend(ncol = nfreq//2, fontsize = 8)

  # brightness temperature vs limb darkening
  ax = axs[0,1]
  # average over time and walker, excluding burn-in
  tb_avg = mean(tb[:,nburn:,:], axis = (1,2))
  tb_std = std(tb[:,nburn:,:], axis = (1,2))
  ld_avg = mean(ld[:,nburn:,:], axis = (1,2))
  ld_std = std(ld[:,nburn:,:], axis = (1,2))
  for i in range(nfreq):
    ax.errorbar(ld_avg[i], tb_avg[i], xerr = ld_std[i], yerr = tb_std[i], color = 'C%d' % (i+1))

  # plot truth or observation
  if args['truth'] != 'none' or obsfile != 'none':
    for i in range(nfreq):
      ax.plot(ld_truth[i], tb_truth[i], '^', ms = 10, alpha = 0.5, color = 'C%d' % (i+1))

  ax.xaxis.tick_top()
  ax.xaxis.set_label_position('top')
  ax.yaxis.tick_right()
  ax.yaxis.set_label_position('right')
  ax.set_xlabel("L$_d$ anomaly (%)")
  ax.set_ylabel("Tb anomaly (K)")
  ax.set_ylim(axs[0,0].get_ylim())

  # variable profile sequence
  # average over walkers
  var_avg = mean(var[:,:,:], axis = 2)

  X, Y = meshgrid(range(nstep), pres)
  ax = axs[1,0]
  h1 = ax.contourf(X, Y, var_avg.T, 20, cmap = 'bwr')

  # add colorbar
  divider = make_axes_locatable(ax)
  cax = divider.append_axes("right", size = "10%", pad = 0., aspect = 40)
  colorbar(h1, cax = cax)

  ax.set_xlim([0, nstep-1])
  ax.set_ylim([pmax, pmin])
  ax.set_yscale('log')
  ax.set_xlabel('MCMC steps')
  ax.set_ylabel('Pressure (bar)', fontsize = 12)

  # average over time and walker, excluding burn-in
  ax = axs[1,1]
  var_std = std(var[nburn:,:,:], axis = (0,2))
  var_avg = mean(var[nburn:,:,:], axis = (0,2))
  if not args['d']:
    ax.plot(var_base, pres, '0.7')
  else :
    ax.plot([0., 0.], [pmax, pmin], 'k')
    if args['truth'] != 'none':
      ax.plot(var_truth, pres, 'C3--')

  ax.plot(var_avg, pres)
  ax.fill_betweenx(pres, var_avg - var_std, var_avg + var_std, alpha = 0.5)

  ax.set_ylim([pmax, pmin])
  ax.set_yscale('log')
  if name == 'tem':
    ax.plot([0., 0.], [pmax, pmin], '0.7', linewidth = 2)
    ax.set_xlabel("T (K)")
    ax.set_xlim([-5., 10.])
  elif name == 'nh3':
    ax.set_xlabel('NH$_3$ mmr (g/kg)')
  ax.set_ylabel('Pressure (bar)', fontsize = 12)
  ax.yaxis.tick_right()
  ax.yaxis.set_label_position('right')

  outname = '%s/%s-%s-profile.png' % (args['dir'], os.path.basename(args['input']), name)
  print('figure saved to %s' % outname)
  savefig('%s' % outname, bbox_inches = 'tight')
  close()

if __name__ == '__main__':
  # number of burn-in steps
  var_list = args['var'].split(',')
  args['nburn'] = int(args['nburn'])

  for var in var_list:
    plot_mcmc_profile(var, int(args['nburn']))
