#! /usr/bin/env python3
import argparse, glob, os
from pylab import *
from netCDF4 import Dataset
from snapy.harp.utils import get_rt_bands

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input',
    required = True,
    choices = [x[:-8] for x in glob.glob('*.nc')] + [x[:-8] for x in glob.glob('data/*.nc')],
    help = 'true default atmospheric profiles'
    )
parser.add_argument('-t', '--truth',
    choices = ['none'] + [x[:-8] for x in glob.glob('*.nc')],
    default = 'none',
    help = 'true atmospheric profiles'
    )
parser.add_argument('--var',
    choices = ['nh3', 'tem'],
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
args = vars(parser.parse_args())
pmin, pmax = float(args['pmin']), float(args['pmax'])

if __name__ == '__main__':
# number of burn-in steps
  nburn = 100
# read atmospheric profiles
  data = Dataset('%s-mcmc.nc' % args['input'])
  if args['var'] == 'tem':
    var = data['temp'][:,:,3,:]
  elif args['var'] == 'nh3':
    var = data['vapor2'][:,:,3,:]*1.E3 # kg/kg -> g/kg
  pres = data['press'][0,:,0,0]/1.E5  # par -> bar
  ang = list(map(int, arccos(data['mu_out'][:])/pi*180.))
  i45 = ang.index(45)

  nstep, nlevel, nwalker = var.shape

# read auxiliary information from the input file
  radio_bands = get_rt_bands('%s.inp' % args['input'])
  freq = radio_bands[:,0]
  nfreq = len(freq)

  tb, ld = zeros((nfreq, nstep, nwalker)), zeros((nfreq, nstep, nwalker))
  for i in range(nfreq):
    tb[i,:,:] = data['b%dtoa' % (i+1,)][:,0,3,:]
    ld[i,:,:] = (data['b%dtoa' % (i+1,)][:,0,3,:] 
               - data['b%dtoa' % (i+1,)][:,i45,3,:])/tb[i,:,:]*100.
# tb0 is the baseline model
  tb0, ld0 = zeros(nfreq), zeros(nfreq)
  for i in range(nfreq):
    tb0[i] = data['b%dtoa' % (i+1,)][0,0,0,0]
    ld0[i] = (data['b%dtoa' % (i+1,)][0,0,0,0] 
            - data['b%dtoa' % (i+1,)][0,i45,0,0])/tb0[i]*100.
# tb is the anomaly with respect to the baseline
  tb -= tb0.reshape(nfreq,1,1)
  ld -= ld0.reshape(nfreq,1,1)

# mean of all walkers
  var_base = var[0,:,0]
  var_avg = mean(var[:,:,:], axis = 2)

# true ammonia
  if args['truth'] != 'none':
      data = Dataset('%s-main.nc' % args['truth'])
      if args['var'] == 'tem':
        var_truth = data['temp'][0,:,0,0]
      elif args['var'] == 'nh3':
        var_truth = data['vapor2'][0,:,0,0]*1.E3
      tb_truth, ld_truth = zeros(nfreq), zeros(nfreq)
      # tb_truth is the anomaly with respect to the baseline
      for i in range(nfreq):
        tb_truth[i] = data['b%dtoa' % (i+1,)][0,0,0,0]
        ld_truth[i] = (data['b%dtoa' % (i+1,)][0,0,0,0]
                     - data['b%dtoa' % (i+1,)][0,i45,0,0])/tb_truth[i]*100.
        tb_truth[i] -= tb0[i]
        ld_truth[i] -= ld0[i]

  fig, axs = subplots(2, 2, figsize = (12, 10),
    gridspec_kw = {'height_ratios':[1,4], 'width_ratios':[4,1]})
  subplots_adjust(hspace = 0.04, wspace = 0.04)

# brightness temperature anomaly
  ax = axs[0,0]
  tb_avg = mean(tb, axis = 2)
  ax.plot(range(nstep), zeros(nstep), '0.7', linewidth = 2)
  for i in range(nfreq):
    ax.plot(range(nstep), tb_avg[i], label = '%.1f GHz' % freq[i])
  ax.set_xlim([0, nstep-1])
  ax.set_ylabel("Tb' (K)", fontsize = 12)
  ax.xaxis.tick_top()
  ax.xaxis.set_label_position('top')
  ax.set_xlabel('MCMC step')
  ax.legend(ncol = nfreq, fontsize = 12)

# brightness temperature vs limb darkening
  ax = axs[0,1]
  tb_avg = mean(tb, axis = (1,2))
  tb_std = std(tb, axis = (1,2))
  ld_avg = mean(ld, axis = (1,2))
  ld_std = std(ld, axis = (1,2))
#ax.plot([-10, 12], [-10, 12], '0.7', linewidth = 2)
# true tb
  if args['truth'] != 'none':
      for i in range(nfreq):
        ax.plot(ld_truth[i], tb_truth[i], '^', ms = 5, alpha = 0.5, color = 'k')
  for i in range(nfreq):
    ax.errorbar(ld_avg[i], tb_avg[i], xerr = ld_std[i], yerr = tb_std[i])
  ax.xaxis.tick_top()
  ax.xaxis.set_label_position('top')
  ax.yaxis.tick_right()
  ax.yaxis.set_label_position('right')
  ax.set_xlabel("L$_d$ anomaly (%)")
  ax.set_ylabel("Tb anomaly (K)")
  #ax.plot([0, 0], ax.get_ylim(), color = '0.7')
  #ax.plot(ax.get_xlim(), [0, 0], color = '0.7')
  #ax.plot([0, 0], [-10.,10.], color = '0.7')
  #ax.plot([-10.,10.], [0, 0], color = '0.7')
  #ax.set_xlim([-10.,10.])
  ax.set_ylim(axs[0,0].get_ylim())

# variable profile sequence
  X, Y = meshgrid(range(nstep), pres)
  ax = axs[1,0]
  if args['var'] == 'tem':
    ax.contourf(X, Y, (var_avg - var_base).T, 20)
  else:
    ax.contourf(X, Y, var_avg.T, 20)
  ax.set_xlim([0, nstep-1])
  ax.set_ylim([pmax, pmin])
  ax.set_yscale('log')
  ax.set_xlabel('MCMC step')
  ax.set_ylabel('Pressure (bar)', fontsize = 12)

# mean of all profiles
  ax = axs[1,1]
  var_std = std(var[nburn:,:,:], axis = (0,2))
  var_avg = mean(var[nburn:,:,:], axis = (0,2))
  if args['var'] == 'tem':
    var_avg -= var_base
    var_truth -= var_base

  ax.plot(var_avg, pres)
  if args['var'] != 'tem':
    ax.plot(var_base, pres, '0.7')
  if args['truth'] != 'none':
    ax.plot(var_truth, pres, 'C3--')
  ax.fill_betweenx(pres, var_avg - var_std, var_avg + var_std, alpha = 0.5)
  ax.set_ylim([pmax, pmin])
  ax.set_yscale('log')
  if args['var'] == 'tem':
    ax.plot([0., 0.], [pmax, pmin], '0.7', linewidth = 2)
    ax.set_xlabel("T - T$_{ad}$ (K)")
    ax.set_xlim([-5., 5.])
  elif args['var'] == 'nh3':
    ax.set_xlabel('NH$_3$ mmr (g/kg)')
  ax.set_ylabel('Pressure (bar)', fontsize = 12)
  ax.yaxis.tick_right()
  ax.yaxis.set_label_position('right')

  savefig('%s-%s.png' % (os.path.basename(args['input']), args['var']), bbox_inches = 'tight')
