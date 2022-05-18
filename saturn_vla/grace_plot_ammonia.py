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
# read atmospheric profiles
  data = Dataset('%s-mcmc.nc' % args['input'])
  nh3 = data['vapor2'][:,:,3,:]*1.E3 # kg/kg -> g/kg
  pres = data['press'][0,:,0,0]/1.E5  # par -> bar
  ang = list(map(int, arccos(data['mu_out'][:])/pi*180.))
  x1f = data['x1f'][:]
  i45 = ang.index(45)
  nstep, nlevel, nwalker = nh3.shape
  dtau1 = data['b1tau'][0,::-1,0,0]
  dtau2 = data['b2tau'][0,::-1,0,0]
  dtau3 = data['b3tau'][0,::-1,0,0]
  dtau4 = data['b4tau'][0,::-1,0,0]
  dtau5 = data['b5tau'][0,::-1,0,0]
  dtau6 = data['b6tau'][0,::-1,0,0]
  dZ = x1f[1:] - x1f[:-1]
  H0 = 50.E3

# normalized weighting functions

  tau = cumsum(dtau1)[::-1]
  wfunc = exp(-tau)*dtau1[::-1]/dZ*H0
  max_wfunc = max(wfunc)
  norm_wfunc_1 = wfunc/ max_wfunc

  tau = cumsum(dtau2)[::-1]
  wfunc = exp(-tau)*dtau2[::-1]/dZ*H0
  max_wfunc = max(wfunc)
  norm_wfunc_2 = wfunc/ max_wfunc

  tau = cumsum(dtau3)[::-1]
  wfunc = exp(-tau)*dtau3[::-1]/dZ*H0
  max_wfunc = max(wfunc)
  norm_wfunc_3 = wfunc/ max_wfunc

  tau = cumsum(dtau4)[::-1]
  wfunc = exp(-tau)*dtau4[::-1]/dZ*H0
  max_wfunc = max(wfunc)
  norm_wfunc_4 = wfunc/ max_wfunc

  tau = cumsum(dtau5)[::-1]
  wfunc = exp(-tau)*dtau5[::-1]/dZ*H0
  max_wfunc = max(wfunc)
  norm_wfunc_5 = wfunc/ max_wfunc

  tau = cumsum(dtau6)[::-1]
  wfunc = exp(-tau)*dtau6[::-1]/dZ*H0
  max_wfunc = max(wfunc)
  norm_wfunc_6 = wfunc/ max_wfunc

# read auxiliary information from the input file
  radio_bands = get_rt_bands('%s.inp' % args['input'])
  freq = radio_bands[:,0]
  nfreq = len(freq)

  tb, ld = zeros((nfreq, nstep, nwalker)), zeros((nfreq, nstep, nwalker))
  for i in range(nfreq):
    tb[i,:,:] = data['radiance'][:,0,3,:]

# tb0 is the baseline model
  tb0 = zeros(nfreq)
  for i in range(nfreq):
    tb0[i] = data['radiance'][0,0,0,0]
    
# tb is the anomaly with respect to the baseline
  tb -= tb0.reshape(nfreq,1,1)
 
# mean of all walkers
  nh3_base = nh3[0,:,0]
  nh3_avg = mean(nh3[:,:,:], axis = 2)

# true ammonia
  if args['truth'] != 'none':
      data = Dataset('%s-main.nc' % args['truth'])
      nh3_truth = data['vapor2'][0,:,0,0]*1.E3
      tb_truth = zeros(nfreq)
      # tb_truth is the anomaly with respect to the baseline
      for i in range(nfreq):
        tb_truth[i] = data['radiance'][0,0,0,0]
        tb_truth[i] -= tb0[i]
        
  fig, axs = subplots(2, 2, figsize = (12, 10),
    gridspec_kw = {'height_ratios':[1,4], 'width_ratios':[4,1]})
  subplots_adjust(hspace = 0.16, wspace = 0.04)

# brightness temperature anomaly
  ax = axs[0,0]
  tb_avg = mean(tb, axis = 2)
  ax.plot(range(nstep), zeros(nstep), '0.7', linewidth = 2)
  for i in range(nfreq):
    ax.plot(range(nstep), tb_avg[i], label = '%.1f GHz' % freq[i])
  ax.set_xlim([0, nstep-1])
  ax.set_ylabel("Tb' (K)", fontsize = 15)
  ax.xaxis.tick_top()
  ax.xaxis.set_label_position('top')
  ax.set_xlabel('MCMC step')
  ax.legend(ncol = nfreq, fontsize = 8)

# brightness temperature vs limb darkening
  ax = axs[0,1]
  tb_avg = mean(tb, axis = (1,2))
  tb_std = std(tb, axis = (1,2))
  
#ax.plot([-10, 12], [-10, 12], '0.7', linewidth = 2)
# true tb
  if args['truth'] != 'none':
      for i in range(nfreq):
        ax.plot(tb_truth[i], '^', ms = 5, alpha = 0.5, color = 'k')
  for i in range(nfreq):
    ax.errorbar(tb_avg[i], yerr = tb_std[i])
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

# ammonia profile sequence
  X, Y = meshgrid(range(nstep), pres)
  ax = axs[1,0]
  ax.contourf(X, Y, nh3_avg.T, 20)
  ax.set_xlim([0, nstep-1])
  ax.set_ylim([pmax, pmin])
  ax.set_yscale('log')
  ax.set_xlabel('MCMC step')
  ax.set_ylabel('Pressure (bar)', fontsize = 12)

# mean of all profiles
  ax = axs[1,1]
  nh3_std = std(nh3, axis = (0,2))
  nh3_avg = mean(nh3, axis = (0,2))
  ax.plot(nh3_base, pres, '0.7')
  ax.plot(nh3_avg, pres)
  if args['truth'] != 'none':
    ax.plot(nh3_truth, pres, 'C3--')
  ax.fill_betweenx(pres, nh3_avg - nh3_std, nh3_avg + nh3_std, alpha = 0.5)
  ax.set_ylim([pmax, pmin])
  ax.set_yscale('log')
  ax.set_xlabel('NH$_3$ mmr (g/kg)')
  ax.set_ylabel('Pressure (bar)', fontsize = 12)
  ax.yaxis.tick_right()
  ax.yaxis.set_label_position('right')

# normalized weighting func on same plot as mean of all profiles
  ax = ax.twiny()
  ax.plot(norm_wfunc_1, pres, label='2.1 GHz', alpha=0.35)
  ax.plot(norm_wfunc_2, pres, label='4.1 GHz', alpha=0.35)
  ax.plot(norm_wfunc_3, pres, label='10.0 GHz', alpha=0.35)
  ax.plot(norm_wfunc_4, pres, label='15.0 GHz', alpha=0.35)
  ax.plot(norm_wfunc_5, pres, label='22.0 GHz', alpha=0.35)
  ax.plot(norm_wfunc_6, pres, label='44.2 GHz', alpha=0.35)
  ax.legend(fontsize = 8)
  ax.set_xlabel('Normalized Weighting Function', fontsize=12)
  

  savefig('%s-nh3.png' % os.path.basename(args['input']), bbox_inches = 'tight')
