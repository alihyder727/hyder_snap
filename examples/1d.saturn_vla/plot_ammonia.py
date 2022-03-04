#! /usr/bin/env python3
import argparse
from pylab import *
from netCDF4 import Dataset
from snapy.harp.utils import get_rt_bands

parser = argparse.ArgumentParser()
parser.add_argument('-t', '--truth',
    default = 'none',
    help = 'true atmospheric profiles'
    )
args = vars(parser.parse_args())

fname = 'vla_saturn-test'
#fname = 'vla_saturn_real-q54'
#ftruth = 'vla_ideal_saturn-n6'

# read atmospheric profiles
data = Dataset('%s-mcmc.nc' % fname)
nh3 = data['vapor2'][:,:,3,:]*1.E3 # kg/kg -> g/kg
pres = data['press'][0,:,0,0]/1.E5  # par -> bar
ang = list(map(int, arccos(data['mu_out'][:])/pi*180.))
i45 = ang.index(45)
nstep, nlevel, nwalker = nh3.shape

# read auxiliary information from the input file
radio_bands = get_rt_bands('%s.inp' % fname)
freq = radio_bands[:,0]
nfreq = len(freq)

tb = zeros((nfreq, nstep, nwalker))
for i in range(nfreq):
  tb[i,:,:] = data['b%dtoa' % (i+1,)][:,0,3,:]
  ld[i,:,:] = (data['b%dtoa' % (i+1,)][:,0,3,:] - data['b%dtoa
# tb0 is the baseline model
tb0 = zeros(nfreq)
for i in range(nfreq):
  tb0[i] = data['b%dtoa' % (i+1,)][0,0,0,0]
# tb is the anomaly with respect to the baseline
tb -= tb0.reshape(nfreq,1,1)
#tb[:,0,:] = 0.

# mean of all walkers
nh3_base = nh3[0,:,0]
nh3_avg = mean(nh3[:,:,:], axis = 2)

# true ammonia
if args['truth'] != 'none':
    data = Dataset('data/%s-main.nc' % args['truth'])
    nh3_truth = data['vapor2'][0,:,0,0]*1.E3
    tb_truth = zeros(nfreq)
    # tb_truth is the anomaly with respect to the baseline
    for i in range(nfreq):
      tb_truth[i] = data['b%dtoa' % (i+1,)][0,0,0,0] - data['b%dtoa' % (i+1,)][0,0,3,0]

fig, axs = subplots(2, 2, figsize = (12, 10),
  gridspec_kw = {'height_ratios':[1,4], 'width_ratios':[4,1]})
subplots_adjust(hspace = 0.08, wspace = 0.08)

# brightness temperature anomaly
ax = axs[0,0]
tb_avg = mean(tb, axis = 2)
ax.plot(range(nstep), zeros(nstep), '0.7', linewidth = 2)
for i in range(6):
  ax.plot(range(nstep), tb_avg[i], label = '%.1f GHz' % freq[i])
ax.set_xlim([0, nstep])
ax.set_ylabel("Tb' (K)", fontsize = 15)
ax.xaxis.tick_top()
ax.xaxis.set_label_position('top')
ax.set_xlabel('MCMC step')
ax.legend(ncol = nfreq, fontsize = 8)

# brightness temperature vs limb darkening
ax = axs[0,1]
tb_avg = mean(tb, axis = (1,2))
tb_std = std(tb, axis = (1,2))
ax.plot([-10, 12], [-10, 12], '0.7', linewidth = 2)
for i in range(6):
  ax.errorbar(tb_truth[i], tb_avg[i], xerr = 1., yerr = tb_std[i])
ax.set_xlim([0, 12])
ax.set_ylim([0, 12])
ax.xaxis.tick_top()
ax.xaxis.set_label_position('top')
ax.yaxis.tick_right()
ax.yaxis.set_label_position('right')
ax.set_xlabel("Tb' measure (K)")
ax.set_ylabel("Tb' simulate (K)")

X, Y = meshgrid(range(nstep), pres)
ax = axs[1,0]
ax.contourf(X, Y, nh3_avg.T, 20)
ax.set_xlim([0, nstep])
ax.set_ylim([100., 0.2])
ax.set_yscale('log')
ax.set_xlabel('MCMC step')
ax.set_ylabel('Pressure (bar)', fontsize = 12)

ax = axs[1,1]
nh3_std = std(nh3, axis = (0,2))
nh3_avg = mean(nh3, axis = (0,2))
ax.plot(nh3_base, pres, '0.7')
ax.plot(nh3_avg, pres)
ax.plot(nh3_truth, pres, 'C3--')
ax.fill_betweenx(pres, nh3_avg - nh3_std, nh3_avg + nh3_std, alpha = 0.5)
ax.set_ylim([100., 0.2])
ax.set_yscale('log')
ax.set_xlabel('NH$_3$ mmr (g/kg)')
ax.set_ylabel('Pressure (bar)', fontsize = 12)
ax.yaxis.tick_right()
ax.yaxis.set_label_position('right')

savefig('figs/%s-nh3.png' % fname)
