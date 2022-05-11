#! /usr/bin/env python3
import argparse, glob, os
from pylab import *
from snapy.harp.utils import athinput, get_rt_bands
from astropy.io import fits
from plot_mcmc_profile import read_variable_simulate

def plot_variable_pdf(name, var, pres, unit = 'ppmv'):
  nburn = int(args['nburn'])
  levels = list(map(float, args['levels'].split()))

  nstep, nlevel, nwalker = var[nburn:,:,:].shape
  fig, axs = subplots(2, 1, figsize = (6, 8), sharex = True)
  subplots_adjust(hspace = 0.08)

  for i in range(len(levels)):
    j = (abs(pres - levels[i])).argmin()
    cvar = var[nburn:,j,:].flatten()
    hist, edges = histogram(cvar, bins = 40, density = True)
    hist[1:-1] = (hist[1:-1] + hist[:-2] + hist[2:])/3.
    if i < 3:
      ax = axs[0]
    else :
      ax = axs[1]
    ax.hist(cvar, bins = 40, alpha = 0.5, color = 'C%d' % i, 
        density = True, label = '%s bar' % levels[i])
    ax.plot((edges[:-1] + edges[1:])/2., hist, color = 'C%d' % i, linewidth = 2)
    ax.set_ylabel('Probability density')

  axs[0].legend()
  axs[1].legend()
  ax = axs[1]
  if name == 'nh3':
    ax.set_xlabel('NH$_3$ concentration (%s)' % unit)
  else :
    ax.set_xlabel('Temperature (K)')

  outfig = '%s/%s-%s-pdf.png' % (args['dir'], os.path.basename(args['input']), name)
  print('figure saved to %s' % outfig)
  savefig('%s' % outfig, bbox_inches = 'tight')
  close()

def write_variable_pdf(name, var, pres, outfile):
  nburn = int(args['nburn'])
  levels = list(map(float, args['levels'].split()))

  with open(outfile, 'a') as file:
    file.write('# %s statistics:\n' % name.upper())
    if name == 'nh3':
      file.write('# Pressure level (bar), NH3 mean (g/kg), NH3 std (g/kg), NH3 mean (ppmv), NH3 std (ppmv)\n')
    else :
      file.write('# Pressure level (bar), T mean (K), T anomaly (K), T std (K)\n')

    for i in range(len(levels)):
      j = (abs(pres - levels[i])).argmin()
      cvar = var[nburn:,j,:].flatten()

      if name == 'nh3':
        file.write('%10s%10.2f%10.2f%10.2f%10.2f\n' % (levels[i], mean(cvar)/350.*2.7,
          std(cvar)/350.*2.7, mean(cvar), std(cvar)))
      else :
        file.write('%10s%10.2f%10.2f%10.2f\n' % (levels[i], mean(cvar) + var_ad[j],
          mean(cvar), std(cvar)))
    print('statistics written to %s' % outfile)

if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument('-i', '--input',
      required = True,
      help = 'mcmc case to plot'
      )
  parser.add_argument('-d', '--dir',
      default = '.',
      help = 'save directory'
      )
  parser.add_argument('--var',
      default = 'nh3',
      help = 'which variable to plot'
      )
  parser.add_argument('--nburn',
      default= '1000',
      help = 'number of burn-in steps'
      )
  parser.add_argument('--levels',
      default = '100 50 20 10 5 2 1',
      help = 'pressure levels'
      )
  parser.add_argument('--unit',
      default = 'ppmv',
      help = 'ammonia concentration unit'
      )
  args = vars(parser.parse_args())

  var_list = args['var'].split(',')
  outfile = '%s/%s-pdf.txt' % (args['dir'], os.path.basename(args['input']))

  with open(outfile, 'w') as file:
    file.write('# Variability from simulation %s\n' % os.path.basename(args['input']))

  # read atmospheric profiles
  for name in var_list:
    var_ad, var_base, var, pres = read_variable_simulate(name, args['input'])
    if name == 'tem':
      var -= var_ad.reshape((1,len(var_ad),1))
    else :
      var *= 350./2.7
    write_variable_pdf(name, var, pres, outfile)
    plot_variable_pdf(name, var, pres, unit = 'ppmv')
