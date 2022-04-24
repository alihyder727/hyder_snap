#! /usr/bin/env python3
import argparse, glob, os
from pylab import *
from snapy.harp.utils import athinput, get_rt_bands
from astropy.io import fits
#import matplotlib
#matplotlib.use('Agg')

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input',
    required = True,
    help = 'mcmc case to plot'
    )
parser.add_argument('-d', '--dir',
    default = '.',
    help = 'save directory'
    )
parser.add_argument('--nburn',
    default= '1000',
    help = 'number of burn-in steps'
    )
args = vars(parser.parse_args())

if __name__ == '__main__':
  # read auxiliary information from the input file
  radio_bands = get_rt_bands('%s.inp' % args['input'])
  freq = radio_bands[:,0]
  nfreq = len(freq)

  # read simulation
  nburn = int(args['nburn'])
  hdul = fits.open('%s.fits' % args['input'])
  val = hdul[1].data
  nstep, nwalker, nval = val.shape
  val = val[nburn:,:,:].reshape((-1, nval))
  val_avg = mean(val, axis = 0)

  # read observation
  obsfile = athinput('%s.inp' % args['input'])['inversion']['obsfile']
  if obsfile != 'none':
    dirname = os.path.dirname(args['input'])
    if dirname == '': dirname = '.'
    # fit target
    target = genfromtxt('%s/%s' % (dirname, obsfile), skip_header = 2, max_rows = 12)
    tb_truth = target[::2]
    ld_truth = target[1::2]

    # covariance matrix
    icov = genfromtxt('%s/%s' % (dirname, obsfile), skip_header = 14)
    dicov = diag(icov)
    chi2 = (target - val_avg)*(target - val_avg)*diag(icov)

  # plot sample and chi2
  cov = sqrt(1./dicov)
  fig, axs = subplots(2, 3, figsize = (16, 10))
  subplots_adjust(hspace = 0.3)
  for i in range(nfreq):
    ax = axs[i//3, i%3]
    j = 2*i
    hist, X, Y = histogram2d(val[:,j], val[:,j+1], bins = 20)
    hflat = hist.flatten()
    inds = argsort(hflat)[::-1]
    hflat = hflat[inds]
    sm = cumsum(hflat)
    sm /= sm[-1]
    vlist = []
    for v in [0.95, 0.68]:
      vlist.append(hflat[sm <= v][-1])
    X1, Y1 = 0.5*(X[1:] + X[:-1]), 0.5*(Y[1:] + Y[:-1])
    ax.scatter(val[:,j], val[:,j+1], s = 0.5, c = '0.2', alpha = 0.5)
    ax.contour(X1, Y1, hist.T, vlist, linewidths = 2, colors = 'C1')
    ax.errorbar(target[j], target[j+1], xerr = cov[j], yerr = cov[j+1],
        capsize = 5, linewidth = 2)
    ax.set_title('%s GHz, $\chi^2$=%.2f + %.2f' % (freq[i], chi2[j], chi2[j+1]))
    ax.set_xlabel('T$_b$ nadir (K)', fontsize = 12)
    ax.set_ylabel('R$_{45}$ (%)', fontsize = 12)

  outname = '%s/%s-chi2.png' % (args['dir'], os.path.basename(args['input']))
  print('figure saved to %s' % outname)
  savefig('%s' % outname, bbox_inches = 'tight')
