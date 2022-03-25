#! /usr/bin/env python3
from pylab import *
import argparse, os

parser = argparse.ArgumentParser(
    formatter_class = argparse.ArgumentDefaultsHelpFormatter
    )
parser.add_argument('-i', '--input',
    default = 'none',
    help = 'synthetic brightness temperature output file'
    )
parser.add_argument('--ecal',
    default = '0.02',
    help = 'calibration accuracy'
    )
parser.add_argument('--noise',
    default = '0.1',
    help = 'limb darkening noise (%%)'
    )
parser.add_argument('--etb',
    default = '0.5',
    help = 'brightness temperature uncertainty (K)'
    )
parser.add_argument('--tref',
    default = '300.',
    help = 'reference temperature (K)'
    )
parser.add_argument('--eld',
    default = '0.1',
    help = 'limb darkening uncertainty (%%)'
    )
parser.add_argument('-d',
    action = 'store_true',
    default = False,
    help = 'fit differential'
    )
parser.add_argument('--case',
    choices = ['synthetic', 'global', 'zonal', 'grs',
               'PJ1', 'PJ3', 'PJ4', 'PJ5', 'PJ6', 'PJ7'],
    default = 'synthetic',
    )
args = vars(parser.parse_args())

def write_observation_file(outname, tb1, etb, ld1, eld, 
    ecal = 0., tref = 0., tb0 = None, ld0 = None):
  nwave = len(tb1)
  ndim = nwave*2
  cov = zeros((ndim, ndim))
# setup covariance matrix
  if tb0 != None :
    tb1 -= tb0
    for i in range(nwave):
      cov[i*2, i*2] = etb[i]*etb[i]
  else:   
    for i in range(nwave):
      cov[i*2, i*2] = etb[i]*etb[i] + (tb1[i]-tref)*(tb1[i]-tref)*ecal*ecal

  if ld0 != None :
    ld1 -= ld0
    for i in range(nwave):
      cov[i*2+1, i*2+1] = eld[i]*eld[i] + 0.1*0.1
  else:   
    for i in range(nwave):
      cov[i*2+1, i*2+1] = eld[i]*eld[i] + 0.1*0.1
  icov = linalg.inv(cov)

  with open('%s.obs' % outname, 'w') as file:
    file.write('# Observation file for case %s\n' % outname)
    file.write('%-d\n' % ndim)

    # write brightness temperature
    for i in range(nwave):
      file.write('%-.2f\n' % tb1[i])
      file.write('%-.3f\n' % ld1[i])

    # write inverse covariance matrix
    for i in range(ndim):
      for j in range(ndim):
        file.write('%-12.4g' % icov[i,j])
      file.write('\n')
  print('Observation file written to %s.obs' % outname)

def read_synthetic_data(fname):
# read emission angles
  with open('%s.out' % fname, 'r') as file:
    file.readline()
    line = file.readline()
  angles = array(list(map(float, line.split()[3:])))
  angles = list(map(int, angles))
  i45 = angles.index(45) + 1

# read brightness temperatures
  data = genfromtxt('%s.out' % fname)
  rows, cols = data.shape
  data = data.reshape(4, rows//4, cols)
  nwave = len(data[0])

# baseline brightness temperature and limb darkening
  tb0 = data[3,:,1]
  ld0 = (data[3,:,1] - data[3,:,i45])/data[3,:,1]*100.
  print("Baseline brightness temperature")
  print(tb0)
  print("Baseline limb darkening")
  print(ld0)

# new brightness temperature and limb darkening
  tb1 = data[0,:,1]
  ld1 = (data[0,:,1] - data[0,:,i45])/data[0,:,1]*100.
  print("New brightness temperature")
  print(tb1)
  print("New limb darkening")
  print(ld1)
  return tb1, ld1, tb0, ld0

def read_mwr_data():
  data = genfromtxt('data/mwr_zonal_average.txt', max_rows = 12)
  tb0 = data[:6,:]
  ld0 = data[6:,:]

  data = genfromtxt('data/mwr_zonal_average.txt', skip_header = 21)
  nrows, ncols = data.shape
  data = data.reshape(2, 6, nrows//12, -1)
  tb1 = data[0,:,:,1:]
  ld1 = data[1,:,:,1:]
  lat = data[0,0,:,0]

  return tb1, ld1, tb0, ld0, lat

if __name__ == '__main__':
  if args['case'] == 'synthetic':
    if args['input'] == 'none':
      raise ValueError('input file unspecified for synthetic inversion')
    tb1, ld1, tb0, ld0 = read_synthetic_data(args['input'])
  else:
    tb1, ld1, tb0, ld0, lat = read_mwr_data()
    # increase limb darkening uncertainty by a factor of 3
    # ld0[:,1] *= 3.
    # increase CH1 limb darkening uncertainty
    ld0[0,0] += 3.5
    ld0[0,1] *= 30.
    tb0[0,0] += 30.
    tb0[0,1] *= 10.

  zonal_str = 'mwr_PJ1-12_zonal'

# fit differential
  if args['case'] == 'synthetic':
    etb = float(args['etb'])*ones(len(tb1))
    eld = float(args['eld'])*ones(len(tb1))
    if args['d']:
      write_observation_file(args['input'],
          tb1, etb, ld1, eld,
          tref = float(args['tref']), tb0 = tb0, ld0 = ld0)
    else:
      write_observation_file(args['input'], 
          tb1, etb, ld1, eld,
          ecal = float(args['ecal']), tref = float(args['tref']))
  elif args['case'] == 'global':
    write_observation_file('mwr_PJ1-12_%s' % args['case'],
        tb0[:,0], tb0[:,1], ld0[:,0], ld0[:,1],
        ecal = float(args['ecal']), tref = float(args['tref']))
  elif args['case'] == 'zonal':
    nlat = len(lat)
    with open('%s_lats.txt' % zonal_str, 'w') as file:
      for i in range(nlat):
        file.write('%8.2f\n' % lat[i])
      file.write('\n')
    print('Latitude file written to %s_lats.txt' % zonal_str)
    for i in range(nlat):
      write_observation_file('%s_%.2f' % (zonal_str, lat[i]),
          tb1[:,i,0], tb1[:,i,1], ld1[:,i,0], ld1[:,i,1],
          ecal = 0.)
  elif args['case'] == 'grs':
    pass
  else :
    raise ValueError('Unrecognized case %' % args['case'])
