#! /usr/bin/env python3
from pylab import *
import argparse

parser = argparse.ArgumentParser(
    formatter_class = argparse.ArgumentDefaultsHelpFormatter
    )
parser.add_argument('-i', '--input',
  required = True,
  help = 'brightness temperature output file (required)'
  )
parser.add_argument('--ecal',
  default = '0.03',
  help = 'calibration accuracy (%%)'
  )
parser.add_argument('--etb',
  default = '1.',
  help = 'brightness temperature uncertainty, (K)'
  )
parser.add_argument('--eld',
  default = '0.5',
  help = 'limb darkening uncertainty (%%)'
  )
parser.add_argument('-d',
  action = 'store_true',
  default = False,
  help = 'fit differential'
  )
args = vars(parser.parse_args())
etb = float(args['etb'])
eld = float(args['eld'])
ecal = float(args['ecal'])

# read emission angles
with open('%s.out' % args['input'], 'r') as file:
  file.readline()
  line = file.readline()
angles = array(list(map(float, line.split()[3:])))
angles = list(map(int, angles))
i45 = angles.index(45) + 1

# read brightness temperatures
data = genfromtxt('%s.out' % args['input'])
rows, cols = data.shape
data = data.reshape(4, rows//4, cols)
nwave = len(data[0])

print(data)
print(i45)

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

#ndim = nwave*2
ndim = nwave
cov = zeros((ndim, ndim))
# setup covariance matrix
print("New brightness temerature coeffs:")
if args['d']: # fit differential
  tb1 -= tb0
  ld1 -= ld0
  for i in range(nwave):
    # tb
    cov[i, i] = etb*etb
    # ld
    #cov[i*2+1, i*2+1] = eld*eld
  outname = '%s.dobs' % args['input']
else:   
  for i in range(nwave):
    # tb
    cov[i*2,j*2] = etb*etb + (tb1[i]-300.)*(tb1[i]-300.)*ecal*ecal
    # ld
    #cov[i*2+1, j*2+1] = eld*eld
  outname = '%s.obs' % args['input']
icov = linalg.inv(cov)

with open(outname, 'w') as file:
  file.write('# Observation file for case %s\n' % args['input'])
  file.write('%-d\n' % ndim)

  # write brightness temperature
  for i in range(nwave):
    file.write('%-.2f\n' % tb1[i])
    #file.write('%-.3f\n' % ld1[i])

  # write inverse covariance matrix
  for i in range(ndim):
    for j in range(ndim):
      file.write('%-8.4g' % icov[i,j])
    file.write('\n')
print('Observation file written to %s' % outname)
