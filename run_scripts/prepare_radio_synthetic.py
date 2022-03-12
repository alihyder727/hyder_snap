#! /usr/bin/env python3
from pylab import *
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input',
  required = True,
  help = 'brightness temperature output file (required)'
  )
parser.add_argument('--acc',
  default = '0.02',
  help = 'calibration accuracy'
  )
parser.add_argument('--noise',
  default = '0.5',
  help = 'precision noise'
  )
parser.add_argument('-d',
  action = 'store_true',
  default = False,
  help = 'fit differential'
  )
args = vars(parser.parse_args())

# read emission angles
with open('%s.out' % args['input'], 'r') as file:
  file.readline()
  line = file.readline()
angles = array(list(map(float, line.split()[3:])))
angles = angles[where(angles <= 45.)]
mu = cos(angles/180.*pi)
nangle = len(angles)

# read brightness temperatures
data = genfromtxt('%s.out' % args['input'])
rows, cols = data.shape
data = data.reshape(4, rows//4, cols)
nwave = len(data[0])

acc_cal = float(args['acc'])
noise = float(args['noise'])

# fit quadratic coefficient using least square
B = zeros((nangle, 3))

for i in range(nangle):
  B[i,0] = 1.
  B[i,1] = 1. - mu[i]
  B[i,2] = (1. - mu[i])*(1. - mu[i])

# baseline brightness temperature coeff
coeffs = zeros((nwave, 3))
for i in range(nwave):
  Tb = data[3,i,1:1+nangle]
  coeffs[i] = linalg.solve(dot(B.T,B), dot(B.T,Tb))
print("Baseline brightness temperature coeffs:")
print(coeffs)

ndim = nwave*3
cov = zeros((ndim, ndim))
# setup covariance matrix
print("New brightness temerature coeffs:")
if args['d']: # fit differential
  data = genfromtxt('%s.out' % args['input'])
  data = data.reshape(4, rows//4, cols)
  for i in range(nwave):
    Tb = data[0,i,1:1+nangle]
    coeff = linalg.solve(dot(B.T,B), dot(B.T,Tb))
    print(coeff)
    coeffs[i] = coeff - coeffs[i]
  for i in range(ndim):
    for j in range(ndim):
      if i == j:
        cov[i,j] = noise*noise
      else:
        cov[i,j] = 0.
  outname = '%s.dobs' % args['input']
else:   
  for i in range(ndim):
    for j in range(ndim):
      if i == j:
        if i%3 == 0:
          cov[i,j] = (acc_cal*coeffs[i//3,0])**2
        else:
          cov[i,j] = noise**2
      else:
        cov[i,j] = 0.
  outname = '%s.obs' % args['input']
icov = linalg.inv(cov)

with open(outname, 'w') as file:
  file.write('# Observation file for case %s\n' % args['input'])
  file.write('%-d\n' % ndim)
  # write brightness temperature
  for i in range(nwave):
    for j in range(3):
      file.write('%-.2f\n' % coeffs[i,j])
  # write inverse covariance matrix
  for i in range(ndim):
    for j in range(ndim):
      file.write('%-8.2g' % icov[i,j])
    file.write('\n')
print('Observation file written to %s' % outname)
