#! /usr/bin/env python3
from pylab import *
import argparse, os
from snapy.harp.radio_arguments import *
from snapy.harp.radio_model import create_input
from snapy.planet_gravity import saturn_gravity

def write_observation_file(outname, tb, stb, mu, smu,
    tb0 = None, ecal = 0.02):
  nwave = len(tb)
  ndim = nwave
  cov = zeros((ndim, ndim))
# setup covariance matrix
  if tb0 != None :
    tb -= tb0
    for i in range(nwave):
      cov[i,i] = (stb[i])**2 + (tb[i]*smu[i]/mu[i])**2
  else:   
    for i in range(nwave):
      cov[i,i] = (stb[i])**2 + (tb[i]*ecal)**2 + (tb[i]*smu[i]/mu[i])**2

  icov = linalg.inv(cov)

  with open('%s.obs' % outname, 'w') as file:
    file.write('# Observation file for case %s\n' % outname)
    file.write('%-d\n' % ndim)

    # write brightness temperature
    for i in range(nwave):
      file.write('%-.2f\n' % tb[i])

    # write inverse covariance matrix
    for i in range(ndim):
      for j in range(ndim):
        file.write('%-12.4g' % icov[i,j])
      file.write('\n')
  print('Observation file written to %s.obs' % outname)

if __name__ == '__main__':
  data = genfromtxt('./dat/vla_saturn_csxukq.txt')
  rows, cols = data.shape
  data = data.reshape((4,-1,cols))
  tb = data[0,:,3:]
  stb = data[1,:,3:]
  mu = data[2,:,3:]
  smu = data[3,:,3:]
  lat1 = data[0,:,1]
  lat2 = data[0,:,2]

  # band correction factor is the disk temperature ratio
  # between B&A and C&B bands
  iX = 2
  tb[:,iX] *= 1.03
  iU = 3
  tb[:,iU] *= 1.06

  # create observation file
  glat = float(args['glat'])
  ilat = searchsorted(lat1, glat)
  if args['output'] == '':
    outname = args['input'][:-4]
  else:
    outname = args['input'][:-4] + '-' + args['output']
  write_observation_file(outname, tb[ilat], stb[ilat], mu[ilat], smu[ilat],
    tb0 = None, ecal = 0.02)

  # create model input file from template file
  args['obs'] = outname
  args['grav'] = str(saturn_gravity(glat, norm = 'graphic'))

  for j in range(len(mu[ilat])):
    args['b%sdir' % (j+1,)] = '%.2f' % (arccos(mu[ilat,j])/pi*180.,)

  inpfile = create_input(args['input'], args)

