#! /usr/bin/env python3
import argparse
from netCDF4 import Dataset

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input',
    required = True,
    help = 'netcdf file for ck bands'
    )
parser.add_argument('--var',
    default= 'kcoeff',
    help = 'absorbing molecules'
    )
args = vars(parser.parse_args())


data = Dataset(args['input'], 'r')
ng = data['ngauss'][:]
edges = data['bin_edges'][:]
weights = ['%.4f' % x for x in data['weights'][:]]
gpoints = ['%.4f' % x for x in data['samples'][:]]

anames = args['var'].split(',')
fname = args['input'][:-3]

with open('%s_bands.txt' % fname, 'w') as file:
  file.write('<radiation>\n')
  for i in range(1, len(edges)):
    file.write('b%d = %.2f %.2f %d\n' % (i, edges[i-1],  edges[i], ng))
    file.write('b%d' % i + '.gpoints = %s\n' % ','.join(gpoints))
    file.write('b%d' % i + '.weights = %s\n' % ','.join(weights))
    file.write('b%d' % i + '.absorbers = %s\n' % ' '.join(['ck-%s-%d' % (x,i-1) for x in anames]))
    for aname in anames:
      file.write('b%d.ck-%s-%d = %s\n' % (i, aname, i-1, args['input']))
    file.write('\n')

print('ck bands written to %s_bands.txt' % args['input'][:-3])
