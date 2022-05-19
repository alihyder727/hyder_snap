#! /usr/bin/env python3
import matplotlib
matplotlib.use('Agg')
from pylab import *
from astropy.io import fits
from matplotlib.patches import Ellipse
from common import read_info
import argparse, os

## Fit an ellipse, based on the reference 
##  (*) Halir R., Flusser J. 'Numerically Stable Direct Least Squares Fitting of Ellipses'
##  (**) Weisstein, Eric W. "Ellipse." From MathWorld--A Wolfram Web Resource. http://mathworld.wolfram.com/Ellipse.html
## and git repo: https://github.com/bdhammel/least-squares-ellipse-fitting

def fit_ellipse(x, y):
  D1 = vstack([x**2, x*y, y**2]).T
  D2 = vstack([x, y, ones(len(x))]).T

  S1 = dot(D1.T, D1)
  S2 = dot(D1.T, D2)
  S3 = dot(D2.T, D2)

  T = -dot(linalg.inv(S3), S2.T)
  M = S1 + dot(S2, T)
  M = vstack([M[2,:]/2., -M[1,:], M[0,:]/2.])

  eigval, eigvec = linalg.eig(M)

  cond = 4.*eigvec[0,:]*eigvec[2,:] - eigvec[1,:]**2
  a1 = eigvec[:, nonzero(cond > 0)[0]]
  a2 = dot(T, a1)
  a = vstack([a1, a2])
  D = hstack([D1, D2])

  gfit = linalg.norm(dot(D, a), 1)/len(D1)

  return a.flatten(), gfit

def get_parameters(coeff):
  # Eigenvectors are the coefficients of an ellipse in general form
  # a*x^2 + 2*b*x*y + c*y^2 + 2*d*x + 2*f*y + g = 0
  # [eqn. 15) from (**) or (***)
  a = coeff[0]
  b = coeff[1] / 2.
  c = coeff[2]
  d = coeff[3] / 2.
  f = coeff[4] / 2.
  g = coeff[5]

  # Finding center of ellipse [eqn.19 and 20] from (**)
  x0 = (c*d - b*f) / (b**2. - a*c)
  y0 = (a*f - b*d) / (b**2. - a*c)

  # Find the semi-axes lengths [eqn. 21 and 22] from (**)
  numerator = 2 * (a*f**2 + c*d**2 + g*b**2 - 2*b*d*f - a*c*g)
  denominator1 = (b * b - a * c) * (
      (c - a) * np.sqrt(1 + 4*b*b / ((a - c)*(a - c))) - (c + a)
  )
  denominator2 = (b*b - a*c) * (
      (a - c) * np.sqrt(1 + 4*b*b / ((a - c) * (a - c))) - (c + a)
  )
  re = np.sqrt(numerator / denominator1)
  rp = np.sqrt(numerator / denominator2)

  # Angle of counterclockwise rotation of major-axis of ellipse to x-axis
  # [eqn. 23] from (**) or [eqn. 26] from (***).
  phi = .5 * np.arctan((2.*b) / (a - c))

  return x0, y0, re, rp, phi
  

if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument('-d', '--dir',
    default = '.',
    help = 'directory of the fits file'
    )
  parser.add_argument('-f', '--fits',
    default = 'tbd',
    help = 'base name of the fits file without extension'
    )
  parser.add_argument('-v', '--value',
    default = '3.',
    help = 'background pixel values'
    )
  args = vars(parser.parse_args())

  print('Running detect_planet ...')

  # info = read_info('../dat/%s.info' % args['fits'])
  hdu = fits.open('%s/%s.fits' % (args['dir'], args['fits']))[0]
  if (args['fits'][-3:] == 'msk'):
    data = hdu.data[:,:]
  else :
    data = hdu.data[0,0,:,:]

  fig, axs = subplots(1, 3, figsize = (12, 4))

  # plot original image
  ax = axs[0]
  ax.imshow(data, origin = 'lower', cmap = 'Greys_r')

  # plot fitted image
  ax = axs[1]
  ax.imshow(data, origin = 'lower', cmap = 'Greys_r')
  value = float(args['value'])
  cs = ax.contour(data, [value], colors = 'C1', linewidths = 1)

  # plot histgram
  ax = axs[2]
  bins, edges, patches = ax.hist(data.flatten(), 50)
  ax.plot([value, value], [bins.min()/2, bins.max()/2], 'r')
  ax.set_yscale('log')

  seg = cs.collections[0].get_segments()
  seglen = array(map(len, seg))
  ix = seglen.argmax()
  vert = seg[ix]

  coeff, gfit = fit_ellipse(vert[:,0], vert[:,1])
  x0, y0, re, rp, phi = get_parameters(coeff)
  #print(re, rp)

  ellipse = Ellipse((x0, y0), width = re*2, height = rp*2, angle = phi/pi*180.,
    fill = False, color = 'g', ls = ':', lw = 1)
  ax.add_patch(ellipse)

  with open('../dat/%s.info' % args['fits'], 'a') as file:
    #px = data.argmax()
    #ix, iy = px//data.shape[1], px%data.shape[1]
    file.write('# PLANET DETECTION\n')
    file.write('PIXELS (PX)       = %dx%d\n' % data.shape)
    file.write('CENTER (PX)       = %.2f,%.2f\n' % (x0, y0))
    file.write('RE (PX)           = %.2f\n' % re)
    file.write('RP (PX)           = %.2f\n' % rp)
    file.write('INC (DEG)         = %.2f\n' % (phi/pi*180.,))
    file.write('FITTING           = %.2f\n' % gfit)

  os.system('mkdir -p ../img')
  savefig('../img/%s-ellipse.png' % args['fits'])
