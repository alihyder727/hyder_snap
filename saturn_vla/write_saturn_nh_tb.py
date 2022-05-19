#! /usr/bin/env python3
from pylab import *

# S band
data = genfromtxt('../dat/sat-s-witb.dat')
latbins = arange(0., 90., 4.)
latg = data[:,1]
lonw = data[:,3]
mu = data[:,4]
tb = data[:,5]
fname = '../dat/vla_saturn_s.txt'
with open(fname, 'w') as file:
  file.write('# VLA smeared map of saturn, S band, 13.33 - 15.04 cm, 2 - 2.25 GHz\n')
  file.write('#%9s%10s%10s%10s%10s%10s%10s%10s\n' % ('LATG1', 'LATG2', 'TB', 'sTB', 'MU',
      'sMU', 'MUMAX', 'N'))
  i, j = 0, 1
  while i < j:
    if j == len(latbins) - 1:
      break
    ix = where((latg > latbins[i]) & (latg < latbins[j]))[0]
    try :
      lonmin, lonmax = min(lonw[ix]), max(lonw[ix])
    except ValueError :
      j += 1
      continue
    #print('latbin = %.1f - %.1f' % (latbins[i], latbins[j]))
    #print('min lon = %.1f' % lonmin)
    #print('max lon = %.1f' % lonmax)
    #print('num lon = %d' % len(lonw[ix]))
    # Take 1/4 of the longitudes
    lon1 = lonmin + (lonmax - lonmin)/3.
    lon2 = lonmax - (lonmax - lonmin)/3.
    iy = where((latg > latbins[i]) & (latg < latbins[j]) & (lonw > lon1) & (lonw < lon2))[0]
    if len(lonw[iy]) < 4:
      j += 1
      continue
    #print('num lon center = %d' % len(lonw[iy]))
    mu_avg = mean(mu[iy])
    mu_max = max(mu[iy])
    tb_avg = mean(tb[iy])
    mu_std = std(mu[iy])
    tb_std = std(tb[iy])
    #print('tb mean = %.1f, mu mean = %.2f' % (tb_avg, mu_avg))
    #print('tb std = %.1f, mu std = %.2f' % (tb_std, mu_std))
    file.write('%10.1f%10.1f%10.2f%10.2f%10.2f%10.3f%10.3f%10d\n' %
        (latbins[i], latbins[j], tb_avg, tb_std, mu_avg, mu_std, mu_max, len(tb[iy])))
    i = j
    j = i + 1
print('output file written to %s' % fname)

# C band
data = genfromtxt('../dat/sat-c-witb.dat')
latbins = arange(0., 91., 1.)
latg = data[:,1]
lonw = data[:,3]
mu = data[:,4]
tb = data[:,5]
fname = '../dat/vla_saturn_c.txt'
with open(fname, 'w') as file:
  file.write('# VLA smeared map of saturn, C band, 7.08 - 7.53 cm, 4. - 4.25 GHz\n')
  file.write('#%9s%10s%10s%10s%10s%10s%10s\n' % ('LATG1', 'LATG2', 'TB', 'sTB', 'MU', 'sMU', 'N'))
  i, j = 0, 1
  while i < j:
    if j == len(latbins) - 1:
      break
    ix = where((latg > latbins[i]) & (latg < latbins[j]))[0]
    lonmin, lonmax = min(lonw[ix]), max(lonw[ix])
    #print('latbin = %.1f - %.1f' % (latbins[i], latbins[j]))
    #print('min lon = %.1f' % lonmin)
    #print('max lon = %.1f' % lonmax)
    #print('num lon = %d' % len(lonw[ix]))
    # Take 1/4 of the longitudes
    lon1 = lonmin + (lonmax - lonmin)/3.
    lon2 = lonmax - (lonmax - lonmin)/3.
    iy = where((latg > latbins[i]) & (latg < latbins[j]) & (lonw > lon1) & (lonw < lon2))[0]
    if len(lonw[iy]) < 4:
      j += 1
      continue
    #print('num lon center = %d' % len(lonw[iy]))
    mu_avg = mean(mu[iy])
    tb_avg = mean(tb[iy])
    mu_std = std(mu[iy])
    tb_std = std(tb[iy])
    #print('tb mean = %.1f, mu mean = %.2f' % (tb_avg, mu_avg))
    #print('tb std = %.1f, mu std = %.2f' % (tb_std, mu_std))
    file.write('%10.1f%10.1f%10.2f%10.2f%10.2f%10.3f%10d\n' %
        (latbins[i], latbins[j], tb_avg, tb_std, mu_avg, mu_std, len(tb[iy])))
    i = j
    j = i + 1
print('output file written to %s' % fname)

# X band
coeff = genfromtxt('../dat/limb_disk.txt', delimiter = ',')
coeff = coeff[2,1:-1]
data = genfromtxt('../dat/sat-x-ictb.dat')
latbins = arange(0., 90.2, 0.2)
latg = data[:,1]
lonw = data[:,3]
lonw[lonw < 180.] += 360.
mu = data[:,4]
tb = data[:,5]
fname = '../dat/vla_saturn_x.txt'
with open(fname, 'w') as file:
  file.write('# VLA smeared map of saturn, X band, 2.5 - 3.76 cm, 8 - 12 GHz\n')
  file.write('#%9s%10s%10s%10s%10s%10s%10s\n' % ('LATG1', 'LATG2', 'TB', 'sTB', 'MU', 'sMU', 'N'))
  i, j = 0, 1
  while i < j:
    if j == len(latbins) - 1:
      break
    ix = where((latg > latbins[i]) & (latg < latbins[j]))[0]
    lonmin, lonmax = min(lonw[ix]), max(lonw[ix])
    #print('latbin = %.1f - %.1f' % (latbins[i], latbins[j]))
    #print('min lon = %.1f' % lonmin)
    #print('max lon = %.1f' % lonmax)
    #print('num lon = %d' % len(lonw[ix]))
    # Take 1/4 of the longitudes
    lon1 = lonmin + (lonmax - lonmin)/3.
    lon2 = lonmax - (lonmax - lonmin)/3.
    iy = where((latg > latbins[i]) & (latg < latbins[j]) & (lonw > lon1) & (lonw < lon2))[0]
    #print(lonw[iy])
    if len(lonw[iy]) < 4:
      j += 1
      continue
    #print('num lon center = %d' % len(lonw[iy]))
    tb_disk = (coeff[0]*pow(mu[iy], coeff[1]) + coeff[5]*pow(mu[iy], coeff[6]))/2.
    mu_avg = mean(mu[iy])
    tb_avg = mean(tb[iy] + tb_disk)
    mu_std = std(mu[iy])
    tb_std = std(tb[iy] + tb_disk)
    #print('tb mean = %.1f, mu mean = %.2f' % (tb_avg, mu_avg))
    #print('tb std = %.1f, mu std = %.2f' % (tb_std, mu_std))
    file.write('%10.1f%10.1f%10.2f%10.2f%10.2f%10.3f%10d\n' %
        (latbins[i], latbins[j], tb_avg, tb_std, mu_avg, mu_std, len(tb[iy])))
    i = j
    j = i + 1
print('output file written to %s' % fname)

# U band
coeff = genfromtxt('../dat/limb_disk.txt', delimiter = ',')
coeff = coeff[4,1:-1]
data = genfromtxt('../dat/sat-u-ictb.dat')
latbins = arange(0., 90.2, 0.2)
latg = data[:,1]
lonw = data[:,3]
lonw[lonw < 100.] += 360.
mu = data[:,4]
tb = data[:,5]
fname = '../dat/vla_saturn_u.txt'
with open(fname, 'w') as file:
  file.write('# VLA smeared map of saturn, U band, 1.66 - 2.50 cm, 12 - 18 GHz\n')
  file.write('#%9s%10s%10s%10s%10s%10s%10s\n' % ('LATG1', 'LATG2', 'TB', 'sTB', 'MU', 'sMU', 'N'))
  i, j = 0, 1
  while i < j:
    if j == len(latbins) - 1:
      break
    ix = where((latg > latbins[i]) & (latg < latbins[j]))[0]
    try :
      lonmin, lonmax = min(lonw[ix]), max(lonw[ix])
    except ValueError :
      j += 1
      continue
    #print('latbin = %.1f - %.1f' % (latbins[i], latbins[j]))
    #print('min lon = %.1f' % lonmin)
    #print('max lon = %.1f' % lonmax)
    #print('num lon = %d' % len(lonw[ix]))
    # Take 1/4 of the longitudes
    lon1 = lonmin + (lonmax - lonmin)/3.
    lon2 = lonmax - (lonmax - lonmin)/3.
    iy = where((latg > latbins[i]) & (latg < latbins[j]) & (lonw > lon1) & (lonw < lon2))[0]
    #print(lonw[iy])
    if len(lonw[iy]) < 4:
      j += 1
      continue
    #print('num lon center = %d' % len(lonw[iy]))
    tb_disk = (coeff[0]*pow(mu[iy], coeff[1]) + coeff[5]*pow(mu[iy], coeff[6]))/2.
    mu_avg = mean(mu[iy])
    tb_avg = mean(tb[iy] + tb_disk)
    mu_std = std(mu[iy])
    tb_std = std(tb[iy] + tb_disk)
    #print('tb mean = %.1f, mu mean = %.2f' % (tb_avg, mu_avg))
    #print('tb std = %.1f, mu std = %.2f' % (tb_std, mu_std))
    file.write('%10.1f%10.1f%10.2f%10.2f%10.2f%10.3f%10d\n' %
        (latbins[i], latbins[j], tb_avg, tb_std, mu_avg, mu_std, len(tb[iy])))
    i = j
    j = i + 1
print('output file written to %s' % fname)

# K band
data = genfromtxt('../dat/sat-k-witb.dat')
latbins = arange(0., 90.2, 0.2)
latg = data[:,1]
lonw = data[:,3]
mu = data[:,4]
tb = data[:,5]
fname = '../dat/vla_saturn_k.txt'
with open(fname, 'w') as file:
  file.write('# VLA smeared map of saturn, K band, 1.15 - 1.67 cm, 18 - 26 GHz\n')
  file.write('#%9s%10s%10s%10s%10s%10s%10s\n' % ('LATG1', 'LATG2', 'TB', 'sTB', 'MU', 'sMU', 'N'))
  i, j = 0, 1
  while i < j:
    if j == len(latbins) - 1:
      break
    ix = where((latg > latbins[i]) & (latg < latbins[j]))[0]
    try :
      lonmin, lonmax = min(lonw[ix]), max(lonw[ix])
    except ValueError :
      j += 1
      continue
    #print('latbin = %.1f - %.1f' % (latbins[i], latbins[j]))
    #print('min lon = %.1f' % lonmin)
    #print('max lon = %.1f' % lonmax)
    #print('num lon = %d' % len(lonw[ix]))
    # Take 1/4 of the longitudes
    lon1 = lonmin + (lonmax - lonmin)/3.
    lon2 = lonmax - (lonmax - lonmin)/3.
    iy = where((latg > latbins[i]) & (latg < latbins[j]) & (lonw > lon1) & (lonw < lon2))[0]
    #print(lonw[iy])
    if len(lonw[iy]) < 8:
      j += 1
      continue
    #print('num lon center = %d' % len(lonw[iy]))
    mu_avg = mean(mu[iy])
    tb_avg = mean(tb[iy])
    mu_std = std(mu[iy])
    tb_std = std(tb[iy])
    #print('tb mean = %.1f, mu mean = %.2f' % (tb_avg, mu_avg))
    #print('tb std = %.1f, mu std = %.2f' % (tb_std, mu_std))
    file.write('%10.1f%10.1f%10.2f%10.2f%10.2f%10.3f%10d\n' %
        (latbins[i], latbins[j], tb_avg, tb_std, mu_avg, mu_std, len(tb[iy])))
    i = j
    j = i + 1
print('output file written to %s' % fname)

# Q band
data = genfromtxt('../dat/sat-q-witb.dat')
latbins = arange(0., 90.1, 0.1)
latg = data[:,1]
lonw = data[:,3]
lonw[lonw > 250.] -= 360.
mu = data[:,4]
tb = data[:,5]
fname = '../dat/vla_saturn_q.txt'
with open(fname, 'w') as file:
  file.write('# VLA smeared map of saturn, Q band, 0.62 - 0.75 cm, 40 - 48.4 GHz\n')
  file.write('#%9s%10s%10s%10s%10s%10s%10s\n' % ('LATG1', 'LATG2', 'TB', 'sTB', 'MU', 'sMU', 'N'))
  i, j = 0, 1
  while i < j:
    if j == len(latbins) - 1:
      break
    ix = where((latg > latbins[i]) & (latg < latbins[j]))[0]
    try :
      lonmin, lonmax = min(lonw[ix]), max(lonw[ix])
    except ValueError :
      j += 1
      continue
    #print('latbin = %.1f - %.1f' % (latbins[i], latbins[j]))
    #print('min lon = %.1f' % lonmin)
    #print('max lon = %.1f' % lonmax)
    #print('num lon = %d' % len(lonw[ix]))
    # Take 1/4 of the longitudes
    lon1 = lonmin + (lonmax - lonmin)/3.
    lon2 = lonmax - (lonmax - lonmin)/3.
    iy = where((latg > latbins[i]) & (latg < latbins[j]) & (lonw > lon1) & (lonw < lon2))[0]
    #print(lonw[iy])
    if len(lonw[iy]) < 4:
      j += 1
      continue
    #print('num lon center = %d' % len(lonw[iy]))
    mu_avg = mean(mu[iy])
    tb_avg = mean(tb[iy])
    mu_std = std(mu[iy])
    tb_std = std(tb[iy])
    #print('tb mean = %.1f, mu mean = %.2f' % (tb_avg, mu_avg))
    #print('tb std = %.1f, mu std = %.2f' % (tb_std, mu_std))
    file.write('%10.1f%10.1f%10.2f%10.2f%10.2f%10.3f%10d\n' %
        (latbins[i], latbins[j], tb_avg, tb_std, mu_avg, mu_std, len(tb[iy])))
    i = j
    j = i + 1
print('output file written to %s' % fname)
