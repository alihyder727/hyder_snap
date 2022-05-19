#! /usr/bin/env python3
from pylab import *
from common import read_info
from astropy.io import fits
import argparse, warnings

def lat_lon(x,y,ob_lon,ob_lat,np_ang,req,rpol):
  '''Find latitude and longitude on planet given x,y pixel locations and
  planet equatorial and polar radius'''
  np_ang = -np_ang
  x1 = cos(radians(np_ang))*x - sin(radians(np_ang))*y
  y1 = sin(radians(np_ang))*x + cos(radians(np_ang))*y
  olrad = radians(ob_lat)
  
  #set up quadratic equation for ellipsoid
  r2 = (req/rpol)**2
  a = 1 + r2*(tan(olrad))**2 #second order
  b = 2*y1*r2*sin(olrad) / (cos(olrad)**2) #first order
  c = x1**2 + r2*y1**2 / (cos(olrad))**2 - req**2 #constant


  radical = b**2 - 4*a*c
  #will equal nan outside planet since radical < 0
  with warnings.catch_warnings():
    warnings.simplefilter("ignore") #suppresses error for taking sqrt nan
    x3s1=(-b+sqrt(radical))/(2*a)
    x3s2=(-b-sqrt(radical))/(2*a)
  z3s1=(y1+x3s1*sin(olrad))/cos(olrad)
  z3s2=(y1+x3s2*sin(olrad))/cos(olrad)
  odotr1=x3s1*cos(olrad)+z3s1*sin(olrad)
  odotr2=x3s2*cos(olrad)+z3s2*sin(olrad)
  #the two solutions are front and rear intersections with planet
  #only want front intersection
  
  #tricky way of putting all the positive solutions into one array
  with warnings.catch_warnings():
    warnings.simplefilter("ignore") #suppresses error for taking < nan
    odotr2[odotr2 < 0] = nan
    x3s2[odotr2 < 0] = nan
    z3s2[odotr2 < 0] = nan
    odotr1[odotr1 < 0] = odotr2[odotr1 < 0]
    x3s1[odotr1 < 0] = x3s2[odotr1 < 0]
    z3s1[odotr1 < 0] = z3s2[odotr1 < 0]
  
  odotr,x3,z3 = odotr1,x3s1,z3s1
  y3 = x1
  r = sqrt(x3**2 + y3**2 + z3**2)
  
  #lon_w = degrees(arctan(y3/x3)) + ob_lon
  lon_w = degrees(arctan2(x3,y3)-pi/2) + ob_lon 
  with warnings.catch_warnings():
    warnings.simplefilter("ignore") #suppresses error for taking < nan
    lon_w[lon_w < 0] += 360
    lon_w = lon_w%360
  lat_c = degrees(arcsin(z3/r))
  lat_g = degrees(arctan(r2*tan(radians(lat_c))))
  #plt.imshow(lon_w, origin = 'lower left')
  #plt.show()
  return lat_g, lat_c, lon_w

def surface_normal(lat_g, lon_w, ob_lon):
  '''Returns the normal vector to the surface of the planet.
  Take dot product with sub-obs or sub-sun vector to find cosine of emission angle'''
  nx = cos(radians(lat_g))*cos(radians(lon_w-ob_lon))
  ny = cos(radians(lat_g))*sin(radians(lon_w-ob_lon))
  nz = sin(radians(lat_g))
  return array([nx,ny,nz])

def emission_angle(ob_lat, surf_n):
  '''Return the cosine of the emission angle of surface wrt observer'''
  ob = array([cos(radians(ob_lat)),0,sin(radians(ob_lat))])
  return dot(surf_n.T, ob).T

if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument('-f', '--fits',
    default = 'tbd',
    help = 'base name of the fits file without extension'
    )
  args = vars(parser.parse_args())

  print('Running extract_latlon ...')

  info = read_info('../dat/%s.info' % args['fits'])

  xx = arange(info['xsize']) - info['xcenter']
  yy = arange(info['ysize']) - info['ycenter']
  x,y = meshgrid(xx,yy)
  if info['target'].lower() == 'saturn':
    re = 60268/info['scale1']
    rp = 54364/info['scale2']
  else:
    re = info['re']
    rp = info['rp']
  #print(re,rp)
  #exit(1)
  lat_g, lat_c, lon_w = lat_lon(x,y,
    info['oblon'],info['oblat'],0.,re,rp)
  surf_n = surface_normal(lat_g, lon_w, info['oblon'])
  mu = emission_angle(info['oblat'], surf_n)

  tb = fits.open(info['fits'])[0].data[0,0,:,:]
  tb = tb[~isnan(lat_c)]
  mu = mu[~isnan(mu)]
  latg = lat_g[~isnan(lat_g)]
  latc = lat_c[~isnan(lat_c)]
  lonw = lon_w[~isnan(lon_w)]

  with open('../dat/%s.dat' % args['fits'], 'w') as file:
    file.write('#%4s%12s%12s%12s%12s%12s\n' % 
      ('IDX', 'LATG', 'LATC', 'LONW', 'MU', 'TB'))
    for i in range(len(latg)):
      file.write('%5d%12.2f%12.2f%12.2f%12.4f%12.2f\n' %
        (i, latg[i], latc[i], lonw[i], mu[i], tb[i]))

  with open('../dat/%s.info' % args['fits'], 'a') as file:
    file.write('# TB \n')
    file.write('TB_MAX              = %.2f\n' % tb.max())
    file.write('TB_MIN              = %.2f\n' % tb.min())
    file.write('ANGLE_MAX           = %.2f\n' % degrees(arccos(mu.max())))
    file.write('ANGLE_MIN           = %.2f\n' % degrees(arccos(mu.min())))
    file.write('LAT_MAX (GRAPHIC)   = %.2f\n' % latg.max())
    file.write('LAT_MIN (GRAPHIC)   = %.2f\n' % latg.min())
    file.write('LON_MAX (W)         = %.2f\n' % lonw.max())
    file.write('LON_MIN (W)         = %.2f\n' % lonw.min())
