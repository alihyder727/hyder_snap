from numpy import sin,sqrt,pi,tan,arctan

def deg2rad(ang):
  return ang*pi/180.

def rad2deg(rad):
  return rad*180./pi

def graphic2centric(lat, rerp):
  if abs(lat) == 90.:
    return lat
  else:
    return rad2deg(arctan(tan(deg2rad(lat))/rerp**2))

def centric2graphic(lat, rerp):
  if abs(lat) == 90.:
    return lat
  else:
    return rad2deg(arctan(tan(deg2rad(lat))*rerp**2))

def jupiter_gravity(lat, norm = 'centric'):
  if norm == 'centric':
    pclat = lat
  elif norm == 'graphic': 
    pclat = graphic2centric(lat, 1.07)
  else: raise ValueError("norm not recognized")

  grav = 23.3;
  S = sin(pclat * pi / 180.);
  SS = S*S;
  CS = S*sqrt(1 - SS);
  GR = - grav + SS*(-4.26594 + SS*(0.47685 + SS*(-0.100513 + SS*(0.0237067 - 0.00305515*SS))));
  GTH = CS*(-3.42313 +  SS*(0.119119 + SS*(0.00533106 + SS*(-0.00647658 + SS*0.000785945))));
  return sqrt(GR*GR+GTH*GTH);

def saturn_gravity(lat, norm = 'centric'):
  if norm == 'centric':
    pclat = lat
  elif norm == 'graphic':
    pclat = graphic2centric(lat, 1.11)
  else: raise ValueError("norm not recognized")

  grav = 9.06656
  S = sin(pclat * pi / 180.);
  SS = S * S;
  CS = S*sqrt(1 - SS);
  GR = - grav + SS*(-3.59253 + SS*(0.704538 + SS*(-0.260158 + SS*(0.0923098 - SS*0.0166287))));
  GTH = CS*(-2.25384 + SS*(.152112 + SS*(-.0102391 + SS*(-.00714765 + SS*.000865634))));
  return sqrt(GR*GR+GTH*GTH);
