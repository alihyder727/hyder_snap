#ifndef INTERPOLATION_H_
#define INTERPOLATION_H_

struct float_triplet;

int locate(double const *xx, double x, int n);
void interpn(double *val, double const *coor, double const *data, double const *axis, 
  int const *len, int ndim, int nval = 1);
double interp1(double x, double const *data, double const *axis, int len);

void spline(int n, float_triplet *table, double y1_bot, double y1_top);
int find_place_in_table(int n, float_triplet *table, double x, double *dx, int il = -1);
double splint(double xx, float_triplet *table, double dx);

void interpnf(double *val, double const *coor, double const *data, double const *axis, 
  int const *len, int ndim);

#endif
