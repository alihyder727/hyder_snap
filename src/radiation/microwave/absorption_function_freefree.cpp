/**************************************************************************

                    =====  JUNO codes  =====

   Copyright (C) 2019 California Institute of Technology (Caltech)
   All rights reserved.

   Written by Fabiano Oyafuso (fabiano@jpl.nasa.gov)
   Adapted by Cheng Li (cli@gps.caltech.edu) to SNAP structure

**************************************************************************/

// C/C++ headers
#include <cmath>

// cli, 190801
inline double SQUARE(double x) { return x*x; }
inline double CUBE(double x) { return x*x*x; }
static double const pi           = 3.14159265358979;
// end

double absorption_coefficient_freefree_Reference(double freq, double P, double T)
   // this is very crude
{
   int const Z = 1;
   double const fine_struct = 1.0/137.036;
   double const mc2 = 0.511e6; // eV
   double const kT = (T / 300) * .0256; // eV
   double const hbarc = 1973e-8; // eV-cm
   double sigma_T = 6.65e-25;  // cm^2
   //double g_ff = 1; // gaunt factor
   double hbar_omega=4.13566733e-15 * freq * 1e9;  // eV (.6GHz) -- check that freq is in GHz

   double debroglie_thermal = sqrt(2.0*pi*SQUARE(hbarc)/(mc2*kT));
   double E_ionization = 5.986; // eV
   double n_Al = 1.6e-6 * (P*1e5) / (kT * 1.609e-19) * 1e-6;
   
   double ni = sqrt(2 * n_Al / pow(debroglie_thermal,3) * exp(-E_ionization/kT));

   double kludge = 0.5;
   double alpha_ff = kludge * sqrt(32*CUBE(pi)/3.0) * (SQUARE(Z)*fine_struct) * sqrt(mc2/kT) * 
      (ni*CUBE(hbarc/hbar_omega)) * ni*sigma_T * (1.0-exp(-hbar_omega/kT));
  return alpha_ff;
}

// from C. Li (in units of 1/m)
double absorption_coefficient_freefree_Chengli(double freq, double P, double T)
{
   double ion_opacity_scale = 1.0;
   double ion_opacity_cutoff = 1600.;

   double alpha_ff;
   if (T < ion_opacity_cutoff)
      alpha_ff=0;
   else {
      alpha_ff=ion_opacity_scale*exp(-20.957 + 9.51112E-3*(1200.+T-ion_opacity_cutoff));
      // convert to cm^-1
      alpha_ff *= 1e-2;
   }
   return alpha_ff;
}
