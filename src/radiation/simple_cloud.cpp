// C/C++ headers
#include <iostream>
#include <algorithm>
#include <stdexcept>
#include <cstring>
#include <cassert>  // assert

// Athena++ headers
//#include "../math_funcs.hpp"
#include "absorber.hpp"
#include "water_cloud.hpp"
#include "radiation_utils.hpp"  // getPhaseMomentum


// For grey cloud
Real SimpleCloud::AbsorptionCoefficient(Real wave, Real const chem[]) const
{
  //Real result= 1.;
  //Real result= 1.E3*exp(- pow( ((log(chem[IPR])-log(5.E6))/2.), 2)) ;
  //Real result= 1.E1*exp(- pow( ((log(chem[IPR])-log(5.E5))/1.), 2)) ;
  Real csize = 1.E0*1.e-6; // one micron size particle
  Real qext = 1.E0;
  Real crho = 5.E3;
  //std::cout << chem[0] << " " << chem[6] << " " << chem[7] << " " << chem[9] << std::endl;
  return  chem[imol_]*qext/(4./3.*csize*crho);     // -> 1/m
  //return chem[imol_]*result*mixr_;     // -> 1/m
}

Real SimpleCloud::SingleScateringAlbedo(Real wave, Real const chem[]) const
{
  // ssalb
  Real ww = 0.9;

  if (chem[IPR] > 1)
    return ww;
  else
    return 0.;
}


void SimpleCloud::PhaseMomentum(Real wave, Real const chem[], Real *pp, int np) const
{
  Real gg=0.9;

  if (chem[IPR] > 1)
    getPhaseMomentum(0, gg, np, pp);  // 0 for HENYEY_GREENSTEIN
  else
    getPhaseMomentum(0, 0.0, np, pp);  // 0 for HENYEY_GREENSTEIN
}

