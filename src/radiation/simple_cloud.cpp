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
#include "radiation_utils.hpp"  // getPhaseHenyeyGreenstein


// For grey cloud
Real SimpleCloud::getAttenuation(Real wave1, Real wave2, GridData const& gdata) const
{
  //Real result= 1.;
  //Real result= 1.E3*exp(- pow( ((log(q[IPR])-log(5.E6))/2.), 2)) ;
  //Real result= 1.E1*exp(- pow( ((log(q[IPR])-log(5.E5))/1.), 2)) ;
  Real csize = 1.E0*1.e-6; // one micron size particle
  Real qext = 1.E0;
  Real crho = 5.E3;
  //std::cout << q[0] << " " << q[6] << " " << q[7] << " " << q[9] << std::endl;
  return  gdata.q[imol_]*qext/(4./3.*csize*crho);     // -> 1/m
  //return q[imol_]*result*mixr_;     // -> 1/m
}

Real SimpleCloud::getSingleScateringAlbedo(Real wave1, Real wave2, GridData const& gdata) const
{
  // ssalb
  Real ww = 0.9;

  if (gdata.q[IPR] > 1)
    return ww;
  else
    return 0.;
}


void SimpleCloud::getPhaseMomentum(Real *pp, Real wave1, Real wave2, GridData const& gdata, int np) const
{
  Real gg=0.9;

  if (gdata.q[IPR] > 1)
    getPhaseHenyeyGreenstein(pp, 0, gg, np); // 0 for HENYEY_GREENSTEIN
  else
    getPhaseHenyeyGreenstein(pp, 0, 0.0, np); // 0 for HENYEY_GREENSTEIN
}
