// C/C++ headers
#include <iostream>
#include <algorithm>
#include <stdexcept>
#include <cstring>
#include <cmath>
#include <cassert>  // assert

// Athena++ headers
#include "../athena_arrays.hpp"
#include "../globals.hpp"
#include "absorber.hpp"
#include "radiation.hpp"
#include "../mesh/mesh.hpp"
#include "../thermodynamics/thermodynamics.hpp"
#include "freedman_simple.hpp"


FreedmanSimple::FreedmanSimple(RadiationBand *pband, ParameterInput *pin):
  Absorber(pband, "freedman_simple")
{
  char str[80];
  sprintf(str, "%s.%s.scale", pband->myname.c_str(), myname.c_str());
  scale_ = pin->GetOrAddReal("radiation", str, 1.);
}

// xiz semigrey
Real FreedmanSimple::getAttenuation(Real wave1, Real wave2, GridData const& gdata) const
{ 
  static const Real Rgas = 8.314462;
  const Absorber *pabs = this;
  Thermodynamics * pthermo = pabs->pmy_band->pmy_rad->pmy_block->pthermo;
  Real mu = Rgas/pthermo->GetRd();
  Real result;
  Real p = gdata.q[IPR];
  Real T = gdata.q[IDN];

/*xiz semigrey
//Tan and Komacek 2019 simple fit m^2/kg
//    k_VIS= 10.0_dp**(0.0478_dp*Pl10**2 - 0.1366_dp*Pl10 - 3.2095_dp)
//    k_IR = 10.0_dp**(0.0498_dp*Pl10**2 - 0.1329_dp*Pl10 - 2.9457_dp)
  Real logp = log10(p); // Pa
if (wave < 40000.) //for semigrey
  result = pow(10.0, (0.0498*pow(logp,2.) - 0.1329*logp - 2.9457));
else
  result = pow(10.0, (0.0478*pow(logp,2.) - 0.1366*logp - 3.2095));
*/

//Komacek et al. 2017
//if (wave < 40000.) //for semigrey
//  result = 2.28e-6*pow(p, 0.53);
//else
  result = 2.28e-6*pow(p, 0.53); //visible opacity scale in disort
//

  Real dens   = p*mu/(Rgas * T); // kg/m^3
//  if (p > 5e1)
    return scale_*dens*result;     // -> 1/m
//  else
//    return 0.;
}
