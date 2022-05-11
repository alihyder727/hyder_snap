// C/C++ headers
#include <stdexcept>

// Athena++ headers
#include "mwr_absorbers.hpp"
#include "absorption_functions.hpp"

Real MwrAbsorberElectron::getAttenuation(Real wave1, Real wave2,
    GridData const& gdata) const
{
  Real P = gdata.q[IPR]/1.E5; // pa -> bar
  Real T = gdata.q[IDN];

  Real abs;

  if (model_name_ == "Reference")
    abs = attenuation_freefree_Reference(wave1, P, T);
  else if (model_name_ == "ChengLi")
    abs = attenuation_freefree_Chengli(wave1, P, T);
  else // AppletonHartree
    abs = attenuation_appleton_hartree_nomag(wave1, P, T, gdata.s[imol_]);

  return 100.*abs;  // 1/cm -> 1/m
}
