// C/C++ headers
#include <stdexcept>

// Athena++ headers
#include "mwr_absorbers.hpp"
#include "absorption_functions.hpp"

Real MwrAbsorberElectron::Attenuation(Real wave, Real const q[], Real const c[], Real const s[]) const
{
  Real P = q[IPR]/1.E5; // pa -> bar
  Real T = q[IDN];

  Real abs;

  if (model_name_ == "Reference")
    abs = attenuation_freefree_Reference(wave, P, T);
  else if (model_name_ == "ChengLi")
    abs = attenuation_freefree_Chengli(wave, P, T);
  else // AppletonHartree
    abs = attenuation_appleton_hartree_nomag(wave, P, T, s[imol_]);

  return 100.*abs;  // 1/cm -> 1/m
}
