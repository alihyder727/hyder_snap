// C/C++ headers
#include <stdexcept>

// Athena++ headers
#include "mwr_absorbers.hpp"
#include "absorption_functions.hpp"

Real MwrAbsorberElectron::getAttenuation(Real wave1, Real wave2,
    CellVariables const& var) const
{
  Real P = var.q[IPR]/1.E5; // pa -> bar
  Real T = var.q[IDN];

  Real abs;

  if (model_name_ == "Reference")
    abs = attenuation_freefree_Reference(wave1, P, T);
  else if (model_name_ == "ChengLi")
    abs = attenuation_freefree_Chengli(wave1, P, T);
  else // AppletonHartree
    abs = attenuation_appleton_hartree_nomag(wave1, P, T, var.s[imol_]);

  return 100.*abs;  // 1/cm -> 1/m
}
