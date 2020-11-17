// C/C++ headers
#include <stdexcept>

// Athena++ headers
#include "mwr_absorbers.hpp"
#include "absorption_functions.hpp"

Real MwrAbsorberFreeFree::AbsorptionCoefficient(Real wave, Real const prim[]) const
{
  Real P = prim[IPR]/1.E5; // pa -> bar
  Real T = prim[IDN];

  Real abs;

  if (model_name_ == "Reference")
    abs = absorption_coefficient_freefree_Reference(wave, P, T);
  else // Cheng Li
    abs = absorption_coefficient_freefree_Chengli(wave, P, T);

  return 100.*abs;  // 1/cm -> 1/m
}
