// C/C++ headers
#include <stdexcept>

// Athena++ headers
#include "mwr_absorbers.hpp"
#include "absorption_functions.hpp"

MwrAbsorberH2O::MwrAbsorberH2O(RadiationBand *pband, int imol, Real xHe, Real scale):
    Absorber(pband, "mw_H2O", imol), xHe_(xHe), scale_(scale)
{
  std::stringstream msg;

  if ((xHe_ < 0.) || (xHe_ > 1.)) {
    msg << "### FATAL ERROR in MwrAbsorberPH3::MwrAbsorberH2O."
        << std::endl << "Value error in molar mixing ratios";
    ATHENA_ERROR(msg);
  }
}

Real MwrAbsorberH2O::Attenuation(Real wave, Real const q[], Real const c[], Real const s[]) const
{
  Real P = q[IPR]/1.E5; // pa -> bar
  Real T = q[IDN];
  Real xdry = 1.;
  for (int i = 1; i <= NVAPOR; ++i) xdry -= q[i];
  Real XHe = xHe_*xdry;
  Real XH2 = xdry - XHe;
  Real XH2O = q[imol_];

  Real abs;

  if (model_name_ == "deBoer")
    abs = attenuation_H2O_deBoer(wave, P, T, XH2, XHe, XH2O);
  else if (model_name_ == "Waters")
    abs = attenuation_H2O_Waters(wave, P, T, XH2, XHe, XH2O);
  else if (model_name_ == "Goodman")
    abs = attenuation_H2O_Goodman(wave, P, T, XH2, XHe, XH2O);
  else // Karpowicz
    abs = attenuation_H2O_Karpowicz(wave, P, T, XH2, XHe, XH2O, scale_);

  return 100.*abs;  // 1/cm -> 1/m
}
