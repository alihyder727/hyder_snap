// C/C++ headers
#include <stdexcept>

// Athena++ headers
#include "../../math/interpolation.h"
#include "mwr_absorbers.hpp"
#include "absorption_functions.hpp"

MwrAbsorberPH3::MwrAbsorberPH3(RadiationBand *pband, int imol, Real xHe):
    Absorber(pband, "mw_PH3", imol), xHe_(xHe)
{
  std::stringstream msg;

  if ((xHe_ < 0.) || (xHe_ > 1.)) {
    msg << "### FATAL ERROR in MwrAbsorberPH3::MwrAbsorberPH3."
        << std::endl << "Value error in molar mixing ratios";
    throw std::runtime_error(msg.str().c_str());
  }
}

MwrAbsorberPH3::MwrAbsorberPH3(RadiationBand *pband, Real xHe, Real *xPH3, Real *pres, int np):
    Absorber(pband, "mw_PH3", IDN), xHe_(xHe)
{
  std::stringstream msg;

  for (int i = 0; i < np; ++i) {
    if ((xPH3[i] < 0.) || (xPH3[i] > 1.)) {
      msg << "### FATAL ERROR in MwrAbsorberNH3::MwrAbsorberNH3."
          << std::endl << "Value error in molar mixing ratios";
      throw std::runtime_error(msg.str().c_str());
    }
    ref_xph3_.push_back(xPH3[i]);
    ref_pres_.push_back(pres[i]);
  }
}

Real MwrAbsorberPH3::getAttenuation(Real wave1, Real wave2,
    CellVariables const& var) const
{
  Real P = var.q[IPR]/1.E5; // pa -> bar
  Real T = var.q[IDN];
  Real xdry = 1.;
  for (int i = 1; i <= NVAPOR; ++i) xdry -= var.q[i];
  Real XHe = xHe_*xdry;
  Real XH2, XPH3;

  if (method_ == 1) {
    XPH3 = var.q[imol_];
    XH2 = xdry - XHe;
  } else {  // method_ == 2
    XPH3 = interp1(var.q[IPR], ref_xph3_.data(), ref_pres_.data(), ref_pres_.size())*xdry;;
    XH2 = xdry - XHe - XPH3;
  }

  Real abs;

  if (model_name_ == "Radtran")
    abs = absorption_coefficient_PH3_radtran(wave1, P, T, XH2, XHe, XPH3);
  else // Hoffman
    abs = absorption_coefficient_PH3_Hoffman(wave1, P, T, XH2, XHe, XPH3);

  return 100.*abs;  // 1/cm -> 1/m
}
