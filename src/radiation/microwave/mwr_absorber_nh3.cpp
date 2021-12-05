// C/C++ headers
#include <sstream>
#include <stdexcept>

// Athena++ headers
#include "../../math/interpolation.h"
#include "mwr_absorbers.hpp"
#include "absorption_functions.hpp"

MwrAbsorberNH3::MwrAbsorberNH3(RadiationBand *pband, int imol, Real xHe, Real xH2O):
  Absorber(pband, "mw_NH3", imol), method_(1), xHe_(xHe), xH2O_(xH2O)
{
  std::stringstream msg;

  if ((xHe < 0.) || (xH2O < 0.) || (xHe + xH2O > 1.)) {
    msg << "### FATAL ERROR in MwrAbsorberNH3::MwrAbsorberNH3."
        << std::endl << "Value error in molar mixing ratios";
    ATHENA_ERROR(msg);
  }
}

MwrAbsorberNH3::MwrAbsorberNH3(RadiationBand *pband, int imol, Real xHe, Real *xH2O, Real *pres, int np):
  Absorber(pband, "mw_NH3", imol), method_(2), xHe_(xHe)
{
  std::stringstream msg;

  for (int i = 0; i < np; ++i) {
    if ((xH2O[i] < 0.) || (xH2O[i] > 1.)) {
      msg << "### FATAL ERROR in MwrAbsorberNH3::MwrAbsorberNH3."
          << std::endl << "Value error in molar mixing ratios";
      ATHENA_ERROR(msg);
    }
    ref_xh2o_.push_back(xH2O[i]);
    ref_pres_.push_back(pres[i]);
  }
}

MwrAbsorberNH3::MwrAbsorberNH3(RadiationBand *pband, std::vector<int> imols, Real xHe):
  Absorber(pband, "mw_NH3", imols), method_(3), xHe_(xHe)
{
  std::stringstream msg;

  if (imols_.size() != 2) {
    msg << "### FATAL ERROR in MwrAbsorberNH3::MwrAbsorberNH3."
        << std::endl << "Number of dependent molecules is not 2";
    ATHENA_ERROR(msg);
  }

  if ((xHe < 0.) || (xHe > 1.)) {
    msg << "### FATAL ERROR in MwrAbsorberNH3::MwrAbsorberNH3."
        << std::endl << "Value error in molar mixing ratios";
    ATHENA_ERROR(msg);
  }
}

Real MwrAbsorberNH3::Attenuation(Real wave, Real const q[], Real const c[], Real const s[]) const
{
  Real P = q[IPR]/1.E5; // pa -> bar
  Real P_idl = q[IPR]/1.E5; // pa -> bar
  Real T = q[IDN];
  Real xdry = 1.;
  for (int i = 1; i <= NVAPOR; ++i) xdry -= q[i];
  Real XHe = xHe_*xdry;
  Real XH2, XNH3, XH2O;

  if (method_ == 1) {
    XNH3 = q[imol_];
    XH2O = xH2O_*xdry;
    XH2 = xdry - XHe - XH2O;
  } else if (method_ == 2) {
    XNH3 = q[imol_];
    XH2O = interp1(q[IPR], ref_xh2o_.data(), ref_pres_.data(), ref_pres_.size())*xdry;
    XH2 = xdry - XHe - XH2O;
  } else {  // method_ == 3
    XNH3 = q[imols_[0]];
    XH2O = q[imols_[1]];
    XH2 = xdry - XHe;
  }

  Real abs;

  if (model_name_ == "Bellotti16")
    abs = absorption_coefficient_NH3_Bellotti(wave, P, P_idl, T, XH2, XHe, XNH3, XH2O);
  else if (model_name_ == "BellottiSwitch16")
    abs = absorption_coefficient_NH3_Bellotti_switch(wave, P, P_idl, T, XH2, XHe, XNH3, XH2O);
  else if (model_name_ == "Devaraj")
    abs = absorption_coefficient_NH3_Devaraj(wave, P, P_idl, T, XH2, XHe, XNH3, XH2O);
  else if (model_name_ == "Radtran")
    abs = absorption_coefficient_NH3_radtran(wave, P, T, XH2, XHe, XNH3);
  else // Hanley09
    abs = absorption_coefficient_NH3_Hanley(wave, P, P_idl, T, XH2, XHe, XNH3, XH2O);

  return 100.*abs;  // 1/cm -> 1/m
}
