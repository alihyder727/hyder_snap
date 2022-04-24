// C/C++ headers
#include <stdexcept>

// Athena++ headers
#include "../../math/interpolation.h"
#include "mwr_absorbers.hpp"
#include "absorption_functions.hpp"

MwrAbsorberH2S::MwrAbsorberH2S(RadiationBand *pband, int imol, Real xHe):
    Absorber(pband, "mw_H2S", imol), xHe_(xHe)
{
  std::stringstream msg;

  if ((xHe_ < 0.) || (xHe_ > 1.)) {
    msg << "### FATAL ERROR in MwrAbsorberH2S::MwrAbsorberH2S."
        << std::endl << "Value error in molar mixing ratios";
    throw std::runtime_error(msg.str().c_str());
  }
}

MwrAbsorberH2S::MwrAbsorberH2S(RadiationBand *pband, Real xHe, Real *xH2S, Real *pres, int np):
    Absorber(pband, "mw_H2S", IDN), xHe_(xHe)
{
  std::stringstream msg;

  for (int i = 0; i < np; ++i) {
    if ((xH2S[i] < 0.) || (xH2S[i] > 1.)) {
      msg << "### FATAL ERROR in MwrAbsorberNH3::MwrAbsorberNH3."
          << std::endl << "Value error in molar mixing ratios";
      throw std::runtime_error(msg.str().c_str());
    }
    ref_xh2s_.push_back(xH2S[i]);
    ref_pres_.push_back(pres[i]);
  }
}

Real MwrAbsorberH2S::getAttenuation(Real wave1, Real wave2,
    Real const q[], Real const c[], Real const s[]) const
{
  // adapted by cli (Cheng Li), Aug 30
  Real P = q[IPR]/1.E5;
  Real T = q[IDN];
  Real xdry = 1.;
  for (int i = 1; i <= NVAPOR; ++i) xdry -= q[i];
  Real XHe = xHe_*xdry;
  Real XH2, XH2S;

  if (method_ == 1) {
    XH2S = q[imol_];
    XH2 = xdry - XHe;
  } else {  // method_ == 2
    XH2S = interp1(q[IPR], ref_xh2s_.data(), ref_pres_.data(), ref_pres_.size())*xdry;;
    XH2 = xdry - XHe - XH2S;
  }

  // 1/cm -> 1/m
  return 100.*absorption_coefficient_H2S(wave1, P, T, XH2, XHe, XH2S);
}
