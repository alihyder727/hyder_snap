// C/C++ headers
#include <sstream>
#include <stdexcept>

// Athena++ headers
#include "mwr_absorbers.hpp"
#include "absorption_functions.hpp"

MwrAbsorberCIA::MwrAbsorberCIA(RadiationBand *pband, Real xHe, Real xCH4, Real fequal):
    Absorber(pband, "mw_CIA", IDN), xHe_(xHe), xCH4_(xCH4), fequal_(fequal)
{
  std::stringstream msg;

  if ((xHe_ < 0.) || (xCH4_ < 0.) || (xHe_ + xCH4_ > 1.)) {
    msg << "### FATAL ERROR in MwrAbsorberCIA::MwrAbsorberCIA."
        << std::endl << "Value error in molar mixing ratios";
    throw std::runtime_error(msg.str().c_str());
  }
}

Real MwrAbsorberCIA::AbsorptionCoefficient(Real wave, Real const prim[]) const
{
  Real P = prim[IPR]/1.E5;  // pa -> bar
  Real T = prim[IDN];
  Real xdry = 1.;
  for (int i = 1; i < NMASS; ++i) xdry -= prim[i];
  Real XHe = xHe_*xdry;
  Real XCH4 = xCH4_*xdry;
  Real XH2 = (1. - xHe_ - xCH4_)*xdry;

  return 100.*absorption_coefficient_CIA(wave, P, T, XH2, XHe, XCH4, fequal_);
}
