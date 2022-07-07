// C/C++ headers
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>

// Athena++ headers
#include "../math/interpolation.h"
#include "../utils/utils.hpp"
#include "hydrogen_cia.hpp"

void XizH2H2CIA::loadCoefficient(std::string fname, int bid)
{
  std::stringstream msg;
  if (!FileExists(fname)) {
    msg << "### FATAL ERROR in XizH2H2CIA::loadCoefficient."
        << std::endl << "File '" << fname << "' " << "does not exist.";
    throw std::runtime_error(msg.str().c_str());
  }
  len_[1] = GetNumCols(fname) - 1;
  len_[0] = GetNumRows(fname) - 1;
  
  std::ifstream infile(fname.c_str(), std::ios::in);
  axis_.resize(len_[0] + len_[1]);
  kcoeff_.resize(len_[0]*len_[1]);
  Real junk;
  if (infile.is_open()) {
    infile >> junk;
    for (int j = 0; j < len_[1]; j++) {
      infile >> axis_[len_[0] + j];
    }
    for (int k = 0; k < len_[0]; k++) {
      infile >> axis_[k];
      for (int j = 0; j < len_[1]; j++)
        infile >> kcoeff_[k*len_[1] + j];
    }
    infile.close();
  } else {
    msg << "### FATAL ERROR in XizH2H2CIA::loadCoefficient: ";
    msg << "Cannot open file: " << fname << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
}

Real XizH2H2CIA::getAttenuation(Real wave1, Real wave2,
    CellVariables const& var) const
{
  static const Real kBoltz = 1.3806504E-23;
  static const Real Lo = 2.68719E25;

  // first axis is wavenumber, second is temperature
  Real val, coord[2] = {wave1, var.q[IDN]};
  interpn(&val, coord, kcoeff_.data(), axis_.data(), len_, 2);

  Real amagat = var.q[IPR]/(kBoltz*var.q[IDN]*Lo);
  Real x0 = 1.;
  if (imol_ == 0) 
    for (int n = 1; n <= NVAPOR; ++n) x0 -= var.q[n];
  else
    x0 = var.q[imol_];

  return 100.*exp(-val)*x0*x0*amagat*amagat*mixr_*mixr_; // 1/cm -> 1/m
}
