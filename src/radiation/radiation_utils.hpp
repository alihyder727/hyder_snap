#ifndef RADIATION_UTILS_HPP
#define RADIATION_UTILS_HPP

// C/C++ header
#include <string>

// Athena++ header
#include "../athena.hpp"


/*void WriteTopFlux(std::string fname) const;
void WriteTopRadiance(std::string fname) const;
void WriteOpticalDepth(std::string fname) const;
void WriteHeatingRate(std::string fname, AthenaArray<Real> const& flux,
      AthenaArray<Real> const& hr, Real const* level); */
void GetPhaseMomentum(int iphas, Real gg, int npmom, Real *pmom);

#endif
