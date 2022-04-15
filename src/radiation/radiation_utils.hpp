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
void getPhaseMomentum(int iphas, Real gg, int npmom, Real *pmom);
void packSpectralProperties(Real *buf, Real const *tau, Real const *ssa, Real const* pmom, int nlayer, int npmom);
void unpackSpectralProperties(Real *tau, Real *ssa, Real *pmom, Real const *buf, int slyr, int npmom, int nblocks, int npmom_max = -1);

#endif
