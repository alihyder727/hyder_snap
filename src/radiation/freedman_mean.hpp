#ifndef FREEDMAN_MEAN_HPP
#define FREEDMAN_MEAN_HPP

// C/C++ header
#include <vector>

// Athena++ header
#include "absorber.hpp"

// Richard S. Freedman 2011. APJS
class FreedmanMean: public Absorber {
public:
  FreedmanMean(RadiationBand *pband):
    Absorber(pband) {}
  virtual ~FreedmanMean() {}
  Real AbsorptionCoefficient(Real wave, Real const prim[]) const;
};

#endif
