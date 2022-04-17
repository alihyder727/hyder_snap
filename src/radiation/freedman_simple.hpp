#ifndef FREEDMAN_SIMPLE_HPP
#define FREEDMAN_SIMPLE_HPP

// C/C++ header
#include <vector>

// Athena++ header
#include "absorber.hpp"

class FreedmanSimple: public Absorber {
public:
  FreedmanSimple(RadiationBand *pband, ParameterInput *pin);
  virtual ~FreedmanSimple() {}
  Real Attenuation(Real wave, Real const q[], Real const c[], Real const s[]) const;

private:
  Real scale_;
};

#endif
