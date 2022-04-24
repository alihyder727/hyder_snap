#ifndef CORRELATEDK_ABSORBER_HPP
#define CORRELATEDK_ABSORBER_HPP

// C/C++ header
#include <vector>

// Athena++ header
#include "absorber.hpp"

class CorrelatedKAbsorber: public Absorber {
public:
  CorrelatedKAbsorber(RadiationBand *pband):
    Absorber(pband) {}
  virtual ~CorrelatedKAbsorber() {}
  void loadCoefficient(std::string fname);
  //Real ckAbsorptionCoefficient(int mw, int mg, Real const prim[]) const;
  Real getAttenuation(Real g1, Real g2, Real const q[], Real const c[], Real const s[]) const;

protected:
  int len_[3];                  /**< length of interpolation axis */
  std::vector<Real> axis_;    /**< interpolation axis */
  std::vector<Real> kcoeff_;    /**< absorption coefficient */
};

#endif
