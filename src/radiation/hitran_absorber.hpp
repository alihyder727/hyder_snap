#ifndef HITRAN_ABSORBER_HPP
#define HITRAN_ABSORBER_HPP

// C/C++ header
#include <vector>

// Athena++ header
#include "absorber.hpp"

class HitranAbsorber: public Absorber {
  friend std::ostream& operator<<(std::ostream &os, HitranAbsorber const& ab);
  //friend HitranAbsorber const& MakeCKAbsorber<>(HitranAbsorber const& albl,
  //  int const *ck_axis, Real const *ck_wave, int nbins);
public:
  HitranAbsorber(RadiationBand *pband, std::string name, int imol, Real mixr = 1.): 
    Absorber(pband, name, imol, mixr) {}
  virtual ~HitranAbsorber() {}
  void LoadCoefficient(std::string fname);
  Real AbsorptionCoefficient(Real wave, Real const prim[]) const;

protected:
  int len_[3];                  /**< length of interpolation axis */
  std::vector<Real> axis_;    /**< interpolation axis */
  std::vector<Real> kcoeff_;    /**< absorption coefficient */
  AthenaArray<Real> refatm_;    /**< reference atmosphere */
  Real RefTemp_(Real pres) const;
};

#endif
