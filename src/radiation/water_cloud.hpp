#ifndef WATER_CLOUD_HPP
#define WATER_CLOUD_HPP

// C/C++ header
#include <vector>

// Athena++ header
#include "absorber.hpp"

class SimpleCloud: public Absorber {
public:
  SimpleCloud(RadiationBand *pband, int id, Real mixr = 1.):
      Absorber(pband, "H2O_l", id, mixr) {}
  Real getAttenuation(Real wave1, Real wave2,
      Real const q[], Real const c[], Real const s[]) const;
  Real getSingleScateringAlbedo(Real wave1, Real wave2, 
      Real const q[], Real const c[], Real const s[]) const;
  void getPhaseMomentum(Real *pp, Real wave1, Real wave2,
      Real const q[], Real const c[], Real const s[], int np) const;
};

class FuWaterLiquidCloud: public Absorber {
public:
  FuWaterLiquidCloud(RadiationBand *pband, int id):
      Absorber(pband, "H2O_l", id) {}
  Real getAttenuation(Real wave1, Real wave2,
      Real const q[], Real const c[], Real const s[]) const;
  Real getSingleScateringAlbedo(Real wave1, Real wave2,
      Real const q[], Real const c[], Real const s[]) const;
  void getPhaseMomentum(Real *pp, Real wave1, Real wave2, 
      Real const q[], Real const c[], Real const s[], int np) const;
};

class FuWaterIceCloud: public Absorber {
public:
  FuWaterIceCloud(RadiationBand *pband, int id):
      Absorber(pband, "H2O_s", id) {}
  Real getAttenuation(Real wave1, Real wave2,
      Real const q[], Real const c[], Real const s[]) const;
  Real getSingleScateringAlbedo(Real wave1, Real wave2,
      Real const q[], Real const c[], Real const s[]) const;
  void getPhaseMomentum(Real *pp, Real wave1, Real wave2, 
      Real const q[], Real const c[], Real const s[], int np) const;
};

class XuWaterIceCloud: public Absorber {
public:
  XuWaterIceCloud(RadiationBand *pband, int id):
    Absorber(pband, "H2O_s", id) {}
  void loadCoefficient(std::string fname);
  Real getAttenuation(Real wave1, Real wave2,
      Real const q[], Real const c[], Real const s[]) const;
  Real getSingleScateringAlbedo(Real wave1, Real wave2,
      Real const q[], Real const c[], Real const s[]) const;
  void getPhaseMomentum(Real *pp, Real wave1, Real wave2,
      Real const q[], Real const c[], Real const s[], int np) const;

protected:
  int len_[2];
  std::vector<Real> axis_;
  std::vector<Real> ssalb_;
  std::vector<Real> gg_;
};

#endif
