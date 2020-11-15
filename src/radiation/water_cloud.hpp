#ifndef WATER_CLOUD_HPP
#define WATER_CLOUD_HPP

// C/C++ header
#include <vector>

// Athena++ header
#include "absorber.hpp"

class FuWaterLiquidCloud: public Absorber {
public:
  FuWaterLiquidCloud(RadiationBand *pband, int id):
    Absorber(pband, "H2O_l", id) {}
  Real AbsorptionCoefficient(Real wave, Real const prim[]) const;
  Real SingleScateringAlbedo(Real wave, Real const prim[]) const;
  void PhaseMomentum(Real wave, Real const prim[], Real *pp, int np) const;
};

class FuWaterIceCloud: public Absorber {
public:
  FuWaterIceCloud(RadiationBand *pband, int id):
    Absorber(pband, "H2O_s", id) {}
  Real AbsorptionCoefficient(Real wave, Real const prim[]) const;
  Real SingleScateringAlbedo(Real wave, Real const prim[]) const;
  void PhaseMomentum(Real wave, Real const prim[], Real *pp, int np) const;
};

class XuWaterIceCloud: public Absorber {
public:
  XuWaterIceCloud(RadiationBand *pband, int id):
    Absorber(pband, "H2O_s", id) {}
  void LoadCoefficient(std::string fname);
  Real AbsorptionCoefficient(Real wave, Real const prim[]) const;
  Real SingleScateringAlbedo(Real wave, Real const prim[]) const;
  void PhaseMomentum(Real wave, Real const prim[], Real *pp, int np) const;

protected:
  int len_[2];
  std::vector<Real> axis_;
  std::vector<Real> ssalb_;
  std::vector<Real> gg_;
};

#endif
