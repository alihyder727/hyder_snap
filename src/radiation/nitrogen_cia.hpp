#ifndef NITROGEN_CIA_HPP
#define NITROGEN_CIA_HPP

// C/C++ header
#include <vector>

// Athena++ header
#include "absober.hpp"

class N2N2CIA: public Absorber {
public:
  N2N2CIA(RadiationBand *pband, int id, Real mixr):
    Absorber("N2-N2", id, mixr) {}
  virtual ~N2N2CIA() {}
  void LoadCoefficient(std::string fname);
  Real AbsorptionCoefficient(Real wave, Real const prim[]) const;

protected:
  int rt_len_[2];
  int fd1_len_[2];
  int fd2_len_[2];
  std::vector<Real> rt_axis_;
  std::vector<Real> fd1_axis_;
  std::vector<Real> fd2_axis_;
  std::vector<Real> rt_;
  std::vector<Real> fd1_;
  std::vector<Real> fd2_;
};

#endif
