#ifndef OXYGEN_CIA_HPP
#define OXYGEN_CIA_HPP

// C/C++ header
#include <vector>

// Athena++ header
#include "../athena.hpp"

class O2O2CIA: public Absorber {
public:
  O2O2CIA(RadiationBand *pband, int id, Real mixr):
    Absorber(pband, "O2-O2", id, mixr) {}
  virtual ~O2O2CIA() {}
  void LoadCoefficient(std::string fname);
  Real AbsorptionCoefficient(Real wave, Real const prim[]) const;

protected:
  int fd_len_[2];
  int a1dg_x3sg00_len_[2];
  int a1dg_x3sg10_len_[2];
  int ab_len_[2];
  int other_len_[2];
  std::vector<Real> fd_axis_;
  std::vector<Real> a1dg_x3sg00_axis_;
  std::vector<Real> a1dg_x3sg10_axis_;
  std::vector<Real> ab_axis_;
  std::vector<Real> other_axis_;
  std::vector<Real> fd_;
  std::vector<Real> a1dg_x3sg00_;
  std::vector<Real> a1dg_x3sg10_;
  std::vector<Real> ab_;
  std::vector<Real> other_;
};

#endif
