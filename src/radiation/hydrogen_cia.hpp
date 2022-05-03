#ifndef HYDROGEN_CIA_HPP
#define HYDROGEN_CIA_HPP

// C/C++ header
#include <iostream>
#include <vector>

// Athena++ header
#include "absorber.hpp"

class XizH2H2CIA: public Absorber {
public:
  XizH2H2CIA(RadiationBand *pband, int id, Real mixr):
    Absorber(pband, "H2-H2", id, mixr) {}
  virtual ~XizH2H2CIA() {}
  void loadCoefficient(std::string fname, int bid = -1);
  Real getAttenuation(Real wave1, Real wave2, GridData const& gdata) const;

protected:
  int len_[2];
  std::vector<Real> axis_;
  std::vector<Real> kcoeff_;
};


class XizH2HeCIA: public Absorber {
public:
  XizH2HeCIA(RadiationBand *pband, int id, Real mixr1, Real mixr2):
    Absorber(pband, "H2-He", id, mixr1), mixr2_(mixr2) {}
  virtual ~XizH2HeCIA() {}
  void loadCoefficient(std::string fname, int bid = -1);
  Real getAttenuation(Real wave1, Real wave2, GridData const& gdata) const;

protected:
  Real mixr2_;
  int len_[2];
  std::vector<Real> axis_;
  std::vector<Real> kcoeff_;
};

class OrtonCIA: public Absorber {
public:
  OrtonCIA(RadiationBand *pband, int id, Real mixr):
    Absorber(pband, "H2-H2", id, mixr) {}
  virtual ~OrtonCIA() {}
  void loadCoefficient(std::string fname, int bid = -1);
  Real getAttenuation(Real wave1, Real wave2, GridData const& gdata) const;

protected:
  int len_[2];
  std::vector<Real> axis_;
  std::vector<Real> kcoeff_;
};

#endif
