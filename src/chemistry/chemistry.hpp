#ifndef CHEMISTRY_HPP
#define CHEMISTRY_HPP

// C/C++ header
#include <vector>

// Athena++ header
#include "../athena.hpp"

class MeshBlock;
class ParameterInput;
template<typename T> class AthenaArray;

class Chemistry {
public:
  // data
  Real last_time, min_dt;

  // functions
  Chemistry(MeshBlock *pmb, ParameterInput *pin);
  ~Chemistry();
  void EvolveOneStep(AthenaArray<Real> &u, Real time, Real dt);
  void AddSedimentationFlux(AthenaArray<Real>& x1flux, AthenaArray<Real> const& wr,
    int k, int j, int il, int iu);
  void ConcentrationCorrection(AthenaArray<Real> &u);
  void AddFrictionalHeating(AthenaArray<Real> &u, AthenaArray<Real> const& w, Real dt);
  
  virtual void SedimentationVelocity(Real vsed[], Real const w[], Real temp = 0.);
  virtual void AssembleReactionMatrix(Real *r0, Real **r1, Real const q[], Real time);

protected:
  MeshBlock *pmy_block_;
  Real **identity_;
  Real *r0_;
  Real **r1_;

  // integration order
  int order_;
  int max_iter_;
  Real gamma_;
  Real ftol_;

  // reaction coefficients 
  std::vector<Real> kc_;

  // default sedimentation velocity
  //Real vsed_default_[NMASS];
};

class Kessler94 : public Chemistry {
public:
  Kessler94(MeshBlock *pmb, ParameterInput *pin):
    Chemistry(pmb, pin) {}
  void AssembleReactionMatrix(Real *r0, Real **r1, Real const q[], Real time);
};

#endif
