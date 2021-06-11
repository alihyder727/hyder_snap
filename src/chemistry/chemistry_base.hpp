#ifndef CHEMISTRY_BASE_HPP
#define CHEMISTRY_BASE_HPP

// C/C++ header
#include <vector>
#include <map>

// Athena++ header
#include "../thermodynamics/thermodynamics.hpp"
#include "../particles/particles.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "../athena.hpp"

// Eigen header files
#include "../math/eigen335/Eigen/Core"

template<typename T>
class ChemistryBase {
public:
// data
  MeshBlock *pmy_block;
  Particles *pmy_part;
  std::string myname;
  ChemistryBase<T> *prev, *next;

// functions
  ChemistryBase(MeshBlock *pmb, ParameterInput *pin):
    pmy_block(pmb), prev(nullptr), next(nullptr) {}

  Real GetCvTotal(Real const q[]) {
    Real cvt = 0, qd = 1;
    for (int i = 1; i <= NVAPOR; ++i) {
      qd -= q[i];
      cvt += q[i]*cv_[i];
    }
    for (int i = 1+NVAPOR; i < cv_.size(); ++i)
      cvt += q[i]*cv_[i];
    cvt = cvt + qd*cv_[0];
    return cvt;
  }

  void IntegrateDense(AthenaArray<Real> &u, AthenaArray<Real> &dc,
    AthenaArray<Real> const& c, Real time, Real dt);

  template<typename D1, typename D2>
  void AssembleReactionMatrix(Eigen::DenseBase<D1>& rate,
    Eigen::DenseBase<D2>& jac, Real const q[], Real time);

  void ApplyChemicalLimits(Real q0[], Real q[]);

protected:
  //! reaction coefficients 
  std::map<std::string, Real> coeffs_;
  std::vector<int> index_;
  //! molar cv ratios
  std::vector<Real> cv_;
  //! internal energy
  std::vector<Real> deltaU_;
  //! stores saturation surplus
  std::vector<Real> dqsat_;
};

#include "integrate_dense.hpp"

#endif
