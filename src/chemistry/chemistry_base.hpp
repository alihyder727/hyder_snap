#ifndef CHEMISTRY_BASE_HPP
#define CHEMISTRY_BASE_HPP

// C/C++ header
#include <vector>
#include <array>
#include <map>

// Athena++ header
#include "../thermodynamics/thermodynamics.hpp"
#include "../particles/particles.hpp"
#include "../parameter_input.hpp"
#include "../athena.hpp"

// Eigen header files
#include "../math/eigen335/Eigen/Core"

class Chemisry;

template<typename T>
class ChemistryBase {
public:
// data
  Chemistry *pmy_chem;
  std::string myname;
  std::string particle_name;
  ChemistryBase<T> *prev, *next;

// functions
  ChemistryBase(Chemistry *pchem, ParameterInput *pin):
    pmy_chem(pchem), prev(nullptr), next(nullptr)
  {}

  void IntegrateDense(AthenaArray<Real> &u, AthenaArray<Real> &c, 
    Real time, Real dt);

  template<typename D1, typename D2>
  void AssembleReactionMatrix(Eigen::DenseBase<D1>& rate,
    Eigen::DenseBase<D2>& jac, Real const q[], Real cv, Real time);

  void ApplyConcentrationLimit(Real q[], Real const q0[]) {}

protected:
  //! reaction coefficients 
  std::map<std::string, Real> coeffs_;
  std::vector<int> index_;
  //! stores saturation vapor surplus
  std::array<Real, 1+NVAPOR> dqsat_;
  //! internal energy
  std::vector<Real> deltaU_;
};

#include "integrate_dense.hpp"

#endif
