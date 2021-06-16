#ifndef CHEMISTRY_BASE_HPP
#define CHEMISTRY_BASE_HPP

// C/C++ header
#include <vector>
#include <array>
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
  std::string myname;
  std::string particle_name;
  ChemistryBase<T> *prev, *next;

// functions
  ChemistryBase(MeshBlock *pmb, ParameterInput *pin):
    pmy_block(pmb), prev(nullptr), next(nullptr)
  {
    int nc1 = pmb->ncells1, nc2 = pmb->ncells2, nc3 = pmb->ncells3;
    cvt_.NewAthenaArray(nc3, nc2, nc1);
    mol_.NewAthenaArray(nc3, nc2, nc1);
  }

  void SetTotalCv(AthenaArray<Real> const& u, Particles *ppart,
    int il, int iu, int jl, int ju, int kl, int ku)
  {
    Thermodynamics *pthermo = pmy_block->pthermo;
    cvt_.ZeroClear();
    mol_.ZeroClear();
    Real gm1 = pmy_block->peos->GetGamma() - 1.;
    Real cvd = pthermo->GetRd()/gm1;
    for (int k = kl; k <= ku; ++k)
      for (int j = jl; j <= ju; ++j)
        for (int i = il; i <= iu; ++i) {
          for (int n = 0; n <= NVAPOR; ++n)
            cvt_(k,j,i) += u(n,k,j,i)*pthermo->GetCvRatio(n)*cvd;
          Particles *p = ppart;
          while (p != nullptr) {
            cvt_(k,j,i) += p->GetTotalCv(k,j,i);
            mol_(k,j,i) += p->GetMolarDensity(k,j,i);
            p = p->next;
          }
        }
  }

  void IntegrateDense(AthenaArray<Real> &u, AthenaArray<Real> &c, 
    Real time, Real dt);

  template<typename D1, typename D2>
  void AssembleReactionMatrix(Eigen::DenseBase<D1>& rate,
    Eigen::DenseBase<D2>& jac, Real const q[], Real cv, Real time);

protected:
  //! total cv, J/(K m^3)
  AthenaArray<Real> cvt_;
  //! molar density, mol/m^3
  AthenaArray<Real> mol_;
  //! reaction coefficients 
  std::map<std::string, Real> coeffs_;
  std::vector<int> qindex_;
  //! stores saturation vapor surplus
  std::array<Real, 1+NVAPOR> dqsat_;
  //! internal energy
  std::vector<Real> deltaU_;
};

#include "integrate_dense.hpp"

#endif
