#ifndef KESSLER94_HPP
#define KESSLER94_HPP

// C/C++ header
#include <sstream>

// Athena++ header
#include "../utils/utils.hpp"
#include "../thermodynamics/thermodynamics.hpp"
#include "../particles/particles.hpp"
#include "chemistry.hpp"
#include "chemistry_base.hpp"
#include "chemistry_solver.hpp"

class Kessler94 : public ChemistryBase<Kessler94> {
public:
// typedefs
  typedef ChemistrySolver<4> Solver;

// data
  Solver solver;

// functions
  Kessler94(Chemistry *pchem, ParameterInput *pin, std::string name) :
    ChemistryBase<Kessler94>(pchem, pin)
  {
    myname = name;
    particle_name = pin->GetString("chemistry", name + ".link_particle");

    coeffs_["condensation"] = pin->GetReal("chemistry", name + ".condensation");
    coeffs_["autoconversion"] = pin->GetReal("chemistry", name + ".autoconversion");
    coeffs_["accretion"] = pin->GetReal("chemistry", name + ".accretion");
    coeffs_["evaporation"] = pin->GetReal("chemistry", name + ".evaporation");

    index_.resize(4);
    index_[0] = IDN;
    index_[1] = pin->GetInteger("chemistry", name + ".link_vapor");
    index_[2] = NHYDRO;
    index_[3] = NHYDRO + 1;

    deltaU_.resize(NHYDRO + 2);
    std::fill(deltaU_.begin(), deltaU_.end(), 0.);
    deltaU_[index_[1]] = pin->GetReal("chemistry", name + ".deltaU");

    Particles *part = pchem->pmy_block->ppart->FindParticle(particle_name);
    Real cfloor_ = 1.;
    for (int n = 0; n < part->u.GetDim4(); ++n)
      cfloor_ = std::min(cfloor_, part->GetDensityFloor()/part->GetMolecularWeight(n));
  }

  template<typename T>
  void ApplySolutionLimit(Real c[], Real const c0[], Eigen::DenseBase<T>& dc, Real cv);

  template<typename D1, typename D2>
  void AssembleReactionMatrix(Eigen::DenseBase<D1>& rate,
    Eigen::DenseBase<D2>& jac, Real const q[], Real cv, Real time);

protected:
  Real cfloor_;
};

template<typename T>
void Kessler94::ApplySolutionLimit(Real c[], Real const c0[], Eigen::DenseBase<T>& dc, Real cv) {
  Thermodynamics *pthermo = pmy_chem->pmy_block->pthermo;
  int iT = index_[0];
  int iqv = index_[1];
  int iqc = index_[2];

  Real dU = deltaU_[iqv] - deltaU_[iqc];
  Real Rv = pthermo->GetRd()/pthermo->GetMassRatio(iqv);

  Real dqsat[1+NVAPOR];
  pthermo->SaturationSurplus(dqsat, c, VariableType::chem);

  if (c[iT] < 1. || dqsat[iqv]*dqsat_[iqv] < 0.) {
    std::memcpy(c, c0, (NHYDRO+2)*sizeof(Real));
    int iter = 0.;
    while (true) {
      Real tem = c[index_[0]];
      Real qs = c[iqv] - dqsat_[iqv];
      //std::cout << "qs = " << qs << std::endl;
      Real lf = pthermo->GetLatent(NVAPOR+iqv,tem) - Rv*tem;
      Real dqsdt = qs/tem*lf/(Rv*tem);
      Real dq = -dqsat_[iqv]/(dqsdt*dU/cv + 1.);
      //std::cout << "dq = " << dq << std::endl;
      for (int n = 0; n < Solver::Size; ++n)
        c[index_[n]] += dc(n)/dc(iqv)*dq;
      pthermo->SaturationSurplus(dqsat_.data(), c, VariableType::chem);
      //std::cout << "tem = " << c[0] << std::endl;
      iter++;
      if (std::abs(c[index_[0]] - tem) < 0.1 || iter > 10.) break;
    }
  }

  if (c[index_[2]] < cfloor_) {
    c[index_[2]] = 0;
    c[index_[1]] += c[index_[2]];
  }

  if (c[index_[3]] < cfloor_) {
    c[index_[3]] = 0;
    c[index_[1]] += c[index_[3]];
  }

  c[index_[1]] = std::max(0., c[index_[1]]);
}

#include "kessler94_impl.hpp"

#endif
