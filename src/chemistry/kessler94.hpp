#ifndef KESSLER94_HPP
#define KESSLER94_HPP

// C/C++ header
#include <sstream>

// Athena++ header
#include "chemistry_base.hpp"
#include "../utils/utils.hpp"

class Kessler94 : public ChemistryBase<Kessler94> {
public:
// typedefs
  typedef ChemistrySolver<4> Solver;

// data
  Solver solver;

// functions
  Kessler94(MeshBlock *pmb, ParameterInput *pin, std::string name) :
    ChemistryBase<Kessler94>(pmb, pin)
  {
    myname = name;

    coeffs_["condensation"] = pin->GetReal("chemistry", name + ".condensation");
    coeffs_["autoconversion"] = pin->GetReal("chemistry", name + ".autoconversion");
    coeffs_["accretion"] = pin->GetReal("chemistry", name + ".accretion");
    coeffs_["evaporation"] = pin->GetReal("chemistry", name + ".evaporation");

    index_.resize(4);
    index_[0] = IDN;
    index_[1] = pin->GetInteger("chemistry", name + ".ivapor");
    index_[2] = NHYDRO;
    index_[3] = NHYDRO + 1;

    Thermodynamics *pthermo = pmb->pthermo;
    Real Rgas = Thermodynamics::Rgas;
    Real cvd = Rgas/(pmb->peos->GetGamma() - 1.);
    Real Rd = pthermo->GetRd();
    for (int i = 0; i <= NVAPOR; ++i) {
      cv_.push_back(pthermo->GetCvRatio(i)*pthermo->GetMassRatio(i)*cvd);
      //deltaU_.push_back(pthermo->GetLatent(i)*pthermo->GetMassRatio(i)*Rgas/Rd);
      deltaU_.push_back(0.);
    }

    int ic = index_[1] + NVAPOR;

    cv_.push_back(pthermo->GetCvRatio(ic)*pthermo->GetMassRatio(ic)*cvd);
    cv_.push_back(pthermo->GetCvRatio(ic)*pthermo->GetMassRatio(ic)*cvd);

    deltaU_.push_back(-pthermo->GetLatent(ic)*pthermo->GetMassRatio(ic)*Rgas/Rd);
    deltaU_.push_back(-pthermo->GetLatent(ic)*pthermo->GetMassRatio(ic)*Rgas/Rd);

    std::string str = pin->GetString("chemistry", name + ".particle");
    pmy_part = pmb->ppart->FindParticle(str);
  }

  void ApplyChemicalLimits(Real q[], Real const q0[]) {
    Thermodynamics *pthermo = pmy_block->pthermo;

    int iT = index_[0];
    int iv = index_[1];
    int ic = index_[2];
    int ip = index_[3];

    Real T0 = q0[iT];
    Real L = deltaU_[iv] - deltaU_[ic];
    Real Rv = pthermo->GetRd()/pthermo->GetMassRatio(iv);
    Real lf = pthermo->GetLatent(ic,T0) - Rv*T0;
    Real qs = q0[iv] - dqsat_[iv];
    Real dqsdt = qs/T0*lf/(Rv*T0);

    if (dqsat_[iv]*(q[iv] - qs) < 0) {
      Real cvt = GetCvTotal(q);
      Real dq = (q[iv] - qs)/(1. + dqsdt*L/cvt);
      Real dqc = dq, dqp = 0.;
      q[iv] -= dq;
      q[ic] += dq;
      if (q[ic] < 0.) {
        dqp = -q[ic];
        dqc = dq - dqp;
        q[ip] += q[ic];
        q[ic] = 0.;
      }
      q[iT] += L*dq/cvt;
    }
  }

  template<typename D1, typename D2>
  void AssembleReactionMatrix(Eigen::DenseBase<D1>& rate,
    Eigen::DenseBase<D2>& jac, Real const q[], Real time);
};

#include "kessler94_impl.hpp"

#endif
