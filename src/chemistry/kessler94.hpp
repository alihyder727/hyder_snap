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
    particle_name = pin->GetString("chemistry", name + ".link_particle");

    coeffs_["condensation"] = pin->GetReal("chemistry", name + ".condensation");
    coeffs_["autoconversion"] = pin->GetReal("chemistry", name + ".autoconversion");
    coeffs_["accretion"] = pin->GetReal("chemistry", name + ".accretion");
    coeffs_["evaporation"] = pin->GetReal("chemistry", name + ".evaporation");

    qindex_.resize(4);
    qindex_[0] = IDN;
    qindex_[1] = pin->GetInteger("chemistry", name + ".link_vapor");
    qindex_[2] = NHYDRO;
    qindex_[3] = NHYDRO + 1;

    deltaU_.resize(NHYDRO + 2);
    std::fill(deltaU_.begin(), deltaU_.end(), 0.);
    deltaU_[qindex_[1]] = pin->GetReal("chemistry", name + ".deltaU");

  }

  void ApplyChemicalLimits(Real q[], Real const q0[], Real cv)
  {
    Thermodynamics *pthermo = pmy_block->pthermo;

    int iT = qindex_[0];
    int iv = qindex_[1];
    int ic = qindex_[2];
    int ip = qindex_[3];

    Real T0 = q0[iT];
    Real L = deltaU_[iv] - deltaU_[ic];
    Real Rv = pthermo->GetRd()/pthermo->GetMassRatio(iv);
    Real lf = pthermo->GetLatent(ic,T0) - Rv*T0;
    Real qs = q0[iv] - dqsat_[iv];
    Real dqsdt = qs/T0*lf/(Rv*T0);

    if (dqsat_[iv]*(q[iv] - qs) < 0) {
      Real dq = (q[iv] - qs)/(1. + dqsdt*L/cv);
      Real dqc = dq, dqp = 0.;
      q[iv] -= dq;
      q[ic] += dq;
      if (q[ic] < 0.) {
        dqp = -q[ic];
        dqc = dq - dqp;
        q[ip] += q[ic];
        q[ic] = 0.;
      }
      q[iT] += L*dq/cv;
    }
  }

  template<typename D1, typename D2>
  void AssembleReactionMatrix(Eigen::DenseBase<D1>& rate,
    Eigen::DenseBase<D2>& jac, Real const q[], Real cv, Real time);
};

#include "kessler94_impl.hpp"

#endif
