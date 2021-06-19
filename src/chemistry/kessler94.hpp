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

    index_.resize(4);
    index_[0] = IDN;
    index_[1] = pin->GetInteger("chemistry", name + ".link_vapor");
    index_[2] = NHYDRO;
    index_[3] = NHYDRO + 1;

    deltaU_.resize(NHYDRO + 2);
    std::fill(deltaU_.begin(), deltaU_.end(), 0.);
    deltaU_[index_[1]] = pin->GetReal("chemistry", name + ".deltaU");

    qfloor_ = pin->GetOrAddReal("chemistry", name + ".qfloor", 0.);
  }

  void ApplyConcentrationLimit(Real c[]) {
    /*if (!(q[index_[0]] > 0.)) {
      for (int n = 0; n < 4; ++n)
        std::cout << q0[n] << " " ;
      std::cout << std::endl;
      for (int n = 0; n < 4; ++n)
        std::cout << q1[n] << " " ;
      std::cout << std::endl;
      for (int n = 0; n < 4; ++n)
        std::cout << q2[n] << " " ;
      std::cout << std::endl;
      for (int n = 0; n < 4; ++n)
        std::cout << q[n] << " " ;
      std::cout << std::endl;
    }*/
    assert(c[index_[0]] > 0.);

    if (c[index_[2]] < qfloor_) {
      c[index_[2]] = 0;
      c[index_[1]] += c[index_[2]];
    }

    if (c[index_[3]] < qfloor_) {
      c[index_[3]] = 0;
      c[index_[1]] += c[index_[3]];
    }

    c[index_[1]] = std::max(0., c[index_[1]]);
  }

  template<typename D1, typename D2>
  void AssembleReactionMatrix(Eigen::DenseBase<D1>& rate,
    Eigen::DenseBase<D2>& jac, Real const q[], Real cv, Real time);

protected:
  Real qfloor_;
};

#include "kessler94_impl.hpp"

#endif
