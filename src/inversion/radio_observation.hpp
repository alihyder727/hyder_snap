#ifndef RADIO_OBSERVATION_HPP_
#define RADIO_OBSERVATION_HPP_

// C/C++ headers
#include <vector>

// Athena++ headers
#include "../athena.hpp"

// Eigen header files
#include "../math/eigen335/Eigen/Core"

class MeshBlock;
class ParameterInput;
class Inversion;

class RadioObservation {
public:
  // data
  Eigen::VectorXd target;
  Eigen::MatrixXd icov;

	std::vector<int> ix;
	std::vector<Real> plevel;

  // functions
  RadioObservation(Inversion *pinvt, ParameterInput *pin);
  void ReadObservationFile(char const *fname);
	Real LogPriorProbability(Real const *TpSample, Real const *XpSample, int nsample) const;
	Real LogPosteriorProbability(Real const *par, Real *val, int ndim, int nvalue, int kwalker) const;

private:
  Inversion *pmy_invt_;
  Real Tstd_, Tlen_, Xstd_, Xlen_, chi_;
  bool fit_differential_;
};

void update_atm_profiles(MeshBlock *pmb, int k,
    Real const *PrSample, Real const *TpSample, Real const *XpSample, int nsample,
		std::vector<int> const& ix, Real Tstd, Real Tlen, Real Xstd, Real Xlen, Real chi = 1.);

void calculate_fit_target(MeshBlock *pmb, Real *val, int nvalue,
    int k, int j, bool differential = false);

#endif
