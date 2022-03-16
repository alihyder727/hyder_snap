/** @file profile_inversion.hpp
 * @brief
 *
 * @author Cheng Li (chengcli@umich.edu)
 * @date Wednesday Mar 16, 2022 07:10:49 EDT
 * @bug No known bugs.
 */

#ifndef PROFILE_INVERSION_HPP
#define PROFILE_INVERSION_HPP

// Athena++ header
#include "inversion.hpp"

class ProfileInversion : public Inversion {
public:
  // functions
  ProfileInversion(MeshBlock *pmb, ParameterInput *pin);
  ~ProfileInversion() {}

	Real LogPriorProbability(Real const *PrSample, Real const* const* XpSample, int nsample,
		Real const *Xstd, Real const *Xlen, Real chi) const;

	Real LogPosteriorProbability(Real const *par, Real *val, int ndim, int nvalue, int kwalker) const;

	void UpdateAtmosphere(void *params, int k, int j = -1) const {
    Real const* const* *XpSample = static_cast<Real const* const* *>(params);
		UpdateProfiles_(k, plevel.data(), (*XpSample), plevel.size(), Xstd_, Xlen_, chi_);
	}

protected:
	// atmospheric functions
	void UpdateProfiles_(int k, Real const *PrSample,
		Real const* const* XpSample, int nsample, Real const *Xstd, Real const *Xlen, Real chi) const;

	// hyper-parameters
	Real chi_;
	Real Xstd_[1+NVAPOR];
	Real Xlen_[1+NVAPOR];

private:
  bool fit_differential_;
};


#endif /* end of include guard PROFILE_INVERSION_HPP */

