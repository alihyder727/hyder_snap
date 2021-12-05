/** @file log_posterior_probability.cpp
 * @brief
 *
 * @author Cheng Li (chengcli@umich.edu)
 * @date Saturday Nov 20, 2021 09:00:13 EST
 * @bug No known bugs.
 */

// C/C++ header
#include <iostream>

// Athena++ header
#include "../mesh/mesh.hpp"
#include "../hydro/hydro.hpp"
#include "../radiation/radiation.hpp"
#include "inversion.hpp"
#include "radio_observation.hpp"

Real RadioObservation::LogPosteriorProbability(Real const *par, Real *val, int ndim, int nvalue) const
{
  MeshBlock *pmb = pmy_invt_->pmy_block;
  Hydro *phydro = pmb->phydro;
  Radiation *prad = pmb->prad;
  int is = pmb->is, js = pmb->js, ks = pmb->ks;
  int ie = pmb->ie, je = pmb->je, ke = pmb->ke;
  int nsample = ndim/ix.size();

  Real *TpSample = new Real [nsample+2];
  Real *XpSample = new Real [ndim+2];
	std::fill(TpSample, TpSample + nsample + 2, 0.);
	std::fill(XpSample, XpSample + ndim + 2, 0.);

  // sample temperature, sample composition #1, sample composition #2, ...
  std::cout << "- Parameters: ";
	for (int i = 0; i < ndim; ++i)
		std::cout << par[i] << " ";
  std::cout << std::endl;

	int ic = 0;
	if (std::find(ix.begin(), ix.end(), 0) != ix.end()) {
    TpSample[0] = 0.;
		for (int i = 1; i <= nsample; ++i)
			TpSample[i] = par[i];
    TpSample[nsample+1] = 0.;
		ic = 1;
	}
	for (std::vector<int>::const_iterator m = ix.begin(); m != ix.end(); ++m) {
		if (*m != 0) {
      XpSample[0] = 0.;
			for (int i = 1; i <= nsample; ++i)
				XpSample[i] = par[ic*nsample+i];
      XpSample[nsample+1] = 0.;
			ic++;
		}
  }

	// update atmosphere based on TpSample and XpSample
  update_atm_profiles(pmb, plevel.data(), TpSample, XpSample, nsample+2,
      ix, Tstd_, Tlen_, Xstd_, Xlen_, chi_);

  // calculate radiation for updated profiles located at j = je
	prad->CalculateRadiances(phydro->w, 0., ks, je, is, ie+1);

	// prior probability
	Real lnprior = LogPriorProbability(TpSample, XpSample, nsample+2);

  // posterior probability
  Eigen::VectorXd misfit(nvalue);

  // calculate model result for profile at j = je
  calculate_fit_target(pmb, val, nvalue, je, fit_differential_);
  for (int m = 0; m < nvalue; ++m)
    misfit(m) = val[m] - target(m);
  Real lnpost = -0.5*misfit.transpose()*icov*misfit;

  // posterior probability
  std::cout << "- Log posterir probability: ";
  std::cout << lnpost << std::endl;

  delete[] TpSample;
  delete[] XpSample;

  return lnprior + lnpost;
}

