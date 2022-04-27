/** @file pi_log_posterior_probability.cpp
 * @brief
 *
 * @author Cheng Li (chengcli@umich.edu)
 * @date Saturday Nov 20, 2021 09:00:13 EST
 * @bug No known bugs.
 */

// C/C++ header
#include <iostream>
#include <sstream>

// Athena++ header
#include "../mesh/mesh.hpp"
#include "../hydro/hydro.hpp"
#include "../radiation/radiation.hpp"
#include "../debugger/debugger.hpp"
#include "../utils/utils.hpp"
#include "profile_inversion.hpp"

Real ProfileInversion::LogPosteriorProbability(Real const *par, Real *val,
	int ndim_, int nvalue_, int kw) const
{
  MeshBlock *pmb = pmy_block;
  Hydro *phydro = pmb->phydro;
  Radiation *prad = pmb->prad;
  int is = pmb->is, js = pmb->js, ks = pmb->ks;
  int ie = pmb->ie, je = pmb->je, ke = pmb->ke;
  std::stringstream &msg = pmb->pdebug->msg;

  // check parameter consistency
  if (ndim != ndim_ || nvalue != nvalue_) {
    msg << "### FATAL ERROR in function ProfileInversion::LogPosteriorProbability"
        << std::endl << "input dimension inconsistent";
    ATHENA_ERROR(msg);
  }

  if (ndim % ix.size() != 0) {
    msg << "### FATAL ERROR in function ProfileInversion::LogPosteriorProbability"
        << std::endl << "inversion dimension (ndim) cannot be divided by number of "
        << "inversion variables";
    ATHENA_ERROR(msg);
  }

  if (ndim/ix.size() != plevel.size() - 2) {
    msg << "### FATAL ERROR in function ProfileInversion::LogPosteriorProbability"
        << std::endl << "inversion dimension (ndim) and inversion levels do not match";
    ATHENA_ERROR(msg);
  }

  // logging
  pmb->pdebug->Call("LogPosteriorProbability");
  msg << "- I am walker " << kw << std::endl;

	Real **XpSample;
	NewCArray(XpSample, 1+NVAPOR, plevel.size());
  std::fill(*XpSample, *XpSample + (1+NVAPOR)*plevel.size(), 0.);

  // sample temperature, sample composition #1, sample composition #2, ...
  msg << "- parameters: ";
  for (int i = 0; i < ndim; ++i)
    msg << par[i] << " ";
  msg << std::endl;

  int ip = 0;
  for (std::vector<int>::const_iterator m = ix.begin(); m != ix.end(); ++m, ++ip) {
		XpSample[*m][0] = 0.;
		for (int i = 1; i <= plevel.size() - 2; ++i)
			XpSample[*m][i] = par[ip*(plevel.size()-2)+i-1];
		XpSample[*m][plevel.size() - 1] = 0.;
  }

  // update atmosphere based on XpSample
  UpdateProfiles_(ks+kw, plevel.data(), XpSample, plevel.size(), Xstd_, Xlen_, chi_);

  // calculate radiation for updated profiles located at j = js+1 ... je
  msg << "- run RT for models 1 to " << je - js << std::endl;
  for (int j = js+1; j <= je; ++j)
    prad->calculateRadiances(phydro->w, 0., ks+kw, j, is, ie+1);

  // prior probability
  Real lnprior = LogPriorProbability(plevel.data(), XpSample, plevel.size(), Xstd_, Xlen_, chi_);

  // posterior probability
  Eigen::VectorXd misfit(nvalue);

  // calculate model result for profile at j = je
  msg << "- calculate output for model " << je - js << std::endl;
  calculate_fit_target(pmb, val, nvalue, ks+kw, je, fit_differential_);
  Real lnpost = 0.;
  if (target.size() > 0) {
    for (int m = 0; m < nvalue; ++m)
      misfit(m) = val[m] - target(m);
    lnpost = -0.5*misfit.transpose()*icov*misfit;
  }

  // posterior probability
  msg << "- log posterir probability = " << lnpost << std::endl;
  pmb->pdebug->Leave();

	FreeCArray(XpSample);

  return lnprior + lnpost;
}

