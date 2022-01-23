/** @file log_posterior_probability.cpp
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
#include "inversion.hpp"
#include "radio_observation.hpp"

Real RadioObservation::LogPosteriorProbability(Real const *par, Real *val, int ndim, int nvalue, int kw) const
{
  MeshBlock *pmb = pmy_invt_->pmy_block;
  Hydro *phydro = pmb->phydro;
  Radiation *prad = pmb->prad;
  int is = pmb->is, js = pmb->js, ks = pmb->ks;
  int ie = pmb->ie, je = pmb->je, ke = pmb->ke;
  std::stringstream &msg = pmb->pdebug->msg;

  // check parameter consistency
  if (ndim != ndim_ || nvalue != nvalue_) {
    msg << "### FATAL ERROR in function RadioObservation::LogPosteriorProbability"
        << std::endl << "input dimension inconsistent";
    ATHENA_ERROR(msg);
  }

  if (ndim % ix.size() != 0) {
    msg << "### FATAL ERROR in function RadioObservation::LogPosteriorProbability"
        << std::endl << "inversion dimension (ndim) cannot be divided by number of "
        << "inversion variables";
    ATHENA_ERROR(msg);
  }

  // logging
  pmb->pdebug->Call("LogPosteriorProbability");
  msg << "- I am walker " << kw << std::endl;

  int nsample = ndim/ix.size();
  Real *TpSample = new Real [nsample+2];
  Real *XpSample = new Real [ix.size()*(nsample+2)];
  std::fill(TpSample, TpSample + (nsample+2), 0.);
  std::fill(XpSample, XpSample + ix.size()*(nsample+2), 0.);

  // sample temperature, sample composition #1, sample composition #2, ...
  msg << "- parameters: ";
  for (int i = 0; i < ndim; ++i)
    msg << par[i] << " ";
  msg << std::endl;

  int ic = 0;
  if (std::find(ix.begin(), ix.end(), 0) != ix.end()) {
    TpSample[0] = 0.;
    for (int i = 1; i <= nsample; ++i)
      TpSample[i] = par[i-1];
    TpSample[nsample+1] = 0.;
    ic = 1;
  }
  for (std::vector<int>::const_iterator m = ix.begin(); m != ix.end(); ++m) {
    if (*m != 0) {
      XpSample[ic*(nsample+2)] = 0.;
      for (int i = 1; i <= nsample; ++i)
        XpSample[ic*(nsample+2) + i] = par[ic*nsample+i-1];
      XpSample[ic*(nsample+2) + nsample+1] = 0.;
      ic++;
    }
  }

  // update atmosphere based on TpSample and XpSample
  update_atm_profiles(pmb, ks+kw, plevel.data(), TpSample, XpSample, nsample+2,
      ix, Tstd_, Tlen_, Xstd_, Xlen_, chi_);

  // calculate radiation for updated profiles located at j = js+1 ... je
  msg << "- run RT for models 1 to " << je - js << std::endl;
  for (int j = js+1; j <= je; ++j)
    prad->CalculateRadiances(phydro->w, 0., ks+kw, j, is, ie+1);

  // prior probability
  Real lnprior = LogPriorProbability(TpSample, XpSample, nsample+2);

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

  delete[] TpSample;
  delete[] XpSample;

  return lnprior + lnpost;
}

