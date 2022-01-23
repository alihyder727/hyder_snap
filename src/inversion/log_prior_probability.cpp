/** @file log_prior_probability.cpp
 * @brief
 *
 * @author Cheng Li (chengcli@umich.edu)
 * @date Saturday Nov 20, 2021 09:21:41 EST
 * @bug No known bugs.
 */

// C/C++ header
#include <iostream>

// Athena++ header
#include "../mesh/mesh.hpp"
#include "../hydro/hydro.hpp"
#include "../debugger/debugger.hpp"
#include "inversion.hpp"
#include "radio_observation.hpp"
#include "gaussian_process.hpp"

Real RadioObservation::LogPriorProbability(Real const *TpSample, Real const *XpSample, int nsample) const
{
  std::stringstream &msg = pmy_invt_->pmy_block->pdebug->msg;
  Hydro *phydro = pmy_invt_->pmy_block->phydro;
  Real *zlev = new Real [nsample];
  Real *stdSample = new Real [nsample];
  Real P0 = phydro->reference_pressure;
  Real H0 = phydro->scale_height;

  for (int i = 0; i < nsample; ++i)
    zlev[i] = -H0*log(plevel[i]/P0);

  for (int i = 0; i < nsample; ++i)
    stdSample[i] = Tstd_*pow(exp(zlev[i]/H0), chi_);
  Real lnp1 = gp_lnprior(SquaredExponential, TpSample, zlev, stdSample, nsample, Tlen_);

  for (int i = 0; i < nsample; ++i)
    stdSample[i] = Xstd_*pow(exp(zlev[i]/H0), chi_);
  Real lnp2 = gp_lnprior(SquaredExponential, XpSample, zlev, stdSample, nsample, Xlen_);

  msg << "- Log prior probability = " << lnp1 + lnp2 << std::endl;

  delete[] zlev;
  delete[] stdSample;

	return lnp1 + lnp2;
}

