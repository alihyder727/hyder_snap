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
#include "inversion.hpp"
#include "radio_observation.hpp"
#include "gaussian_process.hpp"

Real RadioObservation::LogPriorProbability(Real const *TpSample, Real const *XpSample, int nsample) const
{
  Hydro *phydro = pmy_invt_->pmy_block->phydro;
  Real *zlev = new Real [nsample];
  Real *stdSample = new Real [nsample];

  for (int i = 0; i < nsample; ++i)
    zlev[i] = phydro->reference_height - 
			phydro->scale_height*log(plevel[i]/phydro->reference_pressure);

  for (int i = 0; i < nsample; ++i)
    stdSample[i] = Tstd_;
  Real lnp1 = gp_lnprior(SquaredExponential, TpSample, zlev, stdSample, nsample, Tlen_);

  for (int i = 0; i < nsample; ++i)
    stdSample[i] = Xstd_;
  Real lnp2 = gp_lnprior(SquaredExponential, XpSample, zlev, stdSample, nsample, Xlen_);

  std::cout << "- Log prior probability: ";
  std::cout << lnp1 + lnp2 << std::endl;

  delete[] zlev;
  delete[] stdSample;

	return lnp1 + lnp2;
}

