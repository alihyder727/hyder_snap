/** @file pi_log_prior_probability.cpp
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
#include "profile_inversion.hpp"
#include "gaussian_process.hpp"

Real ProfileInversion::LogPriorProbability(Real const *PrSample,
    Real const* const* XpSample, int nsample, Real const *Xstd, Real const *Xlen, Real chi) const
{
  std::stringstream &msg = pmy_block->pdebug->msg;
  Hydro *phydro = pmy_block->phydro;

  Real *zlev = new Real [nsample];
  Real *zstd = new Real [nsample];

  Real P0 = phydro->reference_pressure;
  Real H0 = phydro->scale_height;

  for (int i = 0; i < nsample; ++i)
    zlev[i] = -H0*log(PrSample[i]/P0);

	Real lnprior = 0.;
	for (int n = 0; n <= NVAPOR; ++n) {
		for (int i = 0; i < nsample; ++i)
			zstd[i] = Xstd[n]*pow(exp(zlev[i]/H0), chi);
		lnprior += gp_lnprior(SquaredExponential, XpSample[n], zlev, zstd, nsample, Xlen[n]);
	}

  msg << "- Log prior probability = " << lnprior << std::endl;

  delete[] zlev;
  delete[] zstd;

	return lnprior;
}
