// C/C++ headers
#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <cassert>

// Athena++ headers
#include "../math/core.h"
#include "../utils/utils.hpp"
#include "mcmc.hpp"

template<typename T>
void mcmc_init(
  T const* pobj,
  double **par,
  mcmc_opts *opts,
  mcmc_recs *recs)
{
  int nwalker = recs->nwalker;
  int ndim = recs->ndim;
  int nval = recs->nvalue;

  // make sure that the start points are valid
  for (int k = 0; k < nwalker; ++k) {
    recs->lnp[0][k] = pobj->LogPosteriorProbability(par[k], recs->val[0][k], ndim, nval);

    int niter = 0;
    while (std::isnan(recs->lnp[0][k]) && (niter++ < 10)) {  // if point (k) is invalid
      mcmc_stretch_move(par[k], par, k, nwalker, ndim, opts);
      recs->lnp[0][k] = pobj->LogPosteriorProbability(par[k], recs->val[0][k], ndim, nval);
    }

    if (niter >= 10) {
      std::cerr << "Starting point iteration > 10 times" << std::endl;
      std::exit(1);
    }

    // transfer input parameters to records
    for (int d = 0; d < ndim; ++d)
      recs->par[0][k][d] = par[k][d];

    if (recs->lnp[0][k] > recs->opt_lnp) {
      recs->opt_lnp = recs->lnp[0][k];
      for (int d = 0; d < ndim; ++d)
        recs->opt_par[d] = recs->par[0][k][d];
      for (int d = 0; d < nval; ++d)
        recs->opt_val[d] = recs->val[0][k][d];
    }
    recs->newstate[0][k] = 1;
  }

  recs->cur++;
  recs->accept += nwalker;

  // make initial output
  mcmc_report(opts, recs, "w");
}

template<typename T>
void mcmc_advance(
	T const *pobj,
  mcmc_opts *opts,
  mcmc_recs *recs)
{
  int cur = recs->cur;
  int ndim = recs->ndim;
  int nwalker = recs->nwalker;
  int nval = recs->nvalue;

  double *par = new double [ndim];

  for (int k = 0; k < nwalker; ++k) {
    double zz = mcmc_stretch_move(par, recs->par[cur-1], k, nwalker, ndim, opts);
    //mcmc_walk_move(par_all + k*np, par, np, k, nwalker, zz + k, opts);

    recs->lnp[cur][k] = pobj->LogPosteriorProbability(par, recs->val[cur][k], ndim, nval);

    double lnp0 = recs->lnp[cur-1][k],
           lnp1 = recs->lnp[cur][k],
           pdiff = (ndim - 1.)*log(zz) + lnp1 - lnp0;

    if (pdiff > log(1.*rand()/RAND_MAX)) {  // accept this position
      for (int d = 0; d < ndim; ++d)
        recs->par[cur][k][d] = par[d];

      if (lnp1 > recs->opt_lnp) { // new best position
        recs->opt_lnp = lnp1;
        for (int d = 0; d < ndim; ++d)
          recs->opt_par[d] = recs->par[cur][k][d];
        for (int d = 0; d < nval; ++d)
          recs->opt_val[d] = recs->val[cur][k][d];
      }

      recs->accept++;
      recs->newstate[cur][k] = 1;
    } else {  // do not accept this position
      for (int d = 0; d < ndim; ++d)
        recs->par[cur][k][d] = recs->par[cur-1][k][d];
      for (int d = 0; d < nval; ++d)
        recs->val[cur][k][d] = recs->val[cur-1][k][d];
      recs->lnp[cur][k] = recs->lnp[cur-1][k];
      recs->newstate[cur][k] = 0;
    }
  }

  recs->cur++;
  if (recs->cur % opts->print == 0)
    mcmc_report(opts, recs, "a");
  delete[] par;
}
