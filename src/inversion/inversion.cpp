// C/C++ headers
#include <iostream>
#include <sstream>
#include <string>
#include <cstring>
#include <stdexcept>

// Athena++ headers
#include "../parameter_input.hpp"
#include "../mesh/mesh.hpp"
#include "inversion.hpp"

Inversion::Inversion(MeshBlock *pmb, ParameterInput *pin)
{
  pmy_block_ = pmb;
  method = pin->GetOrAddString("inversion", "method", "none");
  initialized_ = false;
  lnprob_ = NULL;
  obj = NULL;

  if (method == "mcmc") {
    opts_.a = pin->GetOrAddReal("inversion", "stretch", 2.);
    opts_.p = pin->GetOrAddInteger("inversion", "walk", 4);
    opts_.print = pin->GetOrAddInteger("inversion", "print", 100);
#ifdef MPI_PARALLEL
    opts_.mpi_comm = MPI_COMM_WORLD;
#endif

    strcpy(opts_.logfile, pin->GetOrAddString("inversion", "logfile", "inversion.log").c_str());
  }
}

Inversion::~Inversion()
{
  if (initialized_ && method == "mcmc") {
    mcmc_free(&recs_);
  }
}

void __attribute__((weak)) Inversion::Finish() {} // override in pgen

void Inversion::EnrollObjectives(ObjectiveFunction_t lnprob, void *myobj)
{
  lnprob_ = lnprob;
  obj = myobj;
}

void Inversion::Initialize(Real **pos, int nwalker, int ndim, int nvalue)
{
  std::stringstream msg;
  if (lnprob_ == NULL) {
    msg << "### FATAL ERROR in function Inversion::Initialize"
        << std::endl << "Enroll objective function first.";
    ATHENA_ERROR(msg);
  }

  if (method == "mcmc") {
    mcmc_alloc(&recs_, pmy_block_->pmy_mesh->nlim+1, nwalker, ndim, nvalue);
    mcmc_init(lnprob_, pos, &opts_, &recs_, obj);
  }

  initialized_ = true;
}

void Inversion::MakeMCMCOutputs(std::string fname)
{
  mcmc_save_fits(fname.c_str(), &opts_, &recs_);
}

void Inversion::MCMCStep(void *myobj)
{
  if (initialized_ && lnprob_ != NULL)
    mcmc_advance(lnprob_, &opts_, &recs_, myobj);
}

void Inversion::ResetChain()
{
  int nwalker = recs_.nwalker;
  int ndim = recs_.ndim;
  int nvalue = recs_.nvalue;
  int cur = recs_.cur;

  // copy the last state into the first state
  for (int k = 0; k < nwalker; ++k) {
    for (int d = 0; d < ndim; ++d)
      recs_.par[0][k][d] = recs_.par[cur-1][k][d];
    for (int d = 0; d < nvalue; ++d)
      recs_.val[0][k][d] = recs_.val[cur-1][k][d];
    recs_.lnp[0][k] = recs_.lnp[cur-1][k];
    recs_.newstate[0][k] = recs_.newstate[cur-1][k];
  }
  recs_.reset += cur - 1;
  recs_.cur = 1;
}
