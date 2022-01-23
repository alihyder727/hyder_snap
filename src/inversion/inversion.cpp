// C/C++ headers
#include <iostream>
#include <sstream>
#include <string>
#include <cstring>
#include <stdexcept>

// Athena++ headers
#include "../parameter_input.hpp"
#include "../mesh/mesh.hpp"
#include "../debugger/debugger.hpp"
#include "inversion.hpp"
#include "radio_observation.hpp"
#include "mcmc_impl.hpp"

Inversion::Inversion(MeshBlock *pmb, ParameterInput *pin):
  pmy_block(pmb), pradio(nullptr), mcmc_initialized_(false)
{
  pmb->pdebug->Enter("Inversion");
  std::stringstream &msg = pmb->pdebug->msg;
  task = pin->GetOrAddString("inversion", "task", "none");

  opts_.a = pin->GetOrAddReal("inversion", "stretch", 2.);
  opts_.p = pin->GetOrAddInteger("inversion", "walk", 4);
  opts_.print = pin->GetOrAddInteger("inversion", "print", 100);
#ifdef MPI_PARALLEL
  opts_.mpi_comm = MPI_COMM_WORLD;
#endif

  strcpy(opts_.logfile, pin->GetOrAddString("inversion", "logfile", "inversion.log").c_str());

  if (task != "none") {
    msg << "- number of MCMC steps = " << pmb->pmy_mesh->nlim << std::endl;
    if (task == "radio") {
      pradio = new RadioObservation(this, pin);
      mcmc_alloc(&recs_, pmb->pmy_mesh->nlim+1, pradio->nwalker_, pradio->ndim_, pradio->nvalue_);
      mcmc_initialized_ = true;
    } else {
      msg << "### FATAL ERROR in function Inversion::Inversion"
          << std::endl << "Unrecognized inversion task.";
      ATHENA_ERROR(msg);
    }
  }
  pmb->pdebug->Leave();
}

Inversion::~Inversion()
{
  if (pradio != nullptr) delete pradio;
  if (mcmc_initialized_) mcmc_free(&recs_);
}

void Inversion::MakeMCMCOutputs(std::string fname)
{
  // initialize model
  if (recs_.cur == 0)
    mcmc_init(pradio, pradio->init_pos_, &opts_, &recs_);

  std::stringstream msg;
  if (!mcmc_initialized_) {
    msg << "### FATAL ERROR in function Inversion::MakeMCMCOutputs"
        << std::endl << "mcmc chain uninitialized.";
    ATHENA_ERROR(msg);
  }
  mcmc_save_fits(fname.c_str(), &opts_, &recs_);
}

void Inversion::MCMCStep()
{
  // initialize model
  if (recs_.cur == 0)
    mcmc_init(pradio, pradio->init_pos_, &opts_, &recs_);

  std::stringstream msg;
  if (!mcmc_initialized_) {
    msg << "### FATAL ERROR in function Inversion::MCMCStep"
        << std::endl << "mcmc chain uninitialized.";
    ATHENA_ERROR(msg);
  }
  mcmc_advance(pradio, &opts_, &recs_);
}

void Inversion::ResetChain()
{
  std::stringstream msg;
  if (!mcmc_initialized_) {
    msg << "### FATAL ERROR in function Inversion::ResetChain"
        << std::endl << "mcmc chain uninitialized.";
    ATHENA_ERROR(msg);
  }

  int cur = recs_.cur;
  // copy the last state into the first state
  for (int k = 0; k < recs_.nwalker; ++k) {
    for (int d = 0; d < recs_.ndim; ++d)
      recs_.par[0][k][d] = recs_.par[cur-1][k][d];
    for (int d = 0; d < recs_.nvalue; ++d)
      recs_.val[0][k][d] = recs_.val[cur-1][k][d];
    recs_.lnp[0][k] = recs_.lnp[cur-1][k];
    recs_.newstate[0][k] = recs_.newstate[cur-1][k];
  }
  recs_.reset += cur - 1;
  recs_.cur = 1;
}
