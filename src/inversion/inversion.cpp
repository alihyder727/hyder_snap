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
#include "mcmc_impl.hpp"

Inversion::Inversion(MeshBlock *pmb, ParameterInput *pin):
  next(nullptr), prev(nullptr), pmy_block(pmb), mcmc_initialized_(false), init_pos_(nullptr)
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

  strcpy(opts_.logfile, pin->GetOrAddString("inversion", "logfile", task + ".log").c_str());

  std::string obsfile = pin->GetOrAddString("inversion", "obsfile", "none");
  if (obsfile != "none") {
    ReadObservationFile(obsfile.c_str());
    msg << "- target: " << target.transpose() << std::endl
        << "- inverse covariance matrix" << std::endl
        << icov << std::endl;
  }

  pmb->pdebug->Leave();
}

Inversion::~Inversion()
{
  if (mcmc_initialized_) mcmc_free(&recs_);
  if (init_pos_ != nullptr)
    FreeCArray(init_pos_);
}

void Inversion::InitializeChain(int nwalker, int ndim, int nvalue) {
	mcmc_alloc(&recs_, pmy_block->pmy_mesh->nlim+1, nwalker, ndim, nvalue);
	mcmc_initialized_ = true;
}

#define MAX_LINE 512
void Inversion::ReadObservationFile(std::string fname)
{
  std::stringstream msg;
  FILE *fp = fopen(fname.c_str(), "r");
  if (fp == NULL) {
    msg << "### FATAL ERROR in ProfileInversion::ReadObseravtionFile" << std::endl 
        << fname << " cannot be opened.";
    ATHENA_ERROR(msg);
  }
  char line[MAX_LINE], *pl;

  int rows;
  // header
  pl = NextLine(line, MAX_LINE, fp);

  // target values
  sscanf(pl, "%d", &rows);
  target.resize(rows);
  icov.resize(rows, rows);

  for (int i = 0; i < rows; ++i) {
    pl = NextLine(line, MAX_LINE, fp);
    sscanf(pl, "%lf", &target(i));
  }

  // inverse covariance matrix
  for (int i = 0; i < rows; ++i) {
    pl = NextLine(line, MAX_LINE, fp);
    char *p = strtok(pl, " ");
    for (int j = 0; j < rows; ++j) {
      sscanf(p, "%lf", &icov(i,j));
      p = strtok(NULL, " ");
    }
  }

  fclose(fp);
}
#undef MAX_LINE

void Inversion::MakeMCMCOutputs(std::string fname)
{
  // initialize model
  if (recs_.cur == 0)
    mcmc_init(this, init_pos_, &opts_, &recs_);

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
    mcmc_init(this, init_pos_, &opts_, &recs_);

  std::stringstream msg;
  if (!mcmc_initialized_) {
    msg << "### FATAL ERROR in function Inversion::MCMCStep"
        << std::endl << "mcmc chain uninitialized.";
    ATHENA_ERROR(msg);
  }
  mcmc_advance(this, &opts_, &recs_);
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
