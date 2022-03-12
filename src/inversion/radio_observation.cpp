// C/C++ headers
#include <iomanip>
#include <iostream>
#include <cstring>
#include <cassert>

// Athena++ headers
#include "../mesh/mesh.hpp"
#include "../utils/utils.hpp"
#include "../debugger/debugger.hpp"
#include "../radiation/radiation.hpp"
#include "inversion.hpp"
//#include "../math/core.h"
#include "radio_observation.hpp"

RadioObservation::RadioObservation(Inversion *pinvt, ParameterInput *pin):
  pmy_invt_(pinvt)
{
  std::stringstream &msg = pinvt->pmy_block->pdebug->msg;
  pinvt->pmy_block->pdebug->Enter("RadioObservation");
  std::string obsfile = pin->GetOrAddString("inversion", "obsfile", "none");

  if (obsfile != "none") {
    ReadObservationFile(obsfile.c_str());
    msg << "- target: " << target.transpose() << std::endl
        << "- inverse covariance matrix" << std::endl
        << icov << std::endl;
  }

  // T correlation 
  Tstd_ = pin->GetReal("inversion", "Tstd");
  Tlen_ = pin->GetReal("inversion", "Tlen")*1.E3; // km -> m

  // X correlation
  Xstd_ = pin->GetReal("inversion", "Xstd")/1.E3;   // g/kg -> kg/kg
  Xlen_ = pin->GetReal("inversion", "Xlen")*1.E3;   // km -> m

  // power law coefficient
  chi_ = pin->GetOrAddReal("inversion", "chi", 0.1);

  // composition id
  ix = Vectorize<int>(pin->GetString("inversion", "Variables").c_str());

  // Pressure sample
  plevel = Vectorize<Real>(pin->GetString("inversion", "PrSample").c_str());
  int nsample = plevel.size();
  ndim_ = ix.size()*nsample;
  msg << "- number of input dimension = " << ndim_ << std::endl;
  msg << "- inversion pressure levels (bars) = ";
  for (std::vector<Real>::iterator m = plevel.begin(); m != plevel.end(); ++m)
    msg << *m << " ";
  msg << std::endl;
  
  // add boundaries
  Real pmax = pin->GetReal("inversion", "Pmax");
  Real pmin = pin->GetReal("inversion", "Pmin");
  if (pmax < (plevel.front()+1.E-6) || pmin > (plevel.back()-1.E-6)) {
    msg << "### FATAL ERROR in RadioObservation::RadioObservation" << std::endl
        << "Pmax (" << pmax << ")" << " must be greater than the largest value of PrSample" << std::endl
        << "Pmin (" << pmin << ")" << " must be lesser than the smallest value of PrSample";
    ATHENA_ERROR(msg);
  }
  plevel.insert(plevel.begin(), pmax);
  plevel.push_back(pmin);
  msg << "- top boundary = " << pmin << std::endl;
  msg << "- bottom boundary = " << pmax << std::endl;

  for (std::vector<Real>::iterator m = plevel.begin(); m != plevel.end(); ++m)
    (*m) *= 1.E5; // bar -> pa

  // output dimension
  nvalue_ = 3*pinvt->pmy_block->prad->GetNumBands();
  msg << "- number of output dimension = " << nvalue_ << std::endl;

  // initialize random positions
  srand(time(NULL) + Globals::my_rank);
  nwalker_ = pinvt->pmy_block->block_size.nx3;
  msg << "- walkers per block = " << nwalker_ << std::endl;
  msg << "- total number of walkers = " <<  pinvt->pmy_block->pmy_mesh->mesh_size.nx3
      << std::endl;
  if (nwalker_ < 2) {
    msg << "### FATAL ERROR in RadioObservation::RadioObservation"
        << "nwalker (nx3) must be at least " << 2;
    ATHENA_ERROR(msg);
  }

  NewCArray(init_pos_, nwalker_, ndim_);
  for (int n = 0; n < nwalker_; ++n) {
    int ic = 0;
    if (std::find(ix.begin(), ix.end(), 0) != ix.end()) {
      for (int i = 0; i < nsample; ++i)
        init_pos_[n][i] = (1.*rand()/RAND_MAX - 0.5)*Tstd_;
      ic = 1;
    }

    for (std::vector<int>::iterator m = ix.begin(); m != ix.end(); ++m)
      if (*m != 0) {
        for (int i = 0; i < nsample; ++i)
          init_pos_[n][ic*nsample + i] = (1.*rand()/RAND_MAX - 0.5)*Xstd_;
        ic++;
      }
  }

  // fit differential
  fit_differential_ = pin->GetOrAddBoolean("inversion", "differential", false);
  if (fit_differential_)
    msg << "- fit differential" << std::endl;
  pinvt->pmy_block->pdebug->Leave();
}

RadioObservation::~RadioObservation()
{
  FreeCArray(init_pos_);
}

#define MAX_LINE 512
void RadioObservation::ReadObservationFile(char const *fname)
{
  std::stringstream msg;
  FILE *fp = fopen(fname, "r");
  if (fp == NULL) {
    msg << "### FATAL ERROR in RadioObservation::ReadObseravtionFile" << std::endl 
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
