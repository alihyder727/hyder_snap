// C/C++ headers
#include <iomanip>
#include <iostream>
#include <cstring>
#include <cassert>

// Athena++ headers
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../utils/utils.hpp"
#include "../debugger/debugger.hpp"
#include "../radiation/radiation.hpp"
#include "profile_inversion.hpp"

ProfileInversion::ProfileInversion(MeshBlock *pmb, ParameterInput *pin):
  Inversion(pmb, pin)
{
  std::stringstream &msg = pmb->pdebug->msg;
  pmb->pdebug->Enter("ProfileInversion");

  // T correlation 
  Xstd_[0] = pin->GetReal("inversion", "tem.std");
  Xlen_[0] = pin->GetReal("inversion", "tem.corr.km")*1.E3; // km -> m
  msg << "- temperature std = " << Xstd_[0] << std::endl;
  msg << "- temperature correlation length = " << Xlen_[0] << std::endl;

  for (int n = 1; n <= NVAPOR; ++n) {
    Xstd_[n] = pin->GetOrAddReal("inversion", "qvapor" + std::to_string(n) + ".std.gkg", 1.)/1.E3;
    // km -> m
    Xlen_[n] = pin->GetOrAddReal("inversion", "qvapor" + std::to_string(n) + ".corr.km", 1.)*1.E3;
    msg << "- vapor " << n << " std = "  << Xstd_[n] << std::endl;
    msg << "- vapor " << n << " correlation length = " << Xlen_[n] << std::endl;
  }

  // power law coefficient
  chi_ = pin->GetOrAddReal("inversion", "chi", 0.0);

  // composition id
  ix = Vectorize<int>(pin->GetString("inversion", "Variables").c_str());

  // Pressure sample
  plevel = Vectorize<Real>(pin->GetString("inversion", "PrSample").c_str());
  int nsample = plevel.size();

  ndim = ix.size()*nsample;
  msg << "- number of input dimension = " << ndim << std::endl;
  msg << "- inversion pressure levels (bars) = ";
  for (std::vector<Real>::iterator m = plevel.begin(); m != plevel.end(); ++m)
    msg << *m << " ";
  msg << std::endl;
  
  // add boundaries
  Real pmax = pin->GetReal("inversion", "Pmax");
  Real pmin = pin->GetReal("inversion", "Pmin");
  if (pmax < (plevel.front()+1.E-6) || pmin > (plevel.back()-1.E-6)) {
    msg << "### FATAL ERROR in ProfileInversion::ProfileInversion" << std::endl
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

  // fit differential
  fit_differential_ = pin->GetOrAddBoolean("inversion", "differential", false);
  if (fit_differential_)
    msg << "- fit differential" << std::endl;

  // output dimension
  //nvalue = 2*pmb->prad->getNumBands();
  nvalue = pmb->prad->getNumBands();
  msg << "- number of output dimension = " << nvalue << std::endl;

  // number of walkers
  nwalker = pmb->block_size.nx3;
  msg << "- walkers per block = " << nwalker << std::endl;
  msg << "- total number of walkers = " <<  pmb->pmy_mesh->mesh_size.nx3
      << std::endl;
  if ((nwalker < 2) && pmb->pmy_mesh->nlim > 0) {
    msg << "### FATAL ERROR in ProfileInversion::ProfileInversion"
        << "nwalker (nx3) must be at least " << 2;
    ATHENA_ERROR(msg);
  }

  // initialize mcmc chain
  InitializeChain(nwalker, ndim, nvalue);

  // initialize random positions
  msg << "- initialize random positions for walkers" << std::endl;
  srand(time(NULL) + Globals::my_rank);
  NewCArray(init_pos_, nwalker, ndim);
  for (int n = 0; n < nwalker; ++n) {
    int ip = 0;
    for (std::vector<int>::iterator m = ix.begin(); m != ix.end(); ++m, ++ip)
      for (int i = 0; i < nsample; ++i)
        init_pos_[n][ip*nsample + i] = (1.*rand()/RAND_MAX - 0.5)*Xstd_[*m]
          *pow(pmb->phydro->reference_pressure/plevel[i+1], chi_);
  }

  pmb->pdebug->Leave();
}
