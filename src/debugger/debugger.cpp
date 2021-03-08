// C/C++ headers
#include <cstring>
#include <sstream>
#include <algorithm>

// Athena++
#include "debugger.hpp"
#include "../globals.hpp"
#include "../coordinates/coordinates.hpp"

Debugger::Debugger(MeshBlock *pmb):
  pmy_block(pmb), prev(nullptr), next(nullptr), fname_("NULL")
{
  data_.NewAthenaArray(2*NHYDRO, pmb->ncells3, pmb->ncells2, pmb->ncells1);
}

Debugger::~Debugger()
{
  if (prev != nullptr) prev->next = next;
  if (next != nullptr) next->prev = prev;
}

Debugger* Debugger::StartTracking(std::string name)
{
  Debugger *my = this;
  while (my->next != nullptr) my = my->next;
  my->next = new Debugger(pmy_block); 
  my->next->fname_ = name;
  my->next->prev = my;
  my->next->next = nullptr;

  return my->next;
}

void Debugger::Track3D(std::string name, TestFunc_t test, AthenaArray<Real>& var, int n)
{
  std::stringstream msg;
  bool pass = true;
  MeshBlock *pmb = pmy_block;

  int is = pmb->is, js = pmb->js, ks = pmb->ks;
  int ie = pmb->ie, je = pmb->je, ke = pmb->ke;
  int c1, c2, c3;

  is -= NGHOST; ie += NGHOST;
  if (pmb->pmy_mesh->f2) {
    js -= NGHOST; je += NGHOST;
  }
  if (pmb->pmy_mesh->f3) {
    ks -= NGHOST; ke += NGHOST;
  }

  vnames_.push_back(name);

  if (vnames_.size() >= 2*NHYDRO) {
    msg << "### FATAL ERROR in function [Debugger::Track3D]"
        << std::endl << "Maximum number (" 
        << 2*NHYDRO << ") of trakcers reached." << std::endl;
    ATHENA_ERROR(msg);
  }

  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        data_(vnames_.size() - 1,k,j,i) = var(n,k,j,i);
        bool inside = (i >= pmb->is) && (i <= pmb->ie) && (j >= pmb->js) && (j <= pmb->je)
                   && (k >= pmb->ks) && (k <= pmb->ke);
        if (!test(var(n,k,j,i)) && pass && inside) {
          pass = false;
          c1 = i;
          c2 = j;
          c3 = k;
        }
      }

  if (!pass) {
    DumpTracking(name, c1, c2, c3, "w");
    msg << "### FATAL ERROR in function [Debugger::Track3D]"
        << std::endl << "Debugger catched invalid " << name << " at "
        << "(k=" << c3 << ", " << "j=" << c2 << ", " << "i=" << c1 << ") "
        << "on rank " << Globals::my_rank << std::endl;
    ATHENA_ERROR(msg);
  }
}

void Debugger::Track1D(std::string name, TestFunc_t test, AthenaArray<Real>& var, 
  int n, int k, int j)
{
  std::stringstream msg;
  bool pass = true;
  MeshBlock *pmb = pmy_block;

  int is = pmb->is - NGHOST, ie = pmb->ie + NGHOST;
  int c1, c2, c3;

  std::vector<std::string>::iterator it = std::find(vnames_.begin(), vnames_.end(), name);
  if (it == vnames_.end())
    vnames_.push_back(name);
  int m = std::find(vnames_.begin(), vnames_.end(), name) - vnames_.begin();

  if (m >= 2*NHYDRO) {
    msg << "### FATAL ERROR in function [Debugger::Track1D]"
        << std::endl << "Maximum number (" 
        << 2*NHYDRO << ") of trakcers reached." << std::endl;
    ATHENA_ERROR(msg);
  }

  for (int i = is; i <= ie; ++i) {
    data_(m,k,j,i) = var(n,i);
    bool inside = (i >= pmb->is) && (i <= pmb->ie);
    if (!test(var(n,i)) && pass && inside) {
      pass = false;
      c1 = i;
      c2 = j;
      c3 = k;
    }
  }

  if (!pass) {
    DumpTracking(name, c1, c2, c3, "w");
    msg << "### FATAL ERROR in function [Debugger::Track1D]"
        << std::endl << "Debugger catched invalid " << name << " at "
        << "(k=" << c3 << ", " << "j=" << c2 << ", " << "i=" << c1 << ") "
        << "on rank " << Globals::my_rank << std::endl;
    ATHENA_ERROR(msg);
  }
}

void Debugger::DumpTracking(std::string name, int c1, int c2, int c3, char const* mode)
{
  MeshBlock *pmb = pmy_block;
  Coordinates *pcoord = pmb->pcoord;

  int is = c1 - NGHOST, ie = c1 + NGHOST;
  int js, je, ks, ke;
  if (pmb->pmy_mesh->f2) {
    js = c2 - NGHOST;
    je = c2 + NGHOST;
  } else
    js = je = 0;
  if (pmb->pmy_mesh->f3) {
    ks = c3 - NGHOST;
    ke = c3 + NGHOST;
  } else
    ks = ke = 0;

  // loop over all variables
  for (int n = 0; n < vnames_.size(); ++n) {
    // open file for output
    std::string out;
    char blockid[12];
    std::snprintf(blockid, sizeof(blockid), "block%d", pmb->gid);
    out.assign(vnames_[n]);
    out.append(".");
    out.append(blockid);
    out.append(".txt");
    
    FILE *pfile;
    std::stringstream msg;
    if ((pfile = std::fopen(out.c_str(), mode)) == nullptr) {
      msg << "### FATAL ERROR in function [Debugger::DumpTracking]"
          << std::endl << "Output file '" << out << "' could not be opened"
          << std::endl;
      ATHENA_ERROR(msg);
    }

    // file header
    std::fprintf(pfile, "# Debugger dump at time = %e\n", pmb->pmy_mesh->time);
    std::fprintf(pfile, "# Variable = %s\n", vnames_[n].c_str());
    std::fprintf(pfile, "# Subdomain extent = [%.3g, %.3g] x [%.3g, %.3g] x [%.3g, %.3g]\n",
      pcoord->x1f(pmb->is), pcoord->x1f(pmb->ie+1),
      pcoord->x2f(pmb->js), pcoord->x2f(pmb->je+1),
      pcoord->x3f(pmb->ks), pcoord->x3f(pmb->ke+1));
    std::fprintf(pfile, "# Break at (k=%d, j=%d, i=%d) of variable %s\n", c3, c2, c1,
      name.c_str());

    // data table
    std::fprintf(pfile, "# Currently at [%s]\n", fname_.c_str());
    std::fprintf(pfile, "# Previously at [%s]\n", prev->fname_.c_str());
    for (int k = ks; k <= ke; ++k) {
      std::fprintf(pfile, "k=%02d        ", k);
      for (int j = js; j <= je; ++j)
        std::fprintf(pfile, " (k=%02d,j=%02d)", k, j);
      std::fprintf(pfile, "\n");

      for (int i = is; i <= ie; ++i) {
        std::fprintf(pfile, "i=%02d        ", i);
        for (int j = js; j <= je; ++j)
          std::fprintf(pfile, "%12.3g", data_(n,k,j,i));
        std::fprintf(pfile, "\n");
      }
      std::fprintf(pfile, "\n");
    }

    fclose(pfile);
    if (prev != nullptr) prev->DumpTracking(name, c1, c2, c3, "a");
  }
}
