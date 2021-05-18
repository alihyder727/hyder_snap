// C/C++ headers
#include <cstring>
#include <sstream>
#include <stdexcept>

// Athena++ headers
#include "diagnostics.hpp"
#include "../mesh/mesh.hpp"
#include "../coordinates/coordinates.hpp"
#include "../globals.hpp"

Diagnostics::Diagnostics(MeshBlock *pmb, ParameterInput *pin):
  myname("HEAD"), type(""), grid(""), long_name(""), units(""),
  prev(nullptr), next(nullptr), 
  ncycle(0), pmy_block_(pmb)
{
  std::stringstream msg;
  char cstr[80];
  std::string varnames = pin->GetOrAddString("problem", "diagnostics", "");
  std::strcpy(cstr, varnames.c_str());
  char *p = std::strtok(cstr, " ,");
  while (p != NULL) {
    std::string name(p);
    if (name == "div")  // 1.
      AddDiagnostics(Divergence(pmb));
    else if (name == "curl")  // 2.
      AddDiagnostics(Curl(pmb));
    else if (name == "mean")  // 3.
      AddDiagnostics(HydroMean(pmb));
    else if (name == "tempa") // 4.
      AddDiagnostics(TemperatureAnomaly(pmb));
    else if (name == "presa") // 5.
      AddDiagnostics(PressureAnomaly(pmb));
    else if (name == "eddyflux") // 6.
      AddDiagnostics(EddyFlux(pmb));
    else if (name == "flux") // 7.
      AddDiagnostics(TotalFlux(pmb));
    else if (name == "div_h") // 8.
      AddDiagnostics(HorizontalDivergence(pmb));
    else if (name == "b") // 9.
      AddDiagnostics(Buoyancy(pmb));
    else if (name == "radflux") // 10.
      AddDiagnostics(RadiativeFlux(pmb));
    else if (name == "am") // 11.
      AddDiagnostics(AngularMomentum(pmb));
    else if (name == "eke") // 12.
      AddDiagnostics(EddyKineticEnergy(pmb));
    else {
      msg << "### FATAL ERROR in function Diagnostics::Diagnostics"
          << std::endl << "Diagnostic variable " << name << " not defined";
			ATHENA_ERROR(msg);
    }
    p = std::strtok(NULL, " ,");
  }

  if (NGHOST < 2) {
    msg << "### FATAL ERROR in function Diagnostics::Diagnostics"
        << std::endl << "Most diagnostic variables require at least 2 ghost cells";
		ATHENA_ERROR(msg);
  }
}

Diagnostics::Diagnostics(MeshBlock *pmb, std::string name):
  myname(name), prev(nullptr), next(nullptr), ncycle(0), pmy_block_(pmb)
{
  ncells1_ = pmb->block_size.nx1 + 2*(NGHOST);
  ncells2_ = 1; 
  ncells3_ = 1;
  if (pmb->pmy_mesh->f2) ncells2_ = pmb->block_size.nx2 + 2*(NGHOST);
  if (pmb->pmy_mesh->f3) ncells3_ = pmb->block_size.nx3 + 2*(NGHOST);

  x1edge_.NewAthenaArray(ncells1_+1);
  x1edge_p1_.NewAthenaArray(ncells1_);
  x2edge_.NewAthenaArray(ncells1_+1);
  x2edge_p1_.NewAthenaArray(ncells1_);
  x3edge_.NewAthenaArray(ncells1_+1);
  x3edge_p1_.NewAthenaArray(ncells1_);

  x1area_.NewAthenaArray(ncells1_+1);
  x2area_.NewAthenaArray(ncells1_);
  x2area_p1_.NewAthenaArray(ncells1_);
  x3area_.NewAthenaArray(ncells1_);
  x3area_p1_.NewAthenaArray(ncells1_);

  vol_.NewAthenaArray(ncells1_);
  brank_.resize(Globals::nranks);
  color_.resize(Globals::nranks);
}

Diagnostics::~Diagnostics() {
  if (prev != nullptr) prev->next = next;
  if (next != nullptr) next->prev = prev;
}

Diagnostics* Diagnostics::operator[](std::string name)
{
  Diagnostics *p = this;
  while (p != nullptr) {
    if (p->myname == name)
      return p;
    p = p->next;
  }
  return p;
}

void Diagnostics::SetColor(CoordinateDirection dir) {
	MeshBlock *pmb = pmy_block_;
	NeighborBlock bblock, tblock;
	pmb->FindNeighbors(dir, bblock, tblock);

  if (dir == X1DIR) {
    if (pmb->block_size.x1min <= pmb->pmy_mesh->mesh_size.x1min) {
			bblock.snb.gid = -1;
			bblock.snb.rank = -1;
    } 
		if (pmb->block_size.x1max >= pmb->pmy_mesh->mesh_size.x1max) {
			tblock.snb.gid = -1;
			tblock.snb.rank = -1;
		}
  } else if (dir == X2DIR) {
    if (pmb->block_size.x2min <= pmb->pmy_mesh->mesh_size.x2min) {
			bblock.snb.gid = -1;
			bblock.snb.rank = -1;
		}
    if (pmb->block_size.x2max >= pmb->pmy_mesh->mesh_size.x2max) {
			tblock.snb.gid = -1;
			tblock.snb.rank = -1;
		}
  } else { // X3DIR
    if (pmb->block_size.x3min <= pmb->pmy_mesh->mesh_size.x3min) {
			bblock.snb.gid = -1;
			bblock.snb.rank = -1;
		}
    if (pmb->block_size.x3max >= pmb->pmy_mesh->mesh_size.x3max) {
			tblock.snb.gid = -1;
			tblock.snb.rank = -1;
		}
  }

#ifdef MPI_PARALLEL
  MPI_Allgather(&bblock.snb.rank, 1, MPI_INT, brank_.data(), 1, MPI_INT, MPI_COMM_WORLD);
#else
  brank_[0] = -1;
#endif

  int c = 0;
  for (int i = 0; i < Globals::nranks; ++i) {
    if (brank_[i] == -1)
      color_[i] = c++;
    else
      color_[i] = color_[brank_[i]];
  }
}
