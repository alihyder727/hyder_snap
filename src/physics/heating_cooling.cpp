#include "../mesh/mesh.hpp"
#include "../hydro/hydro.hpp"
#include "../thermodynamics/thermodynamics.hpp"
#include "../coordinates/coordinates.hpp"
#include "physics.hpp"

TaskStatus Physics::TopCooling(AthenaArray<Real> &du,
  AthenaArray<Real> const& w, Real time, Real dt)
{
  MeshBlock *pmb = pmy_block;
  Thermodynamics *pthermo = pmb->pthermo;
  int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;

  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j) {
      Real cv = pthermo->GetMeanCv(w.at(k,j,ie));
      du(IEN,k,j,ie) -= dt*dTdt_*w(IDN,k,j,ie)*cv;
    }

  return TaskStatus::success;
}

TaskStatus Physics::BotHeating(AthenaArray<Real> &du,
  AthenaArray<Real> const& w, Real time, Real dt)
{
  MeshBlock *pmb = pmy_block;
	Coordinates *pcoord = pmb->pcoord;
  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int je = pmb->je; int ke = pmb->ke;

  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j) {
      du(IEN,k,j,is) += dt*hflux_*pcoord->GetFace1Area(k,j,is)
				/pcoord->GetCellVolume(k,j,is);
    }

  return TaskStatus::success;
}
