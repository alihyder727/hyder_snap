#include "../mesh/mesh.hpp"
#include "../hydro/hydro.hpp"
#include "../thermodynamics/thermodynamics.hpp"
#include "physics.hpp"

TaskStatus Physics::TopCooling(AthenaArray<Real> &u,
  AthenaArray<Real> const& w, Real time, Real dt)
{
  MeshBlock *pmb = pmy_block;
  Thermodynamics *pthermo = pmb->pthermo;
  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;

  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j) {
      Real cv = pthermo->Cv(w.at(k,j,ie));
      u(IEN,k,j,ie) += dt*Jcool_*w(IDN,k,j,ie)*cv;
    }

  return TaskStatus::success;
}

TaskStatus Physics::BotHeating(AthenaArray<Real> &u,
  AthenaArray<Real> const& w, Real time, Real dt)
{
  MeshBlock *pmb = pmy_block;
  Thermodynamics *pthermo = pmb->pthermo;
  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;

  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j) {
      Real cv = pthermo->Cv(w.at(k,j,is));
      u(IEN,k,j,is) += dt*Jheat_*w(IDN,k,j,is)*cv;
    }

  return TaskStatus::success;
}
