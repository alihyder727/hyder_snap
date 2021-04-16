// C/C++ headers
#include <iostream>

// Athena++ headers
#include "../hydro/implicit/implicit_solver.hpp"
#include "../hydro/hydro.hpp"
#include "../thermodynamics/thermodynamics.hpp"
#include "../chemistry/chemistry.hpp"
#include "../radiation/radiation.hpp"
#include "../physics/physics.hpp"
#include "task_list.hpp"

//! \brief apply implicit correction and physics package
enum TaskStatus TimeIntegratorTaskList::UpdateHydro(MeshBlock *pmb, int stage) {
  Hydro *ph = pmb->phydro;
  Real dt = pmb->pmy_mesh->dt;

  // do implicit coorection at every stage
  if (ph->implicit_flag == 1)
    ph->pimp1->PartialCorrection(ph->du, ph->w, stage_wghts[stage-1].beta*dt);
  else if (ph->implicit_flag == 2)
    ph->pimp2->PartialCorrection(ph->du, ph->w, stage_wghts[stage-1].beta*dt);
  else if (ph->implicit_flag == 3)
    ph->pimp3->PartialCorrection(ph->du, ph->w, stage_wghts[stage-1].beta*dt);
  else { // full implicit
    ph->pimp1->PartialCorrection(ph->du, ph->w, stage_wghts[stage-1].beta*dt);
    ph->pimp2->PartialCorrection(ph->du, ph->w, stage_wghts[stage-1].beta*dt);
    ph->pimp3->PartialCorrection(ph->du, ph->w, stage_wghts[stage-1].beta*dt);
  }
  ph->ImplicitUpdate(ph->du);
  Real wghts[3];
  wghts[0] = 1.;
  wghts[1] = 1.;
  wghts[2] = 0.;
  pmb->WeightedAve(ph->u, ph->du, ph->u2, wghts);

  // do physics package at the last stage
  if (stage == nstages)
    pmb->pphy->ApplyPhysicsPackages(ph->u, ph->w, pmb->pmy_mesh->time, pmb->pmy_mesh->dt);

  return TaskStatus::success;
}

//! \brief integrate chemistry
TaskStatus TimeIntegratorTaskList::IntegrateChemistry(MeshBlock *pmb, int stage) {
  Hydro *ph = pmb->phydro;

  if (stage != nstages) return TaskStatus::next;

  // frictional heating
  pmb->pchem->AddFrictionalHeating(ph->u, ph->w, pmb->pmy_mesh->dt);

  // do slow chemistry
  pmb->pchem->EvolveOneStep(ph->u, pmb->pmy_mesh->time, pmb->pmy_mesh->dt);

  // do fast chemistry
  pmb->pthermo->SaturationAdjustment(ph->u);

  return TaskStatus::success;
}

//! \brief calculate radiation flux
TaskStatus TimeIntegratorTaskList::CalculateRadiationFlux(MeshBlock *pmb, int stage) {
  // only do radiation at last rk step
  if (stage != nstages) return TaskStatus::next;

  Radiation *prad = pmb->prad;
  Hydro *phydro=pmb->phydro;

  if (prad->current > 0.) {
    prad->current -= pmb->pmy_mesh->dt;  // radiation is in cool-down
  } else {
    prad->current = prad->cooldown;
    for (int k = pmb->ks; k <= pmb->ke; ++k)
      for (int j = pmb->js; j <= pmb->je; ++j)
        prad->CalculateFluxes(phydro->w, pmb->pmy_mesh->time, k, j, pmb->is, pmb->ie+1);
  }
  return TaskStatus::success;
}
