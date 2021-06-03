// C/C++ headers
#include <iostream>

// Athena++ headers
#include "../hydro/implicit/implicit_solver.hpp"
#include "../hydro/hydro.hpp"
#include "../thermodynamics/thermodynamics.hpp"
#include "../chemistry/chemistry.hpp"
#include "../radiation/radiation.hpp"
#include "../particles/material_point.hpp"
#include "../particles/particle_buffer.hpp"
#include "../particles/particles.hpp"
#include "task_list.hpp"

//! \brief apply implicit correction and physics package
enum TaskStatus TimeIntegratorTaskList::UpdateHydro(MeshBlock *pmb, int stage) {
  Hydro *ph = pmb->phydro;
  Real dt = pmb->pmy_mesh->dt;

  // do implicit coorection at every stage
  // X3DIR
  if ((ph->implicit_flag & (1<<2)) && (pmb->ncells3 > 1)) {
    ph->pimp->SetDirection(X3DIR);
    if (ph->implicit_flag & (1<<3))
      ph->pimp->FullCorrection(ph->du, ph->w, stage_wghts[stage-1].beta*dt);
    else
      ph->pimp->PartialCorrection(ph->du, ph->w, stage_wghts[stage-1].beta*dt);
  }

  // X2DIR
  if ((ph->implicit_flag & (1<<1)) && (pmb->ncells2 > 1)) {
    ph->pimp->SetDirection(X2DIR);
    if (ph->implicit_flag & (1<<3))
      ph->pimp->FullCorrection(ph->du, ph->w, stage_wghts[stage-1].beta*dt);
    else
      ph->pimp->PartialCorrection(ph->du, ph->w, stage_wghts[stage-1].beta*dt);
  }

  // X1DIR
  if (ph->implicit_flag & 1) {
    ph->pimp->SetDirection(X1DIR);
    if (ph->implicit_flag & (1<<3))
      ph->pimp->FullCorrection(ph->du, ph->w, stage_wghts[stage-1].beta*dt);
    else
      ph->pimp->PartialCorrection(ph->du, ph->w, stage_wghts[stage-1].beta*dt);
  }

  Real wghts[3];
  wghts[0] = 1.;
  wghts[1] = 1.;
  wghts[2] = 0.;
  pmb->WeightedAve(ph->u, ph->du, ph->u2, wghts);

  return TaskStatus::success;
}

//! \brief integrate chemistry
TaskStatus TimeIntegratorTaskList::IntegrateChemistry(MeshBlock *pmb, int stage) {
  Hydro *ph = pmb->phydro;

  if (stage != nstages) return TaskStatus::next;

  // frictional heating
  //pmb->pchem->AddFrictionalHeating(ph->u, ph->w, pmb->pmy_mesh->dt);

  // do slow chemistry
  //pmb->pchem->EvolveOneStep(ph->u, pmb->pmy_mesh->time, pmb->pmy_mesh->dt);

  // do fast chemistry
  //pmb->pthermo->SaturationAdjustment(ph->u);

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

// particle steps
TaskStatus TimeIntegratorTaskList::IntegrateParticles(MeshBlock *pmb, int stage) {
  //! \todo check if it works for vl2 integrator
  Particles *ppar = pmb->ppar;
  while (ppar != nullptr) {
    ppar->ExchangeHydro(pmb->phydro->du, pmb->phydro->w);

    // copy initial state
    if (stage == 1) ppar->mp1 = ppar->mp;
    ppar->TimeIntegrate(ppar->mp, pmb->pmy_mesh->time, pmb->pmy_mesh->dt);

    if (stage > 1) {
      Real ave_wghts[3];
      ave_wghts[0] = stage_wghts[stage-1].gamma_1;
      ave_wghts[1] = stage_wghts[stage-1].gamma_2;
      ave_wghts[2] = stage_wghts[stage-1].gamma_3;
      ppar->WeightedAverage(ppar->mp, ppar->mp1, ave_wghts);
    }

    ppar = ppar->next;
  }

  return TaskStatus::success;
}

TaskStatus TimeIntegratorTaskList::MeshToParticles(MeshBlock *pmb, int stage) {
  if (stage != nstages) return TaskStatus::next;

  Particles *ppar = pmb->ppar;
  while (ppar != nullptr) {
    ppar->Particulate(ppar->dc);
    ppar = ppar->next;
  }

  return TaskStatus::success;
}

TaskStatus TimeIntegratorTaskList::SendParticles(MeshBlock *pmb, int stage) {
  // only do send/recv at last rk step
  if (stage != nstages) return TaskStatus::next;

  Particles *ppar = pmb->ppar;
  while (ppar != nullptr) {
    ppar->ppb->DetachParticle(ppar->mp);
    ppar->ppb->SendParticle();
    ppar = ppar->next;
  }

  return TaskStatus::success;
}

TaskStatus TimeIntegratorTaskList::ReceiveParticles(MeshBlock *pmb, int stage) {
  // only do send/recv at last rk step
  if (stage != nstages) return TaskStatus::next;

  Particles *ppar = pmb->ppar;
  while (ppar != nullptr) {
    ppar->ppb->RecvParticle();
    ppar = ppar->next;
  }

  return TaskStatus::success;
}

TaskStatus TimeIntegratorTaskList::AttachParticles(MeshBlock *pmb, int stage) {
  // only do send/recv at last rk step
  if (stage != nstages) return TaskStatus::next;

  bool ret = true;
  Particles *ppar = pmb->ppar;
  while (ppar != nullptr) {
    ret = ppar->ppb->AttachParticle(ppar->mp) && ret;
    ppar = ppar->next;
  }

  if (ret)
    return TaskStatus::success;
  else
    return TaskStatus::fail;
}

TaskStatus TimeIntegratorTaskList::ParticlesToMesh(MeshBlock *pmb, int stage) {
  // only do send/recv at last rk step
  if (stage != nstages) return TaskStatus::next;

  Particles *ppar = pmb->ppar;
  while (ppar != nullptr) {
    ppar->AggregateMass(ppar->c);
    ppar = ppar->next;
  }

  return TaskStatus::success;
}
