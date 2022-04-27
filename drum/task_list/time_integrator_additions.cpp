// C/C++ headers
#include <iostream>

// Athena++ headers
#include "../hydro/implicit/implicit_solver.hpp"
#include "../hydro/polar_filter/ring_filter.hpp"
#include "../hydro/hydro.hpp"
#include "../thermodynamics/thermodynamics.hpp"
#include "../chemistry/chemistry.hpp"
#include "../radiation/radiation.hpp"
#include "../particles/material_point.hpp"
#include "../particles/particle_buffer.hpp"
#include "../particles/particles.hpp"
#include "task_list.hpp"
#include "../globals.hpp"

//! \brief apply implicit correction and physics package
enum TaskStatus TimeIntegratorTaskList::UpdateHydro(MeshBlock *pmb, int stage) {
  Hydro *ph = pmb->phydro;
  Real dt = pmb->pmy_mesh->dt;

  // polar filter
  //ph->pfilter->ApplyPolarFilter(ph->du);

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

  // polar filter
  //ph->pfilter->ApplyPolarFilter(ph->u);

  return TaskStatus::success;
}

//! \brief calculate radiation flux
TaskStatus TimeIntegratorTaskList::CalculateRadiationFlux(MeshBlock *pmb, int stage) {
  // only do radiation at last rk step
  if (stage != nstages) return TaskStatus::next;
  //if (Globals::my_rank == 0)
  //  std::cout << "I'm calculating radiation flux" << std::endl;

  Radiation *prad = pmb->prad;
  Hydro *phydro=pmb->phydro;

  if (prad->current > 0.) {
    prad->current -= pmb->pmy_mesh->dt;  // radiation is in cool-down
  } else {
    prad->current = prad->cooldown;
    for (int k = pmb->ks; k <= pmb->ke; ++k)
      for (int j = pmb->js; j <= pmb->je; ++j)
        prad->calculateRadiativeFluxes(phydro->w, pmb->pmy_mesh->time, k, j, pmb->is, pmb->ie+1);
  }
  return TaskStatus::success;
}

// particle steps
TaskStatus TimeIntegratorTaskList::IntegrateParticles(MeshBlock *pmb, int stage) {
  //if (Globals::my_rank == 0)
  //  std::cout << "I'm integrating particles" << std::endl;
  //! \todo check if it works for vl2 integrator
  Particles *ppart = pmb->ppart;
  while (ppart != nullptr) {
    // default subroutine only sest particle velocity
    // user defined particles can interact with hydro fields (two-way)
    ppart->ExchangeHydro(ppart->mp, pmb->phydro->du, pmb->phydro->w, pmb->pmy_mesh->dt);

    // copy initial state
    if (stage == 1) ppart->mp1 = ppart->mp;
    ppart->TimeIntegrate(ppart->mp, pmb->pmy_mesh->time, pmb->pmy_mesh->dt);

    if (stage > 1) {
      Real ave_wghts[3];
      ave_wghts[0] = stage_wghts[stage-1].gamma_1;
      ave_wghts[1] = stage_wghts[stage-1].gamma_2;
      ave_wghts[2] = stage_wghts[stage-1].gamma_3;
      ppart->WeightedAverage(ppart->mp, ppart->mp1, ave_wghts);
    }

    ppart = ppart->next;
  }

  return TaskStatus::success;
}

TaskStatus TimeIntegratorTaskList::SendParticles(MeshBlock *pmb, int stage) {
  // only do send/recv at last rk step
  if (stage != nstages) return TaskStatus::next;
  //if (Globals::my_rank == 0)
  //  std::cout << "I'm sending particles" << std::endl;

  Particles *ppart = pmb->ppart;
  while (ppart != nullptr) {
    ppart->ppb->DetachParticle(ppart->mp);
    ppart->ppb->SendParticle();
    ppart = ppart->next;
  }

  return TaskStatus::success;
}

TaskStatus TimeIntegratorTaskList::ReceiveParticles(MeshBlock *pmb, int stage) {
  // only do send/recv at last rk step
  if (stage != nstages) return TaskStatus::next;
  //if (Globals::my_rank == 0)
  //  std::cout << "I'm receiving particles" << std::endl;

  Particles *ppart = pmb->ppart;
  while (ppart != nullptr) {
    ppart->ppb->RecvParticle();
    ppart = ppart->next;
  }

  return TaskStatus::success;
}

TaskStatus TimeIntegratorTaskList::AttachParticles(MeshBlock *pmb, int stage) {
  // only do send/recv at the last rk step
  if (stage != nstages) return TaskStatus::next;
  //if (Globals::my_rank == 0)
  //  std::cout << "I'm attaching particles" << std::endl;

  bool ret = true;
  Particles *ppart = pmb->ppart;
  while (ppart != nullptr) {
    ret = ppart->ppb->AttachParticle(ppart->mp) && ret;
    ppart = ppart->next;
  }

  if (ret)
    return TaskStatus::success;
  else
    return TaskStatus::fail;
}

TaskStatus TimeIntegratorTaskList::ParticlesToMesh(MeshBlock *pmb, int stage) {
  // only do particle/mesh at the last rk step
  if (stage != nstages) return TaskStatus::next;
  //if (Globals::my_rank == 0)
  //  std::cout << "I'm doing particle to mesh" << std::endl;

  Particles *ppart = pmb->ppart;
  while (ppart != nullptr) {
    ppart->AggregateDensity(ppart->u, ppart->mp);
    ppart = ppart->next;
  }

  return TaskStatus::success;
}

//! \brief integrate chemistry
TaskStatus TimeIntegratorTaskList::IntegrateChemistry(MeshBlock *pmb, int stage) {
  // only do chemistry at the last rk step
  if (stage != nstages) return TaskStatus::next;
  //if (Globals::my_rank == 0)
  //  std::cout << "I'm doing chemistry" << std::endl;

  // frictional heating
  //pmb->pchem->AddFrictionalHeating(ph->u, ph->w, pmb->pmy_mesh->dt);

  pmb->pchem->TimeIntegrate(pmb->pmy_mesh->time, pmb->pmy_mesh->dt);

  return TaskStatus::success;
}

TaskStatus TimeIntegratorTaskList::MeshToParticles(MeshBlock *pmb, int stage) {
  // only do mesh/particle at the last rk step
  if (stage != nstages) return TaskStatus::next;
  //  std::cout << "I'm doing mesh to particle" << std::endl;

  Particles *ppart = pmb->ppart;
  while (ppart != nullptr) {
    ppart->Particulate(ppart->mp, ppart->u);
    ppart = ppart->next;
  }

  return TaskStatus::success;
}
