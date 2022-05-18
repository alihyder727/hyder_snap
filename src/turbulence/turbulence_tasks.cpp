//! \file turbulence_tasks.cpp
//! \brief declared in task_list/task_list.hpp

// C/C++ headers
#include <iostream>

// Athena++ headers
#include "../task_list/task_list.hpp"
#include "../mesh/mesh.hpp"
#include "../hydro/hydro.hpp"
#include "../turbulence/turbulence_model.hpp"

TaskStatus TimeIntegratorTaskList::CalculateTurbulenceFlux(MeshBlock *pmb, int stage) {
  //std::cout << "calculate turbulence flux" << std::endl;
  if (stage <= nstages) {
    pmb->pturb->calculateFluxes(pmb->pturb->r, 2);
    return TaskStatus::next;
  }
  return TaskStatus::fail;
}


TaskStatus TimeIntegratorTaskList::SendTurbulenceFlux(MeshBlock *pmb, int stage) {
  //std::cout << "send turbulence flux" << std::endl;
  pmb->pturb->sbvar.SendFluxCorrection();
  return TaskStatus::success;
}


TaskStatus TimeIntegratorTaskList::ReceiveTurbulenceFlux(MeshBlock *pmb, int stage) {
  //std::cout << "receive turbulence flux" << std::endl;
  if (pmb->pturb->sbvar.ReceiveFluxCorrection()) {
    return TaskStatus::next;
  } else {
    return TaskStatus::fail;
  }
}


TaskStatus TimeIntegratorTaskList::IntegrateTurbulence(MeshBlock *pmb, int stage) {
  //std::cout << "integrate turbulence" << std::endl;
  TurbulenceModel *pturb = pmb->pturb;
  if (stage <= nstages) {
    // This time-integrator-specific averaging operation logic is identical to
    // IntegrateHydro, IntegrateField
    Real ave_wghts[3];
    ave_wghts[0] = 1.0;
    ave_wghts[1] = stage_wghts[stage-1].delta;
    ave_wghts[2] = 0.0;
    pmb->WeightedAve(pturb->s1, pturb->s, pturb->s2, ave_wghts);

    ave_wghts[0] = stage_wghts[stage-1].gamma_1;
    ave_wghts[1] = stage_wghts[stage-1].gamma_2;
    ave_wghts[2] = stage_wghts[stage-1].gamma_3;
    if (ave_wghts[0] == 0.0 && ave_wghts[1] == 1.0 && ave_wghts[2] == 0.0)
      pturb->s.SwapAthenaArray(pturb->s1);
    else
      pmb->WeightedAve(pturb->s, pturb->s1, pturb->s2, ave_wghts);

    const Real wght = stage_wghts[stage-1].beta*pmb->pmy_mesh->dt;
    pturb->addFluxDivergence(wght, pturb->s);

    // Hardcode an additional flux divergence weighted average for the penultimate
    // stage of SSPRK(5,4) since it cannot be expressed in a 3S* framework
    if (stage == 4 && integrator == "ssprk5_4") {
      // From Gottlieb (2009), u^(n+1) partial calculation
      ave_wghts[0] = -1.0; // -u^(n) coeff.
      ave_wghts[1] = 0.0;
      ave_wghts[2] = 0.0;
      const Real beta = 0.063692468666290; // F(u^(3)) coeff.
      const Real wght_ssp = beta*pmb->pmy_mesh->dt;
      // writing out to s2 register
      pmb->WeightedAve(pturb->s2, pturb->s1, pturb->s2, ave_wghts);
      pturb->addFluxDivergence(wght_ssp, pturb->s2);
    }

    // turbulence forcing
    Real dt = (stage_wghts[(stage-1)].beta)*(pmb->pmy_mesh->dt);
    pturb->driveTurbulence(pturb->s, pturb->r, pmb->phydro->w, dt);
    return TaskStatus::next;
  }
  return TaskStatus::fail;
}


TaskStatus TimeIntegratorTaskList::SendTurbulence(MeshBlock *pmb, int stage) {
  //std::cout << "send turbulence" << std::endl;
  if (stage <= nstages) {
    // Swap TurbulenceModel quantity in BoundaryVariable interface back to conserved var
    // formulation (also needed in SetBoundariesTurbulence() since the tasks are independent)
    pmb->pturb->sbvar.var_cc = &(pmb->pturb->s);
    pmb->pturb->sbvar.SendBoundaryBuffers();
  } else {
    return TaskStatus::fail;
  }
  return TaskStatus::success;
}


TaskStatus TimeIntegratorTaskList::ReceiveTurbulence(MeshBlock *pmb, int stage) {
  //std::cout << "recv turbulence" << std::endl;
  bool ret;
  if (stage <= nstages) {
    ret = pmb->pturb->sbvar.ReceiveBoundaryBuffers();
  } else {
    return TaskStatus::fail;
  }
  if (ret) {
    return TaskStatus::success;
  } else {
    return TaskStatus::fail;
  }
  return TaskStatus::success;
}


TaskStatus TimeIntegratorTaskList::SetBoundariesTurbulence(MeshBlock *pmb, int stage) {
  //std::cout << "set turbulence boundary" << std::endl;
  if (stage <= nstages) {
    // Set Turbulence quantity in BoundaryVariable interface to cons var formulation
    pmb->pturb->sbvar.var_cc = &(pmb->pturb->s);
    pmb->pturb->sbvar.SetBoundaries();
    return TaskStatus::success;
  }
  return TaskStatus::fail;
}
