// C/C++ headers
#include <stdexcept>
#include <sstream>
#include <iostream>

// Athena++ headers
#include "../task_list/task_manager.hpp"
#include "physics.hpp"

Physics::~Physics() {
  delete ptm;
}

void Physics::ApplyPhysicsPackages(AthenaArray<Real> &u,
  AthenaArray<Real> const& w, Real time, Real dt)
{
  std::stringstream msg;
  int count = 0;
  ptm->Reset();

  while (count < 100) {
    TaskListStatus status = ptm->DoNextJob(u, w, time, dt, packages_);
    if (status == TaskListStatus::complete) break;
    count++;
  }

  if (count >= 100) {
    msg << "### FATAL ERROR in Physics::ApplyPhysicsPackages"
        << std::endl << "Physics Package stuck." << std::endl;
    ATHENA_ERROR(msg);
  }
}

void Physics::DumpRestartData(std::string fname)
{}

void Physics::LoadRestartData(std::string fname)
{}

