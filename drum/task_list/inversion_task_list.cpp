// C/C++ headers
#include <iostream>
#include <sstream>
#include <stdexcept>

// Athena++ headers
#include "inversion_task_list.hpp"
#include "../athena.hpp"
#include "../parameter_input.hpp"
#include "../mesh/mesh.hpp"
#include "../hydro/hydro.hpp"
#include "../bvals/bvals.hpp"
#include "../inversion/inversion.hpp"

using namespace InversionTaskNames;

InversionTaskList::InversionTaskList(ParameterInput *pin, Mesh *pm)
{
  std::cout << std::endl;
  std::cout << "######################################################" << std::endl;
  std::cout << "##                 WELCOME TO HARP                  ##" << std::endl;
  std::cout << "## A High-performance Atmospheric Radiative Package ##" << std::endl;
  std::cout << "######################################################" << std::endl;
  nstages = 1;
  std::string method = pin->GetString("inversion", "method");

  // Now assemble list of tasks for each step of inversion task

  if (method == "LM") {
    AddTask(CALC_GRAD,NONE);
    AddTask(OPTIMIZE,CALC_GRAD);
  } else if (method == "mcmc") {
    AddTask(SAMPLE,NONE);
  }
}

void InversionTaskList::AddTask(const TaskID& id, const TaskID& dep) {
  task_list_[ntasks].task_id=id;
  task_list_[ntasks].dependency=dep;

  if (id == CALC_GRAD) {
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&InversionTaskList::CalculateGradient);
  } else if (id == OPTIMIZE) {
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&InversionTaskList::Optimize);
  } else if (id == SAMPLE) {
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&InversionTaskList::Sample);
  } else {
      std::stringstream msg;
      msg << "### FATAL ERROR in AddTask" << std::endl
          << "Invalid Task is specified" << std::endl;
      ATHENA_ERROR(msg);
  }
  ntasks++;
  return;
}

void InversionTaskList::StartupTaskList(MeshBlock *pmb, int stage) {}

TaskStatus InversionTaskList::CalculateGradient(MeshBlock *pmb, int step) {
  //std::cout << "Calculate gradient" << std::endl;
  return TaskStatus::success;
}

TaskStatus InversionTaskList::Optimize(MeshBlock *pmb, int step) {
  //std::cout << "Optimize" << std::endl;
  return TaskStatus::success;
}

TaskStatus InversionTaskList::Sample(MeshBlock *pmb, int step) {
  Inversion *pinvt = pmb->pinvt;
  pinvt->MCMCStep(pinvt->obj);
  //std::cout << "Sample" << std::endl;
  return TaskStatus::success;
}
