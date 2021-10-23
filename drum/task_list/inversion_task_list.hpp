#ifndef INVERSION_TASK_LIST_HPP_
#define INVERSION_TASK_LIST_HPP_

// Athena++ headers
#include "../athena.hpp"
#include "task_list.hpp"

// forward declarations
class Mesh;
class MeshBlock;
class TaskID;

class InversionTaskList : public TaskList {
public:
  InversionTaskList(ParameterInput *pin, Mesh *pm);
  ~InversionTaskList() {}
  //void AddInversionTask(uint64_t id, uint64_t dep);

  TaskStatus CalculateGradient(MeshBlock *pmb, int step);
  TaskStatus Sample(MeshBlock *pmb, int step);
  TaskStatus Optimize(MeshBlock *pmb, int step);

private:
  void AddTask(const TaskID& id, const TaskID& dep) override;
  void StartupTaskList(MeshBlock *pmb, int stage) override;
};

namespace InversionTaskNames {
  const TaskID NONE(0);
  const TaskID CALC_GRAD(1);
  const TaskID OPTIMIZE(2);
  const TaskID SAMPLE(3);
}

#endif
