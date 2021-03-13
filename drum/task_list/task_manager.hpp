#ifndef TASK_MANAGER_HPP
#define TASK_MANAGER_HPP

// C/C++ header
#include <bitset>
#include <vector>
#include <sstream>

// Athena++ header
#include "task_list.hpp"
#include "../athena.hpp"

//! \brief manage all tasks to do
template<typename T>
class TaskManager {
public:
// functions
  TaskManager(T *pclass): pclass_(pclass),
    all_tasks_(0LL), conflict_tasks_(0LL), finished_tasks_(0LL),
    indx_first_task_(0), num_tasks_left_(0)
  {}

  template<typename Y>
  void AddPackage(Y pkg, char const *name) {
    std::stringstream msg;
    if ((conflict_tasks_ & pkg.id) == pkg.id) {
      msg << "### FATAL ERROR in function TaskManager::AddPackage"
          << std::endl << "Package'" << name << "' "
          << "is incompatible with existing task." << std::endl;
      ATHENA_ERROR(msg);
    } else {
      all_tasks_ |= pkg.id;
      conflict_tasks_ |= pkg.conflict;
    }
  }

  void RemoveTask(uint64_t id) {
    all_tasks_ &= ~id;
  }

  void Reset() {
    // hamming weight
    int ntasks = 0;
    uint64_t tmp = all_tasks_;
    while (tmp > 0) {     // until all bits are zero
      if ((tmp & 1) == 1) // check lower bit
        ntasks++;
      tmp >>= 1;          // shift bits, removing lower bit
    }

    indx_first_task_ = 0;
    num_tasks_left_  = ntasks;
    finished_tasks_  = 0LL;
  }

  bool Unfinished(uint64_t id) const {
    return (finished_tasks_ & id) == 0LL;
  }

  bool DependencyClear(uint64_t dep) const {
    return (finished_tasks_ & dep) == dep;
  }

  bool HasTask(uint64_t id) const {
    return (all_tasks_ & id) == id;
  }

  template<typename Y>
  TaskListStatus DoNextJob(AthenaArray<Real> &u, AthenaArray<Real> const& w,
    Real time, Real dt, std::vector<Y> const &tlist);

private:
// data
  T* pclass_;
  uint64_t all_tasks_;
  uint64_t conflict_tasks_;
  uint64_t finished_tasks_;
  int indx_first_task_;
  int num_tasks_left_;
};

template<typename T>
template<typename Y>
TaskListStatus TaskManager<T>::DoNextJob(AthenaArray<Real> &u, AthenaArray<Real> const& w,
  Real time, Real dt, std::vector<Y> const &tlist)
{
  TaskStatus ret;
  if ((num_tasks_left_ == 0) || (tlist.size() == 0))
    return TaskListStatus::complete;

  std::stringstream msg;
  if (num_tasks_left_ > tlist.size()) {
    msg << "### FATAL ERROR in function TaskManager::DoNextJob"
        << std::endl << num_tasks_left_ << " tasks to do." << std::endl
        << "But only " << tlist.size() << " jobs available." << std::endl;
    ATHENA_ERROR(msg);
  }

  // find next task
  while (!(HasTask(tlist[indx_first_task_].id) && Unfinished(tlist[indx_first_task_].id)))
    indx_first_task_++;

  for (int i = indx_first_task_; i < tlist.size(); i++) {
    Y const& todo = tlist[i];
    // has this task, not done, and denpendency clear
    if (HasTask(todo.id) && Unfinished(todo.id) && DependencyClear(todo.dep)) {
      ret = (pclass_->*todo.Function)(u, w, time, dt);
      if (ret != TaskStatus::fail) { // success
        num_tasks_left_--;
        finished_tasks_ |= todo.id;
        if (num_tasks_left_ == 0) return TaskListStatus::complete;
        if (ret == TaskStatus::next) continue;
        return TaskListStatus::running;
      }
    }
  }
  // there are still tasks to do but nothing can be done now
  return TaskListStatus::stuck;
}

#endif
