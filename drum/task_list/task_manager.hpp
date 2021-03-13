#ifndef TASK_MANAGER_HPP
#define TASK_MANAGER_HPP

// C/C++ header
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
    all_tasks_(0LL), incompatible_tasks_(0LL),
    finished_tasks_(0LL), indx_first_task_(0), indx_last_task_(-1), num_tasks_left_(0)
  {}

  template<typename Y>
  void AddPackage(Y pkg, char const *name) {
    std::stringstream msg;
    if ((incompatible_tasks_ & pkg.id) == pkg.id) {
      msg << "### FATAL ERROR in function TaskManager::AddPackage"
          << std::endl << "Package'" << name << "' "
          << "is incompatible with existing task." << std::endl;
      ATHENA_ERROR(msg);
    } else {
      all_tasks_ |= pkg.id;
      incompatible_tasks_ |= pkg.incompatible;
    }
  }

  void RemoveTask(uint64_t id) {
    all_tasks_ &= ~id;
  }

  void Reset() {
    // hamming weight
    int ntasks = 0;
    while (all_tasks_ > 0) {     // until all bits are zero
      if ((all_tasks_ & 1) == 1) // check lower bit
        ntasks++;
      all_tasks_ >>= 1;          // shift bits, removing lower bit
    }

    indx_first_task_ = 0;
    indx_last_task_  = ntasks - 1;
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
    return all_tasks_ & id;
  }

  template<typename Y>
  TaskListStatus DoNextTask(AthenaArray<Real> &u, AthenaArray<Real> const& w,
    Real time, Real dt, std::vector<Y> const &tlist);

private:
// data
  T* pclass_;
  uint64_t all_tasks_;
  uint64_t incompatible_tasks_;
  uint64_t finished_tasks_;
  int indx_first_task_, indx_last_task_;
  int num_tasks_left_;
};

template<typename T>
template<typename Y>
TaskListStatus TaskManager<T>::DoNextTask(AthenaArray<Real> &u, AthenaArray<Real> const& w,
  Real time, Real dt, std::vector<Y> const &tlist)
{
  bool first_task = true;
  TaskStatus ret;
  if (num_tasks_left_ == 0) return TaskListStatus::complete;

  for (int i=indx_first_task_; i <= indx_last_task_; i++) {
    Y const& todo = tlist[i];
    if (Unfinished(todo.id)) { // task not done
      if (DependencyClear(todo.dependency)) {  // dependency clear
        ret = (pclass_->*todo.Function)(u, w, time, dt);
        if (ret != TaskStatus::fail) { // success
          num_tasks_left_--;
          finished_tasks_ &= todo.id;
          if (first_task) indx_first_task_++;
          if (num_tasks_left_ == 0) return TaskListStatus::complete;
          if (ret == TaskStatus::next) continue;
          return TaskListStatus::running;
        }
      }
    } else if (first_task) { // this task is already done AND it is at the top of the list
      indx_first_task_++;
    }
    first_task = false; // no longer being first task
  }
  // there are still tasks to do but nothing can be done now
  return TaskListStatus::stuck;
}

#endif
