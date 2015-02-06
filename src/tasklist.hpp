#ifndef TASKLIST_HPP
#define TASKLIST_HPP
//======================================================================================
//! \file tasklist.hpp
//  \brief definition of the task list
//======================================================================================

#include "athena.hpp"
#include "mesh.hpp"

class MeshBlock;

typedef bool (*TaskFunc_t)(MeshBlock*, int);

struct Task {
  enum task taskid;
  unsigned long int depend;
  TaskFunc_t TaskFunc;
  int task_arg;
};

class TaskList {
private:
  int ntask;
  Task task[64];
  friend class MeshBlock;
public:
  TaskList() : ntask(0) {};
  ~TaskList() {};
  void AddTask(enum task, unsigned long int dependence);
};

#endif
