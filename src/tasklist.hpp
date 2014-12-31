#ifndef TASKLIST_HPP
#define TASKLIST_HPP
//======================================================================================
//! \file tasklist.hpp
//  \brief definition of the task list
//======================================================================================


struct Task {
public:
  enum task taskid;
  unsigned long int depend;
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


inline void TaskList::AddTask(enum task t, unsigned long int dependence)
{
  task[ntask].taskid=t;
  task[ntask].depend=dependence;
  ntask++;
  return;
}

#endif
