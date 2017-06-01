#ifndef MG_TASK_LIST_HPP
#define MG_TASK_LIST_HPP
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//!   \file mg_task_list.hpp
//    \brief provides functionality to control dynamic execution using tasks

// Athena++ headers
#include "../athena.hpp"
#include "../mesh/mesh.hpp"

// forward declarations
class Mesh;
class MeshBlock;
class TaskList;

// return codes for functions working on individual Tasks and TaskList
enum TaskStatus {TASK_FAIL, TASK_SUCCESS, TASK_NEXT};
enum TaskListStatus {TL_RUNNING, TL_STUCK, TL_COMPLETE, TL_NOTHING_TO_DO};

//----------------------------------------------------------------------------------------
//! \struct IntegratorWeight
//  \brief weights used in time integrator tasks

struct IntegratorWeight {
  Real a,b,c;
};

//----------------------------------------------------------------------------------------
//! \struct MGTask
//  \brief data and function pointer for an individual MGTask

struct MGTask {
  uint64_t task_id;      // encodes step & task using bit positions in HydroTasks
  uint64_t dependency;   // encodes dependencies to other tasks using " " " "
  enum TaskStatus (MultigridTaskList::*TaskFunc)(Multigrid*, int);  // ptr to member function
};


//----------------------------------------------------------------------------------------
//! \class MultigridTaskList
//  \brief data and function definitions for MultigridTaskList class

class MultigridTaskList
{
public:
  MultigridTaskList() {};
  ~MultigridTaskList() {};

  // data
  int ntasks;     // number of tasks in this list

  // functions
  enum TaskListStatus DoAllAvailableTasks(Multigrid *pmg, int step, TaskState &ts);
  void DoTaskListOneSubStep(Mesh *pmesh, int step);
  void ClearTaskList(void);
  void AddTask(uint64_t id, uint64_t dep);

  // functions
  enum TaskStatus MGGravityStartReceive(MeshBlock *pmb, int step);
  enum TaskStatus MGGravityClearBoundary(MeshBlock *pmb, int step);
  enum TaskStatus MGGravitySendBoundary(MeshBlock *pmb, int step);

private:
  MultigridDriver* pmy_mgdriver_;
  struct MGTask task_list_[64];
}

//----------------------------------------------------------------------------------------
// 64-bit integers with "1" in different bit positions used to ID each Multigrid task.

namespace MultigridTaskNames {
  const uint64_t NONE=0;
  const uint64_t MG_STARTRECV1   = 1LL<<0;
  const uint64_t MG_CLEARBND1    = 1LL<<1;
  const uint64_t MG_SENDBND1     = 1LL<<2;
  const uint64_t MG_RECVBND1     = 1LL<<3;
  const uint64_t MG_SMOOTHRED1   = 1LL<<4;
  const uint64_t MG_SMOOTHBLACK1 = 1LL<<5;
  const uint64_t MG_STARTRECV2   = 1LL<<6;
  const uint64_t MG_CLEARBND2    = 1LL<<7;
  const uint64_t MG_SENDBND2     = 1LL<<8;
  const uint64_t MG_RECVBND2     = 1LL<<9;
  const uint64_t MG_SMOOTHRED2   = 1LL<<10;
  const uint64_t MG_SMOOTHBLACK2 = 1LL<<11;
  const uint64_t MG_PHYSBND1     = 1LL<<12;
  const uint64_t MG_PHYSBND2     = 1LL<<13;
  const uint64_t MG_RESTRICT     = 1LL<<14;
  const uint64_t MG_PROLONG      = 1LL<<15;
  const uint64_t MG_FMGPROLONG   = 1LL<<16;
};

#endif // MG_TASK_LIST_HPP
