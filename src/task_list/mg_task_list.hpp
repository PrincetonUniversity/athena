#ifndef TASK_LIST_MG_TASK_LIST_HPP_
#define TASK_LIST_MG_TASK_LIST_HPP_
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
#include "../multigrid/multigrid.hpp"
#include "./task_list.hpp"

// forward declarations
class Mesh;
class MeshBlock;
class Multigrid;
class MultigridDriver;
class MultigridTaskList;

//----------------------------------------------------------------------------------------
//! \struct MGTask
//  \brief data and function pointer for an individual MGTask

struct MGTask {
  uint64_t task_id;      // encodes task using bit positions in MultigridTaskNames
  uint64_t dependency;   // encodes dependencies to other tasks using " " " "
  enum TaskStatus (MultigridTaskList::*TaskFunc)(Multigrid*);  // ptr to a task
};


//----------------------------------------------------------------------------------------
//! \class MultigridTaskList
//  \brief data and function definitions for MultigridTaskList class

class MultigridTaskList {
public:
  explicit MultigridTaskList(MultigridDriver *pmd) : pmy_mgdriver_(pmd) {}
  ~MultigridTaskList() {}

  // data
  int ntasks;     // number of tasks in this list

  // functions
  enum TaskListStatus DoAllAvailableTasks(Multigrid *pmg, TaskState &ts);
  void DoTaskListOneStage(MultigridDriver *pmd);
  void ClearTaskList(void) {ntasks=0;}
  void AddMultigridTask(uint64_t id, uint64_t dep);

  // functions
  enum TaskStatus StartReceive(Multigrid *pmg);
  enum TaskStatus StartReceiveFace(Multigrid *pmg);
  enum TaskStatus ClearBoundary(Multigrid *pmg);
  enum TaskStatus ClearBoundaryFace(Multigrid *pmg);
  enum TaskStatus SendBoundary(Multigrid *pmg);
  enum TaskStatus SendBoundaryFace(Multigrid *pmg);
  enum TaskStatus ReceiveBoundary(Multigrid *pmg);
  enum TaskStatus ReceiveBoundaryFace(Multigrid *pmg);
  enum TaskStatus SmoothRed(Multigrid *pmg);
  enum TaskStatus SmoothBlack(Multigrid *pmg);
  enum TaskStatus PhysicalBoundary(Multigrid *pmg);
  enum TaskStatus Restrict(Multigrid *pmg);
  enum TaskStatus Prolongate(Multigrid *pmg);
  enum TaskStatus FMGProlongate(Multigrid *pmg);

  void SetMGTaskListToFiner(int nsmooth, int ngh, int flag = 0);
  void SetMGTaskListToCoarser(int nsmooth, int ngh);
  void SetMGTaskListFMGProlongate(int flag = 0);

private:
  MultigridDriver* pmy_mgdriver_;
  struct MGTask task_list_[64];
};

//----------------------------------------------------------------------------------------
// 64-bit integers with "1" in different bit positions used to ID each Multigrid task.

namespace MultigridTaskNames {
  const uint64_t NONE           = 0;
  const uint64_t MG_STARTRECV0  = 1LL<<0;
  const uint64_t MG_STARTRECV0F = 1LL<<1;
  const uint64_t MG_STARTRECV1R = 1LL<<2;
  const uint64_t MG_STARTRECV1B = 1LL<<3;
  const uint64_t MG_STARTRECV2R = 1LL<<4;
  const uint64_t MG_STARTRECV2B = 1LL<<5;
  const uint64_t MG_STARTRECVL  = 1LL<<6;
  const uint64_t MG_CLEARBND0   = 1LL<<7;
  const uint64_t MG_CLEARBND0F  = 1LL<<8;
  const uint64_t MG_CLEARBND1R  = 1LL<<9;
  const uint64_t MG_CLEARBND1B  = 1LL<<10;
  const uint64_t MG_CLEARBND2R  = 1LL<<11;
  const uint64_t MG_CLEARBND2B  = 1LL<<12;
  const uint64_t MG_CLEARBNDL   = 1LL<<13;
  const uint64_t MG_SENDBND0    = 1LL<<14;
  const uint64_t MG_SENDBND0F   = 1LL<<15;
  const uint64_t MG_SENDBND1R   = 1LL<<16;
  const uint64_t MG_SENDBND1B   = 1LL<<17;
  const uint64_t MG_SENDBND2R   = 1LL<<18;
  const uint64_t MG_SENDBND2B   = 1LL<<19;
  const uint64_t MG_SENDBNDL    = 1LL<<20;
  const uint64_t MG_RECVBND0    = 1LL<<21;
  const uint64_t MG_RECVBND0F   = 1LL<<22;
  const uint64_t MG_RECVBND1R   = 1LL<<23;
  const uint64_t MG_RECVBND1B   = 1LL<<24;
  const uint64_t MG_RECVBND2R   = 1LL<<25;
  const uint64_t MG_RECVBND2B   = 1LL<<26;
  const uint64_t MG_RECVBNDL    = 1LL<<27;
  const uint64_t MG_SMOOTH1R    = 1LL<<28;
  const uint64_t MG_SMOOTH1B    = 1LL<<29;
  const uint64_t MG_SMOOTH2R    = 1LL<<30;
  const uint64_t MG_SMOOTH2B    = 1LL<<31;
  const uint64_t MG_PHYSBND0    = 1LL<<32;
  const uint64_t MG_PHYSBND1R   = 1LL<<33;
  const uint64_t MG_PHYSBND1B   = 1LL<<34;
  const uint64_t MG_PHYSBND2R   = 1LL<<35;
  const uint64_t MG_PHYSBND2B   = 1LL<<36;
  const uint64_t MG_PHYSBNDL    = 1LL<<37;
  const uint64_t MG_RESTRICT    = 1LL<<38;
  const uint64_t MG_PROLONG     = 1LL<<39;
  const uint64_t MG_FMGPROLONG  = 1LL<<40;
}; // namespace MultigridTaskNames

#endif // TASK_LIST_MG_TASK_LIST_HPP_
