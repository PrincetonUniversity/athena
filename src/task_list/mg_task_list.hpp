#ifndef TASK_LIST_MG_TASK_LIST_HPP_
#define TASK_LIST_MG_TASK_LIST_HPP_
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//!   \file mg_task_list.hpp
//    \brief provides functionality to control dynamic execution using tasks

// C headers

// C++ headers

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
  std::uint64_t task_id;      // encodes task using bit positions in MultigridTaskNames
  std::uint64_t dependency;   // encodes dependencies to other tasks using " " " "
  TaskStatus (MultigridTaskList::*TaskFunc)(Multigrid*);  // ptr to a task
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
  TaskListStatus DoAllAvailableTasks(Multigrid *pmg, TaskStates &ts);
  void DoTaskListOneStage(MultigridDriver *pmd);
  void ClearTaskList() {ntasks=0;}
  void AddMultigridTask(std::uint64_t id, std::uint64_t dep);

  // functions
  TaskStatus StartReceive(Multigrid *pmg);
  TaskStatus StartReceiveFace(Multigrid *pmg);
  TaskStatus ClearBoundary(Multigrid *pmg);
  TaskStatus ClearBoundaryFace(Multigrid *pmg);
  TaskStatus SendBoundary(Multigrid *pmg);
  TaskStatus SendBoundaryFace(Multigrid *pmg);
  TaskStatus ReceiveBoundary(Multigrid *pmg);
  TaskStatus ReceiveBoundaryFace(Multigrid *pmg);
  TaskStatus SmoothRed(Multigrid *pmg);
  TaskStatus SmoothBlack(Multigrid *pmg);
  TaskStatus PhysicalBoundary(Multigrid *pmg);
  TaskStatus Restrict(Multigrid *pmg);
  TaskStatus Prolongate(Multigrid *pmg);
  TaskStatus FMGProlongate(Multigrid *pmg);

  void SetMGTaskListToFiner(int nsmooth, int ngh, int flag = 0);
  void SetMGTaskListToCoarser(int nsmooth, int ngh);
  void SetMGTaskListFMGProlongate(int flag = 0);

 private:
  MultigridDriver* pmy_mgdriver_;
  MGTask task_list_[64];
};

//----------------------------------------------------------------------------------------
// 64-bit integers with "1" in different bit positions used to ID each Multigrid task.

namespace MultigridTaskNames {
const std::uint64_t NONE           = 0ULL;
const std::uint64_t MG_STARTRECV0  = 1ULL<<0;
const std::uint64_t MG_STARTRECV0F = 1ULL<<1;
const std::uint64_t MG_STARTRECV1R = 1ULL<<2;
const std::uint64_t MG_STARTRECV1B = 1ULL<<3;
const std::uint64_t MG_STARTRECV2R = 1ULL<<4;
const std::uint64_t MG_STARTRECV2B = 1ULL<<5;
const std::uint64_t MG_STARTRECVL  = 1ULL<<6;
const std::uint64_t MG_CLEARBND0   = 1ULL<<7;
const std::uint64_t MG_CLEARBND0F  = 1ULL<<8;
const std::uint64_t MG_CLEARBND1R  = 1ULL<<9;
const std::uint64_t MG_CLEARBND1B  = 1ULL<<10;
const std::uint64_t MG_CLEARBND2R  = 1ULL<<11;
const std::uint64_t MG_CLEARBND2B  = 1ULL<<12;
const std::uint64_t MG_CLEARBNDL   = 1ULL<<13;
const std::uint64_t MG_SENDBND0    = 1ULL<<14;
const std::uint64_t MG_SENDBND0F   = 1ULL<<15;
const std::uint64_t MG_SENDBND1R   = 1ULL<<16;
const std::uint64_t MG_SENDBND1B   = 1ULL<<17;
const std::uint64_t MG_SENDBND2R   = 1ULL<<18;
const std::uint64_t MG_SENDBND2B   = 1ULL<<19;
const std::uint64_t MG_SENDBNDL    = 1ULL<<20;
const std::uint64_t MG_RECVBND0    = 1ULL<<21;
const std::uint64_t MG_RECVBND0F   = 1ULL<<22;
const std::uint64_t MG_RECVBND1R   = 1ULL<<23;
const std::uint64_t MG_RECVBND1B   = 1ULL<<24;
const std::uint64_t MG_RECVBND2R   = 1ULL<<25;
const std::uint64_t MG_RECVBND2B   = 1ULL<<26;
const std::uint64_t MG_RECVBNDL    = 1ULL<<27;
const std::uint64_t MG_SMOOTH1R    = 1ULL<<28;
const std::uint64_t MG_SMOOTH1B    = 1ULL<<29;
const std::uint64_t MG_SMOOTH2R    = 1ULL<<30;
const std::uint64_t MG_SMOOTH2B    = 1ULL<<31;
const std::uint64_t MG_PHYSBND0    = 1ULL<<32;
const std::uint64_t MG_PHYSBND1R   = 1ULL<<33;
const std::uint64_t MG_PHYSBND1B   = 1ULL<<34;
const std::uint64_t MG_PHYSBND2R   = 1ULL<<35;
const std::uint64_t MG_PHYSBND2B   = 1ULL<<36;
const std::uint64_t MG_PHYSBNDL    = 1ULL<<37;
const std::uint64_t MG_RESTRICT    = 1ULL<<38;
const std::uint64_t MG_PROLONG     = 1ULL<<39;
const std::uint64_t MG_FMGPROLONG  = 1ULL<<40;
} // namespace MultigridTaskNames

#endif // TASK_LIST_MG_TASK_LIST_HPP_
