#ifndef TASK_LIST_MG_TASK_LIST_HPP_
#define TASK_LIST_MG_TASK_LIST_HPP_
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file mg_task_list.hpp
//! \brief define MultiGridTaskList class

// C headers

// C++ headers

// Athena++ headers
#include "../athena.hpp"
#include "./task_list.hpp"

// forward declarations
class Mesh;
class MeshBlock;
class Multigrid;
class MultigridDriver;
class MultigridTaskList;

//----------------------------------------------------------------------------------------
//! \struct MGTask
//! \brief data and function pointer for an individual MGTask

struct MGTask {
  TaskID task_id;      //!> encodes task using bit positions in MultigridTaskNames
  TaskID dependency;   //!> encodes dependencies to other tasks using MultigridTaskNames
  TaskStatus (MultigridTaskList::*TaskFunc)(Multigrid*);  //!> ptr to a task
};


//----------------------------------------------------------------------------------------
//! \class MultigridTaskList
//! \brief data and function definitions for MultigridTaskList class

class MultigridTaskList {
 public:
  explicit MultigridTaskList(MultigridDriver *pmd) : ntasks(0), pmy_mgdriver_(pmd),
                                                     task_list_{} {}
  // data
  int ntasks;     //!> number of tasks in this list

  // functions
  TaskListStatus DoAllAvailableTasks(Multigrid *pmg, TaskStates &ts);
  void DoTaskListOneStage(MultigridDriver *pmd);
  void ClearTaskList() {ntasks=0;}

  // functions
  TaskStatus StartReceive(Multigrid *pmg);
  TaskStatus StartReceiveFluxCons(Multigrid *pmg);
  TaskStatus StartReceiveForProlongation(Multigrid *pmg);
  TaskStatus ClearBoundary(Multigrid *pmg);
  TaskStatus ClearBoundaryFluxCons(Multigrid *pmg);
  TaskStatus SendBoundary(Multigrid *pmg);
  TaskStatus SendBoundaryFluxCons(Multigrid *pmg);
  TaskStatus SendBoundaryForProlongation(Multigrid *pmg);
  TaskStatus ReceiveBoundary(Multigrid *pmg);
  TaskStatus ReceiveBoundaryFluxCons(Multigrid *pmg);
  TaskStatus ReceiveBoundaryForProlongation(Multigrid *pmg);
  TaskStatus SmoothRed(Multigrid *pmg);
  TaskStatus SmoothBlack(Multigrid *pmg);
  TaskStatus PhysicalBoundary(Multigrid *pmg);
  TaskStatus Restrict(Multigrid *pmg);
  TaskStatus Prolongate(Multigrid *pmg);
  TaskStatus FMGProlongate(Multigrid *pmg);
  TaskStatus ProlongateBoundary(Multigrid *pmg);
  TaskStatus ProlongateBoundaryForProlongation(Multigrid *pmg);
  TaskStatus CalculateFASRHS(Multigrid *pmg);
  TaskStatus StoreOldData(Multigrid *pmg);

  void SetMGTaskListToFiner(int nsmooth, int ngh, int flag = 0);
  void SetMGTaskListToCoarser(int nsmooth, int ngh);
  void SetMGTaskListFMGProlongate(int flag = 0);

 private:
  MultigridDriver* pmy_mgdriver_;
  MGTask task_list_[64*TaskID::kNField_];

  void AddMultigridTask(const TaskID& id, const TaskID& dep);
};

//----------------------------------------------------------------------------------------
//! 64-bit integers with "1" in different bit positions used to ID each Multigrid task.

namespace MultigridTaskNames {
const TaskID NONE(0);
const TaskID MG_STARTRECV0(1);
const TaskID MG_STARTRECV1R(2);
const TaskID MG_STARTRECV1B(3);
const TaskID MG_STARTRECV2R(4);
const TaskID MG_STARTRECV2B(5);
const TaskID MG_STARTRECVP(6);
const TaskID MG_STARTRECVL(7);
const TaskID MG_CLEARBND0(8);
const TaskID MG_CLEARBND1R(9);
const TaskID MG_CLEARBND1B(10);
const TaskID MG_CLEARBND2R(11);
const TaskID MG_CLEARBND2B(12);
const TaskID MG_CLEARBNDP(13);
const TaskID MG_CLEARBNDL(14);
const TaskID MG_SENDBND0(15);
const TaskID MG_SENDBND1R(16);
const TaskID MG_SENDBND1B(17);
const TaskID MG_SENDBND2R(18);
const TaskID MG_SENDBND2B(19);
const TaskID MG_SENDBNDP(20);
const TaskID MG_SENDBNDL(21);
const TaskID MG_RECVBND0(22);
const TaskID MG_RECVBND1R(23);
const TaskID MG_RECVBND1B(24);
const TaskID MG_RECVBND2R(25);
const TaskID MG_RECVBND2B(26);
const TaskID MG_RECVBNDP(27);
const TaskID MG_RECVBNDL(28);
const TaskID MG_PRLNGBNDP(29);
const TaskID MG_PRLNGFC0(30);
const TaskID MG_PRLNGFC1R(31);
const TaskID MG_PRLNGFC1B(32);
const TaskID MG_PRLNGFC2R(33);
const TaskID MG_PRLNGFC2B(34);
const TaskID MG_PRLNGFCL(35);
const TaskID MG_SMOOTH1R(36);
const TaskID MG_SMOOTH1B(37);
const TaskID MG_SMOOTH2R(38);
const TaskID MG_SMOOTH2B(39);
const TaskID MG_PHYSBND0(40);
const TaskID MG_PHYSBND1R(41);
const TaskID MG_PHYSBND1B(42);
const TaskID MG_PHYSBND2R(43);
const TaskID MG_PHYSBND2B(44);
const TaskID MG_PHYSBNDL(45);
const TaskID MG_RESTRICT(46);
const TaskID MG_PROLONG(47);
const TaskID MG_FMGPROLONG(48);
const TaskID MG_CALCFASRHS(49);
} // namespace MultigridTaskNames

#endif // TASK_LIST_MG_TASK_LIST_HPP_
