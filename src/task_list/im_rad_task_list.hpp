#ifndef TASK_LIST_IM_RAD_TASK_LIST_HPP_
#define TASK_LIST_IM_RAD_TASK_LIST_HPP_
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file im_rad_task_list.hpp
//! \brief define im_rad_task_list class
//! \brief This is used to handle the boundary condition

// C headers

// C++ headers

// Athena++ headers
#include "../athena.hpp"
#include "./task_list.hpp"

// forward declarations
class Mesh;
class MeshBlock;
class Radiation;
class IMRadTaskList;
class TaskID;


//----------------------------------------------------------------------------------------
//! \struct IMRadTask
//! \brief data and function pointer for an individual IMRadTask

struct IMRadTask {
  TaskID task_id;      //!> encodes task using bit positions in IMRadTaskNames
  TaskID dependency;   //!> encodes dependencies to other tasks using IMRadTaskNames
  TaskStatus (IMRadTaskList::*TaskFunc)(MeshBlock *); //!> ptr to a task
};

//----------------------------------------------------------------------------------------
//! \class IMRadTaskList
//! \brief data and function definitions for IMRadTaskList class

class IMRadTaskList {
 public:
  IMRadTaskList() : ntasks(0), task_list_{} {} // 2x direct + zero initialization
  virtual ~IMRadTaskList() = default;

  Mesh *pmy_mesh;
  int ntasks;     //!> number of tasks in this list
  Real time, dt;

  TaskListStatus DoAllAvailableTasks(MeshBlock *pmb, TaskStates &ts);
  void DoTaskListOneStage(Real wght);
  TaskStatus PhysicalBoundary(MeshBlock *pmb);
  TaskStatus ProlongateBoundary(MeshBlock *pmb);

 protected:
  IMRadTask task_list_[64*TaskID::kNField_];

 private:
  virtual void AddTask(const TaskID& id, const TaskID& dep) = 0;
  virtual void StartupTaskList(MeshBlock *pmb) = 0;
};

class IMRadITTaskList : public IMRadTaskList {
 public:
  explicit IMRadITTaskList(Mesh *pm);

  TaskStatus ClearRadBoundary(MeshBlock *pmb);
  TaskStatus SendRadBoundary(MeshBlock *pmb);
  TaskStatus ReceiveRadBoundary(MeshBlock *pmb);
  TaskStatus SetRadBoundary(MeshBlock *pmb);
  TaskStatus SendRadBoundaryShear(MeshBlock *pmb);
  TaskStatus ReceiveRadBoundaryShear(MeshBlock *pmb);
  TaskStatus CheckResidual(MeshBlock *pmb);
  TaskStatus AddFluxAndSourceTerms(MeshBlock *pmb);

 private:
  void StartupTaskList(MeshBlock *pmb) override;
  void AddTask(const TaskID& id, const TaskID& dep) override;
};

//----------------------------------------------
//! IMRadHydroTaskList
//! Derived Class to handle Hydro boundary update
class IMRadHydroTaskList : public IMRadTaskList{
 public:
  explicit IMRadHydroTaskList(Mesh *pm);

  TaskStatus ClearHydroBoundary(MeshBlock *pmb);
  TaskStatus SendHydroBoundary(MeshBlock *pmb);
  TaskStatus ReceiveHydroBoundary(MeshBlock *pmb);
  TaskStatus SetHydroBoundary(MeshBlock *pmb);
  TaskStatus SendHydroBoundaryShear(MeshBlock *pmb);
  TaskStatus ReceiveHydroBoundaryShear(MeshBlock *pmb);
  TaskStatus UpdateOpacity(MeshBlock *pmb);
  TaskStatus AddRadSource(MeshBlock *pmb);
  TaskStatus Primitive(MeshBlock *pmb);

 private:
  void StartupTaskList(MeshBlock *pmb) override;
  void AddTask(const TaskID& id, const TaskID& dep) override;
};

// task list for separate compton scattering
// need to update boundary condition
class IMRadComptTaskList : public IMRadTaskList {
 public:
  explicit IMRadComptTaskList(Mesh *pm);

  TaskStatus ClearRadBoundary(MeshBlock *pmb);
  TaskStatus SendRadBoundary(MeshBlock *pmb);
  TaskStatus ReceiveRadBoundary(MeshBlock *pmb);
  TaskStatus SetRadBoundary(MeshBlock *pmb);
  TaskStatus SendRadBoundaryShear(MeshBlock *pmb);
  TaskStatus ReceiveRadBoundaryShear(MeshBlock *pmb);
  TaskStatus CalComptTerms(MeshBlock *pmb);
 private:
  void StartupTaskList(MeshBlock *pmb) override;
  void AddTask(const TaskID& id, const TaskID& dep) override;
};


//----------------------------------------------------------------------------------------
//! 64-bit integers with "1" in different bit positions used to ID each IMRadITT task.

namespace IMRadITTaskNames {
const TaskID NONE(0);
const TaskID CLEAR_RAD(1);     // clear radiation boundary
const TaskID SEND_RAD_BND(2);  // send radiation boundary
const TaskID RECV_RAD_BND(3); // receive radiation boundary
const TaskID SETB_RAD_BND(4); // set radiation physical boundary
const TaskID RAD_PHYS_BND(5);// radiation physical boundary
const TaskID PRLN_RAD_BND(6); // prolongation
const TaskID SEND_RAD_SH(7); // send shearing box boundary
const TaskID RECV_RAD_SH(8); // receive shearing box boundary
const TaskID CHK_RAD_RES(9); // check residual
const TaskID FLX_AND_SRC(10);  // calculate the source term and all others together

} // namespace IMRadITTaskNames

namespace IMRadHydroTaskNames {
const TaskID NONE(0);
const TaskID CLEAR_HYD(1);    // clear radiation boundary
const TaskID SEND_HYD_BND(2); // send radiation boundary
const TaskID RECV_HYD_BND(3); // receive radiation boundary
const TaskID SETB_HYD_BND(4); // set radiation physical boundary
const TaskID HYD_PHYS_BND(5); // radiation physical boundary
const TaskID PRLN_HYD_BND(6); // prolongation
const TaskID SEND_HYD_SH(7);  // send shearing box boundary
const TaskID RECV_HYD_SH(8);  // receive shearing box boundary
const TaskID UPD_OPA(9);      // check residual
const TaskID ADD_RAD_SRC(10); // add radiation source term
const TaskID CONS_TO_PRIM(11);// convert conservative to primitive variables
} // namespace IMRadHydroTaskNames

namespace IMRadComptTaskNames {
const TaskID NONE(0);
const TaskID CLEAR_RAD(1);     // clear radiation boundary
const TaskID SEND_RAD_BND(2);  // send radiation boundary
const TaskID RECV_RAD_BND(3); // receive radiation boundary
const TaskID SETB_RAD_BND(4); // set radiation physical boundary
const TaskID RAD_PHYS_BND(5);// radiation physical boundary
const TaskID PRLN_RAD_BND(6); // prolongation
const TaskID SEND_RAD_SH(7); // send shearing box boundary
const TaskID RECV_RAD_SH(8); // receive shearing box boundary
const TaskID CAL_COMPT(9); // add flux divergence term

} // namespace IMRadComptTaskNames

#endif // TASK_LIST_IM_RAD_TASK_LIST_HPP_
