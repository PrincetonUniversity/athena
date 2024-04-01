#ifndef TASK_LIST_CRDIFFUSION_TASK_LIST_HPP_
#define TASK_LIST_CRDIFFUSION_TASK_LIST_HPP_
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file grav_task_list.hpp
//! \brief define CRDiffusionBoundaryTaskList

// C headers

// C++ headers
#include <cstdint>      // std::uint64_t

// Athena++ headers
#include "../athena.hpp"
#include "task_list.hpp"

// forward declarations
class Mesh;
class MeshBlock;

//----------------------------------------------------------------------------------------
//! \class CRDiffusionBoundaryTaskList
//! \brief data and function definitions for CRDiffusionBoundaryTaskList derived class

class CRDiffusionBoundaryTaskList : public TaskList {
 public:
  CRDiffusionBoundaryTaskList(ParameterInput *pin, Mesh *pm);

  // functions
  TaskStatus ClearCRDiffusionBoundary(MeshBlock *pmb, int stage);
  TaskStatus SendCRDiffusionBoundary(MeshBlock *pmb, int stage);
  TaskStatus ReceiveCRDiffusionBoundary(MeshBlock *pmb, int stage);
  TaskStatus SetCRDiffusionBoundary(MeshBlock *pmb, int stage);
  TaskStatus ProlongateCRDiffusionBoundary(MeshBlock *pmb, int stage);
  TaskStatus PhysicalBoundary(MeshBlock *pmb, int stage);

 private:
  void AddTask(const TaskID& id, const TaskID& dep) override;
  void StartupTaskList(MeshBlock *pmb, int stage) override;
};


//----------------------------------------------------------------------------------------
//! 64-bit integers with "1" in different bit positions used to ID  each hydro task.
namespace CRDiffusionBoundaryTaskNames {
const TaskID NONE(0);
const TaskID CLEAR_CRDIFF(1);

const TaskID SEND_CRDIFF_BND(2);
const TaskID RECV_CRDIFF_BND(3);
const TaskID SETB_CRDIFF_BND(4);
const TaskID PROLONG_CRDIFF_BND(5);
const TaskID CRDIFF_PHYS_BND(6);
}
#endif // TASK_LIST_CRDIFFUSION_TASK_LIST_HPP_
