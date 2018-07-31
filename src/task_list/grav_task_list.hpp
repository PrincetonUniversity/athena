#ifndef TASK_LIST_GRAV_TASK_LIST_HPP_
#define TASK_LIST_GRAV_TASK_LIST_HPP_
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//!   \file grav_task_list.hpp
//    \brief provides functionality to control dynamic execution using tasks

#include <stdint.h>

// Athena++ headers
#include "../athena.hpp"
#include "task_list.hpp"

// forward declarations
class Mesh;
class MeshBlock;

//----------------------------------------------------------------------------------------
//! \class GravitySolverTaskList
//  \brief data and function definitions for GravitySolverTaskList derived class

class GravitySolverTaskList : public TaskList {
public:
  GravitySolverTaskList(ParameterInput *pin, Mesh *pm);
  ~GravitySolverTaskList() {}

  void AddGravitySolverTask(uint64_t id, uint64_t dep);

  // functions
  enum TaskStatus StartGravityReceive(MeshBlock *pmb, int stage);
  enum TaskStatus ClearGravityBoundary(MeshBlock *pmb, int stage);
  enum TaskStatus SendGravityBoundary(MeshBlock *pmb, int stage);
  enum TaskStatus ReceiveGravityBoundary(MeshBlock *pmb, int stage);
  enum TaskStatus PhysicalBoundary(MeshBlock *pmb, int stage);
};


//----------------------------------------------------------------------------------------
// 64-bit integers with "1" in different bit positions used to ID  each hydro task.

namespace GravitySolverTaskNames {
  const uint64_t NONE=0;
  const uint64_t START_GRAV_RECV=1LL<<0;
  const uint64_t CLEAR_GRAV=1LL<<1;

  const uint64_t SEND_GRAV_BND=1LL<<2;
  const uint64_t RECV_GRAV_BND=1LL<<3;
  const uint64_t GRAV_PHYS_BND=1LL<<4;
};

#endif // TASK_LIST_GRAV_TASK_LIST_HPP_
