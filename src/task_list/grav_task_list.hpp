#ifndef TASK_LIST_GRAV_TASK_LIST_HPP_
#define TASK_LIST_GRAV_TASK_LIST_HPP_
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//!   \file grav_task_list.hpp
//    \brief provides functionality to control dynamic execution using tasks

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
//! \class GravitySolverTaskList
//  \brief data and function definitions for GravitySolverTaskList derived class

class GravitySolverTaskList : public TaskList {
 public:
  GravitySolverTaskList(ParameterInput *pin, Mesh *pm);
  ~GravitySolverTaskList() {}

  void AddGravitySolverTask(std::uint64_t id, std::uint64_t dep);

  // functions
  TaskStatus ClearGravityBoundary(MeshBlock *pmb, int stage);
  TaskStatus SendGravityBoundary(MeshBlock *pmb, int stage);
  TaskStatus ReceiveGravityBoundary(MeshBlock *pmb, int stage);
  TaskStatus PhysicalBoundary(MeshBlock *pmb, int stage);

 private:
  void StartupTaskList(MeshBlock *pmb, int stage) override;
};


//----------------------------------------------------------------------------------------
// 64-bit integers with "1" in different bit positions used to ID  each hydro task.

namespace GravitySolverTaskNames {
const std::uint64_t NONE=0ULL;
const std::uint64_t CLEAR_GRAV=1ULL<<0;

const std::uint64_t SEND_GRAV_BND=1ULL<<1;
const std::uint64_t RECV_GRAV_BND=1ULL<<2;
const std::uint64_t GRAV_PHYS_BND=1ULL<<3;
}
#endif // TASK_LIST_GRAV_TASK_LIST_HPP_
