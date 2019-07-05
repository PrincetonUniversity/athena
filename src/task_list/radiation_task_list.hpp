
#ifndef TASK_LIST_RADIATION_TASK_LIST_HPP_
#define TASK_LIST_RADIATION_TASK_LIST_HPP_
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//!   \file radiation_task_list.hpp
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

//--------------------------------------------------------------------------------------
//! \class RadiationIntegratorTaskList
//  \brief data and function definitions for RadiationIntegratorTaskList derived class
//
class RadiationIntegratorTaskList : public TaskList {
 public:
  RadiationIntegratorTaskList(ParameterInput *pin, Mesh *pm);
  std::string integrator;

  //functions
  enum TaskStatus LocalIntegratorJeans(MeshBlock *pmb, int step);
  enum TaskStatus ConstRadiation(MeshBlock *pmb, int step);
 private:
  void StartupTaskList(MeshBlock *pmb, int stage) override;
  void AddTask(const TaskID& id, const TaskID& dep) override;
};

namespace RadiationIntegratorTaskNames {
const TaskID NONE(0);
const TaskID INT_CONST(1); //constant radiation, do nothing
const TaskID INT_LOC_JEANS(2); //local jeans shielding
};

#endif // TASK_LIST_RADIATION_TASK_LIST_HPP_
