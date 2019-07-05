#ifndef TASK_LIST_CHEMISTRY_TASK_LIST_HPP_
#define TASK_LIST_CHEMISTRY_TASK_LIST_HPP_
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//!   \file chemistry_task_list.hpp
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
//! \class ChemistryIntegratorTaskList
//  \brief data and function definitions for ChemistryIntegratorTaskList derived class
//
class ChemistryIntegratorTaskList : public TaskList {
 public:
  ChemistryIntegratorTaskList(ParameterInput *pin, Mesh *pm);

  //functions
  TaskStatus IntegrateSourceTerm(MeshBlock *pmb, int step);

 private:
  void StartupTaskList(MeshBlock *pmb, int stage) override;
  void AddTask(const TaskID& id, const TaskID& dep) override;
};


namespace ChemistryIntegratorTaskNames {
const TaskID NONE(0);
const TaskID INT_CHEM_SRC(1);
};

#endif // TASK_LIST_CHEMISTRY_TASK_LIST_HPP_
