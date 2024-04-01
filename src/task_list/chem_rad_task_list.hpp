
#ifndef TASK_LIST_CHEM_RAD_TASK_LIST_HPP_
#define TASK_LIST_CHEM_RAD_TASK_LIST_HPP_
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//!   \file chem_rad_task_list.hpp
//    \brief provides functionality to control dynamic execution using tasks

// C headers

// C++ headers
#include <cstdint>      // std::uint64_t
#include <string>

// Athena++ headers
#include "../athena.hpp"
#include "task_list.hpp"

// forward declarations
class Mesh;
class MeshBlock;

//--------------------------------------------------------------------------------------
//! \class ChemRadiationIntegratorTaskList
//  \brief data and function definitions for ChemRadiationIntegratorTaskList derived class
//
class ChemRadiationIntegratorTaskList : public TaskList {
  friend class ChemRadIntegrator;
 public:
  ChemRadiationIntegratorTaskList(ParameterInput *pin, Mesh *pm);
  std::string integrator;

  // six-ray
  enum TaskStatus GetColMB_ix1(MeshBlock *pmb, int step);
  enum TaskStatus GetColMB_ox1(MeshBlock *pmb, int step);
  enum TaskStatus GetColMB_ix2(MeshBlock *pmb, int step);
  enum TaskStatus GetColMB_ox2(MeshBlock *pmb, int step);
  enum TaskStatus GetColMB_ix3(MeshBlock *pmb, int step);
  enum TaskStatus GetColMB_ox3(MeshBlock *pmb, int step);
  enum TaskStatus RecvAndSend_ix1(MeshBlock *pmb, int step);
  enum TaskStatus RecvAndSend_ox1(MeshBlock *pmb, int step);
  enum TaskStatus RecvAndSend_ix2(MeshBlock *pmb, int step);
  enum TaskStatus RecvAndSend_ox2(MeshBlock *pmb, int step);
  enum TaskStatus RecvAndSend_ix3(MeshBlock *pmb, int step);
  enum TaskStatus RecvAndSend_ox3(MeshBlock *pmb, int step);
  enum TaskStatus ClearSixrayReceive(MeshBlock *pmb, int step);
  enum TaskStatus UpdateRadiationSixRay(MeshBlock *pmb, int step);
  // constant radiation
  enum TaskStatus ConstRadiation(MeshBlock *pmb, int step);

 private:
  void StartupTaskList(MeshBlock *pmb, int stage) override;
  void AddTask(const TaskID& id, const TaskID& dep) override;
  enum TaskStatus RecvAndSend_direction(MeshBlock *pmb, int step, BoundaryFace direction);
};


namespace ChemRadiationIntegratorTaskNames {
const TaskID NONE(0);
const TaskID INT_CONST(1); // constant radiation, do nothing
// six-ray radiation
const TaskID GET_COL_MB_IX1(2);
const TaskID GET_COL_MB_OX1(3);
const TaskID GET_COL_MB_IX2(4);
const TaskID GET_COL_MB_OX2(5);
const TaskID GET_COL_MB_IX3(6);
const TaskID GET_COL_MB_OX3(7);
const TaskID RECV_SEND_COL_IX1(8);
const TaskID RECV_SEND_COL_OX1(9);
const TaskID RECV_SEND_COL_IX2(10);
const TaskID RECV_SEND_COL_OX2(11);
const TaskID RECV_SEND_COL_IX3(12);
const TaskID RECV_SEND_COL_OX3(13);
const TaskID CLEAR_SIXRAY_RECV(14);
const TaskID UPDATE_RAD(15);
} // namespace ChemRadiationIntegratorTaskNames

#endif // TASK_LIST_CHEM_RAD_TASK_LIST_HPP_
