
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
#include <string>

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
  friend class RadIntegrator;
 public:
  RadiationIntegratorTaskList(ParameterInput *pin, Mesh *pm);
  std::string integrator;

  //six-ray
  enum TaskStatus GetColMB0(MeshBlock *pmb, int step);
  enum TaskStatus GetColMB1(MeshBlock *pmb, int step);
  enum TaskStatus GetColMB2(MeshBlock *pmb, int step);
  enum TaskStatus GetColMB3(MeshBlock *pmb, int step);
  enum TaskStatus GetColMB4(MeshBlock *pmb, int step);
  enum TaskStatus GetColMB5(MeshBlock *pmb, int step);
  enum TaskStatus RecvAndSend0(MeshBlock *pmb, int step);
  enum TaskStatus RecvAndSend1(MeshBlock *pmb, int step);
  enum TaskStatus RecvAndSend2(MeshBlock *pmb, int step);
  enum TaskStatus RecvAndSend3(MeshBlock *pmb, int step);
  enum TaskStatus RecvAndSend4(MeshBlock *pmb, int step);
  enum TaskStatus RecvAndSend5(MeshBlock *pmb, int step);
  enum TaskStatus ClearSixrayReceive(MeshBlock *pmb, int step);
  enum TaskStatus UpdateRadiationSixRay(MeshBlock *pmb, int step);
  //constant radiation
  enum TaskStatus ConstRadiation(MeshBlock *pmb, int step);
 private:
  void StartupTaskList(MeshBlock *pmb, int stage) override;
  void AddTask(const TaskID& id, const TaskID& dep) override;
  enum TaskStatus RecvAndSend_direction(MeshBlock *pmb, int step, int direction);
};

namespace RadiationIntegratorTaskNames {
const TaskID NONE(0);
const TaskID INT_CONST(1); //constant radiation, do nothing
//six-ray radiation
const TaskID GET_COL_MB0(2);
const TaskID GET_COL_MB1(3);
const TaskID GET_COL_MB2(4);
const TaskID GET_COL_MB3(5);
const TaskID GET_COL_MB4(6);
const TaskID GET_COL_MB5(7);
const TaskID RECV_SEND_COL0(8);
const TaskID RECV_SEND_COL1(9);
const TaskID RECV_SEND_COL2(10);
const TaskID RECV_SEND_COL3(11);
const TaskID RECV_SEND_COL4(12);
const TaskID RECV_SEND_COL5(13);
const TaskID CLEAR_SIXRAY_RECV(14);
const TaskID UPDATE_RAD(15);
};

#endif // TASK_LIST_RADIATION_TASK_LIST_HPP_
