//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file radiation_task_list.cpp
//! \brief derived class for radiation integrator task list.

// C/C++ headers
#include <iostream>   // endl
#include <sstream>    // sstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena++ classes headers
#include "../athena.hpp"
#include "../defs.hpp"
#include "../mesh/mesh.hpp"
#include "../radiation/integrators/rad_integrators.hpp"
#include "../radiation/radiation.hpp"
#include "radiation_task_list.hpp"
#include "task_list.hpp"

//--------------------------------------------------------------------------------------
//! RadiationIntegratorTaskList constructor
RadiationIntegratorTaskList::RadiationIntegratorTaskList(ParameterInput *pin, Mesh *pm) {
  integrator = RADIATION_INTEGRATOR;
  // Now assemble list of tasks for each step of chemistry integrator
  {using namespace RadiationIntegratorTaskNames; // NOLINT (build/namespace)
    if (integrator == "six_ray") {
      AddTask(GET_COL_MB0,NONE);
      AddTask(RECV_SEND_COL0,GET_COL_MB0);
      AddTask(GET_COL_MB1,NONE);
      AddTask(RECV_SEND_COL1,GET_COL_MB1);
      AddTask(GET_COL_MB2,NONE);
      AddTask(RECV_SEND_COL2,GET_COL_MB2);
      AddTask(GET_COL_MB3,NONE);
      AddTask(RECV_SEND_COL3,GET_COL_MB3);
      AddTask(GET_COL_MB4,NONE);
      AddTask(RECV_SEND_COL4,GET_COL_MB4);
      AddTask(GET_COL_MB5,NONE);
      AddTask(RECV_SEND_COL5,GET_COL_MB5);
      AddTask(UPDATE_RAD,
          RECV_SEND_COL0|RECV_SEND_COL1|RECV_SEND_COL2|
          RECV_SEND_COL3|RECV_SEND_COL4|RECV_SEND_COL5);
      AddTask(CLEAR_SIXRAY_RECV,UPDATE_RAD);
    } else if (integrator == "const") {
      //do nothing, radiation field constant, remain initial value
      AddTask(INT_CONST,NONE);
    } else {
      std::stringstream msg;
      msg << "### FATAL ERROR in RadiationIntegratorTaskList constructor" << std::endl
        << "integrator=" << integrator << " not valid radiation integrator, "
        << std::endl << "choose from {six_ray, const}" << std::endl;
      ATHENA_ERROR(msg);
    }
  } // end of using namespace block
}

//--------------------------------------------------------------------------------------
//! \fn void RadiationIntegratorTaskList::AddTask(const TaskID& id, const TaskID& dep)
//! \brief Sets id and dependency for "ntask" member of task_list_ array, then iterates
//! value of ntask.
void RadiationIntegratorTaskList::AddTask(const TaskID& id, const TaskID& dep) {
  task_list_[ntasks].task_id=id;
  task_list_[ntasks].dependency=dep;

  using namespace RadiationIntegratorTaskNames; // NOLINT (build/namespace)
  if (id == UPDATE_RAD) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&RadiationIntegratorTaskList::UpdateRadiationSixRay);
  } else if (id == GET_COL_MB0) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&RadiationIntegratorTaskList::GetColMB0);
  } else if (id == GET_COL_MB1) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&RadiationIntegratorTaskList::GetColMB1);
  } else if (id == GET_COL_MB2) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&RadiationIntegratorTaskList::GetColMB2);
  } else if (id == GET_COL_MB3) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&RadiationIntegratorTaskList::GetColMB3);
  } else if (id == GET_COL_MB4) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&RadiationIntegratorTaskList::GetColMB4);
  } else if (id == GET_COL_MB5) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&RadiationIntegratorTaskList::GetColMB5);
  } else if (id == RECV_SEND_COL0) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&RadiationIntegratorTaskList::RecvAndSend0);
  } else if (id == RECV_SEND_COL1) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&RadiationIntegratorTaskList::RecvAndSend1);
  } else if (id == RECV_SEND_COL2) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&RadiationIntegratorTaskList::RecvAndSend2);
  } else if (id == RECV_SEND_COL3) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&RadiationIntegratorTaskList::RecvAndSend3);
  } else if (id == RECV_SEND_COL4) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&RadiationIntegratorTaskList::RecvAndSend4);
  } else if (id == RECV_SEND_COL5) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&RadiationIntegratorTaskList::RecvAndSend5);
  } else if (id == CLEAR_SIXRAY_RECV) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&RadiationIntegratorTaskList::ClearSixrayReceive);
  } else if (id == INT_CONST) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&RadiationIntegratorTaskList::ConstRadiation);
  } else {
    std::stringstream msg;
    msg << "### FATAL ERROR in RadiationIntegratorTaskList::AddTask" << std::endl
      << "Invalid Task is specified" << std::endl;
    ATHENA_ERROR(msg);
  }
  ntasks++;
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationIntegratorTaskList::StartupTaskList(MeshBlock *pmb, int stage)
//! \brief Initialize boundary
void RadiationIntegratorTaskList::StartupTaskList(MeshBlock *pmb, int stage) {
  //TODO(Munan Gong): sixray boundary e.g.
  //for fft gravity: pmb->pgrav->gbvar.StartReceiving(BoundaryCommSubset::all);
  return;
}

//----------------------------------------------------------------------------------------
//! Six-ray

//meshblock column densities
TaskStatus RadiationIntegratorTaskList::GetColMB0(MeshBlock *pmb, int step)
{
#ifdef INCLUDE_CHEMISTRY
  pmb->prad->pradintegrator->GetColMB(0);
#endif
  return TaskStatus::success;
}

TaskStatus RadiationIntegratorTaskList::GetColMB1(MeshBlock *pmb, int step)
{
#ifdef INCLUDE_CHEMISTRY
  pmb->prad->pradintegrator->GetColMB(1);
#endif
  return TaskStatus::success;
}

TaskStatus RadiationIntegratorTaskList::GetColMB2(MeshBlock *pmb, int step)
{
#ifdef INCLUDE_CHEMISTRY
  pmb->prad->pradintegrator->GetColMB(2);
#endif
  return TaskStatus::success;
}

TaskStatus RadiationIntegratorTaskList::GetColMB3(MeshBlock *pmb, int step)
{
#ifdef INCLUDE_CHEMISTRY
  pmb->prad->pradintegrator->GetColMB(3);
#endif
  return TaskStatus::success;
}

TaskStatus RadiationIntegratorTaskList::GetColMB4(MeshBlock *pmb, int step)
{
#ifdef INCLUDE_CHEMISTRY
  pmb->prad->pradintegrator->GetColMB(4);
#endif
  return TaskStatus::success;
}

TaskStatus RadiationIntegratorTaskList::GetColMB5(MeshBlock *pmb, int step)
{
#ifdef INCLUDE_CHEMISTRY
  pmb->prad->pradintegrator->GetColMB(5);
#endif
  return TaskStatus::success;
}

//boundary receive and send

TaskStatus RadiationIntegratorTaskList::RecvAndSend_direction(MeshBlock *pmb,
    int step, int direction)
{
#ifdef INCLUDE_CHEMISTRY
  //TODO (Munan Gong)
#endif
  return TaskStatus::success;
}

TaskStatus RadiationIntegratorTaskList::RecvAndSend0(MeshBlock *pmb, int step){
  return RecvAndSend_direction(pmb, step, 0);
}

TaskStatus RadiationIntegratorTaskList::RecvAndSend1(MeshBlock *pmb, int step){
  return RecvAndSend_direction(pmb, step, 1);
}

TaskStatus RadiationIntegratorTaskList::RecvAndSend2(MeshBlock *pmb, int step){
  return RecvAndSend_direction(pmb, step, 2);
}

TaskStatus RadiationIntegratorTaskList::RecvAndSend3(MeshBlock *pmb, int step){
  return RecvAndSend_direction(pmb, step, 3);
}

TaskStatus RadiationIntegratorTaskList::RecvAndSend4(MeshBlock *pmb, int step){
  return RecvAndSend_direction(pmb, step, 4);
}

TaskStatus RadiationIntegratorTaskList::RecvAndSend5(MeshBlock *pmb, int step){
  return RecvAndSend_direction(pmb, step, 5);
}

TaskStatus RadiationIntegratorTaskList::ClearSixrayReceive(MeshBlock *pmb,
                                                                 int step)
{
#ifdef INCLUDE_CHEMISTRY
  //TODO (Munan Gong)
#endif
  return TaskStatus::success;
}

//update radiation variables
TaskStatus RadiationIntegratorTaskList::UpdateRadiationSixRay(MeshBlock *pmb,
                                                              int step)
{
#ifdef INCLUDE_CHEMISTRY
  for (int i=0; i<6; i++) {
    pmb->prad->pradintegrator->UpdateRadiation(i);
  }
  pmb->prad->pradintegrator->CopyToOutput();
#endif
  return TaskStatus::success;
}


//----------------------------------------------------------------------------------------
//! Trivial constant integrator
TaskStatus RadiationIntegratorTaskList::ConstRadiation(MeshBlock *pmb, int stage) {
#ifdef INCLUDE_CHEMISTRY
  pmb->prad->pradintegrator->CopyToOutput();
#endif
  return TaskStatus::success;
}

