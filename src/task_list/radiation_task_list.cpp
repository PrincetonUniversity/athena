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
      AddTask(GET_COL_MB_IX1,NONE);
      AddTask(RECV_SEND_COL_IX1,GET_COL_MB_IX1);
      AddTask(GET_COL_MB_OX1,NONE);
      AddTask(RECV_SEND_COL_OX1,GET_COL_MB_OX1);
      AddTask(GET_COL_MB_IX2,NONE);
      AddTask(RECV_SEND_COL_IX2,GET_COL_MB_IX2);
      AddTask(GET_COL_MB_OX2,NONE);
      AddTask(RECV_SEND_COL_OX2,GET_COL_MB_OX2);
      AddTask(GET_COL_MB_IX3,NONE);
      AddTask(RECV_SEND_COL_IX3,GET_COL_MB_IX3);
      AddTask(GET_COL_MB_OX3,NONE);
      AddTask(RECV_SEND_COL_OX3,GET_COL_MB_OX3);
      AddTask(UPDATE_RAD,
          RECV_SEND_COL_IX1|RECV_SEND_COL_OX1|RECV_SEND_COL_IX2|
          RECV_SEND_COL_OX2|RECV_SEND_COL_IX3|RECV_SEND_COL_OX3);
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
  } else if (id == GET_COL_MB_IX1) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&RadiationIntegratorTaskList::GetColMB_ix1);
  } else if (id == GET_COL_MB_OX1) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&RadiationIntegratorTaskList::GetColMB_ox1);
  } else if (id == GET_COL_MB_IX2) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&RadiationIntegratorTaskList::GetColMB_ix2);
  } else if (id == GET_COL_MB_OX2) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&RadiationIntegratorTaskList::GetColMB_ox2);
  } else if (id == GET_COL_MB_IX3) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&RadiationIntegratorTaskList::GetColMB_ix3);
  } else if (id == GET_COL_MB_OX3) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&RadiationIntegratorTaskList::GetColMB_ox3);
  } else if (id == RECV_SEND_COL_IX1) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&RadiationIntegratorTaskList::RecvAndSend_ix1);
  } else if (id == RECV_SEND_COL_OX1) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&RadiationIntegratorTaskList::RecvAndSend_ox1);
  } else if (id == RECV_SEND_COL_IX2) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&RadiationIntegratorTaskList::RecvAndSend_ix2);
  } else if (id == RECV_SEND_COL_OX2) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&RadiationIntegratorTaskList::RecvAndSend_ox2);
  } else if (id == RECV_SEND_COL_IX3) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&RadiationIntegratorTaskList::RecvAndSend_ix3);
  } else if (id == RECV_SEND_COL_OX3) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&RadiationIntegratorTaskList::RecvAndSend_ox3);
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
#ifdef INCLUDE_CHEMISTRY
  if (integrator == "six_ray") {
    pmb->prad->pradintegrator->col_bvar.StartReceiving(BoundaryCommSubset::all);
  }
#endif
  return;
}

//----------------------------------------------------------------------------------------
//! Six-ray

//meshblock column densities
TaskStatus RadiationIntegratorTaskList::GetColMB_ix1(MeshBlock *pmb, int step)
{
#ifdef INCLUDE_CHEMISTRY
  pmb->prad->pradintegrator->GetColMB(BoundaryFace::inner_x1);
#endif
  return TaskStatus::success;
}

TaskStatus RadiationIntegratorTaskList::GetColMB_ox1(MeshBlock *pmb, int step)
{
#ifdef INCLUDE_CHEMISTRY
  pmb->prad->pradintegrator->GetColMB(BoundaryFace::outer_x1);
#endif
  return TaskStatus::success;
}

TaskStatus RadiationIntegratorTaskList::GetColMB_ix2(MeshBlock *pmb, int step)
{
#ifdef INCLUDE_CHEMISTRY
  pmb->prad->pradintegrator->GetColMB(BoundaryFace::inner_x2);
#endif
  return TaskStatus::success;
}

TaskStatus RadiationIntegratorTaskList::GetColMB_ox2(MeshBlock *pmb, int step)
{
#ifdef INCLUDE_CHEMISTRY
  pmb->prad->pradintegrator->GetColMB(BoundaryFace::outer_x2);
#endif
  return TaskStatus::success;
}

TaskStatus RadiationIntegratorTaskList::GetColMB_ix3(MeshBlock *pmb, int step)
{
#ifdef INCLUDE_CHEMISTRY
  pmb->prad->pradintegrator->GetColMB(BoundaryFace::inner_x3);
#endif
  return TaskStatus::success;
}

TaskStatus RadiationIntegratorTaskList::GetColMB_ox3(MeshBlock *pmb, int step)
{
#ifdef INCLUDE_CHEMISTRY
  pmb->prad->pradintegrator->GetColMB(BoundaryFace::outer_x3);
#endif
  return TaskStatus::success;
}

//boundary receive and send

TaskStatus RadiationIntegratorTaskList::RecvAndSend_direction(MeshBlock *pmb,
    int step, BoundaryFace direction)
{
#ifdef INCLUDE_CHEMISTRY
  SixRayBoundaryVariable *pbvar = &pmb->prad->pradintegrator->col_bvar;
  BoundaryFace direction_opp = pbvar->GetOppositeBoundaryFace(direction);
  NeighborBlock *pnb = pbvar->GetFaceNeighbor(direction);
  NeighborBlock *pnb_opp = pbvar->GetFaceNeighbor(direction_opp);
  bool ret = true;
  if (pnb == nullptr) {
    if (pnb_opp != nullptr) {
      pbvar->SendSixRayBoundaryBuffers(direction_opp);
    }
  } else {
    ret = pbvar->ReceiveAndSetSixRayBoundaryBuffers(direction);
    if (ret == true) {
      pmb->prad->pradintegrator->UpdateCol(direction);
      if (pnb_opp != nullptr) {
        pbvar->SendSixRayBoundaryBuffers(direction_opp);
      }
    } else {
      return TaskStatus::fail;
    }
  }
#endif
  return TaskStatus::success;
}

TaskStatus RadiationIntegratorTaskList::RecvAndSend_ix1(MeshBlock *pmb, int step){
  return RecvAndSend_direction(pmb, step, BoundaryFace::inner_x1);
}

TaskStatus RadiationIntegratorTaskList::RecvAndSend_ox1(MeshBlock *pmb, int step){
  return RecvAndSend_direction(pmb, step, BoundaryFace::outer_x1);
}

TaskStatus RadiationIntegratorTaskList::RecvAndSend_ix2(MeshBlock *pmb, int step){
  return RecvAndSend_direction(pmb, step, BoundaryFace::inner_x2);
}

TaskStatus RadiationIntegratorTaskList::RecvAndSend_ox2(MeshBlock *pmb, int step){
  return RecvAndSend_direction(pmb, step, BoundaryFace::outer_x2);
}

TaskStatus RadiationIntegratorTaskList::RecvAndSend_ix3(MeshBlock *pmb, int step){
  return RecvAndSend_direction(pmb, step, BoundaryFace::inner_x3);
}

TaskStatus RadiationIntegratorTaskList::RecvAndSend_ox3(MeshBlock *pmb, int step){
  return RecvAndSend_direction(pmb, step, BoundaryFace::outer_x3);
}

TaskStatus RadiationIntegratorTaskList::ClearSixrayReceive(MeshBlock *pmb,
                                                                 int step)
{
#ifdef INCLUDE_CHEMISTRY
  pmb->prad->pradintegrator->col_bvar.ClearBoundary(BoundaryCommSubset::all);
#endif
  return TaskStatus::success;
}

//update radiation variables
TaskStatus RadiationIntegratorTaskList::UpdateRadiationSixRay(MeshBlock *pmb,
                                                              int step)
{
#ifdef INCLUDE_CHEMISTRY
  pmb->prad->pradintegrator->UpdateRadiation();
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

