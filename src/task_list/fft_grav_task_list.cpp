//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file fft_grav_task_list.cpp
//  \brief

// C headers

// C++ headers
#include <iostream>   // endl
#include <sstream>    // sstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena++ headers
#include "../athena.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../gravity/gravity.hpp"
#include "../hydro/hydro.hpp"
#include "../hydro/srcterms/hydro_srcterms.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "fft_grav_task_list.hpp"
#include "task_list.hpp"

//----------------------------------------------------------------------------------------
//  FFTGravitySolverTaskList constructor

FFTGravitySolverTaskList::FFTGravitySolverTaskList(ParameterInput *pin, Mesh *pm) {
  // Now assemble list of tasks for each stage of time integrator
  {using namespace FFTGravitySolverTaskNames; // NOLINT (build/namespace)
    // compute hydro fluxes, integrate hydro variables
    AddTask(SEND_GRAV_BND,NONE);
    AddTask(RECV_GRAV_BND,NONE);
    AddTask(SETB_GRAV_BND,(RECV_GRAV_BND|SEND_GRAV_BND));
    AddTask(GRAV_PHYS_BND,SETB_GRAV_BND);
    AddTask(CLEAR_GRAV, GRAV_PHYS_BND);
  } // end of using namespace block
}

//----------------------------------------------------------------------------------------
//! \fn void FFTGravitySolverTaskList::AddTask(std::uint64_t id, std::uint64_t dep)
//  \brief Sets id and dependency for "ntask" member of task_list_ array, then iterates
//  value of ntask.

void FFTGravitySolverTaskList::AddTask(std::uint64_t id, std::uint64_t dep) {
  task_list_[ntasks].task_id=id;
  task_list_[ntasks].dependency=dep;

  using namespace FFTGravitySolverTaskNames; // NOLINT (build/namespace)
  switch (id) {
    case (CLEAR_GRAV):
      task_list_[ntasks].TaskFunc=
          static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
          (&FFTGravitySolverTaskList::ClearFFTGravityBoundary);
      break;
    case (SEND_GRAV_BND):
      task_list_[ntasks].TaskFunc=
          static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
          (&FFTGravitySolverTaskList::SendFFTGravityBoundary);
      break;
    case (RECV_GRAV_BND):
      task_list_[ntasks].TaskFunc=
          static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
          (&FFTGravitySolverTaskList::ReceiveFFTGravityBoundary);
      break;
    case (SETB_GRAV_BND):
      task_list_[ntasks].TaskFunc=
          static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
          (&FFTGravitySolverTaskList::SetFFTGravityBoundary);
      break;
    case (GRAV_PHYS_BND):
      task_list_[ntasks].TaskFunc=
          static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
          (&FFTGravitySolverTaskList::PhysicalBoundary);
      break;
    default:
      std::stringstream msg;
      msg << "### FATAL ERROR in FFTGravitySolverTaskList::AddTask" << std::endl
          << "Invalid Task "<< id << " is specified" << std::endl;
      ATHENA_ERROR(msg);
  }
  ntasks++;
  return;
}

void FFTGravitySolverTaskList::StartupTaskList(MeshBlock *pmb, int stage) {
  pmb->pgrav->pgbval->StartReceiving(BoundaryCommSubset::all);
  return;
}

TaskStatus FFTGravitySolverTaskList::ClearFFTGravityBoundary(MeshBlock *pmb, int stage) {
  pmb->pgrav->pgbval->ClearBoundary(BoundaryCommSubset::all);
  return TaskStatus::success;
}

TaskStatus FFTGravitySolverTaskList::SendFFTGravityBoundary(MeshBlock *pmb, int stage) {
  pmb->pgrav->pgbval->SendBoundaryBuffers();
  return TaskStatus::success;
}

TaskStatus FFTGravitySolverTaskList::ReceiveFFTGravityBoundary(MeshBlock *pmb,
                                                               int stage) {
  bool ret = pmb->pgrav->pgbval->ReceiveBoundaryBuffers();
  if (ret == false)
    return TaskStatus::fail;
  return TaskStatus::success;
}

TaskStatus FFTGravitySolverTaskList::SetFFTGravityBoundary(MeshBlock *pmb, int stage) {
  pmb->pgrav->pgbval->SetBoundaries();
  return TaskStatus::success;
}

TaskStatus FFTGravitySolverTaskList::PhysicalBoundary(MeshBlock *pmb, int stage) {
  return TaskStatus::next;
}
