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
#include "../bvals/bvals_grav.hpp"
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

FFTGravitySolverTaskList::FFTGravitySolverTaskList(ParameterInput *pin, Mesh *pm)
    : TaskList(pm) {
  // Now assemble list of tasks for each stage of time integrator
  {using namespace FFTGravitySolverTaskNames; // NOLINT (build/namespace)
    // compute hydro fluxes, integrate hydro variables
    AddFFTGravitySolverTask(SEND_GRAV_BND,NONE);
    AddFFTGravitySolverTask(RECV_GRAV_BND,NONE);
    AddFFTGravitySolverTask(GRAV_PHYS_BND,SEND_GRAV_BND|RECV_GRAV_BND);
    AddFFTGravitySolverTask(CLEAR_GRAV, GRAV_PHYS_BND);
  } // end of using namespace block
}

//----------------------------------------------------------------------------------------
//! \fn void FFTGravitySolverTaskList::AddFFTGravitySolverTask(std::uint64_t id,
//                                                       std::uint64_t dep)
//  \brief Sets id and dependency for "ntask" member of task_list_ array, then iterates
//  value of ntask.

void FFTGravitySolverTaskList::AddFFTGravitySolverTask(std::uint64_t id, std::uint64_t dep) {
  task_list_[ntasks].task_id=id;
  task_list_[ntasks].dependency=dep;

  using namespace FFTGravitySolverTaskNames; // NOLINT (build/namespace)
  switch (id) {
    case (CLEAR_GRAV):
      task_list_[ntasks].TaskFunc=
          static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
          (&FFTGravitySolverTaskList::ClearFFTGravityBoundary);
      break;
    case (SEND_GRAV_BND):
      task_list_[ntasks].TaskFunc=
          static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
          (&FFTGravitySolverTaskList::SendFFTGravityBoundary);
      break;
    case (RECV_GRAV_BND):
      task_list_[ntasks].TaskFunc=
          static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
          (&FFTGravitySolverTaskList::ReceiveFFTGravityBoundary);
      break;
    case (GRAV_PHYS_BND):
      task_list_[ntasks].TaskFunc=
          static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
          (&FFTGravitySolverTaskList::PhysicalBoundary);
      break;
    default:
      std::stringstream msg;
      msg << "### FATAL ERROR in AddFFTGravitySolverTask" << std::endl
          << "Invalid Task "<< id << " is specified" << std::endl;
      ATHENA_ERROR(msg);
  }
  ntasks++;
  return;
}


void FFTGravitySolverTaskList::StartupTaskList(MeshBlock *pmb, int stage) {
  pmb->pgbval->StartReceivingFFTGravity();

  return;
}

//----------------------------------------------------------------------------------------
//! \fn
//  \brief

//----------------------------------------------------------------------------------------
// Functions to start/end MPI communication

enum TaskStatus FFTGravitySolverTaskList::ClearFFTGravityBoundary(MeshBlock *pmb,
                                                                  int stage) {
  pmb->pgbval->ClearBoundaryFFTGravity();
  return TASK_SUCCESS;
}

enum TaskStatus FFTGravitySolverTaskList::SendFFTGravityBoundary(MeshBlock *pmb,
                                                                 int stage) {
  if (pmb->pgrav->pgbval->SendFFTGravityBoundaryBuffers(pmb->pgrav->phi)==false)
    return TASK_FAIL;
  return TASK_SUCCESS;
}

enum TaskStatus FFTGravitySolverTaskList::ReceiveFFTGravityBoundary(MeshBlock *pmb,
                                                                    int stage) {
  if (pmb->pgrav->pgbval->ReceiveFFTGravityBoundaryBuffers(pmb->pgrav->phi)==false)
    return TASK_FAIL;
  return TASK_SUCCESS;
}

enum TaskStatus FFTGravitySolverTaskList::PhysicalBoundary(MeshBlock *pmb,
                                                           int stage) {
  pmb->pgrav->pgbval->ApplyPhysicalBoundaries();
  return TASK_NEXT;
}
