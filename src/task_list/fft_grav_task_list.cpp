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

FFTGravitySolverTaskList::FFTGravitySolverTaskList(ParameterInput *pin, Mesh *pm)
    : TaskList(pm) {
  // Now assemble list of tasks for each stage of time integrator
  {using namespace FFTGravitySolverTaskNames; // NOLINT (build/namespace)
    // compute hydro fluxes, integrate hydro variables
    AddTask(SEND_GRAV_BND,NONE);
    AddTask(RECV_GRAV_BND,NONE);
    AddTask(GRAV_PHYS_BND,SEND_GRAV_BND|RECV_GRAV_BND);
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
      msg << "### FATAL ERROR in FFTGravitySolverTaskList::AddTask" << std::endl
          << "Invalid Task "<< id << " is specified" << std::endl;
      ATHENA_ERROR(msg);
  }
  ntasks++;
  return;
}

void FFTGravitySolverTaskList::StartupTaskList(MeshBlock *pmb, int stage) {
  // KGF: BoundaryValues wrapper version called in time_integrator.cpp
  // KGF: what is "time" parameter ever used for in this function?
  // ANSWER: shearing box capabilities, only. Remove this function parameter.
  Real time=0;
  pmb->pgrav->pgbval->StartReceivingAll(time);
  return;
}

enum TaskStatus FFTGravitySolverTaskList::ClearFFTGravityBoundary(MeshBlock *pmb,
                                                                  int stage) {
  // KGF: BoundaryValues wrapper version called in time_integrator.cpp
  pmb->pgrav->pgbval->ClearBoundaryAll();
  return TASK_SUCCESS;
}

enum TaskStatus FFTGravitySolverTaskList::SendFFTGravityBoundary(MeshBlock *pmb,
                                                                 int stage) {
  // KGF: BoundaryBuffer version does not return bool (copied Multigrid implementation)
  pmb->pgrav->pgbval->SendBoundaryBuffers();
  return TASK_SUCCESS;
}

enum TaskStatus FFTGravitySolverTaskList::ReceiveFFTGravityBoundary(MeshBlock *pmb,
                                                                    int stage) {
  if (pmb->pgrav->pgbval->ReceiveBoundaryBuffers() == false)
    return TASK_FAIL;
  return TASK_SUCCESS;
}

enum TaskStatus FFTGravitySolverTaskList::PhysicalBoundary(MeshBlock *pmb,
                                                           int stage) {
  // KGF: FFT self-gravity can only handle periodic boundary conditions

  // KGF: because Gravity() currently does not add the pointer:
  // pgbval = new CellCenteredBoundaryVariable(pmy_block, phi, nullptr);
  // to the BoundaryValues::bvars std::vector, it will probably be safe from:
  //
  // 1) dereferencing the nullptr corresponding to the non-existent fluxes when SMR/AMR is
  // used an automatically tries to ProlongateBoundaries()
  // 2) applying the inherited implementations of BoundaryPhysics functions
  // CellCenteredBoundaryVariable::Outflow...*(), etc.
  // 3) generally messing up the variable's own BoundaryData MPI requests and buffers

  // through unnecessary initialization in Mesh::Initialize() and coupling to the
  // time_integrator.cpp main task list's BoundaryValues function calls

  // KGF: DOES FFT SELF-GRAVITY ACTUALLY NEED "enum BoundaryStatus flag[56]" that was
  // copied from Multigrid's unique "struct MGBoundaryData"?

  // Need to trap errors and incompatibilitis with FFT self-gravity when:
  // 1) any boundary codition runtime flag is non-periodic
  // 2) SMR/AMR is used at all (although we could safely use it only for MHD variables and
  // prevent mesh refinement operations from automatically being applied to gravity phi,
  // in the future)

  // KGF: should this be moved to FFTGravity() (currently default ctor) in
  // gravity/fft_gravity.hpp? What exactly is the Gravity class's shared role between FFT
  // and Multigrid??
  return TASK_NEXT;
}
