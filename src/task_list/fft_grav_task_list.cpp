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
  // KGF: BoundaryValues wrapper version called in time_integrator.cpp
  // KGF: what is "time" parameter ever used for in this function?
  // ANSWER: shearing box capabilities, only. Remove this function parameter.
  pmb->pgrav->pgbval->StartReceiving(BoundaryCommSubset::all);
  return;
}

TaskStatus FFTGravitySolverTaskList::ClearFFTGravityBoundary(MeshBlock *pmb, int stage) {
  // KGF: BoundaryValues wrapper version called in time_integrator.cpp
  pmb->pgrav->pgbval->ClearBoundary(BoundaryCommSubset::all);
  return TaskStatus::success;
}

TaskStatus FFTGravitySolverTaskList::SendFFTGravityBoundary(MeshBlock *pmb, int stage) {
  // KGF: BoundaryBuffer version does not return bool (copied Multigrid implementation)
  pmb->pgrav->pgbval->SendBoundaryBuffers();
  return TaskStatus::success;
}

TaskStatus FFTGravitySolverTaskList::ReceiveFFTGravityBoundary(MeshBlock *pmb,
                                                               int stage) {
  bool ret = pmb->pgrav->pgbval->ReceiveBoundaryBuffers();
  //  std::cout << "ret= " << ret << std::endl;
  if (ret == false)
    return TaskStatus::fail;
  // if (pmb->pgrav->pgbval->ReceiveBoundaryBuffers() == false)
  //   return TaskStatus::fail;
  return TaskStatus::success;
}



TaskStatus FFTGravitySolverTaskList::SetFFTGravityBoundary(MeshBlock *pmb, int stage) {
  pmb->pgrav->pgbval->SetBoundaries();
  return TaskStatus::success;
}

TaskStatus FFTGravitySolverTaskList::PhysicalBoundary(MeshBlock *pmb, int stage) {
  //std::cout << "In FFTGravitySolverTaskList::PhysicalBoundary" << std::endl;

  // KGF: FFT self-gravity can only handle periodic boundary conditions

  // KGF: there is no need for any FFT self-gravity boundary values handling in
  // Mesh::Initialize(), main.cpp, nor time_integrator.cpp. All interactions to the
  // CellCenteredBoundaryVariable object owned by Gravity class occurs in this file.

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

  // BUT, BoundaryVariable::CopyVariableBufferSameProcess() REQUIRES that the calling
  // BoundaryVariable object has a matching pointer in the bvars vector in order to locate
  // the matching BoundaryVariable obejct contained in another MeshBlock object on the
  // same process.

  // KGF: Does FFT self-gravity actually need additional "BoundaryStatus sflag[56]"
  // that was copied from Multigrid's unique "struct MGBoundaryData"?

  // Need to trap errors and incompatibilitis with FFT self-gravity when:
  // 1) any boundary codition runtime flag is non-periodic
  // 2) SMR/AMR is used at all (although we could safely use it only for MHD variables and
  // prevent mesh refinement operations from automatically being applied to gravity phi,
  // in the future)

  // KGF: should this be moved to FFTGravity() (currently default ctor) in
  // gravity/fft_gravity.hpp? What exactly is the Gravity class's shared role between FFT
  // and Multigrid??
  return TaskStatus::next;
}
