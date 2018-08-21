//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file grav_task_list.cpp
//  \brief

// C/C++ headers
#include <iostream>   // endl
#include <sstream>    // sstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena++ classes headers
#include "task_list.hpp"
#include "grav_task_list.hpp"
#include "../athena.hpp"
#include "../parameter_input.hpp"
#include "../mesh/mesh.hpp"
#include "../hydro/hydro.hpp"
#include "../field/field.hpp"
#include "../gravity/gravity.hpp"
#include "../eos/eos.hpp"
#include "../hydro/srcterms/hydro_srcterms.hpp"
#include "../bvals/bvals_grav.hpp"

//----------------------------------------------------------------------------------------
//  GravitySolverTaskList constructor

GravitySolverTaskList::GravitySolverTaskList(ParameterInput *pin, Mesh *pm)
  : TaskList(pm) {

  // Now assemble list of tasks for each stage of time integrator
  {using namespace GravitySolverTaskNames; // NOLINT (build/namespace)
    AddGravitySolverTask(START_GRAV_RECV,NONE);

    // compute hydro fluxes, integrate hydro variables
    AddGravitySolverTask(SEND_GRAV_BND,START_GRAV_RECV);
    AddGravitySolverTask(RECV_GRAV_BND,START_GRAV_RECV);
    AddGravitySolverTask(GRAV_PHYS_BND,SEND_GRAV_BND|RECV_GRAV_BND);
    AddGravitySolverTask(CLEAR_GRAV, GRAV_PHYS_BND);
  } // end of using namespace block
}

//----------------------------------------------------------------------------------------
//! \fn void GravitySolverTaskList::AddGravitySolverTask(uint64_t id, uint64_t dep)
//  \brief Sets id and dependency for "ntask" member of task_list_ array, then iterates
//  value of ntask.

void GravitySolverTaskList::AddGravitySolverTask(uint64_t id, uint64_t dep) {
  task_list_[ntasks].task_id=id;
  task_list_[ntasks].dependency=dep;

  using namespace GravitySolverTaskNames; // NOLINT (build/namespace)
  switch((id)) {
    case (START_GRAV_RECV):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&GravitySolverTaskList::StartGravityReceive);
      break;
    case (CLEAR_GRAV):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&GravitySolverTaskList::ClearGravityBoundary);
      break;
    case (SEND_GRAV_BND):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&GravitySolverTaskList::SendGravityBoundary);
      break;
    case (RECV_GRAV_BND):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&GravitySolverTaskList::ReceiveGravityBoundary);
      break;
    case (GRAV_PHYS_BND):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&GravitySolverTaskList::PhysicalBoundary);
      break;
    default:
      std::stringstream msg;
      msg << "### FATAL ERROR in AddGravitySolverTask" << std::endl
          << "Invalid Task "<< id << " is specified" << std::endl;
      throw std::runtime_error(msg.str().c_str());
  }
  ntasks++;
  return;
}

//----------------------------------------------------------------------------------------
//! \fn
//  \brief

//----------------------------------------------------------------------------------------
// Functions to start/end MPI communication

enum TaskStatus GravitySolverTaskList::StartGravityReceive(MeshBlock *pmb, int stage) {
  pmb->pgbval->StartReceivingGravity();
  return TASK_SUCCESS;
}

enum TaskStatus GravitySolverTaskList::ClearGravityBoundary(MeshBlock *pmb, int stage) {
  pmb->pgbval->ClearBoundaryGravity();
  return TASK_SUCCESS;
}

enum TaskStatus GravitySolverTaskList::SendGravityBoundary(MeshBlock *pmb, int stage) {
  if (pmb->pgbval->SendGravityBoundaryBuffers(pmb->pgrav->phi)==false)
    return TASK_FAIL;
  return TASK_SUCCESS;
}

enum TaskStatus GravitySolverTaskList::ReceiveGravityBoundary(MeshBlock *pmb, int stage) {
  if (pmb->pgbval->ReceiveGravityBoundaryBuffers(pmb->pgrav->phi)==false)
    return TASK_FAIL;
  return TASK_SUCCESS;
}

enum TaskStatus GravitySolverTaskList::PhysicalBoundary(MeshBlock *pmb, int stage) {
  pmb->pgbval->ApplyPhysicalBoundaries();
  return TASK_NEXT;
}
