//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file grav_task_list.cpp
//! \brief function implementation for GravityBoundaryTaskList

// C headers

// C++ headers
#include <iostream>   // endl
#include <sstream>    // sstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena++ headers
#include "../athena.hpp"
#include "../gravity/gravity.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "grav_task_list.hpp"
#include "task_list.hpp"

//----------------------------------------------------------------------------------------
//! GravityBoundaryTaskList constructor

GravityBoundaryTaskList::GravityBoundaryTaskList(ParameterInput *pin, Mesh *pm) {
  {using namespace GravityBoundaryTaskNames; // NOLINT (build/namespace)
    AddTask(SEND_GRAV_BND,NONE);
    AddTask(RECV_GRAV_BND,NONE);
    AddTask(SETB_GRAV_BND,(RECV_GRAV_BND|SEND_GRAV_BND));
    if (pm->multilevel) {
      AddTask(PROLONG_GRAV_BND,SETB_GRAV_BND);
      AddTask(GRAV_PHYS_BND,PROLONG_GRAV_BND);
    } else {
      AddTask(GRAV_PHYS_BND,SETB_GRAV_BND);
    }
    AddTask(CLEAR_GRAV, GRAV_PHYS_BND);
  } // end of using namespace block
}

//----------------------------------------------------------------------------------------
//! \fn void GravityBoundaryTaskList::AddTask(const TaskID& id, const TaskID& dep)
//! \brief Sets id and dependency for "ntask" member of task_list_ array, then iterates
//! value of ntask.

void GravityBoundaryTaskList::AddTask(const TaskID& id, const TaskID& dep) {
  task_list_[ntasks].task_id=id;
  task_list_[ntasks].dependency=dep;

  using namespace GravityBoundaryTaskNames; // NOLINT (build/namespace)
  if (id == CLEAR_GRAV) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&GravityBoundaryTaskList::ClearGravityBoundary);
  } else if (id == SEND_GRAV_BND) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&GravityBoundaryTaskList::SendGravityBoundary);
  } else if (id == RECV_GRAV_BND) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&GravityBoundaryTaskList::ReceiveGravityBoundary);
  } else if (id == SETB_GRAV_BND) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&GravityBoundaryTaskList::SetGravityBoundary);
  } else if (id == PROLONG_GRAV_BND) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&GravityBoundaryTaskList::ProlongateGravityBoundary);
  } else if (id == GRAV_PHYS_BND) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&GravityBoundaryTaskList::PhysicalBoundary);
  } else {
    std::stringstream msg;
    msg << "### FATAL ERROR in GravityBoundaryTaskList::AddTask" << std::endl
        << "Invalid Task is specified" << std::endl;
    ATHENA_ERROR(msg);
  }
  ntasks++;
  return;
}

void GravityBoundaryTaskList::StartupTaskList(MeshBlock *pmb, int stage) {
  pmb->pgrav->gbvar.StartReceiving(BoundaryCommSubset::all);
  return;
}

TaskStatus GravityBoundaryTaskList::ClearGravityBoundary(MeshBlock *pmb, int stage) {
  pmb->pgrav->gbvar.ClearBoundary(BoundaryCommSubset::all);
  return TaskStatus::success;
}

TaskStatus GravityBoundaryTaskList::SendGravityBoundary(MeshBlock *pmb, int stage) {
  if (pmb->pgrav->fill_ghost)
    pmb->pgrav->SaveFaceBoundaries();
  pmb->pgrav->gbvar.SendBoundaryBuffers();
  return TaskStatus::success;
}

TaskStatus GravityBoundaryTaskList::ReceiveGravityBoundary(MeshBlock *pmb,
                                                               int stage) {
  bool ret = pmb->pgrav->gbvar.ReceiveBoundaryBuffers();
  if (!ret)
    return TaskStatus::fail;
  return TaskStatus::success;
}

TaskStatus GravityBoundaryTaskList::SetGravityBoundary(MeshBlock *pmb, int stage) {
  pmb->pgrav->gbvar.SetBoundaries();
  return TaskStatus::success;
}

TaskStatus GravityBoundaryTaskList::ProlongateGravityBoundary(MeshBlock *pmb,
                                                              int stage) {
  pmb->pbval->ProlongateBoundariesPostMG(&(pmb->pgrav->gbvar));
  return TaskStatus::success;
}

TaskStatus GravityBoundaryTaskList::PhysicalBoundary(MeshBlock *pmb, int stage) {
  if (pmb->pgrav->fill_ghost) {
    pmb->pgrav->RestoreFaceBoundaries();
    pmb->pgrav->gbvar.ExpandPhysicalBoundaries();
  }
  return TaskStatus::next;
}
