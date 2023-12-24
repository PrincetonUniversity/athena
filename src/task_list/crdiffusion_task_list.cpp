//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file grav_task_list.cpp
//! \brief function implementation for CRDiffusionBoundaryTaskList

// C headers

// C++ headers
#include <iostream>   // endl
#include <sstream>    // sstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena++ headers
#include "../athena.hpp"
#include "../crdiffusion/crdiffusion.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "crdiffusion_task_list.hpp"
#include "task_list.hpp"

//----------------------------------------------------------------------------------------
//! CRDiffusionBoundaryTaskList constructor

CRDiffusionBoundaryTaskList::CRDiffusionBoundaryTaskList(ParameterInput *pin, Mesh *pm) {
  // Now assemble list of tasks for each stage of time integrator
  {using namespace CRDiffusionBoundaryTaskNames; // NOLINT (build/namespace)
    AddTask(SEND_CRDIFF_BND,NONE);
    AddTask(RECV_CRDIFF_BND,NONE);
    AddTask(SETB_CRDIFF_BND,(RECV_CRDIFF_BND|SEND_CRDIFF_BND));
    if (pm->multilevel) {
      AddTask(PROLONG_CRDIFF_BND,SETB_CRDIFF_BND);
      AddTask(CRDIFF_PHYS_BND,PROLONG_CRDIFF_BND);
    } else {
      AddTask(CRDIFF_PHYS_BND,SETB_CRDIFF_BND);
    }
    AddTask(CLEAR_CRDIFF, CRDIFF_PHYS_BND);
  } // end of using namespace block
}

//----------------------------------------------------------------------------------------
//! \fn void CRDiffusionBoundaryTaskList::AddTask(const TaskID& id, const TaskID& dep)
//! \brief Sets id and dependency for "ntask" member of task_list_ array, then iterates
//! value of ntask.

void CRDiffusionBoundaryTaskList::AddTask(const TaskID& id, const TaskID& dep) {
  task_list_[ntasks].task_id=id;
  task_list_[ntasks].dependency=dep;

  using namespace CRDiffusionBoundaryTaskNames; // NOLINT (build/namespace)
  if (id == CLEAR_CRDIFF) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&CRDiffusionBoundaryTaskList::ClearCRDiffusionBoundary);
  } else if (id == SEND_CRDIFF_BND) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&CRDiffusionBoundaryTaskList::SendCRDiffusionBoundary);
  } else if (id == RECV_CRDIFF_BND) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&CRDiffusionBoundaryTaskList::ReceiveCRDiffusionBoundary);
  } else if (id == SETB_CRDIFF_BND) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&CRDiffusionBoundaryTaskList::SetCRDiffusionBoundary);
  } else if (id == PROLONG_CRDIFF_BND) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&CRDiffusionBoundaryTaskList::ProlongateCRDiffusionBoundary);
  } else if (id == CRDIFF_PHYS_BND) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&CRDiffusionBoundaryTaskList::PhysicalBoundary);
  } else {
    std::stringstream msg;
    msg << "### FATAL ERROR in CRDiffusionBoundaryTaskList::AddTask" << std::endl
        << "Invalid Task is specified" << std::endl;
    ATHENA_ERROR(msg);
  }
  ntasks++;
  return;
}

void CRDiffusionBoundaryTaskList::StartupTaskList(MeshBlock *pmb, int stage) {
  pmb->pcrdiff->crbvar.StartReceiving(BoundaryCommSubset::all);
  return;
}

TaskStatus CRDiffusionBoundaryTaskList::ClearCRDiffusionBoundary(MeshBlock *pmb,
                                                                 int stage) {
  pmb->pcrdiff->crbvar.ClearBoundary(BoundaryCommSubset::all);
  return TaskStatus::success;
}

TaskStatus CRDiffusionBoundaryTaskList::SendCRDiffusionBoundary(MeshBlock *pmb,
                                                                int stage) {
  pmb->pcrdiff->crbvar.SendBoundaryBuffers();
  return TaskStatus::success;
}

TaskStatus CRDiffusionBoundaryTaskList::ReceiveCRDiffusionBoundary(MeshBlock *pmb,
                                                                   int stage) {
  bool ret = pmb->pcrdiff->crbvar.ReceiveBoundaryBuffers();
  if (!ret)
    return TaskStatus::fail;
  return TaskStatus::success;
}

TaskStatus CRDiffusionBoundaryTaskList::SetCRDiffusionBoundary(MeshBlock *pmb,
                                                               int stage) {
  pmb->pcrdiff->crbvar.SetBoundaries();
  return TaskStatus::success;
}

TaskStatus CRDiffusionBoundaryTaskList::ProlongateCRDiffusionBoundary(MeshBlock *pmb,
                                                                      int stage) {
  pmb->pbval->ProlongateBoundariesPostMG(&(pmb->pcrdiff->crbvar));
  return TaskStatus::success;
}

TaskStatus CRDiffusionBoundaryTaskList::PhysicalBoundary(MeshBlock *pmb, int stage) {
  pmb->pcrdiff->crbvar.ExpandPhysicalBoundaries();
  return TaskStatus::next;
}
