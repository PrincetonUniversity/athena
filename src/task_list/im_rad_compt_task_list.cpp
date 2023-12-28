//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file im_rad_compt_task_list.cpp
//! \brief function implementation for implicit radiation with Compton

// C headers

// C++ headers
#include <iostream>   // endl
#include <sstream>    // sstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena++ headers
#include "../athena.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../nr_radiation/implicit/radiation_implicit.hpp"
#include "../nr_radiation/integrators/rad_integrators.hpp"
#include "../nr_radiation/radiation.hpp"
#include "./im_rad_task_list.hpp"

//----------------------------------------------------------------------------------------
//! IMRadTaskList constructor

IMRadComptTaskList::IMRadComptTaskList(Mesh *pm) {
  pmy_mesh = pm;
  {using namespace IMRadComptTaskNames; // NOLINT (build/namespace)
    AddTask(CAL_COMPT,NONE);
    AddTask(SEND_RAD_BND,CAL_COMPT);
    AddTask(RECV_RAD_BND,CAL_COMPT);
    AddTask(SETB_RAD_BND,(RECV_RAD_BND|SEND_RAD_BND));
    if (pm->shear_periodic) {
      AddTask(SEND_RAD_SH,SETB_RAD_BND);
      AddTask(RECV_RAD_SH,SEND_RAD_SH|RECV_RAD_BND);
    }
    TaskID setb = SETB_RAD_BND;
    if (pm->shear_periodic)
      setb=(setb|RECV_RAD_SH);
    if (pm->multilevel) {
      AddTask(PRLN_RAD_BND,setb);
      AddTask(RAD_PHYS_BND,PRLN_RAD_BND);
    } else {
      AddTask(RAD_PHYS_BND,setb);
    }
    AddTask(CLEAR_RAD, RAD_PHYS_BND);
  }
}


void IMRadComptTaskList::AddTask(const TaskID& id, const TaskID& dep) {
  task_list_[ntasks].task_id=id;
  task_list_[ntasks].dependency=dep;

  using namespace IMRadComptTaskNames; // NOLINT (build/namespace)
  if (id == CLEAR_RAD) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (IMRadTaskList::*)(MeshBlock*)>
        (&IMRadComptTaskList::ClearRadBoundary);
  } else if (id == SEND_RAD_BND) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (IMRadTaskList::*)(MeshBlock*)>
        (&IMRadComptTaskList::SendRadBoundary);
  } else if (id == RECV_RAD_BND) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (IMRadTaskList::*)(MeshBlock*)>
        (&IMRadComptTaskList::ReceiveRadBoundary);
  } else if (id == SETB_RAD_BND) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (IMRadTaskList::*)(MeshBlock*)>
        (&IMRadComptTaskList::SetRadBoundary);
  } else if (id == RAD_PHYS_BND) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (IMRadTaskList::*)(MeshBlock*)>
        (&IMRadComptTaskList::PhysicalBoundary);
  } else if (id == PRLN_RAD_BND) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (IMRadTaskList::*)(MeshBlock*)>
        (&IMRadComptTaskList::ProlongateBoundary);
  } else if (id == SEND_RAD_SH) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (IMRadTaskList::*)(MeshBlock*)>
        (&IMRadComptTaskList::SendRadBoundaryShear);
  } else if (id == RECV_RAD_SH) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (IMRadTaskList::*)(MeshBlock*)>
        (&IMRadComptTaskList::ReceiveRadBoundaryShear);
  } else if (id == CAL_COMPT) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (IMRadTaskList::*)(MeshBlock*)>
        (&IMRadComptTaskList::CalComptTerms);
  } else {
    std::stringstream msg;
    msg << "### FATAL ERROR in IMRadTaskList::AddTask" << std::endl
        << "Invalid Task is specified" << std::endl;
    ATHENA_ERROR(msg);
  }
  ntasks++;
  return;
}


TaskStatus IMRadComptTaskList::CalComptTerms(MeshBlock *pmb) {
  NRRadiation *prad = pmb->pnrrad;
  Hydro *ph = pmb->phydro;
  prad->pradintegrator->AddMultiGroupCompt(pmb, dt, ph->u, prad->ir);
  return TaskStatus::success;
}

TaskStatus IMRadComptTaskList::ClearRadBoundary(MeshBlock *pmb) {
  pmb->pnrrad->rad_bvar.ClearBoundary(BoundaryCommSubset::radiation);
  return TaskStatus::success;
}

TaskStatus IMRadComptTaskList::SendRadBoundary(MeshBlock *pmb) {
  pmb->pnrrad->rad_bvar.SendBoundaryBuffers();
  return TaskStatus::success;
}

TaskStatus IMRadComptTaskList::SendRadBoundaryShear(MeshBlock *pmb) {
  pmb->pnrrad->rad_bvar.SendShearingBoxBoundaryBuffers();
  return TaskStatus::success;
}

TaskStatus IMRadComptTaskList::ReceiveRadBoundary(MeshBlock *pmb) {
  bool ret = pmb->pnrrad->rad_bvar.ReceiveBoundaryBuffers();
  if (!ret)
    return TaskStatus::fail;
  return TaskStatus::success;
}

TaskStatus IMRadComptTaskList::ReceiveRadBoundaryShear(MeshBlock *pmb) {
  bool ret = pmb->pnrrad->rad_bvar.ReceiveShearingBoxBoundaryBuffers();
  if (ret) {
    pmb->pnrrad->rad_bvar.SetShearingBoxBoundaryBuffers();
    return TaskStatus::success;
  } else {
    return TaskStatus::fail;
  }
}

TaskStatus IMRadComptTaskList::SetRadBoundary(MeshBlock *pmb) {
  pmb->pnrrad->rad_bvar.SetBoundaries();
  return TaskStatus::success;
}

void IMRadComptTaskList::StartupTaskList(MeshBlock *pmb) {
  pmb->pnrrad->rad_bvar.StartReceiving(BoundaryCommSubset::radiation);
  if (pmy_mesh->shear_periodic) {
    pmb->pnrrad->rad_bvar.StartReceivingShear(BoundaryCommSubset::radiation);
  }
  return;
}
