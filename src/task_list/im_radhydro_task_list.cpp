//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file im_radhydro_task_list.cpp
//! \brief function implementation for radiation hydro task list

// C headers

// C++ headers
#include <iostream>   // endl
#include <sstream>    // sstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena++ headers
#include "../athena.hpp"
#include "../bvals/bvals.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../nr_radiation/implicit/radiation_implicit.hpp"
#include "../nr_radiation/integrators/rad_integrators.hpp"
#include "../nr_radiation/radiation.hpp"
#include "../scalars/scalars.hpp"
#include "./im_rad_task_list.hpp"

//----------------------------------------------------------------------------------------
//! IMRadTaskList constructor

IMRadHydroTaskList::IMRadHydroTaskList(Mesh *pm) {
  pmy_mesh = pm;
  // Now assemble list of tasks for each stage of time integrator

  {using namespace IMRadHydroTaskNames; // NOLINT (build/namespace)
    AddTask(ADD_RAD_SRC,NONE);
    AddTask(SEND_HYD_BND,ADD_RAD_SRC);
    AddTask(RECV_HYD_BND,NONE);
    AddTask(SETB_HYD_BND,(RECV_HYD_BND|SEND_HYD_BND));
    if (pm->shear_periodic) {
      AddTask(SEND_HYD_SH,SETB_HYD_BND);
      AddTask(RECV_HYD_SH,SEND_HYD_SH|RECV_HYD_BND);
    }
    TaskID setb = SETB_HYD_BND;
    if (pm->shear_periodic)
      setb=(setb|RECV_HYD_SH);
    if (pm->multilevel) {
      AddTask(PRLN_HYD_BND,setb);
      AddTask(CONS_TO_PRIM,PRLN_HYD_BND);
    } else {
      AddTask(CONS_TO_PRIM,setb);
    }
    AddTask(HYD_PHYS_BND,CONS_TO_PRIM);
    AddTask(CLEAR_HYD, HYD_PHYS_BND);
    AddTask(UPD_OPA,HYD_PHYS_BND);
  } // end of using namespace block
}

//----------------------------------------------------------------------------------------
//! \fn void FFTGravitySolverTaskList::AddTask(const TaskID& id, const TaskID& dep)
//! \brief Sets id and dependency for "ntask" member of task_list_ array, then iterates
//! value of ntask.

void IMRadHydroTaskList::AddTask(const TaskID& id, const TaskID& dep) {
  task_list_[ntasks].task_id=id;
  task_list_[ntasks].dependency=dep;

  using namespace IMRadHydroTaskNames; // NOLINT (build/namespace)
  if (id == CLEAR_HYD) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (IMRadTaskList::*)(MeshBlock*)>
        (&IMRadHydroTaskList::ClearHydroBoundary);
  } else if (id == SEND_HYD_BND) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (IMRadTaskList::*)(MeshBlock*)>
        (&IMRadHydroTaskList::SendHydroBoundary);
  } else if (id == RECV_HYD_BND) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (IMRadTaskList::*)(MeshBlock*)>
        (&IMRadHydroTaskList::ReceiveHydroBoundary);
  } else if (id == SETB_HYD_BND) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (IMRadTaskList::*)(MeshBlock*)>
        (&IMRadHydroTaskList::SetHydroBoundary);
  } else if (id == HYD_PHYS_BND) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (IMRadTaskList::*)(MeshBlock*)>
        (&IMRadTaskList::PhysicalBoundary);
  } else if (id == PRLN_HYD_BND) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (IMRadTaskList::*)(MeshBlock*)>
        (&IMRadTaskList::ProlongateBoundary);
  } else if (id == SEND_HYD_SH) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (IMRadTaskList::*)(MeshBlock*)>
        (&IMRadHydroTaskList::SendHydroBoundaryShear);
  } else if (id == RECV_HYD_SH) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (IMRadTaskList::*)(MeshBlock*)>
        (&IMRadHydroTaskList::ReceiveHydroBoundaryShear);
  } else if (id == UPD_OPA) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (IMRadTaskList::*)(MeshBlock*)>
        (&IMRadHydroTaskList::UpdateOpacity);
  } else if (id == ADD_RAD_SRC) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (IMRadTaskList::*)(MeshBlock*)>
        (&IMRadHydroTaskList::AddRadSource);
  } else if (id == CONS_TO_PRIM) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (IMRadTaskList::*)(MeshBlock*)>
        (&IMRadHydroTaskList::Primitive);
  } else {
    std::stringstream msg;
    msg << "### FATAL ERROR in IMRadTaskList::AddTask" << std::endl
        << "Invalid Task is specified" << std::endl;
    ATHENA_ERROR(msg);
  }
  ntasks++;
  return;
}

TaskStatus IMRadHydroTaskList::ClearHydroBoundary(MeshBlock *pmb) {
  pmb->phydro->hbvar.ClearBoundary(BoundaryCommSubset::radhydro);
  return TaskStatus::success;
}

TaskStatus IMRadHydroTaskList::SendHydroBoundary(MeshBlock *pmb) {
  pmb->phydro->hbvar.SwapHydroQuantity(pmb->phydro->u, HydroBoundaryQuantity::cons);
  pmb->phydro->hbvar.SendBoundaryBuffers();
  return TaskStatus::success;
}

TaskStatus IMRadHydroTaskList::SendHydroBoundaryShear(MeshBlock *pmb) {
  pmb->phydro->hbvar.SendShearingBoxBoundaryBuffers();
  return TaskStatus::success;
}


TaskStatus IMRadHydroTaskList::ReceiveHydroBoundary(MeshBlock *pmb) {
  bool ret = pmb->phydro->hbvar.ReceiveBoundaryBuffers();
  if (!ret)
    return TaskStatus::fail;
  return TaskStatus::success;
}

TaskStatus IMRadHydroTaskList::ReceiveHydroBoundaryShear(MeshBlock *pmb) {
  bool ret = pmb->phydro->hbvar.ReceiveShearingBoxBoundaryBuffers();
  if (ret) {
    pmb->phydro->hbvar.SetShearingBoxBoundaryBuffers();
    return TaskStatus::success;
  } else {
    return TaskStatus::fail;
  }
}


TaskStatus IMRadHydroTaskList::SetHydroBoundary(MeshBlock *pmb) {
  pmb->phydro->hbvar.SwapHydroQuantity(pmb->phydro->u, HydroBoundaryQuantity::cons);
  pmb->phydro->hbvar.SetBoundaries();
  return TaskStatus::success;
}

TaskStatus IMRadHydroTaskList::Primitive(MeshBlock *pmb) {
  Hydro *ph = pmb->phydro;
  Field *pf = pmb->pfield;
  PassiveScalars *ps = pmb->pscalars;
  BoundaryValues *pbval = pmb->pbval;
  int il = pmb->is, iu = pmb->ie, jl = pmb->js, ju = pmb->je, kl = pmb->ks, ku = pmb->ke;
  if (pbval->nblevel[1][1][0] != -1) il -= NGHOST;
  if (pbval->nblevel[1][1][2] != -1) iu += NGHOST;
  if (pbval->nblevel[1][0][1] != -1) jl -= NGHOST;
  if (pbval->nblevel[1][2][1] != -1) ju += NGHOST;
  if (pbval->nblevel[0][1][1] != -1) kl -= NGHOST;
  if (pbval->nblevel[2][1][1] != -1) ku += NGHOST;
  pmb->peos->ConservedToPrimitive(ph->u, ph->w, pf->b,
                                  ph->w1, pf->bcc, pmb->pcoord,
                                  il, iu, jl, ju, kl, ku);
  if (NSCALARS > 0) {
    // r1/r_old for GR is currently unused:
    pmb->peos->PassiveScalarConservedToPrimitive(ps->s, ph->u, ps->r, ps->r,
                                                 pmb->pcoord, il, iu, jl, ju, kl, ku);
  }
  ph->w.SwapAthenaArray(ph->w1);

  return TaskStatus::success;
}

TaskStatus IMRadHydroTaskList::UpdateOpacity(MeshBlock *pmb) {
  pmb->pnrrad->UpdateOpacity(pmb, pmb->phydro->w);
  return TaskStatus::success;
}

TaskStatus IMRadHydroTaskList::AddRadSource(MeshBlock *pmb) {
  if (pmb->pnrrad->set_source_flag > 0) {
    pmb->pnrrad->pradintegrator->AddSourceTerms(pmb, pmb->phydro->u);
  }
  return TaskStatus::success;
}

void IMRadHydroTaskList::StartupTaskList(MeshBlock *pmb) {
  pmb->phydro->hbvar.StartReceiving(BoundaryCommSubset::radhydro);
  if (pmy_mesh->shear_periodic) {
    pmb->phydro->hbvar.StartReceivingShear(BoundaryCommSubset::radhydro);
  }
  return;
}
