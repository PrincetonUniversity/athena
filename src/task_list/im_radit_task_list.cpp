//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file im_radit_task_list.cpp
//! \brief function implementation for iteration in implicit radiation

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

IMRadITTaskList::IMRadITTaskList(Mesh *pm) {
  pmy_mesh = pm;
  {using namespace IMRadITTaskNames; // NOLINT (build/namespace)
    AddTask(FLX_AND_SRC,NONE);
    AddTask(SEND_RAD_BND,FLX_AND_SRC);
    AddTask(RECV_RAD_BND,FLX_AND_SRC);
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
    // check residual does not need ghost zones
    AddTask(CHK_RAD_RES,FLX_AND_SRC);
  } // end of using namespace block
}


void IMRadITTaskList::AddTask(const TaskID& id, const TaskID& dep) {
  task_list_[ntasks].task_id=id;
  task_list_[ntasks].dependency=dep;

  using namespace IMRadITTaskNames; // NOLINT (build/namespace)
  if (id == CLEAR_RAD) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (IMRadTaskList::*)(MeshBlock*)>
        (&IMRadITTaskList::ClearRadBoundary);
  } else if (id == SEND_RAD_BND) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (IMRadTaskList::*)(MeshBlock*)>
        (&IMRadITTaskList::SendRadBoundary);
  } else if (id == RECV_RAD_BND) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (IMRadTaskList::*)(MeshBlock*)>
        (&IMRadITTaskList::ReceiveRadBoundary);
  } else if (id == SETB_RAD_BND) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (IMRadTaskList::*)(MeshBlock*)>
        (&IMRadITTaskList::SetRadBoundary);
  } else if (id == RAD_PHYS_BND) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (IMRadTaskList::*)(MeshBlock*)>
        (&IMRadITTaskList::PhysicalBoundary);
  } else if (id == PRLN_RAD_BND) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (IMRadTaskList::*)(MeshBlock*)>
        (&IMRadITTaskList::ProlongateBoundary);
  } else if (id == SEND_RAD_SH) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (IMRadTaskList::*)(MeshBlock*)>
        (&IMRadITTaskList::SendRadBoundaryShear);
  } else if (id == RECV_RAD_SH) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (IMRadTaskList::*)(MeshBlock*)>
        (&IMRadITTaskList::ReceiveRadBoundaryShear);
  } else if (id == CHK_RAD_RES) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (IMRadTaskList::*)(MeshBlock*)>
        (&IMRadITTaskList::CheckResidual);
  } else if (id == FLX_AND_SRC) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (IMRadTaskList::*)(MeshBlock*)>
        (&IMRadITTaskList::AddFluxAndSourceTerms);
  } else {
    std::stringstream msg;
    msg << "### FATAL ERROR in IMRadTaskList::AddTask" << std::endl
        << "Invalid Task is specified" << std::endl;
    ATHENA_ERROR(msg);
  }
  ntasks++;
  return;
}

TaskStatus IMRadITTaskList::AddFluxAndSourceTerms(MeshBlock *pmb) {
  NRRadiation *prad = pmb->pnrrad;
  Hydro *ph = pmb->phydro;
  const int &rb_or_not = pmy_mesh->pimrad->rb_or_not;

  int is=pmb->is, ie=pmb->ie;
  int js=pmb->js, je=pmb->je;
  int ks=pmb->ks, ke=pmb->ke;

  if (rb_or_not == 0) {
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie; ++i) {
          prad->pradintegrator->FirstOrderFluxDivergence(k, j, i,
                                                         prad->ir_old);
          // add angular flux
          if (prad->angle_flag == 1) {
            prad->pradintegrator->ImplicitAngularFluxes(k,j,i,prad->ir);
          }
          // add the source terms together
          prad->pradintegrator->CalSourceTerms(pmb, dt, k,j,i, ph->u,
                                               prad->ir1, prad->ir);
        }
      }
    }
  } else if (rb_or_not == 1) {
    // first the red points
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie; ++i) {
          if (((i-is)+(j-js)+(k-ks))%2 == 0) {
            prad->pradintegrator->FirstOrderFluxDivergence(k, j, i,
                                                           prad->ir);
            // add angular flux
            if (prad->angle_flag == 1) {
              prad->pradintegrator->ImplicitAngularFluxes(k,j,i,prad->ir);
            }
            // add the source terms together
            prad->pradintegrator->CalSourceTerms(pmb, dt, k,j,i, ph->u,
                                                 prad->ir1, prad->ir);
          }
        }
      }
    }
  } else if (rb_or_not == 2) {
    // now the black points
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie; ++i) {
          if (((i-is)+(j-js)+(k-ks))%2 == 1) {
            prad->pradintegrator->FirstOrderFluxDivergence(k, j, i,
                                                           prad->ir);

            // add angular flux
            if (prad->angle_flag == 1) {
              prad->pradintegrator->ImplicitAngularFluxes(k,j,i,prad->ir);
            }

            // add the source terms together
            prad->pradintegrator->CalSourceTerms(pmb, dt, k,j,i, ph->u,
                                                 prad->ir1, prad->ir);
          }
        }
      }
    }
  } else {
    std::stringstream msg;
    msg << "### FATAL ERROR in function [Iteration]"
        << std::endl << "red_or_black '" << rb_or_not << "' not allowed!";
    ATHENA_ERROR(msg);
  }
  return TaskStatus::success;
}

TaskStatus IMRadITTaskList::ClearRadBoundary(MeshBlock *pmb) {
  pmb->pnrrad->rad_bvar.ClearBoundary(BoundaryCommSubset::radiation);
  return TaskStatus::success;
}

TaskStatus IMRadITTaskList::SendRadBoundary(MeshBlock *pmb) {
  pmb->pnrrad->rad_bvar.SendBoundaryBuffers();
  return TaskStatus::success;
}

TaskStatus IMRadITTaskList::SendRadBoundaryShear(MeshBlock *pmb) {
  pmb->pnrrad->rad_bvar.SendShearingBoxBoundaryBuffers();
  return TaskStatus::success;
}

TaskStatus IMRadITTaskList::ReceiveRadBoundary(MeshBlock *pmb) {
  bool ret = pmb->pnrrad->rad_bvar.ReceiveBoundaryBuffers();
  if (!ret) {
    return TaskStatus::fail;
  }
  return TaskStatus::success;
}

TaskStatus IMRadITTaskList::ReceiveRadBoundaryShear(MeshBlock *pmb) {
  bool ret = pmb->pnrrad->rad_bvar.ReceiveShearingBoxBoundaryBuffers();
  if (ret) {
    pmb->pnrrad->rad_bvar.SetShearingBoxBoundaryBuffers();
    return TaskStatus::success;
  } else {
    return TaskStatus::fail;
  }
}

TaskStatus IMRadITTaskList::SetRadBoundary(MeshBlock *pmb) {
  pmb->pnrrad->rad_bvar.SetBoundaries();
  return TaskStatus::success;
}

TaskStatus IMRadITTaskList::CheckResidual(MeshBlock *pmb) {
  pmy_mesh->pimrad->CheckResidual(pmb, pmb->pnrrad->ir_old,pmb->pnrrad->ir);
  return TaskStatus::success;
}

void IMRadITTaskList::StartupTaskList(MeshBlock *pmb) {
  pmb->pnrrad->rad_bvar.StartReceiving(BoundaryCommSubset::radiation);
  if (pmy_mesh->shear_periodic) {
    pmb->pnrrad->rad_bvar.StartReceivingShear(BoundaryCommSubset::radiation);
  }
  return;
}
