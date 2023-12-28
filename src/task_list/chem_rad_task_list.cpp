//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file chem_rad_task_list.cpp
//! \brief derived class for radiation integrator task list.

// C/C++ headers
#include <iostream>   // endl
#include <sstream>    // sstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena++ classes headers
#include "../athena.hpp"
#include "../chem_rad/chem_rad.hpp"
#include "../chem_rad/integrators/rad_integrators.hpp"
#include "../defs.hpp"
#include "../mesh/mesh.hpp"
#include "chem_rad_task_list.hpp"
#include "task_list.hpp"

//--------------------------------------------------------------------------------------
//! ChemRadiationIntegratorTaskList constructor
ChemRadiationIntegratorTaskList::ChemRadiationIntegratorTaskList(ParameterInput *pin,
                                                                 Mesh *pm) {
  integrator = CHEMRADIATION_INTEGRATOR;
  // Now assemble list of tasks for each step of chemistry integrator
  {using namespace ChemRadiationIntegratorTaskNames; // NOLINT (build/namespace)
    if (integrator == "six_ray") {
      AddTask(GET_COL_MB_IX1,NONE);
      AddTask(RECV_SEND_COL_IX1,GET_COL_MB_IX1);
      AddTask(GET_COL_MB_OX1,NONE);
      AddTask(RECV_SEND_COL_OX1,GET_COL_MB_OX1);
      AddTask(GET_COL_MB_IX2,NONE);
      AddTask(RECV_SEND_COL_IX2,GET_COL_MB_IX2);
      AddTask(GET_COL_MB_OX2,NONE);
      AddTask(RECV_SEND_COL_OX2,GET_COL_MB_OX2);
      AddTask(GET_COL_MB_IX3,NONE);
      AddTask(RECV_SEND_COL_IX3,GET_COL_MB_IX3);
      AddTask(GET_COL_MB_OX3,NONE);
      AddTask(RECV_SEND_COL_OX3,GET_COL_MB_OX3);
      AddTask(UPDATE_RAD,
          RECV_SEND_COL_IX1|RECV_SEND_COL_OX1|RECV_SEND_COL_IX2|
          RECV_SEND_COL_OX2|RECV_SEND_COL_IX3|RECV_SEND_COL_OX3);
      AddTask(CLEAR_SIXRAY_RECV,UPDATE_RAD);
    } else if (integrator == "const") {
      //do nothing, radiation field constant, remain initial value
      AddTask(INT_CONST,NONE);
    } else {
      std::stringstream msg;
      msg << "### FATAL ERROR in ChemRadiationIntegratorTaskList constructor" << std::endl
        << "integrator=" << integrator << " not valid radiation integrator, "
        << std::endl << "choose from {six_ray, const}" << std::endl;
      ATHENA_ERROR(msg);
    }
  } // end of using namespace block
}

//--------------------------------------------------------------------------------------
//! \fn void ChemRadiationIntegratorTaskList::AddTask(const TaskID& id, const TaskID& dep)
//! \brief Sets id and dependency for "ntask" member of task_list_ array, then iterates
//! value of ntask.
void ChemRadiationIntegratorTaskList::AddTask(const TaskID& id, const TaskID& dep) {
  task_list_[ntasks].task_id=id;
  task_list_[ntasks].dependency=dep;

  using namespace ChemRadiationIntegratorTaskNames; // NOLINT (build/namespace)
  if (id == UPDATE_RAD) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&ChemRadiationIntegratorTaskList::UpdateRadiationSixRay);
  } else if (id == GET_COL_MB_IX1) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&ChemRadiationIntegratorTaskList::GetColMB_ix1);
  } else if (id == GET_COL_MB_OX1) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&ChemRadiationIntegratorTaskList::GetColMB_ox1);
  } else if (id == GET_COL_MB_IX2) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&ChemRadiationIntegratorTaskList::GetColMB_ix2);
  } else if (id == GET_COL_MB_OX2) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&ChemRadiationIntegratorTaskList::GetColMB_ox2);
  } else if (id == GET_COL_MB_IX3) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&ChemRadiationIntegratorTaskList::GetColMB_ix3);
  } else if (id == GET_COL_MB_OX3) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&ChemRadiationIntegratorTaskList::GetColMB_ox3);
  } else if (id == RECV_SEND_COL_IX1) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&ChemRadiationIntegratorTaskList::RecvAndSend_ix1);
  } else if (id == RECV_SEND_COL_OX1) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&ChemRadiationIntegratorTaskList::RecvAndSend_ox1);
  } else if (id == RECV_SEND_COL_IX2) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&ChemRadiationIntegratorTaskList::RecvAndSend_ix2);
  } else if (id == RECV_SEND_COL_OX2) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&ChemRadiationIntegratorTaskList::RecvAndSend_ox2);
  } else if (id == RECV_SEND_COL_IX3) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&ChemRadiationIntegratorTaskList::RecvAndSend_ix3);
  } else if (id == RECV_SEND_COL_OX3) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&ChemRadiationIntegratorTaskList::RecvAndSend_ox3);
  } else if (id == CLEAR_SIXRAY_RECV) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&ChemRadiationIntegratorTaskList::ClearSixrayReceive);
  } else if (id == INT_CONST) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&ChemRadiationIntegratorTaskList::ConstRadiation);
  } else {
    std::stringstream msg;
    msg << "### FATAL ERROR in ChemRadiationIntegratorTaskList::AddTask" << std::endl
      << "Invalid Task is specified" << std::endl;
    ATHENA_ERROR(msg);
  }
  ntasks++;
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void ChemRadiationIntegratorTaskList::StartupTaskList(MeshBlock *pmb, int stage)
//! \brief Initialize boundary
void ChemRadiationIntegratorTaskList::StartupTaskList(MeshBlock *pmb, int stage) {
  if (CHEMISTRY_ENABLED && integrator == "six_ray") {
    pmb->pchemrad->pchemradintegrator->col_bvar.StartReceiving(BoundaryCommSubset::all);
  }
  return;
}

//----------------------------------------------------------------------------------------
//! Six-ray

//meshblock column densities
TaskStatus ChemRadiationIntegratorTaskList::GetColMB_ix1(MeshBlock *pmb, int step) {
  if (CHEMISTRY_ENABLED) {
    pmb->pchemrad->pchemradintegrator->GetColMB(BoundaryFace::inner_x1);
  }
  return TaskStatus::success;
}

TaskStatus ChemRadiationIntegratorTaskList::GetColMB_ox1(MeshBlock *pmb, int step) {
  if (CHEMISTRY_ENABLED) {
    pmb->pchemrad->pchemradintegrator->GetColMB(BoundaryFace::outer_x1);
  }
  return TaskStatus::success;
}

TaskStatus ChemRadiationIntegratorTaskList::GetColMB_ix2(MeshBlock *pmb, int step) {
  if (CHEMISTRY_ENABLED) {
    pmb->pchemrad->pchemradintegrator->GetColMB(BoundaryFace::inner_x2);
  }
  return TaskStatus::success;
}

TaskStatus ChemRadiationIntegratorTaskList::GetColMB_ox2(MeshBlock *pmb, int step) {
  if (CHEMISTRY_ENABLED) {
    pmb->pchemrad->pchemradintegrator->GetColMB(BoundaryFace::outer_x2);
  }
  return TaskStatus::success;
}

TaskStatus ChemRadiationIntegratorTaskList::GetColMB_ix3(MeshBlock *pmb, int step) {
  if (CHEMISTRY_ENABLED) {
    pmb->pchemrad->pchemradintegrator->GetColMB(BoundaryFace::inner_x3);
  }
  return TaskStatus::success;
}

TaskStatus ChemRadiationIntegratorTaskList::GetColMB_ox3(MeshBlock *pmb, int step) {
  if (CHEMISTRY_ENABLED) {
    pmb->pchemrad->pchemradintegrator->GetColMB(BoundaryFace::outer_x3);
  }
  return TaskStatus::success;
}

//boundary receive and send

TaskStatus ChemRadiationIntegratorTaskList::RecvAndSend_direction(MeshBlock *pmb,
    int step, BoundaryFace direction) {
  if (CHEMISTRY_ENABLED) {
    SixRayBoundaryVariable *pbvar = &pmb->pchemrad->pchemradintegrator->col_bvar;
    BoundaryFace direction_opp = pbvar->GetOppositeBoundaryFace(direction);
    NeighborBlock *pnb = pbvar->GetFaceNeighbor(direction);
    NeighborBlock *pnb_opp = pbvar->GetFaceNeighbor(direction_opp);
    bool ret = true;
    if (pnb == nullptr) {
      if (pnb_opp != nullptr) {
        pbvar->SendSixRayBoundaryBuffers(direction_opp);
      }
    } else {
      ret = pbvar->ReceiveAndSetSixRayBoundaryBuffers(direction);
      if (ret == true) {
        pmb->pchemrad->pchemradintegrator->UpdateCol(direction);
        if (pnb_opp != nullptr) {
          pbvar->SendSixRayBoundaryBuffers(direction_opp);
        }
      } else {
        return TaskStatus::fail;
      }
    }
  }
  return TaskStatus::success;
}

TaskStatus ChemRadiationIntegratorTaskList::RecvAndSend_ix1(MeshBlock *pmb, int step) {
  return RecvAndSend_direction(pmb, step, BoundaryFace::inner_x1);
}

TaskStatus ChemRadiationIntegratorTaskList::RecvAndSend_ox1(MeshBlock *pmb, int step) {
  return RecvAndSend_direction(pmb, step, BoundaryFace::outer_x1);
}

TaskStatus ChemRadiationIntegratorTaskList::RecvAndSend_ix2(MeshBlock *pmb, int step) {
  return RecvAndSend_direction(pmb, step, BoundaryFace::inner_x2);
}

TaskStatus ChemRadiationIntegratorTaskList::RecvAndSend_ox2(MeshBlock *pmb, int step) {
  return RecvAndSend_direction(pmb, step, BoundaryFace::outer_x2);
}

TaskStatus ChemRadiationIntegratorTaskList::RecvAndSend_ix3(MeshBlock *pmb, int step) {
  return RecvAndSend_direction(pmb, step, BoundaryFace::inner_x3);
}

TaskStatus ChemRadiationIntegratorTaskList::RecvAndSend_ox3(MeshBlock *pmb, int step) {
  return RecvAndSend_direction(pmb, step, BoundaryFace::outer_x3);
}

TaskStatus ChemRadiationIntegratorTaskList::ClearSixrayReceive(MeshBlock *pmb,
                                                               int step) {
  if (CHEMISTRY_ENABLED) {
    pmb->pchemrad->pchemradintegrator->col_bvar.ClearBoundary(BoundaryCommSubset::all);
  }
  return TaskStatus::success;
}

//update radiation variables
TaskStatus ChemRadiationIntegratorTaskList::UpdateRadiationSixRay(MeshBlock *pmb,
                                                                  int step) {
  if (CHEMISTRY_ENABLED) {
    pmb->pchemrad->pchemradintegrator->UpdateRadiation();
    pmb->pchemrad->pchemradintegrator->CopyToOutput();
  }
  return TaskStatus::success;
}


//----------------------------------------------------------------------------------------
//! Trivial constant integrator
TaskStatus ChemRadiationIntegratorTaskList::ConstRadiation(MeshBlock *pmb, int stage) {
  if (CHEMISTRY_ENABLED) {
    pmb->pchemrad->pchemradintegrator->CopyToOutput();
  }
  return TaskStatus::success;
}
