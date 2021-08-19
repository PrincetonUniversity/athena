//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file radiation_task_list.cpp
//! \brief derived class for radiation integrator task list.

// C/C++ headers
#include <iostream>   // endl
#include <sstream>    // sstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena++ classes headers
#include "../athena.hpp"
#include "../defs.hpp"
#include "../mesh/mesh.hpp"
#include "../radiation/integrators/rad_integrators.hpp"
#include "../radiation/radiation.hpp"
#include "radiation_task_list.hpp"
#include "task_list.hpp"

//--------------------------------------------------------------------------------------
//! RadiationIntegratorTaskList constructor
RadiationIntegratorTaskList::RadiationIntegratorTaskList(ParameterInput *pin, Mesh *pm) {
  integrator = RADIATION_INTEGRATOR;
  // Now assemble list of tasks for each step of chemistry integrator
  {using namespace RadiationIntegratorTaskNames; // NOLINT (build/namespace)
    if (integrator == "loc_jeans") {
      AddTask(INT_LOC_JEANS,NONE);
    } else if (integrator == "const") {
      //do nothing, radiation field constant, remain initial value
      AddTask(INT_CONST,NONE);
    } else {
      std::stringstream msg;
      msg << "### FATAL ERROR in RadiationIntegratorTaskList constructor" << std::endl
        << "integrator=" << integrator << " not valid radiation integrator, "
        << std::endl << "choose from {jeans, const}" << std::endl;
      ATHENA_ERROR(msg);
    }
  } // end of using namespace block
}

//--------------------------------------------------------------------------------------
//! \fn void RadiationIntegratorTaskList::AddTask(const TaskID& id, const TaskID& dep)
//! \brief Sets id and dependency for "ntask" member of task_list_ array, then iterates
//! value of ntask.
void RadiationIntegratorTaskList::AddTask(const TaskID& id, const TaskID& dep) {
  task_list_[ntasks].task_id=id;
  task_list_[ntasks].dependency=dep;

  using namespace RadiationIntegratorTaskNames; // NOLINT (build/namespace)
  if (id == INT_LOC_JEANS) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&RadiationIntegratorTaskList::LocalIntegratorJeans);
  } else if (id == INT_CONST) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&RadiationIntegratorTaskList::ConstRadiation);
  } else {
    std::stringstream msg;
    msg << "### FATAL ERROR in RadiationIntegratorTaskList::AddTask" << std::endl
      << "Invalid Task is specified" << std::endl;
    ATHENA_ERROR(msg);
  }
  ntasks++;
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationIntegratorTaskList::StartupTaskList(MeshBlock *pmb, int stage)
//! \brief Initialize boundary
void RadiationIntegratorTaskList::StartupTaskList(MeshBlock *pmb, int stage) {
  //TODO(Munan Gong): sixray boundary e.g.
  //for fft gravity: pmb->pgrav->gbvar.StartReceiving(BoundaryCommSubset::all);
  return;
}

//----------------------------------------------------------------------------------------
//! Local Jeans integrator
TaskStatus RadiationIntegratorTaskList::LocalIntegratorJeans(MeshBlock *pmb, int stage) {
#ifdef INCLUDE_CHEMISTRY
  pmb->prad->pradintegrator->UpdateRadiation(0);
  pmb->prad->pradintegrator->CopyToOutput();
#endif
  return TaskStatus::success;
}

//----------------------------------------------------------------------------------------
//! Trivial constant integrator
TaskStatus RadiationIntegratorTaskList::ConstRadiation(MeshBlock *pmb, int stage) {
#ifdef INCLUDE_CHEMISTRY
  pmb->prad->pradintegrator->CopyToOutput();
#endif
  return TaskStatus::success;
}

