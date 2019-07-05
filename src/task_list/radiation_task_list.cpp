//======================================================================================
// Athena++ astrophysical MHD code
// Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
//
// This program is free software: you can redistribute and/or modify it under the terms
// of the GNU General Public License (GPL) as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
// PARTICULAR PURPOSE.  See the GNU General Public License for more details.
//
// You should have received a copy of GNU GPL in the file LICENSE included in the code
// distribution.  If not see <http://www.gnu.org/licenses/>.
//======================================================================================
//! \file radiation_task_list.cpp
//  \brief derived class for radiation integrator task list.
//======================================================================================

// C/C++ headers
#include <iostream>   // endl
#include <sstream>    // sstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena++ classes headers
#include "../athena.hpp"
#include "../mesh/mesh.hpp"
#include "../defs.hpp"
#include "../radiation/integrators/rad_integrators.hpp"
#include "../radiation/radiation.hpp"
#include "radiation_task_list.hpp"
#include "task_list.hpp"

//--------------------------------------------------------------------------------------
//  RadiationIntegratorTaskList constructor
RadiationIntegratorTaskList::RadiationIntegratorTaskList(ParameterInput *pin, Mesh *pm)
{
  integrator = RADIATION_INTEGRATOR;
  // Now assemble list of tasks for each step of chemistry integrator
  {using namespace RadiationIntegratorTaskNames;
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
//! \fn
//  \brief Sets id and dependency for "ntask" member of task_list_ array, then iterates
//  value of ntask.  
void RadiationIntegratorTaskList::AddTask(const TaskID& id, const TaskID& dep)
{
  task_list_[ntasks].task_id=id;
  task_list_[ntasks].dependency=dep;

  using namespace RadiationIntegratorTaskNames;
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
void RadiationIntegratorTaskList::StartupTaskList(MeshBlock *pmb, int stage) {
  //TODO: sixray boundary e.g. 
  //for fft gravity: pmb->pgrav->gbvar.StartReceiving(BoundaryCommSubset::all);
  return;
}

TaskStatus RadiationIntegratorTaskList::LocalIntegratorJeans(MeshBlock *pmb, int stage)
{
#ifdef INCLUDE_CHEMISTRY
  pmb->prad->pradintegrator->UpdateRadiation(0);
  pmb->prad->pradintegrator->CopyToOutput();
#endif
  return TaskStatus::success;
}

TaskStatus RadiationIntegratorTaskList::ConstRadiation(MeshBlock *pmb, int stage)
{
#ifdef INCLUDE_CHEMISTRY
  pmb->prad->pradintegrator->CopyToOutput();
#endif
  return TaskStatus::success;
}

