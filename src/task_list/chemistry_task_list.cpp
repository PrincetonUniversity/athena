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
//! \file chemistry_task_list.cpp
//  \brief derived class for chemistry integrator task list.
//======================================================================================

// C/C++ headers
#include <iostream>   // endl
#include <sstream>    // sstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena++ classes headers
#include "../athena.hpp"
#include "../mesh/mesh.hpp"
#include "../scalars/scalars.hpp" 
#ifdef INCLUDE_CHEMISTRY
#include "../chemistry/network/network.hpp" 
#endif
#include "chemistry_task_list.hpp"
#include "task_list.hpp"

//--------------------------------------------------------------------------------------
//  ChemistryIntegratorTaskList constructor
ChemistryIntegratorTaskList::ChemistryIntegratorTaskList(ParameterInput *pin, Mesh *pm)
{
  // Now assemble list of tasks for each step of chemistry integrator
  {using namespace ChemistryIntegratorTaskNames;
    AddTask(INT_CHEM_SRC,NONE);
  } // end of using namespace block
}

//--------------------------------------------------------------------------------------
//! \fn
//  \brief Sets id and dependency for "ntask" member of task_list_ array, then iterates
//  value of ntask.  
void ChemistryIntegratorTaskList::AddTask(const TaskID& id, const TaskID& dep)
{
  task_list_[ntasks].task_id=id;
  task_list_[ntasks].dependency=dep;

  using namespace ChemistryIntegratorTaskNames;
  if (id == INT_CHEM_SRC) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&ChemistryIntegratorTaskList::IntegrateSourceTerm);
  } else {
    std::stringstream msg;
    msg << "### FATAL ERROR ChemistryIntegratorTaskList::AddTask" << std::endl
        << "Invalid Task is specified" << std::endl;
    ATHENA_ERROR(msg);
  }
  ntasks++;
  return;
}

void ChemistryIntegratorTaskList::StartupTaskList(MeshBlock *pmb, int stage) {
  return;
}

TaskStatus ChemistryIntegratorTaskList::IntegrateSourceTerm(MeshBlock *pmb, int stage)
{
#ifdef INCLUDE_CHEMISTRY
  pmb->pscalars->odew.Integrate();
#endif
  return TaskStatus::success;
}

