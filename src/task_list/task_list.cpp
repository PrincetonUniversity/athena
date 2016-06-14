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
//! \file tasklist.hpp
//  \brief task functions
//======================================================================================

// Athena++ classes headers
#include "../athena.hpp"
#include "../globals.hpp"
#include "../mesh/mesh.hpp"
#include "../hydro/hydro.hpp"
#include "../field/field.hpp"
#include "../bvals/bvals.hpp"
#include "../eos/eos.hpp"
#include "../hydro/srcterms/hydro_srcterms.hpp"

// this class header
#include "task_list.hpp"

//--------------------------------------------------------------------------------------
// TaskList constructor

TaskList::TaskList(Mesh *pm)
{
  pmy_mesh_ = pm;
  ntasks = 0;

  // hardwired for VL2 integrator for now
  nloop_over_list = 2;
  CreateTimeIntegrator(pm);

}

// destructor

TaskList::~TaskList()
{
}

//--------------------------------------------------------------------------------------
//! \fn
//  \brief do all possible tasks in this TaskList, return status

enum TaskListStatus TaskList::DoAllTasksPossible(MeshBlock *pmb, int step) {
  int skip=0;
  enum TaskStatus ret;

  if(pmb->num_tasks_left_==0) return TL_NOTHING_TO_DO;

  for(int i=pmb->indx_first_task_; i<ntasks; i++) {
    Task &taski=task_list_[i];

    if((taski.task_id & pmb->finished_tasks) == 0LL) { // task not done
      // check if dependency clear
      if (((taski.dependency & pmb->finished_tasks) == taski.dependency)) {
        ret=taski.TaskFunc(pmb,step);
        if(ret!=TASK_FAIL) { // success
          pmb->num_tasks_left_--;
          pmb->finished_tasks |= taski.task_id;
          if(skip==0) pmb->indx_first_task_++;
          if(pmb->num_tasks_left_==0) return TL_COMPLETE;
          if(ret==TASK_NEXT) continue;
          return TL_RUNNING;
        }
      }
      skip++; // increment number of tasks processed

    } else if(skip==0) // this task is already done AND it is at the top of the list
      pmb->indx_first_task_++;
  }
  return TL_STUCK; // there are still tasks to do but nothing can be done now
}

//--------------------------------------------------------------------------------------
//! \fn void TaskList::ExecuteTaskList(Mesh *pmesh)
//  \brief completes all tasks in this list

void TaskList::ExecuteTaskList(Mesh *pmesh)
{

  for (int step=1; step<=nloop_over_list; ++step) {
    MeshBlock *pmb = pmesh->pblock;
    int nmb_left = pmesh->GetNumMeshBlocksThisRank(Globals::my_rank);
    // initialize, start MPI communications (if needed)
    while (pmb != NULL)  {
      pmb->indx_first_task_ = 0;
      pmb->num_tasks_left_ = ntasks;
      pmb->finished_tasks = 0; // encodes which tasks are done
      pmb->pbval->StartReceivingAll();
      pmb=pmb->next;
    }

    // cycle through all MeshBlocks and perform all tasks possible
    while(nmb_left > 0) {
      pmb = pmesh->pblock;
      while (pmb != NULL)  {
        if (DoAllTasksPossible(pmb,step) == TL_COMPLETE) nmb_left--;
        pmb=pmb->next;
      }
    }

    // clear boundary buffers
    pmb = pmesh->pblock;
    while (pmb != NULL)  {
      pmb->pbval->ClearBoundaryAll();
      pmb=pmb->next;
    }
  }

  return;
}
