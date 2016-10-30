//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file tasklist.cpp
//  \brief functions for TaskList base class

// Athena++ classes headers
#include "../athena.hpp"
#include "../globals.hpp"
#include "../mesh/mesh.hpp"

// this class header
#include "task_list.hpp"

//----------------------------------------------------------------------------------------
// TaskList constructor

TaskList::TaskList(Mesh *pm)
{
  pmy_mesh_ = pm;
  ntasks = 0;
  nsub_steps = 0;
}

// destructor

TaskList::~TaskList()
{
}

//----------------------------------------------------------------------------------------
//! \fn DoAllAvailableTasks
//  \brief do all tasks that can be done (are not waiting for a dependency to be 
//  cleared) in this TaskList, return status.  

enum TaskListStatus TaskList::DoAllAvailableTasks(MeshBlock *pmb, int step) {
  int skip=0;
  enum TaskStatus ret;

  if(pmb->num_tasks_left_==0) return TL_NOTHING_TO_DO;

  for(int i=pmb->indx_first_task_; i<ntasks; i++) {
    Task &taski=task_list_[i];

    if((taski.task_id & pmb->finished_tasks) == 0LL) { // task not done
      // check if dependency clear
      if (((taski.dependency & pmb->finished_tasks) == taski.dependency)) {
        ret=(this->*task_list_[i].TaskFunc)(pmb,step);
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

//----------------------------------------------------------------------------------------
//! \fn void TaskList::DoTaskList(Mesh *pmesh)
//  \brief completes all tasks in this list, will not return until all are tasks done

void TaskList::DoTaskList(Mesh *pmesh)
{

  for (int step=1; step<=nsub_steps; ++step) {
    MeshBlock *pmb = pmesh->pblock;
    int nmb_left = pmesh->GetNumMeshBlocksThisRank(Globals::my_rank);
    // initialize counters stored in each MeshBlock
    while (pmb != NULL)  {
      pmb->indx_first_task_ = 0;
      pmb->num_tasks_left_ = ntasks;
      pmb->finished_tasks = 0LL; // encodes which tasks are done
      pmb=pmb->next;
    }

    // cycle through all MeshBlocks and perform all tasks possible
    while(nmb_left > 0) {
      pmb = pmesh->pblock;
      while (pmb != NULL)  {
        if (DoAllAvailableTasks(pmb,step) == TL_COMPLETE) nmb_left--;
        pmb=pmb->next;
      }
    }

  }

  return;
}
