//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file mg_task_list.cpp
//  \brief functions for TaskList base class

// Athena++ classes headers
#include "../athena.hpp"
#include "../globals.hpp"
#include "../mesh/mesh.hpp"

// this class header
#include "task_list.hpp"

//----------------------------------------------------------------------------------------
//! \fn void MultigridTaskList::DoTaskList(MultigridDriver *pmd)
//  \brief completes all tasks in this list, will not return until all are tasks done

void MultigridTaskList::DoTaskList(MultigridDriver *pmd)
{
  MeshBlock *pmb = pmd->pblock;
  int nmb_left = pmd->GetNumMeshBlocks();

  while (pmb != NULL)  {
    pmb->ts.Reset(ntasks);
    pmb=pmb->next;
  }

  // cycle through all MeshBlocks and perform all tasks possible
  while(nmb_left > 0) {
    pmb = pmd->pblock;
    while (pmb != NULL)  {
      if (DoAllAvailableTasks(pmb,step,pmd->ts_) == TL_COMPLETE) nmb_left--;
      pmb=pmb->next;
    }
  }

  return;
}


//----------------------------------------------------------------------------------------//! \fn
//  \brief Sets id and dependency for "ntask" member of task_list_ array, then iterates
//  value of ntask.  

void MultigridTaskList::AddTask(uint64_t id, uint64_t dep)
{
  task_list_[ntasks].task_id=id;
  task_list_[ntasks].dependency=dep;

  using namespace HydroIntegratorTaskNames;
  switch((id)) {
    case (MGGRAV_STARTRECV):
      task_list_[ntasks].TaskFunc= 
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&MultigridTaskList::MGGravStartReceive);
      break;
    case (MGGRAV_CLEARBND):
      task_list_[ntasks].TaskFunc= 
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&MultigridTaskList::MGGravClearBoundary);
      break;

    default:
      std::stringstream msg;
      msg << "### FATAL ERROR in AddTask" << std::endl
          << "Invalid Task "<< id << " is specified" << std::endl;
      throw std::runtime_error(msg.str().c_str());
  }
  ntasks++;
  return;
}



enum TaskStatus MultigridTaskList::MGGravityStartReceive(MeshBlock *pmb, int step)
{
  int nc=pmb->pmggrav->GetCurrentNumberOfCells();
  pmb->pbval->StartReceivingMultigrid(nc, BND_MGGRAV);
  return TASK_SUCCESS;
}

enum TaskStatus MultigridTaskList::MGGravityClearBoundary(MeshBlock *pmb, int step)
{
  pmb->pbval->ClearBoundaryMultigrid(BND_MGGRAV);
  return TASK_SUCCESS;
}

enum TaskStatus MultigridTaskList::MGGravitySendBoundary(MeshBlock *pmb, int step)
{
  MGGravity *pmg=pmb->pmggrav;
  int nc=pmg->GetCurrentNumberOfCells();
  pmb->pbval->SendMultigridBoundaryBuffers(pmg->GetCurrentData(), nc, BND_MGGRAV);
  return TASK_SUCCESS;
}


enum TaskStatus MultigridTaskList::MGGravityReceiveBoundary(MeshBlock *pmb, int step)
{
  MGGravity *pmg=pmb->pmggrav;
  int nc=pmg->GetCurrentNumberOfCells();
  if(pmb->pbval->ReceiveMultigridBoundaryBuffers(pmg->GetCurrentData(), nc, BND_MGGRAV)
     == true) return TASK_NEXT;
  else return TASK_FAIL;
}
