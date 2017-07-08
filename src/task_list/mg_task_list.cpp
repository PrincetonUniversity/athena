//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file mg_task_list.cpp
//  \brief functions for MultigridTaskList class

// C/C++ headers
#include <iostream>   // endl
#include <sstream>    // sstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena++ classes headers
#include "../athena.hpp"
#include "../globals.hpp"
#include "../mesh/mesh.hpp"
#include "../multigrid/multigrid.hpp"

// this class header
#include "mg_task_list.hpp"

using namespace MultigridTaskNames;

//----------------------------------------------------------------------------------------
//! \fn void MultigridTaskList::DoTaskListOneSubStep(MultigridDriver *pmd)
//  \brief completes all tasks in this list, will not return until all are tasks done

void MultigridTaskList::DoTaskListOneSubStep(MultigridDriver *pmd)
{
  MeshBlock *pmb = pmd->pblock_;
  int nmb_left = pmd->GetNumMeshBlocks();

  while (pmb != NULL)  {
    Multigrid *pmg=pmd->GetMultigridBlock(pmb);
    pmg->ts_.Reset(ntasks);
    pmb=pmb->next;
  }

  // cycle through all MeshBlocks and perform all tasks possible
  while(nmb_left > 0) {
    pmb = pmd->pblock_;
    while (pmb != NULL)  {
      Multigrid *pmg=pmd->GetMultigridBlock(pmb);
      if (DoAllAvailableTasks(pmg, pmg->ts_) == TL_COMPLETE) nmb_left--;
      pmb=pmb->next;
    }
  }

  return;
}


//----------------------------------------------------------------------------------------
//! \fn enum TaskListStatus MultigridTaskList::DoAllAvailableTasks
//  \brief do all tasks that can be done (are not waiting for a dependency to be 
//  cleared) in this TaskList, return status.  

enum TaskListStatus MultigridTaskList::DoAllAvailableTasks(Multigrid *pmg, TaskState &ts)
{
  int skip=0;
  enum TaskStatus ret;

  if(ts.num_tasks_left==0) return TL_NOTHING_TO_DO;

  for(int i=ts.indx_first_task; i<ntasks; i++) {
    MGTask &taski=task_list_[i];

    if((taski.task_id & ts.finished_tasks) == 0LL) { // task not done
      // check if dependency clear
      if (((taski.dependency & ts.finished_tasks) == taski.dependency)) {
        ret=(this->*task_list_[i].TaskFunc)(pmg);
        if(ret!=TASK_FAIL) { // success
          ts.num_tasks_left--;
          ts.finished_tasks |= taski.task_id;
          if(skip==0) ts.indx_first_task++;
          if(ts.num_tasks_left==0) return TL_COMPLETE;
          if(ret==TASK_NEXT) continue;
          return TL_RUNNING;
        }
      }
      skip++; // increment number of tasks processed

    } else if(skip==0) // this task is already done AND it is at the top of the list
      ts.indx_first_task++;
  }
  return TL_STUCK; // there are still tasks to do but nothing can be done now
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridTaskList::AddMultigridTask(uint64_t id, uint64_t dep)
//  \brief Sets id and dependency for "ntask" member of task_list_ array, then iterates
//  value of ntask.  

void MultigridTaskList::AddMultigridTask(uint64_t id, uint64_t dep)
{
  task_list_[ntasks].task_id=id;
  task_list_[ntasks].dependency=dep;

  switch(id) {
    case (MG_STARTRECV1R):
    case (MG_STARTRECV1B):
    case (MG_STARTRECV2R):
    case (MG_STARTRECV2B):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (MultigridTaskList::*)(Multigrid*)>
                                    (&MultigridTaskList::StartReceive);
      break;
    case (MG_CLEARBND1R):
    case (MG_CLEARBND1B):
    case (MG_CLEARBND2R):
    case (MG_CLEARBND2B):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (MultigridTaskList::*)(Multigrid*)>
                                    (&MultigridTaskList::ClearBoundary);
      break;
    case (MG_SENDBND1R):
    case (MG_SENDBND1B):
    case (MG_SENDBND2R):
    case (MG_SENDBND2B):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (MultigridTaskList::*)(Multigrid*)>
                                    (&MultigridTaskList::SendBoundary);
      break;
    case (MG_RECVBND1R):
    case (MG_RECVBND1B):
    case (MG_RECVBND2R):
    case (MG_RECVBND2B):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (MultigridTaskList::*)(Multigrid*)>
                                    (&MultigridTaskList::ReceiveBoundary);
      break;
    case (MG_SMOOTH1R):
    case (MG_SMOOTH2R):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (MultigridTaskList::*)(Multigrid*)>
                                    (&MultigridTaskList::SmoothRed);
      break;
    case (MG_SMOOTH1B):
    case (MG_SMOOTH2B):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (MultigridTaskList::*)(Multigrid*)>
                                    (&MultigridTaskList::SmoothBlack);
      break;
    case (MG_PHYSBND1):
    case (MG_PHYSBND2):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (MultigridTaskList::*)(Multigrid*)>
                                    (&MultigridTaskList::PhysicalBoundary);
      break;
    case (MG_RESTRICT):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (MultigridTaskList::*)(Multigrid*)>
                                    (&MultigridTaskList::Restrict);
      break;
    case (MG_PROLONG):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (MultigridTaskList::*)(Multigrid*)>
                                    (&MultigridTaskList::Prolongate);
      break;
    case (MG_FMGPROLONG):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (MultigridTaskList::*)(Multigrid*)>
                                    (&MultigridTaskList::FMGProlongate);
      break;

    default:
      std::stringstream msg;
      msg << "### FATAL ERROR in AddMultigridTask" << std::endl
          << "Invalid Task "<< id << " is specified" << std::endl;
      throw std::runtime_error(msg.str().c_str());
  }
  ntasks++;
  return;
}

enum TaskStatus MultigridTaskList::StartReceive(Multigrid *pmg)
{
  int nc=pmg->GetCurrentNumberOfCells();
  pmg->pmy_block_->pbval->StartReceivingMultigrid(nc, pmg->btype);
  return TASK_SUCCESS;
}

enum TaskStatus MultigridTaskList::ClearBoundary(Multigrid *pmg)
{
  pmg->pmy_block_->pbval->ClearBoundaryMultigrid(pmg->btype);
  return TASK_NEXT;
}

enum TaskStatus MultigridTaskList::SendBoundary(Multigrid *pmg)
{
  int nc=pmg->GetCurrentNumberOfCells();
  pmg->pmy_block_->pbval->
    SendMultigridBoundaryBuffers(pmg->GetCurrentData(), nc, pmg->btype);
  return TASK_SUCCESS;
}

enum TaskStatus MultigridTaskList::ReceiveBoundary(Multigrid *pmg)
{
  int nc=pmg->GetCurrentNumberOfCells();
  if(pmg->pmy_block_->pbval->
     ReceiveMultigridBoundaryBuffers(pmg->GetCurrentData(), nc, pmg->btype)== false)
    return TASK_FAIL;
  return TASK_NEXT;
}

enum TaskStatus MultigridTaskList::SmoothRed(Multigrid *pmg)
{
  pmg->Smooth(0);
  return TASK_NEXT;
}

enum TaskStatus MultigridTaskList::SmoothBlack(Multigrid *pmg)
{
  pmg->Smooth(1);
  return TASK_NEXT;
}

enum TaskStatus MultigridTaskList::Restrict(Multigrid *pmg)
{
  pmg->Restrict();
  return TASK_NEXT;
}

enum TaskStatus MultigridTaskList::Prolongate(Multigrid *pmg)
{
  pmg->ProlongateAndCorrect();
  return TASK_NEXT;
}

enum TaskStatus MultigridTaskList::FMGProlongate(Multigrid *pmg)
{
  pmg->FMGProlongate();
  return TASK_NEXT;
}

enum TaskStatus MultigridTaskList::PhysicalBoundary(Multigrid *pmg)
{
  return TASK_NEXT;
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridTaskList::SetMGTaskListToFiner(int nsmooth)
//  \brief Sets id and dependency for "ntask" member of task_list_ array, then iterates
//  value of ntask.  

void MultigridTaskList::SetMGTaskListToFiner(int nsmooth)
{
  ntasks=0;
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridTaskList::SetMGTaskListToCoarser(int nsmooth)
//  \brief Sets id and dependency for "ntask" member of task_list_ array, then iterates
//  value of ntask.  

void MultigridTaskList::SetMGTaskListToCoarser(int nsmooth)
{
  ntasks=0;
  if(nsmooth==0) {
    AddMultigridTask(MG_STARTRECV1R, NONE);
  }
  else if(nsmooth==1) {
    AddMultigridTask(MG_STARTRECV1R, NONE);
  }
  else if(nsmooth==2) {
    AddMultigridTask(MG_STARTRECV1R, NONE);
  }
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridTaskList::SetMGTaskListFMGProlongate(void)
//  \brief Sets id and dependency for "ntask" member of task_list_ array, then iterates
//  value of ntask.  

void MultigridTaskList::SetMGTaskListFMGProlongate(void)
{
  ntasks=0;
}
