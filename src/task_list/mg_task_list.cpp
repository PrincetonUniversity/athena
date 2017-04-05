//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file mg_task_list.cpp
//  \brief functions for MultigridTaskList class

// Athena++ classes headers
#include "../athena.hpp"
#include "../globals.hpp"
#include "../mesh/mesh.hpp"
#include "../multigrid/multigrid.hpp"

// this class header
#include "mg_task_list.hpp"

//----------------------------------------------------------------------------------------
//! \fn void MultigridTaskList::DoTaskList(MultigridDriver *pmd)
//  \brief completes all tasks in this list, will not return until all are tasks done

void MultigridTaskList::DoTaskList(MultigridDriver *pmd)
{
  MeshBlock *pmb = pmd->pblock;
  int nmb_left = pmd->GetNumMeshBlocks();

  while (pmb != NULL)  {
    Multigrid *pmg=pmd->GetMultigridBlock(pmb);
    pmg->ts.Reset(ntasks);
    pmb=pmb->next;
  }

  // cycle through all MeshBlocks and perform all tasks possible
  while(nmb_left > 0) {
    pmb = pmd->pblock;
    while (pmb != NULL)  {
      Multigrid *pmg=pmd->GetMultigridBlock(pmb);
      if (DoAllAvailableTasks(pmg,pmg->ts_) == TL_COMPLETE) nmb_left--;
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

  using namespace MultigridTaskNames;
  switch(id) {
    case (MG_STARTRECV):
      task_list_[ntasks].TaskFunc=MultigridTaskList::StartReceive;
      break;
    case (MG_CLEARBND):
      task_list_[ntasks].TaskFunc=MultigridTaskList::ClearBoundary;
      break;
    case (MG_SENDRED):
      task_list_[ntasks].TaskFunc=MultigridTaskList::SendBoundary;
      break;
    case (MG_RECVRED):
      task_list_[ntasks].TaskFunc=MultigridTaskList::ReceiveBoundary;
      break;
    case (MG_SENDBLACK):
      task_list_[ntasks].TaskFunc=MultigridTaskList::SendBoundary;
      break;
    case (MG_RECVBLACK):
      task_list_[ntasks].TaskFunc=MultigridTaskList::ReceiveBoundary;
      break;
    case (MG_SMOOTHRED):
      task_list_[ntasks].TaskFunc=MultigridTaskList::SmoothRed;
      break;
    case (MG_SMOOTHBLACK):
      task_list_[ntasks].TaskFunc=MultigridTaskList::SmoothBlack;
      break;
    case (MG_CALCDEFECT):
      task_list_[ntasks].TaskFunc=MultigridTaskList::CalculateDefect;
      break;
    case (MG_RESTRICT):
      task_list_[ntasks].TaskFunc=MultigridTaskList::Restrict;
      break;
    case (MG_PROLONG):
      task_list_[ntasks].TaskFunc=MultigridTaskList::Prolongate;
      break;
    case (MG_FMGPROLONG):
      task_list_[ntasks].TaskFunc=MultigridTaskList::FMGProlongate;
      break;
    case (MG_PHYSBND):
      task_list_[ntasks].TaskFunc=MultigridTaskList::PhysicalBoundary;
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



enum TaskStatus MultigridTaskList::StartReceive(Multigrid *pmg, int step)
{
  int nc=pmg->GetCurrentNumberOfCells();
  pmg->pmb->pbval->StartReceivingMultigrid(nc, BND_MGGRAV);
  return TASK_SUCCESS;
}

enum TaskStatus MultigridTaskList::ClearBoundary(Multigrid *pmg, int step)
{
  pmg->pmb->pbval->ClearBoundaryMultigrid(BND_MGGRAV);
  return TASK_SUCCESS;
}

enum TaskStatus MultigridTaskList::SendBoundary(Multigrid *pmg, int step)
{
  int nc=pmg->GetCurrentNumberOfCells();
  pmg->pmy_block->pbval->
    SendMultigridBoundaryBuffers(pmg->GetCurrentData(), nc, BND_MGGRAV);
  return TASK_SUCCESS;
}


enum TaskStatus MultigridTaskList::ReceiveBoundary(Multigrid *pmg, int step)
{
  int nc=pmg->GetCurrentNumberOfCells();
  if(pmg->pmy_block->pbval->
     ReceiveMultigridBoundaryBuffers(pmg->GetCurrentData(), nc, BND_MGGRAV)== true)
    return TASK_NEXT;
  else return TASK_FAIL;
}


enum TaskStatus MultigridTaskList::SmoothRed(Multigrid *pmg, int step)
{
  pmg->Smooth(0);
  return TASK_NEXT;
}


enum TaskStatus MultigridTaskList::SmoothBlack(Multigrid *pmg, int step)
{
  pmg->Smooth(1);
  return TASK_NEXT;
}


enum TaskStatus MultigridTaskList::CalculateDefect(Multigrid *pmg, int step)
{
  pmg->CalculateDefect();
  return TASK_SUCCESS;
}


enum TaskStatus MultigridTaskList::Restrict(Multigrid *pmg, int step)
{
  pmg->Restrict();
  pmg->ZeroClearData();
  return TASK_NEXT;
}


enum TaskStatus MultigridTaskList::Prolongate(Multigrid *pmg, int step)
{
  pmg->ProlongateAndCorrect();
  return TASK_NEXT;
}


enum TaskStatus MultigridTaskList::FMGProlongate(Multigrid *pmg, int step)
{
  pmg->FMGProlongate();
  return TASK_NEXT;
}

enum TaskStatus MultigridTaskList::PhysicalBoundary(Multigrid *pmg, int step)
{
  pmg->ApplyPhysicalBoundary();
  return TASK_NEXT;
}

