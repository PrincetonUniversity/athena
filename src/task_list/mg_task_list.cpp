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

using namespace MultigridTaskNames; // NOLINT (build/namespace)

//----------------------------------------------------------------------------------------
//! \fn void MultigridTaskList::DoTaskListOneStage(MultigridDriver *pmd)
//  \brief completes all tasks in this list, will not return until all are tasks done

void MultigridTaskList::DoTaskListOneStage(MultigridDriver *pmd) {
  Multigrid *pmg = pmd->pmg_;
  int nmg_left = pmd->GetNumMultigrids();

  while (pmg != NULL)  {
    pmg->ts_.Reset(ntasks);
    pmg=pmg->next;
  }

  // cycle through all MeshBlocks and perform all tasks possible
  while(nmg_left > 0) {
    pmg = pmd->pmg_;
    while (pmg != NULL)  {
      if (DoAllAvailableTasks(pmg, pmg->ts_) == TL_COMPLETE) nmg_left--;
      pmg=pmg->next;
    }
  }

  return;
}


//----------------------------------------------------------------------------------------
//! \fn enum TaskListStatus MultigridTaskList::DoAllAvailableTasks
//  \brief do all tasks that can be done (are not waiting for a dependency to be
//  cleared) in this TaskList, return status.

enum TaskListStatus MultigridTaskList::DoAllAvailableTasks(Multigrid *pmg,
                                                           TaskState &ts) {
  int skip=0;
  enum TaskStatus ret;

  if (ts.num_tasks_left==0) return TL_NOTHING_TO_DO;

  for (int i=ts.indx_first_task; i<ntasks; i++) {
    MGTask &taski=task_list_[i];

    if ((taski.task_id & ts.finished_tasks) == 0LL) { // task not done
      // check if dependency clear
      if (((taski.dependency & ts.finished_tasks) == taski.dependency)) {
        ret=(this->*task_list_[i].TaskFunc)(pmg);
        if (ret!=TASK_FAIL) { // success
          ts.num_tasks_left--;
          ts.finished_tasks |= taski.task_id;
          if (skip==0) ts.indx_first_task++;
          if (ts.num_tasks_left==0) return TL_COMPLETE;
          if (ret==TASK_NEXT) continue;
          return TL_RUNNING;
        }
      }
      skip++; // increment number of tasks processed

    } else if (skip==0) { // this task is already done AND it is at the top of the list
      ts.indx_first_task++;
    }
  }
  return TL_STUCK; // there are still tasks to do but nothing can be done now
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridTaskList::AddMultigridTask(uint64_t id, uint64_t dep)
//  \brief Sets id and dependency for "ntask" member of task_list_ array, then iterates
//  value of ntask.

void MultigridTaskList::AddMultigridTask(uint64_t id, uint64_t dep) {
  task_list_[ntasks].task_id=id;
  task_list_[ntasks].dependency=dep;

  switch(id) {
    case (MG_STARTRECV0):
    case (MG_STARTRECVL):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (MultigridTaskList::*)(Multigrid*)>
                                    (&MultigridTaskList::StartReceive);
      break;
    case (MG_STARTRECV0F):
    case (MG_STARTRECV1R):
    case (MG_STARTRECV1B):
    case (MG_STARTRECV2R):
    case (MG_STARTRECV2B):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (MultigridTaskList::*)(Multigrid*)>
                                    (&MultigridTaskList::StartReceiveFace);
      break;
    case (MG_CLEARBND0):
    case (MG_CLEARBNDL):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (MultigridTaskList::*)(Multigrid*)>
                                    (&MultigridTaskList::ClearBoundary);
      break;
    case (MG_CLEARBND0F):
    case (MG_CLEARBND1R):
    case (MG_CLEARBND1B):
    case (MG_CLEARBND2R):
    case (MG_CLEARBND2B):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (MultigridTaskList::*)(Multigrid*)>
                                    (&MultigridTaskList::ClearBoundaryFace);
      break;
    case (MG_SENDBND0):
    case (MG_SENDBNDL):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (MultigridTaskList::*)(Multigrid*)>
                                    (&MultigridTaskList::SendBoundary);
      break;
    case (MG_SENDBND0F):
    case (MG_SENDBND1R):
    case (MG_SENDBND1B):
    case (MG_SENDBND2R):
    case (MG_SENDBND2B):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (MultigridTaskList::*)(Multigrid*)>
                                    (&MultigridTaskList::SendBoundaryFace);
      break;
    case (MG_RECVBND0):
    case (MG_RECVBNDL):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (MultigridTaskList::*)(Multigrid*)>
                                    (&MultigridTaskList::ReceiveBoundary);
      break;
    case (MG_RECVBND0F):
    case (MG_RECVBND1R):
    case (MG_RECVBND1B):
    case (MG_RECVBND2R):
    case (MG_RECVBND2B):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (MultigridTaskList::*)(Multigrid*)>
                                    (&MultigridTaskList::ReceiveBoundaryFace);
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
    case (MG_PHYSBND0):
    case (MG_PHYSBND1R):
    case (MG_PHYSBND1B):
    case (MG_PHYSBND2R):
    case (MG_PHYSBND2B):
    case (MG_PHYSBNDL):
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

enum TaskStatus MultigridTaskList::StartReceive(Multigrid *pmg) {
  int nc=pmg->GetCurrentNumberOfCells();
  pmg->pmgbval->StartReceivingMultigrid(nc, pmg->btype);
  return TASK_SUCCESS;
}

enum TaskStatus MultigridTaskList::StartReceiveFace(Multigrid *pmg) {
  int nc=pmg->GetCurrentNumberOfCells();
  pmg->pmgbval->StartReceivingMultigrid(nc, pmg->btypef);
  return TASK_SUCCESS;
}

enum TaskStatus MultigridTaskList::ClearBoundary(Multigrid *pmg) {
  pmg->pmgbval->ClearBoundaryMultigrid(pmg->btype);
  return TASK_NEXT;
}

enum TaskStatus MultigridTaskList::ClearBoundaryFace(Multigrid *pmg) {
  pmg->pmgbval->ClearBoundaryMultigrid(pmg->btypef);
  return TASK_NEXT;
}

enum TaskStatus MultigridTaskList::SendBoundary(Multigrid *pmg) {
  int nc=pmg->GetCurrentNumberOfCells();
  if (pmg->pmgbval->
     SendMultigridBoundaryBuffers(pmg->GetCurrentData(), nc, pmg->btype)==false)
    return TASK_FAIL;
  return TASK_SUCCESS;
}

enum TaskStatus MultigridTaskList::SendBoundaryFace(Multigrid *pmg) {
  int nc=pmg->GetCurrentNumberOfCells();
  if (pmg->pmgbval->
     SendMultigridBoundaryBuffers(pmg->GetCurrentData(), nc, pmg->btypef)==false)
    return TASK_FAIL;
  return TASK_SUCCESS;
}

enum TaskStatus MultigridTaskList::ReceiveBoundary(Multigrid *pmg) {
  int nc=pmg->GetCurrentNumberOfCells();
  if (pmg->pmgbval->
     ReceiveMultigridBoundaryBuffers(pmg->GetCurrentData(), nc, pmg->btype)==false)
    return TASK_FAIL;
  return TASK_NEXT;
}

enum TaskStatus MultigridTaskList::ReceiveBoundaryFace(Multigrid *pmg) {
  int nc=pmg->GetCurrentNumberOfCells();
  if (pmg->pmgbval->
     ReceiveMultigridBoundaryBuffers(pmg->GetCurrentData(), nc, pmg->btypef)==false)
    return TASK_FAIL;
  return TASK_NEXT;
}

enum TaskStatus MultigridTaskList::SmoothRed(Multigrid *pmg) {
  pmg->Smooth(0);
  return TASK_NEXT;
}

enum TaskStatus MultigridTaskList::SmoothBlack(Multigrid *pmg) {
  pmg->Smooth(1);
  return TASK_NEXT;
}

enum TaskStatus MultigridTaskList::Restrict(Multigrid *pmg) {
  pmg->Restrict();
  return TASK_NEXT;
}

enum TaskStatus MultigridTaskList::Prolongate(Multigrid *pmg) {
  pmg->ProlongateAndCorrect();
  return TASK_NEXT;
}

enum TaskStatus MultigridTaskList::FMGProlongate(Multigrid *pmg) {
  pmg->FMGProlongate();
  return TASK_NEXT;
}

enum TaskStatus MultigridTaskList::PhysicalBoundary(Multigrid *pmg) {
  pmg->pmgbval->ApplyPhysicalBoundaries();
  return TASK_NEXT;
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridTaskList::SetMGTaskListToFiner(int nsmooth, int ngh, int flag)
//  \brief Set the task list for prolongation and post smoothing

void MultigridTaskList::SetMGTaskListToFiner(int nsmooth, int ngh, int flag) {
  ClearTaskList();
  // nsmooth==0 should not be used
  if (flag==1) { // first time on the block level
    AddMultigridTask(MG_PROLONG,    NONE);
  } else {
    AddMultigridTask(MG_STARTRECV0, NONE);
    AddMultigridTask(MG_SENDBND0,   MG_STARTRECV0);
    AddMultigridTask(MG_RECVBND0,   MG_STARTRECV0);
    AddMultigridTask(MG_PHYSBND0,   MG_SENDBND0|MG_RECVBND0);
    AddMultigridTask(MG_PROLONG,    MG_PHYSBND0);
    AddMultigridTask(MG_CLEARBND0,  MG_PROLONG);
  }
  if (nsmooth==1) {
    if (flag==1)
      AddMultigridTask(MG_STARTRECV1R, MG_PROLONG);
    else
      AddMultigridTask(MG_STARTRECV1R, MG_CLEARBND0);
    AddMultigridTask(MG_SENDBND1R,   MG_STARTRECV1R);
    AddMultigridTask(MG_RECVBND1R,   MG_STARTRECV1R);
    AddMultigridTask(MG_PHYSBND1R,   MG_SENDBND1R|MG_RECVBND1R);
    AddMultigridTask(MG_SMOOTH1R,    MG_PHYSBND1R);
    AddMultigridTask(MG_CLEARBND1R,  MG_SMOOTH1R);
    AddMultigridTask(MG_STARTRECV1B, MG_CLEARBND1R);
    AddMultigridTask(MG_SENDBND1B,   MG_STARTRECV1B);
    AddMultigridTask(MG_RECVBND1B,   MG_STARTRECV1B);
    AddMultigridTask(MG_PHYSBND1B,   MG_SENDBND1B|MG_RECVBND1B);
    AddMultigridTask(MG_SMOOTH1B,    MG_PHYSBND1B);
    AddMultigridTask(MG_CLEARBND1B,  MG_SMOOTH1B);
  }
  if (nsmooth==2) {
    AddMultigridTask(MG_STARTRECV2R, MG_CLEARBND1B);
    AddMultigridTask(MG_SENDBND2R,   MG_STARTRECV2R);
    AddMultigridTask(MG_RECVBND2R,   MG_STARTRECV2R);
    AddMultigridTask(MG_PHYSBND2R,   MG_SENDBND2R|MG_RECVBND2R);
    AddMultigridTask(MG_SMOOTH2R,    MG_PHYSBND2R);
    AddMultigridTask(MG_CLEARBND2R,  MG_SMOOTH2R);
    AddMultigridTask(MG_STARTRECV2B, MG_CLEARBND2R);
    AddMultigridTask(MG_SENDBND2B,   MG_STARTRECV2B);
    AddMultigridTask(MG_RECVBND2B,   MG_STARTRECV2B);
    AddMultigridTask(MG_PHYSBND2B,   MG_SENDBND2B|MG_RECVBND2B);
    AddMultigridTask(MG_SMOOTH2B,    MG_PHYSBND2B);
    AddMultigridTask(MG_CLEARBND2B,  MG_SMOOTH2B);
  }
  if (flag==2) { // last
    if (nsmooth==1)
      AddMultigridTask(MG_STARTRECVL, MG_CLEARBND1B);
    else if (nsmooth==2)
      AddMultigridTask(MG_STARTRECVL, MG_CLEARBND2B);
    AddMultigridTask(MG_SENDBNDL,   MG_STARTRECVL);
    AddMultigridTask(MG_RECVBNDL,   MG_STARTRECVL);
    AddMultigridTask(MG_PHYSBNDL,   MG_SENDBNDL|MG_RECVBNDL);
    AddMultigridTask(MG_CLEARBNDL,  MG_PHYSBNDL);
  }
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridTaskList::SetMGTaskListToCoarser(int nsmooth, int ngh)
//  \brief Set the task list for pre smoothing and restriction

void MultigridTaskList::SetMGTaskListToCoarser(int nsmooth, int ngh) {
  ClearTaskList();
  if (nsmooth==0) {
    AddMultigridTask(MG_STARTRECV0F, NONE);
    AddMultigridTask(MG_SENDBND0F,   MG_STARTRECV0F);
    AddMultigridTask(MG_RECVBND0F,   MG_STARTRECV0F);
    AddMultigridTask(MG_PHYSBND0,    MG_SENDBND0F|MG_RECVBND0F);
    AddMultigridTask(MG_RESTRICT,    MG_PHYSBND0);
    AddMultigridTask(MG_CLEARBND0F,  MG_RESTRICT);
  } else if (nsmooth==1) {
    AddMultigridTask(MG_STARTRECV1R, NONE);
    AddMultigridTask(MG_SENDBND1R,   MG_STARTRECV1R);
    AddMultigridTask(MG_RECVBND1R,   MG_STARTRECV1R);
    AddMultigridTask(MG_PHYSBND1R,   MG_SENDBND1R|MG_RECVBND1R);
    AddMultigridTask(MG_SMOOTH1R,    MG_PHYSBND1R);
    AddMultigridTask(MG_CLEARBND1R,  MG_SMOOTH1R);
    AddMultigridTask(MG_STARTRECV1B, MG_CLEARBND1R);
    AddMultigridTask(MG_SENDBND1B,   MG_STARTRECV1B);
    AddMultigridTask(MG_RECVBND1B,   MG_STARTRECV1B);
    AddMultigridTask(MG_PHYSBND1B,   MG_SENDBND1B|MG_RECVBND1B);
    AddMultigridTask(MG_SMOOTH1B,    MG_PHYSBND1B);
    AddMultigridTask(MG_CLEARBND1B,  MG_SMOOTH1B);
    AddMultigridTask(MG_STARTRECV0F, MG_CLEARBND1B);
    AddMultigridTask(MG_SENDBND0F,   MG_STARTRECV0F);
    AddMultigridTask(MG_RECVBND0F,   MG_STARTRECV0F);
    AddMultigridTask(MG_PHYSBND0,    MG_SENDBND0F|MG_RECVBND0F);
    AddMultigridTask(MG_RESTRICT,    MG_PHYSBND0);
    AddMultigridTask(MG_CLEARBND0F,  MG_RESTRICT);
  } else if (nsmooth==2) {
    AddMultigridTask(MG_STARTRECV1R, NONE);
    AddMultigridTask(MG_SENDBND1R,   MG_STARTRECV1R);
    AddMultigridTask(MG_RECVBND1R,   MG_STARTRECV1R);
    AddMultigridTask(MG_PHYSBND1R,   MG_SENDBND1R|MG_RECVBND1R);
    AddMultigridTask(MG_SMOOTH1R,    MG_PHYSBND1R);
    AddMultigridTask(MG_CLEARBND1R,  MG_SMOOTH1R);
    AddMultigridTask(MG_STARTRECV1B, MG_CLEARBND1R);
    AddMultigridTask(MG_SENDBND1B,   MG_STARTRECV1B);
    AddMultigridTask(MG_RECVBND1B,   MG_STARTRECV1B);
    AddMultigridTask(MG_PHYSBND1B,   MG_SENDBND1B|MG_RECVBND1B);
    AddMultigridTask(MG_SMOOTH1B,    MG_PHYSBND1B);
    AddMultigridTask(MG_CLEARBND1B,  MG_SMOOTH1B);
    AddMultigridTask(MG_STARTRECV2R, MG_CLEARBND1B);
    AddMultigridTask(MG_SENDBND2R,   MG_STARTRECV2R);
    AddMultigridTask(MG_RECVBND2R,   MG_STARTRECV2R);
    AddMultigridTask(MG_PHYSBND2R,   MG_SENDBND2R|MG_RECVBND2R);
    AddMultigridTask(MG_SMOOTH2R,    MG_PHYSBND2R);
    AddMultigridTask(MG_CLEARBND2R,  MG_SMOOTH2R);
    AddMultigridTask(MG_STARTRECV2B, MG_CLEARBND2R);
    AddMultigridTask(MG_SENDBND2B,   MG_STARTRECV2B);
    AddMultigridTask(MG_RECVBND2B,   MG_STARTRECV2B);
    AddMultigridTask(MG_PHYSBND2B,   MG_SENDBND2B|MG_RECVBND2B);
    AddMultigridTask(MG_SMOOTH2B,    MG_PHYSBND2B);
    AddMultigridTask(MG_CLEARBND2B,  MG_SMOOTH2B);
    AddMultigridTask(MG_STARTRECV0F, MG_CLEARBND2B);
    AddMultigridTask(MG_SENDBND0F,   MG_STARTRECV0F);
    AddMultigridTask(MG_RECVBND0F,   MG_STARTRECV0F);
    AddMultigridTask(MG_PHYSBND0,    MG_SENDBND0F|MG_RECVBND0F);
    AddMultigridTask(MG_RESTRICT,    MG_PHYSBND0);
    AddMultigridTask(MG_CLEARBND0F,  MG_RESTRICT);
  }
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridTaskList::SetMGTaskListFMGProlongate(int flag)
//  \brief Set the task list for FMG prolongation

void MultigridTaskList::SetMGTaskListFMGProlongate(int flag) {
  ClearTaskList();
  if (flag==1) { // first time on the block level
    AddMultigridTask(MG_FMGPROLONG,    NONE);
  } else {
    AddMultigridTask(MG_STARTRECV0, NONE);
    AddMultigridTask(MG_SENDBND0,   MG_STARTRECV0);
    AddMultigridTask(MG_RECVBND0,   MG_STARTRECV0);
    AddMultigridTask(MG_PHYSBND0,   MG_RECVBND0|MG_SENDBND0);
    AddMultigridTask(MG_FMGPROLONG, MG_PHYSBND0);
    AddMultigridTask(MG_CLEARBND0,  MG_FMGPROLONG);
  }
}
