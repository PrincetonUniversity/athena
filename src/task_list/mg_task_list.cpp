//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file mg_task_list.cpp
//  \brief functions for MultigridTaskList class

// C headers

// C++ headers
#include <iostream>   // endl
#include <sstream>    // sstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena++ headers
#include "../athena.hpp"
#include "../globals.hpp"
#include "../mesh/mesh.hpp"
#include "../multigrid/multigrid.hpp"
#include "mg_task_list.hpp"

using namespace MultigridTaskNames; // NOLINT (build/namespace)

//----------------------------------------------------------------------------------------
//! \fn void MultigridTaskList::DoTaskListOneStage(MultigridDriver *pmd)
//  \brief completes all tasks in this list, will not return until all are tasks done

void MultigridTaskList::DoTaskListOneStage(MultigridDriver *pmd) {
  int nmg_left = pmd->GetNumMultigrids();

  for (auto itr = pmd->vmg_.begin(); itr<pmd->vmg_.end(); itr++) {
    Multigrid *pmg = *itr;
    pmg->ts_.Reset(ntasks);
  }

  // cycle through all MeshBlocks and perform all tasks possible
  while (nmg_left > 0) {
    for (auto itr = pmd->vmg_.begin(); itr<pmd->vmg_.end(); itr++) {
      Multigrid *pmg = *itr;
      if (DoAllAvailableTasks(pmg, pmg->ts_) == TaskListStatus::complete) nmg_left--;
    }
  }

  return;
}


//----------------------------------------------------------------------------------------
//! \fn TaskListStatus MultigridTaskList::DoAllAvailableTasks
//  \brief do all tasks that can be done (are not waiting for a dependency to be
//  cleared) in this TaskList, return status.

TaskListStatus MultigridTaskList::DoAllAvailableTasks(Multigrid *pmg, TaskStates &ts) {
  int skip=0;
  TaskStatus ret;

  if (ts.num_tasks_left==0) return TaskListStatus::nothing_to_do;

  for (int i=ts.indx_first_task; i<ntasks; i++) {
    MGTask &taski=task_list_[i];

    if (ts.finished_tasks.IsUnfinished(taski.task_id)) { // task not done
      // check if dependency clear
      if (ts.finished_tasks.CheckDependencies(taski.dependency)) {
        ret=(this->*task_list_[i].TaskFunc)(pmg);
        if (ret!=TaskStatus::fail) { // success
          ts.num_tasks_left--;
          ts.finished_tasks.SetFinished(taski.task_id);
          if (skip==0) ts.indx_first_task++;
          if (ts.num_tasks_left==0) return TaskListStatus::complete;
          if (ret==TaskStatus::next) continue;
          return TaskListStatus::running;
        }
      }
      skip++; // increment number of tasks processed

    } else if (skip==0) { // this task is already done AND it is at the top of the list
      ts.indx_first_task++;
    }
  }
  // there are still tasks to do but nothing can be done now
  return TaskListStatus::stuck;
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridTaskList::AddMultigridTask(const TaskID& id, const TaskID& dep)
//  \brief Sets id and dependency for "ntask" member of task_list_ array, then iterates
//  value of ntask.

void MultigridTaskList::AddMultigridTask(const TaskID& id, const TaskID& dep) {
  task_list_[ntasks].task_id=id;
  task_list_[ntasks].dependency=dep;

  if (id == MG_STARTRECV0 || id == MG_STARTRECVL) {
      task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (MultigridTaskList::*)(Multigrid*)>
        (&MultigridTaskList::StartReceive);
  } else if (id == MG_STARTRECV0F || id == MG_STARTRECV1R || id == MG_STARTRECV1B
          || id == MG_STARTRECV2R || id == MG_STARTRECV2B) {
      task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (MultigridTaskList::*)(Multigrid*)>
        (&MultigridTaskList::StartReceiveFace);
  } else if (id == MG_CLEARBND0 || id == MG_CLEARBNDL) {
      task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (MultigridTaskList::*)(Multigrid*)>
        (&MultigridTaskList::ClearBoundary);
  } else if (id == MG_CLEARBND0F || id == MG_CLEARBND1R || id == MG_CLEARBND1B
          || id == MG_CLEARBND2R || id == MG_CLEARBND2B) {
      task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (MultigridTaskList::*)(Multigrid*)>
        (&MultigridTaskList::ClearBoundaryFace);
  } else if (id == MG_SENDBND0 || id == MG_SENDBNDL) {
      task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (MultigridTaskList::*)(Multigrid*)>
        (&MultigridTaskList::SendBoundary);
  } else if (id == MG_SENDBND0F || id == MG_SENDBND1R || id == MG_SENDBND1B
          || id == MG_SENDBND2R || id == MG_SENDBND2B) {
      task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (MultigridTaskList::*)(Multigrid*)>
        (&MultigridTaskList::SendBoundaryFace);
  } else if (id == MG_RECVBND0 || id == MG_RECVBNDL) {
      task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (MultigridTaskList::*)(Multigrid*)>
        (&MultigridTaskList::ReceiveBoundary);
  } else if (id == MG_RECVBND0F || id == MG_RECVBND1R || id == MG_RECVBND1B
          || id == MG_RECVBND2R || id == MG_RECVBND2B) {
      task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (MultigridTaskList::*)(Multigrid*)>
        (&MultigridTaskList::ReceiveBoundaryFace);
  } else if (id == MG_PRLNGBND0 || id == MG_PRLNGBNDL) {
      task_list_[ntasks].TaskFunc=
          static_cast<TaskStatus (MultigridTaskList::*)(Multigrid*)>
          (&MultigridTaskList::ProlongateBoundary);
  } else if (id == MG_SMOOTH1R || id == MG_SMOOTH2R) {
      task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (MultigridTaskList::*)(Multigrid*)>
        (&MultigridTaskList::SmoothRed);
  } else if (id == MG_SMOOTH1B || id == MG_SMOOTH2B) {
      task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (MultigridTaskList::*)(Multigrid*)>
        (&MultigridTaskList::SmoothBlack);
  } else if (id == MG_PHYSBND0 || id == MG_PHYSBND1R || id == MG_PHYSBND1B
          || id == MG_PHYSBND2R || id == MG_PHYSBND2B || id == MG_PHYSBNDL) {
      task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (MultigridTaskList::*)(Multigrid*)>
        (&MultigridTaskList::PhysicalBoundary);
  } else if (id == MG_RESTRICT) {
      task_list_[ntasks].TaskFunc=
          static_cast<TaskStatus (MultigridTaskList::*)(Multigrid*)>
          (&MultigridTaskList::Restrict);
  } else if (id == MG_PROLONG) {
      task_list_[ntasks].TaskFunc=
          static_cast<TaskStatus (MultigridTaskList::*)(Multigrid*)>
          (&MultigridTaskList::Prolongate);
  } else if (id == MG_FMGPROLONG) {
      task_list_[ntasks].TaskFunc=
          static_cast<TaskStatus (MultigridTaskList::*)(Multigrid*)>
          (&MultigridTaskList::FMGProlongate);
  } else if (id == MG_CALCFASRHS) {
      task_list_[ntasks].TaskFunc=
          static_cast<TaskStatus (MultigridTaskList::*)(Multigrid*)>
          (&MultigridTaskList::CalculateFASRHS);
  } else {
    std::stringstream msg;
    msg << "### FATAL ERROR in AddMultigridTask" << std::endl
        << "Invalid Task is specified" << std::endl;
    ATHENA_ERROR(msg);
  }
  ntasks++;
  return;
}

TaskStatus MultigridTaskList::StartReceive(Multigrid *pmg) {
  pmg->pmgbval->StartReceivingMultigrid(pmg->btype);
  return TaskStatus::success;
}

TaskStatus MultigridTaskList::StartReceiveFace(Multigrid *pmg) {
  pmg->pmgbval->StartReceivingMultigrid(pmg->btypef);
  return TaskStatus::success;
}

TaskStatus MultigridTaskList::ClearBoundary(Multigrid *pmg) {
  pmg->pmgbval->ClearBoundaryMultigrid(pmg->btype);
  return TaskStatus::next;
}

TaskStatus MultigridTaskList::ClearBoundaryFace(Multigrid *pmg) {
  pmg->pmgbval->ClearBoundaryMultigrid(pmg->btypef);
  return TaskStatus::next;
}

TaskStatus MultigridTaskList::SendBoundary(Multigrid *pmg) {
  if (!(pmg->pmgbval->SendMultigridBoundaryBuffers(pmg->btype)))
    return TaskStatus::fail;
  return TaskStatus::success;
}

TaskStatus MultigridTaskList::SendBoundaryFace(Multigrid *pmg) {
  if (!(pmg->pmgbval->SendMultigridBoundaryBuffers(pmg->btypef)))
    return TaskStatus::fail;
  return TaskStatus::success;
}

TaskStatus MultigridTaskList::ReceiveBoundary(Multigrid *pmg) {
  if (!(pmg->pmgbval->ReceiveMultigridBoundaryBuffers(pmg->btype)))
    return TaskStatus::fail;
  return TaskStatus::next;
}

TaskStatus MultigridTaskList::ReceiveBoundaryFace(Multigrid *pmg) {
  if (!(pmg->pmgbval->ReceiveMultigridBoundaryBuffers(pmg->btypef)))
    return TaskStatus::fail;
  return TaskStatus::next;
}

TaskStatus MultigridTaskList::SmoothRed(Multigrid *pmg) {
  pmg->SmoothBlock(0);
  return TaskStatus::next;
}

TaskStatus MultigridTaskList::SmoothBlack(Multigrid *pmg) {
  pmg->SmoothBlock(1);
  return TaskStatus::next;
}

TaskStatus MultigridTaskList::Restrict(Multigrid *pmg) {
  pmg->RestrictBlock();
  return TaskStatus::success;
}

TaskStatus MultigridTaskList::Prolongate(Multigrid *pmg) {
  pmg->ProlongateAndCorrectBlock();
  return TaskStatus::next;
}

TaskStatus MultigridTaskList::FMGProlongate(Multigrid *pmg) {
  pmg->FMGProlongateBlock();
  return TaskStatus::success;
}

TaskStatus MultigridTaskList::CalculateFASRHS(Multigrid *pmg) {
  if (pmy_mgdriver_->current_level_ < pmy_mgdriver_->fmglevel_) {
    pmg->StoreOldData();
    pmg->CalculateFASRHSBlock();
  }
  return TaskStatus::success;
}

TaskStatus MultigridTaskList::PhysicalBoundary(Multigrid *pmg) {
  pmg->pmgbval->ApplyPhysicalBoundaries();
  return TaskStatus::next;
}

TaskStatus MultigridTaskList::ProlongateBoundary(Multigrid *pmg) {
  pmg->pmgbval->ProlongateMultigridBoundaries();
  return TaskStatus::success;
}

//----------------------------------------------------------------------------------------
//! \fn void MultigridTaskList::SetMGTaskListToFiner(int nsmooth, int ngh, int flag)
//  \brief Set the task list for prolongation and post smoothing

void MultigridTaskList::SetMGTaskListToFiner(int nsmooth, int ngh, int flag) {
  ClearTaskList();
  // nsmooth==0 should not be used
  if (flag==1) { // first time on the block level
    AddMultigridTask(MG_PROLONG, NONE);
  } else {
    AddMultigridTask(MG_STARTRECV0, NONE);
    AddMultigridTask(MG_SENDBND0,   MG_STARTRECV0);
    AddMultigridTask(MG_RECVBND0,   MG_STARTRECV0);
    if (pmy_mgdriver_->pmy_mesh_->multilevel) {
      AddMultigridTask(MG_PRLNGBND0, MG_SENDBND0|MG_RECVBND0);
      AddMultigridTask(MG_PHYSBND0,  MG_PRLNGBND0);
    } else {
      AddMultigridTask(MG_PHYSBND0,  MG_SENDBND0|MG_RECVBND0);
    }
    AddMultigridTask(MG_PROLONG,   MG_PHYSBND0);
    AddMultigridTask(MG_CLEARBND0, MG_PROLONG);
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
    AddMultigridTask(MG_SENDBNDL, MG_STARTRECVL);
    AddMultigridTask(MG_RECVBNDL, MG_STARTRECVL);
    if (pmy_mgdriver_->pmy_mesh_->multilevel) {
      AddMultigridTask(MG_PRLNGBNDL, MG_SENDBNDL|MG_RECVBNDL);
      AddMultigridTask(MG_PHYSBNDL,  MG_PRLNGBNDL);
    } else {
      AddMultigridTask(MG_PHYSBNDL,  MG_SENDBNDL|MG_RECVBNDL);
    }
    AddMultigridTask(MG_CLEARBNDL, MG_PHYSBNDL);
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
    if (pmy_mgdriver_->ffas_) {
      AddMultigridTask(MG_CALCFASRHS, MG_PHYSBND0);
      AddMultigridTask(MG_RESTRICT,   MG_CALCFASRHS);
    } else {
      AddMultigridTask(MG_RESTRICT,    MG_PHYSBND0);
    }
    AddMultigridTask(MG_CLEARBND0F,  MG_RESTRICT);
  } else if (nsmooth==1) {
    AddMultigridTask(MG_STARTRECV1R, NONE);
    AddMultigridTask(MG_SENDBND1R,   MG_STARTRECV1R);
    AddMultigridTask(MG_RECVBND1R,   MG_STARTRECV1R);
    AddMultigridTask(MG_PHYSBND1R,   MG_SENDBND1R|MG_RECVBND1R);
    if (pmy_mgdriver_->ffas_) {
      AddMultigridTask(MG_CALCFASRHS,  MG_PHYSBND1R);
      AddMultigridTask(MG_SMOOTH1R,    MG_CALCFASRHS);
    } else {
      AddMultigridTask(MG_SMOOTH1R,    MG_PHYSBND1R);
    }
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
    if (pmy_mgdriver_->ffas_) {
      AddMultigridTask(MG_CALCFASRHS,  MG_PHYSBND1R);
      AddMultigridTask(MG_SMOOTH1R,    MG_CALCFASRHS);
    } else {
      AddMultigridTask(MG_SMOOTH1R,    MG_PHYSBND1R);
    }
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
    if (pmy_mgdriver_->pmy_mesh_->multilevel) {
      AddMultigridTask(MG_PRLNGBND0, MG_SENDBND0|MG_RECVBND0);
      AddMultigridTask(MG_PHYSBND0,  MG_PRLNGBND0);
    } else {
      AddMultigridTask(MG_PHYSBND0,  MG_RECVBND0|MG_SENDBND0);
    }
    AddMultigridTask(MG_FMGPROLONG, MG_PHYSBND0);
    AddMultigridTask(MG_CLEARBND0,  MG_FMGPROLONG);
  }
}
