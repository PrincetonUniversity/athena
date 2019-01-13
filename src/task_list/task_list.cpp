//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file task_list.cpp
//  \brief functions for TaskList base class


// C headers

// C++ headers
//#include <vector>
// used to be needed for vector of pointers in DoTaskListOneStage()

// Athena++ headers
#include "../athena.hpp"
#include "../globals.hpp"
#include "../mesh/mesh.hpp"
#include "task_list.hpp"

#ifdef OPENMP_PARALLEL
#include <omp.h>
#endif

//----------------------------------------------------------------------------------------
// TaskList constructor

TaskList::TaskList(Mesh *pm) {
  pmy_mesh_=pm;
  ntasks = 0;
  nstages = 0;
}

// destructor

TaskList::~TaskList() {
}

//----------------------------------------------------------------------------------------
//! \fn enum TaskListStatus TaskList::DoAllAvailableTasks
//  \brief do all tasks that can be done (are not waiting for a dependency to be
//  cleared) in this TaskList, return status.

enum TaskListStatus TaskList::DoAllAvailableTasks(MeshBlock *pmb, int stage,
                                                  TaskState &ts) {
  int skip=0;
  enum TaskStatus ret;
  if (ts.num_tasks_left==0) return TL_NOTHING_TO_DO;

  for (int i=ts.indx_first_task; i<ntasks; i++) {
    Task &taski=task_list_[i];
    if ((taski.task_id & ts.finished_tasks) == 0ULL) { // task not done
      // check if dependency clear
      if (((taski.dependency & ts.finished_tasks) == taski.dependency)) {
        ret=(this->*task_list_[i].TaskFunc)(pmb, stage);
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
//! \fn void TaskList::DoTaskListOneStage(Mesh *pmesh, int stage)
//  \brief completes all tasks in this list, will not return until all are tasks done

void TaskList::DoTaskListOneStage(Mesh *pmesh, int stage) {
  int nthreads = pmesh->GetNumMeshThreads();
#pragma omp parallel num_threads(nthreads)
  {
    int nmb = pmesh->GetNumMeshBlocksThisRank(Globals::my_rank);
    int tid = 0, tis=0;
    int nmbt = nmb / nthreads;
    int nmbres = nmb % nthreads;
    int nmymb = nmbt;

#ifdef OPENMP_PARALLEL
    // calculate the number and index of the MeshBlocks owned by this thread
    tid = omp_get_thread_num();
    if (tid < nmbres) {
      tis = nmbt * tid + tid;
      nmymb++;
    } else {
      tis = nmbt * tid + nmbres;
    }
#endif

    // initialize the task states, initiate MPI and construct the MeshBlock list
    MeshBlock **pmb_array = new MeshBlock*[nmymb];
    MeshBlock *pmb = pmesh->FindMeshBlock(tis+pmesh->pblock->gid);
    for (int n=0; n < nmymb; ++n) {
      pmb->tasks.Reset(ntasks);
      pmb_array[n] = pmb;
      pmb = pmb->next;
    }

    // cycle through all MeshBlocks and perform all tasks possible
    int nmb_left = nmymb;
    StartupTaskList(pmb_array, nmymb, stage);

#pragma omp barrier
    while (nmb_left > 0) {
      for (int i=0; i<nmymb; ++i) {
        if (DoAllAvailableTasks(pmb_array[i], stage, pmb_array[i]->tasks)
            == TL_COMPLETE) {
          nmb_left--;
        }
      }
    }
    delete [] pmb_array;
  } // omp parallel
  return;
}
