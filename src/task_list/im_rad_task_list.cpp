//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file task_list.cpp
//! \brief functions for TaskList base class


// C headers

// C++ headers
//#include <vector> // formerly needed for vector of MeshBlock ptrs in DoTaskListOneStage

// Athena++ headers
#include "../athena.hpp"
#include "../globals.hpp"
#include "../mesh/mesh.hpp"
#include "../nr_radiation/radiation.hpp"
#include "im_rad_task_list.hpp"

#ifdef OPENMP_PARALLEL
#include <omp.h>
#endif

//----------------------------------------------------------------------------------------
//! \fn TaskListStatus TaskList::DoAllAvailableTasks
//! \brief do all tasks that can be done (are not waiting for a dependency to be
//! cleared) in this TaskList, return status.

TaskListStatus IMRadTaskList::DoAllAvailableTasks(MeshBlock *pmb, TaskStates &ts) {
  int skip = 0;
  TaskStatus ret;
  if (ts.num_tasks_left == 0) return TaskListStatus::nothing_to_do;

  for (int i=ts.indx_first_task; i<ntasks; i++) {
    IMRadTask &taski = task_list_[i];
    if (ts.finished_tasks.IsUnfinished(taski.task_id)) { // task not done
      // check if dependency clear
      if (ts.finished_tasks.CheckDependencies(taski.dependency)) {
        ret = (this->*task_list_[i].TaskFunc)(pmb);
        if (ret != TaskStatus::fail) { // success
          ts.num_tasks_left--;
          ts.finished_tasks.SetFinished(taski.task_id);
          if (skip == 0) ts.indx_first_task++;
          if (ts.num_tasks_left == 0) return TaskListStatus::complete;
          if (ret == TaskStatus::next) continue;
          return TaskListStatus::running;
        }
      }
      skip++; // increment number of tasks processed

    } else if (skip == 0) { // this task is already done AND it is at the top of the list
      ts.indx_first_task++;
    }
  }
  // there are still tasks to do but nothing can be done now
  return TaskListStatus::stuck;
}

//----------------------------------------------------------------------------------------
//! \fn void TaskList::DoTaskListOneStage(Mesh *pmesh, int stage)
//! \brief completes all tasks in this list, will not return until all are tasks done

void IMRadTaskList::DoTaskListOneStage(Real wght) {
  time = pmy_mesh->time + wght;
  dt = wght;
  int nmb = pmy_mesh->nblocal;
  for (int i=0; i<nmb; ++i) {
    pmy_mesh->my_blocks(i)->tasks.Reset(ntasks);
    StartupTaskList(pmy_mesh->my_blocks(i));
  }

  int nmb_left = nmb;
  // cycle through all MeshBlocks and perform all tasks possible
  while (nmb_left > 0) {
    //! \note
    //! KNOWN ISSUE: Workaround for unknown OpenMP race condition. See #183 on GitHub.
    for (int i=0; i<nmb; ++i) {
      if (DoAllAvailableTasks(pmy_mesh->my_blocks(i), pmy_mesh->my_blocks(i)->tasks)
          == TaskListStatus::complete) {
        nmb_left--;
      }
    }
  }
  return;
}

TaskStatus IMRadTaskList::PhysicalBoundary(MeshBlock *pmb) {
  pmb->pnrrad->rad_bvar.var_cc = &(pmb->pnrrad->ir);
  pmb->pbval->ApplyPhysicalBoundaries(time, dt, pmb->pbval->bvars_main_int);
  return TaskStatus::success;
}


TaskStatus IMRadTaskList::ProlongateBoundary(MeshBlock *pmb) {
  pmb->pbval->ProlongateBoundaries(time, dt, pmb->pbval->bvars_main_int);
  return TaskStatus::success;
}


