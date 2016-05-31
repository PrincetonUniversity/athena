#ifndef TASK_LIST_HPP
#define TASK_LIST_HPP
//======================================================================================
// Athena++ astrophysical MHD code
// Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
// See LICENSE file for full public license information.
//======================================================================================
//!   \file tasklist.hpp
//    \brief definition for class TaskList
// Provides functions to control dynamic execution using tasks
//======================================================================================

// Athena++ headers
#include "../athena.hpp"
#include "../mesh/mesh.hpp"

// forward declarations
class Mesh;
class MeshBlock;

// return codes for functions working on individual Tasks and TaskList
enum TaskStatus {TASK_FAIL, TASK_SUCCESS, TASK_NEXT};
enum TaskListStatus {TL_RUNNING, TL_STUCK, TL_COMPLETE, TL_NOTHING_TO_DO};

// definition of TaskFunc_t
typedef enum TaskStatus (*TaskFunc_t)(MeshBlock*, unsigned long int, int);

// 32-bit integers with "1" in different bit positions to label each hydro task.
enum HydroTasks {
  NONE=0,
  CALC_FLX=1L<<0,
  FLX_SEND=1L<<1,
  FLX_RECV=1L<<2,
  CALC_EMF=1L<<3,
  EMF_SEND=1L<<4,
  EMF_RECV=1L<<5,
  HYD_INT =1L<<6,
  HYD_SEND=1L<<7,
  HYD_RECV=1L<<8,
  FLD_INT =1L<<9,
  FLD_SEND=1L<<10,
  FLD_RECV=1L<<11,
  PROLONG =1L<<12,
  CON2PRIM=1L<<13,
  PHY_BVAL=1L<<14,
  USERWORK=1L<<15,
  NEW_DT  =1L<<16,
  AMR_FLAG=1L<<17,
};

//!   \struct Task
//    \brief all data related to an individual Task

struct Task {
  int step_of_task, step_of_depend; 
  unsigned long int task_id;    // encodes step & task using bit positions in HydroTasks
  unsigned long int dependency; // encodes dependencies to other tasks using " " " "
  TaskFunc_t TaskFunc;          // function called by this task
};

//!   \class TaskList
//    \brief data and function definitions for task list

class TaskList {
public:
  TaskList(Mesh *pm);
  ~TaskList();

  // data
  int ntasks;

  // functions
  void AddTask(int st, unsigned long int id, int sd, unsigned long int dep);
  enum TaskListStatus DoOneTask(MeshBlock *pmb);
  void CreateVL2Integrator(Mesh *pm);

private:
  Mesh* pmy_mesh_;
  Task task_list_[64];
};

#endif // TASK_LIST_HPP
