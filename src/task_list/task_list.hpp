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
typedef enum TaskStatus (*TaskFunc_t)(MeshBlock*, int);

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

// 32-bit integers with "1" in different bit positions to label each hydro task.
namespace IntegratorTask {
  const uint64_t NONE=0;
  const uint64_t CALC_HYDFLX=1LL<<0;
  const uint64_t CALC_FLDFLX=1LL<<1;
  const uint64_t CALC_RADFLX=1LL<<2;
  const uint64_t CALC_CHMFLX=1LL<<3;

  const uint64_t ADD_VISCFLX=1LL<<4;
  const uint64_t ADD_HEATFLX=1LL<<5;
  const uint64_t ADD_OHMFLX=1LL<<6;
  const uint64_t ADD_ADFLX=1LL<<7;
  const uint64_t ADD_HALLFLX=1LL<<8;

  const uint64_t SEND_HYDFLX=1LL<<9;
  const uint64_t SEND_FLDFLX=1LL<<10;
  const uint64_t SEND_RADFLX=1LL<<11;
  const uint64_t SEND_CHMFLX=1LL<<12;

  const uint64_t RECV_HYDFLX=1LL<<13;
  const uint64_t RECV_FLDFLX=1LL<<14;
  const uint64_t RECV_RADFLX=1LL<<15;
  const uint64_t RECV_CHMFLX=1LL<<16;

  const uint64_t SRCTERM_HYD=1LL<<17;
  const uint64_t SRCTERM_FLD=1LL<<18;
  const uint64_t SRCTERM_RAD=1LL<<19;
  const uint64_t SRCTERM_CHM=1LL<<20;

  const uint64_t INT_HYD=1LL<<21;
  const uint64_t INT_FLD=1LL<<22;
  const uint64_t INT_RAD=1LL<<23;
  const uint64_t INT_CHM=1LL<<24;

  const uint64_t SEND_HYD=1LL<<25;
  const uint64_t SEND_FLD=1LL<<26;
  const uint64_t SEND_RAD=1LL<<27;
  const uint64_t SEND_CHM=1LL<<28;

  const uint64_t RECV_HYD=1LL<<29;
  const uint64_t RECV_FLD=1LL<<30;
  const uint64_t RECV_RAD=1LL<<31;
  const uint64_t RECV_CHM=1LL<<32;

  const uint64_t PROLONG =1LL<<33;
  const uint64_t CON2PRIM=1LL<<34;
  const uint64_t PHY_BVAL=1LL<<35;
  const uint64_t USERWORK=1LL<<36;
  const uint64_t NEW_DT  =1LL<<37;
  const uint64_t AMR_FLAG=1LL<<38;
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
  void CreateTimeIntegrator(Mesh *pm);

private:
  Mesh* pmy_mesh_;
  Task task_list_[64];
};

#endif // TASK_LIST_HPP
