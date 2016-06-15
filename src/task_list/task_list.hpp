#ifndef TASK_LIST_HPP
#define TASK_LIST_HPP
//======================================================================================
// Athena++ astrophysical MHD code
// Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
// See LICENSE file for full public license information.
//======================================================================================
//!   \file tasklist.hpp
//    \brief provides functionality to control dynamic execution using tasks
//======================================================================================

// Athena++ headers
#include "../athena.hpp"
#include "../mesh/mesh.hpp"

// forward declarations
class Mesh;
class MeshBlock;
class TaskList;

// return codes for functions working on individual Tasks and TaskList
enum TaskStatus {TASK_FAIL, TASK_SUCCESS, TASK_NEXT};
enum TaskListStatus {TL_RUNNING, TL_STUCK, TL_COMPLETE, TL_NOTHING_TO_DO};

//--------------------------------------------------------------------------------------
//! \struct IntegratorWeight
//  \brief weights used in time integrator tasks

struct IntegratorWeight {
  Real a,b,c;
};

//--------------------------------------------------------------------------------------
//! \struct Task
//  \brief data and function pointer for an individual Task

struct Task {
  uint64_t task_id;      // encodes step & task using bit positions in HydroTasks
  uint64_t dependency;   // encodes dependencies to other tasks using " " " "
  enum TaskStatus (TaskList::*TaskFunc)(MeshBlock*, int);  // ptr to member function
};

//--------------------------------------------------------------------------------------
//! \class TaskList
//  \brief data and function definitions for task list base class

class TaskList {
friend class TimeIntegratorTaskList;
public:
  TaskList(Mesh *pm);
  ~TaskList();

  // data
  int ntasks;     // number of tasks in this list
  int nsub_steps; // number of times task list should be repeated per full time step

  // functions
  enum TaskListStatus DoAllAvailableTasks(MeshBlock *pmb, int step);
  void DoTaskList(Mesh *pmesh);

private:
  Mesh* pmy_mesh_;
  struct Task task_list_[64];
};

//--------------------------------------------------------------------------------------
//! \class TimeIntegratorTaskList
//  \brief data and function definitions for TimeIntegratorTaskList derived class

class TimeIntegratorTaskList : public TaskList {
public:
  TimeIntegratorTaskList(ParameterInput *pin, Mesh *pm);
  ~TimeIntegratorTaskList() {};

  // data
  std::string integrator;
  struct IntegratorWeight step_wghts[MAX_NSTEP];

  // functions
  void AddTimeIntegratorTask(uint64_t id, uint64_t dep);

  enum TaskStatus CalculateFluxes(MeshBlock *pmb, int step);
  enum TaskStatus CalculateEMF(MeshBlock *pmb, int step);

  enum TaskStatus FluxCorrectSend(MeshBlock *pmb, int step);
  enum TaskStatus EMFCorrectSend(MeshBlock *pmb, int step);

  enum TaskStatus FluxCorrectReceive(MeshBlock *pmb, int step);
  enum TaskStatus EMFCorrectReceive(MeshBlock *pmb, int step);

  enum TaskStatus HydroIntegrate(MeshBlock *pmb, int step);
  enum TaskStatus FieldIntegrate(MeshBlock *pmb, int step);

  enum TaskStatus HydroSourceTerms(MeshBlock *pmb, int step);

  enum TaskStatus HydroSend(MeshBlock *pmb, int step);
  enum TaskStatus FieldSend(MeshBlock *pmb, int step);

  enum TaskStatus HydroReceive(MeshBlock *pmb, int step);
  enum TaskStatus FieldReceive(MeshBlock *pmb, int step);

  enum TaskStatus Prolongation(MeshBlock *pmb, int step);
  enum TaskStatus Primitives(MeshBlock *pmb, int step);
  enum TaskStatus PhysicalBoundary(MeshBlock *pmb, int step);
  enum TaskStatus UserWork(MeshBlock *pmb, int step);
  enum TaskStatus NewBlockTimeStep(MeshBlock *pmb, int step);
  enum TaskStatus CheckRefinement(MeshBlock *pmb, int step);
};

//--------------------------------------------------------------------------------------
// 64-bit integers with "1" in different bit positions used to ID  each hydro task.

namespace HydroIntegratorTaskNames {
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

#endif // TASK_LIST_HPP
