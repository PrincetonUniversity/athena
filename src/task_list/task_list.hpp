#ifndef TASK_LIST_HPP
#define TASK_LIST_HPP
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//!   \file task_list.hpp
//    \brief provides functionality to control dynamic execution using tasks

#include <stdint.h>

// Athena++ headers
#include "../athena.hpp"

// forward declarations
class Mesh;
class MeshBlock;
class TaskList;
class GravitySolverTaskList;

// return codes for functions working on individual Tasks and TaskList
enum TaskStatus {TASK_FAIL, TASK_SUCCESS, TASK_NEXT};
enum TaskListStatus {TL_RUNNING, TL_STUCK, TL_COMPLETE, TL_NOTHING_TO_DO};

//----------------------------------------------------------------------------------------
//! \struct IntegratorWeight
//  \brief weights used in time integrator tasks

struct IntegratorWeight {
  // 2S or 3S* low-storage RK coefficients, Ketchenson (2010)
  Real delta; // low-storage coefficients to avoid double F() evaluation per substage
  Real gamma_1, gamma_2, gamma_3; // low-storage coeff for weighted ave of registers
  Real beta; // Coefficients from bidiagonal Shu-Osher form Beta matrix, -1 diagonal terms
};

//----------------------------------------------------------------------------------------
//! \struct Task
//  \brief data and function pointer for an individual Task

struct Task {
  uint64_t task_id;      // encodes step & task using bit positions in HydroTasks
  uint64_t dependency;   // encodes dependencies to other tasks using " " " "
  enum TaskStatus (TaskList::*TaskFunc)(MeshBlock*, int);  // ptr to member function
};


//---------------------------------------------------------------------------------------
//! \class TaskState
//  \brief container for task states

class TaskState {
  public:
  uint64_t finished_tasks;
  int indx_first_task, num_tasks_left;
  void Reset(int ntasks) {
    indx_first_task = 0;
    num_tasks_left = ntasks;
    finished_tasks = 0LL;
  };
};


//----------------------------------------------------------------------------------------
//! \class TaskList
//  \brief data and function definitions for task list base class

class TaskList {
friend class TimeIntegratorTaskList;
friend class GravitySolverTaskList;
public:
  TaskList(Mesh *pm);
  virtual ~TaskList();

  // data
  int ntasks;     // number of tasks in this list
  int nsub_steps; // number of times task list should be repeated per full time step

  // functions
  enum TaskListStatus DoAllAvailableTasks(MeshBlock *pmb, int step, TaskState &ts);
  void DoTaskListOneSubstep(Mesh *pmesh, int step);

private:
  Mesh* pmy_mesh_;
  struct Task task_list_[64];
};

//----------------------------------------------------------------------------------------
//! \class TimeIntegratorTaskList
//  \brief data and function definitions for TimeIntegratorTaskList derived class

class TimeIntegratorTaskList : public TaskList {
public:
  TimeIntegratorTaskList(ParameterInput *pin, Mesh *pm);
  ~TimeIntegratorTaskList() {};

  // data
  std::string integrator;
  Real cfl_limit; // dt stability limit for the particular time integrator + spatial order
  struct IntegratorWeight step_wghts[MAX_NSTEP];

  void AddTimeIntegratorTask(uint64_t id, uint64_t dep);

  // functions
  enum TaskStatus StartAllReceive(MeshBlock *pmb, int step);
  enum TaskStatus ClearAllBoundary(MeshBlock *pmb, int step);

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

  enum TaskStatus GravSend(MeshBlock *pmb, int step);
  enum TaskStatus GravReceive(MeshBlock *pmb, int step);
  enum TaskStatus GravSolve(MeshBlock *pmb, int step);
  enum TaskStatus GravFluxCorrection(MeshBlock *pmb, int step);

  enum TaskStatus StartupIntegrator(MeshBlock *pmb, int step);
  enum TaskStatus UpdateTimeStep(MeshBlock *pmb, int step);
};


//----------------------------------------------------------------------------------------
// 64-bit integers with "1" in different bit positions used to ID  each hydro task.

namespace HydroIntegratorTaskNames {
  const uint64_t NONE=0;
  const uint64_t START_ALLRECV=1LL<<0;
  const uint64_t CLEAR_ALLBND=1LL<<1;

  const uint64_t CALC_HYDFLX=1LL<<2;
  const uint64_t CALC_FLDFLX=1LL<<3;
  const uint64_t CALC_RADFLX=1LL<<4;
  const uint64_t CALC_CHMFLX=1LL<<5;

  const uint64_t ADD_VISCFLX=1LL<<6;
  const uint64_t ADD_HEATFLX=1LL<<7;
  const uint64_t ADD_OHMFLX=1LL<<8;
  const uint64_t ADD_ADFLX=1LL<<9;
  const uint64_t ADD_HALLFLX=1LL<<10;

  const uint64_t SEND_HYDFLX=1LL<<11;
  const uint64_t SEND_FLDFLX=1LL<<12;
  const uint64_t SEND_RADFLX=1LL<<13;
  const uint64_t SEND_CHMFLX=1LL<<14;

  const uint64_t RECV_HYDFLX=1LL<<15;
  const uint64_t RECV_FLDFLX=1LL<<16;
  const uint64_t RECV_RADFLX=1LL<<17;
  const uint64_t RECV_CHMFLX=1LL<<18;

  const uint64_t SRCTERM_HYD=1LL<<19;
  const uint64_t SRCTERM_FLD=1LL<<20;
  const uint64_t SRCTERM_RAD=1LL<<21;
  const uint64_t SRCTERM_CHM=1LL<<22;

  const uint64_t INT_HYD=1LL<<23;
  const uint64_t INT_FLD=1LL<<24;
  const uint64_t INT_RAD=1LL<<25;
  const uint64_t INT_CHM=1LL<<26;

  const uint64_t SEND_HYD=1LL<<27;
  const uint64_t SEND_FLD=1LL<<28;
  const uint64_t SEND_RAD=1LL<<29;
  const uint64_t SEND_CHM=1LL<<30;

  const uint64_t RECV_HYD=1LL<<31;
  const uint64_t RECV_FLD=1LL<<32;
  const uint64_t RECV_RAD=1LL<<33;
  const uint64_t RECV_CHM=1LL<<34;

  const uint64_t PROLONG =1LL<<35;
  const uint64_t CON2PRIM=1LL<<36;
  const uint64_t PHY_BVAL=1LL<<37;
  const uint64_t USERWORK=1LL<<38;
  const uint64_t NEW_DT  =1LL<<39;
  const uint64_t AMR_FLAG=1LL<<40;

  const uint64_t SOLV_GRAV=1LL<<41;
  const uint64_t SEND_GRAV=1LL<<42;
  const uint64_t RECV_GRAV=1LL<<43;
  const uint64_t CORR_GFLX=1LL<<44;

  const uint64_t STARTUP_INT=1LL<<45;
  const uint64_t UPDATE_DT=1LL<<46;
};

#endif // TASK_LIST_HPP
