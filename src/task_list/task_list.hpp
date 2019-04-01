#ifndef TASK_LIST_TASK_LIST_HPP_
#define TASK_LIST_TASK_LIST_HPP_
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//!   \file task_list.hpp
//    \brief provides functionality to control dynamic execution using tasks

// C headers

// C++ headers
#include <cstdint>      // std::uint64_t
#include <string>       // std::string

// Athena++ headers
#include "../athena.hpp"

// forward declarations
class Mesh;
class MeshBlock;
class TaskList;
class FFTGravitySolverTaskList;

// TODO(felker): these 4x declarations can be nested in TaskList if MGTaskList is derived

// constants = return codes for functions working on individual Tasks and TaskList
enum class TaskStatus {fail, success, next};
enum class TaskListStatus {running, stuck, complete, nothing_to_do};

//----------------------------------------------------------------------------------------
//! \struct Task
//  \brief data and function pointer for an individual Task

struct Task { // aggregate and POD
  std::uint64_t task_id;    // encodes task with bit positions in HydroIntegratorTaskNames
  std::uint64_t dependency; // encodes dependencies to other tasks using " " " "
  TaskStatus (TaskList::*TaskFunc)(MeshBlock*, int);  // ptr to member function
};

//---------------------------------------------------------------------------------------
//! \struct TaskStates
//  \brief container for task states on a single MeshBlock

struct TaskStates { // aggregate and POD
  std::uint64_t finished_tasks;
  int indx_first_task, num_tasks_left;
  void Reset(int ntasks) {
    indx_first_task = 0;
    num_tasks_left = ntasks;
    finished_tasks = 0ULL;
  }
};

//----------------------------------------------------------------------------------------
//! \class TaskList
//  \brief data and function definitions for task list base class

class TaskList {
 public:
  TaskList() : ntasks(0), nstages(0), task_list_{} {} // 2x direct + zero initialization
  // rule of five:
  virtual ~TaskList() = default;

  // data
  int ntasks;     // number of tasks in this list
  int nstages;    // number of times the tasklist is repeated per each full timestep

  // functions
  TaskListStatus DoAllAvailableTasks(MeshBlock *pmb, int stage, TaskStates &ts);
  void DoTaskListOneStage(Mesh *pmesh, int stage);

 protected:
  Task task_list_[64];

 private:
  virtual void AddTask(std::uint64_t id, std::uint64_t dep) = 0;
  virtual void StartupTaskList(MeshBlock *pmb, int stage) = 0;
};

//----------------------------------------------------------------------------------------
//! \class TimeIntegratorTaskList
//  \brief data and function definitions for TimeIntegratorTaskList derived class

class TimeIntegratorTaskList : public TaskList {
 public:
  TimeIntegratorTaskList(ParameterInput *pin, Mesh *pm);
  ~TimeIntegratorTaskList() {}

  //--------------------------------------------------------------------------------------
  //! \struct IntegratorWeight
  //  \brief weights used in time integrator tasks

  struct IntegratorWeight {
    // 2S or 3S* low-storage RK coefficients, Ketchenson (2010)
    Real delta; // low-storage coefficients to avoid double F() evaluation per substage
    Real gamma_1, gamma_2, gamma_3; // low-storage coeff for weighted ave of registers
    Real beta; // coeff. from bidiagonal Shu-Osher form Beta matrix, -1 diagonal terms
  };

  // data
  std::string integrator;
  Real cfl_limit; // dt stability limit for the particular time integrator + spatial order

  // functions
  TaskStatus ClearAllBoundary(MeshBlock *pmb, int stage);

  TaskStatus HydroFluxCalculate(MeshBlock *pmb, int stage);
  TaskStatus EMFCalculate(MeshBlock *pmb, int stage);

  TaskStatus HydroFluxCorrectSend(MeshBlock *pmb, int stage);
  TaskStatus EMFCorrectSend(MeshBlock *pmb, int stage);

  TaskStatus HydroFluxCorrectReceive(MeshBlock *pmb, int stage);
  TaskStatus EMFCorrectReceive(MeshBlock *pmb, int stage);

  TaskStatus HydroIntegrate(MeshBlock *pmb, int stage);
  TaskStatus FieldIntegrate(MeshBlock *pmb, int stage);

  TaskStatus HydroSourceTerms(MeshBlock *pmb, int stage);

  TaskStatus HydroDiffusion(MeshBlock *pmb, int stage);
  TaskStatus FieldDiffusion(MeshBlock *pmb, int stage);

  TaskStatus HydroSend(MeshBlock *pmb, int stage);
  TaskStatus FieldSend(MeshBlock *pmb, int stage);

  TaskStatus HydroReceive(MeshBlock *pmb, int stage);
  TaskStatus FieldReceive(MeshBlock *pmb, int stage);

  TaskStatus HydroSetBoundaries(MeshBlock *pmb, int stage);
  TaskStatus FieldSetBoundaries(MeshBlock *pmb, int stage);

  TaskStatus HydroShearSend(MeshBlock *pmb, int stage);
  TaskStatus HydroShearReceive(MeshBlock *pmb, int stage);
  TaskStatus FieldShearSend(MeshBlock *pmb, int stage);
  TaskStatus FieldShearReceive(MeshBlock *pmb, int stage);
  TaskStatus EMFShearSend(MeshBlock *pmb, int stage);
  TaskStatus EMFShearReceive(MeshBlock *pmb, int stage);
  TaskStatus EMFShearRemap(MeshBlock *pmb, int stage);

  TaskStatus Prolongation(MeshBlock *pmb, int stage);
  TaskStatus Primitives(MeshBlock *pmb, int stage);
  TaskStatus PhysicalBoundary(MeshBlock *pmb, int stage);
  TaskStatus UserWork(MeshBlock *pmb, int stage);
  TaskStatus NewBlockTimeStep(MeshBlock *pmb, int stage);
  TaskStatus CheckRefinement(MeshBlock *pmb, int stage);

 private:
  IntegratorWeight stage_wghts[MAX_NSTAGE];

  void AddTask(std::uint64_t id, std::uint64_t dep) override;
  void StartupTaskList(MeshBlock *pmb, int stage) override;
};

//----------------------------------------------------------------------------------------
//! \class SuperTimeStepTaskList
//  \brief data and function definitions for SuperTimeStepTaskList derived class

class SuperTimeStepTaskList : public TaskList {
 public:
  SuperTimeStepTaskList(ParameterInput *pin, Mesh *pm, TimeIntegratorTaskList *ptlist);
  ~SuperTimeStepTaskList() {}

  // functions
  TaskStatus HydroFluxCalculate_STS(MeshBlock *pmb, int stage);
  TaskStatus EMFCalculate_STS(MeshBlock *pmb, int stage);

  TaskStatus HydroIntegrate_STS(MeshBlock *pmb, int stage);
  TaskStatus FieldIntegrate_STS(MeshBlock *pmb, int stage);

  TaskStatus PhysicalBoundary_STS(MeshBlock *pmb, int stage);

 private:
  // currently intiialized but unused. May use it for direct calls to TimeIntegrator fns:
  TimeIntegratorTaskList *ptlist_;
  void AddTask(std::uint64_t id, std::uint64_t dep) override;
  void StartupTaskList(MeshBlock *pmb, int stage) override;
};

//----------------------------------------------------------------------------------------
// 64-bit integers with "1" in different bit positions used to ID each hydro task.

namespace HydroIntegratorTaskNames {
const std::uint64_t NONE          = 0ULL;
const std::uint64_t CLEAR_ALLBND  = 1ULL<<1;

const std::uint64_t CALC_HYDFLX = 1ULL<<2;
const std::uint64_t CALC_FLDFLX = 1ULL<<3;
const std::uint64_t CALC_RADFLX = 1ULL<<4;
const std::uint64_t CALC_CHMFLX = 1ULL<<5;

const std::uint64_t SEND_HYDFLX = 1ULL<<6;
const std::uint64_t SEND_FLDFLX = 1ULL<<7;
const std::uint64_t SEND_RADFLX = 1ULL<<8;
const std::uint64_t SEND_CHMFLX = 1ULL<<9;

const std::uint64_t RECV_HYDFLX = 1ULL<<10;
const std::uint64_t RECV_FLDFLX = 1ULL<<11;
const std::uint64_t RECV_RADFLX = 1ULL<<12;
const std::uint64_t RECV_CHMFLX = 1ULL<<13;

const std::uint64_t SRCTERM_HYD = 1ULL<<14;
const std::uint64_t SRCTERM_FLD = 1ULL<<15;
const std::uint64_t SRCTERM_RAD = 1ULL<<16;
const std::uint64_t SRCTERM_CHM = 1ULL<<17;

const std::uint64_t INT_HYD = 1ULL<<18;
const std::uint64_t INT_FLD = 1ULL<<19;
const std::uint64_t INT_RAD = 1ULL<<20;
const std::uint64_t INT_CHM = 1ULL<<21;

const std::uint64_t SEND_HYD = 1ULL<<22;
const std::uint64_t SEND_FLD = 1ULL<<23;
const std::uint64_t SEND_RAD = 1ULL<<24;
const std::uint64_t SEND_CHM = 1ULL<<25;

const std::uint64_t RECV_HYD = 1ULL<<26;
const std::uint64_t RECV_FLD = 1ULL<<27;
const std::uint64_t RECV_RAD = 1ULL<<28;
const std::uint64_t RECV_CHM = 1ULL<<29;

const std::uint64_t SETB_HYD = 1ULL<<30;
const std::uint64_t SETB_FLD = 1ULL<<31;
const std::uint64_t SETB_RAD = 1ULL<<32;
const std::uint64_t SETB_CHM = 1ULL<<33;

const std::uint64_t PROLONG  = 1ULL<<34;
const std::uint64_t CON2PRIM = 1ULL<<35;
const std::uint64_t PHY_BVAL = 1ULL<<36;
const std::uint64_t USERWORK = 1ULL<<37;
const std::uint64_t NEW_DT   = 1ULL<<38;
const std::uint64_t AMR_FLAG = 1ULL<<39;

const std::uint64_t SEND_HYDSH = 1ULL<<40;
const std::uint64_t SEND_EMFSH = 1ULL<<41;
const std::uint64_t SEND_FLDSH = 1ULL<<42;
const std::uint64_t RECV_HYDSH = 1ULL<<43;
const std::uint64_t RECV_EMFSH = 1ULL<<44;
const std::uint64_t RECV_FLDSH = 1ULL<<45;
const std::uint64_t RMAP_EMFSH = 1ULL<<46;

const std::uint64_t DIFFUSE_HYD      = 1ULL<<47;
const std::uint64_t DIFFUSE_FLD      = 1ULL<<48;
const std::uint64_t CALC_DIFFUSIVITY = 1ULL<<49;

const std::uint64_t CALC_SCLRFLX    = 1ULL<<50;
const std::uint64_t SEND_SCLRFLX    = 1ULL<<51;
const std::uint64_t RECV_SCLRFLX    = 1ULL<<52;
const std::uint64_t INT_SCLR       = 1ULL<<53;
const std::uint64_t SEND_SCLR      = 1ULL<<54;
const std::uint64_t RECV_SCLR      = 1ULL<<55;
const std::uint64_t SETB_SCLR      = 1ULL<<56;

} // namespace HydroIntegratorTaskNames
#endif // TASK_LIST_TASK_LIST_HPP_
