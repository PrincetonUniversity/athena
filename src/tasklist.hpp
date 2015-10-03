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
#include "athena.hpp"
#include "mesh.hpp"

// forward declarations
class MeshBlock;

// return codes for functions working on individual Tasks
enum task_status {task_failure, task_success, task_do_next};

// return codes for TaskList
enum tasklist_status {tl_running, tl_stuck, tl_complete, tl_nothing};

// 64-bit pattern that labels each TaskFunction.  Used to express dependencies.
enum task_code {
  none=0,
  hydro_integrate_1=1L<<0,
  calculate_emf_1=1L<<1,
  field_integrate_1=1L<<2,
  hydro_send_1=1L<<3,
  hydro_recv_1=1L<<4,
  flux_correct_send_1=1L<<5,
  flux_correct_recv_1=1L<<6,
  hydro_prolong_1=1L<<7,
  hydro_boundary_1=1L<<8,
  field_send_1=1L<<9,
  field_recv_1=1L<<10,
  emf_correct_send_1=1L<<11,
  emf_correct_recv_1=1L<<12,
  field_prolong_1=1L<<13,
  field_boundary_1=1L<<14,
  primitives_1=1L<<15,

  hydro_integrate_0=1L<<16,
  calculate_emf_0=1L<<17,
  field_integrate_0=1L<<18,
  hydro_send_0=1L<<19,
  hydro_recv_0=1L<<20,
  flux_correct_send_0=1L<<21,
  flux_correct_recv_0=1L<<22,
  hydro_prolong_0=1L<<23,
  hydro_boundary_0=1L<<24,
  field_send_0=1L<<25,
  field_recv_0=1L<<26,
  emf_correct_send_0=1L<<27,
  emf_correct_recv_0=1L<<28,
  field_prolong_0=1L<<29,
  field_boundary_0=1L<<30,
  primitives_0=1L<<31,

  new_blocktimestep=1L<<32
};

// definition of TaskFunc_t
typedef enum task_status (*TaskFunc_t)(MeshBlock*, int);

//!   \struct Task
//    \brief all data related to an individual Task
struct Task {
  enum task_code task_id;        // unique 64-bit pattern that labels task
  unsigned long int dependency;  // dependencies to other tasks encoded via task_id's
  TaskFunc_t TaskFunc;           // function called by this task
  int task_flag;                 // flag passed to task function
};

//!   \class TaskList
//    \brief data and function definitions for task list

class TaskList {
public:
  TaskList() : ntasks_(0) {};
  ~TaskList() {};
  void AddTask(enum task_code, unsigned long int dependency);
private:
  int ntasks_;
  Task task_list_[64];
  friend class MeshBlock;
};

#endif // TASK_LIST_HPP
