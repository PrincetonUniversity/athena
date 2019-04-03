#ifndef TASK_LIST_FFT_GRAV_TASK_LIST_HPP_
#define TASK_LIST_FFT_GRAV_TASK_LIST_HPP_
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//!   \file fft_grav_task_list.hpp
//    \brief

// C headers

// C++ headers
#include <cstdint>      // std::uint64_t

// Athena++ headers
#include "../athena.hpp"
#include "task_list.hpp"

// forward declarations
class Mesh;
class MeshBlock;

//----------------------------------------------------------------------------------------
//! \class FFTGravitySolverTaskList
//  \brief data and function definitions for FFTGravitySolverTaskList derived class

class FFTGravitySolverTaskList : public TaskList {
 public:
  FFTGravitySolverTaskList(ParameterInput *pin, Mesh *pm);
  ~FFTGravitySolverTaskList() {}

  // functions
  TaskStatus ClearFFTGravityBoundary(MeshBlock *pmb, int stage);
  TaskStatus SendFFTGravityBoundary(MeshBlock *pmb, int stage);
  TaskStatus ReceiveFFTGravityBoundary(MeshBlock *pmb, int stage);
  TaskStatus SetFFTGravityBoundary(MeshBlock *pmb, int stage);
  TaskStatus PhysicalBoundary(MeshBlock *pmb, int stage);

 private:
  void AddTask(std::uint64_t id, std::uint64_t dep) override;
  void StartupTaskList(MeshBlock *pmb, int stage) override;
};


//----------------------------------------------------------------------------------------
// 64-bit integers with "1" in different bit positions used to ID  each hydro task.
namespace FFTGravitySolverTaskNames {
const std::uint64_t NONE=0ULL;
const std::uint64_t CLEAR_GRAV=1ULL<<0;

const std::uint64_t SEND_GRAV_BND=1ULL<<1;
const std::uint64_t RECV_GRAV_BND=1ULL<<2;
const std::uint64_t SETB_GRAV_BND=1ULL<<3;
const std::uint64_t GRAV_PHYS_BND=1ULL<<4;
}
#endif // TASK_LIST_FFT_GRAV_TASK_LIST_HPP_
