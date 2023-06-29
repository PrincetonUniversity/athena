#ifndef NR_RADIATION_IMPLICIT_RADIATION_IMPLICIT_HPP_
#define NR_RADIATION_IMPLICIT_RADIATION_IMPLICIT_HPP_
//======================================================================================
// Athena++ astrophysical MHD code
// Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
// See LICENSE file for full public license information.
//======================================================================================
//! \file radiation_implicit.hpp
//  \brief implicit radiation class definitions
//======================================================================================



// C++ headers
#include <cstdint>     // int64_t
#include <functional>  // reference_wrapper
#include <string>
#include <vector>

// Athena++ headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../../task_list/im_rad_task_list.hpp"

class Mesh;
class ParameterInput;
class NRRadiation;
class RadIntegrator;
class TimeIntegratorTaskList;

class IMRadiation {
  friend class NRRadiation;
  friend class RadIntegrator;
 public:
  IMRadiation(Mesh *pm, ParameterInput *pin);
  //  ~Radiation();

  void Iteration(Mesh *pm,
             TimeIntegratorTaskList *ptlist, int stage);
  void CheckResidual(MeshBlock *pmb,
        AthenaArray<Real> &ir_old, AthenaArray<Real> &ir_new);

  Real cfl_rad; // the additional CFL number.
                // Small cfl_rad is good for convergence,
                // but with smaller time step

  int ite_scheme;
  int rb_or_not;

  IMRadITTaskList *pimraditlist;
  IMRadHydroTaskList *pimradhylist;
  IMRadComptTaskList *pimradcomptlist;

  // scheduled relaxation jacobi method
  //   // X.Yang, R.Mittal, JCP 274 (2014) 695-708.
  Real omega;
  int srj_p; // <= 9; 0 means off.
  int srj_q[9]; // repetitions
  Real srj_w[9]; // omegas

  // internal data
  int srj_level, srj_cnt; // internal counts

  // The function pointer for srj parameters
  SRJFunc SetSRJParameters;
  void EnrollSRJFunction(SRJFunc MySRJFunction);


 private:
  Real sum_diff_;
  Real sum_full_;
  int nlimit_;       // threadhold for the number of iterations
  Real error_limit_;
};

#endif // NR_RADIATION_IMPLICIT_RADIATION_IMPLICIT_HPP_
