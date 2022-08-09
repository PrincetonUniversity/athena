#ifndef CHEMISTRY_ODE_WRAPPER_HPP_
#define CHEMISTRY_ODE_WRAPPER_HPP_
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file ode_wrapper.hpp
//! \brief definitions for ode solver classes.

// C headers
// CVODE headers
#include <cvode/cvode.h>            // CVODE solver fcts., consts.
#include <cvode/cvode_direct.h>     // prototype for CVDense
#include <nvector/nvector_serial.h> // N_Vector type
#include <sundials/sundials_dense.h>
#include <sunlinsol/sunlinsol_dense.h> // access to dense SUNLinearSolver
#include <sunmatrix/sunmatrix_dense.h> // access to dense SUNMatrix

// C++ headers

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "network/network.hpp"

class ParameterInput;
class PassiveScalars;

//! \class ODEWrapper
//! \brief Wrapper for ODE solver, CVODE
class ODEWrapper {
 public:
  ODEWrapper(MeshBlock *pmb, ParameterInput *pin);
  ~ODEWrapper();
  void Initialize(ParameterInput *pin);
  void Integrate(const Real tinit, const Real dt);

  void SetInitStep(const Real h_init);
  Real GetLastStep() const;
  Real GetNextStep() const;
  long int GetNsteps() const; // NOLINT (runtime/int)

 private:
  PassiveScalars *pmy_spec_;
  MeshBlock *pmy_block_;
  int dim_; //dimenstion
  Real reltol_;//relative tolerance
  AthenaArray<Real> abstol_;
  SUNContext sunctx_;
  SUNMatrix dense_matrix_;
  SUNLinearSolver dense_ls_;
  void *cvode_mem_;
  N_Vector y_;
  Real *ydata_;
  Real h_init_;
  Real fac_dtmax_;//factor of the max timestep in CVODE relative to the hydrostep
  int output_zone_sec_;

  //CVODE checkflag
  void CheckFlag(const void *flagvalue, const char *funcname,
                 const int opt) const;
};


#endif // CHEMISTRY_ODE_WRAPPER_HPP_
