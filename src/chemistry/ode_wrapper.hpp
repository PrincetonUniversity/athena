#ifndef ODE_WRAPPER_HPP
#define ODE_WRAPPER_HPP
//======================================================================================
// Athena++ astrophysical MHD code
// Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
// See LICENSE file for full public license information.
//======================================================================================
//! \file ode_wrapper.hpp
//  \brief definitions for ode solver classes.
//======================================================================================

// Athena++ classes headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "network/network.hpp" 
//CVODE headers
#include <sundials/sundials_types.h> /* realtype type*/
#include <sundials/sundials_dense.h>
#include <nvector/nvector_serial.h> /* N_Vector type*/
#include <cvode/cvode.h>            /* CVODE solver fcts., consts. */
#include <cvode/cvode_direct.h>       /* prototype for CVDense */
#include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNMatrix            */
#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver      */

class ParameterInput;
class PassiveScalars;

//! \class ODEWrapper
//  \brief Wrapper for ODE solver, CVODE
class ODEWrapper {
public:
  ODEWrapper(MeshBlock *pmb, ParameterInput *pin);
  ~ODEWrapper();
  //Initialize CVODE
  void Initialize(ParameterInput *pin);
  //Update abundance in PassiveScalars over time dt.
  // For each cell:
  // Step 1: Set the radiation field strength in ChemNetwork.
  // Depends on the data structure of radiation field, this can be copying
  // the value from Radiation class to ChemNetwork class, or just pass a pointer.
  //
  // Step 2: re-initialize CVODE with starting time t, and starting abundance
  // y. If x(k, j, i, ispec), we can just pass a pointer to CVODE, otherwise,
  // we need to copy the abundance of PassiveScalars to an array.
  //
  // Step 3: Integration. Update the array of PassiveScalars abundance in that
  // cell over time dt.
  // 
  // Note that this will be not vectorizable(?).
  void Integrate(const Real tinit, const Real dt);

  //solve the chemical abundance to equilibrium. Useful for post-processing.
  void SolveEq();

  void SetInitStep(const Real h_init);
  //Get the last step size
  Real GetLastStep() const;
  //Get the next step size
  Real GetNextStep() const;
  //Get the number of steps between two reinits.
  long int GetNsteps() const;

private:
  PassiveScalars *pmy_spec_;
  MeshBlock *pmy_block_;
  int dim_; //dimenstion
  Real reltol_;//relative tolerance
  AthenaArray<Real> abstol_;
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


#endif // ODE_WRAPPER_HPP
