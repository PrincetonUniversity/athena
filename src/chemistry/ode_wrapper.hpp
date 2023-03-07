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

 private:
  PassiveScalars *pmy_spec_;
  MeshBlock *pmy_block_;
  int dim_; //dimension  of the ODEs
  bool output_zone_sec_; //option to output solver performance
};


#endif // CHEMISTRY_ODE_WRAPPER_HPP_
