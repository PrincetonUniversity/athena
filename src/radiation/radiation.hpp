#ifndef RADIATION_RADIATION_HPP_
#define RADIATION_RADIATION_HPP_
//======================================================================================
// Athena++ astrophysical MHD code
// Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
// See LICENSE file for full public license information.
//======================================================================================
//! \file radiation.hpp
//  \brief definitions for Radiation class
//======================================================================================

//c++ header
#include <string>

//Athena++ classes headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../defs.hpp"
#include "../parameter_input.hpp"

class RadIntegrator;

class Radiation {
 public:
  Radiation(MeshBlock *pmb, ParameterInput *pin);

  MeshBlock* pmy_block;    // ptr to MeshBlock containing this Fluid
  RadIntegrator* pradintegrator; //ptr to radiation integrator
  std::string integrator;

  AthenaArray<Real> ir; // radiation specific intensity
  // average radiation intensity over all angles, for output
  AthenaArray<Real> ir_avg;

  int nang, nfreq, n_fre_ang; // n_fre_ang=nang*nfreq
};

#endif // RADIATION_RADIATION_HPP_
