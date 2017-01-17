#ifndef GRAVITY_HPP
#define GRAVITY_HPP

//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file gravity.hpp
//  \brief defines Gravity class which implements data and functions for gravitational potential

// Athena++ classes headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../parameter_input.hpp"

class MeshBlock;
class ParameterInput;

//! \class Gravity
//  \brief gravitational potential data and functions

class Gravity {
public:
  Gravity(MeshBlock *pmb, ParameterInput *pin);
  ~Gravity();

  MeshBlock* pmy_block;  // ptr to MeshBlock containing this Field

  AthenaArray<Real> phi, phi_old;  // gravitational potential
  Real gconst, four_pi_gconst;
  Real grav_mean_rho;

  void Initialize(ParameterInput *pin);
  void Solver(const AthenaArray<Real> &den);

private:
  // scratch space used to compute fluxes
  AthenaArray<Real> den_;
};
#endif // GRAVITY_HPP
