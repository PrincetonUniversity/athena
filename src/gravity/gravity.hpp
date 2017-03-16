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

class MeshBlock;
class ParameterInput;
class Coordinates;

//! \class Gravity
//  \brief gravitational potential data and functions

class Gravity {
public:
  Gravity(MeshBlock *pmb, ParameterInput *pin);
  ~Gravity();

  MeshBlock* pmy_block;  // ptr to MeshBlock containing this Field

  AthenaArray<Real> phi, phi_old;  // gravitational potential
  AthenaArray<Real> gflx[3]; // gravity tensor
  Real gconst, four_pi_G;
  Real grav_mean_rho;

  void AddGravityFlux(AthenaArray<Real> *flx);
  void CalculateGravityFlux();
  void Initialize(ParameterInput *pin);
  void Solver(const AthenaArray<Real> &u);

private:
};
#endif // GRAVITY_HPP
