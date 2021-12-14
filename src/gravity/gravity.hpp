#ifndef GRAVITY_GRAVITY_HPP_
#define GRAVITY_GRAVITY_HPP_
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file gravity.hpp
//! \brief defines Gravity class which implements data and functions for gravitational
//!        potential. Shared by both Multigrid and FFT schemes for self-gravity.

// C headers

// C++ headers
#include <utility>
#include <vector>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../bvals/bvals.hpp"
#include "../bvals/cc/bvals_cc.hpp"

class MeshBlock;
class ParameterInput;
class Coordinates;
class GravityBoundaryValues;
class MGGravity;
class MGGRavityDriver;

//! \class Gravity
//! \brief gravitational potential data and functions

class Gravity {
 public:
  Gravity(MeshBlock *pmb, ParameterInput *pin);

  MeshBlock* pmy_block;  // ptr to MeshBlock containing this Field
  AthenaArray<Real> src;   // shallow copy of accumulated source
  AthenaArray<Real> phi;   // gravitational potential
  AthenaArray<Real> empty_flux[3];
  Real gconst, four_pi_G;
  bool srcterm;

  // TODO(felker): consider creating a CellCentered.. derived class, and changing to
  //GravityBoundaryVariable *pgbval;
  CellCenteredBoundaryVariable gbvar;

  void Initialize(ParameterInput *pin);
  void Solver(const AthenaArray<Real> &u);

  void EnrollSource(AthenaArray<Real> &arr, int idx);
  void UnenrollSource(AthenaArray<Real> &arr, int idx);
  void ComputeSource();

  friend class MGGravityDriver;

 private:
  bool gravity_tensor_momentum_;
  bool gravity_tensor_energy_;

  // actual copy of accumulated source and (pointers, index) pair to enrolled sources
  AthenaArray<Real> accumulated_src_;
  std::vector<std::pair<AthenaArray<Real>*, int>> enrolled_src_;
};

#endif // GRAVITY_GRAVITY_HPP_
