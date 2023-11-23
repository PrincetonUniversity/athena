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
class MGGravityDriver;

//! \class Gravity
//! \brief gravitational potential data and functions

class Gravity {
 public:
  Gravity(MeshBlock *pmb, ParameterInput *pin);
  ~Gravity();

  MeshBlock* pmy_block;  // ptr to MeshBlock containing this Field
  AthenaArray<Real> phi, coarse_phi;   // gravitational potential
  AthenaArray<Real> def;   // defect from the Multigrid solver
  AthenaArray<Real> empty_flux[3];
  Real four_pi_G;
  bool output_defect;
  bool fill_ghost;

  // TODO(felker): consider creating a CellCentered.. derived class, and changing to
  //GravityBoundaryVariable *pgbval;
  CellCenteredBoundaryVariable gbvar;

  void SaveFaceBoundaries();
  void RestoreFaceBoundaries();

  friend class MGGravityDriver;

 private:
  MGGravity *pmg;
  AthenaArray<Real> fbuf_[6];
};

#endif // GRAVITY_GRAVITY_HPP_
