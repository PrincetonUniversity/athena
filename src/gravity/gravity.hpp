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
#include "../coordinates/coordinates.hpp"
#include "../task_list/task_list.hpp"

class MeshBlock;
class ParameterInput;
class Hydro;

//! \class Gravity
//  \brief gravitational potential data and functions

class Gravity {
friend class Hydro;
public:
  Gravity(MeshBlock *pmb, ParameterInput *pin);
  ~Gravity();

  MeshBlock* pmy_block;  // ptr to MeshBlock containing this Field

  AthenaArray<Real> phi;  // gravitational potential

  void GravitySolver(AthenaArray<Real> &den, AthenaArray<Real> &phi, Coordinates *pco, 
       int is, int ie, int js, int je, int ks, int ke);

private:
  // scratch space used to compute fluxes
  AthenaArray<Real> den_;
};
#endif // GRAVITY_HPP
