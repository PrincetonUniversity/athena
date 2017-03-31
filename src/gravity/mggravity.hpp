#ifndef MGGRAVITY_HPP
#define MGGRAVITY_HPP

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
#include "../multigrid/multigrid.hpp"

class MeshBlock;
class ParameterInput;
class Coordinates;
class Multigrid

//! \class MGGravity
//  \brief Multigrid gravity Solver

class MGGravity : public Multigrid {
public:
  MGGravity(MeshBlock *pmb) : Multigrid(pmb, 1) {};
  ~MGGravity() {};
  void Smooth(int lev, int color);
  void CalculateResidual(int lev);

private:
  const Real omega = 1.15;
};
#endif // MGGRAVITY_HPP
