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
//  \brief Multigrid gravity solver for each block

class MGGravity : public Multigrid {
public:
  MGGravity(int nx, int ny, int nz, int ngh, Real dx)
    : Multigrid(1, nx, ny, nz, ngh, dx) {};
  ~MGGravity() {};
  void Smooth(int lev, int color);
  void CalculateResidual(int lev);

private:
  const Real omega_ = 1.15;
};


//! \class GravityDriver
//  \brief Multigrid gravity solver

class GravityDriver {
public:
  GravityDriver();
  ~GravityDriver();

private:
  MGGravity *mgroot_;
  int umode_, mode_;
};

#endif // MGGRAVITY_HPP
