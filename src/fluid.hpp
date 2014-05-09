#ifndef FLUID_HPP
#define FLUID_HPP
//======================================================================================
/* Athena++ astrophysical MHD code
 * Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
 * See LICENSE file for full public license information.
 *====================================================================================*/
/*! \file fluid.hpp
 *  \brief defines derived class Fluid, based on base class Mesh
 *  contains data structures and functions for a Fluid stored on the Mesh
 *====================================================================================*/

class ParameterInput;
class BoundaryConditions;

//! \class Fluid
//  \brief fluid data and functions

class Fluid {
public:
  Fluid(ParameterInput *pin, Block *pb);
  ~Fluid();

  AthenaArray<Real> u,w;     // conserved and primitive variables
  AthenaArray<Real> u1_,w1_; // conserved and primitive variables at the half-time step
  Real time, dt;

  Block* pmy_block;          // pointer to parent Block of this Fluid
  Real GetGamma() const { return gamma_; }

  BoundaryConditions *pbvals;
  void Problem(ParameterInput *pin);

private:
  Real gamma_;               // ratio of specific heats

  AthenaArray<Real> wl_,wr_,flx_;

};
#endif
