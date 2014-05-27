#ifndef BOUNDARY_CONDITIONS_HPP
#define BOUNDARY_CONDITIONS_HPP
//======================================================================================
/* Athena++ astrophysical MHD code
 * Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
 * See LICENSE file for full public license information.
 *====================================================================================*/
/*! \file bvals.hpp
 *  \brief defines class FluidBoundaryConditions
 *  Contains data structures and functions related to BCs for the fluid
 *====================================================================================*/

class ParameterInput;
class Fluid;

//! \class FluidBoundaryConditions
//  \brief BCs data and functions for fluid

class FluidBoundaryConditions {
public:
  FluidBoundaryConditions(Block *pb);
  ~FluidBoundaryConditions();

  void ApplyBoundaryConditions(AthenaArray<Real> &a);

private:
  Block *pparent_block;  // ptr to parent Block

// function pointers, set in Init function based on parameters in input file

  void (*FluidInnerX1_) (Block *pb, AthenaArray<Real> &a);
  void (*FluidOuterX1_) (Block *pb, AthenaArray<Real> &a);
  void (*FluidInnerX2_) (Block *pb, AthenaArray<Real> &a);
  void (*FluidOuterX2_) (Block *pb, AthenaArray<Real> &a);
  void (*FluidInnerX3_) (Block *pb, AthenaArray<Real> &a);
  void (*FluidOuterX3_) (Block *pb, AthenaArray<Real> &a);

};
#endif
