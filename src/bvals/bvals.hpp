#ifndef BOUNDARY_CONDITIONS_HPP
#define BOUNDARY_CONDITIONS_HPP
//======================================================================================
/* Athena++ astrophysical MHD code
 * Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
 * See LICENSE file for full public license information.
 *====================================================================================*/
/*! \file boundary_conditions.hpp
 *  \brief defines class FluidBoundaryConditions
 *  Contains data structures and functions related to BCs for the fluid
 *====================================================================================*/

class ParameterInput;
class Fluid;

//! \class BoundaryConditions
//  \brief BCs data and functions

class FluidBoundaryConditions {
public:
  FluidBoundaryConditions(Fluid *pf);
  ~FluidBoundaryConditions();

// functions to initialize BC function pointers, and to apply BCs at each edge of Block
  void InitBoundaryConditions(ParameterInput *pin);
  void ApplyBoundaryConditions(AthenaArray<Real> &a);

private:
  Fluid *pparent_fluid;  // ptr to parent Fluid

// function pointers, set in Init function based on parameters in input file

  void (*FluidInnerX1_) (Fluid *pmy_fluid, AthenaArray<Real> &a);
  void (*FluidOuterX1_) (Fluid *pmy_fluid, AthenaArray<Real> &a);
  void (*FluidInnerX2_) (Fluid *pmy_fluid, AthenaArray<Real> &a);
  void (*FluidOuterX2_) (Fluid *pmy_fluid, AthenaArray<Real> &a);
  void (*FluidInnerX3_) (Fluid *pmy_fluid, AthenaArray<Real> &a);
  void (*FluidOuterX3_) (Fluid *pmy_fluid, AthenaArray<Real> &a);

};
#endif
