#ifndef BOUNDARY_CONDITIONS_HPP
#define BOUNDARY_CONDITIONS_HPP
//======================================================================================
/* Athena++ astrophysical MHD code
 * Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
 * See LICENSE file for full public license information.
 *====================================================================================*/
/*! \file boundary_conditions.hpp
 *  \brief defines class BoundaryConditions
 *  Contains data structures and functions related to BCs for the fluid
 *====================================================================================*/

class ParameterInput;
class Fluid;

//! \class BoundaryConditions
//  \brief BCs data and functions

class BoundaryConditions {
public:
  BoundaryConditions(ParameterInput *pin, Fluid *pf);
  ~BoundaryConditions();

  void SetFluidBoundaryValues(Fluid *pf);

private:
  Fluid *pmy_fluid;

// function pointers, set in constructor based on parameters in input file

  void (*FluidInnerX1_) (Fluid *pmy_fluid);
  void (*FluidOuterX1_) (Fluid *pmy_fluid);
  void (*FluidInnerX2_) (Fluid *pmy_fluid);
  void (*FluidOuterX2_) (Fluid *pmy_fluid);
  void (*FluidInnerX3_) (Fluid *pmy_fluid);
  void (*FluidOuterX3_) (Fluid *pmy_fluid);

};
#endif
