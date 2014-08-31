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

// Athena headers
#include "../../athena.hpp"         // Real
#include "../../athena_arrays.hpp"  // AthenaArray

// Declarations
class MeshBlock;
class Fluid;

//! \class FluidBoundaryConditions
//  \brief BCs data and functions for fluid

class FluidBoundaryConditions {
public:
  FluidBoundaryConditions(Fluid *pf);
  ~FluidBoundaryConditions();

  void ApplyBoundaryConditions(AthenaArray<Real> &a);

private:
  Fluid *pmy_fluid;  // ptr to Fluid containing this FluidBoundaryConditions

// function pointers, set in Init function based on parameters in input file

  void (*FluidInnerX1_) (MeshBlock *pmb, AthenaArray<Real> &a);
  void (*FluidOuterX1_) (MeshBlock *pmb, AthenaArray<Real> &a);
  void (*FluidInnerX2_) (MeshBlock *pmb, AthenaArray<Real> &a);
  void (*FluidOuterX2_) (MeshBlock *pmb, AthenaArray<Real> &a);
  void (*FluidInnerX3_) (MeshBlock *pmb, AthenaArray<Real> &a);
  void (*FluidOuterX3_) (MeshBlock *pmb, AthenaArray<Real> &a);

};

//-------------------- prototypes for all BC functions ---------------------------------

  void ReflectInnerX1(MeshBlock *pmb, AthenaArray<Real> &a);
  void ReflectOuterX1(MeshBlock *pmb, AthenaArray<Real> &a);
  void ReflectInnerX2(MeshBlock *pmb, AthenaArray<Real> &a);
  void ReflectOuterX2(MeshBlock *pmb, AthenaArray<Real> &a);
  void ReflectInnerX3(MeshBlock *pmb, AthenaArray<Real> &a);
  void ReflectOuterX3(MeshBlock *pmb, AthenaArray<Real> &a);

  void OutflowInnerX1(MeshBlock *pmb, AthenaArray<Real> &a);
  void OutflowOuterX1(MeshBlock *pmb, AthenaArray<Real> &a);
  void OutflowInnerX2(MeshBlock *pmb, AthenaArray<Real> &a);
  void OutflowOuterX2(MeshBlock *pmb, AthenaArray<Real> &a);
  void OutflowInnerX3(MeshBlock *pmb, AthenaArray<Real> &a);
  void OutflowOuterX3(MeshBlock *pmb, AthenaArray<Real> &a);

  void PeriodicInnerX1(MeshBlock *pmb, AthenaArray<Real> &a);
  void PeriodicOuterX1(MeshBlock *pmb, AthenaArray<Real> &a);
  void PeriodicInnerX2(MeshBlock *pmb, AthenaArray<Real> &a);
  void PeriodicOuterX2(MeshBlock *pmb, AthenaArray<Real> &a);
  void PeriodicInnerX3(MeshBlock *pmb, AthenaArray<Real> &a);
  void PeriodicOuterX3(MeshBlock *pmb, AthenaArray<Real> &a);
#endif
