#ifndef BOUNDARY_VALUES_HPP
#define BOUNDARY_VALUES_HPP
//======================================================================================
/* Athena++ astrophysical MHD code
 * Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
 * See LICENSE file for full public license information.
 *====================================================================================*/
/*! \file bvals.hpp
 *  \brief defines BoundaryValues class used for setting BCs on all data types
 *====================================================================================*/

// Athena headers
#include "../athena.hpp"         // Real
#include "../athena_arrays.hpp"  // AthenaArray

// forward declarations
class MeshBlock;
class Fluid;
class ParameterInput;

enum EdgeNames {inner_x1, outer_x1, inner_x2, outer_x2, inner_x3, outer_x3};
typedef void (*BValFluid_t)(MeshBlock *pmb, AthenaArray<Real> &a);

//! \class BoundaryValues
//  \brief BVals data and functions

class BoundaryValues {
public:
  BoundaryValues(MeshBlock *pmb, ParameterInput *pin);
  ~BoundaryValues();

  void ApplyBVals(AthenaArray<Real> &a);
  void EnrollBoundaryFunction(enum EdgeNames edge, BValFluid_t my_bc);

private:
  MeshBlock *pmy_mblock_;  // ptr to MeshBlock containing this BVals

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
