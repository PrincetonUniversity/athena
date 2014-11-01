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
struct InterfaceBField;

//-------------------- prototypes for all BC functions ---------------------------------

void ReflectInnerX1(MeshBlock *pmb, AthenaArray<Real> &a);
void ReflectInnerX2(MeshBlock *pmb, AthenaArray<Real> &a);
void ReflectInnerX3(MeshBlock *pmb, AthenaArray<Real> &a);
void ReflectOuterX1(MeshBlock *pmb, AthenaArray<Real> &a);
void ReflectOuterX2(MeshBlock *pmb, AthenaArray<Real> &a);
void ReflectOuterX3(MeshBlock *pmb, AthenaArray<Real> &a);

void ReflectInnerX1(MeshBlock *pmb, InterfaceBField &a);
void ReflectInnerX2(MeshBlock *pmb, InterfaceBField &a);
void ReflectInnerX3(MeshBlock *pmb, InterfaceBField &a);
void ReflectOuterX1(MeshBlock *pmb, InterfaceBField &a);
void ReflectOuterX2(MeshBlock *pmb, InterfaceBField &a);
void ReflectOuterX3(MeshBlock *pmb, InterfaceBField &a);

void OutflowInnerX1(MeshBlock *pmb, AthenaArray<Real> &a);
void OutflowInnerX2(MeshBlock *pmb, AthenaArray<Real> &a);
void OutflowInnerX3(MeshBlock *pmb, AthenaArray<Real> &a);
void OutflowOuterX1(MeshBlock *pmb, AthenaArray<Real> &a);
void OutflowOuterX2(MeshBlock *pmb, AthenaArray<Real> &a);
void OutflowOuterX3(MeshBlock *pmb, AthenaArray<Real> &a);

void OutflowInnerX1(MeshBlock *pmb, InterfaceBField &a);
void OutflowInnerX2(MeshBlock *pmb, InterfaceBField &a);
void OutflowInnerX3(MeshBlock *pmb, InterfaceBField &a);
void OutflowOuterX1(MeshBlock *pmb, InterfaceBField &a);
void OutflowOuterX2(MeshBlock *pmb, InterfaceBField &a);
void OutflowOuterX3(MeshBlock *pmb, InterfaceBField &a);

void PeriodicInnerX1(MeshBlock *pmb, AthenaArray<Real> &a);
void PeriodicInnerX2(MeshBlock *pmb, AthenaArray<Real> &a);
void PeriodicInnerX3(MeshBlock *pmb, AthenaArray<Real> &a);
void PeriodicOuterX1(MeshBlock *pmb, AthenaArray<Real> &a);
void PeriodicOuterX2(MeshBlock *pmb, AthenaArray<Real> &a);
void PeriodicOuterX3(MeshBlock *pmb, AthenaArray<Real> &a);

void PeriodicInnerX1(MeshBlock *pmb, InterfaceBField &a);
void PeriodicInnerX2(MeshBlock *pmb, InterfaceBField &a);
void PeriodicInnerX3(MeshBlock *pmb, InterfaceBField &a);
void PeriodicOuterX1(MeshBlock *pmb, InterfaceBField &a);
void PeriodicOuterX2(MeshBlock *pmb, InterfaceBField &a);
void PeriodicOuterX3(MeshBlock *pmb, InterfaceBField &a);

enum EdgeNames {inner_x1, outer_x1, inner_x2, outer_x2, inner_x3, outer_x3};
typedef void (*BValFluid_t)(MeshBlock *pmb, AthenaArray<Real> &a);
typedef void (*BValBField_t)(MeshBlock *pmb, InterfaceBField &a);

//! \class BoundaryValues
//  \brief BVals data and functions

class BoundaryValues {
public:
  BoundaryValues(MeshBlock *pmb, ParameterInput *pin);
  ~BoundaryValues();

  template<typename T> void ApplyBVals(T &input);
  void EnrollFluidBValFunction (enum EdgeNames edge, BValFluid_t  my_bc);
  void EnrollBFieldBValFunction(enum EdgeNames edge, BValBField_t my_bc);

private:
  MeshBlock *pmy_mblock_;  // ptr to MeshBlock containing this BVals

  template<typename T> void BValsInnerX1_(MeshBlock* pmb, T &input);
  template<typename T> void BValsInnerX2_(MeshBlock* pmb, T &input);
  template<typename T> void BValsInnerX3_(MeshBlock* pmb, T &input);
  template<typename T> void BValsOuterX1_(MeshBlock* pmb, T &input);
  template<typename T> void BValsOuterX2_(MeshBlock* pmb, T &input);
  template<typename T> void BValsOuterX3_(MeshBlock* pmb, T &input);

  BValFluid_t FluidInnerX1_, FluidInnerX2_, FluidInnerX3_;
  BValFluid_t FluidOuterX1_, FluidOuterX2_, FluidOuterX3_;
  BValBField_t BFieldInnerX1_, BFieldInnerX2_, BFieldInnerX3_;
  BValBField_t BFieldOuterX1_, BFieldOuterX2_, BFieldOuterX3_;
};
#endif
