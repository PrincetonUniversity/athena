#ifndef BFIELD_BOUNDARY_CONDITIONS_HPP
#define BFIELD_BOUNDARY_CONDITIONS_HPP
//======================================================================================
/* Athena++ astrophysical MHD code
 * Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
 * See LICENSE file for full public license information.
 *====================================================================================*/
/*! \file bfield_bvals.hpp
 *  \brief defines class FieldBCs
 *  Contains data structures and functions related to BCs for the magnetic field
 *====================================================================================*/

// Athena headers
#include "../../athena.hpp"         // Real
#include "../../athena_arrays.hpp"  // AthenaArray

// forward declarations
class MeshBlock;
class Fluid;
class ParameterInput;

//! \class BFieldBCs
//  \brief BCs data and functions for magnetic field

class BFieldBCs {
public:
  BFieldBCs(BField *pb, ParameterInput *pin);
  ~BFieldBCs();

  void ApplyBFieldBCs(AthenaArray<Real> &a);
  void ApplyEFieldBCs(AthenaArray<Real> &a);

private:
  BField *pmy_bfield;  // ptr to BField containing this BFieldBCs

// function pointers, set in constructor based on parameters in input file

  void (*BFieldInnerX1_) (MeshBlock *pmb, AthenaArray<Real> &a);
  void (*BFieldOuterX1_) (MeshBlock *pmb, AthenaArray<Real> &a);
  void (*BFieldInnerX2_) (MeshBlock *pmb, AthenaArray<Real> &a);
  void (*BFieldOuterX2_) (MeshBlock *pmb, AthenaArray<Real> &a);
  void (*BFieldInnerX3_) (MeshBlock *pmb, AthenaArray<Real> &a);
  void (*BFieldOuterX3_) (MeshBlock *pmb, AthenaArray<Real> &a);

};

//-------------------- prototypes for all BC functions ---------------------------------

  void ReflectBField_ix1(MeshBlock *pmb, AthenaArray<Real> &a);
  void ReflectBField_ox1(MeshBlock *pmb, AthenaArray<Real> &a);
  void ReflectBField_ix2(MeshBlock *pmb, AthenaArray<Real> &a);
  void ReflectBField_ox2(MeshBlock *pmb, AthenaArray<Real> &a);
  void ReflectBField_ix3(MeshBlock *pmb, AthenaArray<Real> &a);
  void ReflectBField_ox3(MeshBlock *pmb, AthenaArray<Real> &a);

  void OutflowBField_ix1(MeshBlock *pmb, AthenaArray<Real> &a);
  void OutflowBField_ox1(MeshBlock *pmb, AthenaArray<Real> &a);
  void OutflowBField_ix2(MeshBlock *pmb, AthenaArray<Real> &a);
  void OutflowBField_ox2(MeshBlock *pmb, AthenaArray<Real> &a);
  void OutflowBField_ix3(MeshBlock *pmb, AthenaArray<Real> &a);
  void OutflowBField_ox3(MeshBlock *pmb, AthenaArray<Real> &a);

  void PeriodicBField_ix1(MeshBlock *pmb, AthenaArray<Real> &a);
  void PeriodicBField_ox1(MeshBlock *pmb, AthenaArray<Real> &a);
  void PeriodicBField_ix2(MeshBlock *pmb, AthenaArray<Real> &a);
  void PeriodicBField_ox2(MeshBlock *pmb, AthenaArray<Real> &a);
  void PeriodicBField_ix3(MeshBlock *pmb, AthenaArray<Real> &a);
  void PeriodicBField_ox3(MeshBlock *pmb, AthenaArray<Real> &a);
#endif
