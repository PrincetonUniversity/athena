#ifndef FIELD_HPP
#define FIELD_HPP
//======================================================================================
// Athena++ astrophysical MHD code
// Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
// See LICENSE file for full public license information.
//======================================================================================
//! \file field.hpp
//  \brief defines Field class which implements data and functions for E/B fields
//======================================================================================

// Athena headers
#include "../athena.hpp"         // Real
#include "../athena_arrays.hpp"  // AthenaArray

typedef struct InterfaceField {
  AthenaArray<Real> x1,x2,x3;
} InterfaceField;

typedef struct EdgeField {
  AthenaArray<Real> x1,x2,x3;
} EdgeField;

class MeshBlock;
class ParameterInput;
class Fluid;
//class BFieldIntegrator;

//! \class Field
//  \brief electric and magnetic field data and functions

class Field {
friend class Fluid;
public:
  Field(MeshBlock *pmb, ParameterInput *pin);
  ~Field();

  InterfaceField bi;     // interface magnetic fields
  InterfaceField bi1;    // interface magnetic fields at intermediate step
  AthenaArray<Real> bc;  // cell-centered fields
  AthenaArray<Real> bc1; // cell-centered fields at intermediate step

  InterfaceField ei;     // interface electric fields
  EdgeField emf;         // edge electric fields used to update B using CT

//  BFieldIntegrator *pb_integrator;   // integration algorithm (CT)

private:
  MeshBlock* pmy_mblock_;    // ptr to MeshBlock containing this Field
};
#endif
