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

typedef struct InterfaceBField {
  AthenaArray<Real> x1,x2,x3;
} InterfaceBField;

class MeshBlock;
class ParameterInput;
//class BFieldIntegrator;

//! \class Field
//  \brief electric and magnetic field data and functions

class Field {
friend class FluidIntegrator;
public:
  Field(MeshBlock *pmb, ParameterInput *pin);
  ~Field();

  InterfaceBField bi;    // interface magnetic fields
  InterfaceBField bi1;   // interface magnetic fields at intermediate step
  AthenaArray<Real> bc;  // cell-centered fields
  AthenaArray<Real> bc1; // cell-centered fields at intermediate step

//  BFieldIntegrator *pb_integrator;   // integration algorithm (CT)

private:
  MeshBlock* pmy_mblock_;    // ptr to MeshBlock containing this Field
};
#endif
