#ifndef BFIELD_HPP
#define BFIELD_HPP
//======================================================================================
/* Athena++ astrophysical MHD code
 * Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
 * See LICENSE file for full public license information.
 *====================================================================================*/
/*! \file bfield.hpp
 *  \brief defines BField class which implements data and functions for magnetic field
 *====================================================================================*/

// Athena headers
#include "../athena.hpp"         // Real
#include "../athena_arrays.hpp"  // AthenaArray

class MeshBlock;
class ParameterInput;
//class BFieldIntegrator;
//class BFieldBoundaryConditions;

//! \class BField
//  \brief magnetic field data and functions

class BField {
friend class BFieldIntegrator;
public:
  BField(MeshBlock *pmb, ParameterInput *pin);
  ~BField();

  AthenaArray<Real> b1i,b2i,b3i;     // interface magnetic fields
  AthenaArray<Real> b1i1,b2i1,b3i1;  // interface magnetic fields at intermediate step

//  BFieldIntegrator *pb_integrator;   // integration algorithm (CT)
//  BFieldBoundaryConditions *pb_bcs;  // boundary conditions

private:
  MeshBlock* pmy_block_;    // ptr to MeshBlock containing this Fluid
};
#endif
