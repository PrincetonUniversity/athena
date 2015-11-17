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

// Athena++ classes headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp" // Coordinates

//typedef struct InterfaceField {
//  AthenaArray<Real> x1f,x2f,x3f;
//} InterfaceField;

class MeshBlock;
class ParameterInput;
class Hydro;
class FieldIntegrator;

//! \class Field
//  \brief electric and magnetic field data and functions

class Field {
friend class Hydro;
public:
  Field(MeshBlock *pmb, ParameterInput *pin);
  ~Field();
  void CalculateCellCenteredField(const InterfaceField &bf, AthenaArray<Real> &bc,
       Coordinates *pco, int is, int ie, int js, int je, int ks, int ke);

  MeshBlock* pmy_mblock;  // ptr to MeshBlock containing this Field

  InterfaceField b;       // face-centered magnetic fields
  InterfaceField b1;      // face-centered magnetic fields at intermediate step
  AthenaArray<Real> bcc;  // cell-centered magnetic fields
  AthenaArray<Real> bcc1; // cell-centered magnetic fields at intermediate step

  EdgeField e;         // edge-centered electric fields used in CT
  InterfaceField ei;   // face-centered electric fields (e.g. from Riemann solver)
  InterfaceField wght; // weights used to integrate E to corner using GS algorithm

  FieldIntegrator *pint;  // integration algorithm (CT)

private:
};
#endif // FIELD_HPP
