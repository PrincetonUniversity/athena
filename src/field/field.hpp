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
  void CalculateCellCenteredField(const FaceField &bf, AthenaArray<Real> &bc,
       Coordinates *pco, int is, int ie, int js, int je, int ks, int ke);

  MeshBlock* pmy_mblock;  // ptr to MeshBlock containing this Field

  FaceField b;       // face-centered magnetic fields
  FaceField b1;      // face-centered magnetic fields at intermediate step
  AthenaArray<Real> bcc;  // cell-centered magnetic fields
  AthenaArray<Real> bcc1; // cell-centered magnetic fields at intermediate step

  EdgeField e;         // edge-centered electric fields used in CT
  FaceField ei;   // face-centered electric fields (e.g. from Riemann solver)
  FaceField wght; // weights used to integrate E to corner using GS algorithm

  FieldIntegrator *pintegrator;  // integration algorithm (CT)

  void CopyOrAverageField(FaceField &a, FaceField &b, FaceField &c, Real factor);

private:
};
#endif // FIELD_HPP
