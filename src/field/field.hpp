#ifndef FIELD_FIELD_HPP_
#define FIELD_FIELD_HPP_
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file field.hpp
//  \brief defines Field class which implements data and functions for E/B fields

// Athena++ classes headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../task_list/task_list.hpp"

class MeshBlock;
class ParameterInput;
class Hydro;
class FieldDiffusion;

//! \class Field
//  \brief electric and magnetic field data and functions

class Field {
friend class Hydro;
public:
  Field(MeshBlock *pmb, ParameterInput *pin);
  ~Field();

  MeshBlock* pmy_block;  // ptr to MeshBlock containing this Field
  FieldDiffusion *pfdif;

  // face-centered magnetic fields
  FaceField b;       // time-integrator memory register #1
  FaceField b1;      // time-integrator memory register #2
  FaceField b2;      // time-integrator memory register #3
  // (no more than MAX_NREGISTER allowed)

  // cell-centered magnetic fields
  AthenaArray<Real> bcc;  // time-integrator memory register #1

  EdgeField e;    // edge-centered electric fields used in CT
  FaceField wght; // weights used to integrate E to corner using GS algorithm
  AthenaArray<Real> e2_x1f, e3_x1f; // electric fields at x1-face from Riemann solver
  AthenaArray<Real> e1_x2f, e3_x2f; // electric fields at x2-face from Riemann solver
  AthenaArray<Real> e1_x3f, e2_x3f; // electric fields at x3-face from Riemann solver

  void CalculateCellCenteredField(const FaceField &bf, AthenaArray<Real> &bc,
       Coordinates *pco, int is, int ie, int js, int je, int ks, int ke);
  void CT(const Real wght, FaceField &b_out);
  void WeightedAveB(FaceField &b_out, FaceField &b_in1, FaceField &b_in2,
       const Real wght[3]);
  void ComputeCornerE(AthenaArray<Real> &w, AthenaArray<Real> &bcc);

private:
  // scratch space used to compute fluxes
  AthenaArray<Real> cc_e_;
  AthenaArray<Real> face_area_, edge_length_, edge_length_p1_;
  AthenaArray<Real> g_, gi_;  // only used in GR
};
#endif // FIELD_FIELD_HPP_
