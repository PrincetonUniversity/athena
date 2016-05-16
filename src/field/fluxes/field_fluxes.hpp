#ifndef FIELD_FLUXES_HPP
#define FIELD_FLUXES_HPP
//======================================================================================
// Athena++ astrophysical MHD code
// Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
// See LICENSE file for full public license information.
//======================================================================================
//! \file field_fluxes.hpp
//  \brief defines class FieldFluxes; data and functions for fluxes of B-field
//======================================================================================

// Athena headers
#include "../../athena.hpp"         // Real
#include "../../athena_arrays.hpp"  // AthenaArray

// Forward declarations
class Field;
class ParameterInput;
class MeshBlock;

//! \class FieldFluxes
//  \brief member functions implement various flux functions for the B-field

class FieldFluxes {
public:
  FieldFluxes(Field *pfd, ParameterInput *pin);
  ~FieldFluxes();

  Field *pmy_field;  // ptr to Field containing this FieldFluxes

  void CT(MeshBlock *pmb, FaceField &b, AthenaArray<Real> &w,
    AthenaArray<Real> &bcc, const int step);

  void ComputeCornerE(MeshBlock *pmb, AthenaArray<Real> &w, AthenaArray<Real> &bcc);

private:
  // scratch space used to compute fluxes
  AthenaArray<Real> cc_e_;
  AthenaArray<Real> face_area_, edge_length_, edge_length_p1_;
  AthenaArray<Real> g_, gi_;  // only used in GR
};
#endif
