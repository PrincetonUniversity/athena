#ifndef FIELD_INTEGRATOR_HPP
#define FIELD_INTEGRATOR_HPP
//======================================================================================
// Athena++ astrophysical MHD code
// Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
// See LICENSE file for full public license information.
//======================================================================================
//! \file field_integrator.hpp
//  \brief defines class FieldIntegrator for evolving magnetic field
//======================================================================================

// Athena headers
#include "../../athena.hpp"         // Real
#include "../../athena_arrays.hpp"  // AthenaArray

// Forward declarations
class Field;
class ParameterInput;
class MeshBlock;

//! \class FieldIntegrator
//  \brief member functions implement various integration algorithms for the field

class FieldIntegrator {
public:
  FieldIntegrator(Field *pfd, ParameterInput *pin);
  ~FieldIntegrator();

  Field *pmy_field;  // ptr to Field containing this FieldIntegrator

  void CT(MeshBlock *pmb, InterfaceField &b, AthenaArray<Real> &w,
    AthenaArray<Real> &bcc, const int step);
  void ComputeCornerE(MeshBlock *pmb, AthenaArray<Real> &w, AthenaArray<Real> &bcc);
  void BoundaryValuesFaceCenteredE1(MeshBlock *pmb);
  void BoundaryValuesFaceCenteredE2(MeshBlock *pmb);
  void BoundaryValuesFaceCenteredE3(MeshBlock *pmb);

private:
// scratch space used in integrator
  AthenaArray<Real> cc_e1_, cc_e2_, cc_e3_;
  AthenaArray<Real> face_area_, edge_length_, edge_lengthp1_;
};
#endif
