#ifndef COORDINATES_HPP
#define COORDINATES_HPP
//======================================================================================
/* Athena++ astrophysical MHD code
 * Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
 * See LICENSE file for full public license information.
 *====================================================================================*/
/*! \file coordinates.hpp
 *  \brief defines Coordinates class used to compute/store geometrical factors (areas,
 *  volumes, source terms) related to a Mesh
 *====================================================================================*/

// Athena headers
#include "../athena.hpp"  // macros, Real
#include "../athena_arrays.hpp"  // AthenaArray

// Declarations
class MeshBlock;
class ParameterInput;

//! \class Coordinates
//  \brief coordinate data and functions

class Coordinates {
public:
  Coordinates(MeshBlock *pmb, ParameterInput *pin);
  ~Coordinates();

  MeshBlock *pmy_block;  // ptr to MeshBlock containing this Coordinates

  void Area1Face(
    const int k, const int j, const int il, const int iu, AthenaArray<Real> *parea);
  void Area2Face(
    const int k, const int j, const int il, const int iu, AthenaArray<Real> *parea);
  void Area3Face(
    const int k, const int j, const int il, const int iu, AthenaArray<Real> *parea);
  void CellVolume(
    const int k, const int j, const int il, const int iu, AthenaArray<Real> *pvol);
  void CoordinateSourceTerms(
    const int k, const int j, AthenaArray<Real> &prim, AthenaArray<Real> &src);

  void CellMetric(const int k, const int j, AthenaArray<Real> &g,
      AthenaArray<Real> &g_inv);
  void PrimToLocal1(const int k, const int j, AthenaArray<Real> &prim);
  void PrimToLocal2(const int k, const int j, AthenaArray<Real> &prim);
  void PrimToLocal3(const int k, const int j, AthenaArray<Real> &prim);
  void FluxToGlobal1(const int k, const int j, AthenaArray<Real> &flux);
  void FluxToGlobal2(const int k, const int j, AthenaArray<Real> &flux);
  void FluxToGlobal3(const int k, const int j, AthenaArray<Real> &flux);
  void PrimToCons(AthenaArray<Real> &prim, AthenaArray<Real> &cons);

// these are scratch arrays used by integrators and allocated in this class
  AthenaArray<Real> face_area, cell_volume;

private:
  AthenaArray<Real> face1_area_i_, face1_area_j_;
  AthenaArray<Real> face2_area_i_, face2_area_j_;
  AthenaArray<Real> face3_area_i_, face3_area_j_;
  AthenaArray<Real> src_terms_i_,  src_terms_j_;
  AthenaArray<Real> src_terms_i1_, src_terms_i2_, src_terms_i3_, src_terms_i4_;
  AthenaArray<Real> src_terms_j1_, src_terms_j2_, src_terms_j3_;
  AthenaArray<Real> volume_i_,     volume_j_;
  AthenaArray<Real> metric_cell_i1_, metric_cell_i2_, metric_cell_i3_, metric_cell_i4_,
      metric_cell_i5_, metric_cell_i6_;
  AthenaArray<Real> metric_cell_j1_, metric_cell_j2_;
  AthenaArray<Real> metric_face1_i1_, metric_face1_i2_, metric_face1_i3_;
  AthenaArray<Real> metric_face1_j1_;
  AthenaArray<Real> metric_face2_i1_, metric_face2_i2_, metric_face2_i3_;
  AthenaArray<Real> metric_face2_j1_;
  AthenaArray<Real> metric_face3_i1_, metric_face3_i2_, metric_face3_i3_;
  AthenaArray<Real> metric_face3_j1_;
  AthenaArray<Real> trans_face1_i1_, trans_face1_i2_, trans_face1_i3_, trans_face1_i4_;
  AthenaArray<Real> trans_face1_j1_, trans_face1_j2_;
  AthenaArray<Real> trans_face2_i1_, trans_face2_i2_, trans_face2_i3_, trans_face2_i4_;
  AthenaArray<Real> trans_face2_j1_, trans_face2_j2_;
  AthenaArray<Real> trans_face3_i1_, trans_face3_i2_, trans_face3_i3_, trans_face3_i4_;
  AthenaArray<Real> trans_face3_j1_, trans_face3_j2_;
};
#endif
