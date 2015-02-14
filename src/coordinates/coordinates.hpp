#ifndef COORDINATES_HPP
#define COORDINATES_HPP
//======================================================================================
// Athena++ astrophysical MHD code
// Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
// See LICENSE file for full public license information.
//======================================================================================
//! \file coordinates.hpp
//  \brief defines Coordinates class used to compute/store geometrical factors (areas,
//  volumes, source terms) related to a Mesh
//======================================================================================

// Athena headers
#include "../athena.hpp"         // macros, Real
#include "../athena_arrays.hpp"  // AthenaArray

// forward declarations
class MeshBlock;
class ParameterInput;

//! \class Coordinates
//  \brief coordinate data and functions

class Coordinates {
public:
  friend class FluidSourceTerms;
  Coordinates(MeshBlock *pmb, ParameterInput *pin);
  ~Coordinates();

  MeshBlock *pmy_block;  // ptr to MeshBlock containing this Coordinates

// functions to compute length of edges
  void Edge1Length(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &len);
  void Edge2Length(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &len);
  void Edge3Length(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &len);

// functions to compute physical width at cell center
  Real CenterWidth1(const int k, const int j, const int i);
  Real CenterWidth2(const int k, const int j, const int i);
  Real CenterWidth3(const int k, const int j, const int i);

// functions to compute area of faces
  void Face1Area(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &area);
  void Face2Area(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &area);
  void Face3Area(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &area);

// function to compute volume of cells
  void CellVolume(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &vol);

// function to compute geometrical source terms
  void CoordSrcTermsX1(const int k, const int j, const Real dt,
    const AthenaArray<Real> &flx,
    const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc, AthenaArray<Real> &u);

  void CoordSrcTermsX2(const int k, const int j, const Real dt,
    const AthenaArray<Real> &flx,  const AthenaArray<Real> &flx_m1,
    const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc, AthenaArray<Real> &u);

  void CoordSrcTermsX3(const int k, const int j, const Real dt,
    const AthenaArray<Real> &flx,  const AthenaArray<Real> &flx_m1,
    const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc, AthenaArray<Real> &u);

  // Functions for use in general relativity
  void CellMetric(const int k,const int j, AthenaArray<Real> &g, AthenaArray<Real> &gi);
  void Face1Metric(const int k, const int j, AthenaArray<Real> &g,
      AthenaArray<Real> &g_inv);
  void Face2Metric(const int k, const int j, AthenaArray<Real> &g,
      AthenaArray<Real> &g_inv);
  void Face3Metric(const int k, const int j, AthenaArray<Real> &g,
      AthenaArray<Real> &g_inv);
  void PrimToLocal1(const int k, const int j, const AthenaArray<Real> &b1_vals,
      AthenaArray<Real> &prim_left, AthenaArray<Real> &prim_right,
      AthenaArray<Real> &bx);
  void PrimToLocal2(const int k, const int j, const AthenaArray<Real> &b2_vals,
      AthenaArray<Real> &prim_left, AthenaArray<Real> &prim_right,
      AthenaArray<Real> &by);
  void PrimToLocal3(const int k, const int j, const AthenaArray<Real> &b3_vals,
      AthenaArray<Real> &prim_left, AthenaArray<Real> &prim_right,
      AthenaArray<Real> &bz);
  void FluxToGlobal1(const int k, const int j, AthenaArray<Real> &flux);
  void FluxToGlobal2(const int k, const int j, AthenaArray<Real> &flux);
  void FluxToGlobal3(const int k, const int j, AthenaArray<Real> &flux);
  void PrimToCons(
      const AthenaArray<Real> &prim, const AthenaArray<Real> &b, Real gamma_adi_red,
      AthenaArray<Real> &cons);
  Real DistanceBetweenPoints(Real a1, Real a2, Real a3, Real bx, Real by, Real bz);
  void TransformVectorCell(Real at, Real ax, Real ay, Real az, int k, int j, int i,
      Real *a0, Real *a1, Real *a2, Real *a3);
  void TransformVectorFace1(Real at, Real ax, Real ay, Real az, int k, int j, int i,
      Real *a0, Real *a1, Real *a2, Real *a3);
  void TransformVectorFace2(Real at, Real ax, Real ay, Real az, int k, int j, int i,
      Real *a0, Real *a1, Real *a2, Real *a3);
  void TransformVectorFace3(Real at, Real ax, Real ay, Real az, int k, int j, int i,
      Real *a0, Real *a1, Real *a2, Real *a3);

private:
// scratch arrays containing precomputed factors used by functions in this class
  AthenaArray<Real> coord_area1_i_, coord_area2_i_, coord_area3_i_;
  AthenaArray<Real> coord_area1_j_, coord_area2_j_, coord_area3_j_;
  AthenaArray<Real> coord_vol_i_,   coord_vol_j_;
  AthenaArray<Real> coord_src1_i_,  coord_src2_i_;
  AthenaArray<Real> coord_src1_j_,  coord_src2_j_;
  AthenaArray<Real> edge1_length_i_, edge1_length_j_;
  AthenaArray<Real> edge2_length_i_, edge2_length_j_;
  AthenaArray<Real> edge3_length_i_, edge3_length_j_;
  AthenaArray<Real> cell_width1_i_;
  AthenaArray<Real> src_terms_i1_, src_terms_i2_, src_terms_i3_, src_terms_i4_;
  AthenaArray<Real> src_terms_j1_, src_terms_j2_, src_terms_j3_;
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
