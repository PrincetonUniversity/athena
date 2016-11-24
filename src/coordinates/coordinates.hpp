#ifndef COORDINATES_HPP
#define COORDINATES_HPP
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file coordinates.hpp
//  \brief defines abstract base and derived classes for coordinates.  These classes
//  provide data and functions to compute/store coordinate positions and spacing, as well
//  as geometrical factors (areas, volumes, coordinate source terms) for various 
//  coordinate systems.

// Athena++ classes headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../mesh/mesh.hpp"
#include <iostream>

// forward declarations
class MeshBlock;
class ParameterInput;

//----------------------------------------------------------------------------------------
//! \class Coordinates
//  \brief abstract base class for all coordinate derived classes

class Coordinates {
public:
//  friend class HydroSourceTerms;
  Coordinates(MeshBlock *pmb, ParameterInput *pin, bool flag = false);
  virtual ~Coordinates();

  // data
  MeshBlock *pmy_block;  // ptr to MeshBlock containing this Coordinates
  AthenaArray<Real> dx1f, dx2f, dx3f, x1f, x2f, x3f; // face   spacing and positions
  AthenaArray<Real> dx1v, dx2v, dx3v, x1v, x2v, x3v; // volume spacing and positions
  AthenaArray<Real> x1s2, x1s3, x2s1, x2s3, x3s1, x3s2; // area averaged posn for AMR

  // functions...
  // ...to compute length of edges
  virtual void Edge1Length(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &len);
  virtual void Edge2Length(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &len);
  virtual void Edge3Length(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &len);
  virtual Real GetEdge1Length(const int k, const int j, const int i);
  virtual Real GetEdge2Length(const int k, const int j, const int i);
  virtual Real GetEdge3Length(const int k, const int j, const int i);

  // ...to compute physical width at cell center
  virtual Real CenterWidth1(const int k, const int j, const int i);
  virtual Real CenterWidth2(const int k, const int j, const int i);
  virtual Real CenterWidth3(const int k, const int j, const int i);

  // ...to compute area of faces
  virtual void Face1Area(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &area);
  virtual void Face2Area(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &area);
  virtual void Face3Area(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &area);
  virtual Real GetFace1Area(const int k, const int j, const int i);
  virtual Real GetFace2Area(const int k, const int j, const int i);
  virtual Real GetFace3Area(const int k, const int j, const int i);

  // ...to compute volume of cells
  virtual void CellVolume(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &vol);
  virtual Real GetCellVolume(const int k, const int j, const int i);

  // ...to compute geometrical source terms
  virtual void CoordSrcTerms(const Real dt, const AthenaArray<Real> *flux,
    const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc, AthenaArray<Real> &u);

  // ...to determine if index is a pole
  bool IsPole(int j);

#if GENERAL_RELATIVITY==1
  // In GR, functions...
  // ...to return private variables
  Real GetMass() const {return bh_mass_;}
  Real GetSpin() const {return bh_spin_;}

  // ...to compute metric
  virtual void CellMetric(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &g, AthenaArray<Real> &gi);
  virtual void Face1Metric(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &g, AthenaArray<Real> &g_inv);
  virtual void Face2Metric(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &g, AthenaArray<Real> &g_inv);
  virtual void Face3Metric(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &g, AthenaArray<Real> &g_inv);

  // ...to transform primitives to locally flat space
  virtual void PrimToLocal1(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &b1_vals, AthenaArray<Real> &prim_left,
    AthenaArray<Real> &prim_right, AthenaArray<Real> &bx);
  virtual void PrimToLocal2(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &b2_vals, AthenaArray<Real> &prim_left,
    AthenaArray<Real> &prim_right, AthenaArray<Real> &bx);
  virtual void PrimToLocal3(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &b3_vals, AthenaArray<Real> &prim_left,
    AthenaArray<Real> &prim_right, AthenaArray<Real> &bx);

  // ...to transform fluxes in locally flat space to global frame
  virtual void FluxToGlobal1(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &cons, const AthenaArray<Real> &bx,
    AthenaArray<Real> &flux);
  virtual void FluxToGlobal2(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &cons, const AthenaArray<Real> &bx,
    AthenaArray<Real> &flux);
  virtual void FluxToGlobal3(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &cons, const AthenaArray<Real> &bx,
    AthenaArray<Real> &flux);

  // ...to raise (lower) covariant (contravariant) components of a vector
  virtual void RaiseVectorCell(Real a_0, Real a_1, Real a_2, Real a_3,int k,int j,int i,
    Real *pa0, Real *pa1, Real *pa2, Real *pa3);
  virtual void LowerVectorCell(Real a0, Real a1, Real a2, Real a3, int k, int j, int i,
    Real *pa_0, Real *pa_1, Real *pa_2, Real *pa_3);

  // ...to convert from global to Boyer-Lindquist coordinates (useful for pgens)
  void GetBoyerLindquistCoordinates(Real x1, Real x2, Real x3,
    Real *pr, Real *ptheta, Real *pphi);

  // ...to transform 4-vector from Boyer-Lindquist to global coordinates (for pgens)
  virtual void TransformVectorCell(Real at, Real ax, Real ay, Real az,int k,int j,int i,
    Real *a0, Real *a1, Real *a2, Real *a3);
  virtual void TransformVectorFace1(Real at, Real ax, Real ay, Real az,int k,int j,int i,
    Real *a0, Real *a1, Real *a2, Real *a3);
  virtual void TransformVectorFace2(Real at, Real ax, Real ay, Real az,int k,int j,int i,
    Real *a0, Real *a1, Real *a2, Real *a3);
  virtual void TransformVectorFace3(Real at, Real ax, Real ay, Real az,int k,int j,int i,
    Real *a0, Real *a1, Real *a2, Real *a3);
#endif // GENERAL_RELATIVITY

protected:
  bool coarse_flag;  // true if this coordinate object is parent (coarse) mesh in AMR

  // Scratch arrays for coordinate factors
  // Format: coord_<type>[<direction>]_<index>[<count>]_
  //   type: vol[ume], area, etc.
  //   direction: 1/2/3 depending on which face, edge, etc. is in play
  //   index: i/j/k indicating which coordinates index array
  //   count: 1/2/... in case multiple arrays are needed for different terms
  AthenaArray<Real> coord_vol_i_, coord_vol_i1_, coord_vol_i2_;
  AthenaArray<Real> coord_vol_j_, coord_vol_j1_, coord_vol_j2_;
  AthenaArray<Real> coord_vol_k1_;
  AthenaArray<Real> coord_area1_i_, coord_area1_i1_;
  AthenaArray<Real> coord_area1_j_, coord_area1_j1_, coord_area1_j2_;
  AthenaArray<Real> coord_area1_k1_;
  AthenaArray<Real> coord_area2_i_, coord_area2_i1_, coord_area2_i2_;
  AthenaArray<Real> coord_area2_j_, coord_area2_j1_, coord_area2_j2_;
  AthenaArray<Real> coord_area2_k1_;
  AthenaArray<Real> coord_area3_i_, coord_area3_i1_, coord_area3_i2_;
  AthenaArray<Real> coord_area3_j1_, coord_area3_j2_;
  AthenaArray<Real> coord_len1_i1_, coord_len1_i2_;
  AthenaArray<Real> coord_len1_j1_, coord_len1_j2_;
  AthenaArray<Real> coord_len2_i1_;
  AthenaArray<Real> coord_len2_j1_, coord_len2_j2_;
  AthenaArray<Real> coord_len3_i1_;
  AthenaArray<Real> coord_len3_j1_, coord_len3_j2_;
  AthenaArray<Real> coord_len3_k1_;
  AthenaArray<Real> coord_width1_i1_;
  AthenaArray<Real> coord_width2_i1_;
  AthenaArray<Real> coord_width2_j1_;
  AthenaArray<Real> coord_width3_j1_, coord_width3_j2_, coord_width3_j3_;
  AthenaArray<Real> coord_width3_k1_;
  AthenaArray<Real> coord_width3_ji1_;
  AthenaArray<Real> coord_src_i1_, coord_src_i2_, coord_src_i3_, coord_src_i4_;
  AthenaArray<Real> coord_src_j1_, coord_src_j2_, coord_src_j3_;
  AthenaArray<Real> coord_src1_i_;
  AthenaArray<Real> coord_src1_j_;
  AthenaArray<Real> coord_src2_i_;
  AthenaArray<Real> coord_src2_j_;
  AthenaArray<Real> coord_src3_j_;

  // Scratch arrays for physical source terms
  AthenaArray<Real> phy_src1_i_, phy_src2_i_;

#if GENERAL_RELATIVITY==1
  // GR-specific scratch arrays
  AthenaArray<Real> metric_cell_i1_, metric_cell_i2_;
  AthenaArray<Real> metric_cell_j1_, metric_cell_j2_;
  AthenaArray<Real> metric_face1_i1_, metric_face1_i2_;
  AthenaArray<Real> metric_face1_j1_, metric_face1_j2_;
  AthenaArray<Real> metric_face2_i1_, metric_face2_i2_;
  AthenaArray<Real> metric_face2_j1_, metric_face2_j2_;
  AthenaArray<Real> metric_face3_i1_, metric_face3_i2_;
  AthenaArray<Real> metric_face3_j1_, metric_face3_j2_;
  AthenaArray<Real> trans_face1_i1_, trans_face1_i2_;
  AthenaArray<Real> trans_face1_j1_;
  AthenaArray<Real> trans_face1_ji1_, trans_face1_ji2_, trans_face1_ji3_,
      trans_face1_ji4_, trans_face1_ji5_, trans_face1_ji6_, trans_face1_ji7_;
  AthenaArray<Real> trans_face2_i1_, trans_face2_i2_;
  AthenaArray<Real> trans_face2_j1_;
  AthenaArray<Real> trans_face2_ji1_, trans_face2_ji2_, trans_face2_ji3_,
      trans_face2_ji4_, trans_face2_ji5_, trans_face2_ji6_;
  AthenaArray<Real> trans_face3_i1_, trans_face3_i2_;
  AthenaArray<Real> trans_face3_j1_;
  AthenaArray<Real> trans_face3_ji1_, trans_face3_ji2_, trans_face3_ji3_,
      trans_face3_ji4_, trans_face3_ji5_, trans_face3_ji6_;
  AthenaArray<Real> g_, gi_;

  Real bh_mass_; 
  Real bh_spin_; 

#endif // GENERAL_RELATIVITY
};

//----------------------------------------------------------------------------------------
//! \class Cartesian
//  \brief derived class for Cartesian coordinates.  None of the virtual funcs
//  in the Coordinates abstract base class need to be over-written.

class Cartesian : public Coordinates {
public:
  Cartesian(MeshBlock *pmb, ParameterInput *pin, bool flag);
  ~Cartesian();
};

//----------------------------------------------------------------------------------------
//! \class Cylindrical
//  \brief derived class for Cylindrical coordinates.  Some of the length, area,
//  and volume functions in the Coordinates abstract base class are over-written.

class Cylindrical : public Coordinates {
public:
  Cylindrical(MeshBlock *pmb, ParameterInput *pin, bool flag);
  ~Cylindrical();

  // functions...
  // ...to compute length of edges
  void Edge2Length(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &len);
  Real GetEdge2Length(const int k, const int j, const int i);

  // ...to compute physical width at cell center
  Real CenterWidth2(const int k, const int j, const int i);

  // ...to compute area of faces
  void Face1Area(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &area);
  void Face3Area(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &area);
  Real GetFace1Area(const int k, const int j, const int i);
  Real GetFace3Area(const int k, const int j, const int i);

  // ...to compute volumes of cells
  void CellVolume(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &vol);
  Real GetCellVolume(const int k, const int j, const int i);

  // ...to compute geometrical source terms
  void CoordSrcTerms(const Real dt, const AthenaArray<Real> *flux,
    const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc, AthenaArray<Real> &u);
};

//----------------------------------------------------------------------------------------
//! \class SphericalPolar
//  \brief derived class for spherical polar coordinates.  Many of the length, area,
//  and volume functions in the Coordinates abstract base class are over-written.

class SphericalPolar : public Coordinates {
public:
  SphericalPolar(MeshBlock *pmb, ParameterInput *pin, bool flag);
  ~SphericalPolar();

  // functions...
  // ...to compute length of edges
  void Edge2Length(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &len);
  void Edge3Length(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &len);
  Real GetEdge2Length(const int k, const int j, const int i);
  Real GetEdge3Length(const int k, const int j, const int i);

  // ...to compute physical width at cell center
  Real CenterWidth2(const int k, const int j, const int i);
  Real CenterWidth3(const int k, const int j, const int i);

  // ...to compute area of faces
  void Face1Area(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &area);
  void Face2Area(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &area);
  void Face3Area(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &area);
  Real GetFace1Area(const int k, const int j, const int i);
  Real GetFace2Area(const int k, const int j, const int i);
  Real GetFace3Area(const int k, const int j, const int i);

  // ...to compute volumes of cells
  void CellVolume(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &vol);
  Real GetCellVolume(const int k, const int j, const int i);

  // ...to compute geometrical source terms
  void CoordSrcTerms(const Real dt, const AthenaArray<Real> *flux,
    const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc, AthenaArray<Real> &u);
};

//----------------------------------------------------------------------------------------
//! \class Minkowski
//  \brief derived class for Minkowski (flat) spacetime and Cartesian coordinates in GR.
//  None of the length, area, and volume functions in the abstract base class need to be
//  overwritten, but all the metric and transforms functions are.

class Minkowski : public Coordinates {
public:
  Minkowski(MeshBlock *pmb, ParameterInput *pin, bool flag);
  ~Minkowski();

#if GENERAL_RELATIVITY==1
  // In GR, functions...
  // ...to compute metric
  void CellMetric(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &g, AthenaArray<Real> &gi);
  void Face1Metric(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &g, AthenaArray<Real> &g_inv);
  void Face2Metric(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &g, AthenaArray<Real> &g_inv);
  void Face3Metric(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &g, AthenaArray<Real> &g_inv);

  // ...to transform primitives to locally flat space
  void PrimToLocal1(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &b1_vals, AthenaArray<Real> &prim_left,
    AthenaArray<Real> &prim_right, AthenaArray<Real> &bx);
  void PrimToLocal2(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &b2_vals, AthenaArray<Real> &prim_left,
    AthenaArray<Real> &prim_right, AthenaArray<Real> &bx);
  void PrimToLocal3(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &b3_vals, AthenaArray<Real> &prim_left,
    AthenaArray<Real> &prim_right, AthenaArray<Real> &bx);

  // ...to transform fluxes in locally flat space to global frame
  void FluxToGlobal1(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &cons, const AthenaArray<Real> &bx,
    AthenaArray<Real> &flux);
  void FluxToGlobal2(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &cons, const AthenaArray<Real> &bx,
    AthenaArray<Real> &flux);
  void FluxToGlobal3(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &cons, const AthenaArray<Real> &bx,
    AthenaArray<Real> &flux);

  // .. to transform 4-vector from Boyer-Lindquist to global coordinates
  void TransformVectorCell(Real at, Real ax, Real ay, Real az,int k,int j,int i,
    Real *a0, Real *a1, Real *a2, Real *a3);
  void TransformVectorFace1(Real at, Real ax, Real ay, Real az,int k,int j,int i,
    Real *a0, Real *a1, Real *a2, Real *a3);
  void TransformVectorFace2(Real at, Real ax, Real ay, Real az,int k,int j,int i,
    Real *a0, Real *a1, Real *a2, Real *a3);
  void TransformVectorFace3(Real at, Real ax, Real ay, Real az,int k,int j,int i,
    Real *a0, Real *a1, Real *a2, Real *a3);

  // for raising (lowering) covariant (contravariant) components of a vector
  void RaiseVectorCell(Real a_0, Real a_1, Real a_2, Real a_3,int k,int j,int i,
    Real *pa0, Real *pa1, Real *pa2, Real *pa3);
  void LowerVectorCell(Real a0, Real a1, Real a2, Real a3, int k, int j, int i,
    Real *pa_0, Real *pa_1, Real *pa_2, Real *pa_3);
#endif // GENERAL_RELATIVITY
};

//----------------------------------------------------------------------------------------
//! \class Schwarzschild
//  \brief derived class for Schwarzschild spacetime and spherical polar coordinates in GR
//  Nearly every function in the abstract base class need to be overwritten.

class Schwarzschild : public Coordinates {
public:
  Schwarzschild(MeshBlock *pmb, ParameterInput *pin, bool flag);
  ~Schwarzschild();

  // functions...
  // ...to compute length of edges
  void Edge1Length(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &len);
  void Edge2Length(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &len);
  void Edge3Length(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &len);
  Real GetEdge1Length(const int k, const int j, const int i);
  Real GetEdge2Length(const int k, const int j, const int i);
  Real GetEdge3Length(const int k, const int j, const int i);

  // ...to compute physical width at cell center
  Real CenterWidth1(const int k, const int j, const int i);
  Real CenterWidth2(const int k, const int j, const int i);
  Real CenterWidth3(const int k, const int j, const int i);

  // ...to compute area of faces
  void Face1Area(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &area);
  void Face2Area(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &area);
  void Face3Area(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &area);
  Real GetFace1Area(const int k, const int j, const int i);
  Real GetFace2Area(const int k, const int j, const int i);
  Real GetFace3Area(const int k, const int j, const int i);

  // ...to compute volumes of cells
  void CellVolume(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &vol);
  Real GetCellVolume(const int k, const int j, const int i);

  // ...to compute geometrical source terms
  void CoordSrcTerms(const Real dt, const AthenaArray<Real> *flux,
    const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc, AthenaArray<Real> &u);

#if GENERAL_RELATIVITY==1
  // In GR, functions...
  // ...to compute metric
  void CellMetric(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &g, AthenaArray<Real> &gi);
  void Face1Metric(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &g, AthenaArray<Real> &g_inv);
  void Face2Metric(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &g, AthenaArray<Real> &g_inv);
  void Face3Metric(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &g, AthenaArray<Real> &g_inv);

  // ...to transform primitives to locally flat space
  void PrimToLocal1(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &b1_vals, AthenaArray<Real> &prim_left,
    AthenaArray<Real> &prim_right, AthenaArray<Real> &bx);
  void PrimToLocal2(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &b2_vals, AthenaArray<Real> &prim_left,
    AthenaArray<Real> &prim_right, AthenaArray<Real> &bx);
  void PrimToLocal3(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &b3_vals, AthenaArray<Real> &prim_left,
    AthenaArray<Real> &prim_right, AthenaArray<Real> &bx);

  // ...to transform fluxes in locally flat space to global frame
  void FluxToGlobal1(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &cons, const AthenaArray<Real> &bx,
    AthenaArray<Real> &flux);
  void FluxToGlobal2(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &cons, const AthenaArray<Real> &bx,
    AthenaArray<Real> &flux);
  void FluxToGlobal3(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &cons, const AthenaArray<Real> &bx,
    AthenaArray<Real> &flux);

  // .. to transform 4-vector from Boyer-Lindquist to global coordinates
  void TransformVectorCell(Real at, Real ax, Real ay, Real az,int k,int j,int i,
    Real *a0, Real *a1, Real *a2, Real *a3);
  void TransformVectorFace1(Real at, Real ax, Real ay, Real az,int k,int j,int i,
    Real *a0, Real *a1, Real *a2, Real *a3);
  void TransformVectorFace2(Real at, Real ax, Real ay, Real az,int k,int j,int i,
    Real *a0, Real *a1, Real *a2, Real *a3);
  void TransformVectorFace3(Real at, Real ax, Real ay, Real az,int k,int j,int i,
    Real *a0, Real *a1, Real *a2, Real *a3);

  // for raising (lowering) covariant (contravariant) components of a vector
  void RaiseVectorCell(Real a_0, Real a_1, Real a_2, Real a_3,int k,int j,int i,
    Real *pa0, Real *pa1, Real *pa2, Real *pa3);
  void LowerVectorCell(Real a0, Real a1, Real a2, Real a3, int k, int j, int i,
    Real *pa_0, Real *pa_1, Real *pa_2, Real *pa_3);
#endif // GENERAL_RELATIVITY
};

//----------------------------------------------------------------------------------------
//! \class KerrSchild
//  \brief derived class for Kerr spacetime and Kerr-Schild coordinates in GR.
//  Nearly every function in the abstract base class need to be overwritten.

class KerrSchild : public Coordinates {
public:
  KerrSchild(MeshBlock *pmb, ParameterInput *pin, bool flag);
  ~KerrSchild();

  // functions...
  // ...to compute length of edges
  void Edge1Length(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &len);
  void Edge2Length(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &len);
  void Edge3Length(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &len);
  Real GetEdge1Length(const int k, const int j, const int i);
  Real GetEdge2Length(const int k, const int j, const int i);
  Real GetEdge3Length(const int k, const int j, const int i);

  // ...to compute physical width at cell center
  Real CenterWidth1(const int k, const int j, const int i);
  Real CenterWidth2(const int k, const int j, const int i);
  Real CenterWidth3(const int k, const int j, const int i);

  // ...to compute area of faces
  void Face1Area(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &area);
  void Face2Area(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &area);
  void Face3Area(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &area);
  Real GetFace1Area(const int k, const int j, const int i);
  Real GetFace2Area(const int k, const int j, const int i);
  Real GetFace3Area(const int k, const int j, const int i);

  // ...to compute volumes of cells
  void CellVolume(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &vol);
  Real GetCellVolume(const int k, const int j, const int i);

  // ...to compute geometrical source terms
  void CoordSrcTerms(const Real dt, const AthenaArray<Real> *flux,
    const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc, AthenaArray<Real> &u);

#if GENERAL_RELATIVITY==1
  // In GR, functions...
  // ...to compute metric
  void CellMetric(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &g, AthenaArray<Real> &gi);
  void Face1Metric(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &g, AthenaArray<Real> &g_inv);
  void Face2Metric(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &g, AthenaArray<Real> &g_inv);
  void Face3Metric(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &g, AthenaArray<Real> &g_inv);

  // ...to transform primitives to locally flat space
  void PrimToLocal1(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &b1_vals, AthenaArray<Real> &prim_left,
    AthenaArray<Real> &prim_right, AthenaArray<Real> &bx);
  void PrimToLocal2(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &b2_vals, AthenaArray<Real> &prim_left,
    AthenaArray<Real> &prim_right, AthenaArray<Real> &bx);
  void PrimToLocal3(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &b3_vals, AthenaArray<Real> &prim_left,
    AthenaArray<Real> &prim_right, AthenaArray<Real> &bx);

  // ...to transform fluxes in locally flat space to global frame
  void FluxToGlobal1(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &cons, const AthenaArray<Real> &bx,
    AthenaArray<Real> &flux);
  void FluxToGlobal2(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &cons, const AthenaArray<Real> &bx,
    AthenaArray<Real> &flux);
  void FluxToGlobal3(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &cons, const AthenaArray<Real> &bx,
    AthenaArray<Real> &flux);

  // .. to transform 4-vector from Boyer-Lindquist to global coordinates
  void TransformVectorCell(Real at, Real ax, Real ay, Real az,int k,int j,int i,
    Real *a0, Real *a1, Real *a2, Real *a3);
  void TransformVectorFace1(Real at, Real ax, Real ay, Real az,int k,int j,int i,
    Real *a0, Real *a1, Real *a2, Real *a3);
  void TransformVectorFace2(Real at, Real ax, Real ay, Real az,int k,int j,int i,
    Real *a0, Real *a1, Real *a2, Real *a3);
  void TransformVectorFace3(Real at, Real ax, Real ay, Real az,int k,int j,int i,
    Real *a0, Real *a1, Real *a2, Real *a3);

  // for raising (lowering) covariant (contravariant) components of a vector
  void RaiseVectorCell(Real a_0, Real a_1, Real a_2, Real a_3,int k,int j,int i,
    Real *pa0, Real *pa1, Real *pa2, Real *pa3);
  void LowerVectorCell(Real a0, Real a1, Real a2, Real a3, int k, int j, int i,
    Real *pa_0, Real *pa_1, Real *pa_2, Real *pa_3);
#endif // GENERAL_RELATIVITY
};

#endif // COORDINATES_HPP
