#ifndef COORDINATES_HPP
#define COORDINATES_HPP
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file coordinates.hpp
//  \brief defines abstract base class Coordinates containing data and functions used by
//  all coordinate derived classes.  The Coordinates class is used to compute/store
//  geometrical factors (areas, volumes, coordinate source terms) related to a Mesh. The
//  GR coordinates class is derived from this ABC, and contains additional data/fns

// Athena++ classes headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../mesh/mesh.hpp"
#include <iostream>

// forward declarations
class MeshBlock;
class ParameterInput;

//! \class Coordinates
//  \brief abstract base class for coordinate data and functions

class Coordinates {
public:
//  friend class HydroSourceTerms;
  Coordinates(MeshBlock *pmb, ParameterInput *pin, int flag = 0);
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

  // ...to compute metric (no-op except in GR)
  virtual void CellMetric(const int, const int, const int, const int, AthenaArray<Real> &,
        AthenaArray<Real> &) {return;}

  // ...to compute volume of cells (pure virtual)
  virtual void CellVolume(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &vol)=0;
  virtual Real GetCellVolume(const int k, const int j, const int i)=0;

  // ...to compute geometrical source terms (pure virtual)
  virtual void CoordSrcTerms(const Real dt, const AthenaArray<Real> *flux,
    const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc, AthenaArray<Real> &u)=0;

protected:
  int cflag;
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
};

//----------------------------------------------------------------------------------------
//! \class Cartesian
//  \brief derived Coordinates class for Cartesian coordinates.  None of the virtual funcs
//  in the abstract base class are over-written.

class Cartesian : public Coordinates {
public:
  Cartesian(MeshBlock *pmb, ParameterInput *pin, int flag);
  ~Cartesian();

  // functions...
  // ...to compute volumes of cells
  void CellVolume(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &vol);
  Real GetCellVolume(const int k, const int j, const int i);

  // ...to compute geometrical source terms
  void CoordSrcTerms(const Real dt, const AthenaArray<Real> *flux,
    const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc, AthenaArray<Real> &u);
};

//----------------------------------------------------------------------------------------
//! \class Cylindrical
//  \brief derived Coordinates class for Cylindrical coordinates.  Some of the length
//  and area functions in the abstract base class are over-written.

class Cylindrical : public Coordinates {
public:
  Cylindrical(MeshBlock *pmb, ParameterInput *pin, int flag);
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
//  \brief derived Coordinates class for spherical polar coordinates.  Some of the length
//  and area functions in the abstract base class are over-written.

class SphericalPolar : public Coordinates {
public:
  SphericalPolar(MeshBlock *pmb, ParameterInput *pin, int flag);
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

  // ...to compute geometrical source terms (pure virtual)
  void CoordSrcTerms(const Real dt, const AthenaArray<Real> *flux,
    const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc, AthenaArray<Real> &u);
};

#endif // COORDINATES_HPP
