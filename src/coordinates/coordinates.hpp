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
class Block;

//! \class Coordinates
//  \brief coordinate data and functions

class Coordinates {
public:
  Coordinates(Block *pb);
  ~Coordinates();

  Block *pparent_block;

  void Area1Face(
    const int k, const int j, const int il, const int iu, AthenaArray<Real> &area);
  void Area2Face(
    const int k, const int j, const int il, const int iu, AthenaArray<Real> &area);
  void Area3Face(
    const int k, const int j, const int il, const int iu, AthenaArray<Real> &area);
  void CellVolume(
    const int k, const int j, const int il, const int iu, AthenaArray<Real> &area);
  void CoordinateSourceTerms(
    const int k, const int j, AthenaArray<Real> &prim, AthenaArray<Real> &src);

  void CellMetric(const int k, const int j, AthenaArray<Real> &g,
      AthenaArray<Real> &g_inv);
  void PrimToLocal(const int k, const int j, const int il, const int iu, const int ivx,
      AthenaArray<Real> &prim);
  void FluxToGlobal(const int k, const int j, const int il, const int iu, const int ivx,
      AthenaArray<Real> &flux);

  AthenaArray<Real> face_area, cell_volume;

private:
  AthenaArray<Real> face1_area_i_, face1_area_j_;
  AthenaArray<Real> face2_area_i_, face2_area_j_;
  AthenaArray<Real> face3_area_i_, face3_area_j_;
  AthenaArray<Real> src_terms_i_,  src_terms_j_;
  AthenaArray<Real> volume_i_,     volume_j_;
};
#endif
