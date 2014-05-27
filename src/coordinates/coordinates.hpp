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

//! \class Coordinates
//  \brief coordinate data and functions

namespace COORDINATE_SYSTEM {

class Coordinates {
public:
  Coordinates(Block *pb);
  ~Coordinates();

  Block *pparent_block;

  void Area1Face(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &area);
  void Area2Face(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &area);
  void Area3Face(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &area);
  void CellVolume(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &area);
  void SourceTerms(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &src);

  AthenaArray<Real> face_area, cell_volume;
};
} // end namespace COORDINATE_SYSTEM
#endif
