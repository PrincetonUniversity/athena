#ifndef GEOMETRY_HPP
#define GEOMETRY_HPP
//======================================================================================
/* Athena++ astrophysical MHD code
 * Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
 * See LICENSE file for full public license information.
 *====================================================================================*/
/*! \file geometry.hpp
 *  \brief defines Geometry class used to compute/store geometrical factors (areas,
 *  volumes, source terms) related to a Mesh
 *====================================================================================*/

//! \class Geometry
//  \brief geometry data and functions

class Geometry {
public:
  Geometry(Block *pb);
  ~Geometry();

  Block *pmy_block;

  void Area1Face(const int k, const int j, const int il, const int iu,
       AthenaArray<Real> &area);
  void Area2Face(const int k, const int j, const int il, const int iu,
       AthenaArray<Real> &area);
  void Area3Face(const int k, const int j, const int il, const int iu,
       AthenaArray<Real> &area);
  void VolumeOfCell(const int k, const int j, const int il, const int iu,
       AthenaArray<Real> &area);
  void SourceTerms(const int k, const int j, const int il, const int iu,
       AthenaArray<Real> &src);

  AthenaArray<Real> face_area, cell_volume;
};
#endif
