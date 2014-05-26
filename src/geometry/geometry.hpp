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

//namespace COORDINATE_SYSTEM {

class Geometry {
public:
  Geometry(Block *pb);
  ~Geometry();

  Block *pparent_block;

  AthenaArray<Real> dx1v, dx2v, dx3v, x1v, x2v, x3v;

  void InitGeometryFactors(ParameterInput *pin);
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
//} // end namespace COORDINATE_SYSTEM
#endif
