//======================================================================================
/* Athena++ astrophysical MHD code
 * Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
 *
 * This program is free software: you can redistribute and/or modify it under the terms
 * of the GNU General Public License (GPL) as published by the Free Software Foundation,
 * either version 3 of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
 * PARTICULAR PURPOSE.  See the GNU General Public License for more details.
 *
 * You should have received a copy of GNU GPL in the file LICENSE included in
 * the code distribution.  If not see <http://www.gnu.org/licenses/>.
 *====================================================================================*/

#include <iostream>
#include <string>
#include <math.h>
#include <float.h>

#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../parameter_input.hpp"
#include "../mesh.hpp"
#include "coordinates.hpp"

//======================================================================================
//! \file cartesian.cpp
//  \brief implements functions in class Coordinates for Cartesian coordinates
//======================================================================================

namespace cartesian_coordinates {

// constructor

Coordinates::Coordinates(Block *pb)
{
  pparent_block = pb;

  int is = pb->is; int js = pb->js; int ks = pb->ks;
  int ie = pb->ie; int je = pb->je; int ke = pb->ke;

  int ncells1 = pb->block_size.nx1 + 2*(NGHOST);
  int ncells2 = 1, ncells3 = 1;
  if (pb->block_size.nx2 > 1) ncells2 = pb->block_size.nx2 + 2*(NGHOST);
  if (pb->block_size.nx3 > 1) ncells3 = pb->block_size.nx3 + 2*(NGHOST);
  face_area.NewAthenaArray(ncells1);
  cell_volume.NewAthenaArray(ncells1);

// initialize volume-averaged positions and spacing
// x1-direction

  for (int i=is-(NGHOST); i<=ie+(NGHOST); ++i) {
    pb->x1v(i) = 0.5*(pb->x1f(i+1) + pb->x1f(i));
  }
  for (int i=is-(NGHOST)+1; i<=ie+(NGHOST); ++i) {
    pb->dx1v(i) = pb->x1v(i) - pb->x1v(i-1);
  }

// x2-direction

  if (ncells2 == 1) {
    pb->x2v(js) = 0.5*(pb->x2f(js+1) + pb->x2f(js));
    pb->dx2v(js) = pb->dx2f(js);
  } else {
    for (int j=js-(NGHOST); j<=je+(NGHOST); ++j) {
      pb->x2v(j) = 0.5*(pb->x2f(j+1) + pb->x2f(j));
    }
    for (int j=js-(NGHOST)+1; j<=je+(NGHOST); ++j) {
      pb->dx2v(j) = pb->x2v(j) - pb->x2v(j-1);
    }
  }

// x3-direction

  if (ncells3 == 1) {
    pb->x3v(ks) = 0.5*(pb->x3f(ks+1) + pb->x3f(ks));
    pb->dx3v(ks) = pb->dx3f(ks);
  } else {
    for (int k=ks-(NGHOST); k<=ke+(NGHOST); ++k) {
      pb->x3v(k) = 0.5*(pb->x3f(k+1) + pb->x3f(k));
    }
    for (int k=ks-(NGHOST)+1; k<=ke+(NGHOST); ++k) {
      pb->dx3v(k) = pb->x3v(k) - pb->x3v(k-1);
    }
  }
}

// destructor

Coordinates::~Coordinates()
{
  face_area.DeleteAthenaArray();
  cell_volume.DeleteAthenaArray();
}

//--------------------------------------------------------------------------------------
// \!fn 
// \brief

void Coordinates::Area1Face(const int k, const int j, const int il, const int iu,
  AthenaArray<Real> &area)
{
#pragma simd
  for (int i=il; i<=iu; ++i){
    Real& area_i = area(i);
    area_i = (pparent_block->dx2f(j))*(pparent_block->dx3f(k));
  }
  return;
}

void Coordinates::Area2Face(const int k, const int j, const int il, const int iu,
  AthenaArray<Real> &area)
{
  AthenaArray<Real> dx1f = pparent_block->dx1f.ShallowCopy();
#pragma simd
  for (int i=il; i<=iu; ++i){
    Real& area_i = area(i);
    Real& dx1_i  = dx1f(i);
    area_i = dx1_i*(pparent_block->dx3f(k));
  }
  return;
}

void Coordinates::Area3Face(const int k, const int j, const int il, const int iu,
  AthenaArray<Real> &area)
{
  AthenaArray<Real> dx1f = pparent_block->dx1f.ShallowCopy();
#pragma simd
  for (int i=il; i<=iu; ++i){
    Real& area_i = area(i);
    Real& dx1_i  = dx1f(i);
    area_i = dx1_i*(pparent_block->dx2f(j));
  }
  return;
}

//--------------------------------------------------------------------------------------
// \!fn 
// \brief

void Coordinates::CellVolume(const int k, const int j, const int il, const int iu,
  AthenaArray<Real> &vol)
{
  AthenaArray<Real> dx1f = pparent_block->dx1f.ShallowCopy();
#pragma simd
  for (int i=il; i<=iu; ++i){
    Real& vol_i = vol(i);
    Real& dx1_i = dx1f(i);
    vol_i = dx1_i*(pparent_block->dx2f(j))*(pparent_block->dx3f(k));
  }
  return;
}

//--------------------------------------------------------------------------------------
// \!fn 
// \brief

void Coordinates::SourceTerms(const int k, const int j, const int il, const int iu, 
  AthenaArray<Real> &src)
{
  return;
}

} // end cartesian_ccordinates namespace
