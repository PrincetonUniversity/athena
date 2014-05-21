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
#include "../mesh.hpp"
#include "geometry.hpp"

//======================================================================================
//! \file cartesian.cpp
//  \brief implements geometry functions for Cartesian coordinates (the trivial case)
//======================================================================================


//--------------------------------------------------------------------------------------
// \!fn 
// \brief

void Geometry::Area1Face(const int k, const int j, const int il, const int iu, 
  AthenaArray<Real> &area)
{
#pragma simd
  for (int i=il; i<=iu; ++i){
    Real& area_i = area(i);

    area_i = (pmy_block->dx2v(j))*(pmy_block->dx3v(k));
  }
  return;
}

void Geometry::Area2Face(const int k, const int j, const int il, const int iu, 
  AthenaArray<Real> &area)
{
  AthenaArray<Real> dx1v = pmy_block->dx1v.ShallowCopy();
#pragma simd
  for (int i=il; i<=iu; ++i){
    Real& area_i = area(i);
    Real& dx1_i  = dx1v(i);

    area_i = dx1_i*(pmy_block->dx3v(k));
  }
  return;
}

void Geometry::Area3Face(const int k, const int j, const int il, const int iu, 
  AthenaArray<Real> &area)
{
  AthenaArray<Real> dx1v = pmy_block->dx1v.ShallowCopy();
#pragma simd
  for (int i=il; i<=iu; ++i){
    Real& area_i = area(i);
    Real& dx1_i  = dx1v(i);

    area_i = dx1_i*(pmy_block->dx2v(j));
  }
  return;
}

//--------------------------------------------------------------------------------------
// \!fn 
// \brief

void Geometry::VolumeOfCell(const int k, const int j, const int il, const int iu, 
  AthenaArray<Real> &vol)
{
  AthenaArray<Real> dx1v = pmy_block->dx1v.ShallowCopy();
#pragma simd
  for (int i=il; i<=iu; ++i){
    Real& vol_i = vol(i);
    Real& dx1_i = dx1v(i);

    vol_i = dx1_i*(pmy_block->dx2v(j))*(pmy_block->dx3v(k));
  }
  return;
}

//--------------------------------------------------------------------------------------
// \!fn 
// \brief

void Geometry::SourceTerms(const int k, const int j, const int il, const int iu, 
  AthenaArray<Real> &src)
{
  return;
}
