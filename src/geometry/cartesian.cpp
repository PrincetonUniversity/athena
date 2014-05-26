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
#include "geometry.hpp"

//======================================================================================
//! \file cartesian.cpp
//  \brief implements geometry functions for Cartesian coordinates (the trivial case)
//======================================================================================

//namespace cartesian_coordinates {

//--------------------------------------------------------------------------------------
// \!fn 
// \brief

void Geometry::InitGeometryFactors(ParameterInput *pin)
{
  int is = pparent_block->is; int js = pparent_block->js; int ks = pparent_block->ks;
  int ie = pparent_block->ie; int je = pparent_block->je; int ke = pparent_block->ke;

// initialize volume-averaged positions and spacing
// x1-direction

  for (int i=is-(NGHOST); i<=ie+(NGHOST); ++i) {
    x1v(i) = (pparent_block->x1f(i+1) - pparent_block->x1f(i))/2.0;
  }
  for (int i=is-(NGHOST)+1; i<=ie+(NGHOST); ++i) {
    dx1v(i) = x1v(i) - x1v(i-1);
  }
}

//--------------------------------------------------------------------------------------
// \!fn 
// \brief

void Geometry::Area1Face(const int k, const int j, const int il, const int iu, 
  AthenaArray<Real> &area)
{
#pragma simd
  for (int i=il; i<=iu; ++i){
    Real& area_i = area(i);
    area_i = (pparent_block->dx2f(j))*(pparent_block->dx3f(k));
  }
  return;
}

void Geometry::Area2Face(const int k, const int j, const int il, const int iu, 
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

void Geometry::Area3Face(const int k, const int j, const int il, const int iu, 
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

void Geometry::CellVolume(const int k, const int j, const int il, const int iu, 
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

void Geometry::SourceTerms(const int k, const int j, const int il, const int iu, 
  AthenaArray<Real> &src)
{
  return;
}

//} // end cartesian_ccordinates namespace
