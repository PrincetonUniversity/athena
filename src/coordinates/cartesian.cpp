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
 * You should have received a copy of GNU GPL in the file LICENSE included in the code
 * distribution.  If not see <http://www.gnu.org/licenses/>.
 *====================================================================================*/

// Primary header
#include "coordinates.hpp"

// Athena++ headers
#include "../athena.hpp"          // macros, Real
#include "../athena_arrays.hpp"   // AthenaArray
#include "../parameter_input.hpp" // ParameterInput
#include "../mesh.hpp"            // MeshBlock

//======================================================================================
//! \file cartesian.cpp
//  \brief implements functions in class Coordinates for Cartesian coordinates
//======================================================================================

//--------------------------------------------------------------------------------------
// Coordinates constructor

Coordinates::Coordinates(MeshBlock *pmb, ParameterInput *pin)
{
  pmy_block = pmb;
  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;

// initialize volume-averaged positions and spacing
// x1-direction: x1v = dx/2

  for (int i=is-(NGHOST); i<=ie+(NGHOST); ++i) {
    pmb->x1v(i) = 0.5*(pmb->x1f(i+1) + pmb->x1f(i));
  }
  for (int i=is-(NGHOST); i<=ie+(NGHOST)-1; ++i) {
    pmb->dx1v(i) = pmb->x1v(i+1) - pmb->x1v(i);
  }

// x2-direction: x2v = dy/2

  if (pmb->block_size.nx2 == 1) {
    pmb->x2v(js) = 0.5*(pmb->x2f(js+1) + pmb->x2f(js));
    pmb->dx2v(js) = pmb->dx2f(js);
  } else {
    for (int j=js-(NGHOST); j<=je+(NGHOST); ++j) {
      pmb->x2v(j) = 0.5*(pmb->x2f(j+1) + pmb->x2f(j));
    }
    for (int j=js-(NGHOST); j<=je+(NGHOST)-1; ++j) {
      pmb->dx2v(j) = pmb->x2v(j+1) - pmb->x2v(j);
    }
  }

// x3-direction: x3v = dz/2

  if (pmb->block_size.nx3 == 1) {
    pmb->x3v(ks) = 0.5*(pmb->x3f(ks+1) + pmb->x3f(ks));
    pmb->dx3v(ks) = pmb->dx3f(ks);
  } else {
    for (int k=ks-(NGHOST); k<=ke+(NGHOST); ++k) {
      pmb->x3v(k) = 0.5*(pmb->x3f(k+1) + pmb->x3f(k));
    }
    for (int k=ks-(NGHOST); k<=ke+(NGHOST)-1; ++k) {
      pmb->dx3v(k) = pmb->x3v(k+1) - pmb->x3v(k);
    }
  }
}

// destructor

Coordinates::~Coordinates()
{
}

//--------------------------------------------------------------------------------------
// Edge Length functions


//--------------------------------------------------------------------------------------
// Face Area functions

/* \!fn void Coordinates::Face1Area(const int k,const int j, const int il, const int iu,
      AthenaArray<Real> &area)
 * \brief  functions to compute area of cell faces in each direction    */

void Coordinates::Face1Area(const int k, const int j, const int il, const int iu,
  AthenaArray<Real> *parea)
{
#pragma simd
  for (int i=il; i<=iu; ++i){
    Real& area_i = (*parea)(i);
    area_i = (pmy_block->dx2f(j))*(pmy_block->dx3f(k));  // dy*dz
  }
  return;
}

void Coordinates::Face2Area(const int k, const int j, const int il, const int iu,
  AthenaArray<Real> *parea)
{
  AthenaArray<Real> dx1f = pmy_block->dx1f.ShallowCopy();
#pragma simd
  for (int i=il; i<=iu; ++i){
    Real& area_i = (*parea)(i);
    Real& dx1_i  = dx1f(i);
    area_i = dx1_i*(pmy_block->dx3f(k));  // dx*dz
  }
  return;
}

void Coordinates::Face3Area(const int k, const int j, const int il, const int iu,
  AthenaArray<Real> *parea)
{
  AthenaArray<Real> dx1f = pmy_block->dx1f.ShallowCopy();
#pragma simd
  for (int i=il; i<=iu; ++i){
    Real& area_i = (*parea)(i);
    Real& dx1_i  = dx1f(i);
    area_i = dx1_i*(pmy_block->dx2f(j));  // dx*dy
  }
  return;
}

//--------------------------------------------------------------------------------------
// Cell Volume function

/* \!fn void Coordinates::CellVolume(const int k,const int j,const int il, const int iu,
 *   AthenaArray<Real> &vol)
 * \brief function to compute cell volume    */

void Coordinates::CellVolume(const int k, const int j, const int il, const int iu,
  AthenaArray<Real> *pvol)
{
  AthenaArray<Real> dx1f = pmy_block->dx1f.ShallowCopy();
#pragma simd
  for (int i=il; i<=iu; ++i){
    Real& vol_i = (*pvol)(i);
    Real& dx1_i = dx1f(i);
    vol_i = dx1_i*(pmy_block->dx2f(j))*(pmy_block->dx3f(k));  // dx*dy*dz
  }
  return;
}

//--------------------------------------------------------------------------------------
// Width and Distance functions

Real Coordinates::CellPhysicalWidth1(const int k, const int j, const int i)
{
  return (pmy_block->dx1f(i));
}

Real Coordinates::CellPhysicalWidth2(const int k, const int j, const int i)
{
  return (pmy_block->dx2f(j));
}

Real Coordinates::CellPhysicalWidth3(const int k, const int j, const int i)
{
  return (pmy_block->dx3f(k));
}

ThreeVector Coordinates::VectorBetweenPoints(const ThreeVector p1, const ThreeVector p2)
{
  ThreeVector r;
  r.x1 = p1.x1 - p2.x1;
  r.x2 = p1.x2 - p2.x2;
  r.x3 = p1.x3 - p2.x3;
  return r;
}

//--------------------------------------------------------------------------------------
/* \!fn void Coordinates::CoordinateSourceTerms(
 *   const int k, const int j, AthenaArray<Real> &prim, AthenaArray<Real> &src)
 * \brief function to compute coordinate source terms (no-op function for cartesian)  */

void Coordinates::CoordinateSourceTerms(const Real dt, const AthenaArray<Real> &prim,
  AthenaArray<Real> &src)
{
  return;
}
