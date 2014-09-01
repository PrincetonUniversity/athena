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

// Athena headers
#include "../athena.hpp"         // macros, Real
#include "../athena_arrays.hpp"  // AthenaArray
#include "../mesh.hpp"           // MeshBlock

//======================================================================================
//! \file cartesian.cpp
//  \brief implements functions in class Coordinates for Cartesian coordinates
//======================================================================================

// Coordinates constructor

Coordinates::Coordinates(MeshBlock *pmb, ParameterInput *pin)
{
  pmy_block = pmb;
  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;

// initialize volume-averaged positions and spacing
// x1-direction

  for (int i=is-(NGHOST); i<=ie+(NGHOST); ++i) {
    pmb->x1v(i) = 0.5*(pmb->x1f(i+1) + pmb->x1f(i));
  }
  for (int i=is-(NGHOST); i<=ie+(NGHOST)-1; ++i) {
    pmb->dx1v(i) = pmb->x1v(i+1) - pmb->x1v(i);
  }

// x2-direction

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

// x3-direction

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

// Allocate memory for scratch arrays used in integrator, and internal scratch arrays
// For cartesian coordinates, no local scratch arrays are needed

  int ncells1 = pmb->block_size.nx1 + 2*(NGHOST);
  int ncells2 = 1;
  if (pmb->block_size.nx2 > 1) ncells2 = pmb->block_size.nx2 + 2*(NGHOST);

  face_area.NewAthenaArray(ncells1);
  cell_volume.NewAthenaArray(ncells1);
}

// destructor

Coordinates::~Coordinates()
{
  pmy_block = NULL; // MeshBlock destructor will free this memory
  face_area.DeleteAthenaArray();
  cell_volume.DeleteAthenaArray();
}

//--------------------------------------------------------------------------------------
/* \!fn void Coordinates::Area1Face(const int k, const int j, const int il,
 *  const int iu, AthenaArray<Real> &area)
 * \brief  functions to compute area at each face of a grid cell    */

void Coordinates::Area1Face(const int k, const int j, const int il, const int iu,
  AthenaArray<Real> &area)
{
#pragma simd
  for (int i=il; i<=iu; ++i){
    Real& area_i = area(i);
    area_i = (pmy_block->dx2f(j))*(pmy_block->dx3f(k));
  }
  return;
}

void Coordinates::Area2Face(const int k, const int j, const int il, const int iu,
  AthenaArray<Real> &area)
{
  AthenaArray<Real> dx1f = pmy_block->dx1f.ShallowCopy();
#pragma simd
  for (int i=il; i<=iu; ++i){
    Real& area_i = area(i);
    Real& dx1_i  = dx1f(i);
    area_i = dx1_i*(pmy_block->dx3f(k));
  }
  return;
}

void Coordinates::Area3Face(const int k, const int j, const int il, const int iu,
  AthenaArray<Real> &area)
{
  AthenaArray<Real> dx1f = pmy_block->dx1f.ShallowCopy();
#pragma simd
  for (int i=il; i<=iu; ++i){
    Real& area_i = area(i);
    Real& dx1_i  = dx1f(i);
    area_i = dx1_i*(pmy_block->dx2f(j));
  }
  return;
}

//--------------------------------------------------------------------------------------
/* \!fn void Coordinates::CellVolume(const int k,const int j,const int il, const int iu,
 *   AthenaArray<Real> &vol)
 * \brief function to compute cell volume    */

void Coordinates::CellVolume(const int k, const int j, const int il, const int iu,
  AthenaArray<Real> &vol)
{
  AthenaArray<Real> dx1f = pmy_block->dx1f.ShallowCopy();
#pragma simd
  for (int i=il; i<=iu; ++i){
    Real& vol_i = vol(i);
    Real& dx1_i = dx1f(i);
    vol_i = dx1_i*(pmy_block->dx2f(j))*(pmy_block->dx3f(k));
  }
  return;
}

//--------------------------------------------------------------------------------------
/* \!fn void Coordinates::CoordinateSourceTerms(
 *   const int k, const int j, AthenaArray<Real> &prim, AthenaArray<Real> &src)
 * \brief function to compute source terms associated with geometry */

void Coordinates::CoordinateSourceTerms(
  const int k, const int j, AthenaArray<Real> &prim, AthenaArray<Real> &src)
{
  return;
}
