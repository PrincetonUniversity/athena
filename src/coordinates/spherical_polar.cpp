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
//! \file spherical_polar.cpp
//  \brief implements functions in class Coordinates for spherical polar coordinates
//======================================================================================

// constructor

Coordinates::Coordinates(Block *pb)
{
  pmy_block = pb;
  int is = pb->is; int js = pb->js; int ks = pb->ks;
  int ie = pb->ie; int je = pb->je; int ke = pb->ke;

// initialize volume-averaged positions and spacing
// x1-direction

  for (int i=is-(NGHOST); i<=ie+(NGHOST); ++i) {
    pb->x1v(i) = 0.75*(pow(pb->x1f(i+1),4) - pow(pb->x1f(i),4))
                     /(pow(pb->x1f(i+1),3) - pow(pb->x1f(i),3));
  }
  for (int i=is-(NGHOST); i<=ie+(NGHOST)-1; ++i) {
    pb->dx1v(i) = pb->x1v(i+1) - pb->x1v(i);
  }

// x2-direction

  if (pb->block_size.nx2 == 1) {
    pb->x2v(js) = 0.5*(pb->x2f(js+1) + pb->x2f(js));
    pb->dx2v(js) = pb->dx2f(js);
  } else {
    for (int j=js-(NGHOST); j<=je+(NGHOST); ++j) {
      pb->x2v(j) = ((sin(pb->x2f(j+1)) - pb->x2f(j+1)*cos(pb->x2f(j+1))) -
                    (sin(pb->x2f(j  )) - pb->x2f(j  )*cos(pb->x2f(j  ))))/
                         (cos(pb->x2f(j  )) - cos(pb->x2f(j+1)));
    }
    for (int j=js-(NGHOST); j<=je+(NGHOST)-1; ++j) {
      pb->dx2v(j) = pb->x2v(j+1) - pb->x2v(j);
    }
  }

// x3-direction

  if (pb->block_size.nx3 == 1) {
    pb->x3v(ks) = 0.5*(pb->x3f(ks+1) + pb->x3f(ks));
    pb->dx3v(ks) = pb->dx3f(ks);
  } else {
    for (int k=ks-(NGHOST); k<=ke+(NGHOST); ++k) {
      pb->x3v(k) = 0.5*(pb->x3f(k+1) + pb->x3f(k));
    }
    for (int k=ks-(NGHOST); k<=ke+(NGHOST)-1; ++k) {
      pb->dx3v(k) = pb->x3v(k+1) - pb->x3v(k);
    }
  }

// Allocate memory for scratch arrays used in integrator, and internal scratch arrays
// Allocate only those local scratch arrays needed for spherical polar coordinates

  int ncells1 = pb->block_size.nx1 + 2*(NGHOST);
  int ncells2 = 1;
  if (pb->block_size.nx2 > 1) ncells2 = pb->block_size.nx2 + 2*(NGHOST);

  face_area.NewAthenaArray(ncells1);   // scratch used in integrator
  cell_volume.NewAthenaArray(ncells1); // scratch used in integrator

  face1_area_i_.NewAthenaArray(ncells1);
  face2_area_i_.NewAthenaArray(ncells1);
  face3_area_i_.NewAthenaArray(ncells1);
  src_terms_i_.NewAthenaArray(ncells1);
  volume_i_.NewAthenaArray(ncells1);

  face1_area_j_.NewAthenaArray(ncells2);
  face2_area_j_.NewAthenaArray(ncells2);
  src_terms_j_.NewAthenaArray(ncells2);
  volume_j_.NewAthenaArray(ncells2);

// compute constant factors used to compute face-areas and cell volumes and store in
// local scratch arrays.  This helps improve performance.

#pragma simd
  for (int i=is-(NGHOST); i<=ie+(NGHOST); ++i){
    face1_area_i_(i) = pb->x1f(i)*pb->x1f(i);
    face2_area_i_(i) = 0.5*(pb->x1f(i+1)*pb->x1f(i+1) - pb->x1f(i)*pb->x1f(i));
    face3_area_i_(i) = 0.5*(pb->x1f(i+1)*pb->x1f(i+1) - pb->x1f(i)*pb->x1f(i));
    volume_i_(i)     = (1.0/3.0)*(pow(pb->x1f(i+1),3) - pow(pb->x1f(i),3));
    src_terms_i_(i)  = face2_area_i_(i)/volume_i_(i);
  }
  if (pb->block_size.nx2 > 1) {
#pragma simd
    for (int j=js-(NGHOST); j<=je+(NGHOST); ++j){
      face1_area_j_(j) = cos(pb->x2f(j)) - cos(pb->x2f(j+1));
      face2_area_j_(j) = sin(pb->x2f(j));
      volume_j_(j)     = cos(pb->x2f(j)) - cos(pb->x2f(j+1));
      src_terms_j_(j)  = (sin(pb->x2f(j+1)) - sin(pb->x2f(j)))/volume_j_(j);
    }
  } else {
    face1_area_j_(js) = cos(pb->x2f(js)) - cos(pb->x2f(js+1));
    face2_area_j_(js) = sin(pb->x2f(js));
    volume_j_(js)     = cos(pb->x2f(js)) - cos(pb->x2f(js+1));
    src_terms_j_(js)  = (sin(pb->x2f(js)) - sin(pb->x2f(js)))/volume_j_(js);
  }

}

// destructor

Coordinates::~Coordinates()
{
// delete scratch arrays used in integrator, and local arrays used internally
  face_area.DeleteAthenaArray();
  cell_volume.DeleteAthenaArray();

  face1_area_i_.DeleteAthenaArray();
  face2_area_i_.DeleteAthenaArray();
  face3_area_i_.DeleteAthenaArray();
  src_terms_i_.DeleteAthenaArray();
  volume_i_.DeleteAthenaArray();

  face1_area_j_.DeleteAthenaArray();
  face2_area_j_.DeleteAthenaArray();
  src_terms_j_.DeleteAthenaArray();
  volume_j_.DeleteAthenaArray();
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
    area_i = face1_area_i_(i)*face1_area_j_(j)*(pmy_block->dx3f(k));
  }
  return;
}

void Coordinates::Area2Face(const int k, const int j, const int il, const int iu,
  AthenaArray<Real> &area)
{
#pragma simd
  for (int i=il; i<=iu; ++i){
    Real& area_i = area(i);
    area_i = face2_area_i_(i)*face2_area_j_(j)*(pmy_block->dx3f(k));
  }
  return;
}

void Coordinates::Area3Face(const int k, const int j, const int il, const int iu,
  AthenaArray<Real> &area)
{
#pragma simd
  for (int i=il; i<=iu; ++i){
    Real& area_i = area(i);
    area_i = face3_area_i_(i)*(pmy_block->dx2f(j));
  }
  return;
}

//--------------------------------------------------------------------------------------
// \!fn 
// \brief

void Coordinates::CellVolume(const int k, const int j, const int il, const int iu,
  AthenaArray<Real> &vol)
{
#pragma simd
  for (int i=il; i<=iu; ++i){
    Real& vol_i = vol(i);
    vol_i = volume_i_(i)*volume_j_(j)*(pmy_block->dx3f(k));
  }
  return;
}

//--------------------------------------------------------------------------------------
// \!fn 
// \brief

void Coordinates::CoordinateSourceTerms(
  const int k, const int j, AthenaArray<Real> &prim, AthenaArray<Real> &src)
{
#pragma simd
  for (int i=(pmy_block->is); i<=(pmy_block->ie); ++i) {
    Real m_tt = prim(IDN,k,j,i)*prim(IM2,k,j,i)*prim(IM2,k,j,i) + prim(IEN,k,j,i);
    Real m_pp = prim(IDN,k,j,i)*prim(IM3,k,j,i)*prim(IM3,k,j,i) + prim(IEN,k,j,i);
    src(IM1,i) = src_terms_i_(i)*(m_tt + m_pp);

    Real m_tr = prim(IDN,k,j,i)*prim(IM2,k,j,i)*prim(IM1,k,j,i);
    src(IM2,i) = (-1.0)*src_terms_i_(i)*(m_tr - src_terms_j_(j)*m_pp);

    Real m_pr = prim(IDN,k,j,i)*prim(IM3,k,j,i)*prim(IM1,k,j,i);
    Real m_pt = prim(IDN,k,j,i)*prim(IM3,k,j,i)*prim(IM2,k,j,i);
    src(IM3,i) = (-1.0)*src_terms_i_(i)*(m_pr + src_terms_j_(j)*m_pt);
  }
  return;
}
