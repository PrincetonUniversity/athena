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

// C headers
#include <math.h>  // pow, trig functions

// Athena++ headers
#include "../athena.hpp"          // macros, Real
#include "../athena_arrays.hpp"   // AthenaArray
#include "../parameter_input.hpp" // ParameterInput
#include "../mesh.hpp"            // MeshBlock

//======================================================================================
//! \file spherical_polar.cpp
//  \brief implements functions in class Coordinates for 3D (r-theta-phi) spherical
//    polar coordinates
//======================================================================================

//--------------------------------------------------------------------------------------
// Coordinates constructor

Coordinates::Coordinates(MeshBlock *pmb, ParameterInput *pin)
{
  pmy_block = pmb;
  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;

// initialize volume-averaged positions and spacing
// x1-direction: x1v = (\int r dV / \int dV) = d(r^4/4)/d(r^3/3)

  for (int i=is-(NGHOST); i<=ie+(NGHOST); ++i) {
    pmb->x1v(i) = 0.75*(pow(pmb->x1f(i+1),4) - pow(pmb->x1f(i),4))
                     /(pow(pmb->x1f(i+1),3) - pow(pmb->x1f(i),3));
  }
  for (int i=is-(NGHOST); i<=ie+(NGHOST)-1; ++i) {
    pmb->dx1v(i) = pmb->x1v(i+1) - pmb->x1v(i);
  }

// x2-direction: x2v = (\int sin[theta] theta dV / \int dV) =
//   d(sin[theta] - theta cos[theta])/d(-cos[theta])

  if (pmb->block_size.nx2 == 1) {
    pmb->x2v(js) = 0.5*(pmb->x2f(js+1) + pmb->x2f(js));
    pmb->dx2v(js) = pmb->dx2f(js);
  } else {
    for (int j=js-(NGHOST); j<=je+(NGHOST); ++j) {
      pmb->x2v(j) = ((sin(pmb->x2f(j+1)) - pmb->x2f(j+1)*cos(pmb->x2f(j+1))) -
                    (sin(pmb->x2f(j  )) - pmb->x2f(j  )*cos(pmb->x2f(j  ))))/
                         (cos(pmb->x2f(j  )) - cos(pmb->x2f(j+1)));
    }
    for (int j=js-(NGHOST); j<=je+(NGHOST)-1; ++j) {
      pmb->dx2v(j) = pmb->x2v(j+1) - pmb->x2v(j);
    }
  }

// x3-direction: x3v = (\int phi dV / \int dV) = dphi/2

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
// Allocate only those local scratch arrays needed for spherical polar coordinates

  int ncells1 = pmb->block_size.nx1 + 2*(NGHOST);
  face_area.NewAthenaArray(ncells1);   // scratch used in integrator
  cell_volume.NewAthenaArray(ncells1); // scratch used in integrator
  face1_area_i_.NewAthenaArray(ncells1);
  face2_area_i_.NewAthenaArray(ncells1);
  face3_area_i_.NewAthenaArray(ncells1);
  src_terms_i_.NewAthenaArray(ncells1);
  volume_i_.NewAthenaArray(ncells1);

  int ncells2 = 1;
  if (pmb->block_size.nx2 > 1) ncells2 = pmb->block_size.nx2 + 2*(NGHOST);
  face1_area_j_.NewAthenaArray(ncells2);
  face2_area_j_.NewAthenaArray(ncells2);
  src_terms_j_.NewAthenaArray(ncells2);
  volume_j_.NewAthenaArray(ncells2);

// compute constant factors used to compute face-areas and cell volumes and store in
// local scratch arrays.  This helps improve performance.

#pragma simd
  for (int i=is-(NGHOST); i<=ie+(NGHOST); ++i){
    face1_area_i_(i) = pmb->x1f(i)*pmb->x1f(i);
    face2_area_i_(i) = 0.5*(pmb->x1f(i+1)*pmb->x1f(i+1) - pmb->x1f(i)*pmb->x1f(i));
    face3_area_i_(i) = 0.5*(pmb->x1f(i+1)*pmb->x1f(i+1) - pmb->x1f(i)*pmb->x1f(i));
    volume_i_(i)     = (1.0/3.0)*(pow(pmb->x1f(i+1),3) - pow(pmb->x1f(i),3));
    src_terms_i_(i)  = face2_area_i_(i)/volume_i_(i);
  }
  if (pmb->block_size.nx2 > 1) {
#pragma simd
    for (int j=js-(NGHOST); j<=je+(NGHOST); ++j){
      face1_area_j_(j) = cos(pmb->x2f(j)) - cos(pmb->x2f(j+1));
      face2_area_j_(j) = sin(pmb->x2f(j));
      volume_j_(j)     = cos(pmb->x2f(j)) - cos(pmb->x2f(j+1));
      src_terms_j_(j)  = (sin(pmb->x2f(j+1)) - sin(pmb->x2f(j)))/volume_j_(j);
    }
  } else {
    face1_area_j_(js) = cos(pmb->x2f(js)) - cos(pmb->x2f(js+1));
    face2_area_j_(js) = sin(pmb->x2f(js));
    volume_j_(js)     = cos(pmb->x2f(js)) - cos(pmb->x2f(js+1));
    src_terms_j_(js)  = (sin(pmb->x2f(js+1)) - sin(pmb->x2f(js)))/volume_j_(js);
  }

}

// destructor

Coordinates::~Coordinates()
{
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
// \!fn void Coordinates::Area1Face(const int k,const int j, const int il, const int iu,
//        AthenaArray<Real> &area)
// \brief functions to compute area of cell faces in each direction

void Coordinates::Area1Face(const int k, const int j, const int il, const int iu,
  AthenaArray<Real> &area)
{
// area1 = r^2 sin[theta] dtheta dphi = r^2 d(-cos[theta]) dphi
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
// area2 = dr r sin[theta] dphi = d(r^2/2) sin[theta] dphi
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
// area3 = dr r dtheta = d(r^2/2) dtheta
#pragma simd
  for (int i=il; i<=iu; ++i){
    Real& area_i = area(i);
    area_i = face3_area_i_(i)*(pmy_block->dx2f(j));
  }
  return;
}

//--------------------------------------------------------------------------------------
// \!fn void Coordinates::CellVolume(const int k,const int j,const int il, const int iu,
//        AthenaArray<Real> &vol)
// \brief function to compute cell volume

void Coordinates::CellVolume(const int k, const int j, const int il, const int iu,
  AthenaArray<Real> &vol)
{
// volume = r^2 sin(theta) dr dtheta dphi = d(r^3/3) d(-cos theta) dphi
#pragma simd
  for (int i=il; i<=iu; ++i){
    Real& vol_i = vol(i);
    vol_i = volume_i_(i)*volume_j_(j)*(pmy_block->dx3f(k));
  }
  return;
}

//--------------------------------------------------------------------------------------
// \!fn void Coordinates::CoordinateSourceTerms(const int k, const int j,
//        AthenaArray<Real> &prim, AthenaArray<Real> &src)
// \brief function to compute coordinate source term

void Coordinates::CoordinateSourceTerms(const int k, const int j,
  AthenaArray<Real> &prim, AthenaArray<Real> &src)
{
// src_1 = < M_{theta theta} + M_{phi phi} ><1/r> = <M_{tt} + M_{pp}> d(r^2/3)d(r^3/3)
// src_2 = -< M_{theta r} - cot[theta]M_{phi phi} ><1/r> 
//         = (<M_{pp}> d(sin[theta])/d(-cos[theta]) - M_{tr}>) d(r^2/3)d(r^3/3)
// src_3 = -< M_{phi r} + cot[theta]M_{phi theta} ><1/r> 
//         = -(<M_{pr} + M_{pt} d(sin[theta])/d(-cos[theta]) >) d(r^2/3)d(r^3/3)
#pragma simd
  for (int i=(pmy_block->is); i<=(pmy_block->ie); ++i) {
    Real m_tt = prim(IDN,k,j,i)*prim(IM2,k,j,i)*prim(IM2,k,j,i) + prim(IEN,k,j,i);
    Real m_pp = prim(IDN,k,j,i)*prim(IM3,k,j,i)*prim(IM3,k,j,i) + prim(IEN,k,j,i);
    src(IM1,i) = src_terms_i_(i)*(m_tt + m_pp);

    Real m_tr = prim(IDN,k,j,i)*prim(IM2,k,j,i)*prim(IM1,k,j,i);
    src(IM2,i) = src_terms_i_(i)*(src_terms_j_(j)*m_pp - m_tr);

    Real m_pr = prim(IDN,k,j,i)*prim(IM3,k,j,i)*prim(IM1,k,j,i);
    Real m_pt = prim(IDN,k,j,i)*prim(IM3,k,j,i)*prim(IM2,k,j,i);
    src(IM3,i) = (-1.0)*src_terms_i_(i)*(m_pr + src_terms_j_(j)*m_pt);
  }
  return;
}
